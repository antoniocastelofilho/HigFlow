# example2d_DynamicMeshAdapt.TestInitialization

Teste de adaptação dinâmica de malha com inicialização por interface analítica.

## Problemas corrigidos

### 1. Adaptação inicial (step 0) com domínios particionados inconsistentes

**Sintoma:** A malha refinada no passo 0 não era aplicada corretamente — a simulação usava a malha original sem refinamento ou apresentava comportamento inconsistente.

**Causa:** O bloco de adaptação do step 0 fazia `sd_add_higtree(ns->sdp, root)` para trocar a árvore da malha **depois** que `create_initialize_all_domains()` já havia criado os domínios particionados (`psdp`, `psfdu`, `psdmult`) em cima da malha original. Os domínios particionados continuavam referenciando a árvore velha — a árvore refinada ficava "órfã" nos domínios serial e os particionados operavam na topologia errada.

**Solução:** Usar a mesma abordagem da adaptação dinâmica (step % 5): criar um solver novo `ns2`, construir a árvore adaptada via `higflow_make_adapted_tree_params` (que clona a árvore velha e refina na cópia), adicionar nos domínios seriais de `ns2`, particionar `ns2` com `partition_graph` vazio (válido apenas em serial), criar stencils, e trocar `ns = ns2`. Desta forma DPs, stencils e domínios particionados são criados do zero em cima da árvore já refinada.

```c
higflow_solver *ns2 = higflow_create();
higflow_load_data_file_names(argc, argv, ns2);
higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
higflow_set_external_functions(ns2, ...);
higflow_create_domain(ns2, cache, order_center);
higflow_create_domain_multiphase(ns2, cache, order_center, ...);

hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);
sd_add_higtree(ns2->sdp, root);
sd_add_higtree(ns2->sdF, root);
sd_add_higtree(ns2->ed.mult.sdmult, root);

partition_graph *pg = pg_create(MPI_COMM_WORLD);
pg_set_fringe_size(pg, 5);
load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
lb_destroy(lb);
higflow_create_partitioned_domain(ns2, pg, order_center);
higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);
higflow_create_stencil(ns2);
higflow_create_stencil_multiphase(ns2);

ns2->par = ns->par;
ns2->contr = ns->contr;
ns = ns2;
```

**Arquivo:** `ns-example-2d.c`, bloco `ADAPT INICIAL BASEADO NA INTERFACE ANALÍTICA`.

### 2. `num_levels` manual e vetor `thresholds` duplicado

**Sintoma:** O número de níveis de refinamento era definido manualmente (`num_levels = 2`) enquanto o vetor `thresholds` tinha 3 valores — código morto e risco de inconsistência.

**Causa:** Dois call sites da função `higflow_make_adapted_tree_params` definiam vetores `thresholds` diferentes, e `num_levels` era passado como parâmetro separado, sem relação com o tamanho real do vetor.

**Solução:** Os thresholds foram unificados em um único array global `REFINE_THRESHOLDS[]` em `ns-example-2d.c`, terminado com sentinela `-1.0`. A função `higflow_make_adapted_tree_params` passou a calcular `num_levels` internamente percorrendo o array até a sentinela, eliminando o parâmetro redundante.

**Arquivos:** `ns-example-2d.c` (array global), `mesh_adapt_function.c` (função adaptada).

### 3. `partition_graph` vazio nas adaptações (impede paralelismo)

**Sintoma:** A simulação só funcionava em serial (1 processo MPI). Em paralelo, os resultados eram incorretos.

**Causa:** Ambos os blocos de adaptação (step 0 e step % 5) criavam um `partition_graph` novo e vazio com `pg_create(MPI_COMM_WORLD)` e descartavam o `load_balancer` sem usá-lo. Esse pg vazio não tinha informação de vizinhos MPI. Quando `dp_sync()` era chamado, ele consultava `pg->filtered_neighbors` (vazio) e não enviava/recebia nada — células fantasmas (ghost/fringe) nunca eram sincronizadas.

**Tentativa 1 — Reuso do `partition_graph` velho:** `psd_get_partition_graph(ns->psdp)`. Falhou porque o pg armazena ponteiros para as árvores antigas em `pg->tree_props` e `pg->neighbors`. Após `hig_clone` + `hig_refine_uniform`, os ponteiros das novas árvores são diferentes. Quando `psd_set_local_domain` tenta consultar `pg->tree_props` com os novos ponteiros, não encontra nada (`NULL`). E `_create_filtered_neighbors` tenta casar os ponteiros de árvores vizinhas antigas contra a `tree_map` nova — não casa nada, produzindo `filtered_neighbors` vazio.

**Tentativa 2 — Refino in-place da árvore original:** `higflow_refine_tree_inplace` modifica a árvore diretamente sem clonar, preservando os ponteiros. Funciona para o step 0 mas quebra no step 5 porque a árvore já refinada in-place, quando clonada pelo `higflow_make_adapted_tree_params`, produz uma árvore cuja interpolação para a malha nova corrompe a memória (heap corruption).

**Tentativa 3 — Load balancer com árvores clonadas:** Usar `lb_add_input_tree` + `lb_calc_partition` para reparticionar as árvores adaptadas. Falhou com `MPI_ERR_TRUNCATE` porque o load balancer espera árvores completas, não parciais (cada rank só tem seu pedaço local do domínio).

**Solução final (atual):** `partition_graph` vazio (`pg_create` + `pg_set_fringe_size` sem `lb_calc_partition`). Válido apenas em serial (1 MPI rank). Em paralelo, `dp_sync` não funciona porque não há vizinhos configurados. A implementação de paralelismo na adaptação dinâmica de malha ainda está pendente.

```c
partition_graph *pg = pg_create(MPI_COMM_WORLD);
pg_set_fringe_size(pg, 5);
load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
lb_destroy(lb);
```

**Arquivo:** `ns-example-2d.c`, ambos os blocos de adaptação.

### 4. Filtro de `fracvol` com `printf` por célula

**Sintoma:** Milhares de linhas "Step 5" no log da simulação, poluindo a saída e degradando performance.

**Causa:** O filtro empírico de `fracvol` na função `higflow_interpolate_viscosity` imprimia `printf("Step %d\n", ns->par.step)` para cada célula da malha nova durante a interpolação.

**Solução:** Mantido conforme decisão do usuário (o filtro é necessário para o sharpening da fração de volume, mas o `printf` pode ser removido se o volume de saída for problemático).

**Arquivo:** `ns-example-2d.c`, função `higflow_interpolate_viscosity`.

### 5. Velocidade da tampa deslizante (lid-driven cavity)

**Mudança:** A velocidade da parede superior (bc1) foi dobrada de `u=1.0` para `u=2.0` para aumentar a intensidade do escoamento na cavidade. A condição inicial da velocidade também foi atualizada de `u=y` para `u=2.0*y` para manter consistência com a condição de contorno.

**Arquivos:** `ns-user-functions-newtonian-gn.c`, funções `get_boundary_velocity` (case 1) e `get_velocity`.

## Status atual

- **Serial (NP=1):** Funcionando. Adaptação em step 0 e step 5 com conservação de massa ~0.05%.
- **Paralelo (NP>1):** Não funciona. O `partition_graph` vazio não possui vizinhos MPI configurados, impedindo `dp_sync` de comunicar células fantasmas entre ranks. Pendente de implementação futura.

## Otimizações futuras (pendentes)

### 6. Otimizar a busca de sementes no coarsening (O(N²) → O(N))

**Problema:** Em `mesh_adapt_function.c`, o coarsening (linhas 169-251) faz para cada célula pai uma verificação de distância contra TODAS as sementes da interface:

```c
// Para cada child de cada parent, busca distância mínima para TODAS as sementes
for(int s=0; s<seed_count; s++) {
    real d2 = dist_sq_c(center, seeds[s].center);
    if (d2 < min_d2) min_d2 = d2;
}
```

Para 1000 sementes × 500 pais × 4 filhos × 4 passes = **8M de cálculos de distância**. Complexidade O(seed_count × parent_count) = O(N²) na região da interface.

**Solução proposta:** Spatial hashing (grid de bins 2D):

```c
// Divide o domínio em bins de tamanho max(thresholds) + coarsen_hys
// Cada bin armazena índices das sementes naquela região
// Para cada child, consulta apenas sementes no mesmo bin + bins vizinhos (9 bins no máximo)
// Reduz de O(N²) para O(N) com overhead de memória O(seed_count)
```

**Arquivo:** `mesh_adapt_function.c`, função `higflow_make_adapted_tree_params`, etapa de coarsening.

### 7. Fundir as 8 interpolações do AMR em 1 único loop

**Problema:** No bloco AMR (`ns-example-2d.c:828-842`), 8 funções de interpolação iteram a malha completa sequencialmente:

1. `higflow_interpolate_velocity` — itera todos os facets
2. `higflow_interpolate_viscosity` — itera todas as células (2× `compute_value_at_point`)
3. `higflow_interpolate_density` — itera todas as células
4. `higflow_interpolate_pressure` — itera todas as células
5. `higflow_interpolate_bc_for_velocity` — itera células de contorno
6. `higflow_interpolate_bc_for_pressure` — itera células de contorno

Cada função abre seu próprio `celliterator`/`facetiterator`, faz `mp_lookup` por célula, e fecha o iterador. O overhead de criar/destruir iteradores e a falta de localidade de cache tornam isso ~30% mais lento que uma única passagem.

**Solução proposta:** Função `higflow_interpolate_all` que faz uma única iteração sobre a malha nova:

```c
void higflow_interpolate_all(higflow_solver *ns, higflow_solver *ns2) {
    // Única iteração sobre células da malha nova
    higcit_celliterator *it = sd_get_domain_celliterator(sdm2);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp2, hig_get_cid(c));
        Point ccenter;
        hig_get_center(c, ccenter);

        // Interpolar TUDO de uma vez
        real fracvol = compute_value_at_point(sdm1, ccenter, ccenter, 1.0,
                                               ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        real visc    = compute_value_at_point(sdm1, ccenter, ccenter, 1.0,
                                               ns->ed.mult.dpvisc, ns->ed.mult.stn);
        // ... densidade, pressão, etc.

        // Escrever tudo
        dp_set_value(ns2->ed.mult.dpfracvol, clid, fracvol);
        dp_set_value(ns2->ed.mult.dpvisc, clid, visc);
        // ...
    }
    // Um único dp_sync por propriedade no final
}
```

**Ganho esperado:** ~30% de redução no tempo do bloco AMR.

**Arquivos:** `ns-example-2d.c`, funções `higflow_interpolate_*`, bloco AMR.

### 8. AMR incremental (estilo Basilisk)

**Problema:** Atualmente o AMR **destrói e recria o solver inteiro** a cada adaptação. Isso inclui:

- `higflow_destroy(ns)` — libera domínios, stencils, DPs, KSP do PETSc
- `higflow_create()` — aloca novo solver
- `higflow_create_domain()` + `higflow_create_domain_multiphase()` — novos domínios
- `higflow_create_partitioned_domain()` — particionamento
- `higflow_create_stencil()` — stencils
- `higflow_create_distributed_properties()` — novos DPs
- `higflow_create_solver()` — novo KSP do PETSc
- 8 interpolações full-mesh para transferir dados

Custo total ≈ **5-10× o de um step normal do solver**. No Basilisk, o AMR custa tipicamente <5% do passo de integração.

**Solução proposta estilo Basilisk:**

1. **Refinar/coarsen in-place** — usar `hig_refine_uniform` / `hig_merge_children` diretamente na árvore ativa sem clonar. A função `higflow_refine_tree_inplace` em `mesh_adapt_function.c` já existe parcialmente implementada mas não é usada.

2. **Atualizar só células afetadas** — após refinar/coarsen, apenas as novas células e as removidas precisam de inicialização. Transferir valores por interpolação local (aresta a aresta), não full-mesh.

3. **AMR todo step mas barato** — em vez de adaptar 100% da interface a cada N steps, adaptar 1-2% das células por step (as que cruzaram o limiar de refinamento). O custo amortizado é desprezível e a malha fica sempre ótima.

4. **Reusar KSP do PETSc** — o `KSP` pode ser reinicializado (`KSPSetOperators`) em vez de destruído/recriado. As matrizes do PETSc podem ser re-preenchidas com nova topologia sem realocar o objeto KSP.

5. **Critério local wavelet** (como Basilisk) — em vez de busca global por sementes, usar o gradiente da fração volumétrica como estimador de erro por célula:

```c
// Para cada célula: se |∇f| > threshold e nível < max_level → refinar
// Se |∇f| < threshold/coarsen_factor e nível > min_level → coarsen
```

**Ganho esperado:** Redução de ~80% no custo do AMR (de 5-10× steps para <1× step).

**Arquivos:** `mesh_adapt_function.c`, `ns-example-2d.c`.

### 9. (Bônus) Corrigir `create_velocity_copy` — ponteiro não inicializado

**Problema:** Em `utilities-2d.c:1031`, a função `create_velocity_copy` declara `distributed_property **u_copy` sem `malloc`, e depois escreve em `u_copy[dim]`. Isso é **undefined behavior** (stack corruption).

```c
distributed_property **create_velocity_copy(higflow_solver *ns){
    distributed_property **u_copy;  // NUNCA alocado!
    for(int dim=0; dim<DIM; dim++){
        u_copy[dim] = psfd_create_property(ns->psfdu[dim]);  // escreve em lixo
    }
    return u_copy;
}
```

**Solução:** Adicionar `malloc`:

```c
distributed_property **u_copy = malloc(DIM * sizeof(distributed_property*));
```

**Arquivo:** `utilities-2d.c`, função `create_velocity_copy`.
