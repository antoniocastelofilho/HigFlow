# O instantâneo em produção

Proposta de projeto, não de implementação. O objetivo é decidir **onde o
`hig_mesh_snapshot` vive, quem o produz e o que o invalida** antes de encostar
nos sítios de consulta — porque varrer sem isso produz código que fica verde e
errado, e cinco dos arquivos com mais sítios são solvers que nenhum exemplo
alcança.

Tudo abaixo foi medido em `higflow/src` hoje, 20 de setembro de 2026, com o
`src_hugo` excluído da contagem (ele não entra na compilação).

---

## 1. O que a medição mudou no plano

### A adaptação existe, e reconstrói em vez de mutar

Eu havia afirmado aqui que **não há adaptação de malha em tempo de execução**.
Estava errado, e o erro foi de medição: grepei `higflow/examples`, diretório que
não existe — os exemplos são `higflow/example2d_*` e `example3d_*`. O
`2>/dev/null` engoliu o *No such file or directory* e eu li zero resultados como
ausência de adaptação.

Medido de novo, nos 14 exemplos:

```
higflow/src                       0 sítios de hig_refine_uniform / hig_split
example2d_ElectroOsmotic          2 sítios
example2d_DynamicMeshAdapt        2 sítios
```

O que a leitura de cada um mostra é que **a premissa de projeto sobrevive, por um
motivo diferente do que eu tinha escrito**:

- `example2d_ElectroOsmotic`, sítio 1: refina um `hig_clone` da raiz, usado só
  para escrever VTK e destruído em seguida. A malha viva não é tocada.
- `example2d_ElectroOsmotic`, sítio 2: está em `adapt_mesh_and_update_solver`,
  que refina a malha viva — e cuja **única chamada está comentada** (linha 364).
- `example2d_DynamicMeshAdapt`: o caminho vivo é `higflow_rebuild_with_amr`, que
  constrói uma árvore adaptada nova, chama `lb_calc_partition` e **reconstrói** o
  solver sobre ela.

Ou seja: a adaptação em tempo de execução **existe**, e ela funciona
reconstruindo o domínio, não mutando o que está vivo. Nenhum caminho alcançável
refina uma árvore que já pertence a um `sim_domain`.

Para o instantâneo isso é melhor do que a premissa errada que eu tinha: como a
adaptação passa por `psd_create` / `psd_synced_mapper` outra vez, **o instantâneo
é reproduzido sozinho**. Continua não havendo política de invalidação a desenhar.

**Mas existe uma mina, e ela é o motivo de o detector não ser supérfluo.** Duas
funções de refino *no lugar* estão escritas e dormentes:

```
higflow_refine_tree_inplace      definida, nunca chamada
adapt_mesh_and_update_solver     definida, chamada comentada
```

Nenhuma das duas passa pela API do domínio, então a guarda de
`sd_add_higtree` não as vê. Se alguém as acordar, o instantâneo fica velho **em
silêncio** e os laços migrados passam a ler a célula errada. Quem detecta esse
caso é o `sd_snapshot_is_current`, e quem o acordar deve chamar
`sd_compute_snapshot` logo depois — está anotado nos dois lugares.

### A lacuna da franja está fechada, e o que ela revelou

O `hms_from_domain` indexa por `mp_lookup(sd_get_domain_mapper(sd),
hig_get_cid(c))` e rejeita índice fora de `[0, n)`, onde `n` exclui a franja
pela cláusula C6. Essa afirmação só estava exercitada em `np=1`, e com a franja
montada à mão pelo próprio teste.

`test-mesh-snapshot-parallel` fecha isso: monta o domínio pelo caminho de
produção (`lb_calc_partition`, `psd_create`, `psd_synced_mapper`) e afirma, por
redução global, em np = 1, 2 e 3 e nos dois DIM. **A convenção se sustenta** —
os locais ficam em `[0, n)` e a franja vem depois, em todas as configurações.

Mas o exercício revelou algo que muda como se deve ler essa garantia.

**A identidade de índice é verdadeira por construção.** O `_psd_setmapper`
numera os locais com
`mp_assign_from_celliterator(m, sd_get_domain_celliterator(sd), 0)` — o *mesmo*
iterador que o `hms_from_domain` percorre. O id do mapeador é, por definição, o
contador do percurso. Medido: trocar `mp_lookup(m, hig_get_cid(c))` por
`n_visitadas++` dentro do `hms_from_domain` deixa a suíte **inteira verde**, nos
dois DIM e nos três np. Não é oráculo fraco — os dois são a mesma função.

O que resta a proteger, portanto, não é a aritmética do instantâneo: é a
**convenção de numeração da produção**. E aí o teste discrimina. Numerando a
franja primeiro no `_psd_setmapper`, o teste serial passa nos seus 3 casos (ele
imita a convenção à mão) e o paralelo reprova em np=2 e np=3, nos dois DIM.

Isso tem consequência de projeto: o passo 3 da seção 5 — o oráculo diferencial —
**não pode** comparar o instantâneo contra o mesmo percurso que o gerou, porque
isso é tautologia. Ele tem de alcançar a célula por outra via (busca por ponto)
e conferir que o id leva à linha certa.

## 2. O terreno

```
laços de célula   (for it = sd_get_domain_celliterator)      233
sítios de faceta  (sfd_get_domain_facetiterator)              97

sd_get_domain_celliterator   255      hig_get_center   250
hig_get_cid                  207      hig_get_delta    144
sfd_get_stencil              133      sd_get_stencil    67
sd_get_cell_with_point        18
                                      TOTAL          1074
```

O laço canônico — e são 233 variações dele — é exatamente a forma do
instantâneo:

```c
for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
    hig_cell *c = higcit_getcell(it);
    int clid = mp_lookup(m, hig_get_cid(c));   //  ->  i
    Point center; hig_get_center(c, center);   //  ->  s->center[i*DIM + d]
    real val = ns->problem->pressure(center, ns->par.t);
    dp_set_value(ns->dpp, clid, val);
}
```

Dos 306 `mp_lookup`, **183 são literalmente o par `mp_lookup(m,
hig_get_cid(c))`**, e `psd_get_local_id` não é usado em lugar nenhum. O
mapeador é sempre o mesmo objeto que o `hms_from_domain` já consulta
(`sd_get_domain_mapper`), então a identidade de índice vale por construção, não
por coincidência — uma vez confirmada com franja.

---

## 3. Onde o instantâneo vive

**Recomendação: no `sim_domain`, produzido ao fim do `psd_synced_mapper`.** Não
no `higflow_solver`, e não na primeira leitura.

### Por que no domínio

O argumento é de contagem. O solver declara mais de vinte campos de domínio —
`sdp`, `sdF`, `sdED`, `sdmult`, `sdVisc`, `sdSBnA`, `sdSBnB`, `sdphi`, os quatro
de eletro-osmose, mais `sfdu[DIM]`, `sfdF[DIM]`, `sfdEOFeo[DIM]`, cada um com
seu par `psd`/`psfd`. Pôr um instantâneo ao lado de cada um duplica essa lista,
e cada laço migrado passa a ter de escolher o instantâneo certo para o domínio
certo. **Emparelhamento errado não dá erro de compilação e dá resultado
plausível.** Pendurado no domínio, esse erro deixa de ser possível de cometer.

O ganho colateral é o da fronteira: quem consome deixa de saber qual backend
produziu a malha.

### Por que ansioso, e não na primeira leitura

Uma versão anterior desta seção propunha `sd_snapshot(sd)` produzindo na
primeira chamada. **Está errado**, e pelo mesmo motivo que o teste de franja
expôs: o instantâneo depende do mapeador estar atribuído, e quem o atribui é o
`_psd_setmapper`, dentro do `psd_synced_mapper`. Produção preguiçosa deixa essa
ordem implícita — chamada antes, ela devolve um instantâneo indexado por lixo, e
o sintoma não é falha, são números errados.

A medição mostra que há um ponto único e não ambíguo para produzi-lo:

```
psd_create                                   19
psd_synced_mapper                            19
mp_assign_from_celliterator em higflow/src    1
```

Um para um, e o solver não numera domínio nenhum por conta própria — delega.
Então o instantâneo nasce onde o mapeador acabou de ficar pronto, e a ordem
deixa de ser uma coisa que alguém precise lembrar.

Junto vai o guarda da seção 1: refino depois da montagem tem de falhar alto, e
não devolver geometria velha em silêncio.

### O custo, e o que mudaria a decisão

O `sim_domain` passa a carregar estado derivado, e é a HiGTree que passa a
mantê-lo. Isso só é aceitável porque **a malha é imutável depois da montagem**.
Se a adaptação em tempo de execução entrar algum dia, esta é a primeira decisão
a rever.

A alternativa legítima é manter a HiGTree sem memória — uma biblioteca de malha
que responde e não lembra. Não é errada; o preço dela são os vinte e poucos
campos novos no solver e o emparelhamento manual, que é pagar no lugar onde o
erro é silencioso.

## 4. O que ainda não existe

**~~O instantâneo de facetas~~ — feito.** `hig_facet_snapshot` guarda a **caixa
da célula** da faceta, mais `dim` e `dir` por faceta; centro e tamanho saem daí
pelas mesmas contas do `hig_get_facet_center` e do `hig_get_facet_delta`.

Guardar a caixa não foi simetria com o caso das células — foi medição. Dos 97
laços de faceta, 95 leem apenas centro, tamanho e id; os outros 2 chamam o
buffer de resíduo, que precisa da caixa da **célula**. Fosse só pelos 95,
guardar centro e tamanho prontos bastaria.

Produção ansiosa no `psfd_synced_mapper`, oráculo diferencial próprio
(`sfd_snapshot_verify`, localizando por `sfd_get_facet_with_point`), e a mesma
flag `HIGTREE_VERIFY_SNAPSHOT` o liga.

O `sfd_get_stencil` (133 sítios) continua fora: é topologia, não leitura, e fica
atrás do backend como as consultas de ponto.

**A caixa da célula no instantâneo — ~~decisão nova~~ feita.** O instantâneo
guardava centro e delta. Agora guarda **a caixa** (`low` e `high`), e centro e
delta saem dela por derivação, com exatamente as mesmas contas do
`hig_get_center` e do `hig_get_delta`.

O motivo é que a caixa é o primário: é o que a célula guarda, e os outros dois
são derivados dela. Guardando a caixa, toda leitura sai bit a bit igual à de
hoje; guardando centro e delta, reconstruir `low = centro − delta/2` é exato em
álgebra e **não é exato em ponto flutuante**.

Isso não era preciosismo, e agora está **medido** em vez de afirmado. Os 19
laços cobertos do `hig-flow-io.c` usam `c->lowpoint` e `c->highpoint` como
*pontos de interpolação*, e o resultado vai para o VTK que a suíte compara. A
perda da reconstrução:

```
centro − delta/2 ≠ low      1 de 3 células por direção na malha não diádica
                            ~0,6% de 200 mil caixas aleatórias
```

O caso `reconstruir_o_canto_nao_seria_exato` afirma essa perda, e é o único
lugar que guarda a decisão: se alguém "simplificar" o instantâneo de volta para
centro e delta, ele é o que acusa.

Uma coisa que eu **não** consegui prender, e vale registrada: a escolha entre
formas algebricamente equivalentes do centro. `(lo+hi)/2`, `lo/2+hi/2` e
`lo+(hi−lo)/2` dão o mesmo bit em 200 mil caixas aleatórias — dividir por dois é
exato em binário. Não há ali o que discriminar.

Custo: dois arranjos de `n * DIM` em vez de dois (mesma memória), e o produtor
do t8code passou a preencher a caixa a partir do centroide e do nível.

## 4b. Produção por rank — tentada, e o que ela ensinou

**Não entrou**, e o motivo é uma restrição da HiGTree que eu não conhecia.

A ideia era: a floresta nasce em `COMM_WORLD`, o t8code já a reparte, e cada
rank materializa só os seus elementos. O conjunto local de um rank é um
intervalo da curva de preenchimento, e uma árvore hig é uma **caixa** refinada —
então eu ia representar a região com uma árvore **com buracos**, usando
`hig_refine_empty` (que cria o vetor de filhos zerado) e preenchendo só os
compartimentos do rank.

**Árvore com buracos não é navegável.** `hig_get_cell_coords_of_point`
desreferencia `cell->children[...]` para ler as caixas dos filhos; com filhos
nulos, é SEGV. E `hig_get_cell_with_point`, que a chama, ainda faz
`cell = cell->children[p]` sem testar nulo. Toda a localização por ponto —
inclusive a do fechamento de contorno e a do oráculo — depende dessas duas.

**E a minha premissa sobre o `lbal` estava errada.** Eu tinha escrito que ele
produz árvores com buracos. Ele não produz: `hig_refine_empty` ali serve para
*alocar o vetor*, e em seguida o laço preenche **todos** os `tree_size` filhos.
Cada árvore de saída é uma caixa **completa** — e é por isso que um domínio tem
várias (medido: até três por rank com np=3), com um passo de fusão entre elas.

**O que a próxima tentativa tem de fazer:** decompor o conjunto de elementos do
rank em **caixas completas** e emitir uma árvore por caixa. Isso é um problema
de cobertura por retângulos sobre os índices possuídos, não uma adaptação do que
já existe.

*Verificado durante a tentativa, e vale guardar:* o caminho local funciona —
em np=2 cada rank materializou 35 elementos e viu 8 ghosts. O que quebrou foi a
navegação, não a produção. E pedir `set_adapt`, `set_ghost` e `set_partition` no
mesmo commit do t8code derruba com np>1 e funciona em np=1, que é o modo como
esse tipo de defeito se esconde.

## 5. Ordem, e o que dá rede

A restrição real não é tamanho, é cobertura. A suíte de exemplos roda 33 casos;
migrar onde ela não chega é mexer sem rede.

1. **Fechar a lacuna da franja** — `test-mesh-snapshot` com `nps=(1,2,3)`.
   Barato, e é o que torna tudo abaixo confiável.
2. **A produção no `psd_synced_mapper`**, com o guarda contra refino
   pós-montagem. Nenhum
   sítio migrado ainda; a suíte inteira tem de continuar idêntica.
3. **O oráculo diferencial** — ~~a fazer~~ **feito**. `sd_snapshot_verify`
   localiza cada célula **por ponto** e exige que ela volte com o mesmo id e a
   mesma geometria. Com `HIGTREE_VERIFY_SNAPSHOT` no ambiente, todo domínio que
   o programa montar passa por ele, e divergência aborta. Os 33 casos da suíte
   de exemplos passam com o modo ligado.

   Ele **não** percorre o iterador, e isso não é detalhe: o instantâneo foi
   preenchido pelo iterador indexado pelo mapeador, e o mapeador saiu desse mesmo
   iterador. Conferir um contra o outro mede a aritmética contra ela própria.
4. **Migrar os laços cobertos por exemplo**, um arquivo por vez, com a suíte
   entre cada um.
5. **Os solvers sem exemplo, por último** — ou não. Ver a decisão D3.

O caminho de árvore **não sai**. O instantâneo é aditivo: um laço migrado lê do
arranjo, os demais seguem lendo da árvore, e os dois descrevem a mesma malha
enquanto o oráculo do item 3 passar.

---

## 6. O que é decisão sua

**D1 — Onde o instantâneo vive, e quando nasce.** Recomendo no domínio,
produzido ao fim do `psd_synced_mapper`, pelo argumento da seção 3. A
alternativa é no solver, que mantém o `sim_domain` sem estado derivado ao
preço de vinte e poucos campos novos e do emparelhamento manual.

**D2 — Os cinco solvers que nenhum exemplo alcança.** Migrá-los sem rede,
deixá-los por último, ou deixá-los como estão. Eu não migraria sem antes haver
exemplo que os exercite — mas isso os deixa permanentemente no caminho antigo,
e é uma decisão de projeto, não minha.

**D3 — A classificação na interface entre árvores — decidido.** Vale a
semântica do t8code: *sobre o contorno* significa que o domínio **termina** ali.
Interface entre blocos do mesmo domínio não é contorno, porque há malha dos dois
lados.

A investigação mudou o caráter da decisão. Não eram duas convenções
defensáveis — o `cell_find_in_center` para no **primeiro** higtree cuja caixa
contém o ponto (`break`, "assume domains do not overlap") e nunca pergunta se
outro bloco continua o domínio. Isso faz a resposta depender da **ordem das
árvores**, que é o que a cláusula C8 proíbe para localização.

E não é hipotético: medido, com `np=3` um domínio do `example2d_Newt` chega a
**três** higtrees, então as interfaces internas existem em toda execução
paralela.

O efeito hoje é **benigno**, e rastreei até o fim: `ON_BOUNDARY` faz o
`get_stencil` procurar condição de contorno — Dirichlet no centro, Neumann,
Dirichlet — e numa interface interna nenhuma casa, então cai no caminho normal
de interpolação. O preço é busca desperdiçada na camada de facetas daquele
plano.

**O que não fica garantido** é que nenhuma condição registrada coincida com uma
interface interna; se coincidir, seria aplicada num ponto interno. Corrigir
exige perguntar se o domínio *continua* além do ponto — uma sonda com epsilon —
e isso merece verificação própria em vez de vir de carona nesta decisão. Fica
escrito na C14 como desvio conhecido do MTree, não como convenção alternativa.
