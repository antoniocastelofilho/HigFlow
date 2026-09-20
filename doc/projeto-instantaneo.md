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

**O instantâneo de facetas.** São 97 sítios, com mapeador próprio
(`sfd_get_domain_mapper`) e uma direção por faceta. O tipo de célula não serve:
faceta não tem `delta` nos DIM eixos do mesmo jeito, e o `sfd_get_stencil` (133
sítios) é topologia, não leitura — fica atrás do backend, como as consultas de
ponto.

Proponho tratá-lo **depois** do de células, e como peça separada: o de células
cobre 233 laços com um tipo que já existe e já tem teste; o de facetas é tipo
novo, e misturar os dois faz a primeira metade esperar pela segunda.

---

## 5. Ordem, e o que dá rede

A restrição real não é tamanho, é cobertura. A suíte de exemplos roda 33 casos;
migrar onde ela não chega é mexer sem rede.

1. **Fechar a lacuna da franja** — `test-mesh-snapshot` com `nps=(1,2,3)`.
   Barato, e é o que torna tudo abaixo confiável.
2. **A produção no `psd_synced_mapper`**, com o guarda contra refino
   pós-montagem. Nenhum
   sítio migrado ainda; a suíte inteira tem de continuar idêntica.
3. **Um oráculo diferencial**: um modo que percorre o domínio pelos dois
   caminhos — árvore e instantâneo — e afirma que centro, delta e id coincidem
   célula a célula, para *todo* domínio que o exemplo construiu. É o que permite
   migrar um laço sem depender do resultado físico para saber se ele está certo.
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

**D3 — A classificação na interface entre árvores** (pendente de ontem).
HiGTree diz `ON_BOUNDARY`, t8code diz dentro. Não bloqueia nada acima: o
instantâneo carrega geometria, não classificação. Mas bloqueia o fechamento de
contorno pelo t8code, que vem depois.
