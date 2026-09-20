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

### O problema mais difícil que eu havia previsto não existe

Eu vinha dizendo que o instantâneo exigiria invalidação a cada adaptação de
malha. **Não há adaptação de malha em tempo de execução.**

```
hig_refine_uniform / hig_split em higflow/src   0 sítios
hig_refine_uniform / hig_split nos exemplos     0 sítios
```

As únicas chamadas de refino vivem dentro da própria HiGTree, e todas no
caminho de **montagem**: leitura de arquivo AMR (`higtree-io.c`), construção de
árvore de contorno (`domain.c:2532`), distribuição inicial
(`higtree-parallel.c`). Depois que `psd_create` retorna, a malha é imutável até
o fim da simulação.

Isso remove a parte cara do projeto. O instantâneo é **produzido uma vez por
domínio, ao fim da montagem, e nunca invalidado.** Não há política de
invalidação a desenhar, não há risco de instantâneo velho, não há custo por
passo de tempo.

Fica um guarda a escrever, não uma política: se algum dia alguém chamar refino
depois da montagem, isso tem de falhar alto, e não silenciosamente devolver
geometria velha.

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

**Recomendação: no `sim_domain`, produzido ao fim da montagem, alcançado por um
acessor.** Não no `higflow_solver`.

O motivo é contagem. O solver declara mais de vinte campos de domínio
(`sdp`, `sdF`, `sdED`, `sdmult`, `sdVisc`, `sdSBnA`, `sdSBnB`, `sdphi`,
`sdEOphi`, `sdEOpsi`, `sdEOnplus`, `sdEOnminus`, `sfdu[DIM]`, `sfdF[DIM]`,
`sfdEOFeo[DIM]`, cada um com seu par `psd`/`psfd`). Pôr um instantâneo ao lado
de cada um duplica essa lista, e cada sítio migrado passa a ter de escolher o
instantâneo certo para o domínio certo — um emparelhamento errado não dá erro
de compilação e dá resultado plausível.

Pendurado no domínio:

```c
const hig_mesh_snapshot *sd_snapshot(sim_domain *sd);   // produz na 1a chamada
```

o emparelhamento deixa de ser possível de errar, o solver não muda de forma, e
o backend que produziu a malha fica irrelevante para quem consome — que é o
ponto da fronteira.

*O que isso custa:* o `sim_domain` passa a ter estado derivado. Aceitável
porque a malha é imutável; seria inaceitável se ela adaptasse.

---

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
2. **O acessor no domínio**, com o guarda contra refino pós-montagem. Nenhum
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

**D1 — Onde o instantâneo vive.** Recomendo no domínio, pelo argumento da
seção 3. A alternativa é no solver, que mantém o `sim_domain` sem estado
derivado ao preço de vinte e poucos campos novos e do emparelhamento manual.

**D2 — Os cinco solvers que nenhum exemplo alcança.** Migrá-los sem rede,
deixá-los por último, ou deixá-los como estão. Eu não migraria sem antes haver
exemplo que os exercite — mas isso os deixa permanentemente no caminho antigo,
e é uma decisão de projeto, não minha.

**D3 — A classificação na interface entre árvores** (pendente de ontem).
HiGTree diz `ON_BOUNDARY`, t8code diz dentro. Não bloqueia nada acima: o
instantâneo carrega geometria, não classificação. Mas bloqueia o fechamento de
contorno pelo t8code, que vem depois.
