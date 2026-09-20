#ifndef HIG_MESH_SNAPSHOT_H
#define HIG_MESH_SNAPSHOT_H

// =============================================================================
// A FRONTEIRA DAS CONSULTAS
//
// Este tipo existe para tirar a navegacao de arvore de dentro do laco quente.
//
// MEDIDO em higflow/src: das 907 chamadas de consulta que estao dentro de laco
// sobre celulas, 777 -- 86% -- sao LEITURA DE DADO:
//
//     sd_get_domain_celliterator   237
//     hig_get_center               228
//     hig_get_cid                  190
//     hig_get_delta                122
//
// Nenhuma delas precisa da arvore.  Precisam de centro, tamanho e indice, por
// celula local.  Postos num arranjo plano preenchido UMA VEZ por producao, o
// laco quente vira leitura de memoria contigua: sem indirecao, sem despacho, e
// -- o que importa para a substituicao -- sem depender de QUAL backend produziu
// a malha.
//
// As outras 130 chamadas (14%) sao topologia de verdade -- `sfd_get_stencil` e
// `sd_get_cell_with_point` -- e continuam atras do backend.  Nao se tenta
// achatar o que e' consulta.
//
// POR QUE ISSO E' A FRONTEIRA, e nao uma otimizacao.  O adaptador do t8code
// hoje materializa a floresta dele como octree de ponteiros a cada producao,
// O(folhas), so' para que as consultas possam le-la.  Com o instantaneo, um
// backend preenche estes arranjos DIRETO da sua propria estrutura, e a
// materializacao deixa de ser necessaria.  E' a diferenca entre traduzir a cada
// passo e produzir uma vez.
//
// O INDICE E' A IDENTIDADE.  `center[i*DIM + d]` e' a coordenada d do centro da
// celula local i, e i e' o mesmo id que `mp_lookup(m, hig_get_cid(c))` devolve.
// Isso apaga as 190 chamadas de identificador do laco: o id e' o proprio
// contador.
// =============================================================================

#include "types.h"
#include "domain.h"

#ifdef __cplusplus
extern "C" {
#endif

//! O instantaneo guarda A CAIXA, e nao centro e tamanho.
//!
//! POR QUE A CAIXA E' O PRIMARIO.  E' o que a celula guarda: `hig_get_center` e
//! `hig_get_delta` DERIVAM dela, por `(low+high)/2` e `high-low`.  Guardando a
//! caixa e derivando com as mesmas contas, toda leitura sai bit a bit igual a'
//! de hoje.  O contrario nao vale: reconstruir `low = centro - delta/2` e' exato
//! em algebra e NAO em ponto flutuante.
//!
//! Isso nao e' preciosismo.  19 lacos cobertos do `hig-flow-io.c` usam
//! `c->lowpoint` e `c->highpoint` como PONTOS DE INTERPOLACAO, e o resultado vai
//! para o VTK que a suite compara com referencia.  Com centro e delta guardados,
//! migrar esses lacos mudaria a saida no ultimo bit; com a caixa guardada, nao
//! muda nada.
typedef struct hig_mesh_snapshot {
    int   n;          //!< numero de celulas LOCAIS (a franja nao entra, ver C6)
    int   dim;        //!< DIM com que foi construido, para conferencia
    real *low;        //!< n * DIM, indexado por [i*DIM + d]
    real *high;       //!< n * DIM, indexado por [i*DIM + d]
} hig_mesh_snapshot;

//! \brief Preenche um instantaneo a partir de um sim_domain (backend MTree).
//!
//! Percorre o iterador de dominio -- so' as celulas locais, como manda a
//! clausula C6 -- e indexa por `mp_lookup(m, hig_get_cid(c))`, de modo que o
//! indice do arranjo E' o id local.
hig_mesh_snapshot *hms_from_domain(sim_domain *sd);

//! \brief O centro da celula local `i`, no formato que os lacos ja' usam.
//!
//! Existe para que a migracao de um laco troque UMA linha --
//! `hig_get_center(c, center)` por `hms_center(hms, i, center)` -- em vez de
//! espalhar aritmetica de indice por 250 sitios.
static inline void hms_center(const hig_mesh_snapshot *s, int i, Point p) {
    // A MESMA conta do `hig_get_center`: POINT_ADD seguido de POINT_DIV_SCALAR.
    // Escrita de outro jeito -- `low/2 + high/2`, por exemplo -- ela deixaria de
    // ser bit a bit igual, e o ponto todo de guardar a caixa se perderia.
    for (int d = 0; d < DIM; d++) p[d] = (s->low[i * DIM + d] + s->high[i * DIM + d]) / 2.0;
}

//! \brief O canto inferior da celula local `i`, exatamente como a celula o guarda.
static inline void hms_low(const hig_mesh_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->low[i * DIM + d];
}

//! \brief O canto superior da celula local `i`.  Ver `hms_low`.
static inline void hms_high(const hig_mesh_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->high[i * DIM + d];
}

//! \brief O tamanho da celula local `i`.  Ver `hms_center`.
static inline void hms_delta(const hig_mesh_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->high[i * DIM + d] - s->low[i * DIM + d];
}

//! \brief Aloca um instantaneo vazio para `n` celulas.  Para um backend que
//! preenche os arranjos direto, sem passar por sim_domain.
hig_mesh_snapshot *hms_create(int n);

void hms_destroy(hig_mesh_snapshot *s);

// =============================================================================
// O INSTANTANEO DE FACETAS
//
// Mesma fronteira, do outro lado: 97 lacos de higflow/src percorrem facetas, e
// MEDIDO, 95 deles leem apenas centro, tamanho e id.  Os outros 2 chamam o buffer
// de residuo, que precisa da CAIXA DA CELULA da faceta.
//
// POR ISSO ELE GUARDA A CAIXA DA CELULA, e nao o centro da faceta ja' pronto.
// Nao e' simetria com o instantaneo de celulas -- e' a medicao acima.  Com a
// caixa, `hig_get_facet_center` e `hig_get_facet_delta` sao reproduzidos pelas
// MESMAS contas, e o buffer de residuo tambem e' atendido.
//
// Uma faceta e' (celula, dim, dir): o centro e' o centro da celula com a
// coordenada `dim` trocada pelo canto de baixo (dir=0) ou de cima (dir=1), e o
// tamanho e' o tamanho da CELULA.  `dim` e `dir` sao guardados por faceta em vez
// de um `dim` unico para o dominio: um campo por dominio exigiria supor que todas
// as facetas dele compartilham a direcao, e a suposicao nao paga o byte que
// economiza.
// =============================================================================

typedef struct hig_facet_snapshot {
    int   n;            //!< facetas LOCAIS
    real *low;          //!< n * DIM, caixa da CELULA da faceta
    real *high;         //!< n * DIM
    signed char *dim;   //!< direcao da faceta, por faceta
    signed char *dir;   //!< 0 = face de baixo, 1 = face de cima
} hig_facet_snapshot;

hig_facet_snapshot *hfs_create(int n);
void hfs_destroy(hig_facet_snapshot *s);

//! \brief Preenche a partir de um sim_facet_domain (backend MTree).
//!
//! Indexa por `mp_lookup(sfd_get_domain_mapper(sfd), hig_get_fid(f))`, de modo
//! que o indice do arranjo E' o id local da faceta -- a mesma construcao do
//! instantaneo de celulas, e pela mesma razao.
hig_facet_snapshot *hfs_from_facet_domain(sim_facet_domain *sfd);

//! \brief O tamanho da faceta `i`.  E' o tamanho da CELULA dela, como manda o
//! `hig_get_facet_delta`.
static inline void hfs_delta(const hig_facet_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->high[i * DIM + d] - s->low[i * DIM + d];
}

//! \brief O centro da faceta `i`, pelas mesmas contas do `hig_get_facet_center`.
static inline void hfs_center(const hig_facet_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = (s->low[i * DIM + d] + s->high[i * DIM + d]) / 2.0;
    const int k = s->dim[i];
    p[k] = s->dir[i] ? s->high[i * DIM + k] : s->low[i * DIM + k];
}

//! \brief O canto inferior da CELULA da faceta `i`.  Para o buffer de residuo.
static inline void hfs_cell_low(const hig_facet_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->low[i * DIM + d];
}

//! \brief O canto superior da CELULA da faceta `i`.
static inline void hfs_cell_high(const hig_facet_snapshot *s, int i, Point p) {
    for (int d = 0; d < DIM; d++) p[d] = s->high[i * DIM + d];
}

//! \brief Compara dois instantaneos como CONJUNTOS de celulas.
//!
//! Dois backends podem numerar as celulas em ordens diferentes -- e vao, porque
//! a ordem e' consequencia de como cada um percorre --, entao comparar posicao a
//! posicao acusaria diferenca onde nao ha'.  Esta funcao casa por centro, dentro
//! de `tol`, e verifica que os tamanhos batem.
//!
//! Devolve 0 se forem o mesmo conjunto.  Caso contrario devolve o numero de
//! celulas de `a` sem par em `b`, e escreve em `detalhe` (se nao for NULL, com
//! ao menos 256 bytes) a primeira divergencia encontrada.
int hms_mesmo_conjunto(const hig_mesh_snapshot *a, const hig_mesh_snapshot *b,
                       real tol, char *detalhe);

#ifdef __cplusplus
}
#endif

#endif
