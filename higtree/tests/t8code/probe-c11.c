// SONDA: o t8code cumpre a clausula C11 do contrato de Mesh?
//
// C11 diz: "O estencil ATRAVESSA salto de nivel maior que 2:1."  A contribuicao de
// Sousa et al. (2019) e' minimos quadrados moveis em arvore NAO GRADUADA, entao isso
// e' requisito da interface, nao detalhe de implementacao.  E' a clausula de maior
// risco para o t8code, cujo balanceamento e' opcional mas cujas rotinas de vizinhanca
// e de ghost foram construidas em torno do caso balanceado.
//
// A pergunta se parte em duas, e as duas sao medidas aqui:
//
//   Q1  o t8code REPRESENTA uma floresta com salto 4:1, sem forcar balanceamento?
//   Q2  ele REPORTA os vizinhos de face atraves desse salto?
//
// Se Q1 ou Q2 falhar, a premissa de que t8code e MTree podem ser pares
// intercambiaveis precisa ser revista ANTES da integracao.
//
// MEDIDO em 19/09/2026, t8code v4.0.0-26.08, 2-D, np=1:
//
//     Q1  niveis presentes: min=2 max=4  -> razao de tamanho 4
//         REPRESENTA salto 4:1 sem balanceamento
//     Q2  elemento nivel 2, face 0: 4 vizinhos de nivel 4
//         REPORTA vizinhos atraves do salto
//
// Os 4 vizinhos sao o esperado: a face de um elemento de nivel 2, em 2-D, e'
// subdividida 2^2 = 4 vezes por dois niveis de refino.  A PREMISSA DO PROJETO
// SOBREVIVE AO SEU TESTE DE MAIOR RISCO.
//
// ISTO NAO E' O ADAPTADOR.  Responde a pergunta de viabilidade -- a unica capaz
// de invalidar a premissa -- contra a API do t8code diretamente.  O adaptador
// tem de PRODUZIR as estruturas que as consultas da HiGTree leem, e e' ai' que
// mora o resto do trabalho.
//
// TRES RESTRICOES DE INTEGRACAO DESCOBERTAS AO COMPILAR ISTO:
//
//   1. o t8code exige C++20 (`cxx_std_20` no alvo exportado).  A HiGTree
//      constroi com `-std=gnu++17`.  Ou a unidade do adaptador compila em C++20
//      separada, ou a arvore inteira sobe de padrao.
//   2. o consumidor TEM de definir -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1.  Sem
//      isso o t8.h para com "MPI configured differently in t8code and libsc",
//      que aponta para o libsc e nao para a definicao que falta.
//   3. `t8_cmesh_init(&cmesh)` e' obrigatorio ANTES de `t8_cmesh_new_hypercube`.
//      Sem ele o segfault sai dentro de t8_cmesh_set_tree_class, tres quadros
//      abaixo, sem mencionar a inicializacao.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <stdio.h>
#include <stdlib.h>

// Refina, sem balancear, apenas o elemento que contem o canto (0,0,...): duas vezes
// seguidas, o que produz salto de QUATRO para um contra a vizinha imediata.
static int adapt_canto (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                        const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                        const t8_scheme_c *scheme, const int is_family,
                        const int num_elements, t8_element_t *elements[])
{
    double centro[3];
    t8_forest_element_centroid (forest_from, which_tree, elements[0], centro);
    // canto inferior: refina so' quem esta' no primeiro oitavo em todas as direcoes
    if (centro[0] < 0.25 && centro[1] < 0.25) return 1;
    return 0;
}

int main (int argc, char **argv)
{
    int mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
    t8_init (SC_LP_ERROR);

    const int NIVEL = 2;            // 4 por direcao
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);       // obrigatorio ANTES do gerador
    t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
    const t8_scheme_c *scheme = t8_scheme_new_default ();

    t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, NIVEL, 0, sc_MPI_COMM_WORLD);

    // Dois passos de refino SEM balancear -> 4:1
    for (int passo = 0; passo < 2; passo++) {
        t8_forest_t novo;
        t8_forest_init (&novo);
        t8_forest_set_adapt (novo, f, adapt_canto, 0);
        t8_forest_set_ghost (novo, 1, T8_GHOST_FACES);
        t8_forest_commit (novo);
        f = novo;
    }

    const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
    int nivel_min = 99, nivel_max = -1;
    for (t8_locidx_t it = 0; it < ntrees; it++) {
        const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
        for (t8_locidx_t ie = 0; ie < ne; ie++) {
            const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
            const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
            const int lv = scheme->element_get_level (ec, e);
            if (lv < nivel_min) nivel_min = lv;
            if (lv > nivel_max) nivel_max = lv;
        }
    }
    printf ("  Q1  niveis presentes: min=%d max=%d  -> razao de tamanho 2^%d = %d\n",
            nivel_min, nivel_max, nivel_max - nivel_min, 1 << (nivel_max - nivel_min));
    printf ("  Q1  %s\n", (nivel_max - nivel_min >= 2)
            ? "REPRESENTA salto 4:1 sem balanceamento" : "NAO chegou a 4:1");

    // Q2: procura um elemento GROSSO que tenha vizinho FINO 4x menor, e conta
    // quantos vizinhos o t8code reporta atraves daquela face.
    int achou = 0, viz_max = 0;
    for (t8_locidx_t it = 0; it < ntrees && !achou; it++) {
        const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
        const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
        for (t8_locidx_t ie = 0; ie < ne && !achou; ie++) {
            const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
            const int lv = scheme->element_get_level (ec, e);
            if (lv != nivel_min) continue;
            const int nfaces = scheme->element_get_num_faces (ec, e);
            for (int fa = 0; fa < nfaces; fa++) {
                t8_element_t **vz; int *dual; int nviz = 0;
                t8_locidx_t *idx; t8_eclass_t vec;
                t8_forest_leaf_face_neighbors (f, it, e, (const t8_element_t ***) &vz,
                                               fa, &dual, &nviz, &idx, &vec);
                if (nviz > 0) {
                    const int lvv = scheme->element_get_level (vec, vz[0]);
                    if (lvv - lv >= 2) {
                        printf ("  Q2  elemento nivel %d, face %d: %d vizinho(s) de "
                                "nivel %d  -> salto de %dx atravessado\n",
                                lv, fa, nviz, lvv, 1 << (lvv - lv));
                        achou = 1; viz_max = nviz;
                    }
                    T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
                    if (achou) break;
                } else if (nviz == 0) {
                    T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
                }
            }
        }
    }
    printf ("  Q2  %s\n", achou
            ? "REPORTA vizinhos atraves do salto 4:1"
            : "NAO achou vizinhanca atraves de salto >= 4:1");
    (void) viz_max;

    t8_forest_unref (&f);
    sc_finalize ();
    return sc_MPI_Finalize ();
}
