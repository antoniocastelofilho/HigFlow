// O instantaneo sob particionamento real: np = 1, 2 e 3.
//
// POR QUE ESTE TESTE EXISTE, e por que o test-mesh-snapshot nao bastava.  O
// `hms_from_domain` afirma que O INDICE E' A IDENTIDADE -- `center[i*DIM+d]` e' a
// celula de id local i -- e e' essa afirmacao que apaga `hig_get_cid` dos lacos
// quentes de higflow/src.  Ela depende de o mapeador dar aos locais exatamente
// [0, n) e a' franja os ids SEGUINTES.
//
// Ate' aqui isso so' fora exercitado EM SERIE, e com a franja montada a mao: o
// test-mesh-snapshot roda em np=1 e chama `sd_add_fringe_higtree` ele proprio,
// numerando a franja depois dos locais porque ele escolhe assim.  Ou seja, o teste
// confirmava a convencao que ele mesmo tinha imposto.  Quem impoe a convencao de
// verdade e' o `psd_synced_mapper`, que so' roda sob MPI -- e e' quem o solver usa.
//
// A DIFERENCA IMPORTA porque o modo de falhar e' silencioso.  Se a franja recebesse
// id dentro de [0, n), o `hms_from_domain` nao acusaria nada: ele sobrescreveria a
// linha de uma celula local com a geometria de uma celula de franja, e o
// instantaneo sairia com o tamanho certo e uma celula errada dentro.  O laco
// migrado leria a celula errada sem erro nenhum.
//
// O QUE E' AFIRMADO, tudo por REDUCAO GLOBAL e nao por rank -- um teste paralelo em
// que cada rank conclui esconde a falha de um rank no meio das linhas dos outros:
//
//   instantaneo_existe_em_todo_rank      `hms_from_domain` nao devolve NULL em
//       rank nenhum.  Ele devolve NULL quando o id cai fora de [0, n), entao este
//       caso ja' e' a primeira forma da afirmacao -- reduzida por MIN.
//
//   tamanho_e_o_numero_de_celulas_locais o instantaneo tem exatamente as celulas
//       que o iterador de dominio ve' (C6), em cada rank, e a soma global delas e'
//       o dominio inteiro.  Pega franja que entrou E celula local que ficou de
//       fora, que contagem por rank nao distingue.
//
//   indice_e_o_id_local_com_franja       para toda celula local, o id que o
//       mapeador devolve indexa a PROPRIA celula.  E' o test-mesh-snapshot, agora
//       com franja de verdade em volta.
//
//   franja_nao_ocupa_indice_local        para toda celula de franja, o id NAO cai
//       em [0, n).  E' o caso que discrimina: os tres acima passam mesmo que uma
//       celula de franja tenha sobrescrito uma local, porque o que eles conferem e'
//       o conjunto que o ITERADOR ve', e ele nao ve' a franja.
//
// POR QUE A IDENTIDADE DE INDICE E' VERDADEIRA POR CONSTRUCAO, e por que isso
// NAO torna o teste dispensavel.  O `_psd_setmapper` numera os locais com
// `mp_assign_from_celliterator(m, sd_get_domain_celliterator(sd), 0)` -- o MESMO
// iterador que o `hms_from_domain` percorre.  O id que o mapeador devolve e', por
// definicao, o contador do percurso.  Nao ha' como um divergir do outro enquanto
// os dois sairem dali.
//
// Consequencia, MEDIDA e nao suposta: trocar `mp_lookup(m, hig_get_cid(c))` por
// `n_visitadas++` dentro do `hms_from_domain` deixa a suite inteira VERDE, nos dois
// DIM e nos tres np.  Nao e' falha de oraculo -- e' que os dois sao a mesma funcao.
// Qualquer teste que tente distinguir um do outro esta' medindo nada.
//
// O que este teste afirma, entao, e' a CONVENCAO DE NUMERACAO da producao, e nao a
// aritmetica do instantaneo.  E ai' ele discrimina.  VALIDADO NO SENTIDO INVERSO
// fazendo o `_psd_setmapper` numerar a FRANJA PRIMEIRO:
//
//     test-mesh-snapshot (serial)   3 de 3 passam -- ele numera a franja a mao e
//                                   nunca toca no `_psd_setmapper`
//     test-mesh-snapshot-parallel   np=1  passa   (nao ha' franja)
//                                   np=2  REPROVA em 2D e em 3D
//                                   np=3  REPROVA em 2D e em 3D
//
// E' essa a lacuna que este arquivo fecha: a convencao de que os locais ficam em
// [0, n) e a franja vem depois estava sendo IMITADA pelo teste serial, nao
// exercitada.  Quem a impoe e' o `_psd_setmapper`, e ate' aqui nada falhava se ela
// mudasse.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "hig-mesh-snapshot.h"
#include "utils.h"
#include "testing.h"

#define NC       8        // celulas por direcao na raiz
#define FRINGE   2

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    // Mesmo caminho de PRODUCAO do test-fringe-parallel: lb_create,
    // lb_add_input_tree, lb_calc_partition, psd_create, psd_synced_mapper.  E' o
    // caminho que o solver percorre, e o unico em que a franja e a numeracao dela
    // sao produzidas por quem as produz de verdade.
    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, FRINGE);

    load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
    if (rank == 0) {
        Point lo, hi;
        POINT_ASSIGN_SCALAR(lo, 0.0);
        POINT_ASSIGN_SCALAR(hi, 1.0);
        hig_cell *raiz = hig_create_root(lo, hi);
        int nc[DIM];
        for (int d = 0; d < DIM; d++) nc[d] = NC;
        hig_refine_uniform(raiz, nc);
        lb_add_input_tree(lb, raiz, true, 0);
    }
    lb_calc_partition(lb, pg);

    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);
    for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i) {
        sd_add_higtree(sd, lb_get_local_tree(lb, i, NULL));
    }
    lb_destroy(lb);

    psim_domain *psd = psd_create(sd, pg);
    psd_synced_mapper(psd);
    mp_mapper *m = sd_get_domain_mapper(sd);

    // ----------------------------------------------- o que este rank enxerga
    hig_mesh_snapshot *s = hms_from_domain(sd);

    long n_local = 0;
    higcit_celliterator *cit;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        n_local++;
    }
    higcit_destroy(cit);

    // Divergencias entre o instantaneo e a arvore, celula local a celula local.
    long fora_de_lugar = 0;
    if (s != NULL) {
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            hig_cell *c = higcit_getcell(cit);
            const int i = mp_lookup(m, hig_get_cid(c));
            if (i < 0 || i >= s->n) { fora_de_lugar++; continue; }
            Point ce, de;
            hig_get_center(c, ce);
            hig_get_delta(c, de);
            for (int d = 0; d < DIM; d++) {
                Point sa_c, sa_d;
                hms_center(s, i, sa_c);
                hms_delta(s, i, sa_d);
                if (fabs(sa_c[d] - ce[d]) > 1e-15 ||
                    fabs(sa_d[d] - de[d]) > 1e-15) {
                    fora_de_lugar++;
                    break;
                }
            }
        }
        higcit_destroy(cit);
    }

    // A franja invadiu a faixa dos locais?  Percorre as arvores de franja e
    // pergunta o id de cada folha.  `psd_synced_mapper` deve te-las numerado a
    // partir de n_local; qualquer id em [0, n) e' uma linha local sobrescrita.
    long franja_invasora = 0, franja_vista = 0;
    if (s != NULL) {
        for (unsigned k = 0; k < sd_get_num_fringe_higtrees(sd); ++k) {
            cit = higcit_create_all_leaves(sd_get_fringe_higtree(sd, k));
            for (; !higcit_isfinished(cit); higcit_nextcell(cit)) {
                const int i = mp_lookup(m, hig_get_cid(higcit_getcell(cit)));
                franja_vista++;
                if (i >= 0 && i < s->n) franja_invasora++;
            }
            higcit_destroy(cit);
        }
    }

    const long tem_instantaneo = (s != NULL) ? 1 : 0;
    const long n_snapshot      = (s != NULL) ? (long) s->n : -1;
    const long tamanho_bate    = (n_snapshot == n_local) ? 1 : 0;

    // ------------------------------------------------------------ reducoes
    long instantaneo_min, tamanho_min, fora_max, invasora_max, franja_max;
    long n_global;
    MPI_Allreduce(&tem_instantaneo, &instantaneo_min, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&tamanho_bate,    &tamanho_min,     1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&fora_de_lugar,   &fora_max,        1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&franja_invasora, &invasora_max,    1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&franja_vista,    &franja_max,      1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&n_snapshot,      &n_global,        1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    long esperado = 1;
    for (int d = 0; d < DIM; d++) esperado *= NC;

    // So' o rank 0 conclui.  Ninguem chama MPI_Finalize: o higtree_initialize
    // registra o PetscFinalize no atexit, e ele ja' o chama.
    if (rank != 0) { if (s != NULL) hms_destroy(s); return 0; }

    t_case("instantaneo_existe_em_todo_rank");
    T_CHECK_MSG(instantaneo_min == 1,
        "hms_from_domain devolveu NULL em ao menos um rank (np=%d) -- algum id "
        "local caiu fora de [0, n)", ntasks);

    t_case("tamanho_e_o_numero_de_celulas_locais");
    T_CHECK_MSG(tamanho_min == 1,
        "em ao menos um rank o instantaneo nao tem o mesmo numero de celulas que "
        "o iterador de dominio (np=%d): a franja entrou, ou celula local ficou de "
        "fora", ntasks);
    T_CHECK_MSG(n_global == esperado,
        "soma global do tamanho dos instantaneos = %ld, esperado %ld (np=%d)",
        n_global, esperado, ntasks);

    t_case("indice_e_o_id_local_com_franja");
    T_CHECK_MSG(fora_max == 0,
        "ate' %ld celula(s) locais num mesmo rank cujo id nao indexa a propria "
        "posicao no instantaneo (np=%d)", fora_max, ntasks);

    t_case("franja_nao_ocupa_indice_local");
    if (ntasks == 1) {
        T_CHECK_MSG(franja_max == 0,
            "com np=1 nao deveria haver franja, e ha' ate' %ld celula(s)", franja_max);
    } else {
        T_CHECK_MSG(franja_max > 0,
            "com np=%d todo rank deveria ver franja, e o maximo visto foi %ld -- "
            "sem franja este caso nao afirma nada", ntasks, franja_max);
        T_CHECK_MSG(invasora_max == 0,
            "ate' %ld celula(s) de franja num mesmo rank receberam id em [0, n) e "
            "sobrescreveram linha de celula local (np=%d)", invasora_max, ntasks);
    }

    if (s != NULL) hms_destroy(s);
    return t_end();
}
