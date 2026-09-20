// O instantaneo pendurado no dominio: quando nasce, e o que acontece se a malha
// mudar depois.
//
// A decisao de projeto (doc/projeto-instantaneo.md, secao 3) e' que o instantaneo
// vive no `sim_domain` e e' produzido ANSIOSO, ao fim do `psd_synced_mapper`.  Isso
// so' se sustenta porque a malha e' IMUTAVEL depois da montagem -- medido: nao ha'
// `hig_refine_uniform` nem `hig_split` em higflow/src nem nos exemplos.
//
// "E' imutavel" nao e' garantia, e' observacao sobre o codigo de hoje.  O que a
// transforma em garantia sao os dois mecanismos abaixo, e este arquivo existe para
// afirmar os dois -- inclusive o que ABORTA, que sem teste seria so' uma intencao
// escrita em comentario.
//
//   nulo_antes_de_produzir          `sd_get_snapshot` devolve NULL enquanto
//       ninguem produziu.  Sem isto, "nao produzido" e "produzido vazio" seriam a
//       mesma coisa para quem le.
//
//   producao_serial_bate_com_a_arvore  o caminho SEM MPI produz, e produz certo.
//       E' a clausula P3 aplicada ao instantaneo: interface que so' funciona
//       acoplada ao particionamento e' a mais dificil de substituir, e substituir
//       e' o objetivo.
//
//   detector_acusa_refino_no_lugar  refinar uma celula que ja' esta' no dominio
//       NAO passa pela API do dominio, entao nenhuma guarda o intercepta.  O
//       `sd_snapshot_is_current` e' o detector para esse caso -- O(n), para teste
//       e depuracao, nao para laco quente.
//
//   guarda_aborta_se_a_malha_mudar_depois  acrescentar arvore depois de produzir
//       PARA o programa.  Verificado em processo filho: o pai confere que o filho
//       morreu de SIGABRT.  Afirmar isso de dentro do proprio processo seria
//       impossivel -- e' justamente o processo que deixaria de existir.
//
// POR QUE A GUARDA NAO E' `assert`.  A higtree compila com -DNDEBUG no modo
// otimizado, que e' o modo que se roda; `assert` sumiria exatamente ali.  Um guarda
// que desaparece em release nao e' guarda.  Por isso ela e' `fprintf` + `abort`, e
// por isso o caso acima mede o sinal e nao a mensagem.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "hig-mesh-snapshot.h"
#include "utils.h"
#include "testing.h"

#define NC 4

// Dominio montado EM SERIE: sem lb, sem pg, sem psd.  O mapeador e' atribuido a
// mao, que e' o que a montagem serial faz.
static sim_domain *monta_serial(hig_cell **raiz_out) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    hig_refine_uniform(raiz, nc);

    sim_domain *sd = sd_create(NULL);
    sd_add_higtree(sd, raiz);
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *it = sd_get_domain_celliterator(sd);
    mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);

    if (raiz_out != NULL) *raiz_out = raiz;
    return sd;
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    // ------------------------------------------------------------------
    t_case("nulo_antes_de_produzir");
    {
        sim_domain *sd = monta_serial(NULL);
        T_CHECK_MSG(sd_get_snapshot(sd) == NULL,
            "o dominio recem-montado ja' trouxe instantaneo -- alguem o produziu "
            "em lugar que nao o `sd_compute_snapshot`");
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("producao_serial_bate_com_a_arvore");
    {
        sim_domain *sd = monta_serial(NULL);
        sd_compute_snapshot(sd);
        const hig_mesh_snapshot *s = sd_get_snapshot(sd);

        T_CHECK_MSG(s != NULL, "sd_compute_snapshot nao pendurou o instantaneo");

        if (s != NULL) {
            mp_mapper *m = sd_get_domain_mapper(sd);
            int visitadas = 0, divergentes = 0;
            higcit_celliterator *it;
            for (it = sd_get_domain_celliterator(sd); !higcit_isfinished(it);
                 higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                const int i = mp_lookup(m, hig_get_cid(c));
                visitadas++;
                if (i < 0 || i >= s->n) { divergentes++; continue; }
                Point ce, de;
                hig_get_center(c, ce);
                hig_get_delta(c, de);
                for (int d = 0; d < DIM; d++) {
                    if (fabs(s->center[i * DIM + d] - ce[d]) > 1e-15 ||
                        fabs(s->delta[i * DIM + d]  - de[d]) > 1e-15) {
                        divergentes++;
                        break;
                    }
                }
            }
            higcit_destroy(it);
            T_CHECK_MSG(visitadas == s->n,
                "o iterador visitou %d celulas e o instantaneo tem %d",
                visitadas, s->n);
            T_CHECK_MSG(divergentes == 0,
                "%d celula(s) divergem entre o instantaneo serial e a arvore",
                divergentes);
            T_CHECK_MSG(sd_snapshot_is_current(sd),
                "o detector diz que o instantaneo recem-produzido ja' esta' velho");
        }
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("detector_acusa_refino_no_lugar");
    {
        hig_cell *raiz = NULL;
        sim_domain *sd = monta_serial(&raiz);
        sd_compute_snapshot(sd);
        T_CHECK_MSG(sd_snapshot_is_current(sd),
            "o detector ja' reprovava antes de a malha mudar -- ele nao mede o "
            "que este caso precisa medir");

        // Refino EM CIMA de uma celula que ja' esta' no dominio.  Nao passa por
        // `sd_add_higtree`, entao a guarda nao ve'; quem tem de ver e' o detector.
        Point p;
        POINT_ASSIGN_SCALAR(p, 0.375);
        hig_cell *c = hig_get_cell_with_point(raiz, p);
        T_CHECK_MSG(c != NULL, "nao achei a celula para refinar");
        if (c != NULL) {
            int n2[DIM];
            for (int d = 0; d < DIM; d++) n2[d] = 2;
            hig_refine_uniform(c, n2);
            T_CHECK_MSG(!sd_snapshot_is_current(sd),
                "a malha ganhou celulas e o detector continua dizendo que o "
                "instantaneo esta' em dia");
        }
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("guarda_aborta_se_a_malha_mudar_depois");
    {
        fflush(NULL);
        pid_t filho = fork();
        if (filho == 0) {
            // O filho deve morrer aqui.  A mensagem da guarda vai para /dev/null:
            // o que se afirma e' o sinal, e a saida suja confundiria quem le o log.
            freopen("/dev/null", "w", stderr);
            sim_domain *sd = monta_serial(NULL);
            sd_compute_snapshot(sd);

            Point lo, hi;
            POINT_ASSIGN_SCALAR(lo, 1.0);
            POINT_ASSIGN_SCALAR(hi, 2.0);
            hig_cell *outra = hig_create_root(lo, hi);
            int nc[DIM];
            for (int d = 0; d < DIM; d++) nc[d] = NC;
            hig_refine_uniform(outra, nc);
            sd_add_higtree(sd, outra);      // <-- tem de abortar

            _exit(0);                       // se chegou aqui, a guarda nao existe
        }
        T_CHECK_MSG(filho > 0, "fork falhou");
        if (filho > 0) {
            int st = 0;
            waitpid(filho, &st, 0);
            T_CHECK_MSG(WIFSIGNALED(st) && WTERMSIG(st) == SIGABRT,
                "acrescentar arvore depois de produzir o instantaneo deveria "
                "abortar; o filho terminou %s (codigo %d, sinal %d)",
                WIFEXITED(st) ? "normalmente" : "por outro sinal",
                WIFEXITED(st) ? WEXITSTATUS(st) : -1,
                WIFSIGNALED(st) ? WTERMSIG(st) : -1);
        }
    }

    return t_end();
}
