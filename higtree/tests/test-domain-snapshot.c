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

// Dominio com coordenadas NAO DIADICAS.  Em [0,1] dividido por potencia de dois,
// todo centro e' exato em binario e `(a+b)/2` coincide com `a/2 + b/2` -- medido:
// com a malha diadica, sabotar a formula deixava o caso VERDE.  Comecando em 0,1
// com passo 0,2 nenhuma coordenada e' representavel, e as duas formas divergem.
static sim_domain *monta_nao_diadico(void) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.1);
    POINT_ASSIGN_SCALAR(hi, 0.7);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = 3;
    hig_refine_uniform(raiz, nc);

    sim_domain *sd = sd_create(NULL);
    sd_add_higtree(sd, raiz);
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *it = sd_get_domain_celliterator(sd);
    mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);
    return sd;
}

// Dominio de FACETAS montado em serie, no mesmo espirito do de celulas.
static sim_facet_domain *monta_facetas(int dim, hig_cell **raiz_out) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.1);
    POINT_ASSIGN_SCALAR(hi, 0.7);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = 3;
    hig_refine_uniform(raiz, nc);

    sim_facet_domain *sfd = sfd_create(NULL, dim);
    sfd_add_higtree(sfd, raiz);
    sfd_adjust_facet_ids(sfd);
    sfd_compute_sfbi(sfd);
    mp_mapper *mf = sfd_get_domain_mapper(sfd);
    higfit_facetiterator *fit = sfd_get_domain_facetiterator(sfd);
    mp_assign_from_facetiterator(mf, fit, 0);
    higfit_destroy(fit);
    if (raiz_out) *raiz_out = raiz;
    return sfd;
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
                    Point sa_c, sa_d;
                    hms_center(s, i, sa_c);
                    hms_delta(s, i, sa_d);
                    if (fabs(sa_c[d] - ce[d]) > 1e-15 ||
                        fabs(sa_d[d] - de[d]) > 1e-15) {
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

    // ------------------------------------------------------------------
    t_case("derivacao_e_bit_a_bit_igual_a_arvore");
    {
        // A JUSTIFICATIVA DE GUARDAR A CAIXA depende disto, e nada menos.  Os
        // lacos migrados alimentam `compute_value_at_point` e
        // `compute_facet_value_at_point` com centro, delta e cantos; se o
        // instantaneo devolvesse valores que diferem no ultimo bit, a saida VTK
        // mudaria e a referencia so' nao acusaria por causa da tolerancia.
        //
        // Por isso a comparacao aqui e' `!=` e nao tolerancia: o instantaneo
        // COPIA `c->lowpoint`, e centro e delta saem das MESMAS contas do
        // `hig_get_center` e do `hig_get_delta`.  Nao ha' arredondamento a
        // tolerar -- ou e' identico, ou a premissa caiu.
        //
        // O QUE ESTE CASO NAO PRENDE, e vale escrito para ninguem confiar demais
        // nele: a escolha entre formas algebricamente equivalentes do centro.
        // MEDIDO -- `(lo+hi)/2`, `lo/2+hi/2` e `lo+(hi-lo)/2` dao o MESMO bit em
        // 200 mil caixas aleatorias.  Dividir por dois e' exato em binario, entao
        // nao ha' o que discriminar ali.  Quem sustenta a decisao de guardar a
        // caixa e' o caso seguinte.
        sim_domain *sd = monta_nao_diadico();
        sd_compute_snapshot(sd);
        const hig_mesh_snapshot *s = sd_get_snapshot(sd);
        mp_mapper *m = sd_get_domain_mapper(sd);

        int difs = 0;
        char primeira[256]; primeira[0] = '\0';
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sd); !higcit_isfinished(it);
             higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            const int i = mp_lookup(m, hig_get_cid(c));
            Point a_lo, a_hi, a_ce, a_de, t_lo, t_hi, t_ce, t_de;
            hms_low(s, i, a_lo);     hig_get_lowpoint(c, t_lo);
            hms_high(s, i, a_hi);    hig_get_highpoint(c, t_hi);
            hms_center(s, i, a_ce);  hig_get_center(c, t_ce);
            hms_delta(s, i, a_de);   hig_get_delta(c, t_de);
            for (int d = 0; d < DIM; d++) {
                if (a_lo[d] != t_lo[d] || a_hi[d] != t_hi[d] ||
                    a_ce[d] != t_ce[d] || a_de[d] != t_de[d]) {
                    if (!difs) {
                        snprintf(primeira, sizeof primeira,
                            "celula %d, direcao %d: centro %.17g contra %.17g, "
                            "delta %.17g contra %.17g", i, d,
                            (double) a_ce[d], (double) t_ce[d],
                            (double) a_de[d], (double) t_de[d]);
                    }
                    difs++;
                    break;
                }
            }
        }
        higcit_destroy(it);
        T_CHECK_MSG(difs == 0,
            "%d celula(s) em que o instantaneo nao devolve exatamente o que a "
            "arvore devolve.  %s", difs, primeira);
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("reconstruir_o_canto_nao_seria_exato");
    {
        // ESTE CASO GUARDA A DECISAO DE PROJETO, e e' o unico que a guarda.
        //
        // Se alguem "simplificar" o instantaneo de volta para centro e delta, os
        // lacos do `hig-flow-io.c` que hoje leem `hms_low`/`hms_high` passariam a
        // reconstruir o canto por `centro - delta/2`.  Isso e' exato em algebra e
        // NAO em ponto flutuante: medido, difere em ~0,6% de caixas aleatorias e
        // em 1 das 3 celulas por direcao desta malha.  A diferenca entraria em
        // `compute_facet_value_at_point` e sairia no VTK.
        //
        // Entao o caso AFIRMA a perda: se a reconstrucao passasse a ser exata
        // para toda celula desta malha, ele falha -- nao porque algo quebrou, mas
        // porque a razao de guardar a caixa deixou de valer aqui e o teste
        // precisa ser refeito em malha que a exponha.
        sim_domain *sd = monta_nao_diadico();
        sd_compute_snapshot(sd);
        const hig_mesh_snapshot *s = sd_get_snapshot(sd);

        int perdas = 0;
        for (int i = 0; i < s->n; i++) {
            Point lo, ce, de;
            hms_low(s, i, lo);
            hms_center(s, i, ce);
            hms_delta(s, i, de);
            for (int d = 0; d < DIM; d++) {
                if (ce[d] - de[d] / 2.0 != lo[d]) { perdas++; break; }
            }
        }
        T_CHECK_MSG(perdas > 0,
            "em nenhuma das %d celulas a reconstrucao `centro - delta/2` perdeu "
            "bit -- esta malha deixou de sustentar a decisao de guardar a caixa",
            s->n);
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("facetas_derivacao_e_bit_a_bit_igual_a_arvore");
    {
        // O instantaneo de facetas guarda a CAIXA DA CELULA, e centro e tamanho
        // saem dela pelas mesmas contas do `hig_get_facet_center` e do
        // `hig_get_facet_delta`.  Como no caso das celulas, a comparacao e' `!=`:
        // nao ha' arredondamento a tolerar.
        int difs_tot = 0;
        char primeira[256]; primeira[0] = '\0';
        for (int dim = 0; dim < DIM; dim++) {
            sim_facet_domain *sfd = monta_facetas(dim, NULL);
            sfd_compute_snapshot(sfd);
            const hig_facet_snapshot *s = sfd_get_snapshot(sfd);
            mp_mapper *mf = sfd_get_domain_mapper(sfd);

            higfit_facetiterator *fit;
            for (fit = sfd_get_domain_facetiterator(sfd); !higfit_isfinished(fit);
                 higfit_nextfacet(fit)) {
                hig_facet *f = higfit_getfacet(fit);
                const int i = mp_lookup(mf, hig_get_fid(f));
                Point a_ce, a_de, t_ce, t_de;
                hfs_center(s, i, a_ce);
                hfs_delta(s, i, a_de);
                hig_get_facet_center(f, t_ce);
                hig_get_facet_delta(f, t_de);
                for (int d = 0; d < DIM; d++) {
                    if (a_ce[d] != t_ce[d] || a_de[d] != t_de[d]) {
                        if (!difs_tot) {
                            snprintf(primeira, sizeof primeira,
                                "dim %d, faceta %d, direcao %d: centro %.17g "
                                "contra %.17g", dim, i, d,
                                (double) a_ce[d], (double) t_ce[d]);
                        }
                        difs_tot++;
                        break;
                    }
                }
            }
            higfit_destroy(fit);
            sfd_destroy(sfd);
        }
        T_CHECK_MSG(difs_tot == 0,
            "%d faceta(s) em que o instantaneo nao devolve exatamente o que a "
            "arvore devolve.  %s", difs_tot, primeira);
    }

    // ------------------------------------------------------------------
    t_case("facetas_oraculo_aprova_e_acusa");
    {
        sim_facet_domain *sfd = monta_facetas(0, NULL);
        sfd_compute_snapshot(sfd);
        char detalhe[256];

        const int ok = sfd_snapshot_verify(sfd, detalhe, sizeof detalhe);
        T_CHECK_MSG(ok == 0,
            "o oraculo de facetas reprovou um instantaneo recem-produzido: "
            "%d linha(s).  %s", ok, detalhe);

        // Corrompe de proposito: sem isto, "aprova" nao diz nada.
        hig_facet_snapshot *s = (hig_facet_snapshot *) sfd_get_snapshot(sfd);
        if (s != NULL && s->n >= 2) {
            for (int d = 0; d < DIM; d++) {
                real t = s->low[d];
                s->low[d] = s->low[(s->n - 1) * DIM + d];
                s->low[(s->n - 1) * DIM + d] = t;
                t = s->high[d];
                s->high[d] = s->high[(s->n - 1) * DIM + d];
                s->high[(s->n - 1) * DIM + d] = t;
            }
            const int ruins = sfd_snapshot_verify(sfd, detalhe, sizeof detalhe);
            T_CHECK_MSG(ruins >= 2,
                "troquei duas linhas e o oraculo de facetas acusou %d "
                "divergencia(s) -- deveria acusar as duas", ruins);
        }
        sfd_destroy(sfd);
    }

    return t_end();
}
