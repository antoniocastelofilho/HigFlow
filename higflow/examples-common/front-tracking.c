// Adaptador entre o solver e o front-tracking -- o analogo de
// examples-common/fronteira-imersa.c, para a INTERFACE (tensao superficial), nao
// para o corpo rigido.  Compilado POR EXEMPLO (2D), fora da lib: o nucleo
// hig-flow-front-tracking.{c,h} nao toca o solver, e este arquivo faz a ponte.
//
// FASE B2 (gota estatica, lei de Laplace): a forca de tensao superficial
// F = sigma*kappa*n dos marcadores e' espalhada em dpFU (o campo de fonte que a
// equacao de momento le'), ANTES do preditor -- mesmo gancho que o corpo rigido
// usa.  Reusa a busca de faceta escalonada do corpo rigido (fi_suporte_facetas),
// para nao duplicar a parte delicada.
//
// SERIAL POR ENQUANTO: o reduce de franja (dono<-franja) NAO esta' aqui.  Em np=1
// nao ha' franja e o espalhamento e' completo; o paralelo e' trabalho da fase de
// escala (o mesmo padrao do baseline VOF, que ficou serial primeiro).

#include "hig-flow-kernel.h"
#include "hig-flow-fronteira-imersa.h"   // fi_suporte_facetas, fi_suporte_capacidade
#include "hig-flow-front-tracking.h"
#include "hig-mesh-snapshot.h"           // hig_facet_snapshot, hfs_center (para o oraculo)
#include "hig-flow-eval.h"               // compute_center_p_left/right
#include "hig-flow-discret.h"            // compute_dpdx_at_point
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    ft_frente *frente;
    real       sigma;
    int        advecta;      // 0 = frente fixa (Laplace puro); 1 = advecta
    int        balanceado;   // 1 = forca como sigma*kappa*grad(H) balanceado
    int        passo;        // contador proprio, para o diagnostico periodico
    Point      centro_ini;   // centroide inicial, para medir deriva
    real       area_ini;     // area inicial, para medir conservacao
} ft_ctx;

// Tamanho de celula na posicao X (h do nucleo ali).  Malha uniforme: delta[0].
static real _h_em(sim_facet_domain *sfd, const Point X)
{
    hig_cell *c = sfd_get_cell_with_point(sfd, (real *) X);
    if (c == NULL) return 0.0;
    Point d; hig_get_delta(c, d);
    return d[0];
}

// Espalha F = sigma*kappa*n na malha escalonada (dpFU), ACUMULANDO.  Serial.
//
//   contribuicao a' faceta = F_dim * ds * (produto phi) / h^DIM
//
// Conservacao (o mesmo oraculo do standalone): SUM_faceta valor*h^DIM = SUM_k
// F_k*ds_k, porque o produto phi soma 1 no suporte.  Para curva fechada isso
// tende a zero (Laplace: forca liquida nula).
// Regra de posse do espalhamento, por ambiente.  Existe para MEDIR: as tres
// variantes rodam no mesmo binario e o oraculo de conservacao global as separa.
static int _posse_regra(void)
{
    static int r = -1;
    if (r < 0) { const char *s = getenv("FT_POSSE"); r = (s != NULL) ? atoi(s) : 0; }
    return r;
}

static void ft_espalha_tensao_solver(ft_frente *frente, real sigma,
                                     sim_facet_domain *sfd[DIM],
                                     distributed_property *dpF[DIM])
{
    int n = ft_num(frente);
    Point *pos = (Point *) malloc((size_t) n * sizeof(Point));
    Point *F   = (Point *) malloc((size_t) n * sizeof(Point));
    real  *ds  = (real  *) malloc((size_t) n * sizeof(real));
    ft_forcas_tensao(frente, sigma, pos, F, ds);

    int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    for (int k = 0; k < n; k++) {
        real h = _h_em(sfd[0], pos[k]);
        if (h <= 0.0) continue;               // marcador fora do dominio local
        real hd = 1.0;
        for (int d = 0; d < DIM; d++) hd *= h;

        for (int dim = 0; dim < DIM; dim++) {
            int m = fi_suporte_facetas(sfd[dim], dim, pos[k], h, lids, pesos, capac);
            real esc = F[k][dim] * ds[k] / hd;
            // REGRA DE POSSE (FT_POSSE), para medir qual da' conservacao exata:
            //   0 = acumula em tudo que encontra (inclusive franja) + dp_sync
            //       -- o estado atual: a contribuicao em franja se perde
            //   1 = acumula em tudo + reducao franja->dono
            //       -- MEDIDO pior, porque marcador da costura conta duas vezes
            //   2 = acumula SO' em faceta PROPRIA + dp_sync
            //       -- a frente e' replicada, entao o dono de uma faceta tem
            //          todos os marcadores e calcula a contribuicao dela sozinho,
            //          sem comunicacao nenhuma
            for (int i = 0; i < m; i++) {
                if (_posse_regra() == 2 && lids[i] >= fi_dp_num_proprias(dpF[dim]))
                    continue;                      // faceta de franja: nao e' minha
                dp_add_value(dpF[dim], lids[i], esc * pesos[i]);
            }
        }
    }

    // ORACULO DE CONSERVACAO, uma vez (FT_DIAG_CONSERVA): a soma da forca na
    // grade tem de igualar a dos marcadores -- SUM_faceta valor*h^DIM == SUM_k
    // F_k*ds_k -- porque o nucleo de Roma soma 1.  E' o teste do espalhamento
    // escalonado, independente das condicoes de contorno.  So' vale em serial
    // (sem franja) e antes do dp_sync.
    static int ja_conferiu = 0;
    if (!ja_conferiu && getenv("FT_DIAG_CONSERVA") != NULL) {
        ja_conferiu = 1;
        // Malha uniforme: h^DIM constante, tirado de um marcador dentro do dominio.
        real h = 0.0;
        for (int k = 0; k < n && h <= 0.0; k++) h = _h_em(sfd[0], pos[k]);
        real hd = 1.0; for (int d = 0; d < DIM; d++) hd *= h;
        for (int dim = 0; dim < DIM; dim++) {
            real lado_marc = 0.0;
            for (int k = 0; k < n; k++) lado_marc += F[k][dim] * ds[k];
            real lado_grade = 0.0;
            const hig_facet_snapshot *hfs = sfd_get_snapshot(sfd[dim]);
            for (int flid = 0; flid < hfs->n; flid++)
                lado_grade += dp_get_value(dpF[dim], flid) * hd;
            fprintf(stderr, "FT conserva dim=%d: grade=%.6e  marcadores=%.6e  "
                    "erro=%.3e\n", dim, lado_grade, lado_marc,
                    fabs(lado_grade - lado_marc));
        }
    }

    // ORACULO DE CONSERVACAO GLOBAL (FT_DIAG_CONSERVA_PAR), que e' o que separa
    // as regras de posse.  O oraculo antigo soma as facetas do dominio local,
    // franja inclusive, e por isso so' vale em serial.  Este soma apenas as
    // PROPRIAS e reduz entre os ranks: e' a forca que de fato existe na malha.
    // Tem de igualar SUM_k F_k*ds_k, que e' replicada e igual em todo rank.
    if (getenv("FT_DIAG_CONSERVA_PAR") != NULL) {
        static int ja = 0;
        if (!ja) {
            ja = 1;
            real h = 0.0;
            for (int k = 0; k < n && h <= 0.0; k++) h = _h_em(sfd[0], pos[k]);
            real hd = 1.0; for (int d = 0; d < DIM; d++) hd *= h;
            int rk = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
            for (int dim = 0; dim < DIM; dim++) {
                real marc = 0.0;
                for (int k = 0; k < n; k++) marc += F[k][dim] * ds[k];
                real loc = 0.0;
                const int nprop = fi_dp_num_proprias(dpF[dim]);
                for (int flid = 0; flid < nprop; flid++)
                    loc += dp_get_value(dpF[dim], flid) * hd;
                real glob = 0.0;
                MPI_Allreduce(&loc, &glob, 1, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
                if (rk == 0) {
                    int nt = 1; MPI_Comm_size(MPI_COMM_WORLD, &nt);
                    fprintf(stderr, "POSSE=%d np=%d dim=%d: grade(propria,global)="
                            "%+.8e  marcadores=%+.8e  erro=%.3e\n",
                            _posse_regra(), nt, dim, (double) glob, (double) marc,
                            (double) fabs(glob - marc));
                }
            }
        }
    }

    free(lids); free(pesos);
    free(pos); free(F); free(ds);

    // FRANJA: AINDA SEM CONSERTO, E A TENTATIVA OBVIA PIORA.  MEDIDO.
    //
    // O laco acima acumula com dp_add_value nas facetas que este rank encontra, e
    // `fi_suporte_facetas` busca no dominio local, que INCLUI FRANJA.  Uma
    // contribuicao que caia em faceta de franja nunca chega ao dono, e `dp_sync`
    // nao resolve -- ele manda dono->franja e SOBRESCREVE.
    //
    // A tentativa obvia era chamar `fi_reduz_franja`, a MESMA reducao que o
    // espalhamento do corpo rigido usa.  MEDIDA: a discrepancia do np=2 contra o
    // serial PIOROU de 1,5e-11 para 8,3e-6 -- pior que o defeito original.
    //
    // A razao e' que a reducao supoe que cada contribuicao foi feita UMA VEZ.  No
    // corpo rigido isso vale: o DMSwarm reparte os marcadores por posse
    // exclusiva.  Aqui a frente e' REPLICADA e a franja cria sobreposicao -- um
    // marcador na costura e' visto pelos DOIS ranks e espalhado pelos dois.  Antes
    // da reducao, a contribuicao do nao-dono caia na copia de franja e era
    // descartada pelo dp_sync: acidentalmente menos errada.  Com a reducao, ela e'
    // somada no dono e aquele marcador conta duas vezes.
    //
    // O conserto pede a mesma coisa que a adveccao pediu -- uma regra de POSSE,
    // decidida por medida e nao por analogia com o corpo rigido, cujo modelo de
    // distribuicao de marcadores e' outro.  Ate' entao, o espalhamento do FT
    // continua correto so' em np=1.
    if (_posse_regra() == 1)
        for (int dim = 0; dim < DIM; dim++) fi_reduz_franja(dpF[dim], MPI_COMM_WORLD);
    else
        for (int dim = 0; dim < DIM; dim++) dp_sync(dpF[dim]);
}

// --------------------------------------------------------------------------
// ADVECCAO ACOPLADA: os marcadores andam com o fluido.
//
// E' a mudanca 1 das tres do projeto, e o que separa o front-tracking da
// fronteira imersa rigida: onde o corpo rigido INTERPOLA u para calcular a forca
// que impoe u=0, aqui interpola u para MOVER o marcador.  Mesma interpolacao,
// destino oposto.
// --------------------------------------------------------------------------

// Contexto do callback: os dominios e os campos de velocidade da malha.
typedef struct {
    sim_facet_domain     **sfd;   // ns->sfdu
    distributed_property **dpu;   // ns->dpu
} interp_ctx;

// Diagnostico da PARTICAO DA UNIDADE (FT_DIAG_INTERP).  Os pesos de Roma tem de
// somar 1 em cada marcador; se o suporte estiver incompleto -- marcador perto da
// fronteira, ou suporte atravessando nivel de refino -- a soma cai abaixo de 1 e
// a velocidade interpolada sai PEQUENA DEMAIS, em silencio.  Foi esta checagem
// que pegou a forca "exatamente pela metade" no corpo rigido: ela compara com um
// valor absoluto conhecido, coisa que a conservacao nao faz.
static real _pior_desvio_unidade = 0.0;

// Maior |u| interpolado nos marcadores.  NAO e' enfeite: uma gota estatica que
// fica parada e' exatamente o que se veria se a interpolacao devolvesse ZERO em
// silencio -- area perfeita, deriva nula, circularidade constante.  Este numero
// separa "adveccao funcionando, velocidade genuinamente minuscula" de "adveccao
// que nao faz nada".  Comparar com o Vmax do campo euleriano.
static real _maior_u_marcador = 0.0;

// u(X) interpolado da malha escalonada.  O campo fica CONGELADO em u^n durante o
// passo (nao existe u em t+dt/2), entao `t` e' ignorado: o RK2 do ft_advecta
// vira avaliacao de ponto medio no ESPACO, e o esquema segue de primeira ordem
// no tempo -- coerente com o resto do acoplamento explicito.
// ORACULO PARA DECIDIR O DESENHO DA ADVECCAO PARALELA.
//
// Ha' duas formas de tornar a adveccao correta em np>1, e elas NAO sao
// equivalentes:
//
//   (A) cada rank interpola o que enxerga e os SOMATORIOS PARCIAIS sao reduzidos
//       entre os ranks.  Simples e barato.  So' esta' certo se a posse das
//       facetas do suporte for uma PARTICAO -- se dois ranks encontrarem a mesma
//       faceta (um como dona, outro como copia de franja), a soma conta duas
//       vezes.
//
//   (B) reunir as posicoes, cada rank interpola apenas os marcadores cujo
//       suporte ele tem INTEIRO, e as velocidades sao reunidas.  Correto por
//       construcao, mais caro, e falha se algum marcador nao for de ninguem.
//
// Argumentar qual e' o certo e' desnecessario: a particao da unidade mede.
// `_suporte` localiza facetas no dominio LOCAL, que inclui franja, entao a
// pergunta e' factual -- a soma GLOBAL dos pesos por marcador da' 1, mais que 1,
// ou menos que 1?
//
//   == 1   a posse e' particao: (A) esta' exatamente certa, e vence por ser mais
//          simples e mais barata.
//    > 1   ha' duplicacao de franja: (A) precisa de filtro de posse, ou (B).
//    < 1   ha' ponto de suporte que nao e' de ninguem: (B) tambem falharia ali,
//          e o problema e' a largura da franja, nao o desenho da reducao.
//
// Uma reducao por marcador seria inviavel em producao; aqui e' UMA reducao do
// vetor inteiro, so' quando FT_DIAG_UNIDADE_PAR esta' no ambiente.
// FILTRO DE POSSE: a interpolacao de velocidade correta em np>1.
//
// O DEFEITO QUE ISTO SUBSTITUI.  A versao por marcador (`_campo_da_malha`,
// removida quando esta passou a ser a unica em uso) desistia com
// `if (h <= 0.0) return;`, devolvendo velocidade ZERO para todo marcador fora do
// dominio do rank.  A frente e' replicada, entao cada rank advectava a sua copia
// com a velocidade que enxergava e congelava o resto: MEDIDO em np=2, 109 dos 252
// marcadores (43%) nao se moveram, numa faixa contigua -- o pedaco da bolha do
// outro rank.  Sem erro, sem aviso, codigo de saida zero.
//
// POR QUE O FILTRO, E NAO UMA SOMA.  Tres desenhos foram considerados e a
// PARTICAO DA UNIDADE decidiu entre eles, por medida e nao por argumento:
//
//   somar os parciais            duplica.  `_suporte` acha facetas no dominio
//                                local, que inclui FRANJA, e 172 das 504
//                                entradas eram encontradas COMPLETAS pelos dois
//                                ranks: a soma daria exatamente o dobro.
//   normalizar pela soma dos     nao serve.  12 entradas tinham um rank COMPLETO
//   pesos                        e outro PARCIAL; dividir mistura a interpolacao
//                                boa com uma estimativa parcial diferente.
//   filtro de posse (este)       so' contribui o rank cuja soma local de pesos e'
//                                1, isto e', que tem o suporte INTEIRO.  A media
//                                entre os ranks completos e' exata: eles somam
//                                sobre o MESMO conjunto de facetas.
//
// E a medida mostrou que ele e' sempre aplicavel aqui: `parcial orfa` = 0 e
// `sem ninguem` = 0, ou seja TODA entrada tem ao menos um rank completo (172 tem
// dois, 332 tem um).  Onde isso falhar nao existe resposta certa a dar, e por
// isso ABORTA -- devolver zero seria repor o defeito original com outra roupa.
static void _campo_da_malha_lote(const Point *x, int n, real t, void *vctx,
                                 Point *u)
{
    interp_ctx *ic = (interp_ctx *) vctx;
    (void) t;

    const int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    // val_loc[k*DIM+d] = interpolacao deste rank, SE ele tem o suporte inteiro.
    // cnt_loc[k*DIM+d] = 1 nesse caso, 0 caso contrario.
    real *val_loc = (real *) calloc((size_t) n * DIM, sizeof(real));
    real *cnt_loc = (real *) calloc((size_t) n * DIM, sizeof(real));

    for (int k = 0; k < n; k++) {
        const real h = _h_em(ic->sfd[0], x[k]);
        if (h <= 0.0) continue;            // este rank nao enxerga: contribui 0
        for (int dim = 0; dim < DIM; dim++) {
            const int m = fi_suporte_facetas(ic->sfd[dim], dim, x[k], h,
                                             lids, pesos, capac);
            real val = 0.0, soma = 0.0;
            for (int i = 0; i < m; i++) {
                val  += dp_get_value(ic->dpu[dim], lids[i]) * pesos[i];
                soma += pesos[i];
            }
            // SUPORTE INTEIRO?  E' esta comparacao que faz o filtro.  A tolerancia
            // e' folgada de proposito: os pesos de Roma somam 1 em aritmetica
            // exata, e o desvio observado em serie e' da ordem de 1e-14.
            if (fabs(soma - 1.0) < 1e-9) {
                val_loc[k*DIM + dim] = val;
                cnt_loc[k*DIM + dim] = 1.0;
            }
            const real desvio = fabs(soma - 1.0);
            if (desvio > _pior_desvio_unidade) _pior_desvio_unidade = desvio;
            if (fabs(val) > _maior_u_marcador) _maior_u_marcador = fabs(val);
        }
    }

    int ntasks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    real *val_g, *cnt_g;
    if (ntasks == 1) {
        // Serial: a reducao e' identidade.  Evita-la mantem o caminho de np=1
        // bit a bit o mesmo, que e' a condicao para as corridas ja' feitas
        // continuarem comparaveis.
        val_g = val_loc; cnt_g = cnt_loc;
    } else {
        val_g = (real *) malloc((size_t) n * DIM * sizeof(real));
        cnt_g = (real *) malloc((size_t) n * DIM * sizeof(real));
        MPI_Allreduce(val_loc, val_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(cnt_loc, cnt_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    }

    for (int k = 0; k < n; k++)
        for (int dim = 0; dim < DIM; dim++) {
            const real c = cnt_g[k*DIM + dim];
            if (c < 0.5) {
                // NENHUM rank tem o suporte inteiro deste marcador nesta direcao.
                // Nao ha' velocidade correta a devolver.  Abortar e' a unica
                // resposta honesta: zero aqui e' exatamente o defeito que este
                // codigo veio consertar, e ele nao se anuncia.
                int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if (rank == 0)
                    fprintf(stderr,
                        "_campo_da_malha_lote: marcador %d em (%.6f,%.6f), "
                        "direcao %d: NENHUM rank tem o suporte inteiro do nucleo. "
                        "A franja e' estreita demais para esta particao -- "
                        "ver a subsecao do paralelo no relatorio.\n",
                        k, (double) x[k][0], (double) x[k][1], dim);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            // Media entre os ranks completos.  Eles somam sobre o MESMO conjunto
            // de facetas, entao a media e' o proprio valor; dividir e' o que
            // desfaz a contagem multipla da franja.
            u[k][dim] = val_g[k*DIM + dim] / c;
        }

    free(lids); free(pesos);
    if (ntasks > 1) { free(val_g); free(cnt_g); }
    free(val_loc); free(cnt_loc);
}

static void _diag_unidade_paralela(ft_frente *frente, sim_facet_domain **sfd,
                                   int passo)
{
    const int n = ft_num(frente);
    Point *pos = (Point *) malloc((size_t) n * sizeof(Point));
    ft_posicoes(frente, pos);

    const int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    // soma_loc[k*DIM+dim] = soma dos pesos que ESTE rank encontra.
    real *soma_loc = (real *) calloc((size_t) n * DIM, sizeof(real));
    real *cont_loc = (real *) calloc((size_t) n * DIM, sizeof(real));

    for (int k = 0; k < n; k++) {
        real h = _h_em(sfd[0], pos[k]);
        if (h <= 0.0) continue;               // este rank nao enxerga: contribui 0
        for (int dim = 0; dim < DIM; dim++) {
            int m = fi_suporte_facetas(sfd[dim], dim, pos[k], h, lids, pesos, capac);
            real sw = 0.0;
            for (int i = 0; i < m; i++) sw += pesos[i];
            soma_loc[k*DIM + dim] = sw;
            cont_loc[k*DIM + dim] = (real) m;
        }
    }

    real *soma_g = (real *) malloc((size_t) n * DIM * sizeof(real));
    real *cont_g = (real *) malloc((size_t) n * DIM * sizeof(real));
    MPI_Allreduce(soma_loc, soma_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(cont_loc, cont_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);

    real menor = 1e300, maior = -1e300;
    int  n_zero = 0, n_acima = 0, n_abaixo = 0;
    real cmin = 1e300, cmax = -1e300;
    for (int k = 0; k < n * DIM; k++) {
        if (soma_g[k] < menor) menor = soma_g[k];
        if (soma_g[k] > maior) maior = soma_g[k];
        if (cont_g[k] < cmin) cmin = cont_g[k];
        if (cont_g[k] > cmax) cmax = cont_g[k];
        if (soma_g[k] == 0.0)            n_zero++;
        else if (soma_g[k] > 1.0 + 1e-9) n_acima++;
        else if (soma_g[k] < 1.0 - 1e-9) n_abaixo++;
    }
    // AS SOMAS LOCAIS, POR RANK.  A soma GLOBAL nao distingue "um rank com o
    // suporte inteiro" de "dois ranks com 0,7 e 0,3" -- as duas dao 1.  A
    // distincao decide se normalizar pela soma dos pesos e' exato:
    //
    //   todo rank com 0 ou 1  -> normalizar e' EXATO (u = SUM valor / SUM peso),
    //                            tanto na particao quanto na duplicacao
    //   algum rank PARCIAL    -> normalizar mistura interpolacao completa com
    //                            parcial e sai errado; e' preciso filtro de posse
    int loc_zero = 0, loc_um = 0, loc_parcial = 0;
    real pmin = 1e300, pmax = -1e300;
    for (int k = 0; k < n * DIM; k++) {
        const real v = soma_loc[k];
        if (v == 0.0)                   loc_zero++;
        else if (fabs(v - 1.0) < 1e-9)  loc_um++;
        else {
            loc_parcial++;
            if (v < pmin) pmin = v;
            if (v > pmax) pmax = v;
        }
    }

    // A PERGUNTA QUE DECIDE: alguma entrada tem um rank COMPLETO e outro PARCIAL?
    //
    // Parciais por si nao condenam normalizar pela soma dos pesos.  O que
    // condena e' a COEXISTENCIA na mesma entrada:
    //
    //   dois completos (1+1)                  -> 2u/2 = u          correto
    //   dois parciais que particionam (0,6+0,4) -> u/1 = u          correto
    //   um completo + um parcial (1 + 0,6)    -> mistura u com uma
    //                                            estimativa parcial  ERRADO
    //
    // Entao classifica-se cada entrada por QUANTOS ranks a tem completa e
    // quantos a tem parcial.  (A contagem "acima de 1" da soma global NAO
    // responde isso: 1+1 e 1+0,6 caem no mesmo balde.)
    real *comp_loc = (real *) calloc((size_t) n * DIM, sizeof(real));
    real *parc_loc = (real *) calloc((size_t) n * DIM, sizeof(real));
    for (int k = 0; k < n * DIM; k++) {
        const real v = soma_loc[k];
        if (v == 0.0)                  continue;
        if (fabs(v - 1.0) < 1e-9)      comp_loc[k] = 1.0;
        else                           parc_loc[k] = 1.0;
    }
    real *comp_g = (real *) malloc((size_t) n * DIM * sizeof(real));
    real *parc_g = (real *) malloc((size_t) n * DIM * sizeof(real));
    MPI_Allreduce(comp_loc, comp_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(parc_loc, parc_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);

    int b_misto = 0, b_dup = 0, b_unico = 0, b_part = 0, b_orfa = 0, b_nada = 0;
    for (int k = 0; k < n * DIM; k++) {
        const int c = (int) (comp_g[k] + 0.5), pa = (int) (parc_g[k] + 0.5);
        if      (c >= 1 && pa >= 1) b_misto++;   // <-- mata a normalizacao
        else if (c >= 2)            b_dup++;     // duplicacao exata
        else if (c == 1)            b_unico++;   // um dono so': limpo
        else if (pa >= 2)           b_part++;    // parciais que particionam
        else if (pa == 1)           b_orfa++;    // parcial sozinha: soma < 1
        else                        b_nada++;    // ninguem enxerga
    }

    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Cada rank imprime a SUA distribuicao: com poucos ranks isso e' mais
    // informativo que uma reducao, porque mostra o desequilibrio entre eles.
    fprintf(stderr, "UNIDADE-LOC passo %d rank %d: de %d entradas -- "
            "%d zeradas, %d com soma 1, %d PARCIAIS%s\n",
            passo, rank, n * DIM, loc_zero, loc_um, loc_parcial,
            loc_parcial > 0 ? "" : "  (nenhuma parcial)");
    if (loc_parcial > 0)
        fprintf(stderr, "UNIDADE-LOC passo %d rank %d: parciais entre "
                "%.12f e %.12f\n", passo, rank, (double) pmin, (double) pmax);

    if (rank == 0) {
        int ntasks = 1; MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
        fprintf(stderr, "UNIDADE-PAR passo %d np=%d: soma global dos pesos por "
                "marcador -- min %.12f  max %.12f | de %d entradas: "
                "%d zeradas, %d acima de 1, %d abaixo de 1 | "
                "pontos de suporte por entrada: %.0f a %.0f\n",
                passo, ntasks, (double) menor, (double) maior, n * DIM,
                n_zero, n_acima, n_abaixo, (double) cmin, (double) cmax);
        fprintf(stderr, "UNIDADE-CLS passo %d np=%d: por entrada -- "
                "MISTO(completo+parcial) %d | duplicada %d | dono unico %d | "
                "parciais que particionam %d | parcial orfa %d | sem ninguem %d"
                "   ==> normalizar pela soma dos pesos %s\n",
                passo, ntasks, b_misto, b_dup, b_unico, b_part, b_orfa, b_nada,
                (b_misto == 0 && b_orfa == 0) ? "E' EXATO" : "NAO SERVE");
    }
    free(comp_loc); free(parc_loc); free(comp_g); free(parc_g);
    free(pos); free(lids); free(pesos);
    free(soma_loc); free(cont_loc); free(soma_g); free(cont_g);
}

// Metricas de forma da gota, para o oraculo do B2 dinamico: area, deriva do
// centroide e circularidade (4*pi*A/P^2, que vale 1 para o circulo).  Uma gota
// estatica em equilibrio tem de manter as tres.
static void _metricas(const ft_frente *f, const Point centro_ini,
                      real *area, real *deriva, real *circ, real *defor)
{
    int n = ft_num(f);
    Point *p = (Point *) malloc((size_t) n * sizeof(Point));
    ft_posicoes(f, p);
    real cx = 0.0, cy = 0.0;
    for (int i = 0; i < n; i++) { cx += p[i][0]; cy += p[i][1]; }
    cx /= n; cy /= n;
    *deriva = sqrt((cx - centro_ini[0]) * (cx - centro_ini[0])
                 + (cy - centro_ini[1]) * (cy - centro_ini[1]));
    *area = ft_area(f);
    real P = ft_perimetro(f);
    *circ = (P > 0.0) ? 4.0 * M_PI * (*area) / (P * P) : 0.0;
    // DEFORMACAO COM SINAL, pelos semieixos equivalentes (momentos de area).
    // Troca de sinal a cada meio periodo -- a circularidade, sempre <= 1, nao
    // trocaria, e oscilaria no dobro da frequencia.
    real sa, sb;
    ft_semieixos(f, &sa, &sb);
    *defor = (sa + sb > 0.0) ? (sa - sb) / (sa + sb) : 0.0;
    free(p);
}

// ---------------------------------------------------------------------------
// FORCA BALANCEADA (FT_BALANCEADO=1) -- o teste da hipotese do balanco.
//
// A forma padrao espalha sigma*kappa*n pelo nucleo de Roma direto nas facetas.
// Esse operador NAO e' o gradiente discreto de escalar nenhum que a projecao
// produza, e o descasamento entre ele e o gradiente de pressao aparece como
// corrente parasita -- e' a hipotese levantada pela comparacao com o VOF.
//
// Aqui a forca e' montada como
//
//     F = sigma * kappa * grad(H),
//
// com H a fracao de area da gota por celula (da GEOMETRIA da frente, via
// ft_area_na_caixa) e `grad` o MESMO operador que higflow_final_velocity usa
// para corrigir a velocidade:
//
//     compute_center_p_left/right  +  compute_dpdx_at_point
//
// Com kappa constante -- que e' o caso do circulo -- isso e' exatamente
// grad(sigma*kappa*H), isto e', a forca PERTENCE a' imagem do gradiente
// discreto, e a pressao pode cancela-la termo a termo.  Se a hipotese estiver
// certa, as correntes parasitas caem sem que o salto de pressao piore.
// Ver Francois et al. (2006).
// ---------------------------------------------------------------------------

// Campo H por celula, preso ao dominio corrente.  Criado sob demanda.
static distributed_property *_dpH = NULL;

static void _preenche_indicadora(higflow_solver *ns, ft_frente *frente)
{
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    mp_mapper  *mp  = sd_get_domain_mapper(sdp);
    if (_dpH == NULL) _dpH = psd_create_property(ns->psdp);

    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point lo, hi;
        hig_get_lowpoint(c, lo);
        hig_get_highpoint(c, hi);
        real vol = 1.0;
        for (int d = 0; d < DIM; d++) vol *= (hi[d] - lo[d]);
        real a = ft_area_na_caixa(frente, lo, hi);
        dp_set_value(_dpH, clid, (vol > 0.0) ? a / vol : 0.0);
    }
    higcit_destroy(it);
    dp_sync(_dpH);
}

// Curvatura interpolada dos marcadores em `x`, por MEDIA PONDERADA com o nucleo
// de Roma de largura `h`.  So' e' consultada onde grad(H) e' nao nulo, a menos
// de uma celula da frente.
//
// A VERSAO ANTERIOR USAVA O MARCADOR MAIS PROXIMO, e isso custou uma corrida.
// No circulo era exato (kappa constante) e o teste passou espetacularmente.  Na
// ELIPSE, com kappa variavel, o campo do vizinho mais proximo e' descontinuo:
// salta quando a atribuicao troca de marcador.  Enquanto a gota ainda estava
// longe do circulo a forca fisica dominava e nada aparecia; assim que ela ficou
// quase circular (passo ~900, circ 0,9887) o artefato passou a dominar, a
// circularidade REVERTEU e a corrida foi a instabilidade.  A descontinuidade
// alimenta a adveccao, que move a frente, que muda a atribuicao -- realimentacao.
//
// Fora do suporte de qualquer marcador devolve 0; ali grad(H) tambem e' nulo.
static real _kappa_interp(const Point *pos, const Point *forca, int n,
                          real sigma, const Point x, real h)
{
    real num = 0.0, den = 0.0;
    for (int k = 0; k < n; k++) {
        real w = ft_delta_roma((pos[k][0] - x[0]) / h)
               * ft_delta_roma((pos[k][1] - x[1]) / h);
        if (w == 0.0) continue;
        // |forca| = sigma*kappa, por construcao de ft_forcas_tensao.
        real fx = forca[k][0], fy = forca[k][1];
        num += w * sqrt(fx * fx + fy * fy) / sigma;
        den += w;
    }
    return (den > 0.0) ? num / den : 0.0;
}

static void ft_espalha_tensao_balanceada(higflow_solver *ns, ft_frente *frente,
                                         real sigma,
                                         sim_facet_domain *sfd[DIM],
                                         distributed_property *dpF[DIM])
{
    _preenche_indicadora(ns, frente);

    int n = ft_num(frente);
    Point *pos = (Point *) malloc((size_t) n * sizeof(Point));
    Point *F   = (Point *) malloc((size_t) n * sizeof(Point));
    real  *ds  = (real  *) malloc((size_t) n * sizeof(real));
    ft_forcas_tensao(frente, sigma, pos, F, ds);

    sim_domain *sdp = psd_get_local_domain(ns->psdp);

    for (int dim = 0; dim < DIM; dim++) {
        mp_mapper *mu = sfd_get_domain_mapper(sfd[dim]);
        (void) mu;
        const hig_facet_snapshot *hfs = sfd_get_snapshot(sfd[dim]);
        for (int flid = 0; flid < hfs->n; flid++) {
            Point fcenter, fdelta;
            hfs_center(hfs, flid, fcenter);
            hfs_delta(hfs, flid, fdelta);

            // O MESMO par de amostragem que higflow_final_velocity usa para a
            // pressao -- e' o que torna a forca cancelavel pelo gradiente.
            real Hl = compute_center_p_left (sdp, fcenter, fdelta, dim, 0.5,
                                             _dpH, ns->stn);
            real Hr = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5,
                                             _dpH, ns->stn);
            real dHdx = compute_dpdx_at_point(fdelta, dim, 0.5, Hl, Hr);
            if (dHdx == 0.0) continue;          // longe da interface

            real hcel = _h_em(sfd[dim], fcenter);
            if (hcel <= 0.0) hcel = fdelta[dim];
            real kappa = _kappa_interp(pos, F, n, sigma, fcenter, hcel);
            dp_add_value(dpF[dim], flid, sigma * kappa * dHdx);
        }
    }

    free(pos); free(F); free(ds);
    for (int dim = 0; dim < DIM; dim++) dp_sync(dpF[dim]);
}


// ---------------------------------------------------------------------------
// FASE B4: a frente alimenta o CAMPO DE FRACAO do caminho multifasico.
//
// E' a terceira mudanca estrutural do projeto -- as propriedades saltam --, e a
// forma escolhida reusa tudo que o VOF ja' tem: em vez de advectar a fracao com
// PLIC, ela e' RECALCULADA da geometria da frente a cada passo, e dai' para
// frente rho(x), mu(x) e o momento de coeficiente variavel sao o mesmo codigo
// que o VOF usa.  A comparacao entre os dois passa a isolar exatamente a
// representacao da interface.
//
// Convencao: fracvol = 1 DENTRO da gota.  ft_area_na_caixa devolve a area da
// frente dentro da celula; dividida pelo volume da celula da' a fracao.
// ---------------------------------------------------------------------------
static void _preenche_fracvol(higflow_solver *ns, ft_frente *frente)
{
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper  *mp  = sd_get_domain_mapper(sdm);
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point lo, hi;
        hig_get_lowpoint(c, lo);
        hig_get_highpoint(c, hi);
        real vol = 1.0;
        for (int d = 0; d < DIM; d++) vol *= (hi[d] - lo[d]);
        real a = ft_area_na_caixa(frente, lo, hi);
        real f = (vol > 0.0) ? a / vol : 0.0;
        if (f < 0.0) f = 0.0;
        if (f > 1.0) f = 1.0;
        dp_set_value(ns->ed.mult.dpfracvol, clid, f);
    }
    higcit_destroy(it);
    dp_sync(ns->ed.mult.dpfracvol);
}

// Gancho: roda antes do preditor.  Nesse ponto `ns->dpu` e' u^n -- a velocidade
// final, JA' PROJETADA (discretamente livre de divergencia), do passo anterior.
// Entao a ordem e': move a frente com u^n, faz a cirurgia, e espalha a forca de
// tensao superficial nas posicoes NOVAS, que e' o que o preditor vai ver.
//
// Adveccao com campo livre de divergencia preserva a area fechada -- e' o
// oraculo mais barato de que a interpolacao esta' certa.
static void _aplica_tensao(higflow_solver *ns, void *vctx)
{
    ft_ctx *ctx = (ft_ctx *) vctx;

    if (ctx->advecta) {
        if (getenv("FT_DIAG_UNIDADE_PAR") != NULL)
            _diag_unidade_paralela(ctx->frente, ns->sfdu, ctx->passo);
        interp_ctx ic = { ns->sfdu, ns->dpu };
        // Adveccao em LOTE, com o filtro de posse: em np>1 a velocidade de um
        // marcador nao e' calculavel localmente, e uma reducao por marcador seria
        // inviavel.  Em np=1 o caminho e' bit a bit o mesmo de antes.
        ft_advecta_lote(ctx->frente, _campo_da_malha_lote, &ic, ns->par.t, ns->par.dt);
        ft_cirurgia(ctx->frente);

        static int cada = 0;
        if (cada == 0) { const char *c = getenv("FT_DIAG_CADA");
                         cada = (c != NULL) ? atoi(c) : 50; if (cada < 1) cada = 50; }
        if (getenv("FT_DIAG_INTERP") != NULL && ctx->passo % cada == 0) {
            real area, deriva, circ, defor;
            _metricas(ctx->frente, ctx->centro_ini, &area, &deriva, &circ, &defor);
            fprintf(stderr, "FT passo %5d: n=%4d  area=%.8f (dA/A=%.2e)  "
                    "deriva=%.3e  circ=%.6f  D=%+.6e  unidade_pior=%.2e  "
                    "max|u_marc|=%.3e\n",
                    ctx->passo, ft_num(ctx->frente), (double) area,
                    (double) fabs(area - ctx->area_ini) / ctx->area_ini,
                    (double) deriva, (double) circ, (double) defor,
                    (double) _pior_desvio_unidade,
                    (double) _maior_u_marcador);
        }
        // Despejo das POSICOES da frente (FT_DUMP_FRENTE=prefixo): um arquivo
        // por quadro com x y por marcador, fechando no primeiro.  E' o que as
        // figuras do relatorio desenham -- a forma, nao so' o numero.
        { const char *pre = getenv("FT_DUMP_FRENTE");
          static int cadaq = 0;
          if (cadaq == 0) { const char *c = getenv("FT_DUMP_CADA");
                            cadaq = (c != NULL) ? atoi(c) : 100; if (cadaq < 1) cadaq = 100; }
          if (pre != NULL && ctx->passo % cadaq == 0) {
              // Nome pelo PASSO, e nao por um contador de quadros: numa retomada
              // o contador privado voltaria a zero e sobrescreveria os quadros
              // anteriores a' queda.  O passo tambem da' o instante (t=passo*dt).
              char nome[512];
              snprintf(nome, sizeof nome, "%s_%07d.dat", pre, ctx->passo);
              FILE *fp = fopen(nome, "w");
              if (fp != NULL) {
                  int nn = ft_num(ctx->frente);
                  Point *pp = (Point *) malloc((size_t) nn * sizeof(Point));
                  ft_posicoes(ctx->frente, pp);
                  fprintf(fp, "x y\n");
                  for (int i = 0; i <= nn; i++)
                      fprintf(fp, "%.8f %.8f\n", pp[i % nn][0], pp[i % nn][1]);
                  fclose(fp); free(pp);
              }
          } }
        ctx->passo++;
    }

    // MULTIFASICO: a frente manda na fracao, e as propriedades saem dela.
    // Tem de acontecer ANTES de higflow_compute_viscosity/density_multiphase,
    // que e' exatamente onde o gancho foi posto no passo multifasico.
    if (ns->contr.flowtype == MULTIPHASE)
        _preenche_fracvol(ns, ctx->frente);

    if (ctx->balanceado)
        ft_espalha_tensao_balanceada(ns, ctx->frente, ctx->sigma,
                                     ns->sfdF, ns->dpFU);
    else
        ft_espalha_tensao_solver(ctx->frente, ctx->sigma, ns->sfdF, ns->dpFU);
}

//! Instala o front-tracking (tensao superficial) no solver.  `sigma` e' o
//! coeficiente; a frente e' o corpo ja' criado (ft_cria_circulo etc.).
//! extern "C": compilado como C++, mas o exemplo o declara com ligacao C.
// Passo em que a corrida COMECA.  Numa retomada nao e' zero, e isso importa para
// o despejo da frente: um contador privado reiniciado em zero sobrescreveria os
// arquivos dos quadros ja' escritos antes da queda.
static int _passo_inicial = 0;

extern "C" void front_tracking_passo_inicial(int passo)
{
    _passo_inicial = (passo > 0) ? passo : 0;
}

extern "C" void front_tracking_instala(higflow_solver *ns, ft_frente *frente,
                                       real sigma)
{
    ft_ctx *ctx = (ft_ctx *) malloc(sizeof *ctx);
    ctx->frente = frente;
    ctx->sigma  = sigma;
    const char *sa = getenv("FT_ADVECTA");
    ctx->advecta = (sa != NULL) ? atoi(sa) : 0;   // padrao: frente fixa (Laplace)
    ctx->passo   = _passo_inicial;
    { const char *s = getenv("FT_BALANCEADO"); ctx->balanceado = (s != NULL) ? atoi(s) : 0; }

    // Referencias para as metricas do oraculo: area e centroide iniciais.
    {
        int n = ft_num(frente);
        Point *p = (Point *) malloc((size_t) n * sizeof(Point));
        ft_posicoes(frente, p);
        real cx = 0.0, cy = 0.0;
        for (int i = 0; i < n; i++) { cx += p[i][0]; cy += p[i][1]; }
        ctx->centro_ini[0] = cx / n;
        ctx->centro_ini[1] = cy / n;
        for (int d = 2; d < DIM; d++) ctx->centro_ini[d] = 0.0;
        ctx->area_ini = ft_area(frente);
        free(p);
    }

    higflow_set_fronteira_imersa(ns, _aplica_tensao, ctx);
}

// ---------------------------------------------------------------------------
// FASE B4: a malha adaptada, pelo caminho VERIFICADO do remalhamento.
//
// Escreve um .amr MULTINIVEL a partir da geometria da frente, no mesmo formato
// que o modo criterio do example2d_SchaeferTurek usa -- e entao
// higflow_reconstroi_dominio o rele' pelo caminho normal de arranque.  E' a
// maquina que a fronteira imersa verificou (F1-F3, t8code e mtree, consenso de
// franja, transferencia por posicao), e nao a do example2d_DynamicMeshAdapt, que
// a propria suite marca quebrada em np=2.
//
// O CRITERIO E' O MESMO DO VOF: celula com 0<fracvol<1 e' SEMENTE de interface, e
// uma celula recebe o nivel L se estiver a menos de `limiar[L-1]` da semente mais
// proxima.  Aqui o fracvol da malha base sai da geometria da frente
// (ft_area_na_caixa) em vez de ser advectado -- mesma definicao, mesma conta.
//
// Os limiares vem da REGRA DAS CELULAS MINIMAS e sao DERIVADOS, nao digitados:
// o nivel mais fino cobre ao menos `cel_min` celulas finas de cada lado da
// interface, somada a margem que a discretizacao exige (a semente e' o CENTRO de
// uma celula base, a ate' h/2 da interface real).
// ---------------------------------------------------------------------------
extern "C" void ft_escreve_amr_criterio(higflow_solver *ns, ft_frente *frente,
                                        real lx, real ly, int nx, int ny,
                                        int niveis, int cel_min,
                                        const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny;
    const real h_base = (hx < hy) ? hx : hy;

    // Limiares derivados da regra (mesma formula de malha_adapt_limiares).
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    // Um h_base INTEIRO: o criterio decide por celula base, entao a borda da
    // banda e' incerta em h_base/2 pelo teste no centro E a celula se estende
    // outro h_base/2 alem.  Medido: com 0,5*h_base a banda dava 4,02 celulas.
    const real margem = h_base + h_fino;
    const real banda_fina = cel_min * h_fino + margem;
    real *thr = (real *) malloc((size_t) niveis * sizeof(real));
    for (int l = niveis - 1; l >= 0; l--) {
        thr[l] = banda_fina;
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        for (int k = 0; k < niveis - 1 - l; k++) thr[l] += cel_min * h_nivel;
    }

    // Sementes: centros das celulas base cortadas pela frente.
    const long ncel = (long) nx * ny;
    Point *sem = (Point *) malloc((size_t) ncel * sizeof(Point));
    long nsem = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            real lo[DIM], hi[DIM];
            lo[0] = i * hx;  hi[0] = (i + 1) * hx;
            lo[1] = j * hy;  hi[1] = (j + 1) * hy;
            for (int d = 2; d < DIM; d++) { lo[d] = 0.0; hi[d] = 1.0; }
            real f = ft_area_na_caixa(frente, lo, hi) / (hx * hy);
            if (f > 0.001 && f < 0.999) {
                sem[nsem][0] = lo[0] + 0.5 * hx;
                sem[nsem][1] = lo[1] + 0.5 * hy;
                nsem++;
            }
        }
    }

    // Tabela de nivel por celula base.
    signed char *tab = (signed char *) calloc((size_t) ncel, 1);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            real cx = (i + 0.5) * hx, cy = (j + 0.5) * hy;
            real d2 = 1e300;
            for (long k = 0; k < nsem; k++) {
                real a = cx - sem[k][0], b = cy - sem[k][1];
                real s = a * a + b * b;
                if (s < d2) d2 = s;
            }
            int nivel = 0;
            for (int l = niveis - 1; l >= 0; l--)
                if (d2 <= thr[l] * thr[l]) { nivel = l + 1; break; }
            tab[i + (long) j * nx] = (signed char) nivel;
        }
    }

    // ---------------------------------------------------------------
    // A ESTEIRA.  A banda de interface sozinha nao cobre o rastro que a
    // bolha deixa atras de si, e e' la' que a vorticidade vive.  O criterio
    // vira o MAXIMO entre os dois: banda geometrica (a regra das celulas
    // minimas) e vorticidade acima de limiar.  E' o criterio hibrido que a
    // F3 da fronteira imersa usou.
    //
    // HISTERESE (FT_VORT_HISTERESE): quem JA' estava refinado so' engrossa se
    // |omega| cair abaixo da METADE do limiar.  Sem isso a F3 mediu FLAPPING --
    // celulas cruzando o limiar seco a cada ciclo, a malha derivando
    // (67,9k -> 71,1k -> 64,9k) e o arrasto sujo.  Sair da malha fina tem de
    // ser mais dificil que entrar.
    // ---------------------------------------------------------------
    static signed char *tab_ant = NULL;
    static long ncel_ant = 0;
    {
        const char *e1 = getenv("FT_VORT1"), *e2 = getenv("FT_VORT2");
        const real v1 = (e1 != NULL) ? atof(e1) : 0.0;
        const real v2 = (e2 != NULL) ? atof(e2) : 0.0;
        const int  hist = (getenv("FT_VORT_HISTERESE") != NULL);
        // ns==NULL no ARRANQUE: a malha inicial sai so' da geometria da frente,
        // porque ainda nao ha' campo de velocidade para medir esteira.
        if (v1 > 0.0 && ns != NULL) {
            sim_domain *sdp = psd_get_local_domain(ns->psdp);
            sim_facet_domain *sfdu[DIM];
            for (int d = 0; d < DIM; d++)
                sfdu[d] = psfd_get_local_domain(ns->psfdu[d]);
            const hig_mesh_snapshot *hms = sd_get_snapshot(sdp);
            for (int clid = 0; clid < hms->n; clid++) {
                Point cc, cd;
                hms_center(hms, clid, cc);
                hms_delta(hms, clid, cd);
                Point pp;
                POINT_ASSIGN(pp, cc); pp[0] = cc[0] + cd[0];
                real vr = compute_facet_value_at_point(sfdu[1], cc, pp, 1.0, ns->dpu[1], ns->stn);
                pp[0] = cc[0] - cd[0];
                real vl = compute_facet_value_at_point(sfdu[1], cc, pp, 1.0, ns->dpu[1], ns->stn);
                POINT_ASSIGN(pp, cc); pp[1] = cc[1] + cd[1];
                real ut = compute_facet_value_at_point(sfdu[0], cc, pp, 1.0, ns->dpu[0], ns->stn);
                pp[1] = cc[1] - cd[1];
                real ub = compute_facet_value_at_point(sfdu[0], cc, pp, 1.0, ns->dpu[0], ns->stn);
                real om = fabs((vr - vl) / (2.0*cd[0]) - (ut - ub) / (2.0*cd[1]));
                int i = (int)(cc[0] / hx), j = (int)(cc[1] / hy);
                if (i < 0) i = 0; if (i >= nx) i = nx - 1;
                if (j < 0) j = 0; if (j >= ny) j = ny - 1;
                long q = i + (long) j * nx;
                signed char alvo = 0;
                if (v2 > 0.0 && om > v2) alvo = (signed char) niveis;
                else if (om > v1)        alvo = 1;
                if (hist && tab_ant != NULL && ncel_ant == ncel && tab_ant[q] > alvo) {
                    // manter exige so' METADE do limiar
                    if (om > 0.5 * v1) alvo = tab_ant[q];
                }
                if (alvo > tab[q]) tab[q] = alvo;   // MAXIMO com a banda
            }
        }
        if (ncel_ant != ncel) { free(tab_ant); tab_ant = NULL; }
        if (tab_ant == NULL) tab_ant = (signed char *) malloc((size_t) ncel);
        memcpy(tab_ant, tab, (size_t) ncel);
        ncel_ant = ncel;
    }

    // Escrita do .amr multinivel (mesmo formato do modo criterio).
    //
    // SO' O RANK 0 ESCREVE, E TODOS ESPERAM.  As sementes aqui NAO precisam de
    // reuniao entre ranks -- ao contrario do criterio do VOF, que as colhe do
    // dominio local, estas saem da FRENTE, que e' replicada e portanto completa
    // em todo rank.  O que precisa de conserto e' so' a escrita: sem a guarda
    // todos truncam o mesmo arquivo ao mesmo tempo, e sem a barreira alguem pode
    // rele-lo (em higflow_reconstroi_dominio, logo adiante) antes de estar pronto.
    int meurank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &meurank);
    if (meurank != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        free(tab); free(sem); free(thr);
        return;
    }
    FILE *f = fopen(caminho, "w");
    if (f == NULL) { perror(caminho); MPI_Barrier(MPI_COMM_WORLD);
                     free(tab); free(sem); free(thr); return; }
    fprintf(f, "0.0 %.10g 0.0 %.10g\n", (double) lx, (double) ly);
    long cont[8] = {0};
    for (long q = 0; q < ncel; q++)
        for (int l = 1; l <= niveis; l++) if (tab[q] >= l) cont[l]++;
    int niv_esc = 1;
    for (int l = 1; l <= niveis; l++) if (cont[l] > 0) niv_esc = l + 1;
    fprintf(f, "%d\n", niv_esc);
    // INDICE BASE 1, como os .amr simples do repositorio ("1 1 61 61").  Com
    // base 0 a contagem de celulas finas saia CERTA e as posicoes ERRADAS --
    // deslocadas de uma celula do nivel --, e a banda medida caia de 5 para 1,5
    // celulas.  Contagem certa com posicao errada e' o modo de falha que nao
    // aparece em nenhum total.
    fprintf(f, "%.10g %.10g\n1\n1 1 %d %d\n", (double) hx, (double) hy, nx, ny);
    for (int l = 1; l < niv_esc; l++) {
        int fator = 1 << l;
        fprintf(f, "%.10g %.10g\n%ld\n",
                (double) (hx / fator), (double) (hy / fator), cont[l]);
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                if (tab[i + (long) j * nx] >= l)
                    fprintf(f, "%d %d %d %d\n", fator*i + 1, fator*j + 1, fator, fator);
    }
    fclose(f);
    fprintf(stderr, "FT amr: %ld sementes, niveis=%d, limiares=", nsem, niv_esc-1);
    for (int l = 0; l < niveis; l++) fprintf(stderr, " %.5f", (double) thr[l]);
    fprintf(stderr, "  (>=%d celulas finas por lado)\n", cel_min);
    free(tab); free(sem); free(thr);
    MPI_Barrier(MPI_COMM_WORLD);   // o arquivo esta' completo no retorno
}
