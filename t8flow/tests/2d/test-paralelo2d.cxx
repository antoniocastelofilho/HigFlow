// Arreio do PARALELO: o que o caminho serial nao exercita.
//
// POR QUE UM ARREIO SEPARADO, e nao casos espalhados nos outros: todos os casos
// daqui precisam de np >= 2 para medir algo, e em np=1 quase toda afirmacao
// deste arquivo e' trivialmente verdadeira.  Juntos, eles podem ser pedidos em
// np=2 e np=3 pelo arreio de execucao sem arrastar a suite inteira; espalhados,
// seriam casos que passam sem medir nada na maior parte das execucoes.
//
// A REGRA DESTE ARQUIVO: todo caso que afirma uma propriedade do paralelo tem de
// RECUSAR a execucao serial, ou dizer explicitamente que em np=1 nao mede nada.
// Um caso que passa em np=1 e nunca roda em np>1 e' pior que caso ausente: ele
// consta como cobertura.
//
// Cobre as clausulas X1 a X7 do contrato.

#include "../testing.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "estencil2d.hxx"
#include "malha2d.hxx"
#include "manufaturada2d.hxx"
#include "poisson2d.hxx"

using namespace t8flow2d;

namespace {

int NP = 1;
int EU = 0;

//! Soma um `long` entre todos os processos.
long soma(long v)
{
  if (NP <= 1) return v;
  long r = v;
  sc_MPI_Allreduce(&v, &r, 1, sc_MPI_LONG, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  return r;
}

double soma_d(double v)
{
  if (NP <= 1) return v;
  double r = v;
  sc_MPI_Allreduce(&v, &r, 1, sc_MPI_DOUBLE, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  return r;
}

//! Resolve o Poisson e devolve o erro -- o mesmo programa em qualquer np.
struct Caso {
  double e_inf = 0.0;
  double e_l2 = 0.0;
  double resid_op = 0.0;
  long celulas = 0;
  long truncados = 0;
};

Caso roda(const char *nome_sol, int nivel, double fracao, Esquema esq, int ordem = 3)
{
  const Manufaturada &sol = manufaturada(nome_sol);
  Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  if (fracao > 0.0) m.refina_aleatorio(fracao, 2026, nivel + 1);

  ContornoPorRegiao bc;
  bc.padrao(Contorno::dirichlet([&sol](const Ponto2 &p) { return sol.u(p); }));
  Poisson po(m, bc);
  po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
  po.esquema(esq);
  if (esq == Esquema::MlsGauss) po.mls(ordem, 60);
  po.resolve(Motor::Petsc, 1e-14);

  Caso c;
  c.e_inf = po.erro_inf([&sol](const Ponto2 &p) { return sol.u(p); });
  c.e_l2 = po.erro_l2([&sol](const Ponto2 &p) { return sol.u(p); });
  c.resid_op = po.residuo_da_exata([&sol](const Ponto2 &p) { return sol.u(p); });
  c.celulas = m.num_folhas_globais();
  c.truncados = (esq == Esquema::MlsGauss) ? po.truncados_pela_franja() : 0;
  return c;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &NP);
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &EU);

  if (EU == 0) printf("  (np = %d)\n", NP);

  // -------------------------------------------------------------------------
  // X1 -- o id do vizinho de franja
  // -------------------------------------------------------------------------

  t_case("o_id_do_vizinho_de_franja_nao_e_deslocamento_mais_indice");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);

    // TODO ID DE VIZINHO TEM DE EXISTIR na faixa global, e tem de ser conhecido
    // por este processo -- nas locais ou na franja.  A soma
    // `primeiro_id() + indice` aplicada a um indice de franja produz um id que
    // EXISTE (e' uma folha local deste processo), logo este caso nao a pegaria
    // sozinho; o que a pega e' a contagem de vizinhos distintos, abaixo.
    const long ng = m.num_folhas_globais();
    long fora = 0, desconhecidos = 0;
    std::set<long> vizinhos_vistos;

    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        for (int fa = 0; fa < 4; fa++)
          for (long v : m.ids_vizinhos(b, i, fa)) {
            if (v < 0 || v >= ng) fora++;
            if (m.indice_do_id(v) < 0) desconhecidos++;
            vizinhos_vistos.insert(v);
          }

    T_CHECK_MSG(fora == 0, "%ld ids de vizinho fora da faixa global [0, %ld)", fora, ng);
    T_CHECK_MSG(desconhecidos == 0,
                "%ld ids de vizinho que este processo nao conhece -- nem local nem "
                "na franja; o estencil os descartaria em silencio", desconhecidos);

    // O ORACULO QUE DISCRIMINA: a uniao dos vizinhos vistos por todos os
    // processos tem de cobrir TODAS as folhas da malha.  Com a soma errada, os
    // ids de franja viram ids locais repetidos, a uniao encolhe, e a contagem
    // acusa.  Em np=1 a uniao e' a malha inteira trivialmente, e o caso diz isso
    // em vez de fingir que mediu.
    const long distintos = soma(static_cast<long>(vizinhos_vistos.size()));
    if (NP > 1) {
      T_CHECK_MSG(distintos >= ng,
                  "a uniao dos vizinhos vistos pelos %d processos tem %ld ids "
                  "distintos (com repeticao entre processos) para %ld folhas -- se "
                  "os ids de franja estivessem sendo calculados por deslocamento, "
                  "colidiriam com os locais e a soma seria menor", NP, distintos, ng);
      if (EU == 0) printf("  (X1: %ld ids de vizinho vistos, %ld folhas globais)\n",
                          distintos, ng);
    } else {
      T_CHECK_MSG(m.num_fantasmas() == 0,
                  "np=1 nao deveria ter franja, e tem %ld folhas nela",
                  m.num_fantasmas());
      if (EU == 0) printf("  (X1: np=1, sem franja -- este caso NAO mede o paralelo)\n");
    }
  }

  // -------------------------------------------------------------------------
  // X2 -- as camadas da franja
  // -------------------------------------------------------------------------

  t_case("a_franja_tem_as_camadas_pedidas_e_exigir_e_idempotente");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    T_CHECK_MSG(m.camadas_de_franja() == 1,
                "a franja nasce com uma camada (a nativa do t8code), e nasceu com %d",
                m.camadas_de_franja());

    const long f1 = m.num_fantasmas();
    m.exige_camadas(3);
    const long f3 = m.num_fantasmas();
    T_CHECK_MSG(m.camadas_de_franja() == 3, "pedi 3 camadas e a malha diz %d",
                m.camadas_de_franja());

    // IDEMPOTENTE: pedir menos do que ja' existe nao desfaz.  Desfazer seria o
    // comportamento mais plausivel de uma implementacao descuidada, e quebraria
    // quem pediu mais antes.
    m.exige_camadas(2);
    T_CHECK_MSG(m.camadas_de_franja() == 3,
                "pedir 2 depois de 3 derrubou para %d -- exigir nao pode desfazer",
                m.camadas_de_franja());
    T_CHECK_MSG(m.num_fantasmas() == f3,
                "a franja mudou de %ld para %ld ao pedir MENOS camadas",
                f3, m.num_fantasmas());

    if (NP > 1) {
      T_CHECK_MSG(f3 > f1,
                  "tres camadas (%ld folhas) tinham de trazer mais que uma (%ld)",
                  f3, f1);
      if (EU == 0) printf("  (X2: franja de %ld folhas com 1 camada, %ld com 3)\n", f1, f3);
    } else {
      T_CHECK_MSG(f1 == 0 && f3 == 0,
                  "np=1 nao tem franja, e tem %ld com 1 camada e %ld com 3", f1, f3);
      if (EU == 0) printf("  (X2: np=1, sem franja -- este caso NAO mede o paralelo)\n");
    }
  }

  // -------------------------------------------------------------------------
  // X3 -- camadas = aneis + 1
  // -------------------------------------------------------------------------

  t_case("com_camadas_de_menos_o_suporte_trunca_e_com_as_certas_nao");
  {
    // O CONTRASTE E' O CASO.  Com `aneis` camadas o suporte trunca na fronteira
    // de particao; com `aneis + 1` nao.  A diferenca de um foi medida, nao
    // deduzida: o ponto de quadratura esta' SOBRE a faceta, logo a celula que o
    // contem pode ser de franja, e dela o suporte ainda precisa abrir os aneis.
    const int ANEIS = 2;

    Malha curta({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    curta.exige_camadas(ANEIS);              // DE MENOS, de proposito
    ConstrutorEstencil cc(curta, 3, 60, ANEIS);
    cc.dirichlet([](const Ponto2 &) { return 0.0; });

    Malha certa({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    certa.exige_camadas(ANEIS + 1);
    ConstrutorEstencil ce(certa, 3, 60, ANEIS);
    ce.dirichlet([](const Ponto2 &) { return 0.0; });

    // Consulta nos CENTROS DAS FACETAS das folhas locais, que e' onde o esquema
    // consulta -- e onde a celula que contem o ponto pode ser a de la'.
    long pedidos = 0;
    for (long b = 0; b < certa.num_blocos(); b++)
      for (long i = 0; i < certa.num_folhas_do_bloco(b); i++) {
        const Ponto2 c = certa.centroide(b, i);
        const Ponto2 h = certa.lados_da_folha(b, i);
        for (int d = 0; d < 2; d++)
          for (int lado = 0; lado < 2; lado++) {
            Ponto2 x = c;
            x[d] += (lado == 0 ? -1.0 : 1.0) * h[d] / 2.0;
            try { ce.suporte_em(x); } catch (const std::exception &) { /* conta abaixo */ }
            try { cc.suporte_em(x); } catch (const std::exception &) { }
            pedidos++;
          }
      }

    const long t_curta = soma(cc.truncados_pela_franja());
    const long t_certa = soma(ce.truncados_pela_franja());

    T_CHECK_MSG(t_certa == 0,
                "com %d camadas (aneis+1) ainda houve %ld suportes truncados -- os "
                "pesos perto da fronteira de particao nao seriam os do serial",
                ANEIS + 1, t_certa);
    // A REDUCAO VEM ANTES DE QUALQUER `if` POR RANK, SEMPRE.
    //
    // `soma` e' coletiva: se apenas o rank 0 a chama, os outros nunca entram no
    // `Allreduce` e o programa TRAVA -- nao falha, trava.  Num arreio de testes
    // isso aparece como teste LENTO e nao como teste errado, e foi o que
    // aconteceu aqui: duas execucoes de dez minutos antes de ser notado.
    const long total_pedidos = soma(pedidos);
    if (NP > 1) {
      T_CONTRASTE(t_certa + 1, t_curta + 1, 2.0,
                  "truncamentos com aneis camadas contra aneis+1");
      if (EU == 0)
        printf("  (X3: %ld consultas; truncados: %ld com %d camadas, %ld com %d)\n",
               total_pedidos, t_curta, ANEIS, t_certa, ANEIS + 1);
    } else if (EU == 0) {
      printf("  (X3: np=1, sem franja -- este caso NAO mede o paralelo)\n");
    }
  }

  // -------------------------------------------------------------------------
  // X4 -- borda da franja nao e' borda do dominio
  // -------------------------------------------------------------------------

  t_case("a_borda_da_franja_nao_e_confundida_com_borda_do_dominio");
  {
    // AS DUAS SITUACOES DAO LISTA DE VIZINHOS VAZIA, e confundi-las faria o
    // esquema aplicar condicao de contorno ao longo de toda a fronteira de
    // particao -- no MEIO do dominio, e em silencio.
    //
    // O ORACULO NAO PRECISA DE NUMERO NENHUM, e a primeira versao deste caso
    // precisava: ela previa 68 paredes numa malha 16x16 e mediu 248, porque o
    // estencil de dois aneis alcanca celulas a distancia 2 do contorno e cada uma
    // ve as paredes dela.  O numero previsto estava errado e o caso reprovou --
    // corretamente, por motivo errado.
    //
    // O criterio certo e' GEOMETRICO e local: todo ponto de parede que um
    // estencil registra tem de cair na fronteira do QUADRADO.  Uma faceta de
    // fronteira de PARTICAO tratada como contorno produziria um ponto de parede
    // no interior, e o caso o localiza -- em qualquer np, sem referencia externa.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    m.exige_camadas(3);
    ConstrutorEstencil ce(m, 3, 60, 2);
    ce.dirichlet([](const Ponto2 &) { return 0.0; });

    long paredes = 0, no_interior = 0;
    Ponto2 pior{-1.0, -1.0};
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const ConstrutorEstencil::Suporte s = ce.suporte_em(m.centroide(b, i));
        for (size_t k = 0; k < s.ids.size(); k++) {
          if (s.ids[k] >= 0) continue;          // id negativo marca a parede
          paredes++;
          const Ponto2 &x = s.pontos[k];
          const bool na_fronteira =
              std::fabs(x[0]) < 1e-12 || std::fabs(x[0] - 1.0) < 1e-12 ||
              std::fabs(x[1]) < 1e-12 || std::fabs(x[1] - 1.0) < 1e-12;
          if (!na_fronteira) { no_interior++; pior = x; }
        }
      }

    const long total = soma(paredes);
    const long ruins = soma(no_interior);
    T_CHECK_MSG(ruins == 0,
                "%ld pontos de parede caíram no INTERIOR do dominio (np=%d); um "
                "deles em (%.4f, %.4f).  Isso e' a borda da FRANJA sendo tratada "
                "como borda do DOMINIO", ruins, NP, pior[0], pior[1]);
    T_CHECK_MSG(total > 0,
                "nenhuma parede registrada -- o caso nao mediu nada, porque o "
                "fechamento de contorno nao foi exercitado");
    if (EU == 0)
      printf("  (X4: %ld paredes vistas em np=%d, %ld no interior)\n", total, NP, ruins);
  }

  // -------------------------------------------------------------------------
  // X5, X6 -- as reducoes, e a solucao independente de np
  // -------------------------------------------------------------------------

  t_case("as_normas_e_o_residuo_sao_os_mesmos_em_qualquer_np");
  {
    // GRAU 1 COM FLUXO DE SEGUNDA ORDEM: o esquema e' exato, logo o erro e' de
    // maquina em qualquer np -- e o residuo do OPERADOR tambem, que e' o que
    // separa o esquema do solver.
    //
    // Sem a reducao, cada processo devolveria o pior erro da SUA particao.  Isso
    // da' um numero MENOR que o verdadeiro, e num relatorio de verificacao e'
    // pior que um erro grande: e' um erro grande que se apresenta como pequeno.
    const Caso c = roda("grau1", 4, 0.0, Esquema::Fluxo2);
    T_BELOW(c.e_inf, 1e-10, "erro inf do grau1 com fluxo2");
    T_BELOW(c.resid_op, 1e-10, "residuo do operador do grau1 com fluxo2");
    T_CHECK_MSG(c.celulas == 256,
                "contagem GLOBAL de celulas deu %ld em np=%d; o esperado e' 256 "
                "independente da particao", c.celulas, NP);
    if (EU == 0)
      printf("  (X5: grau1 fluxo2, np=%d: inf %.3e, residuo do operador %.3e)\n",
             NP, c.e_inf, c.resid_op);
  }

  t_case("a_solucao_e_a_mesma_em_qualquer_np_na_malha_uniforme");
  {
    const Caso c = roda("grau3", 4, 0.0, Esquema::MlsGauss);
    T_BELOW(c.e_inf, 1e-10, "erro inf do grau3 com MlsGauss em malha uniforme");
    T_BELOW(c.e_l2, 1e-10, "erro L2 do grau3 com MlsGauss em malha uniforme");
    T_CHECK_MSG(c.truncados == 0,
                "a franja truncou %ld suportes -- com truncamento os pesos nao "
                "sao os do serial, e a solucao passa a DEPENDER da particao",
                c.truncados);
    T_CHECK_MSG(c.celulas == 256, "celulas globais %ld, esperado 256", c.celulas);
    if (EU == 0)
      printf("  (X6 uniforme: np=%d, inf %.3e, L2 %.3e, truncados %ld)\n",
             NP, c.e_inf, c.e_l2, c.truncados);
  }

  t_case("a_solucao_e_a_mesma_em_qualquer_np_na_malha_adaptada");
  {
    // O CASO QUE IMPORTA: e' na malha adaptada que a transicao 2:1 cruza
    // fronteira de particao, e e' ali que a sub-faceta, a correcao transversal e
    // a franja tem de concordar todas de uma vez.
    const Caso c = roda("grau2", 4, 0.4, Esquema::MlsGauss);
    T_BELOW(c.e_inf, 1e-10, "erro inf do grau2 com MlsGauss em malha ADAPTADA");
    T_BELOW(c.e_l2, 1e-10, "erro L2 do grau2 com MlsGauss em malha ADAPTADA");
    T_CHECK_MSG(c.truncados == 0, "a franja truncou %ld suportes na malha adaptada",
                c.truncados);
    if (EU == 0)
      printf("  (X6 adaptada: np=%d, %ld celulas, inf %.3e, truncados %ld)\n",
             NP, c.celulas, c.e_inf, c.truncados);
  }

  // -------------------------------------------------------------------------
  // X7 -- a reparticao
  // -------------------------------------------------------------------------

  t_case("a_reparticao_preserva_a_malha_e_refaz_a_franja");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    m.refina_aleatorio(0.4, 2026, 5);
    m.exige_camadas(3);

    // O QUE A REPARTICAO NAO PODE MUDAR: a malha.  O que ela muda e' quem fica
    // com o que -- e isso o oraculo nao afirma, porque a decisao e' da curva de
    // preenchimento do t8code.
    const long ng = m.num_folhas_globais();
    const int nmin = m.nivel_minimo(), nmax = m.nivel_maximo();
    const long transicao = m.num_folhas_em_transicao();
    double area = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) area += m.area(b, i);
    const double area_antes = soma_d(area);
    const long locais_antes = m.num_folhas();

    m.reparticiona();

    double area2 = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) area2 += m.area(b, i);

    T_CHECK_MSG(m.num_folhas_globais() == ng,
                "a reparticao mudou a contagem global de %ld para %ld",
                ng, m.num_folhas_globais());
    const double area_depois = soma_d(area2);
    T_NEAR(area_depois, area_antes, 1e-12, "area total depois da reparticao");
    T_CHECK_MSG(m.nivel_minimo() == nmin && m.nivel_maximo() == nmax,
                "niveis mudaram de [%d, %d] para [%d, %d]",
                nmin, nmax, m.nivel_minimo(), m.nivel_maximo());
    T_CHECK_MSG(m.num_folhas_em_transicao() == transicao,
                "folhas em transicao mudaram de %ld para %ld",
                transicao, m.num_folhas_em_transicao());

    // A FRANJA FOI REFEITA.  Esquecer isso deixaria a franja velha apontando para
    // folhas que mudaram de dono -- e isso nao falha, apenas responde errado.  O
    // oraculo: todo id de vizinho continua conhecido por este processo.
    long desconhecidos = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        for (int fa = 0; fa < 4; fa++)
          for (long v : m.ids_vizinhos(b, i, fa))
            if (m.indice_do_id(v) < 0) desconhecidos++;
    // A REDUCAO FICA NUMA VARIAVEL, e nao nos argumentos da macro: os argumentos
    // de formatacao so' sao avaliados quando a condicao FALHA, logo uma reducao
    // ali seria chamada por uns ranks e nao por outros -- e travaria.
    const long tot_desc = soma(desconhecidos);
    T_CHECK_MSG(tot_desc == 0,
                "%ld ids de vizinho desconhecidos depois da reparticao -- a franja "
                "nao foi refeita", tot_desc);

    if (EU == 0)
      printf("  (X7: %ld folhas globais, %ld locais antes e %ld depois, area %.6f)\n",
             ng, locais_antes, m.num_folhas(), area_antes);
  }

  // -------------------------------------------------------------------------
  // X8, X9 -- a saturacao da carga, e a ordem dos ids de franja
  // -------------------------------------------------------------------------

  t_case("a_carga_da_troca_nao_saturou");
  {
    // O canal e' um array de itens de tamanho FIXO, e a vizinhanca publicada por
    // folha e' cortada quando nao cabe.  Cortar e' inevitavel; cortar em silencio
    // nao e' -- foi a saturacao invisivel que fez quatro camadas de franja ainda
    // deixarem niveis ausentes, e mandou a investigacao atras da camada errada.
    //
    // As malhas da suite tem de caber.  Se um dia nao couberem, este caso diz
    // isso em vez de deixar o sintoma aparecer como erro de truncamento.
    Malha u({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    u.exige_camadas(3);
    const long sat_u = soma(u.saturacoes_da_franja());
    T_CHECK_MSG(sat_u == 0,
                "a carga saturou %ld vezes na malha UNIFORME de nivel 4 com 3 "
                "camadas", sat_u);

    Malha a({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    a.refina_aleatorio(0.4, 2026, 5);
    a.exige_camadas(3);
    const long sat = soma(a.saturacoes_da_franja());
    T_CHECK_MSG(sat == 0,
                "a carga saturou %ld vezes na malha ADAPTADA -- em malha adaptada a "
                "vizinhanca de uma celula grossa tem muitos vizinhos finos, e e' "
                "aqui que a carga estoura primeiro", sat);
    if (EU == 0) printf("  (X8: saturacoes %ld, np=%d)\n", sat, NP);
  }

  t_case("o_mapa_de_id_para_geometria_nao_depende_da_particao");
  {
    // X9.  O ORACULO QUE PEGOU O DEFEITO MAIS CARO desta frente.
    //
    // A documentacao do t8code e' explicita: `leaf_face_neighbors` numera
    // `0 .. num_local-1` para as locais e `num_local .. num_local+num_ghosts-1`
    // para as da camada 1.  Ordenar a franja por id -- determinista, e aparentemente
    // a escolha limpa -- faz o indice do terceiro fantasma designar OUTRA folha.
    //
    // O que o caso afirma: a correspondencia `id -> (centro, nivel)` nao depende
    // da particao.  Ela e' uma propriedade da MALHA, e um id que muda de lugar
    // conforme o numero de processos nao e' um id.
    //
    // E a conferencia e' LOCAL, sem precisar comparar com uma execucao serial: a
    // posicao e o nivel determinam o id na curva de preenchimento, logo o centro
    // reconstruido do id tem de cair dentro da celula que o id designa.  Um id
    // trocado aponta para outra celula, e a verificacao falha.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 4}}, sc_MPI_COMM_WORLD);
    m.refina_aleatorio(0.4, 2026, 5);
    m.exige_camadas(3);

    long errados = 0;
    long id_ruim = -1;
    for (long x = 0; x < m.num_folhas() + m.num_fantasmas(); x++) {
      const long id = m.id_global_do_indice(x);
      const Ponto2 c = m.centroide_do_indice(x);

      // `id_do_ponto` percorre locais e franja pela regra semiaberta e devolve o
      // id da celula que CONTEM o ponto.  Para o centro de uma celula, isso tem
      // de ser o id dela -- e e' um caminho de codigo independente de
      // `id_global_do_indice`, que e' o que faz disto um oraculo e nao uma
      // tautologia.
      if (m.id_do_ponto(c) != id) { errados++; id_ruim = id; }
    }
    const long tot_err = soma(errados);
    T_CHECK_MSG(tot_err == 0,
                "%ld folhas cujo id nao corresponde a' geometria (uma delas id %ld) "
                "com np=%d -- o indice de franja esta' designando outra folha",
                tot_err, id_ruim, NP);

    // E O CONTRASTE: a malha tem de TER franja em np>1, senao o caso acima
    // percorreu apenas folhas locais e nao mediu a franja.
    const long tot_franja = soma(m.num_fantasmas());
    const long tot_folhas = soma(m.num_folhas() + m.num_fantasmas());
    if (NP > 1) {
      T_CHECK_MSG(tot_franja > 0,
                  "nenhuma folha de franja em np=%d -- o caso nao mediu nada", NP);
      if (EU == 0)
        printf("  (X9: %ld folhas conferidas, %ld na franja)\n",
               tot_folhas, tot_franja);
    } else if (EU == 0) {
      printf("  (X9: np=1, sem franja -- este caso NAO mede o paralelo)\n");
    }
  }

  t_case("o_desbalanceamento_e_medido_e_a_reparticao_o_zera");
  {
    // O QUE ESTE CASO AFIRMA NAO E' QUE A PARTICAO ESTEJA BOA, e sim que ela
    // seja MENSURAVEL -- porque as medidas deste projeto nao fecham e a decisao
    // de reparticionar depende delas.
    //
    // O t8code ja' tenta balancear (`set_balance` com `no_repartition = 0`, e a
    // documentacao promete "load-balanced per default"), mas:
    //
    //     refino em CAIXA, nivel 6, np 2/4/8     1,00 (perfeito)
    //     refino ALEATORIO, nivel 6, np 8        1,11
    //     `bench/escalabilidade`, outra config   2,8
    //
    // Por isso a reparticao e' EXPLICITA: um segundo mecanismo automatico em
    // cima de um primeiro que nao se entende seria redundante num caso e
    // insuficiente noutro.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 5}}, sc_MPI_COMM_WORLD);
    m.refina_aleatorio(0.4, 2026, 6);

    const double antes = m.desbalanceamento();
    T_CHECK_MSG(antes >= 1.0,
                "desbalanceamento %.3f e' menor que 1 -- a razao max/min nao pode "
                "ser", antes);

    m.reparticiona();
    const double depois = m.desbalanceamento();

    // A REPARTICAO LEVA A 1,0 EXATAMENTE quando as folhas dividem pelo numero de
    // processos, e a no maximo `1 + 1/minimo` quando nao dividem -- a curva de
    // preenchimento distribui o resto uma folha por processo.
    const long por_processo = m.num_folhas_globais() / NP;
    const double teto = 1.0 + (por_processo > 0 ? 1.0 / por_processo : 1.0);
    T_CHECK_MSG(depois <= teto + 1e-12,
                "depois de reparticionar o desbalanceamento e' %.6f, acima do teto "
                "%.6f que a divisao inteira admite", depois, teto);
    T_CHECK_MSG(depois <= antes + 1e-12,
                "reparticionar PIOROU o balanco: %.3f -> %.3f", antes, depois);

    if (NP > 1) {
      if (EU == 0)
        printf("  (Y2: desbalanceamento %.4f -> %.4f, teto %.4f, np=%d)\n",
               antes, depois, teto, NP);
    } else {
      T_CHECK_MSG(antes == 1.0 && depois == 1.0,
                  "np=1 deveria dar desbalanceamento 1,0 exato; deu %.3f e %.3f",
                  antes, depois);
      if (EU == 0)
        printf("  (Y2: np=1, desbalanceamento 1,0 por construcao -- este caso NAO "
               "mede o paralelo)\n");
    }
  }

  return t_end();
}
