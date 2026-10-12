// AS TRES CONDICOES DE FRONTEIRA: Dirichlet, Neumann e Robin.
//
// POR CELULAS FANTASMAS para Neumann e Robin, e amostrando a parede para
// Dirichlet -- a razao de cada escolha esta' no `contorno2d.hxx`.  O que este
// arquivo afirma e' que as tres funcionam, e o oraculo e' o mesmo em todas:
//
//   REPRODUCAO DE UM CAMPO ANALITICO ATE' A PAREDE.  Dado um campo `u` conhecido,
//   a condicao e' alimentada com o dado EXATO que `u` satisfaz naquela parede --
//   `u` para Dirichlet, `du/dn` para Neumann, `alfa*u + beta*du/dn` para Robin --
//   e o estencil tem de reproduzir `u` inclusive nas celulas da borda, onde parte
//   do suporte vem do fechamento.  Se a conta de (A,B) estiver errada, ou o sinal
//   da normal trocado, ou o `h` errado, aqui aparece.
//
// E O CAMPO VARIA NAS DUAS DIRECOES, com coeficientes distintos.  Com campo
// constante ou unidirecional, sinal de normal trocado e `h` errado ficam
// invisiveis -- e' a mesma razao pela qual nove defeitos de fechamento
// sobreviveram cinco anos na HiGTree.
//
// O CONTRASTE AQUI NAO E' ENTRE ORDENS, e sim entre DADO CERTO E DADO ERRADO: a
// mesma malha, o mesmo ajuste, alimentando a condicao com a derivada de OUTRO
// campo.  Se o resultado nao piorar, o fechamento nao esta' sendo usado -- e o
// verde do caso de cima nao significaria nada.
#include "../testing.h"

#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "estencil2d.hxx"
#include "malha2d.hxx"

using t8flow2d::Bloco;
using t8flow2d::ConstrutorEstencil;
using t8flow2d::Contorno;
using t8flow2d::Estencil;
using t8flow2d::Malha;
using t8flow2d::Mls;
using t8flow2d::Ponto2;
using t8flow2d::TipoContorno;

// O campo de teste e as suas derivadas, analiticos.  Grau 1 para o ajuste de
// ordem 1 reproduzir exatamente; coeficientes distintos por direcao.
static double u(const Ponto2 &p)     { return 1.0 + 2.0 * p[0] - 3.0 * p[1]; }
static double dudx(const Ponto2 &)   { return 2.0; }
static double dudy(const Ponto2 &)   { return -3.0; }

//! Outro campo, para o contraste do dado errado.
static double v(const Ponto2 &p)     { return 0.5 - 1.0 * p[0] + 7.0 * p[1]; }

//! Normal EXTERIOR no ponto da parede do quadrado unitario.
static void normal_em(const Ponto2 &p, double n[2])
{
  const double eps = 1e-9;
  n[0] = n[1] = 0.0;
  if (std::fabs(p[0]) < eps) n[0] = -1.0;
  else if (std::fabs(p[0] - 1.0) < eps) n[0] = 1.0;
  else if (std::fabs(p[1]) < eps) n[1] = -1.0;
  else if (std::fabs(p[1] - 1.0) < eps) n[1] = 1.0;
}

//! du/dn no ponto da parede, para o campo `f` de derivadas dadas.
static double dudn(const Ponto2 &p, double (*fx)(const Ponto2 &),
                   double (*fy)(const Ponto2 &))
{
  double n[2];
  normal_em(p, n);
  return fx(p) * n[0] + fy(p) * n[1];
}

static double aplica(const Estencil &e, const Malha &m, double (*f)(const Ponto2 &))
{
  double s = e.rhs;
  for (size_t k = 0; k < e.ids.size(); k++)
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        if (m.id_global(b, i) == e.ids[k]) { s += e.pesos[k] * f(m.centroide(b, i)); break; }
  return s;
}

//! Pior erro de reproducao de `f` perto da BORDA -- onde o fechamento entra.
//
// O ALVO NAO E' O CENTRO DA CELULA, e isso NAO e' detalhe: no proprio centro o
// ajuste movel e' DEGENERADO.  O peso do MLS e' `1/dist^2`, logo o ponto
// coincidente recebe `1/EPSMACH` = 1e12 e domina tudo -- MEDIDO: peso 1,0 na
// propria celula, 6,4e-11 nas outras, e `rhs` de -3,2e-11.  O estencil devolve o
// valor daquela celula e o contorno NAO INFLUENCIA NADA.
//
// A primeira versao deste arquivo media exatamente ali, e os casos passavam com
// 1e-16 sem exercitar o fechamento.  Quem denunciou foi o CONTRASTE: alimentar a
// condicao com a derivada de outro campo piorava o erro para so' 1e-10, numero
// pequeno demais para um dado grosseiramente errado.
//
// Daqui para a frente o alvo e' deslocado do centro -- que e' tambem o uso REAL do
// estencil: interpolar em faceta, em canto, em ponto de projecao.  No proprio no',
// o estencil e' a identidade, e isso esta' afirmado num caso proprio.
static double pior_perto_da_borda(const Malha &m, ConstrutorEstencil &ce,
                                  double (*f)(const Ponto2 &), long *quantas,
                                  long *com_rhs)
{
  double pior = 0.0;
  *quantas = 0;
  *com_rhs = 0;
  for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
    const Ponto2 c = m.centroide(0, i);
    const bool borda = (c[0] < 0.2 || c[0] > 0.8 || c[1] < 0.2 || c[1] > 0.8);
    if (!borda) continue;
    const Ponto2 h = m.lados_da_folha(0, i);
    // deslocamento assimetrico, dentro da celula, e distinto por direcao
    const Ponto2 x = {c[0] + 0.24 * h[0], c[1] + 0.16 * h[1]};
    const Estencil e = ce.em(x);
    (*quantas)++;
    if (e.rhs != 0.0) (*com_rhs)++;
    pior = std::fmax(pior, std::fabs(aplica(e, m, f) - f(x)));
  }
  return pior;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, 3}};
  Malha m(blocos, sc_MPI_COMM_SELF);
  const int ORDEM = 1;
  const int MAXP = 4 * Mls::min_pontos(ORDEM);

  // ---- a conta de (A, B), conferida contra a algebra -----------------------

  t_case("os_coeficientes_afins_batem_com_a_algebra");
  {
    const Ponto2 par = {0.0, 0.5};
    const double h = 0.125;
    double A = 0.0, B = 0.0;

    // NEUMANN: u_f = u_c + h*g  ->  A = 1, B = h*g
    Contorno::neumann([](const Ponto2 &) { return 3.0; }).afim(par, h, &A, &B);
    T_NEAR(A, 1.0, 1e-14, "A de Neumann");
    T_NEAR(B, h * 3.0, 1e-14, "B de Neumann");

    // DIRICHLET pela formula e' a reflexao ANTISSIMETRICA: A = -1, B = 2g.
    // (No estencil, Dirichlet amostra a parede em vez de refletir -- refletir
    // amplifica o erro do interior.  A formula vale como verificacao da algebra.)
    Contorno::dirichlet([](const Ponto2 &) { return 5.0; }).afim(par, h, &A, &B);
    T_NEAR(A, -1.0, 1e-14, "A de Dirichlet (limite beta=0)");
    T_NEAR(B, 10.0, 1e-14, "B de Dirichlet (limite beta=0)");

    // ROBIN com alfa=2, beta=1, g=4:  den = alfa*h/2 + beta
    const double al = 2.0, be = 1.0, g = 4.0;
    Contorno::robin(al, be, [](const Ponto2 &) { return 4.0; }).afim(par, h, &A, &B);
    const double den = al * h / 2.0 + be;
    T_NEAR(A, (be - al * h / 2.0) / den, 1e-14, "A de Robin");
    T_NEAR(B, g * h / den, 1e-14, "B de Robin");
  }

  t_case("condicao_degenerada_e_recusada");
  {
    bool lancou = false;
    try { Contorno::robin(0.0, 0.0, [](const Ponto2 &) { return 1.0; }); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "Robin com alfa=beta=0 passou sem reclamar");

    // E o denominador anulando com h: alfa*h/2 + beta == 0 exige beta negativo.
    lancou = false;
    try {
      double A = 0.0, B = 0.0;
      Contorno::robin(2.0, -0.125, [](const Ponto2 &) { return 1.0; })
          .afim({0.0, 0.5}, 0.125, &A, &B);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "denominador nulo em Robin passou sem reclamar");

    lancou = false;
    try {
      double A = 0.0, B = 0.0;
      Contorno::neumann([](const Ponto2 &) { return 1.0; }).afim({0.0, 0.5}, 0.0, &A, &B);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "h = 0 passou sem reclamar");
  }

  // ---- DIRICHLET ------------------------------------------------------------

  t_case("dirichlet_reproduz_o_campo_ate_a_parede");
  double pior_dir = 0.0;
  {
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::dirichlet([](const Ponto2 &p) { return u(p); }));
    long n = 0, r = 0;
    pior_dir = pior_perto_da_borda(m, ce, u, &n, &r);
    T_CHECK_MSG(n > 0 && r > 0, "%ld celulas de borda, %ld com rhs", n, r);
    T_NEAR(pior_dir, 0.0, 1e-10, "pior erro na borda com Dirichlet exato");
  }

  // ---- NEUMANN --------------------------------------------------------------

  t_case("neumann_reproduz_o_campo_ate_a_parede_pelo_fantasma");
  double pior_neu = 0.0;
  {
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::neumann([](const Ponto2 &p) { return dudn(p, dudx, dudy); }));
    long n = 0, r = 0;
    pior_neu = pior_perto_da_borda(m, ce, u, &n, &r);
    T_CHECK_MSG(n > 0 && r > 0, "%ld celulas de borda, %ld com rhs", n, r);
    T_NEAR(pior_neu, 0.0, 1e-10, "pior erro na borda com Neumann exato");
  }

  t_case("o_suporte_de_neumann_usa_fantasma_e_nao_amostra_a_parede");
  {
    // A distincao e' observavel, e tem de ser: Dirichlet amostra a PAREDE,
    // Neumann poe um FANTASMA do lado de fora.  Se os dois fizessem o mesmo, a
    // conta de (A,B) nunca seria exercitada.
    ConstrutorEstencil cd(m, ORDEM, MAXP);
    cd.contorno(Contorno::dirichlet([](const Ponto2 &p) { return u(p); }));
    ConstrutorEstencil cn(m, ORDEM, MAXP);
    cn.contorno(Contorno::neumann([](const Ponto2 &p) { return dudn(p, dudx, dudy); }));

    const Ponto2 canto = m.centroide(0, 0);     // celula do canto inferior esquerdo
    const ConstrutorEstencil::Suporte sd = cd.suporte_em(canto);
    const ConstrutorEstencil::Suporte sn = cn.suporte_em(canto);

    T_CHECK_MSG(sd.paredes > 0 && sd.fantasmas == 0,
                "Dirichlet: %ld paredes, %ld fantasmas (esperava paredes e nenhum fantasma)",
                sd.paredes, sd.fantasmas);
    T_CHECK_MSG(sn.fantasmas > 0 && sn.paredes == 0,
                "Neumann: %ld fantasmas, %ld paredes (esperava fantasmas e nenhuma parede)",
                sn.fantasmas, sn.paredes);

    // E o fantasma esta' FORA do dominio -- se estivesse dentro, seria uma celula
    // comum e a condicao nao estaria sendo imposta.
    long fora = 0;
    for (size_t k = 0; k < sn.pontos.size(); k++) {
      const Ponto2 &p = sn.pontos[k];
      if (p[0] < 0.0 || p[0] > 1.0 || p[1] < 0.0 || p[1] > 1.0) fora++;
    }
    T_CHECK_MSG(fora == sn.fantasmas,
                "%ld pontos fora do dominio para %ld fantasmas -- o fantasma tem de "
                "estar do lado de fora", fora, sn.fantasmas);
  }

  // O CONTRASTE DO NEUMANN: dado ERRADO tem de piorar.  Sem isto, "o erro e'
  // pequeno" nao distingue fechamento correto de fechamento ignorado.
  t_case("neumann_com_a_derivada_de_outro_campo_piora");
  {
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::neumann([](const Ponto2 &p) {
      return dudn(p, [](const Ponto2 &) { return -1.0; },    // derivadas de `v`
                     [](const Ponto2 &) { return 7.0; });
    }));
    long n = 0, r = 0;
    const double pior = pior_perto_da_borda(m, ce, u, &n, &r);
    T_CONTRASTE(pior_neu > 1e-16 ? pior_neu : 1e-16, pior, 1e3,
                "Neumann com a derivada certa contra a de outro campo");
    // E O ERRO ERRADO TEM DE SER GRANDE EM VALOR ABSOLUTO, nao so' maior.  Um
    // dado grosseiramente errado produzindo erro de 1e-10 significa que o
    // fechamento quase nao influencia -- foi exatamente esse numero que denunciou
    // o alvo degenerado no centro da celula.
    T_CHECK_MSG(pior > 1e-3,
                "derivada de outro campo deu erro de so' %.3e: o fechamento quase "
                "nao influencia o resultado, logo o caso de cima nao afirma o que "
                "parece afirmar", pior);
    printf("  (Neumann: derivada certa %.3e, derivada de outro campo %.3e)\n",
           pior_neu, pior);
  }

  // ---- ROBIN ----------------------------------------------------------------

  t_case("robin_reproduz_o_campo_ate_a_parede");
  double pior_rob = 0.0;
  {
    const double al = 2.0, be = 1.0;
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::robin(al, be, [al, be](const Ponto2 &p) {
      return al * u(p) + be * dudn(p, dudx, dudy);   // o dado EXATO que `u` satisfaz
    }));
    long n = 0, r = 0;
    pior_rob = pior_perto_da_borda(m, ce, u, &n, &r);
    T_CHECK_MSG(n > 0 && r > 0, "%ld celulas de borda, %ld com rhs", n, r);
    T_NEAR(pior_rob, 0.0, 1e-10, "pior erro na borda com Robin exato");
  }

  t_case("robin_com_beta_zero_concorda_com_dirichlet_na_reproducao");
  {
    // Robin com beta=0 E' Dirichlet (alfa*u = g).  Os CAMINHOS diferem -- um
    // reflete, o outro amostra --, mas a reproducao tem de valer nos dois.  Este
    // caso amarra a formula ao caso conhecido: se (A,B) estivesse errado, o
    // limite nao fecharia.
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::robin(1.0, 0.0, [](const Ponto2 &p) { return u(p); }));
    long n = 0, r = 0;
    const double pior = pior_perto_da_borda(m, ce, u, &n, &r);
    T_NEAR(pior, 0.0, 1e-10, "pior erro na borda com Robin de beta=0");
    printf("  (Dirichlet %.3e, Neumann %.3e, Robin %.3e, Robin beta=0 %.3e)\n",
           pior_dir, pior_neu, pior_rob, pior);
  }

  t_case("robin_com_alfa_e_beta_errados_piora");
  {
    // Dado certo, COEFICIENTES errados: a condicao passa a ser outra, e a
    // reproducao tem de cair.  Pega o caso em que (A,B) ignora alfa ou beta.
    const double al = 2.0, be = 1.0;
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::robin(5.0, 0.3, [al, be](const Ponto2 &p) {
      return al * u(p) + be * dudn(p, dudx, dudy);   // dado de OUTRA condicao
    }));
    long n = 0, r = 0;
    const double pior = pior_perto_da_borda(m, ce, u, &n, &r);
    T_CONTRASTE(pior_rob > 1e-16 ? pior_rob : 1e-16, pior, 1e3,
                "Robin com coeficientes certos contra errados");
    T_CHECK_MSG(pior > 1e-3,
                "coeficientes errados deram erro de so' %.3e -- o fechamento quase "
                "nao influencia", pior);
  }

  // A DEGENERACAO NO NO', afirmada em vez de evitada em silencio.  E' propriedade
  // do ajuste MOVEL, nao defeito: no ponto coincidente o peso `1/dist^2` vira
  // `1/EPSMACH` e domina.  Afirmar isto aqui e' o que impede alguem de voltar a
  // medir fechamento no centro da celula e achar que mediu.
  t_case("no_proprio_centro_da_celula_o_estencil_e_a_identidade");
  {
    ConstrutorEstencil ce(m, ORDEM, MAXP);
    ce.contorno(Contorno::neumann([](const Ponto2 &p) { return dudn(p, dudx, dudy); }));
    const Ponto2 c = m.centroide(0, 0);
    const Estencil e = ce.em(c);
    const long id = m.id_global(0, 0);

    double peso_dele = 0.0, peso_outros = 0.0;
    for (size_t k = 0; k < e.ids.size(); k++)
      if (e.ids[k] == id) peso_dele += e.pesos[k];
      else peso_outros += std::fabs(e.pesos[k]);

    T_NEAR(peso_dele, 1.0, 1e-8, "peso da propria celula no seu centro");
    T_CHECK_MSG(peso_outros < 1e-8,
                "os outros pesos somam %.3e no centro da celula -- esperava ~0, "
                "porque o ponto coincidente domina o ajuste movel", peso_outros);
    T_CHECK_MSG(std::fabs(e.rhs) < 1e-8,
                "rhs de %.3e no centro da celula -- esperava ~0, e e' por isso que "
                "medir fechamento ali nao mede nada", e.rhs);
  }

  t_case("sem_condicao_instalada_o_suporte_curto_e_recusado");
  {
    ConstrutorEstencil ce(m, 2, 40, 1);
    bool lancou = false;
    try { ce.em(m.centroide(0, 0)); }
    catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "suporte curto sem condicao passou sem reclamar");
  }

  return t_end();
}
