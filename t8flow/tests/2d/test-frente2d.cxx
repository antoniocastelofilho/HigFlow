// A FRENTE (fase B1): advecao cinematica, cirurgia, e as formas fechadas.
//
// OS ORACULOS SAO OS QUE A SECAO COMPANHEIRA LISTA para o lado bidimensional:
// curvatura contra circulo e elipse, e a reversao de adveccao do vortice unico
// de LeVeque.  Nenhum deles depende do escoamento -- e' o ponto de B1: o que
// se mede aqui e' o transporte da frente e mais nada.
//
// A CURVATURA DO CIRCULO E' ORACULO DE EXATIDAO, nao de ordem: ela sai do
// CIRCUNRAIO de tres marcadores, e para tres pontos sobre um circulo o
// circunraio E' o raio -- com tres marcadores ou com trezentos.  Uma diferenca
// finita parametrica daria so' segunda ordem, e o teste teria de aceitar erro.
//
// E CADA AFIRMACAO DE EXATIDAO TEM O PAR QUE NAO E' EXATO: a elipse, onde a
// curvatura varia e a formula converge em vez de acertar.  Sem o par, "o erro
// deu 1e-16" nao distingue a formula certa de uma que devolve uma constante.

#include "../testing.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "frente2d.hxx"

using namespace t8flow2d;

namespace {
constexpr double PI = 3.14159265358979323846;

//! O VORTICE UNICO DE LeVEQUE, com reversao em `T`.
//
// E' solenoidal por construcao (vem de uma funcao de corrente) e, por causa do
// `cos(pi t / T)`, inverte na metade: em `t = T` a frente tem de voltar ao
// circulo inicial.  O que ele mede e' a deformacao REVERSIVEL -- o erro que
// sobra e' do integrador e da cirurgia, nao da fisica.
Ponto2 leveque(const Ponto2 &x, double t, double T)
{
  const double c = std::cos(PI * t / T);
  const double sx = std::sin(PI * x[0]), sy = std::sin(PI * x[1]);
  return {-sx * sx * std::sin(2 * PI * x[1]) * c,
          std::sin(2 * PI * x[0]) * sy * sy * c};
}
}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("a_curvatura_do_circulo_e_exata_em_qualquer_resolucao");
  {
    // EXATA, e e' por isso que o teste varre de 3 a 400 marcadores: se fosse
    // aproximacao de segunda ordem, o caso com 3 falharia.
    for (int n : {3, 5, 12, 97, 400}) {
      const double R = 0.37;
      Frente f = Frente::circulo({0.5, 0.5}, R, n);
      double pior = 0.0;
      for (long k = 0; k < f.num_marcadores(); k++)
        pior = std::max(pior, std::fabs(f.curvatura(k) - 1.0 / R));
      // A TOLERANCIA SEGUE O CONDICIONAMENTO, e nao e' um numero escolhido.
      //
      // `kappa = 2 cruz / (a b c)` com `cruz ~ a b sin(theta)`: quando os tres
      // marcadores ficam quase colineares (`theta = 2 pi / n` pequeno), o
      // produto vetorial cancela e o erro relativo cresce como `1/theta^2`.
      // E' propriedade da formula, nao defeito -- mas fixar `1e-12` faria o
      // caso com 400 marcadores falhar por aritmetica, e nao por erro de
      // metodo.  Medido: 9,6e-12 em n=400, contra o limite 2,4e-10 abaixo.
      const double th = 2.0 * PI / n;
      const double limite = 100.0 * 2.22e-16 * (1.0 / R) / (th * th);
      T_CHECK_MSG(pior < limite,
                  "com %d marcadores a curvatura errou %.3e (limite de "
                  "condicionamento %.3e, exata = %.6f)", n, pior, limite, 1.0 / R);
    }

    // O SINAL: positiva para a frente no sentido anti-horario, que e' como ela
    // e' montada.  Sem o sinal, uma reentrancia seria indistinguivel de uma
    // saliencia e a tensao superficial empurraria para o mesmo lado nas duas.
    Frente f = Frente::circulo({0.0, 0.0}, 1.0, 40);
    T_CHECK_MSG(f.curvatura(0) > 0.0, "a curvatura saiu negativa num circulo");
  }

  t_case("na_elipse_a_curvatura_converge_em_vez_de_acertar");
  {
    // O PAR QUE NAO E' EXATO.  Sem ele, o caso acima passa numa formula que
    // devolve uma constante qualquer.
    const double a = 0.4, b = 0.2;
    auto exata = [a, b](double t) {
      const double s = std::sin(t), c = std::cos(t);
      const double den = std::pow(a * a * s * s + b * b * c * c, 1.5);
      return a * b / den;
    };
    double erro[2];
    const int ns[2] = {64, 256};
    for (int q = 0; q < 2; q++) {
      Frente f = Frente::elipse({0.0, 0.0}, a, b, ns[q]);
      double pior = 0.0;
      for (long k = 0; k < f.num_marcadores(); k++) {
        const double t = 2.0 * PI * k / ns[q];
        pior = std::max(pior, std::fabs(f.curvatura(k) - exata(t)));
      }
      erro[q] = pior;
    }
    T_CHECK_MSG(erro[0] > 1e-6,
                "a elipse deu erro %.3e -- se fosse exata aqui tambem, o caso "
                "do circulo nao afirmaria nada", erro[0]);
    const double ordem = std::log2(erro[0] / erro[1]);
    T_CHECK_MSG(ordem > 1.7,
                "a curvatura da elipse convergiu em ordem %.2f (%.3e -> %.3e)",
                ordem, erro[0], erro[1]);
  }

  t_case("a_area_e_o_centroide_saem_da_regiao_e_nao_dos_marcadores");
  {
    // A AREA DO POLIGONO E' EXATA pelo cadarco; a do CIRCULO nao, e o erro tem
    // forma fechada -- e' ela que o teste cobra, nao uma tolerancia escolhida.
    const double R = 0.3;
    for (int n : {8, 64, 512}) {
      Frente f = Frente::circulo({0.25, -0.1}, R, n);
      const double exata = 0.5 * n * R * R * std::sin(2.0 * PI / n);
      T_NEAR(f.area(), exata, 1e-13, "a area do poligono inscrito");
      // E o centroide volta ao centro, por simetria.
      const Ponto2 c = f.centroide();
      T_NEAR(c[0], 0.25, 1e-12, "o centroide em x");
      T_NEAR(c[1], -0.1, 1e-12, "o centroide em y");
    }

    // O CENTROIDE E' DA REGIAO, e o caso que distingue e' o de marcadores NAO
    // equidistantes: ali a media dos vertices se desloca e o da regiao nao.
    std::vector<Ponto2> m;
    const int n = 64;
    for (int k = 0; k < n; k++) {
      // Densidade dobrada na metade direita: a media dos vertices puxa para la'.
      const double t = (k < n / 2) ? (PI * k / (n / 2))
                                   : (PI + PI * (k - n / 2) / (n / 2));
      const double u = (k < n / 2) ? t : PI + 0.5 * (t - PI);
      m.push_back({std::cos(u), std::sin(u)});
    }
    Frente g(m);
    double mx = 0.0;
    for (const Ponto2 &p : g.marcadores()) mx += p[0];
    mx /= n;
    const Ponto2 cg = g.centroide();
    T_CHECK_MSG(std::fabs(mx) > 1e-3,
                "a media dos vertices nao se deslocou (%.3e) -- o caso nao "
                "distingue as duas definicoes", mx);
    T_CHECK_MSG(std::fabs(cg[0]) < std::fabs(mx) / 2.0,
                "o centroide da regiao (%.3e) nao ficou mais perto do centro "
                "que a media dos vertices (%.3e)", cg[0], mx);
  }

  t_case("a_cirurgia_recusa_limites_que_a_fariam_oscilar");
  {
    Frente f = Frente::circulo({0.5, 0.5}, 0.2, 40);
    // `l_min >= l_max/2`: partir produz duas arestas abaixo do minimo, e a
    // cirurgia seguinte as funde de volta -- oscilacao a cada passo, sem nada
    // acusar.
    bool lancou = false;
    std::string msg;
    try { f.cirurgia(0.05, 0.08); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "limites que fariam a frente oscilar foram aceitos");
    T_CHECK_MSG(msg.find("oscilaria") != std::string::npos,
                "a mensagem nao diz o que aconteceria: \"%s\"", msg.c_str());

    // CONTRASTE: limites validos passam, e nao mexem numa frente ja' regular.
    Frente g = Frente::circulo({0.5, 0.5}, 0.2, 40);
    const double L = g.perimetro() / 40.0;
    const Frente::Cirurgia c = g.cirurgia(0.4 * L, 2.0 * L);
    T_CHECK_MSG(c.partidas == 0 && c.fundidas == 0,
                "a cirurgia mexeu numa frente regular: %ld partidas, %ld fusoes",
                c.partidas, c.fundidas);
    T_CHECK_MSG(g.num_marcadores() == 40, "o numero de marcadores mudou");
  }

  t_case("a_cirurgia_mantem_a_resolucao_sob_estiramento");
  {
    // SEM CIRURGIA A FRENTE MORRE em qualquer campo com estiramento.  Este
    // caso mede os dois lados: sem ela, a maior aresta explode; com ela, fica
    // limitada.
    const double T = 2.0, dt = 2e-3;
    auto campo = [T](const Ponto2 &x, double t) { return leveque(x, t, T); };

    auto maior_aresta = [](const Frente &f) {
      double L = 0.0;
      for (long k = 0; k < f.num_marcadores(); k++) {
        const Ponto2 &a = f.marcador(k);
        const Ponto2 &b = f.marcador(k + 1);
        L = std::max(L, std::hypot(b[0] - a[0], b[1] - a[1]));
      }
      return L;
    };

    Frente sem = Frente::circulo({0.5, 0.75}, 0.15, 100);
    Frente com = Frente::circulo({0.5, 0.75}, 0.15, 100);
    const double L0 = com.perimetro() / 100.0;

    // O QUE SE AFIRMA NAO E' UM CONTRASTE NUM INSTANTE, E' A TENDENCIA.
    //
    // Um contraste num instante so' depende de quanto o campo esticou ate'
    // ali, e o limiar viraria numero escolhido.  O que a cirurgia promete e'
    // que a maior aresta fica LIMITADA por `l_max` enquanto, sem ela, cresce
    // MONOTONAMENTE -- e isso nao depende de quanto o campo estica.
    //
    // Medido, em multiplos de `L0`: sem cirurgia 2,02 / 3,05 / 3,61 / 3,77 nos
    // quartos de `T/2`; com cirurgia 1,50 / 1,49 / 1,50 / 1,49, com os
    // marcadores indo de 100 a 203.
    const long npassos = static_cast<long>(T / 2 / dt);
    const long marco = npassos / 4;
    double am_sem[4], am_com[4];
    int q = 0;
    for (long k = 0; k < npassos; k++) {
      sem.avanca(campo, k * dt, dt);
      com.avanca(campo, k * dt, dt);
      com.cirurgia(0.4 * L0, 1.5 * L0);
      if ((k + 1) % marco == 0 && q < 4) {
        am_sem[q] = maior_aresta(sem) / L0;
        am_com[q] = maior_aresta(com) / L0;
        q++;
      }
    }
    T_CHECK_MSG(q == 4, "so' %d marcos medidos -- o oraculo nao percorreu a serie", q);

    // COM CIRURGIA: limitada por `l_max`, em TODOS os marcos.
    for (int i = 0; i < 4; i++)
      T_CHECK_MSG(am_com[i] <= 1.51,
                  "no marco %d a maior aresta foi %.3f L0, acima de l_max = 1.5",
                  i, am_com[i]);

    // SEM CIRURGIA: cresce monotonamente.  E' o par que torna o acima uma
    // afirmacao -- sem ele, "ficou limitada" nao distingue a cirurgia
    // funcionando de um campo que nao estica.
    for (int i = 1; i < 4; i++)
      T_CHECK_MSG(am_sem[i] > am_sem[i - 1],
                  "sem cirurgia a maior aresta NAO cresceu do marco %d (%.3f) "
                  "para o %d (%.3f) -- o campo nao estica e o caso nao afirma "
                  "nada", i - 1, am_sem[i - 1], i, am_sem[i]);
    T_CHECK_MSG(am_sem[3] > 2.0 * am_com[3],
                "sem cirurgia a maior aresta chegou a %.3f L0 contra %.3f com "
                "-- separacao insuficiente", am_sem[3], am_com[3]);
    T_CHECK_MSG(com.num_marcadores() > 100,
                "a cirurgia nao acrescentou marcador nenhum (%ld)",
                com.num_marcadores());
  }

  t_case("a_reversao_de_leveque_devolve_o_circulo");
  {
    // O ORACULO DE B1.  O campo e' solenoidal e inverte em `T/2`: em `t = T` a
    // frente tem de voltar ao circulo inicial.  O que sobra e' erro do
    // integrador e da cirurgia -- nao ha' fisica nenhuma para explicar
    // diferenca.
    const double T = 2.0, R = 0.15;
    const Ponto2 c0{0.5, 0.75};
    auto campo = [T](const Ponto2 &x, double t) { return leveque(x, t, T); };

    // O QUE A MEDIDA MOSTROU, e refutou a expectativa com que este caso
    // nasceu: o desvio NAO cai com o passo quando ha' cirurgia.
    //
    //     dt        sem cirurgia    com cirurgia
    //     8e-03      6,05e-06        9,91e-05
    //     4e-03      7,56e-07        1,01e-04
    //     2e-03      9,45e-08        1,01e-04
    //     1e-03      1,18e-08        1,01e-04
    //
    // Sem cirurgia converge com razao 8,0 a cada bisseccao -- TERCEIRA ordem,
    // e nao segunda: a reversao e' simetrica no tempo e o termo de segunda
    // ordem do ponto medio cancela entre as duas metades.
    //
    // Com cirurgia ha' um PISO independente do passo.  Nao e' defeito: partir
    // e fundir arestas e' irreversivel -- os marcadores criados no caminho de
    // ida nao desaparecem na volta --, e o piso e' o preco de manter a
    // resolucao da frente.  O que se afirma aqui e' justamente isso, e nao
    // "converge": afirmar convergencia seria afirmar o que nao acontece.
    double sem_cir[2], com_cir[2];
    const double dts[2] = {4e-3, 1e-3};
    for (int q = 0; q < 2; q++)
      for (int cir = 0; cir < 2; cir++) {
        Frente f = Frente::circulo(c0, R, 200);
        const double L0 = f.perimetro() / 200.0;
        const long npassos = static_cast<long>(std::llround(T / dts[q]));
        for (long k = 0; k < npassos; k++) {
          f.avanca(campo, k * dts[q], dts[q]);
          if (cir) f.cirurgia(0.4 * L0, 1.5 * L0);
        }
        double pior = 0.0;
        for (const Ponto2 &p : f.marcadores())
          pior = std::max(pior,
                          std::fabs(std::hypot(p[0] - c0[0], p[1] - c0[1]) - R));
        (cir ? com_cir : sem_cir)[q] = pior;
      }

    // SEM CIRURGIA o integrador manda, e CONVERGE.
    const double ordem = std::log2(sem_cir[0] / sem_cir[1]) / 2.0;
    T_CHECK_MSG(ordem > 2.5,
                "sem cirurgia a reversao convergiu em ordem %.2f (%.3e -> %.3e)"
                " -- esperada 3, pela simetria temporal da reversao",
                ordem, sem_cir[0], sem_cir[1]);

    // COM CIRURGIA ha' um PISO, e ele NAO cai com o passo.  E' afirmacao
    // positiva: o desvio tem de ficar praticamente igual.
    const double razao = com_cir[0] / com_cir[1];
    T_CHECK_MSG(razao > 0.9 && razao < 1.1,
                "com cirurgia o desvio mudou por %.3f ao quadruplicar o numero "
                "de passos (%.3e -> %.3e) -- se cai, o piso nao e' da cirurgia",
                razao, com_cir[0], com_cir[1]);

    // E O PISO E' PEQUENO: menos de 0,1% do raio.  Sem este limite, "nao cai
    // com o passo" seria compativel com a frente ter se desfeito.
    T_BELOW(com_cir[1], 1e-3 * R,
            "o piso da reversao com cirurgia, em fracao do raio");
  }

  t_case("o_rk2_bate_o_euler_na_rotacao_pura");
  {
    // ROTACAO PURA: o circulo tem de ficar circulo.  Com Euler o marcador sai
    // pela tangente e a gota INCHA -- deriva radial de `O(dt)` por passo, que
    // se acumula.  E' o caso que justifica o RK2 ser o padrao.
    auto giro = [](const Ponto2 &x, double) {
      return Ponto2{-x[1], x[0]};
    };
    const double dt = 1e-3;
    const long n = 2000;             // duas voltas inteiras
    double desvio[2];
    for (int q = 0; q < 2; q++) {
      Frente f = Frente::circulo({0.0, 0.0}, 1.0, 120);
      for (long k = 0; k < n; k++) f.avanca(giro, k * dt, dt, q == 1);
      double pior = 0.0;
      for (const Ponto2 &p : f.marcadores())
        pior = std::max(pior, std::fabs(std::hypot(p[0], p[1]) - 1.0));
      desvio[q] = pior;
    }
    T_CONTRASTE(desvio[1], desvio[0], 50.0,
                "a deriva radial do RK2 contra a do Euler em duas voltas");
    T_BELOW(desvio[1], 1e-5, "a deriva radial do RK2 em duas voltas");
  }

  return t_end();
}
