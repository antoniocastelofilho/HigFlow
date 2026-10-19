// O NAVIER-STOKES 3D: o que existe hoje e' a DIVERGENCIA, e e' o que se mede.
//
// O ORACULO E' DE EXATIDAO, nao de ordem.  A divergencia de volumes finitos e'
// a regra do ponto medio sobre a superficie, EXATA para campo linear.  Logo um
// campo linear solenoidal tem de dar ZERO DE MAQUINA, e um de divergencia
// conhecida tem de dar o valor exato -- em qualquer resolucao.  Ordem certa
// com constante errada existe; exatidao nao.
//
// E CADA AFIRMACAO DE ZERO VEM COM O PAR QUE NAO E' ZERO.  Sem ele, "a
// divergencia deu 1e-16" nao distingue um operador correto de um que devolve
// zero sempre -- e um operador que descarta faces devolve justamente isso em
// metade dos casos.
//
// O CASO QUE SEPARA AREA DE COMPRIMENTO e' o de blocos NAO CUBICOS, e ele e'
// especifico de tres dimensoes.  Numa celula cubica, `area/volume` e
// `comprimento/area` dao o MESMO numero -- trocar a medida da faceta por
// engano passaria despercebido.  Com lados diferentes, nao passa.

#include "../testing.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "campos3d.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"
#include "manufaturada-ns3d.hxx"
#include "navierstokes3d.hxx"

using namespace t8flow3d;

namespace {

//! Contorno minimo: so' para construir o objeto e medir o operador.
void arma(ContornoPorRegiao *b0, ContornoPorRegiao *b1, ContornoPorRegiao *b2,
          ContornoPorRegiao *bp, double h = 0.0)
{
  auto zero = [](const Ponto3 &) { return 0.0; };
  b0->padrao(Contorno::dirichlet(zero));
  b1->padrao(Contorno::dirichlet(zero));
  b2->padrao(Contorno::dirichlet(zero));
  // UM PONTO DE DIRICHLET NA PRESSAO: o solver recusa toda-Neumann, e com
  // razao.  A regiao usa o `h` da parede, e nao um valor digitado -- com `h`
  // errado ela pega zero facetas (e volta a recusa) ou duas (e o Poisson fica
  // sobredeterminado).
  if (h > 0.0)
    bp->registra("pino",
                 [h](const Ponto3 &x) {
                   return std::fabs(x[2]) < 1e-9 && x[0] < h && x[1] < h;
                 },
                 Contorno::dirichlet(zero));
  bp->padrao(Contorno::neumann(zero));
}

//! Escreve um campo analitico nas facetas.
void poe(const DominioFacetas &fd, CampoVetorial *u,
         const std::function<double(const Ponto3 &, int)> &f)
{
  for (int d = 0; d < 3; d++)
    for (long k = 0; k < fd.num_facetas(d); k++)
      u->em(d, k) = f(fd.centro(d, k), d);
}

//! O maximo da divergencia nas celulas que NAO encostam no contorno.
//
// SO' AS INTERIORES, e a razao e' de formulacao: na celula de borda a faceta
// externa existe e carrega o valor do campo analitico, entao o balanco ali
// mede o campo, nao o operador.  O que se verifica aqui e' o operador.
double div_interior(const Malha3 &m, const NavierStokes &ns, double h)
{
  double pior = 0.0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto3 c = m.centroide(b, i);
      bool borda = false;
      for (int d = 0; d < 3; d++)
        if (c[d] < h || c[d] > 1.0 - h) borda = true;
      if (borda) continue;
      pior = std::max(pior, std::fabs(ns.divergencia_em(b, i)));
    }
  return pior;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Os campos do arquivo, e o que cada um afirma.
  auto solenoidal = [](const Ponto3 &x, int d) {
    // div = 1 + 1 - 2 = 0
    return (d == 0) ? x[0] : (d == 1 ? x[1] : -2.0 * x[2]);
  };
  auto com_fonte = [](const Ponto3 &x, int d) {
    return (d == 0) ? 3.0 * x[0] : 0.0;        // div = 3, exatamente
  };
  auto cisalhante = [](const Ponto3 &x, int d) {
    // div = 0, mas cada componente VARIA ao longo da face que a carrega --
    // o que pega um operador que amostra no lugar errado.
    return (d == 0) ? x[1] : (d == 1 ? x[2] : x[0]);
  };

  t_case("a_divergencia_e_exata_em_malha_uniforme");
  {
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    poe(fd, &ns.velocidade(), solenoidal);
    const double z1 = div_interior(m, ns, h);
    poe(fd, &ns.velocidade(), cisalhante);
    const double z2 = div_interior(m, ns, h);

    // O PAR QUE NAO E' ZERO, e sem ele os dois acima nao afirmam nada.
    poe(fd, &ns.velocidade(), com_fonte);
    const double v3 = div_interior(m, ns, h);

    T_BELOW(z1, 1e-13, "campo solenoidal: a divergencia tem de ser zero");
    T_BELOW(z2, 1e-13, "campo cisalhante: a divergencia tem de ser zero");
    T_NEAR(v3, 3.0, 1e-12,
           "campo de divergencia 3: o operador tem de devolver 3");
  }

  t_case("a_medida_da_faceta_e_uma_AREA_e_o_caso_cubico_nao_mostra");
  {
    // BLOCOS DE LADOS DIFERENTES.  Numa celula cubica `area/volume` e
    // `comprimento/area` dao o mesmo numero, entao trocar a medida da faceta
    // por engano passaria.  Aqui o bloco e' 1 x 2 x 4 coberto por celulas
    // cubicas de lado `h`, e as tres direcoes tem contagens diferentes --
    // qualquer confusao entre as medidas separa.
    //
    // A folha do t8code e' CUBICA, entao o bloco nao-cubico vem de varios
    // blocos lado a lado, e nao de uma celula esticada.
    std::vector<Bloco> blocos;
    for (int jz = 0; jz < 4; jz++)
      for (int jy = 0; jy < 2; jy++)
        blocos.push_back({{0.0, 1.0 * jy, 1.0 * jz},
                          {1.0, 1.0 * (jy + 1), 1.0 * (jz + 1)}, 2});
    Malha3 m(blocos, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, 1.0 / 4.0);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    // Um campo com divergencia conhecida e DIFERENTE em cada direcao: se as
    // tres contribuicoes entrassem com peso errado, o total nao fecharia.
    auto tres_pesos = [](const Ponto3 &x, int d) {
      return (d == 0) ? 2.0 * x[0] : (d == 1 ? 5.0 * x[1] : -1.0 * x[2]);
    };
    poe(fd, &ns.velocidade(), tres_pesos);

    double pior = 0.0;
    long contadas = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        // A celula interior AO BLOCO basta: a faceta existe de qualquer modo,
        // e o campo e' analitico em todo o dominio.
        pior = std::max(pior, std::fabs(ns.divergencia_em(b, i) - 6.0));
        contadas++;
      }
    T_CHECK_MSG(contadas > 100, "so' %ld celulas -- o caso nao mede nada",
                contadas);
    T_BELOW(pior, 1e-12,
            "divergencia 2 + 5 - 1 = 6 em blocos de lados diferentes");
  }

  t_case("a_transicao_2_1_e_o_paralelo_sao_RECUSADOS_na_construcao");
  {
    // RECUSAR NA CONSTRUCAO, e nao no primeiro passo: quem monta uma malha
    // adaptada e so' descobre no passo mil ja' gastou a montagem, e pior,
    // pode ter visto numeros no caminho.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.0, 0.0, 0.0}, {0.4, 0.4, 0.4}, 3);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, 1.0 / 4.0);

    T_CHECK_MSG(m.num_folhas_em_transicao() > 0,
                "a malha do caso nao tem transicao 2:1 -- a recusa abaixo nao "
                "seria exercitada");

    bool lancou = false;
    std::string msg;
    try { NavierStokes ns(m, fd, b0, b1, b2, bp); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "a malha em transicao 2:1 foi aceita");
    // A MENSAGEM TEM DE DIZER O QUE E' DIFERENTE EM 3D, e nao so' "nao
    // suportado": quatro facetas por face, e duas direcoes transversais por
    // faceta em vez de uma.
    //
    // ESTE CASO JA' EXIGIU "OITO", porque a mensagem dizia que a media da
    // viscosidade tocava oito celulas.  Nao toca: o ponto fica sobre uma
    // ARESTA, tocada por quatro, como a quina em 2D.  O caso passava e
    // sustentava a afirmacao errada -- um teste so' verifica o que lhe pedem.
    T_CHECK_MSG(msg.find("QUATRO") != std::string::npos &&
                    msg.find("DUAS direcoes") != std::string::npos,
                "a recusa nao diz por que 3D nao e' 2D com um indice a mais: "
                "\"%s\"", msg.c_str());

    // CONTRASTE: sem transicao, constroi.
    Malha3 m2({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    DominioFacetas fd2(m2);
    bool ok = true;
    try { NavierStokes ns2(m2, fd2, b0, b1, b2, bp); }
    catch (const std::exception &) { ok = false; }
    T_CHECK_MSG(ok, "a malha uniforme foi recusada");
  }

  // =========================================================================
  t_case("a_projecao_zera_a_divergencia_de_um_campo_qualquer");
  {
    // O ORACULO NAO PRECISA DE FISICA.  Dado um campo de velocidade QUALQUER,
    // projeta-lo tem de deixar divergencia zero de maquina -- e isso exercita
    // a cadeia inteira: divergencia, montagem do Poisson, solucao, gradiente.
    //
    // E' O PAR `div`/`grad` QUE FECHA, nao cada um.  Se a projecao usasse um
    // gradiente diferente do que a divergencia do Poisson supos, o campo
    // corrigido nao teria divergencia nula e NADA DIRIA: o solver convergiria
    // igual, e o defeito so' apareceria na fisica, muito depois.
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);
    ns.passo_de_tempo(1e-3);

    // CAMPO COM DIVERGENCIA DE VERDADE, e nulo nas facetas de contorno -- ali
    // a velocidade vem da condicao e a projecao nao a toca.
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) {
        const Ponto3 x = fd.centro(d, f);
        bool borda = (std::fabs(x[d]) < 1e-12 || std::fabs(x[d] - 1.0) < 1e-12);
        ns.velocidade().em(d, f) =
            borda ? 0.0
                  : 0.3 + 1.7 * x[0] - 0.9 * x[1] + 2.1 * x[0] * x[2]
                        + 0.6 * x[1] * x[1];
      }

    const double antes = ns.divergencia_maxima();
    ns.projeta_campo();
    const double depois = ns.divergencia_maxima();

    // O CONTRASTE QUE DA' SENTIDO AO ZERO: sem ele, "a divergencia deu 1e-14"
    // nao distingue a projecao de um campo que ja' era solenoidal.
    T_CHECK_MSG(antes > 1.0,
                "o campo de partida ja' tinha divergencia quase nula (%.3e) "
                "-- a projecao nao teria o que fazer e o caso seria vazio",
                antes);
    T_BELOW(depois, 1e-11,
            "depois da projecao a divergencia tem de ser zero de maquina");
    T_CHECK_MSG(ns.pressao_convergiu(), "o Poisson nao convergiu");
  }

  t_case("a_divergencia_e_o_gradiente_sao_adjuntos_em_tres_dimensoes");
  {
    // AFIRMA SOBRE O PAR, e nao sobre cada operador -- ver a clausula O1 no
    // lado 2D.  Em tres dimensoes o peso em que a adjuncao vale e'
    // `AREA x distancia`, e nao `comprimento x distancia`: usar o peso de 2D
    // aqui nao fecharia nem em malha uniforme.
    const int nv = 3;
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, 1.0 / (1 << nv));
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    poe(fd, &ns.velocidade(), [](const Ponto3 &x, int d) {
      return (d == 0) ? 0.7 + 1.3 * x[0] - 0.4 * x[1] + 2.0 * x[0] * x[2]
             : (d == 1) ? -0.2 + 0.9 * x[0] * x[0] + 1.1 * x[1]
                        : 0.5 - 0.8 * x[2] + 1.4 * x[1] * x[2];
    });
    // Zera o contorno: ali a velocidade nao e' incognita.
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++)
        if (ns.volume_dual_da_faceta(d, f) == 0.0)
          ns.velocidade().em(d, f) = 0.0;

    CampoEscalar q(m);
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        q.em(m.id_global(b, i)) =
            0.3 + 1.7 * c[0] - 0.8 * c[1] + 2.3 * c[0] * c[1] + 0.6 * c[2] * c[2];
      }

    double lado_d = 0.0, lado_g = 0.0, escala = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 hh = m.lados_da_folha(b, i);
        const double v = hh[0] * hh[1] * hh[2];
        const double t = ns.divergencia_em(b, i) * q.em(m.id_global(b, i)) * v;
        lado_d += t;
        escala += std::fabs(t);
      }
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) {
        const double w = ns.volume_dual_da_faceta(d, f);
        if (w == 0.0) continue;
        const double t = ns.velocidade().em(d, f) * ns.gradiente_em(q, d, f) * w;
        lado_g += t;
        escala += std::fabs(t);
      }

    T_CHECK_MSG(escala > 1e-3,
                "as formas bilineares deram quase zero (%.3e) -- os campos "
                "nao exercitam nada", escala);
    const double rel = std::fabs(lado_d + lado_g) / escala;
    T_BELOW(rel, 1e-13,
            "D e G nao sao adjuntos em malha uniforme tridimensional");
  }

  // =========================================================================
  t_case("com_mu_uniforme_o_termo_viscoso_e_o_laplaciano_EXATO");
  {
    // O ORACULO E' DE EXATIDAO.  A diferenca de tres pontos e' exata para
    // QUADRATICA, logo com `mu` uniforme o termo viscoso tem de dar
    // `mu * lap(u)` em precisao de maquina -- em qualquer resolucao.
    //
    // Isto verifica a estrutura de fluxo nas TRES direcoes, e e' o caso que
    // pegaria uma direcao esquecida: em duas dimensoes o laco tinha dois
    // termos, aqui tem tres.
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);
    const double RE = 2.5;
    ns.reynolds(RE);

    // Quadratica com as TRES segundas derivadas diferentes de zero e
    // DIFERENTES ENTRE SI: `lap = 2a + 2b + 2c`.  Com coeficientes iguais,
    // somar duas direcoes em vez de tres daria um numero proporcional ao
    // certo, e o erro passaria por constante.
    const double A = 1.7, B = -0.9, C = 2.3;
    poe(fd, &ns.velocidade(), [&](const Ponto3 &x, int d) {
      if (d != 0) return 0.0;
      return A * x[0] * x[0] + B * x[1] * x[1] + C * x[2] * x[2];
    });
    const double lap = 2.0 * (A + B + C);
    const double esperado = lap / RE;

    double pior = 0.0;
    long contadas = 0;
    for (long f = 0; f < fd.num_facetas(0); f++) {
      const Ponto3 x = fd.centro(0, f);
      bool borda = false;
      for (int d = 0; d < 3; d++)
        if (x[d] < h || x[d] > 1.0 - h) borda = true;
      if (borda) continue;
      pior = std::max(pior, std::fabs(ns.tensao_viscosa(0, f) - esperado));
      contadas++;
    }
    T_CHECK_MSG(contadas > 50, "so' %ld facetas interiores", contadas);
    T_BELOW(pior, 1e-11,
            "com mu uniforme o termo tem de dar mu*lap(u) exatamente");

    // E O TRANSPOSTO E' IDENTICAMENTE NULO com `mu` uniforme -- e' o que faz
    // o caso monofasico reproduzir o laplaciano bit a bit.
    double pior_t = 0.0;
    for (long f = 0; f < fd.num_facetas(0); f++)
      pior_t = std::max(pior_t, std::fabs(ns.tensao_viscosa_transposto(0, f)));
    T_BELOW(pior_t, 0.0,
            "sem campo de viscosidade o transposto tem de ser EXATAMENTE zero");
  }

  t_case("mu_na_aresta_e_a_media_harmonica_de_QUATRO_celulas");
  {
    // QUATRO, E NAO OITO.  O ponto onde `mu` e' avaliado fica sobre dois
    // planos de face com a terceira coordenada num centro de celula -- isso
    // e' uma ARESTA, tocada por quatro celulas.  Oito seria um VERTICE, que
    // este termo nao usa.  A mensagem de recusa e a clausula do contrato ja'
    // disseram oito, e estavam erradas.
    //
    // E HARMONICA, NAO ARITMETICA: o caso exige que as duas DIFIRAM, senao
    // nao distinguiria uma da outra.
    const int nv = 2;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    // `mu` com contraste forte, para harmonica e aritmetica se separarem.
    CampoEscalar mu(m);
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        mu.em(m.id_global(b, i)) = (c[0] < 0.5) ? 1.0 : 10.0;
      }
    ns.viscosidade(&mu);

    // Uma aresta no interior, sobre o plano `x = 0,5` e um plano de face em y.
    const Ponto3 aresta{0.5, 0.5, 0.5 + h / 2.0};
    // As quatro celulas que a tocam: duas com mu=1 e duas com mu=10.
    const double harm = 4.0 / (1.0 / 1.0 + 1.0 / 1.0 + 1.0 / 10.0 + 1.0 / 10.0);
    const double arit = (1.0 + 1.0 + 10.0 + 10.0) / 4.0;
    const double obtido = ns.mu_na_aresta(aresta, 0, 1);

    T_NEAR(obtido, harm, 1e-12,
           "mu na aresta tem de ser a media HARMONICA das quatro celulas");
    T_CHECK_MSG(std::fabs(harm - arit) > 1.0,
                "harmonica (%.4f) e aritmetica (%.4f) nao se distinguem neste "
                "caso -- ele nao verifica qual das duas e'", harm, arit);
  }

  t_case("o_transposto_existe_e_so_com_mu_variavel");
  {
    // O TERMO QUE EU NAO TINHA NO LADO 2D, e cuja falta custou circularidade
    // 0,8361 contra 0,9012 na bolha de Hysing.  Aqui so' se verifica que ele
    // existe, e' nulo com `mu` uniforme e NAO nulo com `mu` variavel.
    //
    // O CONTRASTE E' O CASO: um transposto que devolvesse zero sempre
    // passaria na metade de cima.
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    // Campo com derivadas CRUZADAS nao nulas: `d_x u_y` e `d_x u_z`.
    poe(fd, &ns.velocidade(), [](const Ponto3 &x, int d) {
      return (d == 0) ? 0.3 * x[1] : (d == 1 ? 1.9 * x[0] : -1.1 * x[0]);
    });

    CampoEscalar mu(m);
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        mu.em(m.id_global(b, i)) = 1.0 + 2.0 * c[0] + 3.0 * c[1] + 0.5 * c[2];
      }
    ns.viscosidade(&mu);

    double pior = 0.0;
    for (long f = 0; f < fd.num_facetas(0); f++) {
      const Ponto3 x = fd.centro(0, f);
      bool borda = false;
      for (int d = 0; d < 3; d++)
        if (x[d] < 2 * h || x[d] > 1.0 - 2 * h) borda = true;
      if (borda) continue;
      pior = std::max(pior, std::fabs(ns.tensao_viscosa_transposto(0, f)));
    }
    T_CHECK_MSG(pior > 0.1,
                "com mu variavel e derivadas cruzadas nao nulas o transposto "
                "deu %.3e -- um termo que devolve zero sempre passaria no "
                "caso de mu uniforme", pior);
  }

  // =========================================================================
  t_case("o_convectivo_e_EXATO_para_campo_uniforme_e_linear");
  {
    // DOIS ORACULOS DE FORMA FECHADA, e o segundo distingue os esquemas.
    //
    //   campo UNIFORME `u = (a,b,c)`:  div(uu) = 0 em QUALQUER esquema.
    //   campo `u_x = A x`:             div(uu) = 2 A^2 x.
    //
    // O PRIMEIRO passa com o termo inteiro trocado por zero -- ele verifica a
    // estrutura, nao o valor.  O segundo nao passa.
    //
    // E O SEGUNDO NAO VALE PARA O UPWIND, e a primeira versao deste caso
    // afirmava que valia para os quatro.  O upwind toma `phi` na FACETA a
    // montante (`A x`), e nao no CENTRO da celula (`A(x +- h/2)`) onde fica a
    // face do volume de controle.  Fazendo a conta:
    //
    //     V+ phi+ - V- phi-  =  A(x+h/2) A x - A(x-h/2) A x  =  A^2 x h
    //
    // dividido por `h` da' `A^2 x` -- EXATAMENTE METADE do certo.  Isso e'
    // oraculo melhor do que "nao e' exato": o upwind tem de dar a metade, e
    // o caso exige isso.  Central, CUBISTA e QUICK sao exatos (em NVD o campo
    // linear da' `fc = 1/2`, e os tres devolvem `A(x + h/2)`).
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    const Convectivo esquemas[] = {Convectivo::Central, Convectivo::Upwind1,
                                   Convectivo::Cubista, Convectivo::Quick};
    const char *nomes[] = {"central", "upwind1", "cubista", "quick"};

    for (int k = 0; k < 4; k++) {
      ns.convectivo(esquemas[k]);
      const bool exato = (esquemas[k] != Convectivo::Upwind1);

      // --- campo uniforme: div(uu) = 0 --------------------------------
      poe(fd, &ns.velocidade(), [](const Ponto3 &, int d) {
        return (d == 0) ? 1.3 : (d == 1 ? -0.7 : 2.1);
      });
      double pior0 = 0.0;
      for (int d = 0; d < 3; d++)
        for (long f = 0; f < fd.num_facetas(d); f++) {
          const Ponto3 x = fd.centro(d, f);
          bool borda = false;
          for (int q = 0; q < 3; q++)
            if (x[q] < 2 * h || x[q] > 1.0 - 2 * h) borda = true;
          if (borda) continue;
          pior0 = std::max(pior0, std::fabs(ns.conveccao(d, f)));
        }
      T_BELOW(pior0, 1e-12,
              "campo uniforme: div(uu) tem de ser zero em qualquer esquema");

      // --- `u_x = A x`: div(uu) = 2 A^2 x ------------------------------
      const double A = 1.4;
      poe(fd, &ns.velocidade(),
          [&](const Ponto3 &x, int d) { return (d == 0) ? A * x[0] : 0.0; });
      ns.zera_rebaixamentos();
      double pior1 = 0.0;
      long contadas = 0;
      for (long f = 0; f < fd.num_facetas(0); f++) {
        const Ponto3 x = fd.centro(0, f);
        bool borda = false;
        for (int q = 0; q < 3; q++)
          if (x[q] < 2 * h || x[q] > 1.0 - 2 * h) borda = true;
        if (borda) continue;
        // O upwind tem de dar METADE, e os outros o valor cheio.
        const double esperado = (exato ? 2.0 : 1.0) * A * A * x[0];
        pior1 = std::max(pior1, std::fabs(ns.conveccao(0, f) - esperado));
        contadas++;
      }
      T_CHECK_MSG(contadas > 20, "%s: so' %ld facetas interiores", nomes[k],
                  contadas);
      T_BELOW(pior1, 1e-11,
              exato ? "campo linear: div(uu) tem de dar 2 A^2 x"
                    : "campo linear: o upwind tem de dar EXATAMENTE A^2 x, "
                      "metade do certo");
    }
  }

  t_case("os_esquemas_DISCRIMINAM_num_campo_de_grau_maior");
  {
    // SE OS QUATRO DAO SEMPRE O MESMO, o caso acima nao esta' verificando
    // esquema nenhum -- esta' verificando a estrutura, e qualquer
    // interpolacao serviria.  Este caso existe para mostrar que eles
    // DIFEREM onde tem de diferir: num campo que nenhum deles reproduz.
    //
    // E' a licao do canal de Poiseuille no lado 2D, onde o convectivo se
    // anulava identicamente e os tres esquemas davam o mesmo numero -- de
    // modo que "o CUBISTA esta' correto" nao era afirmacao verificada.
    const int nv = 3;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);

    auto mede = [&](Convectivo c) {
      ns.convectivo(c);
      // Campo CUBICO e sem simetria, nas tres componentes: grau acima do que
      // qualquer um dos esquemas reproduz.
      poe(fd, &ns.velocidade(), [](const Ponto3 &x, int d) {
        const double t = x[0] * x[0] * x[1] + 0.7 * x[2] * x[2] * x[2]
                         - 0.3 * x[0] * x[1] * x[2];
        return (d == 0) ? t : (d == 1 ? 0.6 * t - 0.2 : -0.4 * t + 0.1);
      });
      double s = 0.0;
      for (long f = 0; f < fd.num_facetas(0); f++) {
        const Ponto3 x = fd.centro(0, f);
        bool borda = false;
        for (int q = 0; q < 3; q++)
          if (x[q] < 2 * h || x[q] > 1.0 - 2 * h) borda = true;
        if (borda) continue;
        s += std::fabs(ns.conveccao(0, f));
      }
      return s;
    };

    const double up = mede(Convectivo::Upwind1);
    const double ce = mede(Convectivo::Central);
    const double cu = mede(Convectivo::Cubista);
    const double qu = mede(Convectivo::Quick);

    T_CHECK_MSG(up > 1e-6, "o campo nao gera conveccao (%.3e)", up);
    T_CHECK_MSG(std::fabs(up - ce) / up > 1e-6,
                "upwind e central deram o mesmo (%.6g, %.6g) -- o caso nao "
                "discrimina esquema", up, ce);
    T_CHECK_MSG(std::fabs(cu - ce) / up > 1e-6,
                "CUBISTA e central deram o mesmo (%.6g, %.6g)", cu, ce);
    T_CHECK_MSG(std::fabs(cu - qu) / up > 1e-9,
                "CUBISTA e QUICK deram o mesmo (%.6g, %.6g) -- eles so' "
                "coincidem na faixa central do NVD, e o campo tem de sair "
                "dela em algum lugar", cu, qu);
  }

  t_case("o_rebaixamento_do_convectivo_e_CONTADO");
  {
    // UM ESQUEMA QUE SE ANUNCIA DE SEGUNDA ORDEM e cai para primeira em parte
    // do dominio entrega menos do que diz.  Perto da parede a queda e'
    // inevitavel -- o ponto a montante cai fora --, mas tem de ser VISIVEL.
    //
    // IMPORTA MAIS EM TRES DIMENSOES: sao tres direcoes por faceta em vez de
    // duas, logo mais chances de o montante faltar.  O caso exige que a conta
    // seja nao nula varrendo o dominio inteiro, e ZERO varrendo so' o miolo.
    //
    // NIVEL 4 E NAO 3: com margem de `3h` o miolo do nivel 3 tem 36 facetas,
    // e uma amostra desse tamanho nao sustenta "zero rebaixamentos no miolo".
    const int nv = 4;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);
    ns.convectivo(Convectivo::Cubista);
    poe(fd, &ns.velocidade(), [](const Ponto3 &x, int d) {
      return (d == 0) ? 1.0 + x[1] : (d == 1 ? 0.5 * x[0] : 0.2 * x[2]);
    });

    ns.zera_rebaixamentos();
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) ns.conveccao(d, f);
    const long com_borda = ns.rebaixamentos_convectivos();

    ns.zera_rebaixamentos();
    long no_miolo = 0;
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) {
        const Ponto3 x = fd.centro(d, f);
        bool borda = false;
        for (int q = 0; q < 3; q++)
          if (x[q] < 3 * h || x[q] > 1.0 - 3 * h) borda = true;
        if (!borda) { ns.conveccao(d, f); no_miolo++; }
      }
    const long rebaixos_miolo = ns.rebaixamentos_convectivos();

    T_CHECK_MSG(com_borda > 0,
                "nenhum rebaixamento no dominio inteiro -- junto a' parede o "
                "ponto a montante TEM de faltar, logo o contador nao esta' "
                "contando");
    T_CHECK_MSG(no_miolo > 50, "so' %ld facetas no miolo", no_miolo);
    T_CHECK_MSG(rebaixos_miolo == 0,
                "%ld rebaixamentos NO MIOLO do dominio -- ali o ponto a "
                "montante existe, e cair para primeira ordem seria defeito e "
                "nao geometria", rebaixos_miolo);
  }

  t_case("as_derivadas_a_mao_batem_com_diferencas_finitas");
  {
    // A GUARDA CONTRA ESCORREGAO DE DERIVACAO.  O termo fonte da solucao
    // manufaturada e' montado de DOZE formulas escritas a mao -- nove
    // `du_i/dx_j`, tres `dp/dx_i` --, e uma delas errada produz um `f` errado,
    // que produz uma solucao numerica que converge BONITO para a coisa errada.
    //
    // ISTO NAO E' CALCULAR `f` NUMERICAMENTE.  O codigo usa a formula a mao; o
    // caso compara a formula com uma diferenca central do proprio `u`.  Se as
    // duas concordam ate' o erro de truncamento da diferenca, a formula esta'
    // certa.  Calcular `f` por diferencas tornaria o oraculo uma comparacao
    // entre dois erros meus.
    ManufaturadaNS3 ms;
    ms.re = 1.7;
    const double t = 0.3, eps = 1e-5;
    // O erro da diferenca central e' `O(eps^2)` vezes a terceira derivada, que
    // aqui e' `O(pi^3)`: cerca de 3e-9.  A tolerancia tem de ser maior que
    // isso e muito menor que a escala das derivadas, que e' `O(pi^2) ~ 10`.
    const double tol = 1e-6;

    double pior_u = 0.0, pior_p = 0.0, escala = 0.0;
    const Ponto3 pontos[] = {{0.23, 0.41, 0.67}, {0.80, 0.15, 0.37},
                             {0.51, 0.74, 0.09}, {0.33, 0.33, 0.33}};
    for (const Ponto3 &x : pontos)
      for (int j = 0; j < 3; j++) {
        Ponto3 xp = x, xm = x;
        xp[j] += eps;
        xm[j] -= eps;
        for (int i = 0; i < 3; i++) {
          const double fd = (ms.u(xp, i, t) - ms.u(xm, i, t)) / (2.0 * eps);
          pior_u = std::max(pior_u, std::fabs(fd - ms.dudx(x, i, j, t)));
          escala = std::max(escala, std::fabs(fd));
        }
        const double fdp = (ms.p(xp) - ms.p(xm)) / (2.0 * eps);
        pior_p = std::max(pior_p, std::fabs(fdp - ms.dpdx(x, j)));
      }
    T_CHECK_MSG(escala > 1.0, "as derivadas sao quase nulas (%.3e) nos pontos "
                "escolhidos -- a comparacao nao mede nada", escala);
    T_BELOW(pior_u, tol, "alguma das nove du_i/dx_j nao bate com a diferenca");
    T_BELOW(pior_p, tol, "alguma das tres dp/dx_i nao bate com a diferenca");

    // E O LAPLACIANO, que vale `-3 pi^2 u` por construcao.
    double pior_l = 0.0;
    for (const Ponto3 &x : pontos)
      for (int i = 0; i < 3; i++) {
        double soma = 0.0;
        for (int j = 0; j < 3; j++) {
          Ponto3 xp = x, xm = x;
          xp[j] += eps;
          xm[j] -= eps;
          soma += (ms.u(xp, i, t) - 2.0 * ms.u(x, i, t) + ms.u(xm, i, t)) /
                  (eps * eps);
        }
        pior_l = std::max(pior_l, std::fabs(soma - ms.lap(x, i, t)));
      }
    // A diferenca de segunda ordem tem erro `O(eps^-2)` de arredondamento com
    // `eps` pequeno; a tolerancia reflete isso e nao a formula.
    T_BELOW(pior_l, 1e-2, "lap(u) nao bate com a diferenca de segunda ordem");

    // E A SOLUCAO E' SOLENOIDAL -- por construcao, mas conferido.
    double pior_d = 0.0;
    for (const Ponto3 &x : pontos) {
      double div = 0.0;
      for (int j = 0; j < 3; j++) div += ms.dudx(x, j, j, t);
      pior_d = std::max(pior_d, std::fabs(div));
    }
    T_BELOW(pior_d, 1e-12, "a manufaturada nao e' solenoidal");
  }

  t_case("cada_operador_reproduz_o_termo_da_manufaturada");
  {
    // A CAMADA QUE FALTAVA entre os casos unitarios e o MMS completo.
    //
    // Os operadores tem casos proprios com oraculo de exatidao, e o MMS mede a
    // ordem do passo inteiro.  Entre os dois nao havia nada -- e quando a
    // ordem do MMS saiu errada, nao havia como dizer QUAL operador discorda da
    // manufaturada sem instrumentar o exemplo.
    //
    // Aqui cada um e' comparado com o termo analitico correspondente, no mesmo
    // campo.  O oraculo nao e' exatidao: o campo e' trigonometrico e nenhum
    // esquema o reproduz.  E' CONSISTENCIA -- o erro tem de ser da ordem de
    // `h^2` vezes a escala do termo, e nao da ordem do termo.
    ManufaturadaNS3 ms;
    ms.re = 10.0;
    const int nv = 4;
    const double h = 1.0 / (1 << nv);
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, b2, bp;
    arma(&b0, &b1, &b2, &bp, h);
    NavierStokes ns(m, fd, b0, b1, b2, bp);
    ns.reynolds(ms.re);
    ns.convectivo(Convectivo::Central);

    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++)
        ns.velocidade().em(d, f) = ms.u(fd.centro(d, f), d, 0.0);

    double pior_v = 0.0, esc_v = 0.0, pior_c = 0.0, esc_c = 0.0;
    long contadas = 0;
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) {
        const Ponto3 x = fd.centro(d, f);
        bool borda = false;
        for (int q = 0; q < 3; q++)
          if (x[q] < 2 * h || x[q] > 1.0 - 2 * h) borda = true;
        if (borda) continue;
        contadas++;

        const double visc_exato = ms.lap(x, d, 0.0) / ms.re;
        pior_v = std::max(pior_v,
                          std::fabs(ns.tensao_viscosa(d, f) - visc_exato));
        esc_v = std::max(esc_v, std::fabs(visc_exato));

        double conv_exato = 0.0;
        for (int j = 0; j < 3; j++)
          conv_exato += ms.u(x, j, 0.0) * ms.dudx(x, d, j, 0.0);
        pior_c = std::max(pior_c, std::fabs(ns.conveccao(d, f) - conv_exato));
        esc_c = std::max(esc_c, std::fabs(conv_exato));
      }

    T_CHECK_MSG(contadas > 100, "so' %ld facetas interiores", contadas);
    T_CHECK_MSG(esc_v > 1.0 && esc_c > 1.0,
                "os termos sao quase nulos (%.3g, %.3g) -- nada medido",
                esc_v, esc_c);
    // `h = 1/16`, logo `h^2 = 4e-3`.  Erro de segunda ordem da' alguns por
    // cento da escala; erro da ORDEM da escala e' inconsistencia.
    T_CHECK_MSG(pior_v < 0.10 * esc_v,
                "o termo VISCOSO discorda da manufaturada em %.3g, escala "
                "%.3g (%.1f%%) -- com h = 1/16 isso nao e' discretizacao",
                pior_v, esc_v, 100.0 * pior_v / esc_v);
    // O LIMIAR DO CONVECTIVO E' 5% DEPOIS DA CORRECAO DA INTERPOLACAO.  Com
    // `u_no_ponto` interpolando so' ao longo da propria direcao ele dava
    // 19,6%, e 15% teria sido limiar frouxo o bastante para aceitar o
    // defeito -- o que decidiu nao foi o limiar, foi medir a ORDEM.
    T_CHECK_MSG(pior_c < 0.05 * esc_c,
                "o termo CONVECTIVO discorda da manufaturada em %.3g, escala "
                "%.3g (%.1f%%)", pior_c, esc_c, 100.0 * pior_c / esc_c);

    // E O GRADIENTE DE PRESSAO.
    CampoEscalar pp(m);
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        pp.em(m.id_global(b, i)) = ms.p(m.centroide(b, i));
    double pior_g = 0.0, esc_g = 0.0;
    for (int d = 0; d < 3; d++)
      for (long f = 0; f < fd.num_facetas(d); f++) {
        const Ponto3 x = fd.centro(d, f);
        bool borda = false;
        for (int q = 0; q < 3; q++)
          if (x[q] < 2 * h || x[q] > 1.0 - 2 * h) borda = true;
        if (borda) continue;
        const double g = ms.dpdx(x, d);
        pior_g = std::max(pior_g, std::fabs(ns.gradiente_em(pp, d, f) - g));
        esc_g = std::max(esc_g, std::fabs(g));
      }
    T_CHECK_MSG(pior_g < 0.05 * esc_g,
                "o GRADIENTE discorda da manufaturada em %.3g, escala %.3g",
                pior_g, esc_g);
  }

  return t_end();
}
