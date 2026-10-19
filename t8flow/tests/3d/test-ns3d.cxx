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
#include "navierstokes3d.hxx"

using namespace t8flow3d;

namespace {

//! Contorno minimo: so' para construir o objeto e medir o operador.
void arma(ContornoPorRegiao *b0, ContornoPorRegiao *b1, ContornoPorRegiao *b2,
          ContornoPorRegiao *bp)
{
  auto zero = [](const Ponto3 &) { return 0.0; };
  b0->padrao(Contorno::dirichlet(zero));
  b1->padrao(Contorno::dirichlet(zero));
  b2->padrao(Contorno::dirichlet(zero));
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
    arma(&b0, &b1, &b2, &bp);
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
    arma(&b0, &b1, &b2, &bp);
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
    arma(&b0, &b1, &b2, &bp);

    T_CHECK_MSG(m.num_folhas_em_transicao() > 0,
                "a malha do caso nao tem transicao 2:1 -- a recusa abaixo nao "
                "seria exercitada");

    bool lancou = false;
    std::string msg;
    try { NavierStokes ns(m, fd, b0, b1, b2, bp); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "a malha em transicao 2:1 foi aceita");
    // A MENSAGEM TEM DE DIZER O QUE E' DIFERENTE EM 3D, e nao so' "nao
    // suportado": quatro facetas por face e oito celulas na quina.
    T_CHECK_MSG(msg.find("QUATRO") != std::string::npos &&
                    msg.find("OITO") != std::string::npos,
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

  return t_end();
}
