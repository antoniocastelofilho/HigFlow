// O estencil em tres dimensoes: suporte por aneis, fechamento de contorno.
//
// A PERGUNTA QUE ESTE ARREIO RESPONDE, e que a medida dos aneis deixou aberta:
// o suporte que o crescimento por faces entrega em 3D SUSTENTA o ajuste?
//
// Medido antes: a ordem 3 pede 20 monomios, e tres aneis entregam 20 pontos no
// canto -- matriz normal quadrada, sem redundancia nenhuma.  Suporte sem folga
// e' suporte que rebaixa, e o contador de rebaixamento e' quem diz.

#include "../testing.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "contorno3d.hxx"
#include "estencil3d.hxx"
#include "malha3d.hxx"

using namespace t8flow3d;

namespace {

double grau1(const Ponto3 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 1.5 * p[2]; }

double grau2(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return grau1(p) + 0.5 * x * x - 0.6 * y * y + 0.9 * z * z
         + 0.25 * x * y - 0.4 * x * z + 0.7 * y * z;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("o_suporte_cresce_por_aneis_e_alcanca_o_minimo_da_ordem");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);

    // ORDEM 2 (10 monomios) e ORDEM 3 (20).  O teto de aneis e' generoso, e o
    // que o caso mede e' se o suporte ALCANCA o minimo -- nao se ele para cedo.
    for (int ord : {2, 3}) {
      ConstrutorEstencil ce(m, ord, 80, 6);
      ce.dirichlet([](const Ponto3 &) { return 0.0; });

      long menor = 1 << 30, maior = 0, soma = 0, n = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          const ConstrutorEstencil::Suporte s = ce.suporte_em(m.centroide(b, i));
          const long t = static_cast<long>(s.pontos.size());
          menor = std::min(menor, t);
          maior = std::max(maior, t);
          soma += t;
          n++;
        }
      T_CHECK_MSG(menor >= Mls3::min_pontos(ord),
                  "ordem %d: o menor suporte tem %ld pontos, abaixo do minimo %d",
                  ord, menor, Mls3::min_pontos(ord));
      printf("  (ordem %d: suporte de %ld a %ld pontos, media %.1f; minimo %d)\n",
             ord, menor, maior, double(soma) / n, Mls3::min_pontos(ord));
    }
  }

  t_case("a_caixa_sustenta_a_ordem_tres_onde_a_cruz_nao_sustenta");
  {
    // O CASO QUE DECIDIU A ESTRATEGIA DE SUPORTE, e o par de medidas e' o caso
    // inteiro: as duas na MESMA malha, mesma ordem, mesmo teto de pontos.
    //
    // O crescimento por aneis de face entrega um suporte em CRUZ.  Os pontos
    // bastam em NUMERO -- media de 24,4 para 20 monomios -- e nao em
    // DISTRIBUICAO: quase tudo alinhado aos eixos, e os monomios mistos `xy`,
    // `xz`, `yz` e `xyz` ficam mal determinados.  Em duas dimensoes isso nao
    // aparece, porque ha' menos termos mistos.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);

    long reb[2] = {0, 0}, consultas = 0;
    const Suporte3 estrategias[2] = {Suporte3::Aneis, Suporte3::Caixa};
    for (int e = 0; e < 2; e++) {
      ConstrutorEstencil ce(m, 3, 80, 6);
      ce.estrategia(estrategias[e]);
      ce.dirichlet([&](const Ponto3 &p) { return grau2(p); });
      ce.mls().zera_contadores();

      consultas = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          // ALVO DESLOCADO DO CENTRO: no centro o ajuste e' degenerado -- peso 1
          // na propria celula -- e o fechamento nao seria exercitado.  Em duas
          // dimensoes isso fez um caso inteiro medir nada.
          const Ponto3 c = m.centroide(b, i);
          const Ponto3 h = m.lados_da_folha(b, i);
          ce.derivada_em({c[0] + 0.3 * h[0], c[1], c[2]}, 0);
          consultas++;
        }
      reb[e] = ce.mls().rebaixamentos();
    }

    printf("  (ordem 3, %ld consultas: aneis rebaixam %ld (%.1f%%), "
           "caixa rebaixa %ld (%.1f%%))\n", consultas, reb[0],
           100.0 * reb[0] / consultas, reb[1], 100.0 * reb[1] / consultas);

    T_CHECK_MSG(reb[1] == 0,
                "a CAIXA rebaixou %ld de %ld ajustes (%.1f%%) -- ela foi escolhida "
                "exatamente por nao rebaixar", reb[1], consultas,
                100.0 * reb[1] / consultas);

    // O CONTRASTE E' A OUTRA METADE: sem ele, "a caixa nao rebaixa" nao
    // distingue uma boa distribuicao de um contador quebrado.
    T_CHECK_MSG(reb[0] > consultas / 2,
                "os ANEIS rebaixaram so' %ld de %ld -- se a cruz sustentasse a "
                "ordem 3, a caixa nao teria razao de existir e este caso nao "
                "estaria medindo o que diz", reb[0], consultas);
  }

  t_case("o_estencil_reproduz_o_campo_do_grau_da_ordem");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    ConstrutorEstencil ce(m, 2, 80, 6);
    ce.dirichlet([&](const Ponto3 &p) { return grau2(p); });

    // O ORACULO COMPLETO: montar o estencil, aplicar ao campo amostrado nos
    // centros, e comparar com o valor exato no alvo.  Isso exercita o suporte,
    // os pesos E o fechamento de contorno de uma vez.
    double pior = 0.0;
    Ponto3 onde{-1, -1, -1};
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        const Ponto3 h = m.lados_da_folha(b, i);
        const Ponto3 alvo{c[0] + 0.3 * h[0], c[1] - 0.2 * h[1], c[2] + 0.1 * h[2]};

        const Estencil e = ce.em(alvo);
        double v = e.rhs;
        for (size_t k = 0; k < e.ids.size(); k++) {
          const long x = m.indice_do_id(e.ids[k]);
          if (x < 0) continue;
          v += e.pesos[k] * grau2(m.centroide_do_indice(x));
        }
        const double err = std::fabs(v - grau2(alvo));
        if (err > pior) { pior = err; onde = alvo; }
      }
    T_CHECK_MSG(pior < 1e-9,
                "pior erro %.3e em (%.3f, %.3f, %.3f) ao reproduzir um campo de "
                "grau 2 com ajuste de ordem 2", pior, onde[0], onde[1], onde[2]);
    printf("  (reproducao de grau 2: pior erro %.3e)\n", pior);
  }

  t_case("o_fechamento_de_contorno_entra_nas_seis_faces");
  {
    // EM TRES DIMENSOES SAO SEIS PAREDES, e o erro mais facil de cometer ao
    // portar de 2D e' esquecer `z-` e `z+`.  O oraculo e' geometrico: todo
    // ponto de parede que um estencil registra tem de cair na fronteira do
    // CUBO, e as seis faces tem de aparecer.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    ConstrutorEstencil ce(m, 2, 80, 6);
    ce.dirichlet([](const Ponto3 &) { return 1.0; });

    long paredes = 0, no_interior = 0;
    bool vista[6] = {false, false, false, false, false, false};
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const ConstrutorEstencil::Suporte s = ce.suporte_em(m.centroide(b, i));
        for (size_t k = 0; k < s.ids.size(); k++) {
          if (s.ids[k] >= 0) continue;        // id negativo marca a parede
          paredes++;
          const Ponto3 &x = s.pontos[k];
          bool na_fronteira = false;
          for (int d = 0; d < 3; d++) {
            if (std::fabs(x[d]) < 1e-12) { na_fronteira = true; vista[2 * d] = true; }
            if (std::fabs(x[d] - 1.0) < 1e-12) { na_fronteira = true; vista[2 * d + 1] = true; }
          }
          if (!na_fronteira) no_interior++;
        }
      }

    T_CHECK_MSG(paredes > 0, "nenhuma parede registrada -- o fechamento nao foi "
                "exercitado e o caso nao mediu nada");
    T_CHECK_MSG(no_interior == 0,
                "%ld pontos de parede no INTERIOR do cubo", no_interior);
    for (int f = 0; f < 6; f++)
      T_CHECK_MSG(vista[f],
                  "a face %d (%s%c) nunca apareceu -- em tres dimensoes o erro "
                  "mais facil e' esquecer z- e z+", f,
                  (f % 2) ? "+" : "-", "xxyyzz"[f]);
    printf("  (%ld paredes, as seis faces vistas)\n", paredes);
  }

  t_case("o_suporte_de_neumann_usa_fantasma_e_o_de_dirichlet_amostra_a_parede");
  {
    // OS DOIS TRATAMENTOS SAO DIFERENTES, e a diferenca e' deliberada:
    //
    //   DIRICHLET amostra a PAREDE -- o valor ali e' conhecido e nao depende de
    //   incognita nenhuma.  Refletir (A=-1) amplificaria o erro do interior;
    //   NEUMANN usa FANTASMA -- um ponto do lado de fora cujo valor e' afim na
    //   incognita da celula.
    //
    // O contador de cada um esta' no suporte, e o caso afirma que o tipo da
    // condicao decide qual aparece.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);

    ConstrutorEstencil cd(m, 2, 80, 6);
    cd.dirichlet([](const Ponto3 &) { return 1.0; });
    const ConstrutorEstencil::Suporte sd = cd.suporte_em(m.centroide(0, 0));

    ContornoPorRegiao bn;
    bn.padrao(Contorno::neumann([](const Ponto3 &) { return 0.5; }));
    ConstrutorEstencil cn(m, 2, 80, 6);
    cn.contorno(bn);
    const ConstrutorEstencil::Suporte sn = cn.suporte_em(m.centroide(0, 0));

    T_CHECK_MSG(sd.paredes > 0 && sd.fantasmas == 0,
                "Dirichlet devia dar paredes e nenhum fantasma; deu %ld e %ld",
                sd.paredes, sd.fantasmas);
    T_CHECK_MSG(sn.fantasmas > 0 && sn.paredes == 0,
                "Neumann devia dar fantasmas e nenhuma parede; deu %ld e %ld",
                sn.fantasmas, sn.paredes);
    printf("  (Dirichlet: %ld paredes; Neumann: %ld fantasmas)\n",
           sd.paredes, sn.fantasmas);
  }

  t_case("ponto_fora_do_dominio_e_recusado");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    ConstrutorEstencil ce(m, 2, 80, 6);
    ce.dirichlet([](const Ponto3 &) { return 0.0; });

    bool lancou = false;
    try { ce.em({5.0, 0.5, 0.5}); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "um ponto fora do dominio recebeu estencil");
  }

  t_case("o_ponto_fora_da_franja_e_contado_e_o_fora_do_dominio_nao");
  {
    // AS DUAS RESPOSTAS NEGATIVAS DE `id_do_ponto` SAO DIFERENTES, e confundi-las
    // esconde a que importa.  Fora do DOMINIO e' legitimo -- o contorno trata.
    // Dentro do dominio e nao achado significa que falta FRANJA, e um suporte
    // truncado da' pesos diferentes dos do serial.
    //
    // ESTE CASO EXISTE PORQUE A PRIMEIRA VERSAO DESCARTAVA OS DOIS EM SILENCIO,
    // e a medida de "quantas camadas bastam" dava ZERO truncamento com qualquer
    // numero -- um resultado bom demais, que e' sintoma e nao conquista.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    ConstrutorEstencil ce(m, 2, 200, 6);
    ce.dirichlet([](const Ponto3 &) { return 0.0; });
    ce.zera_contadores_de_franja();

    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        ce.suporte_em(m.centroide(b, i));

    // EM SERIAL NAO HA' FRANJA, logo nenhum ponto pode ser truncado por ela: os
    // pontos da caixa que nao existem sao os de FORA do dominio, e esses o
    // contorno trata.  Qualquer truncamento aqui seria o contador confundindo as
    // duas situacoes.
    T_CHECK_MSG(ce.truncados_pela_franja() == 0,
                "%ld suportes truncados em SERIAL, onde nao ha' franja -- o "
                "contador esta' marcando os pontos fora do DOMINIO",
                ce.truncados_pela_franja());

    // E a malha sabe distinguir, que e' de onde a distincao vem.
    T_CHECK_MSG(m.dentro_do_dominio({0.5, 0.5, 0.5}), "o centro nao esta' no dominio?");
    T_CHECK_MSG(!m.dentro_do_dominio({5.0, 0.5, 0.5}), "um ponto distante entrou");
    T_CHECK_MSG(!m.dentro_do_dominio({0.5, 0.5, -0.1}),
                "um ponto logo abaixo de `z=0` entrou -- em tres dimensoes a face "
                "que se esquece ao portar de 2D e' justamente essa");
  }

  return t_end();
}
