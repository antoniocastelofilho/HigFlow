// A BOLHA ASCENDENTE DE HYSING: o degrau B3, e o caso de comparacao.
//
// E' o unico caso deste projeto cujo oraculo NAO e' de forma fechada.  Ele e'
// um BENCHMARK PUBLICADO -- Hysing et al. --, e tem uma propriedade que o
// torna utilizavel mesmo assim: ele publica series de TRES codigos
// independentes, de modo que a discordancia ENTRE as referencias esta'
// quantificada.
//
// ESSA DISCORDANCIA E' O LIMITE DE RESOLUCAO DO ORACULO.  Pedir a este sistema
// que faca melhor que ela seria pedir mais precisao do que a referencia possui
// -- e um resultado "melhor que a referencia" aqui nao seria uma conquista,
// seria sinal de que alguma coisa foi ajustada para dar aquele numero.
//
//     grandeza                 TP2D     FreeLIFE   MooNMD    consenso
//     circularidade minima    0,9013    0,9011     0,9013    0,9012 +- 0,0001
//     instante do minimo      1,9041    1,8750     1,9000    ~1,90
//     velocidade maxima       0,2417    0,2421     0,2417    0,2417 - 0,2421
//     instante do maximo      0,9213    0,9313     0,9239    ~0,92
//     centro de massa em t=3  1,0813    1,0799     1,0817    1,0799 - 1,0817
//
// AS TRES GRANDEZAS SAIEM DA FRACAO VOLUMETRICA, e isso nao e' detalhe de
// implementacao: e' o que torna a comparacao com o VOF possivel.  O VOF
// advecta esse campo; o front-tracking o escreve da geometria da frente.  A
// medicao e' a MESMA nos dois, e qualquer diferenca que sobreviva pertence ao
// metodo e nao a' forma de medir.
//
// A CIRCULARIDADE E' `2 sqrt(pi A) / P`, e nao `4 pi A / P^2`.  As duas sao
// chamadas de circularidade na literatura, e uma e' o QUADRADO da outra.
// Comparar sem converter da' 0,7987 contra um publicado de 0,9013 -- quando o
// numero comparavel e' `sqrt(0,7987) = 0,8937`.  A discrepancia parece erro de
// metodo e e' erro de definicao.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <memory>
#include <vector>

#include "acoplamento2d.hxx"
#include "arranque.hxx"
#include "campos2d.hxx"
#include "criterio2d.hxx"
#include "contorno2d.hxx"
#include "frente2d.hxx"
#include "malha2d.hxx"
#include "navierstokes2d.hxx"
#include "saida_vtk2d.hxx"

using namespace t8flow2d;

namespace {
constexpr double PI = 3.14159265358979323846;

struct Medidas {
  double area = 0.0, circ = 0.0, yc = 0.0, vc = 0.0;
};

//! AS TRES GRANDEZAS, todas da fracao volumetrica (menos o perimetro).
Medidas mede(const Malha &m, const DominioFacetas &fd, const CampoEscalar &c,
             const CampoVetorial &u, const Frente &f)
{
  Medidas r;
  double soma = 0.0, somay = 0.0, somav = 0.0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const long id = m.id_global(b, i);
      const Ponto2 cc = m.centroide(b, i);
      const Ponto2 hh = m.lados_da_folha(b, i);
      const double w = c.em(id) * hh[0] * hh[1];
      if (w <= 0.0) continue;

      // A VELOCIDADE VERTICAL NO CENTRO, media das duas facetas `y` da celula.
      const long fm = fd.id_faceta(b, i, 1, 0), fp = fd.id_faceta(b, i, 1, 1);
      double v = 0.0;
      int n = 0;
      if (fm >= 0) { v += u.em(1, fm); n++; }
      if (fp >= 0) { v += u.em(1, fp); n++; }
      if (n > 0) v /= n;

      soma += w;
      somay += cc[1] * w;
      somav += v * w;
    }
  r.area = soma;
  if (soma > 0.0) { r.yc = somay / soma; r.vc = somav / soma; }

  // O PERIMETRO vem da GEOMETRIA da frente.  E' o unico que nao sai da
  // fracao, e o VOF o tiraria do comprimento das retas reconstruidas -- e'
  // onde as duas medicoes nao podem ser identicas, e isso fica dito.
  const double P = f.perimetro();
  if (P > 0.0) r.circ = 2.0 * std::sqrt(PI * r.area) / P;
  return r;
}
}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0, np = 1;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &np);

  // OS PARAMETROS DO CASO 1 DE HYSING, na forma DIMENSIONAL.  Ela e'
  // inequivoca: `Re = 35` e `Eo = 10` saem deles, e nao o contrario.
  double rho_f = 1000.0, rho_b = 100.0, mu_f = 10.0, mu_b = 1.0;
  double sigma = 24.5, g = 0.98, R = 0.25, T = 3.0, dt = 0.0;
  int nivel = 5, tsur_cada = 20, base = 1, amr = 0, amr_cada = 20,
      amr_celmin = 5;
  std::string nome_tempo = "semi_implicit_euler";
  double tsur_rlx = 0.5, tsur_lim = 0.05, ds_sobre_h = 1.0;
  double amr_largura = -1.0;   // negativo: deriva de R mais abaixo
  long cada_vtk = 0, cada_med = 10;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--base=", 7)) base = std::atoi(a + 7);
    else if (!std::strncmp(a, "--ds-sobre-h=", 13)) ds_sobre_h = std::atof(a + 13);
    else if (!std::strncmp(a, "--tempo=", 8)) nome_tempo = a + 8;
    else if (!std::strncmp(a, "--amr=", 6)) amr = std::atoi(a + 6);
    else if (!std::strncmp(a, "--amr-cada=", 11)) amr_cada = std::atoi(a + 11);
    else if (!std::strncmp(a, "--amr-celmin=", 13)) amr_celmin = std::atoi(a + 13);
    else if (!std::strncmp(a, "--amr-largura=", 14)) amr_largura = std::atof(a + 14);
    else if (!std::strncmp(a, "--T=", 4)) T = std::atof(a + 4);
    else if (!std::strncmp(a, "--dt=", 5)) dt = std::atof(a + 5);
    else if (!std::strncmp(a, "--vtk=", 6)) cada_vtk = std::atol(a + 6);
    else if (!std::strncmp(a, "--medir=", 8)) cada_med = std::atol(a + 8);
    else if (!std::strncmp(a, "--tsur-cada=", 12)) tsur_cada = std::atoi(a + 12);
    else if (!std::strncmp(a, "--tsur-rlx=", 11)) tsur_rlx = std::atof(a + 11);
    else if (!std::strncmp(a, "--tsur-limiar=", 14)) tsur_lim = std::atof(a + 14);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf("uso: ./exemplos/hysing/bolha [opcoes]\n\n"
                    "  --nivel=N       refino dentro de cada bloco\n"
                    "  --base=N        base de N x 2N blocos de lado 1/N;\n"
                    "                  h = 1/(N*2^nivel).  --base=5 --nivel=5\n"
                    "                  da' h = 0,00625, o h_min do HiGFlow\n"
                    "  --ds-sobre-h=R  espacamento de marcador em unidades de h;\n"
                    "                  1,0 e' o FT_DS_SOBRE_H do HiGFlow\n"
                    "  --tempo=NOME    explicit_euler ou semi_implicit_euler\n"
                    "                  (padrao: o do HiGFlow neste caso,\n"
                    "                  semi_implicit_euler)\n"
                    "  --amr=N         niveis de refino em volta da frente\n"
                    "                  (0 desliga).  --base=5 --nivel=3\n"
                    "                  --amr=2 e' a malha do HiGFlow:\n"
                    "                  base 40x80 com h_min = 0,00625\n"
                    "  --amr-cada=N    remalha a cada N passos (FT_AMR_CADA)\n"
                    "  --amr-celmin=N  celulas finas minimas de cada lado\n"
                    "                  da frente (FT_AMR_CELMIN)\n"
                    "  --amr-largura=L largura FISICA da banda fina.  O padrao\n"
                    "                  e' 0,2 R; com 0 a banda volta a ser so'\n"
                    "                  contagem de celulas, que ANTICONVERGE\n"
                    "  --T, --dt       instante final e passo\n"
                    "  --vtk=N         grava VTK a cada N passos\n"
                    "  --medir=N       mede as grandezas a cada N passos\n"
                    "  --tsur-cada=N   aplica o TSUR a cada N passos (0 desliga)\n"
                    "  --tsur-rlx=R    relaxacao; 0,5 e' a conta em 2D\n"
                    "  --tsur-limiar=L inclinacao minima, em fracao da meia-aresta\n");
      return 0;
    }
  }
  if (np > 1) {
    if (rank == 0) std::printf("ESTE CASO AINDA E' SERIAL (%d pedidos).\n", np);
    return 2;
  }

  // O DOMINIO E' 1 x 2, COBERTO POR BLOCOS QUADRADOS.  A folha do t8code e'
  // quadrada: um bloco 1x2 com refino uniforme daria celula 1:2, e o caso
  // deixaria de ser o de Hysing.
  //
  // POR QUE A BASE E' PARAMETRO, e nao fixa em 1x2 blocos.  O refino dentro do
  // bloco dobra, logo `h = 1/2^nivel` so' alcanca potencias de dois -- e o
  // `h_min` do HiGFlow neste caso e' 0,00625 = 1/160, que NAO e' uma delas
  // (160 = 2^5 * 5).  Com uma base de 5 x 10 blocos e nivel 5 sai exatamente
  // 1/(5*32) = 0,00625.  Sem isso a comparacao com o HiGFlow teria de aceitar
  // malha diferente, e a diferenca de resultado ficaria sem dono.
  if (base < 1) {
    if (rank == 0) std::printf("--base tem de ser pelo menos 1 (%d)\n", base);
    return 2;
  }
  const double lb = 1.0 / base;
  std::vector<Bloco> blocos;
  for (int jy = 0; jy < 2 * base; jy++)
    for (int ix = 0; ix < base; ix++)
      blocos.push_back({{ix * lb, jy * lb}, {(ix + 1) * lb, (jy + 1) * lb},
                        nivel});
  // DUAS ESCALAS, e confundi-las ja' me custou uma corrida: `h_base` e' a
  // celula da base (onde a tabela de AMR decide) e `h` e' a MENOR do dominio
  // (de onde sai o passo e a contagem de marcadores).  Tirar o `dt` da base
  // num dominio refinado aperta o passo por `4^amr` tarde demais: diverge
  // calada, com NaN centenas de passos depois.
  const double h_base = lb / (1 << nivel);
  const double h = h_base / (1 << amr);
  const int nx_base = base * (1 << nivel);   // 1 / h_base
  const int ny_base = 2 * nx_base;           // o dominio e' 1 x 2

  // O PASSO, PELO LIMITE MAIS APERTADO DOS QUE VALEM, e impresso para que a
  // escolha fique a' vista em vez de embutida:
  //   viscoso    dt < rho h^2 / (4 mu)   -- o menor `rho/mu` dos dois fluidos
  //   capilar    dt < sqrt((rho1+rho2) h^3 / (4 pi sigma))
  //   convectivo dt < h / |u|
  //
  // O VISCOSO SO' VALE NO EXPLICITO.  Com o termo `div(mu grad u)` na matriz
  // ele deixa de limitar, e manter o limite seria pagar por implicitacao sem
  // receber nada.
  //
  // O QUE SOBRA E' O CAPILAR, e isso nao se remove implicitando a difusao: ele
  // e' a onda capilar na interface, propriedade da tensao superficial.  Neste
  // caso, em `h = 0,00625`, viscoso 9,77e-4 contra capilar 9,34e-4 -- o capilar
  // ja' mandava, e o semi-implicito NAO acelera nada aqui.  Ele entra por
  // CORRESPONDENCIA com o HiGFlow, que usa `semi_implicit_euler`, e nao por
  // velocidade.  Em `h` mais fino a ordem se inverte (viscoso cai com `h^2`, o
  // capilar com `h^{3/2}`) e ai' o ganho aparece.
  const Tempo tempo = tempo_do_nome(nome_tempo);
  const double dt_visc = std::min(rho_b / mu_b, rho_f / mu_f) * h * h / 4.0;
  const double dt_cap = std::sqrt((rho_f + rho_b) * h * h * h / (4.0 * PI * sigma));
  const double dt_conv = h / 0.3;
  const bool visc_limita = (tempo == Tempo::EulerExplicito);
  if (dt <= 0.0)
    dt = 0.4 * (visc_limita ? std::min({dt_visc, dt_cap, dt_conv})
                            : std::min(dt_cap, dt_conv));
  const long passos = static_cast<long>(std::llround(T / dt));

  // A BOLHA: fracao UM dentro dela.  Quatro marcadores por celula cortada.
  const Ponto2 centro{0.5, 0.5};
  // A CONTAGEM E' PAR, E ISSO NAO E' ARREDONDAMENTO: os marcadores ficam em
  // `theta_k = 2 pi k / n`, e a reflexao do caso em `x = 0,5` leva
  // `theta -> pi - theta`, ou seja `k -> n/2 - k`.  Com `n` IMPAR esse indice
  // nao existe e o conjunto inicial NAO e' simetrico -- semente de assimetria
  // em t = 0, muito acima do arredondamento, que o escoamento depois amplifica.
  // Medido com n = 201 (o que o nivel 5 dava): o centro de area saiu de 0 e
  // chegou a +2,0e-2 em t = 3, 8% do raio, crescendo de forma monotona.
  int nmarc = std::max(32, static_cast<int>(2.0 * PI * R / (ds_sobre_h * h)));
  if (nmarc % 2) nmarc++;
  Frente frente = Frente::circulo(centro, R, nmarc);
  const double ds0 = frente.perimetro() / nmarc;

  // ========================================================================
  // A MONTAGEM, E O REMALHAMENTO COMO "RESET"
  // ========================================================================
  //
  // Tudo o que depende da MALHA vive junto num bloco que se joga fora e se
  // refaz: malha, facetas, campos, contornos e solver.  E' o que o HiGFlow faz
  // em `higflow_reconstroi_dominio`, e o que foi pedido: quando a malha muda,
  // recalculam-se o solver, os estencis e os identificadores.
  //
  // POR QUE REFAZER E NAO ADAPTAR NO LUGAR.  `Malha::refina_por_criterio` so'
  // sabe PEDIR MAIS nivel -- engrossar exige decidir o que fazer com o campo
  // dos filhos, e isso e' do remalhamento.  Como a frente SOBE, a banda fina
  // que ela deixa para tras tem de voltar a ser grossa; sem engrossar, em
  // poucas dezenas de remalhamentos a malha inteira estaria fina e o AMR nao
  // teria servido para nada.  Refazer da base resolve os dois sentidos com um
  // caminho so'.
  //
  // A ORDEM DE DECLARACAO IMPORTA: `ns` guarda PONTEIROS para os contornos e
  // referencias para malha e facetas, e e' destruido ANTES deles porque a
  // destruicao segue a ordem inversa da declaracao.  Trocar a ordem aqui da'
  // uso depois de liberado, que nao falha na hora.
  struct Montagem {
    std::unique_ptr<Malha> m;
    std::unique_ptr<DominioFacetas> fd;
    std::unique_ptr<ContornoPorRegiao> bcu0, bcu1, bcp;
    std::unique_ptr<CampoEscalar> frac, rho, mu;
    std::unique_ptr<CampoVetorial> forca;
    std::unique_ptr<NavierStokes> ns;
  };

  const double E = 1e-9;
  auto zero = [](const Ponto2 &) { return 0.0; };
  auto em_cima_ou_baixo = [E](const Ponto2 &x) {
    return std::fabs(x[1]) < E || std::fabs(x[1] - 2.0) < E;
  };
  auto nos_lados = [E](const Ponto2 &x) {
    return std::fabs(x[0]) < E || std::fabs(x[0] - 1.0) < E;
  };

  // OS LIMIARES SAO DERIVADOS, nao digitados: a regra das celulas minimas E a
  // largura fisica da banda, o maior dos dois.
  //
  // O PADRAO SAI DO RAIO DA BOLHA, e e' aqui que essa conta pertence -- o
  // criterio nao sabe de que forma vem a semente, mas o exemplo sabe o seu `R`.
  //
  // E O VALOR E' MEDIDO, NAO ESCOLHIDO.  Varredura em `h_min = 0,0125`, contra
  // o uniforme na mesma `h_min` (`v_max` 0,2431) e a referencia
  // (0,2417-0,2421):
  //
  //     L        banda     celulas   v_max     y_c(3)
  //     0,05     0,50 R      3704    0,2262    1,0646
  //     0,125    0,75 R      4886    0,2395    1,0694
  //     0,25     1,25 R      6545    0,2429    1,0751
  //     uniforme   --       12800    0,2431    1,0764
  //
  // Convergencia MONOTONA para o uniforme com a banda alargando.  Em `L = R` o
  // adaptativo o reproduz a 0,08% com METADE das celulas.
  //
  // EU HAVIA POSTO `0,2 R` PRIMEIRO, por corresponder a's 8 celulas finas que o
  // HiGFlow resolveu em `h_min = 0,00625`.  A medida mostrou que e' estreito
  // demais: da' -6,5% em `v_max`.  Contar celulas da corrida dele nao substitui
  // medir a largura que o ESCOAMENTO pede, e as duas nao coincidiram.
  if (amr_largura < 0.0) amr_largura = 1.0 * R;
  const std::vector<double> thr =
      (amr > 0) ? amr_limiares(h_base, amr, amr_celmin, amr_largura)
                : std::vector<double>{};

  auto monta = [&](const std::vector<signed char> *tab) {
    auto M = std::make_unique<Montagem>();
    M->m = std::make_unique<Malha>(blocos, sc_MPI_COMM_WORLD);
    if (tab != nullptr) {
      const Ponto2 lo{0.0, 0.0};
      M->m->refina_por_criterio([&](const Ponto2 &x) {
        return nivel + amr_nivel_em(*tab, lo, h_base, h_base, nx_base, ny_base,
                                    x);
      });
    }
    M->fd = std::make_unique<DominioFacetas>(*M->m);

    // CONTORNO DE HYSING: nao escorregamento em cima e embaixo, DESLIZAMENTO
    // LIVRE nos lados.  Deslizamento livre e' normal nula com tangencial de
    // gradiente nulo -- tratar os lados como parede mudaria o caso.
    M->bcu0 = std::make_unique<ContornoPorRegiao>();
    M->bcu1 = std::make_unique<ContornoPorRegiao>();
    M->bcp = std::make_unique<ContornoPorRegiao>();
    M->bcu0->registra("topo_fundo", em_cima_ou_baixo, Contorno::dirichlet(zero));
    M->bcu0->registra("lados", nos_lados, Contorno::dirichlet(zero));
    M->bcu1->registra("topo_fundo", em_cima_ou_baixo, Contorno::dirichlet(zero));
    M->bcu1->registra("lados", nos_lados, Contorno::neumann(zero));

    // O PINO DA PRESSAO SAI DA MALHA, e nao de um `h` digitado.  A regiao tem
    // de conter EXATAMENTE um centro de faceta de parede: com dois, dois pontos
    // ficam presos e o Poisson resolve um problema sobredeterminado; com zero,
    // a pressao volta a ser toda-Neumann e o solver recusa.  Com AMR nao ha'
    // um `h` de parede conhecido de antemao -- o refino pode chegar ali --,
    // entao procura-se a celula de parede mais a' esquerda e crava-se nela.
    double px = 0.0, ph = h_base, py_min = 1e30;
    for (long b = 0; b < M->m->num_blocos(); b++)
      for (long i = 0; i < M->m->num_folhas_do_bloco(b); i++) {
        const Ponto2 c = M->m->centroide(b, i);
        const Ponto2 hh = M->m->lados_da_folha(b, i);
        if (c[1] - hh[1] / 2.0 > E) continue;          // nao encosta no fundo
        if (c[0] < py_min) { py_min = c[0]; px = c[0]; ph = hh[0]; }
      }
    M->bcp->registra("pino",
                     [E, px, ph](const Ponto2 &x) {
                       return std::fabs(x[1]) < E &&
                              std::fabs(x[0] - px) < ph / 4.0;
                     },
                     Contorno::dirichlet(zero));
    M->bcp->padrao(Contorno::neumann(zero));

    M->ns = std::make_unique<NavierStokes>(*M->m, *M->fd, *M->bcu0, *M->bcu1,
                                           *M->bcp);
    M->ns->tempo(tempo);
    M->ns->passo_de_tempo(dt);
    M->ns->gravidade({0.0, -g});

    M->frac = std::make_unique<CampoEscalar>(*M->m);
    M->rho = std::make_unique<CampoEscalar>(*M->m);
    M->mu = std::make_unique<CampoEscalar>(*M->m);
    M->forca = std::make_unique<CampoVetorial>(*M->fd);
    M->ns->densidade(M->rho.get());
    M->ns->viscosidade(M->mu.get());
    M->ns->forca_de_faceta(M->forca.get());
    return M;
  };

  // A TABELA DE NIVEIS, DA FRENTE.  As sementes sao os MARCADORES -- o
  // analogo direto da celula cortada no VOF -- e a dilatacao faz o resto.
  auto tabela_da_frente = [&](const Frente &f) {
    return amr_tabela(f.marcadores(), {0.0, 0.0}, h_base, h_base, nx_base,
                      ny_base, thr);
  };

  // A TRANSFERENCIA.  So' dois campos atravessam: `p` e `u`.
  //
  // `frac`, `rho`, `mu` e a forca NAO atravessam -- sao reescritos da GEOMETRIA
  // DA FRENTE a cada passo, e interpola-los seria introduzir erro onde ha'
  // fonte exata.  Essa e' a vantagem do front-tracking aqui: a interface nao
  // vive na malha, entao remalhar nao a difunde.
  //
  // A VELOCIDADE TRANSFERIDA NAO E' DISCRETAMENTE SOLENOIDAL, e isso e'
  // esperado: interpolar nao preserva a divergencia discreta.  Nao precisa de
  // passe de projecao proprio -- o Poisson do passo SEGUINTE resolve com
  // `div(u*)` no lado direito, e a sujeira entra ali e sai.  O medidor `div`
  // acusa se nao sair.
  auto transfere = [&](const Montagem &velho, Montagem *novo) {
    for (long b = 0; b < novo->m->num_blocos(); b++)
      for (long i = 0; i < novo->m->num_folhas_do_bloco(b); i++) {
        const long id = novo->m->id_global(b, i);
        novo->ns->pressao().em(id) =
            velho.ns->pressao().em_ponto(novo->m->centroide(b, i));
      }
    for (int d = 0; d < 2; d++)
      for (long f = 0; f < novo->fd->num_facetas(d); f++) {
        bool ok = false;
        const double v =
            velho.ns->velocidade_em(d, novo->fd->centro(d, f), &ok);
        novo->ns->velocidade().em(d, f) = ok ? v : 0.0;
      }
    novo->ns->instante(velho.ns->instante());
  };

  // ---- O QUE A TRANSFERENCIA FAZ COM A ENERGIA -------------------------
  //
  // MEDIDA, E NAO INFERIDA DA TENDENCIA.  O par de controle mostrou que o
  // adaptativo e o uniforme discordam 7% em `v_max` com a MESMA `h_min` e o
  // mesmo `dt`, e que o erro cresce com o NUMERO de remalhamentos (136 ->
  // -6,5%, 391 -> -10,5%) em vez de cair com `h`.  Isso aponta a
  // transferencia, mas apontar nao e' medir.
  //
  // O PESO E' O VOLUME DE CONTROLE DUAL da faceta, que e' zero no contorno --
  // ali a velocidade nao e' incognita e nao deve entrar na soma.  `rho` fica
  // DE FORA de proposito: o que se quer e' o salto RELATIVO da mesma soma
  // antes e depois, e a frente nao se move entre os dois, entao incluir `rho`
  // mudaria a escala e nao o salto.
  auto energia = [](const Montagem &M) {
    double e = 0.0;
    for (int d = 0; d < 2; d++)
      for (long f = 0; f < M.fd->num_facetas(d); f++) {
        const double w = M.ns->volume_dual_da_faceta(d, f);
        if (w == 0.0) continue;
        const double u = M.ns->velocidade().em(d, f);
        e += u * u * w;
      }
    return e;
  };
  double perda_total = 0.0, perda_pior = 0.0;

  std::vector<signed char> tab;
  if (amr > 0) tab = tabela_da_frente(frente);
  std::unique_ptr<Montagem> M = monta(amr > 0 ? &tab : nullptr);
  long remalhas = 0;

  if (rank == 0) {
    std::printf("caso        bolha ascendente de Hysing, caso 1\n");
    std::printf("dominio     [0,1]x[0,2], %d x %d blocos, %ld celulas, "
                "h = %.5g\n", base, 2 * base, M->m->num_folhas_globais(), h);
    std::printf("bolha       R = %.3g em (%.2f, %.2f), %d marcadores "
                "(ds/h = %.3g)\n", R, centro[0], centro[1], nmarc, ds_sobre_h);
    std::printf("fluidos     rho %.0f/%.0f   mu %.3g/%.3g   sigma %.4g   g %.4g\n",
                rho_f, rho_b, mu_f, mu_b, sigma, g);
    std::printf("            Re = %.3g   Eo = %.3g\n",
                rho_f * std::sqrt(g) * std::pow(2 * R, 1.5) / mu_f,
                rho_f * g * (2 * R) * (2 * R) / sigma);
    std::printf("dt          %.5g  (visc %.3g, capilar %.3g, convect %.3g)\n",
                dt, dt_visc, dt_cap, dt_conv);
    std::printf("T           %.4g em %ld passos\n", T, passos);
    std::printf("tempo       %s%s\n", nome_tempo.c_str(),
                visc_limita ? "" : "   (o limite viscoso nao se aplica)");
    if (amr > 0)
      std::printf("AMR         %d niveis, a cada %d passos, banda fina %.4g "
                  "(= %.2f R), piso %d celulas\n",
                  amr, amr_cada, thr.empty() ? 0.0 : thr.back(),
                  (thr.empty() ? 0.0 : thr.back()) / R, amr_celmin);
    std::printf("TSUR        a cada %d passos, rlx %.3g, limiar %.3g\n\n",
                tsur_cada, tsur_rlx, tsur_lim);
    // OS MEDIDORES DO ESCOAMENTO vao na MESMA linha das grandezas, e nao num
    // relatorio ao final: uma corrida de milhares de passos que degenera no
    // meio tem de ser vista ENQUANTO degenera.  Cada um tem valor esperado, e
    // sair dele e' sinal antes de ser erro:
    //
    //   CFL      acima de 1 e' instabilidade convectiva a caminho;
    //   div      divergencia maxima; se sobe, a projecao parou de fechar;
    //   its      iteracoes do Poisson; se disparam, o sistema esta' mal
    //            condicionado -- tipicamente a razao de densidade;
    //   marc     marcadores; crescimento continuo e' frente esticando;
    //   area/A0  a bolha NAO muda de area.  Este e' o unico destes que e'
    //            erro, e nao sintoma;
    //   x-0,5    QUEBRA DE SIMETRIA.  O caso e' simetrico em x = 0,5 -- dado
    //            inicial, contorno e gravidade, todos.  Qualquer valor aqui e'
    //            numerico, e se CRESCE de forma monotona ha' semente
    //            estrutural sendo amplificada.  Ficou visivel porque uma
    //            animacao mostrou a bolha torta e o medidor nao dizia nada.
    std::printf("%10s %12s %12s %12s %10s %6s %8s %8s %7s %5s %9s\n",
                "t", "circular.", "y_centro", "v_ascensao", "area/A0", "marc",
                "|u|max", "CFL", "div", "its", "x-0,5");
    std::fflush(stdout);
  }

  // ---- a saida de visualizacao -------------------------------------------
  //
  // TODAS AS PROPRIEDADES, e a interface como poligonal em arquivo proprio.
  // Sem a fracao, a densidade e a viscosidade no mesmo quadro nao da' para
  // ver se a interface que a frente desenha e' a mesma que o campo ve' -- e
  // essa e' justamente a conferencia que uma imagem resolve na hora.
  // O ESCRITOR NASCE A CADA QUADRO, e nao uma vez: ele guarda REFERENCIA para
  // a malha, e depois de um remalhamento a malha de antes nao existe mais.
  // Guarda-lo fora do laco daria uso depois de liberado -- que escreve arquivo
  // com numero, e nao falha.
  auto grava = [&](int quadro) {
    t8flow::SaidaVtk2d vtk(*M->m);
    vtk.limpa();
    vtk.vetorial("velocidade", (*M->fd), M->ns->velocidade());
    vtk.escalar("pressao", M->ns->pressao());
    vtk.escalar("fracao", (*M->frac));
    vtk.escalar("densidade", (*M->rho));
    vtk.escalar("viscosidade", (*M->mu));
    vtk.vetorial("forca_tensao", (*M->fd), (*M->forca));
    vtk.fronteira("interface", {frente.marcadores()});
    vtk.escreve("hysing", quadro);
  };

  escreve_fracao(frente, (*M->m), M->frac.get());
  interpola_propriedade((*M->m), (*M->frac), rho_b, rho_f, M->rho.get());
  interpola_propriedade((*M->m), (*M->frac), mu_b, mu_f, M->mu.get());
  const double A0 = mede((*M->m), (*M->fd), (*M->frac), M->ns->velocidade(), frente).area;
  int quadro = 0;
  if (cada_vtk > 0) grava(quadro++);

  long tot_part = 0, tot_fund = 0, tot_tsur = 0;
  double circ_min = 1e300, t_circ = 0.0, v_max = -1e300, t_vmax = 0.0;
  Medidas fim;

  for (long k = 0; k < passos; k++) {
    // A ORDEM DO PASSO, e ela nao e' livre: a fracao vem da GEOMETRIA da
    // frente, e as propriedades vem da fracao.  Atualizar a frente antes de
    // escrever a fracao deixaria o campo um passo adiantado em relacao a'
    // forca que o moveu.
    escreve_fracao(frente, (*M->m), M->frac.get());
    interpola_propriedade((*M->m), (*M->frac), rho_b, rho_f, M->rho.get());
    interpola_propriedade((*M->m), (*M->frac), mu_b, mu_f, M->mu.get());
    forca_tensao(frente, (*M->m), (*M->fd), sigma, M->forca.get());

    M->ns->passo();
    if (!M->ns->pressao_convergiu()) {
      if (rank == 0)
        std::printf("PARANDO no passo %ld: o Poisson nao convergiu\n", k + 1);
      return 1;
    }

    // A FRENTE SEGUE O CAMPO.  A interpolacao e' a ADJUNTA da distribuicao da
    // forca -- mesmo nucleo --, que e' o que impede a frente de ganhar ou
    // perder energia na travessia.
    frente.move(velocidade_nos_marcadores(frente, (*M->m), (*M->fd), M->ns->velocidade()), dt);
    const Frente::Cirurgia cir = frente.cirurgia(0.4 * ds0, 1.5 * ds0);
    tot_part += cir.partidas;
    tot_fund += cir.fundidas;

    // O TSUR VEM ANTES de a fracao do proximo passo ser escrita, e por isso
    // esta' aqui e nao depois da medida: suavizar depois deixaria o campo um
    // passo atrasado em relacao a' superficie que o gerou.
    if (tsur_cada > 0 && (k + 1) % tsur_cada == 0)
      tot_tsur += frente.tsur(tsur_rlx, tsur_lim);

    if ((k + 1) % cada_med == 0 || k + 1 == passos) {
      escreve_fracao(frente, (*M->m), M->frac.get());
      const Medidas md = mede((*M->m), (*M->fd), (*M->frac), M->ns->velocidade(), frente);
      fim = md;
      if (md.circ < circ_min) { circ_min = md.circ; t_circ = M->ns->instante(); }
      if (md.vc > v_max) { v_max = md.vc; t_vmax = M->ns->instante(); }
      double umax = 0.0;
      for (int d = 0; d < 2; d++)
        for (long q = 0; q < M->fd->num_facetas(d); q++)
          umax = std::max(umax, std::fabs(M->ns->velocidade().em(d, q)));
      const double cfl = umax * dt / h;
      const double dv = M->ns->divergencia_maxima();

      if (rank == 0) {
        std::printf("%10.4f %12.6f %12.6f %12.6f %10.6f %6ld %8.4f %8.4f "
                    "%7.0e %5d %9.2e%s\n",
                    M->ns->instante(), md.circ, md.yc, md.vc, md.area / A0,
                    frente.num_marcadores(), umax, cfl, dv,
                    M->ns->iteracoes_da_pressao(), frente.centroide()[0] - 0.5,
                    (cfl > 1.0) ? "  CFL>1" : (dv > 1e-8 ? "  DIV!" : ""));
        std::fflush(stdout);
      }
    }
    if (cada_vtk > 0 && (k + 1) % cada_vtk == 0) grava(quadro++);

    // A MEDIDA E O VTK VEM ANTES DO REMALHAMENTO, e a ordem custou um
    // diagnostico.  Com o remalhamento antes, a medida lia uma montagem
    // RECEM-CRIADA que nunca deu passo: `its` saia 0 -- o Poisson dela nao
    // tinha rodado -- e `div` saia 4e-02, que e' a divergencia do campo
    // TRANSFERIDO, nao a do campo resolvido.  Interpolar nao preserva
    // divergencia discreta; quem a remove e' o Poisson do passo SEGUINTE.
    //
    // Nao era defeito do solver, mas tambem nao era ruido: com essa ordem o
    // medidor `div` deixava de responder a pergunta para a qual existe -- se a
    // projecao continua fechando depois de um remalhamento --, porque passava
    // a mostrar sempre o campo de antes dela.
    // ---- O REMALHAMENTO -----------------------------------------------
    //
    // VEM DEPOIS DE A FRENTE SE MOVER, e nao antes: a banda tem de cobrir
    // onde a interface ESTA', e remalhar antes do movimento a deixaria um
    // passo atrasada -- com a frente saindo da banda fina pela borda de cima
    // justamente no sentido em que ela sobe.
    //
    // A TABELA NOVA SO' SE COMPARA COM A VELHA; se derem iguais, a malha
    // seria a mesma e o remalhamento e' PULADO.  Nao e' so' economia: a
    // transferencia e' identidade quando a malha nao muda, mas pula-la deixa
    // o campo intocado por construcao em vez de por argumento.
    if (amr > 0 && amr_cada > 0 && (k + 1) % amr_cada == 0) {
      std::vector<signed char> nova = tabela_da_frente(frente);
      if (nova != tab) {
        const double e_antes = energia(*M);
        std::unique_ptr<Montagem> N = monta(&nova);
        transfere(*M, N.get());
        const double e_depois = energia(*N);
        if (e_antes > 0.0) {
          const double rel = (e_depois - e_antes) / e_antes;
          perda_total += rel;
          if (std::fabs(rel) > std::fabs(perda_pior)) perda_pior = rel;
        }
        M = std::move(N);
        tab = std::move(nova);
        remalhas++;
      }
    }

  }
  grava(quadro);

  if (rank == 0) {
    std::printf("\n%-28s %12s   %-22s\n", "grandeza", "medido", "referencia");
    std::printf("%-28s %12.4f   %-22s\n", "circularidade minima", circ_min,
                "0,9012 +- 0,0001");
    std::printf("%-28s %12.4f   %-22s\n", "instante do minimo", t_circ, "~1,90");
    std::printf("%-28s %12.4f   %-22s\n", "velocidade maxima", v_max,
                "0,2417 - 0,2421");
    std::printf("%-28s %12.4f   %-22s\n", "instante do maximo", t_vmax, "~0,92");
    std::printf("%-28s %12.4f   %-22s\n", "centro de massa no fim", fim.yc,
                "1,0799 - 1,0817 (t=3)");
    std::printf("%-28s %12.6f   %-22s\n", "area / area inicial", fim.area / A0,
                "1 (a bolha nao muda de area)");
    std::printf("\nmanutencao da frente: %ld arestas partidas, %ld fundidas,"
                " %ld marcadores movidos pelo TSUR\n",
                tot_part, tot_fund, tot_tsur);
    // O REMALHAMENTO E' CONTADO E DITO, e isso nao e' estatistica.  Uma
    // corrida com `--amr` em que a conta sai ZERO rodou com malha ESTATICA --
    // a banda fina ficou onde a bolha estava no inicio e a interface saiu dela
    // --, e o resultado seria plausivel e errado.  Ja' aconteceu neste
    // sistema: o FT do HiGFlow rodou 6479 passos com a frente FIXA e nada no
    // log dizia.
    if (amr > 0)
      std::printf("remalhamentos: %ld em %ld passos (a cada %d), malha final "
                  "com %ld celulas\n",
                  remalhas, passos, amr_cada, M->m->num_folhas_globais());
      // O NUMERO QUE DECIDE SE A TRANSFERENCIA E' A CAUSA.  Soma dos saltos
      // relativos de energia, e o pior salto.  Se a soma for negativa e da
      // ordem do deficit de `v_max`, a transferencia e' o mecanismo; se for
      // ruido em torno de zero, a causa esta' noutro lugar e eu estava
      // olhando para o lado errado.
      std::printf("energia na transferencia: soma dos saltos %+.4e, "
                  "pior salto %+.3e\n", perda_total, perda_pior);
    std::printf(
        "\nLEITURA: a discordancia ENTRE as tres referencias publicadas e' o\n"
        "limite de resolucao deste oraculo.  Bater melhor que ela nao seria\n"
        "conquista -- seria sinal de que algo foi ajustado para dar o numero.\n");
  }
  return 0;
}
