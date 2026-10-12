// Quanto custa a franja de K camadas.
//
// A PERGUNTA, dita com precisao: a franja de K camadas foi construida para que o
// esquema de alta ordem funcionasse em paralelo (ver as clausulas X1..X9), e o
// custo dela nunca foi medido.  A medida de escalabilidade que o projeto tem --
// 1,6x em np=4 -- PRECEDE a franja, logo nao vale mais.
//
// O QUE ESTE ARREIO MEDE, e por que cada coisa:
//
//   t_constroi    a malha com a franja NATIVA do t8code (uma camada).  E' a
//                 linha de base: o custo que existiria sem este trabalho;
//   t_franja(K)   `exige_camadas(K)`, que refaz a franja inteira com K camadas.
//                 Refaz e nao estende, de proposito -- estender exigiria estado
//                 incremental, e estado incremental de franja e' a classe de
//                 defeito que esta frente passou o dia consertando;
//   t_poisson     resolver o sistema.  E' a REFERENCIA que da' sentido ao resto:
//                 um custo de franja de 50 ms nao diz nada; 50 ms contra 2 s de
//                 solver diz;
//   t_reparticao  `reparticiona()`, que o refino local torna necessario.
//
// E O VOLUME PUBLICADO, sempre ao lado do tempo.  A carga por folha e' `(2K+1)^2`
// entradas, logo o volume cresce com o QUADRADO de K; tempo sem volume e' numero
// sem mecanismo -- nao se sabe se dobrou por comunicacao, por memoria, ou por ter
// medido outra coisa.
//
// DISCIPLINA DE MEDIDA, e cada item ja' produziu um resultado falso neste projeto:
//
//   1. o tempo e' o do rank MAIS LENTO (`MPI_MAX`), nao a media.  O relogio de
//      parede de um programa paralelo e' o do ultimo a terminar;
//   2. barreira ANTES de cada cronometragem, senao a medida inclui o atraso de
//      quem chegou tarde na etapa anterior;
//   3. repeticoes com o MENOR valor valendo.  A maquina tem vizinhos -- outras
//      corridas --, e o minimo e' menos sensivel a eles que a media;
//   4. a primeira repeticao e' DESCARTADA.  A primeira toca memoria fria e paga
//      alocacao de buffers que as seguintes reusam, e incluir isso mede o
//      aquecimento em vez do regime;
//   5. a malha e' reconstruida a cada repeticao, porque `exige_camadas` e'
//      idempotente: medir a segunda chamada mediria um `return` imediato.
//
// Uso:  mpirun -np N ./bench/franja [--nivel=L] [--camadas=K] [--fracao=F] [--reps=R]

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <mpi.h>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "manufaturada2d.hxx"
#include "poisson2d.hxx"

using namespace t8flow2d;

namespace {

int g_np = 1;
int g_rank = 0;

//! Imprime uma medida com as contagens de carga ao lado.
//
// `locmin` e `locmax` vao juntos porque o desequilibrio da particao e' o que
// explica o tempo do rank mais lento: 2,8x de desequilibrio foi medido depois de
// refino local, e sozinho explicava a perda de escalabilidade.
void relata(const char *etapa, int nivel, int camadas, double fracao, long globais,
            long locais, long franja, long bytes, double seg)
{
  long lmin = 0, lmax = 0, fmax = 0, btot = 0;
  double tmax = 0.0;
  MPI_Allreduce(&locais, &lmin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&locais, &lmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&franja, &fmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&bytes, &btot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&seg, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if (g_rank == 0) {
    printf("medida %s nivel=%d np=%d camadas=%d fracao=%.2f globais=%ld "
           "locmin=%ld locmax=%ld franjamax=%ld bytes=%ld seg=%.6f\n",
           etapa, nivel, g_np, camadas, fracao, globais, lmin, lmax, fmax, btot, tmax);
    fflush(stdout);
  }
}

//! Constroi a malha do caso: uniforme ou com refino aleatorio.
Malha *faz_malha(int nivel, double fracao)
{
  Malha *m = new Malha({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  if (fracao > 0.0) m->refina_aleatorio(fracao, 2026, nivel + 1);
  return m;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &g_np);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

  int nivel = 6;
  int camadas_max = 4;
  double fracao = 0.0;
  int reps = 4;
  bool com_poisson = true;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--camadas=", 10)) camadas_max = std::atoi(a + 10);
    else if (!std::strncmp(a, "--fracao=", 9)) fracao = std::atof(a + 9);
    else if (!std::strncmp(a, "--reps=", 7)) reps = std::atoi(a + 7);
    else if (!std::strcmp(a, "--sem-poisson")) com_poisson = false;
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (g_rank == 0)
        std::printf("uso: mpirun -np N ./bench/franja [--nivel=L] [--camadas=K]\n"
                    "                [--fracao=F] [--reps=R] [--sem-poisson]\n");
      return 0;
    }
  }
  if (reps < 2) reps = 2;   // uma e' descartada; precisa de pelo menos duas

  if (g_rank == 0) {
    std::printf("config nivel=%d np=%d camadas_max=%d fracao=%.2f reps=%d\n",
                nivel, g_np, camadas_max, fracao, reps);
    std::fflush(stdout);
  }

  // ---- LINHA DE BASE: a malha com a franja nativa, de uma camada ------------
  {
    double melhor = 1e300;
    long glob = 0, loc = 0, fr = 0, by = 0;
    for (int r = 0; r < reps; r++) {
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      Malha *m = faz_malha(nivel, fracao);
      const double dt = MPI_Wtime() - t0;
      glob = m->num_folhas_globais();
      loc = m->num_folhas();
      fr = m->num_fantasmas();
      by = m->bytes_publicados_na_franja();
      delete m;
      if (r > 0 && dt < melhor) melhor = dt;   // a primeira e' aquecimento
    }
    relata("constroi", nivel, 1, fracao, glob, loc, fr, by, melhor);
  }

  // ---- O CUSTO DE CADA NUMERO DE CAMADAS -----------------------------------
  //
  // `k = 1` entra na varredura de proposito, e e' o CONTRASTE: `exige_camadas(1)`
  // nao faz rodada nenhuma e tem de sair em tempo desprezivel.  Se sair com
  // tempo, a medida nao esta' medindo as rodadas -- esta' medindo outra coisa.
  for (int k = 1; k <= camadas_max; k++) {
    double melhor = 1e300;
    long glob = 0, loc = 0, fr = 0, by = 0, sat = 0;
    for (int r = 0; r < reps; r++) {
      Malha *m = faz_malha(nivel, fracao);   // fora do cronometro
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      m->exige_camadas(k);
      const double dt = MPI_Wtime() - t0;
      glob = m->num_folhas_globais();
      loc = m->num_folhas();
      fr = m->num_fantasmas();
      by = m->bytes_publicados_na_franja();
      sat = m->saturacoes_da_franja();
      delete m;
      if (r > 0 && dt < melhor) melhor = dt;
    }
    relata("franja", nivel, k, fracao, glob, loc, fr, by, melhor);
    if (sat > 0 && g_rank == 0) {
      std::printf("AVISO camadas=%d: a carga SATUROU %ld vezes -- o tempo abaixo "
                  "e' de uma troca que publicou menos do que conhecia\n", k, sat);
      std::fflush(stdout);
    }
  }

  // ---- A REPARTICAO --------------------------------------------------------
  {
    double melhor = 1e300;
    long glob = 0, loc = 0, fr = 0;
    for (int r = 0; r < reps; r++) {
      Malha *m = faz_malha(nivel, fracao);
      m->exige_camadas(3);
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      m->reparticiona();
      const double dt = MPI_Wtime() - t0;
      glob = m->num_folhas_globais();
      loc = m->num_folhas();
      fr = m->num_fantasmas();
      delete m;
      if (r > 0 && dt < melhor) melhor = dt;
    }
    relata("reparticao", nivel, 3, fracao, glob, loc, fr, 0, melhor);
  }

  // ---- A REFERENCIA: resolver o Poisson ------------------------------------
  //
  // E' O NUMERO QUE DA' SENTIDO AOS OUTROS.  Sem ele, "a franja levou 50 ms" nao
  // e' resultado: nao se sabe se e' desprezivel ou dominante.  Dois esquemas,
  // porque eles pedem franjas diferentes -- o de fluxo pede duas camadas, o MLS
  // de ordem 3 pede tres -- e portanto a fracao e' diferente em cada um.
  if (com_poisson) {
    const Manufaturada &sol = manufaturada("grau2");
    const struct { const char *nome; Esquema e; int ordem; } casos[] = {
      {"poisson-fluxo2", Esquema::Fluxo2, 0},
      {"poisson-mlsgauss", Esquema::MlsGauss, 3},
    };
    for (const auto &c : casos) {
      double melhor = 1e300;
      long glob = 0, loc = 0, fr = 0;
      for (int r = 0; r < reps; r++) {
        Malha *m = faz_malha(nivel, fracao);
        ContornoPorRegiao bc;
        bc.padrao(Contorno::dirichlet([&sol](const Ponto2 &p) { return sol.u(p); }));
        Poisson po(*m, bc);
        po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
        po.esquema(c.e);
        if (c.ordem > 0) po.mls(c.ordem, 60);

        MPI_Barrier(MPI_COMM_WORLD);
        const double t0 = MPI_Wtime();
        po.resolve(Motor::Petsc, 1e-10);
        const double dt = MPI_Wtime() - t0;

        glob = m->num_folhas_globais();
        loc = m->num_folhas();
        fr = m->num_fantasmas();
        delete m;
        if (r > 0 && dt < melhor) melhor = dt;
      }
      relata(c.nome, nivel, 0, fracao, glob, loc, fr, 0, melhor);
    }
  }

  return 0;
}
