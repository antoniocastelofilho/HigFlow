// Quantos nao-zeros por linha cada esquema produz.
//
// EXISTE POR CAUSA DE UMA HIPOTESE QUE ELE REFUTOU.  O ganho superlinear do MLS
// de np=1 para np=8 (39x POR ITERACAO) foi atribuido, provisoriamente, a' matriz
// nao caber no cache -- com a estimativa de 300 nao-zeros por linha.  Medido:
// sao 19,6.  A matriz do nivel 7 tem 3,85 MB e CABE no L3 de 16 MiB desta
// maquina, logo a explicacao nao se sustenta como estava.
//
// E' LENTO DE PROPOSITO SIMPLES: `linha_montada` remonta o indice a cada
// chamada, o que o torna O(n^2).  Serve para malhas pequenas, e a densidade por
// linha NAO depende do tamanho da malha -- medir no nivel 5 e extrapolar e'
// legitimo, medir no nivel 7 custa minutos sem acrescentar nada.
//
// Uso:  mpirun -np 1 ./bench/densidade [nivel]

#include <cstdio>
#include <map>
#include "arranque.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "manufaturada2d.hxx"
#include "poisson2d.hxx"
using namespace t8flow2d;
int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  const int nivel = (argc > 1) ? atoi(argv[1]) : 7;
  const Manufaturada &sol = manufaturada("grau2");
  for (int caso = 0; caso < 2; caso++) {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
    ContornoPorRegiao bc;
    bc.padrao(Contorno::dirichlet([&sol](const Ponto2 &p) { return sol.u(p); }));
    Poisson po(m, bc);
    po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
    if (caso == 1) { po.esquema(Esquema::MlsGauss); po.mls(3, 60); }
    long nnz = 0, linhas = 0, maior = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        double r = 0.0;
        const long n = static_cast<long>(po.linha_montada(b, i, &r).size());
        nnz += n; linhas++; if (n > maior) maior = n;
      }
    // 12 bytes por nao-zero: 8 do valor (double) e 4 do indice de coluna (int),
    // que e' o que o formato CSR do PETSc guarda.
    const double mb = nnz * 12.0 / (1024.0 * 1024.0);
    std::printf("%-9s nivel=%d linhas=%ld nnz=%ld media=%.1f maior=%ld matriz=%.1f MB\n",
                (caso == 0) ? "fluxo2" : "mlsgauss", nivel, linhas, nnz,
                double(nnz) / linhas, maior, mb);
  }
  return 0;
}
