// O t8code sob particionamento real: np = 1, 2 e 3 -- C13 e P1 a P4.
//
// As sete clausulas ja' verificadas para o t8code rodam em UM processo.  Esta
// familia e' a que exige a floresta distribuida, e e' onde a substituicao do
// particionador de fato acontece: o MTree reparte pelo `lbal.c` com Zoltan, o
// t8code por curva de preenchimento de espaco.  Os dois tem de cumprir o mesmo
// contrato por caminhos que nao se parecem.
//
// Espelha o test-fringe-parallel, que faz o mesmo para o MTree, inclusive na
// forma: tudo por REDUCAO GLOBAL, com o rank 0 concluindo.  Um teste paralelo em
// que cada rank imprime o seu veredicto esconde a falha de um rank no meio das
// linhas dos outros.
//
//   P1  particao_cobre_o_dominio_uma_vez   soma global dos volumes = 1 e soma
//       global das celulas = total.  Pega celula perdida E celula contada duas
//       vezes, que contagem por rank nao distingue.
//   P2  ghost_existe_quando_ha_vizinho     com np>1 TODO rank tem ghost; com
//       np=1 nenhum.  Reduzido por MIN e por MAX.
//   P3  monta_em_serie                     em np=1 a floresta se monta e tem o
//       dominio inteiro.  Para o MTree isto precisou de correcao (o sfbi so' era
//       preenchido pelo caminho particionado); para o t8code e' nativo.
//   P4  valor_nao_depende_da_particao      campo linear num ponto FIXO da' o
//       mesmo valor em qualquer np -- o oraculo que separa "particiona
//       diferente" de "particiona errado".
//   C13 suporte_alcanca_o_ghost            com np>1, em todo rank existe celula
//       cujo vizinho de face esta' no ghost.  E' o criterio do que a franja
//       entrega, nao do tamanho dela.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "coord.h"
#include "utils.h"
#include "testing.h"
#include "t8code/t8-partition.h"

#define NC 8

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    t8_particao_monta();

    const long   n_local   = t8_particao_num_locais();
    const double vol_local = t8_particao_volume_local();
    const long   tem_ghost = (t8_particao_num_ghosts() > 0) ? 1 : 0;
    const long   alcanca   = t8_particao_suporte_alcanca_ghost();

    // Ponto ESTRITAMENTE dentro de uma celula, em todas as direcoes: sobre uma
    // face, dois ranks o reivindicariam e a reducao por soma somaria duas vezes
    // -- foi o que aconteceu no teste de franja do MTree, medido em np=3.
    Point alvo;
    POINT_ASSIGN_SCALAR(alvo, 0.59375);       // centro 0,5625 mais 1/4 de celula
    double v_local = 0.0;
    const long possui = t8_particao_interpola(alvo, &v_local) ? 1 : 0;
    if(!possui) v_local = 0.0;

    long   n_global, ghost_min, ghost_max, alcanca_min, donos;
    double vol_global, v_global;
    MPI_Allreduce(&n_local,   &n_global,    1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol_local, &vol_global,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_ghost, &ghost_min,   1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_ghost, &ghost_max,   1, MPI_LONG,   MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&alcanca,   &alcanca_min, 1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&possui,    &donos,       1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&v_local,   &v_global,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    long esperado = 1;
    for(int d = 0; d < DIM; d++) esperado *= NC;

    // So' o rank 0 conclui.  E ninguem chama MPI_Finalize: o higtree_initialize
    // registra o PetscFinalize no atexit, e ele ja' o chama.
    if(rank != 0) { t8_particao_encerra(); return 0; }

    t_case("particao_cobre_o_dominio_uma_vez");
    T_CHECK_MSG(n_global == esperado,
        "soma global das celulas locais = %ld, esperado %ld (np=%d): celula "
        "perdida ou contada duas vezes", n_global, esperado, ntasks);
    T_NEAR(vol_global, 1.0, 1.0e-12, "soma global dos volumes locais");

    t_case("ghost_existe_quando_ha_vizinho");
    if(ntasks == 1) {
        T_CHECK_MSG(ghost_max == 0,
            "com np=1 nenhum rank deveria ter ghost, e algum tem");
    } else {
        T_CHECK_MSG(ghost_min == 1,
            "com np=%d TODO rank deveria ter ghost, e ao menos um nao tem",
            ntasks);
    }

    t_case("monta_em_serie");
    if(ntasks == 1) {
        T_CHECK_MSG(n_local == esperado,
            "em np=1 o rank unico deveria ter o dominio inteiro (%ld celulas) e "
            "tem %ld", esperado, n_local);
    } else {
        T_CHECK_MSG(n_local > 0 && n_local < esperado,
            "com np=%d o rank 0 deveria ter uma PARTE do dominio, e tem %ld de "
            "%ld", ntasks, n_local, esperado);
    }

    t_case("suporte_alcanca_o_ghost");
    if(ntasks == 1) {
        T_CHECK_MSG(alcanca_min == 0,
            "com np=1 nenhum suporte deveria alcancar ghost");
    } else {
        T_CHECK_MSG(alcanca_min == 1,
            "com np=%d todo rank deveria ter celula cujo vizinho de face esta' "
            "no ghost, e ao menos um nao tem", ntasks);
    }

    t_case("valor_nao_depende_da_particao");
    T_CHECK_MSG(donos == 1,
        "o ponto alvo deveria pertencer a exatamente um rank, e pertence a %ld",
        donos);
    // O centroide da celula que contem 0,59375 e' 0,5625 em cada direcao.
    {
        double c = 0.5625;
        double esperado_v = 1.0 + 2.0 * c + 3.0 * c;
        T_NEAR(v_global, esperado_v, 1.0e-12,
               "campo linear no ponto fixo, independente da particao");
    }

    t8_particao_encerra();
    return t_end();
}
