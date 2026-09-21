// Sonda: o PetscSF acumula na direcao que falta -- FRANJA -> DONO?
//
// E' o buraco que a escolha por estrutura distribuida reabre.  O espalhamento do
// delta regularizado escreve em facetas de OUTRO rank; o dp_sync da HiGTree
// (pdomain.c:1308) manda dono->franja e SOBRESCREVE, entao a contribuicao
// escrita numa franja seria descartada, calada.
//
// Aqui se imita a situacao: cada rank tem entradas proprias e UMA entrada de
// franja que pertence ao vizinho, escreve nas duas, e o PetscSFReduce com
// MPI_SUM deve somar a franja no dono.
#include <petsc.h>
#include <petscsf.h>

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    PetscMPIInt rank, size;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    const PetscInt nproprias = 4;          // "facetas que este rank possui"
    const PetscInt nfranja   = 1;          // uma facadeta de franja, do vizinho

    // A franja deste rank aponta para a entrada 0 do rank seguinte -- que e'
    // como o gid_map da HiGTree ja' descreve dono e posicao.
    PetscSFNode franja[1];
    franja[0].rank  = (rank + 1) % size;
    franja[0].index = 0;

    PetscSF sf;
    PetscCall(PetscSFCreate(PETSC_COMM_WORLD, &sf));
    PetscCall(PetscSFSetGraph(sf, nproprias, nfranja, NULL, PETSC_COPY_VALUES,
                              franja, PETSC_COPY_VALUES));
    PetscCall(PetscSFSetUp(sf));

    PetscReal dono[4], folha[1];
    for (PetscInt i = 0; i < nproprias; i++) dono[i] = 1.0;   // contribuicao local
    folha[0] = 10.0;                                          // contribuicao no vizinho

    PetscCall(PetscSFReduceBegin(sf, MPIU_REAL, folha, dono, MPI_SUM));
    PetscCall(PetscSFReduceEnd  (sf, MPIU_REAL, folha, dono, MPI_SUM));

    // A entrada 0 deve ter 1 (propria) + 10 (do rank anterior) = 11.
    // As outras, 1.  Se o PetscSF sobrescrevesse em vez de somar, daria 10.
    const PetscBool ok = (PetscBool)(PetscAbsReal(dono[0] - 11.0) < 1e-12 &&
                                     PetscAbsReal(dono[1] -  1.0) < 1e-12);

    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,
        "  rank %d: entrada 0 = %.1f (esperado 11 = 1 propria + 10 do rank %d), "
        "entrada 1 = %.1f   %s\n",
        rank, (double) dono[0], (rank - 1 + size) % size, (double) dono[1],
        ok ? "ACUMULOU" : "NAO ACUMULOU"));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

    PetscCall(PetscSFDestroy(&sf));
    PetscCall(PetscFinalize());
    return 0;
}
