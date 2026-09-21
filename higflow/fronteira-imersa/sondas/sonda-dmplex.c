// Sonda: o DMPlex distribui uma SUPERFICIE TRIANGULADA imersa em 3D -- topologia
// de dimensao 2 num espaco de dimensao 3?  E' o caso da malha lagrangeana 3D, e
// e' o que a estrutura da HiGTree nao sabe representar.
#include <petsc.h>
#include <petscdmplex.h>

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    PetscMPIInt rank, size;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    // Um octaedro: 8 triangulos, 6 vertices, fechado -- superficie de corpo.
    const PetscInt  ncel = 8, nver = 6, cantos = 3;
    PetscInt cels[24] = {0,1,2, 0,2,3, 0,3,4, 0,4,1,
                         5,2,1, 5,3,2, 5,4,3, 5,1,4};
    PetscReal coords[18] = { 0, 0, 1,   1, 0, 0,   0, 1, 0,
                            -1, 0, 0,   0,-1, 0,   0, 0,-1};

    DM dm = NULL;
    // topologia 2D (triangulos) mergulhada em espaco 3D
    PetscCall(DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, 2,
              rank == 0 ? ncel : 0, rank == 0 ? nver : 0, cantos,
              PETSC_TRUE, cels, 3, coords, &dm));

    PetscInt ini, fim;
    PetscCall(DMPlexGetHeightStratum(dm, 0, &ini, &fim));
    const PetscInt cel_antes = fim - ini;
    PetscCall(DMPlexGetDepthStratum(dm, 0, &ini, &fim));
    const PetscInt ver_antes = fim - ini;

    DM dmd = NULL;
    PetscSF sf = NULL;
    PetscCall(DMPlexDistribute(dm, 0, &sf, &dmd));
    if (dmd) { PetscCall(DMDestroy(&dm)); dm = dmd; }

    PetscCall(DMPlexGetHeightStratum(dm, 0, &ini, &fim));
    const PetscInt cel_dep = fim - ini;
    PetscCall(DMPlexGetDepthStratum(dm, 0, &ini, &fim));
    const PetscInt ver_dep = fim - ini;

    // A vizinhanca sobreviveu?  Area e normal do vertice dependem dela.
    PetscInt sup_max = 0;
    PetscCall(DMPlexGetDepthStratum(dm, 0, &ini, &fim));
    for (PetscInt v = ini; v < fim; v++) {
        PetscInt n; PetscCall(DMPlexGetSupportSize(dm, v, &n));
        if (n > sup_max) sup_max = n;
    }

    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,
        "  rank %d: antes %d cel / %d ver   depois %d cel / %d ver   suporte max de vertice %d\n",
        rank, (int) cel_antes, (int) ver_antes, (int) cel_dep, (int) ver_dep, (int) sup_max));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

    if (sf) PetscCall(PetscSFDestroy(&sf));
    PetscCall(DMDestroy(&dm));
    PetscCall(PetscFinalize());
    return 0;
}
