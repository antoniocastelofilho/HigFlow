// Sonda: quanto halo o DMPlexDistribute entrega por nivel de sobreposicao?
//
// ATENCAO AO QUE ISTO MEDE.  A sobreposicao e' da malha LAGRANGEANA -- a
// superficie do corpo.  Ela da' a cada rank triangulos vizinhos dos seus, que e'
// o que area dual e normal do vertice exigem na fronteira da particao.  NAO tem
// relacao com o suporte do delta regularizado, que e' euleriano e e' coberto
// pela franja da HiGTree mais o PetscSF (ver sonda-petscsf.c).
#include <petsc.h>
#include <petscdmplex.h>

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    PetscMPIInt rank, size;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    for (PetscInt sobrep = 0; sobrep <= 2; sobrep++) {
        // Esfera triangulada: superficie de dimensao 2 em espaco 3, que e' a
        // forma de um corpo.  Refinada, para a particao ter interior.
        DM dm;
        PetscCall(DMPlexCreateSphereMesh(PETSC_COMM_WORLD, 2, PETSC_TRUE, 1.0, &dm));
        PetscCall(DMSetFromOptions(dm));

        DM dmd = NULL;
        PetscSF sf = NULL;
        PetscCall(DMPlexDistribute(dm, sobrep, &sf, &dmd));
        if (dmd) { PetscCall(DMDestroy(&dm)); dm = dmd; }

        PetscInt ci, cf, vi, vf;
        PetscCall(DMPlexGetHeightStratum(dm, 0, &ci, &cf));   // celulas (triangulos)
        PetscCall(DMPlexGetDepthStratum (dm, 0, &vi, &vf));   // vertices

        // Fantasmas: as FOLHAS do pointSF sao os pontos que pertencem a outro rank.
        PetscSF psf;
        PetscCall(DMGetPointSF(dm, &psf));
        PetscInt nraiz, nfolha;
        const PetscInt *folhas = NULL;
        PetscCall(PetscSFGetGraph(psf, &nraiz, &nfolha, &folhas, NULL));

        PetscInt cel_fant = 0, ver_fant = 0;
        for (PetscInt i = 0; i < nfolha; i++) {
            const PetscInt p = folhas ? folhas[i] : i;
            if (p >= ci && p < cf) cel_fant++;
            if (p >= vi && p < vf) ver_fant++;
        }

        PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "  sobrep %d | rank %d | celulas %3d (%2d fantasmas) | vertices %3d (%2d fantasmas)\n",
            (int) sobrep, rank, (int)(cf - ci), (int) cel_fant,
            (int)(vf - vi), (int) ver_fant));
        PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

        if (sf) PetscCall(PetscSFDestroy(&sf));
        PetscCall(DMDestroy(&dm));
    }

    PetscCall(PetscFinalize());
    return 0;
}
