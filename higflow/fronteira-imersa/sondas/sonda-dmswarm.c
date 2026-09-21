// Sonda: o DMSwarm desta instalacao do PETSc sabe MIGRAR particula entre ranks?
// E' a capacidade que a malha lagrangeana de corpo movel exige e que o
// distributed_property da HiGTree nao tem.
#include <petsc.h>
#include <petscdmswarm.h>

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    PetscMPIInt rank, size;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    DM swarm;
    PetscCall(DMCreate(PETSC_COMM_WORLD, &swarm));
    PetscCall(DMSetType(swarm, DMSWARM));
    PetscCall(DMSetDimension(swarm, 3));
    PetscCall(DMSwarmSetType(swarm, DMSWARM_BASIC));

    // Campos por marcador: forca e peso -- o que a fronteira imersa carrega.
    PetscCall(DMSwarmRegisterPetscDatatypeField(swarm, "forca", 3, PETSC_REAL));
    PetscCall(DMSwarmRegisterPetscDatatypeField(swarm, "peso",  1, PETSC_REAL));
    PetscCall(DMSwarmFinalizeFieldRegister(swarm));

    const PetscInt nlocal = 4;
    PetscCall(DMSwarmSetLocalSizes(swarm, nlocal, 4));

    PetscReal *forca;
    PetscCall(DMSwarmGetField(swarm, "forca", NULL, NULL, (void **) &forca));
    for (PetscInt i = 0; i < nlocal; i++)
        for (int d = 0; d < 3; d++) forca[3*i + d] = 100.0*rank + i + 0.1*d;
    PetscCall(DMSwarmRestoreField(swarm, "forca", NULL, NULL, (void **) &forca));

    PetscInt antes;
    PetscCall(DMSwarmGetLocalSize(swarm, &antes));

    // Manda todo marcador para o rank seguinte -- e' a migracao que interessa.
    PetscInt *destino;
    PetscCall(DMSwarmGetField(swarm, DMSwarmField_rank, NULL, NULL, (void **) &destino));
    for (PetscInt i = 0; i < antes; i++) destino[i] = (rank + 1) % size;
    PetscCall(DMSwarmRestoreField(swarm, DMSwarmField_rank, NULL, NULL, (void **) &destino));

    PetscCall(DMSwarmMigrate(swarm, PETSC_TRUE));

    PetscInt depois;
    PetscCall(DMSwarmGetLocalSize(swarm, &depois));

    // O conteudo sobreviveu a' viagem?
    PetscCall(DMSwarmGetField(swarm, "forca", NULL, NULL, (void **) &forca));
    const PetscReal esperado = 100.0 * ((rank - 1 + size) % size);
    PetscBool ok = PETSC_TRUE;
    for (PetscInt i = 0; i < depois; i++)
        if (forca[3*i] < esperado - 0.5 || forca[3*i] > esperado + 3.5) ok = PETSC_FALSE;
    PetscCall(DMSwarmRestoreField(swarm, "forca", NULL, NULL, (void **) &forca));

    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,
        "  rank %d: antes %d, depois %d, carga do rank %d %s\n",
        rank, (int) antes, (int) depois, (int)((rank - 1 + size) % size),
        ok ? "INTACTA" : "CORROMPIDA"));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));

    PetscCall(DMDestroy(&swarm));
    PetscCall(PetscFinalize());
    return 0;
}
