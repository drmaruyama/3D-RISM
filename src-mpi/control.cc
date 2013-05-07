#include "control.h"

void Control :: setup_mpi() {
  MPI_Bcast(&convergence, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxstep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ksave, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
