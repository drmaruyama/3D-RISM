#include "solute.h"

void Solute :: init(int n) {
  num = n;
  q = new double[num];
  sig = new double[num];
  eps = new double[num];
  r = new double[num * 3];
}

void Solute :: setup_mpi() {
  MPI_Bcast(q, num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(sig, num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(eps, num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(r, num * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
