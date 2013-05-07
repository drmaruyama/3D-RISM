#include <iostream>
#include <fstream>
#include <string>

#include "solvent.h"

void Solvent :: setup_mpi(int myrank) {
  void alloc3D (vector <vector <double * > > &, int, int, int);

  MPI_Bcast(&temper, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xikt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&natv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntab, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (myrank != 0) {
    multv = new int[natv];
    rhov = new double[natv];
    qv = new double[natv];
    sigv = new double[natv];
    epsv = new double[natv];
    ttab = new double[ntab];
    alloc3D (xvv, natv, natv, ntab);
  }

  MPI_Bcast(multv, natv, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rhov, natv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(qv, natv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(sigv, natv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(epsv, natv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ttab, ntab, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int iv2 = 0; iv2 < natv; ++iv2) {
    for (int iv1 = 0; iv1 < natv; ++iv1) {
      MPI_Bcast(xvv[iv2][iv1], ntab, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }
}
