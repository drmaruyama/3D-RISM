#include <iostream>
#include <fstream>
#include <mpi.h>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  int procs, myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  RISM3D * system;

  system = new RISM3D;

  if (argc == 1) {
    cout << "No parameter file!" << endl;
    return (1);
  } 

  if (argc > 2) {
    cout << "Too much arguments!" << endl;
    return(1);
  }
  system -> initialize(argv[1]);
  system -> iterate();
  system -> output();    

  MPI_Finalize();
  return(0);
}
