#include <iostream>
#include <fstream>
#include "rism3d.h"

int main (int argc, char * argv[]) {
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

  //  MPI_Init(&argc, &argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  system -> initialize(argv[1]);
  system -> iterate();
  system -> output();    

  MPI_Finalize();
  return(0);
}
