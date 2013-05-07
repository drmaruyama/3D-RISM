#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "version.h"

void RISM3D :: read_input (char inputfile[]) {
  int num;

  if (myrank == 0) {
    ifstream in_file;
    in_file.open (inputfile);

    cout << "reading input data file:  " << inputfile << endl;

    string check;
    in_file >> outlist >> co -> ksave >> check;
    if (check != version) {
      cout << "This input file is for old version." << endl;
      exit (1);
    }

    string closure;
    in_file >> closure;
    if (closure == "KH") {
      clos = 0;
    } else if (closure == "HNC") {
      clos = 1;
    } else {
      cout << "3D-RISM: unexpected closure switch " << endl;
      exit(1);
    }
    in_file >> fsolvent;
    in_file >> co -> convergence >> co -> maxstep;
    in_file >> ma -> count >> ma -> m >> ma -> mp;
    in_file >> ce -> box[0] >> ce -> box[1] >> ce -> box[2];
    in_file >> ce -> grid[0] >> ce -> grid[1] >> ce -> grid[2];

    in_file >> num;
    su -> init(num);

    for (int iu = 0; iu < su -> num; ++iu) {
      int n = iu * 3;
      in_file >> su -> q[iu] >> su -> sig[iu] >> su -> eps[iu]
	      >> su -> r[n] >> su -> r[n + 1]
	      >> su -> r[n + 2];
    }

    in_file.close ();
  }
  MPI_Bcast(&clos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myrank != 0) su -> init(num);
  co -> setup_mpi();
  ma -> setup_mpi();
  ce -> setup();
  su -> setup_mpi();
}
