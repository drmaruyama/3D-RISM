#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "version.h"

void RISM3D :: read_input (char inputfile[]) {

  ifstream in_file;
  in_file.open (inputfile);

  cout << "reading input data file:  " << inputfile << endl;

  string check;
  in_file >> outlist >> co -> ksave >> check;
  if (check != version) {
    cout << "This input file is for old version." << endl;
    exit (1);
  }

  in_file >> closure;
  in_file >> fsolvent;
  in_file >> co -> convergence >> co -> maxstep;
  in_file >> ma -> count >> ma -> m >> ma -> mp;
  in_file >> ce -> box[0] >> ce -> box[1] >> ce -> box[2];
  in_file >> ce -> grid[0] >> ce -> grid[1] >> ce -> grid[2];

  int num;
  in_file >> num;

  ce -> setup();
  su -> init(num);

  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    in_file >> su -> q[iu] >> su -> sig[iu] >> su -> eps[iu]
	    >> su -> r[num] >> su -> r[num + 1]
	    >> su -> r[num + 2];
  }

  in_file.close ();
}