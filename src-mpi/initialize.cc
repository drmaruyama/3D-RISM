#include "rism3d.h"

void RISM3D :: initialize(char inputfile[]) {
  set_mpi();
  read_input(inputfile);
  set_fname(inputfile);
  initialize_g();
  set_solvent();
  cal_potential();
} 
