#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);

  int flag = 0;
  if (myrank == 0) {
    if (outlist.find("m") != -1) flag += 1;
    if (outlist.find("d") != -1) flag += 2;
    if (outlist.find("c") != -1) flag += 4;
    if (outlist.find("g") != -1) flag += 8;
    if (outlist.find("h") != -1) flag += 16;
  }

  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ((flag & 1) == 1) {
    double pmv = cal_pmv();
    double * xmu = new double[sv -> natv];
    cal_exchem(xmu);
    if (myrank == 0) {
      output_xmu(xmu, pmv);
    }
    delete[] xmu;
  }

  if ((flag & 2) == 2) {
    double * du;
    if (myrank == 0) du = new double[su -> num * 3];
    cal_grad(du);
    if (myrank == 0) {
      output_grad(du);
      delete[] du;
    }
  }

  if ((flag & 4) == 4) {
    output_cuv();
  }

  if ((flag & 8) == 8) {
    output_guv();
  }

  if ((flag & 16) == 16) {
    output_huv();
  }
} 
