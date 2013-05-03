#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);
  
  if (outlist.find("m") != -1) {
    double pmv = cal_pmv();
    valarray <double> xmu = cal_exchem();
    output_xmu(xmu, pmv);
  }

  if (outlist.find("d") != -1) {
    double * du;
    du = new double[su -> num * 3];
    cal_grad(du);
    output_grad(du);
    delete[] du;
  }

  if (outlist.find("c") != -1) {
    output_cuv();
  }

  if (outlist.find("g") != -1) {
    output_guv();
  }

  if (outlist.find("h") != -1) {
    output_huv();
  }
} 
