#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);
  
  if (outlist.find("m") != string::npos) {
    double pmv = cal_pmv();
    valarray <double> xmu = cal_exchem();
    double dft = cal_exnew();
    output_xmu(xmu, dft, pmv);
  }

  if (outlist.find("d") != string::npos) {
    double * du;
    du = new double[su -> num * 3];
    cal_grad(du);
    output_grad(du);
    delete[] du;
  }

  if (outlist.find("c") != string::npos) {
    output_cuv();
  }

  if (outlist.find("g") != string::npos) {
    output_guv();
  }

  if (outlist.find("h") != string::npos) {
    output_huv();
  }
} 
