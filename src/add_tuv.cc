#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"

void RISM3D :: add_tuv (double cuf) {

  for (int iv = 0; iv < sv -> natv; ++iv) {
    double q = sv -> qv[iv] * cuf;
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      tuv[iv][ig] += q * fr[ig];
    }
  }
} 
