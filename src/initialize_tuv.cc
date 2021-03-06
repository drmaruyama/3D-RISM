#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"

void RISM3D :: initialize_tuv (double cf) {

  cout << "synthesizing initial estimate for Tuv ..." << endl ;
  
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double q = sv -> qv[iv] * cf;
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      tuv[iv][ig] = q * fr[ig];
    }
  }
} 
