#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"

void RISM3D :: initialize_tuv () {

  cout << "synthesizing initial estimate for Tuv ..." << endl ;
  
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      tuv[iv][ig] = sv -> qv[iv] * fr[ig];
    }
  }
} 
