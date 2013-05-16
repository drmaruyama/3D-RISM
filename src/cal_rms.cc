#include "rism3d.h"

double RISM3D :: cal_rms () {
  double rms = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for reduction(+: rms)
    for (int ig = 0; ig < ce -> ngrid; ++ig){
      rms += tuvdif[iv][ig] * tuvdif[iv][ig];
    }
  }
  rms = sqrt (rms / (ce -> ngrid * sv -> natv));
  return rms;
}
