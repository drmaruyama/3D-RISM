#include "rism3d.h"

void RISM3D :: cal_sum (double & ovl00, double & sum) {
  double ovl = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for reduction(+: ovl)
    for (int ig = 0; ig < ce -> ngrid; ++ig){
      ovl += tuvdif[iv][ig] * tuvdif[iv][ig];
    }
  }
  ovl00 = ovl;
  sum = sqrt (ovl00 / (ce -> ngrid * sv -> natv));
}
