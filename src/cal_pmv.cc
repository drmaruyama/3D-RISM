#include <iostream>
#include "rism3d.h"

double RISM3D :: cal_pmv () {

  double cuv = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double cuv0 = 0.0;
#pragma omp parallel for reduction(+: cuv0)
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      cuv0 += (huv[iv][ig].real() - tuv[iv][ig]);
    }
    cuv += cuv0 * sv -> rhov[iv];
  }
  cuv = cuv * ce -> dv;
  double pmv = sv -> xikt * (1.0 - cuv);

  return pmv;
}
