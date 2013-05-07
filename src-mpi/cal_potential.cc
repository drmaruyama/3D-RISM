#include "rism3d.h"

void RISM3D :: cal_potential() {
  double * euv;
  cal_LJ();
  cal_Coulomb(euv);
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> mgrid; ++ig) {
      uuv[iv][ig] += sv -> qv[iv] * euv[ig];
    }
  }
  delete[] euv;
}
