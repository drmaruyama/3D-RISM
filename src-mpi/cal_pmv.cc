#include <iostream>
#include "rism3d.h"

double RISM3D :: cal_pmv () {

  double cuv = 0.0;
  double cuva = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for reduction(+: cuv)
    for (int ig = 0; ig < ce -> mgrid; ++ig) {
      cuv += (huv[iv][ig].real() - tuv[iv][ig]) * sv -> rhov[iv];
    }
  }
  MPI_Reduce(&cuv, &cuva, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  cuv = cuva * ce -> dv;
  double pmv = sv -> xikt * (1.0 - cuv);

  return pmv;
}
