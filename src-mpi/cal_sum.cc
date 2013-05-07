#include "rism3d.h"

void RISM3D :: cal_sum (double & ovl00, double & sum) {
  double ovl = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for reduction(+: ovl)
    for (int ig = 0; ig < ce -> mgrid; ++ig){
      ovl += tuvdif[iv][ig] * tuvdif[iv][ig];
    }
  }
  MPI_Allreduce(&ovl, &ovl00, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum = sqrt (ovl00 / (ce -> ngrid * sv -> natv));
}
