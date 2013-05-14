#include "rism3d.h"

double RISM3D :: cal_rms () {
  double rms = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for reduction(+: rms)
    for (int ig = 0; ig < ce -> mgrid; ++ig){
      rms += tuvdif[iv][ig] * tuvdif[iv][ig];
    }
  }
  double rms00;
  MPI_Allreduce(&rms, &rms00, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rms = sqrt (rms00 / (ce -> ngrid * sv -> natv));
  return rms;
}
