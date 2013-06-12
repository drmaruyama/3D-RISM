#include <iostream>
#include "rism3d.h"

void RISM3D :: cal_exchem (double * & xmu) {

  double dxmu = 0.0;
  if (clos == 0) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double dxmuv = 0.0;
#pragma omp parallel for reduction(+: dxmuv)
      for (int ig = 0; ig < ce -> mgrid; ++ig) {
	double h = huv[iv][ig].real();
	double cuv = h - tuv[iv][ig];
	if (h > 0.0 ) {
	  dxmuv -= (h * 0.5 + 1.0) * cuv;
	} else {
	  dxmuv += h * 0.5 * tuv[iv][ig] - cuv;
	}
      }
      xmu[iv] = dxmuv;
    }
  } else if (clos == 1) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double dxmuv = 0.0 ;
#pragma omp parallel for reduction(+: dxmuv)
      for (int ig = 0; ig < ce -> mgrid; ++ig) {
	double h = huv[iv][ig].real();
	dxmuv += h * 0.5 * tuv[iv][ig] + tuv[iv][ig] - h;
      }
      xmu[iv] = dxmuv;
    }
  } 
  double * xmu2;
  if (myrank == 0) xmu2 = new double[sv -> natv];
  MPI_Reduce(xmu, xmu2, sv -> natv, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      xmu[iv] = xmu2[iv] * ce -> dv * sv -> rhov[iv];
    }
    delete[] xmu2;
  }
} 
