#include <iostream>
#include <fstream>
#include "alloc.h"
#include "rism3d.h"

double RISM3D :: cal_exnew () {
  double excp_fd(double, double);

  double xmu1 = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double tmp = 0.0;
#pragma omp parallel for reduction(+: tmp)
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      double h = huv[0][ig].real();
      tmp -= h; 
    }
    xmu1 += tmp * sv -> rhov[iv];
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(huv[iv], - 1);
  }

  vector <complex <double> *> work;
  alloc2D (work, sv -> natv, ce -> ngrid);

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      work[iv][ig] = complex <double> (0.0, 0.0);
    }
  }

  for (int iv2 = 0; iv2 < sv -> natv; ++iv2) {	
    for (int iv1 = 0; iv1 < sv -> natv; ++iv1) {
#pragma omp parallel for
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	work[iv2][ig] += sv -> cvva[iv2][iv1][indga[ig]] 
	  * sv -> rhov[iv1] * huv[iv1][ig];
      }
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(work[iv], 1);
    fft -> execute(huv[iv], 1);
  }

  double xmu2 = 0.0;
  double xmu3 = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double dxmu2 = 0.0;
    double dxmu3 = 0.0;
#pragma omp parallel for reduction(+: dxmu2, dxmu3)
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      double ch = work[iv][ig].real();
      dxmu2 += ch;
      dxmu3 += ch * huv[iv][ig].real();
    }
    xmu2 += dxmu2 * sv -> rhov[iv];
    xmu3 += dxmu3 * sv -> rhov[iv];
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(huv[iv], - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
//      work[iv][ig] = sv -> wfka[iv][indga[ig]] * sv -> rhov[iv] * huv[iv][ig];
      work[iv][ig] = sv -> wfka[iv][indga[ig]] * sv -> rhov[0] * huv[iv][ig];
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(huv[iv], 1);
    fft -> execute(work[iv], 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for 
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
//      work[iv][ig] += sv -> rhov[iv];
      work[iv][ig] += sv -> rhov[0];
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    double bexcp_fd = excp_fd(sv -> pfhs[iv], sv -> rhov[0]);
    double excp0_st = sv -> rhov[0] * bexcp_fd * sv -> wfk0[iv];
    double a = sv -> pfhs[iv] / sv -> rhov[0];
#pragma omp parallel for 
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      double dens2 = work[iv][ig].real();
      double pfhs2 = a * dens2;
      double fdr = excp_fd(pfhs2, dens2);
//      rfd[iv][ig] = sv -> rhov[iv] * (guv[iv][ig].real() * fdr - bexcp_fd);
      work[iv][ig] = sv -> rhov[0] 
	* ((huv[iv][ig].real() + 1.0) * fdr - bexcp_fd);
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(work[iv], - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for 
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      work[iv][ig] *= sv -> wfka[iv][indga[ig]];
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(work[iv], 1);
  }

  double eda = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    double bexcp_fd = excp_fd(sv -> pfhs[iv], sv -> rhov[0]);
    double excp0_st = sv -> rhov[0] * bexcp_fd * sv -> wfk0[iv];
    double deda = 0.0;
#pragma omp parallel for reduction(+: deda)
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      deda += (work[iv][ig].real() + excp0_st) * (huv[iv][ig].real() + 1.0) 
	- excp0_st;
    }
    eda += deda * sv -> rhov[iv];
  }

  double xmu = (xmu1 + xmu2 + 0.5 * xmu3 - eda) * ce -> dv;

  return xmu;
} 
