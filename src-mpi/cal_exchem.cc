#include <iostream>
#include "rism3d.h"

valarray <double> RISM3D :: cal_exchem () {

  double dxmu = 0.0;
  valarray <double> xmu;
  xmu.resize(sv -> natv);
  if (closure == "KH") {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double dxmuv = 0.0;
#pragma omp parallel for reduction(+: dxmuv)
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	double h = huv[iv][ig].real();
	double cuv = h - tuv[iv][ig];
	if (h > 0.0 ) {
	  dxmuv -= (h * 0.5 + 1.0) * cuv;
	} else {
	  dxmuv += h * 0.5 * tuv[iv][ig] - cuv;
	}
      }
      xmu[iv] = sv -> rhov[iv] * dxmuv;
    }
  } else if (closure == "HNC") {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double dxmuv = 0.0 ;
#pragma omp parallel for reduction(+: dxmuv)
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	double h = huv[iv][ig].real();
	dxmuv += h * 0.5 * tuv[iv][ig] + tuv[iv][ig] - h;
      }
      xmu[iv] = sv -> rhov[iv] * dxmuv;
    }
  } else {
    cout << "EXCHEM: unexpected closure switch " << closure << endl ;
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    xmu[iv] = xmu[iv] * ce -> dv;
  }
  return xmu;
} 
