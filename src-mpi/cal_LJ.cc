#include <iostream>
#include "rism3d.h"

void RISM3D :: cal_LJ() {
  void alloc2D (vector <double *> &, int, int);

  const double cut = 1.0e-2;
  const double cut2 = cut * cut;

  if (myrank == 0) {
    cout << "tabulating solute Lennard-Jones potential ..." << endl;
  }

  alloc2D (siguv, sv -> natv, su -> num);
  alloc2D (epsuv, sv -> natv, su -> num);

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int iu = 0; iu < su -> num; ++iu) {
      siguv[iv][iu] = (su -> sig[iu] + sv -> sigv[iv]) * 0.5;
      epsuv[iv][iu] = sqrt (su -> eps[iu] * sv -> epsv[iv]) * kcal2J;
    }
  }

  alloc2D (uuv, sv -> natv, ce -> mgrid);
    
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> mgrid; ++ig) {	
      uuv[iv][ig] = 0.0;
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int igz = ce -> zstart; igz < ce -> zend; ++igz) {
      double rz = (igz - ce -> grid[2] / 2) * ce -> dr[2];
      for (int igy = 0; igy < ce -> grid[1]; ++igy) {
	double ry = (igy - ce -> grid[1] / 2) * ce -> dr[1];
	for (int igx = 0; igx < ce -> grid[0]; ++igx) {
	  double rx = (igx - ce -> grid[0] / 2) * ce -> dr[0];
	  int ig = igx + igy * ce -> grid[0] 
	    + (igz - ce -> zstart) * ce -> grid[0] * ce -> grid[1];
	  for (int iu = 0; iu < su -> num; ++iu) {
	    int num = iu * 3;
	    double dx = rx - su -> r[num];
	    double dy = ry - su -> r[num + 1];
	    double dz = rz - su -> r[num + 2];

	    double r2 = dx * dx + dy * dy + dz * dz;
	    
	    if (r2 < cut2) r2 = cut2;

	    double irs2 = siguv[iv][iu] * siguv[iv][iu] / r2;
	    
	    double irs6 = irs2 * irs2 * irs2;
	    uuv[iv][ig] += epsuv[iv][iu] * 4.0 * irs6 * (irs6 - 1.0);
	  }
	}
      }
    }
  }

  double iKbT = 1.0 / (avogadoro * boltzmann * sv -> temper);
  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> mgrid; ++ig) {
      uuv[iv][ig] *= iKbT;
    }
  }
}
