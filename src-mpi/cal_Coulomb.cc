#include <iostream>
#include "rism3d.h"

void RISM3D :: cal_Coulomb (double * & euv) {

  euv = new double[ce -> mgrid];
  fr = new double[ce -> mgrid];
  fk = new complex <double>[ce -> mgrid];

  if (myrank == 0) {
    cout << "synthesizing solute Coulomb potential ..." << endl;
  }

#pragma omp parallel for
  for (int ig = 0; ig < ce -> mgrid; ++ig){
    euv[ig] = 0.0;
    fr[ig] = 0.0;
    fk[ig] = complex<double>(0.0, 0.0);
  }

#pragma omp parallel for
  for (int igz = ce -> zstart; igz < ce -> zend; ++igz) {
    double rz = (igz - ce -> grid[2] / 2) * ce -> dr[2];
    for (int igy = 0; igy < ce -> grid[1]; ++igy) {
      double ry = (igy - ce -> grid[1] / 2) * ce -> dr[1];
      for (int igx = 0; igx < ce -> grid[0]; ++igx) {
	double rx = (igx - ce -> grid[0] / 2) * ce -> dr[0];
	int ig = igx + igy * ce -> grid[0] 
	  + (igz - ce -> zstart) * ce -> grid[0] * ce -> grid[1];
	double rk2 = gv[ig][0] * gv[ig][0] + gv[ig][1] * gv[ig][1]
	  + gv[ig][2] * gv[ig][2];
	double irk4 = 1.0 / (rk2 * (rk2 + 1.0));

	for (int iu = 0; iu < su -> num; ++iu) {
	  int num = iu * 3;
	  double delx = rx - su -> r[num];
	  double dely = ry - su -> r[num + 1];
	  double delz = rz - su -> r[num + 2];
	  double ra = sqrt(delx * delx + dely * dely + delz * delz);
	  double ruk = gv[ig][0] * su -> r[num] + gv[ig][1] * su -> r[num + 1]
	    + gv[ig][2] * su -> r[num + 2];
	  if (ra >= 1.0e-5) {
	    euv[ig] += su -> q[iu] / ra;
	    fr[ig] += su -> q[iu] / ra * (1 - exp(- ra));
	  } else {
	    fr[ig] += su -> q[iu];
	  }
	  fk[ig] += 
	    complex<double> (4.0 * M_PI * su -> q[iu] * irk4, 0.0) 
	    * exp(complex<double> (0.0, - ruk));
	}
      }
    }
  }

  double ubeta = hartree * bohr / (boltzmann * sv -> temper);
  complex <double> cubeta = complex <double> (ubeta, 0.0);

#pragma omp parallel for
  for (int ig = 0; ig < ce -> mgrid; ++ig){
    euv[ig] = euv[ig] * ubeta;
    fr[ig] = fr[ig] * ubeta;
    fk[ig] = fk[ig] * cubeta;
  }
} 
