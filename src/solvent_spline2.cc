#include <iostream>
#include <fstream>
#include <string>
#include "alloc.h"
#include "solvent.h"

void Solvent :: spline2 (vector <double> & ga, int * & indga,
			int nga, int ngrid) {
  void spline (double * &, double * &, int, vector <double *> &);
  double splint (double * &, double * &, vector <double *> &, int, double);

  if (ga[nga - 1] > ttab2[ntab - 1]) {
    cout << "insufficient maximal T tabulated" << endl;
    exit (1);
  }

  vector <double *> chsa;
  alloc2D(chsa, natv, nga);
  alloc2D(wfka, natv, nga);

  int ntab1 = 0;
  for (int n = 0; n < ntab; ++n) {
    if (ttab2[n] < ga[0]) 
      ntab1 = n;
  }

  int ntab2 = ntab - 1;
  for (int n = ntab - 1; n > 0; --n) {
    if (ttab2[n] > ga[nga - 1]) 
      ntab2 = n;
  }

  int np = ntab2 - ntab1 + 1;

  double * x = new double[np];
  double * y = new double[np];
  vector <double *> coe;
  alloc2D(coe, 3, np);

  for (int n = 0; n < np; ++n) {
    x[n] = ttab2[ntab1 + n] ;
  }

  for (int iv = 0; iv < natv; ++iv) {
#pragma omp parallel for
    for (int n = 0; n < np; ++n) {
      y[n] = chs[iv][ntab1 + n] ;
    }
    spline(x, y, np, coe) ;
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      chsa[iv][i] = splint(x, y, coe, np, ga[i]);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
#pragma omp parallel for
    for (int n = 0; n < np; ++n) {
      y[n] = wfk[iv][ntab1 + n] ;
    }
    spline(x, y, np, coe) ;
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      wfka[iv][i] = splint(x, y, coe, np, ga[i]);
    }
  }

  for (int iv = 0; iv < natv; ++iv) {
#pragma omp parallel for
    for (int i = 0; i < nga; ++i) {
      cvva[iv][iv][i] -= chsa[iv][i] * rhov[0] / rhov[iv];
    }
  }

  dealloc2D(coe);
  dealloc2D(chsa);
  delete[] x;
  delete[] y;
  delete[] ttab2;
}
