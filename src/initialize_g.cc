#include <iostream>
#include <fstream>
#include "rism3d.h"

using namespace std;

void RISM3D :: initialize_g() {
  void alloc2D (vector <double *> &, int, int);
  void index(double * &, int * &, int);

  int ngx = ce -> grid[0];
  int ngy = ce -> grid[1];
  int ngz = ce -> grid[2];

  alloc2D(gv, ce -> ngrid, 3);
  indga = new int[ce -> ngrid];
  double * g2 = new double[ce -> ngrid];
  int * indg2 = new int[ce -> ngrid];


#pragma omp parallel for
  for (int igz = 0; igz < ngz; ++igz) {
    double lgz = igz - ngz / 2 + 0.5;
    for (int igy = 0; igy < ngy; ++igy) {
      double lgy = igy - ngy / 2 + 0.5;
      for (int igx = 0; igx < ngx; ++igx) {
	double lgx = igx - ngx / 2 + 0.5;
	int igk = igx + igy * ngx + igz * ngy * ngx;

	double gx = 2.0 * M_PI / ce -> box[0] * lgx;
	double gy = 2.0 * M_PI / ce -> box[1] * lgy;
	double gz = 2.0 * M_PI / ce -> box[2] * lgz;
	gv[igk][0] = gx;
	gv[igk][1] = gy;
	gv[igk][2] = gz;
	g2[igk] = gx * gx + gy * gy + gz * gz;
      }
    }
  }

  index(g2, indg2, ce -> ngrid);

  double ga2o = - 1.0;
  nga = 0;

  for (int igk = 0; igk < ce -> ngrid; ++igk) {
    int igs = indg2[igk];
    double ga2 = g2[igs];
    if (ga2 > ga2o) {
      ++nga;
      ga . push_back (sqrt(ga2));
      ga2o = ga2;
    }
    indga[igs] = nga - 1;
  }
  delete[] g2;
  delete[] indg2;
} 
