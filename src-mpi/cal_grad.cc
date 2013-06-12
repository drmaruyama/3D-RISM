#include "rism3d.h"

void RISM3D :: cal_grad(double * & du) {

  const double cc = hartree * bohr * avogadoro;

  double * du2;
  du2 = new double[su -> num * 3];

#pragma omp parallel for
  for (int iu = 0; iu < su -> num * 3; ++iu) {
    du2[iu] = 0.0;
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int iz = ce -> zstart; iz < ce -> zend; ++iz) { 
      double z = (iz - ce -> grid[2] / 2) * ce -> dr[2];
      for (int iy = 0; iy < ce -> grid[1]; ++iy) { 
	double y = (iy - ce -> grid[1] / 2) * ce -> dr[1];
	for (int ix = 0; ix < ce -> grid[0]; ++ix) { 
	  double x = (ix - ce -> grid[0] / 2) * ce -> dr[0];
	  int ig = ix + iy * ce -> grid[0] 
	    + (iz - ce -> zstart) * ce -> grid[0] * ce -> grid[1];
#pragma omp parallel for 
	  for (int iu = 0; iu < su -> num; ++iu) {
	    int num = iu * 3;
	    double dx = x - su -> r[num];
	    double dy = y - su -> r[num + 1];
	    double dz = z - su -> r[num + 2];

	    double r2 = dx * dx + dy * dy + dz * dz;
	    double r1 = sqrt(r2);

	    if (r1 >= siguv[iv][iu] * 0.5) {
	      double rs2i = (siguv[iv][iu] * siguv[iv][iu]) / r2;
	      double rs6i = rs2i * rs2i * rs2i;
	      double gr = guv[iv][ig].real() * sv -> rhov[iv];
	      double ulj = epsuv[iv][iu] * 24.0 * rs6i / r2 
		* (2.0 * rs6i - 1.0) * gr;
	      double uco = su -> q[iu] * sv -> qv[iv] / (r2 * r1) * cc * gr;
	      du2[iu * 3] += (ulj + uco) * dx;
	      du2[iu * 3 + 1] += (ulj + uco) * dy;
	      du2[iu * 3 + 2] += (ulj + uco) * dz;
	    }
	  }
	}
      }
    }
  }
  MPI_Reduce(du2, du, su -> num * 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  delete[] du2;
}
