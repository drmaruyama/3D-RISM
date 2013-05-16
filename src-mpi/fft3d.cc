#include <iostream>
#include "fft3d.h"

void FFT3D :: initialize (double * & box, int * & ng3, int zs, int ze) {
  fftw_init_threads();
  fftw_mpi_init();
  int omp_num;
#pragma omp parallel
  omp_num = omp_get_num_threads();
  fftw_plan_with_nthreads(omp_num);

  double dkx = M_PI / ng3[0];
  double dky = M_PI / ng3[1];
  double dkz = M_PI / ng3[2];
  ngr = ng3[0] * ng3[1] * ng3[2];
  mgr = ng3[0] * ng3[1] * (ze - zs);
  volf = box[0] * box[1] * box[2] / ngr;
  volb = 1.0 / (box[0] * box[1] * box[2]);
  kf = new complex<double>[mgr];
  kb = new complex<double>[mgr];
  ir = new int[mgr];
  in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * mgr);
  out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * mgr);
  pf = fftw_mpi_plan_dft_3d(ng3[2], ng3[1], ng3[0], in, out, MPI_COMM_WORLD,
			    FFTW_FORWARD, FFTW_MEASURE);
  pb = fftw_mpi_plan_dft_3d(ng3[2], ng3[1], ng3[0], in, out, MPI_COMM_WORLD,
			    FFTW_BACKWARD, FFTW_MEASURE);
#pragma omp parallel for 
  for (int iz = zs; iz < ze; ++iz) {
    double lgz = iz - ng3[2] / 2;
    for (int iy = 0; iy < ng3[1]; ++iy) {
      double lgy = iy - ng3[1] / 2;
      for (int ix = 0; ix < ng3[0]; ++ix) {
	double lgx = ix - ng3[0] / 2;
	int ig = ix + iy * ng3[0] + (iz - zs) * ng3[1] * ng3[0];
	double dkr = dkx * lgx + dky * lgy + dkz * lgz;
	complex<double> cdkr = complex <double> (0.0, dkr);
	kf[ig] = exp(- cdkr);
	kb[ig] = exp(cdkr);
	if ((ix + iy + iz) % 2 == 0) {
	  ir[ig] = 1;
	} else {
	  ir[ig] = -1;
	}
      }
    }
  }
}

void FFT3D :: execute (complex <double> * & data, int key) {
  if (key == - 1) {
#pragma omp parallel for
    for (int ig = 0; ig < mgr; ++ig) {
      complex <double> tmp = data[ig] * kf[ig];
      in[ig][0] = tmp.real() * ir[ig];
      in[ig][1] = tmp.imag() * ir[ig];
    }

    fftw_execute(pf);

#pragma omp parallel for
    for (int ig = 0; ig < mgr; ++ig) {
      data[ig] = complex <double> (out[ig][0], out[ig][1]);
      data[ig] *= volf * ir[ig];
    }
  } else {
#pragma omp parallel for
    for (int ig = 0; ig < mgr; ++ig) {
      in[ig][0] = data[ig].real() * ir[ig];
      in[ig][1] = data[ig].imag() * ir[ig];
    }

    fftw_execute(pb);

#pragma omp parallel for
    for (int ig = 0; ig < mgr; ++ig) {
      data[ig] = complex <double> (out[ig][0], out[ig][1]);
      data[ig] *= volb * ir[ig] * kb[ig];
    }
  }
}
