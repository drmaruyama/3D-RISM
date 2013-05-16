#include "fft3d.h"

void FFT3D :: initialize (double * & box, int * & ng3) {
  fftw_init_threads();
  int omp_num;
#pragma omp parallel
  omp_num = omp_get_num_threads();

  fftw_plan_with_nthreads(omp_num);

  in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) 
				    * ng3[0] * ng3[1] * ng3[2]);
  out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) 
				     * ng3[0] * ng3[1] * ng3[2]);
  pf = fftw_plan_dft_3d(ng3[2], ng3[1], ng3[0], in, out, 
			FFTW_FORWARD, FFTW_MEASURE);
  pb = fftw_plan_dft_3d(ng3[2], ng3[1], ng3[0], in, out, 
			FFTW_BACKWARD, FFTW_MEASURE);

  double dkx = M_PI / ng3[0];
  double dky = M_PI / ng3[1];
  double dkz = M_PI / ng3[2];
  ngr = ng3[0] * ng3[1] * ng3[2];
  volf = box[0] * box[1] * box[2] / ngr;
  volb = 1.0 / (box[0] * box[1] * box[2]);
  kf = new complex<double>[ngr];
  kb = new complex<double>[ngr];
  ir = new int[ngr];
#pragma omp parallel for 
  for (int iz = 0; iz < ng3[2]; ++iz) {
    double lgz = iz - ng3[2] / 2;
    for (int iy = 0; iy < ng3[1]; ++iy) {
      double lgy = iy - ng3[1] / 2;
      for (int ix = 0; ix < ng3[0]; ++ix) {
	double lgx = ix - ng3[0] / 2;
	int ig = ix + iy * ng3[0] + iz * ng3[1] * ng3[0];
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
    for (int ig = 0; ig < ngr; ++ig) {
      complex <double> tmp = data[ig] * kf[ig];
      in[ig][0] = tmp.real() * ir[ig];
      in[ig][1] = tmp.imag() * ir[ig];
    }

    fftw_execute(pf);

#pragma omp parallel for
    for (int ig = 0; ig < ngr; ++ig) {
      data[ig] = complex <double> (out[ig][0], out[ig][1]);
      data[ig] *= volf * ir[ig];
    }
  } else {
#pragma omp parallel for
    for (int ig = 0; ig < ngr; ++ig) {
      in[ig][0] = data[ig].real() * ir[ig];
      in[ig][1] = data[ig].imag() * ir[ig];
    }

    fftw_execute(pb);

#pragma omp parallel for
    for (int ig = 0; ig < ngr; ++ig) {
      data[ig] = complex <double> (out[ig][0], out[ig][1]);
      data[ig] *= volb * ir[ig] * kb[ig];
    }
  }
}
