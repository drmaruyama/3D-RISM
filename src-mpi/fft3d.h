#ifndef FFT3D_H
#define FFT3D_H
#include <valarray>
#include <complex>
#include <omp.h>
#include <fftw3-mpi.h>
using namespace std;

class FFT3D {
 public:
  FFT3D () {}
  ~FFT3D () {fftw_destroy_plan(pf);
             fftw_destroy_plan(pb);
             fftw_free(in);
             fftw_free(out);
	     fftw_mpi_cleanup();}
  void initialize (double * &, int * &, int, int);
  void execute (complex <double> * &, int);
 private:
  fftw_complex *in, *out;
  fftw_plan pf, pb;
  complex <double> * kf;
  complex <double> * kb;
  int * ir;
  double volf, volb;
  int ngr;
  int mgr;
};

#endif 
