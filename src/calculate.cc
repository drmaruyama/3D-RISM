#include "rism3d.h"

void RISM3D :: calculate (double cf) {

  if (clos == 0) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double q = sv -> qv[iv] * cf;
#pragma omp parallel for
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	double earg = - uuv[iv][ig] - euv[ig] * q + tuv[iv][ig];
	if (earg >= 0.0) {
	  tuvdif[iv][ig] = 1.0 + earg;
	} else {
	  tuvdif[iv][ig] = exp(earg);
	}	
      }
    }
  } else if (clos == 1) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      double q = sv -> qv[iv] * cf;
#pragma omp parallel for
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	tuvdif[iv][ig] = exp(- uuv[iv][ig] - euv[ig] * q + tuv[iv][ig]);
      }
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    double q = sv -> qv[iv] * cf;
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      guv[iv][ig] = tuvdif[iv][ig] - 1.0 - tuv[iv][ig] + q * fr[ig];
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(guv[iv], - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    double q = sv -> qv[iv] * cf;
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      guv[iv][ig] -= q * fk[ig];
      huv[iv][ig] = complex <double> (0.0, 0.0);
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      huv[iv][ig] = complex <double> (0.0, 0.0);
    }
  }

  for (int iv2 = 0; iv2 < sv -> natv; ++iv2) {	
    for (int iv1 = 0; iv1 < sv -> natv; ++iv1) {
#pragma omp parallel for
      for (int ig = 0; ig < ce -> ngrid; ++ig) {
	huv[iv1][ig] += guv[iv2][ig] * sv -> xvva[iv1][iv2][ig];
      }
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(huv[iv], 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      guv[iv][ig] = complex <double> (tuvdif[iv][ig], 0.0);
    }
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
#pragma omp parallel for
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      tuvdif[iv][ig] = huv[iv][ig].real() + 1.0 - guv[iv][ig].real();
    }
  }
} 
