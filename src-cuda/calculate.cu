#include <iostream>
#include "rism3d.h"

void RISM3D :: calculate () {
  __global__ void kh(double * dtr, double * dt, double * du);
  __global__ void hnc(double * dtr, double * dt, double * du);
  __global__ void trm1mt(double2 * dguv, double * dtr, double * dt,
                         double * dfr, double qv);
  __global__ void mqvfk(double2 * dguv, double2 * dfk, double qv);
  __global__ void oz(double2 * dhuv, double2 * dguv, double * dx, int natv);
  __global__ void tr(double2 * dguv, double * dtr, double2 * dhuv);

  int ng = ce -> ngrid;

  if (clos == 0) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      kh <<< g, b >>> (dtr + (iv * ng), dt + (iv * ng), du + (iv * ng));
    }
  } else if (clos == 1) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      hnc <<< g, b >>> (dtr + (iv * ng), dt + (iv * ng), du + (iv * ng));
    }
  } 

  for (int iv = 0; iv < sv -> natv; ++iv) {
    trm1mt <<< g, b >>> (dguv + (iv * ng), dtr + (iv * ng),
			  dt + (iv * ng), dfr, sv -> qv[iv]);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dguv + (iv * ng), - 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    mqvfk <<< g, b >>> (dguv + (iv * ng), dfk, sv -> qv[iv]);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    oz <<< g, b >>> (dhuv + (iv * ng), dguv,
		      sv -> dx + (iv * sv -> natv * ng), sv -> natv);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    fft -> execute(dhuv + (iv * ng), 1);
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    tr <<< g, b >>> (dguv + (iv * ng), dtr + (iv * ng), dhuv + (iv * ng));
  }
} 
