#ifndef SOLUTE_H
#define SOLUTE_H
#include <stdlib.h>

class Solute {
 public:
  Solute(){}
  ~Solute(){delete[] q, sig, eps, r;}
  void init(int);
  void setup_cuda();
  double * q;
  double * sig;
  double * eps;
  double * r;
  double * dq;
  double3 * dr;
  int num;
};

#endif
