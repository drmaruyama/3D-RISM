#ifndef SOLUTE_H
#define SOLUTE_H
#include <stdlib.h>
#include <mpi.h>

class Solute {
 public:
  Solute(){}
  ~Solute(){delete[] q, sig, eps, r;}
  void init(int);
  void setup_mpi();
  double * q;
  double * sig;
  double * eps;
  double * r;
  int num;
};

#endif
