#ifndef SOLUTE_H
#define SOLUTE_H

class Solute {
 public:
  Solute(){}
  ~Solute(){delete[] q, sig, eps, r;}
  void init(int);
  double * q;
  double * sig;
  double * eps;
  double * r;
  int num;
};

#endif
