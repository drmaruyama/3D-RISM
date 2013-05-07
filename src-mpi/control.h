#ifndef Control_H
#define Control_H
#include <mpi.h>

class Control {
public:
  Control() {}
  ~Control() {}
  void setup_mpi();
  double convergence;
  int maxstep;
  int ksave;
};

#endif // Control_H
