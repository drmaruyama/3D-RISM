#ifndef Cell_H
#define Cell_H
#include <fftw3-mpi.h>
using namespace std;

class Cell {
public:
  Cell() {box = new double[3]; dr = new double[3]; grid = new int[3];}
  ~Cell() {delete[] box, dr, grid;}
  void setup();
  void setup_mpi();
  double * box;
  double * dr;
  int * grid;
  double volume, dv;
  int ngrid;
  int mgrid, zstart, zend;
};

#endif // Cell_H
