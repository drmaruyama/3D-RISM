#ifndef AN2_H
#define AN2_H

#include <valarray>
#include <vector>
#include <mpi.h>
using namespace std;

class AN2 {
public:
  AN2 () {}
  ~AN2 ();
  void initialize (double * &, int, int);
  void calculate (vector <double *> &, vector <double *> &);
  void setup_mpi();
  double m;
  double mp;
  int count;
private:
  void cal_theta (vector <double *> &, vector <double *> &);
  void newt (vector <double *> &, vector <double *> &);
  void newt0 (vector <double *> &, vector <double *> &);
  vector <double *> tp;
  vector <double *> rp;
  double * irho;
  double * a;
  double * c;
  double * x;
  double s0;
  double s1;
  double s2;
  double theta;
  int niv;
  int ngrid;
  int binary;
};

#endif // AN2_H
