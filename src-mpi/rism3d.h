#ifndef RISM3D_H
#define RISM3D_H

#include <complex>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <valarray>
#include <vector>
#include "physical.h"
#include "cell.h"
#include "control.h"
#include "solute.h"
#include "solvent.h"
#include "anderson.h"
#include "fft3d.h"
using namespace std;

class RISM3D {
public:
  RISM3D () {ce = new Cell; co = new Control; su = new Solute;
    sv = new Solvent; ma = new AN2; fft = new FFT3D;}
  ~RISM3D () {delete ce, co, su, sv;} 
  void initialize (char[]);
  void iterate ();    
  void output ();
private:
  void cal_Coulomb (double * &);
  valarray <double> cal_exchem ();
  void cal_grad (double * &);
  void cal_LJ ();
  double cal_pmv ();
  void cal_potential ();
  void cal_sum (double &, double &);
  void calculate ();
  void initialize_g ();
  void initialize_tuv ();
  void output_cuv ();
  void output_grad (double * &);
  void output_guv ();
  void output_huv ();
  void output_xmu (valarray <double> &, double);
  void read_input (char[]);
  void read_tuv ();
  void set_fname (char[]);
  void set_solvent ();
  void write_tuv ();

  vector <complex <double> *> guv;
  vector <complex <double> *> huv;
  vector <double *> tuv;
  vector <double *> tuvdif;
  vector <double *> gv;
  vector <double *> uuv;
  vector <double *> siguv;
  vector <double *> epsuv;
  double * fr;
  complex <double> * fk;
  vector <double> ga;
  int * indga;

  string outlist;
  string closure;
  string fsolvent;
  string fname;
  int nga;

  Cell * ce;
  Control * co;
  Solute * su;
  Solvent * sv;
  AN2 * ma;
  FFT3D * fft;
};

#endif
