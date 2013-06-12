#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: iterate() {
  void alloc2D (vector <double *> &, int, int);
  void calloc2D (vector <complex <double> *> &, int, int);

  calloc2D (guv, sv -> natv, ce -> mgrid);
  calloc2D (huv, sv -> natv, ce -> mgrid);
  alloc2D (tuv, sv -> natv, ce -> mgrid);
  alloc2D (tuvdif, sv -> natv, ce -> mgrid);

  ifstream in_file ;
  in_file.open((fname + exttuv).c_str());
  bool saved = in_file.is_open();
  in_file.close();

  if (saved) {
    read_tuv();
  } else {
    initialize_tuv();
  }
  ma -> initialize (sv -> rhov, ce -> mgrid, sv -> natv);
  fft -> initialize (ce -> box, ce -> grid, ce -> zstart, ce -> zend);

  if (myrank == 0) cout << "relaxing 3D UV RISM:" << endl;
  bool conver = false;
  for (int istep = 1; istep <= co -> maxstep; ++istep) {
    calculate();
    double rms = cal_rms ();
    if (rms <= co -> convergence) {
      conver = true;
    } else {
      ma -> calculate (tuv, tuvdif);
    }
    if (myrank == 0) {
      cout << " Step = " << istep << " Reside = " << rms << endl;
    }
    if (co -> ksave > 0 && istep % co -> ksave == 0) {
      write_tuv();
    }
    if (conver) {
      if (co -> ksave != 0) write_tuv();
      break;
    }
  }
  if (!conver) {
    if (myrank == 0) {
      cout << "3D UV RISM: reached limit # of relaxation steps: "
	   << co -> maxstep << endl;
    }
  }
  delete ma;
  delete fft;
} 
