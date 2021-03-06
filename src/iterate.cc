#include <iostream>
#include <fstream>
#include "alloc.h"
#include "rism3d.h"
#include "extension.h"

void RISM3D :: iterate(int cu) {
  double cf, cuf;

  alloc2D (guv, sv -> natv, ce -> ngrid);
  alloc2D (huv, sv -> natv, ce -> ngrid);
  alloc2D (tuv, sv -> natv, ce -> ngrid);
  alloc2D (tuvdif, sv -> natv, ce -> ngrid);

  ma -> initialize (sv -> rhov, ce -> ngrid, sv -> natv);
  fft -> initialize (ce -> box, ce -> grid);

  ifstream in_file ;
  in_file.open((fname + exttuv).c_str());
  bool saved = in_file.is_open();
  in_file.close();

  if (saved) {
    read_tuv();
    cu = 0;
    cf = 1.0;
  } else if (cu == 0) {
    cf = 1.0;
    initialize_tuv(cf);
  } else {
    cuf = 1.0 / cu;
    cf = 0.0;
    initialize_tuv(cf);
  }

  for (int c = 0; c <= cu; ++c) {
    if (c > 0) {
      cf += cuf;
      if (cf > 1.0) cf = 1.0;
      add_tuv(cuf);
    }
    cout << "relaxing 3D UV RISM: Charge Up Factor = " << cf << endl;
    bool conver = false;
    for (int istep = 1; istep <= co -> maxstep; ++istep) {
      calculate(cf);
      double rms = cal_rms ();
      if (rms <= co -> convergence) {
	conver = true;
      } else {
	ma -> calculate (tuv, tuvdif);
      }
      cout << " Step = " << istep << " Reside = " << rms << endl;
      if (co -> ksave > 0 && istep % co -> ksave == 0 && c == cu) {
	write_tuv();
      }
      if (conver) {
	if (co -> ksave != 0 && c == cu) write_tuv();
	break;
      }
    }
    if (!conver) {
      cout << "3D UV RISM: reached limit # of relaxation steps: "
	   << co -> maxstep << endl ;
    }
  }
} 
