#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: iterate() {
  void alloc2D (vector <double *> &, int, int);
  void calloc2D (vector <complex <double> *> &, int, int);

  calloc2D (guv, sv -> natv, ce -> ngrid);
  calloc2D (huv, sv -> natv, ce -> ngrid);
  alloc2D (tuv, sv -> natv, ce -> ngrid);

  cudaMalloc(&dguv, ce -> ngrid * sv -> natv * sizeof(double2));
  cudaMalloc(&dhuv, ce -> ngrid * sv -> natv * sizeof(double2));
  cudaMalloc(&dt, ce -> ngrid * sv -> natv * sizeof(double));
  cudaMalloc(&dtr, ce -> ngrid * sv -> natv * sizeof(double));
  cudaMalloc(&ds, ce -> grid[1] * ce -> grid[2] * sizeof(double));

  ifstream in_file ;
  in_file.open((fname + exttuv).c_str());
  bool saved = in_file.is_open();
  in_file.close();

  if (saved) {
    read_tuv();
  } else {
    initialize_tuv();
  }

  ma -> initialize (ce, sv);
  fft -> initialize (ce);

  cout << "relaxing 3D UV RISM:" << endl;
  bool conver = false;
  for (int istep = 1; istep <= co -> maxstep; ++istep) {
    calculate();
    double rms = cal_rms ();
    if (rms <= co -> convergence) {
      conver = true;
    } else {
      ma -> calculate (dt, dtr);
    }
    cout << " Step = " << istep << " Reside = " << rms << endl;
    if (co -> ksave > 0 && istep % co -> ksave == 0) {
      write_tuv();
    }
    if (conver) {
      if (co -> ksave != 0) {
	write_tuv();
      }
      break;
    }
  }
  if (!conver) {
    cout << "3D UV RISM: reached limit # of relaxation steps: "
	 << co -> maxstep << endl;
  }
  for (int iv = 0; iv < sv -> natv; ++iv) {
    cudaMemcpyAsync(huv[iv], dhuv + (iv * ce -> ngrid), 
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault);
    cudaMemcpyAsync(guv[iv], dguv + (iv * ce -> ngrid), 
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault);
  }
  delete ma;
  delete fft;
} 
