#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_xmu(valarray <double> & xmu, double dft, double pmv) {
    
  ofstream out_file;
  out_file.open((fname + extxmu).c_str());

  double ibeta = avogadoro * boltzmann * sv -> temper / kcal2J;

  out_file << "SFE (DFT) = " << ibeta * dft
	   << " (kcal/mol)" << endl;
  out_file << endl;

  out_file << "SFE (SC) = " << ibeta * xmu.sum() 
	   << " (kcal/mol)" << endl;
  for (int iv = 0; iv < xmu.size(); ++iv) {
    out_file << "  " << iv << " -----> " << ibeta * xmu[iv] << endl;
  }
  out_file << endl;

  out_file << "PMV  = " << pmv << " (cc/mol)" << endl;

  out_file.close();
} 
