#include <iostream>
#include <fstream>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  RISM3D * system;
  int dn;

  system = new RISM3D;

  if (argc == 1) {
    cout << "No parameter file!" << endl ;
    return (1) ;
  }
  if (argc == 2) {
    cout << "Set device 0" << endl ;
    dn = 0;
  }
  if (argc == 3) {
    dn = atoi(argv[2]);
    cout << "Set device " << dn << endl ;
  }
  if (argc > 3) {
    cout << "Too much arguments!" << endl ;
    return (1) ;
  }

  cudaSetDevice(dn);
  system -> initialize(argv[1]);
  system -> iterate();
  system -> output();    

  return(0);
}
