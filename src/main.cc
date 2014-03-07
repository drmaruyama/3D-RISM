#include <iostream>
#include <fstream>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  RISM3D * system;
  int ch;
  int cu = 0;

  system = new RISM3D;

  while ((ch = getopt(argc, argv, "c:")) != -1) {
    switch (ch) {
    case 'c':
      cu = atoi(optarg);
      break;
    }
  }

  if (argc == 1) {
    cout << "No parameter file!" << endl;
    return (1);
  } 

  if (cu > 0) cout << "Charge up " << cu << endl;
  system -> initialize(argv[optind]);
  system -> iterate(cu);
  system -> output();    

  return(0);
}
