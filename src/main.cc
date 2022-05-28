#include <iostream>
#include <fstream>
#include <unistd.h>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  RISM3D * system;
  int ch;
  int cu = 0;
  string input;
  string structure;

  system = new RISM3D;

  while ((ch = getopt(argc, argv, "c:i:s:")) != -1) {
    switch (ch) {
    case 'c':
      cu = atoi(optarg);
      break;
    case 'i':
      input = optarg;
      break;
    case 's':
      structure = optarg;
      break;
    }
  }

  if (input.empty() || structure.empty()) {
    if (argv[optind] == NULL) {
      cout << "No input file!" << endl;
      return (1);
    }
    input = argv[optind];
  }

  cout << "Input     : " << input << endl;
  if (!structure.empty()) {
    cout << "Structure : " << structure << endl;
  }

  if (cu > 0) cout << "Charge up " << cu << endl;
  system -> initialize(input, structure);
  system -> iterate(cu);
  system -> output();    

  return(0);
}
