#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_solvent () {
  if (myrank == 0) sv -> read(fsolvent);
  sv -> spline(ga, indga, nga, ce -> ngrid);
}
