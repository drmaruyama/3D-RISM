#include "cell.h"

void Cell :: setup() {
  volume = box[0] * box[1] * box[2];
  ngrid = grid[0] * grid[1] * grid[2];
  dv = volume / ngrid;
  dr[0] = box[0] / grid[0];
  dr[1] = box[1] / grid[1];
  dr[2] = box[2] / grid[2];
}
