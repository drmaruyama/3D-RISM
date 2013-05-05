#include "cell.h"

void Cell :: setup() {
  MPI_Bcast(box, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(grid, 3, MPI_INT, 0, MPI_COMM_WORLD);

  volume = box[0] * box[1] * box[2];
  ngrid = grid[0] * grid[1] * grid[2];
  dv = volume / ngrid;
  dr[0] = box[0] / grid[0];
  dr[1] = box[1] / grid[1];
  dr[2] = box[2] / grid[2];
  setup_mpi();
}

void Cell :: setup_mpi() {
  ptrdiff_t n0, n0s;

  mgrid = (int) fftw_mpi_local_size_3d(grid[2], grid[1], grid[0],
				       MPI_COMM_WORLD, & n0, & n0s);
  zstart = (int) n0s;
  zend   = zstart + (int) n0;
}
