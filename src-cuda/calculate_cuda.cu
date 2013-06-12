__global__ void kh(double * dtr, double * dt, double * du) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  double earg = - du[ip] + dt[ip];
  if (earg >= 0.0) {
    dtr[ip] = 1.0 + earg;
  } else {
    dtr[ip] = exp(earg);
  }
}

__global__ void hnc(double * dtr, double * dt, double * du) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  dtr[ip] = exp(- du[ip] + dt[ip]);
}

__global__ void trm1mt(double2 * dguv, double * dtr, double * dt,
		       double * dfr, double qv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  dguv[ip].x = dtr[ip] - 1.0 - dt[ip] + qv * dfr[ip];
  dguv[ip].y = 0.0;
}

__global__ void pqvfr(double2 * dguv, double * dfr, double qv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  dguv[ip].x += qv * dfr[ip];
}

__global__ void mqvfk(double2 * dguv, double2 * dfk, double qv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  dguv[ip].x -= qv * dfk[ip].x;
  dguv[ip].y -= qv * dfk[ip].y;
}

__global__ void oz(double2 * dhuv, double2 * dguv, double * dx, int natv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  unsigned int ngr = blockDim.x * gridDim.x * gridDim.y;
  double hr = 0.0;
  double hi = 0.0;
  for (unsigned int iv = 0; iv < natv; ++iv) {
    unsigned int i = ip + iv * ngr;
    hr += dguv[i].x * dx[i];
    hi += dguv[i].y * dx[i];
  }
  dhuv[ip].x = hr;
  dhuv[ip].y = hi;
}

__global__ void tr(double2 * dguv, double * dtr, double2 * dhuv) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x 
    + blockIdx.y * blockDim.x * gridDim.x;
  dguv[ip].x = dtr[ip];
  dtr[ip] = dhuv[ip].x + 1.0 - dguv[ip].x;
  //  dguv[ip].y = 0.0;
}
