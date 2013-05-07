#include "anderson.h"

AN2 :: ~AN2 () {
  void dealloc2D (vector <double *> &);
  delete[] a, c, x;
  dealloc2D (tp);
  dealloc2D (rp);
}


void AN2 :: setup_mpi() {
  MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void AN2 :: initialize (double * & rhov, int r, int v) {
  void alloc2D (vector <double *> &, int, int);

  ngrid = r ;
  niv = v ;
  binary = 0;
  alloc2D(tp, niv * 2, ngrid);
  alloc2D(rp, niv * 2, ngrid);
  a = new double[4];
  c = new double[2];
  x = new double[2];
  irho = new double[niv];
  for (int iv = 0; iv < niv; ++iv) {
    irho[iv] = 1.0 / rhov[iv];
  }
}


void AN2 :: calculate (vector <double *> & t, vector <double *> & r) {
  void leg(double * a, double * x, double * b, int n, int nn);

  if (count > 0) {
    --count;
    newt0(t, r);
  } else {
    cal_theta(t, r);
    a[2] = a[1];
    leg(a, x, c, 2, 2);
    s1 = x[0] + mp * ((binary == 0)? 0 : 1);
    s2 = x[1] + mp * ((binary == 0)? 1 : 0);
    newt(t, r);
  }
  binary = (binary == 0)? 1 : 0;
}


void AN2 :: cal_theta (vector <double *> & t, vector <double *> & r) {
  double s0, s1, s2, s3, s4;
  double ss0, ss1, ss2, ss3, ss4;
  c[0] = c[1] = a[0] = a[1] = a[3] = 0.0;
  for (int iv = 0; iv < niv; ++iv) {
    s0 = s1 = s2 = s3 = s4 = 0.0;
#pragma omp parallel for reduction(+: s0, s1, s2, s3, s4)
    for (int ig = 0; ig < ngrid; ++ig) {
      double t1 = r[iv][ig] - rp[iv][ig];
      double t2 = r[iv][ig] - rp[iv + niv][ig];
      s0 += r[iv][ig] * t1;
      s1 += r[iv][ig] * t2;
      s2 += t1 * t1;
      s3 += t1 * t2;
      s4 += t2 * t2;
    }
    MPI_Allreduce(&s0, &ss0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&s1, &ss1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&s2, &ss2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&s3, &ss3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&s4, &ss4, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    c[0] += ss0 * irho[iv];
    c[1] += ss1 * irho[iv];
    a[0] += ss2 * irho[iv];
    a[1] += ss3 * irho[iv];
    a[3] += ss4 * irho[iv];
  }
}


void AN2 :: newt (vector <double *> & t, vector <double *> & r) {
  for (int iv = 0; iv < niv; ++iv) {
    int iv2 = iv + niv;
    int biv = iv + binary * niv;
#pragma omp parallel for 
    for (int ig = 0; ig < ngrid; ++ig) {
      double u = t[iv][ig] + s1 * (tp[iv][ig] - t[iv][ig]) 
	+ s2 * (tp[iv2][ig] - t[iv][ig]);
      double v = t[iv][ig] + r[iv][ig] 
	+ s1 * (tp[iv][ig] + rp[iv][ig]- t[iv][ig] - r[iv][ig])
	+ s2 * (tp[iv2][ig] + rp[iv2][ig]- t[iv][ig] - r[iv][ig]);
      tp[biv][ig] = t[iv][ig];
      rp[biv][ig] = r[iv][ig];
      t[iv][ig] = u + m * (v - u);
    }
  }
}

void AN2 :: newt0 (vector <double *> & t, vector <double *> & r) {
  for (int iv = 0; iv < niv; ++iv) {
    int biv = iv + binary * niv;
#pragma omp parallel for
    for (int ig = 0; ig < ngrid; ++ig) {
      tp[biv][ig] = t[iv][ig];
      rp[biv][ig] = r[iv][ig];
      t[iv][ig] += r[iv][ig];
    }
  }
}

void leg(double * a, double * x, double * b, int n, int nn) {
  int i, j, k, k1;
  double s, p, q, r;
  double *ai, *ak;
  for (k = 0, ak = a; k < n -1; ++k, ak += nn) {
    k1 = k + 1;
    p = ak[k];
    for (j = k1; j < n; ++j)
      ak[j] /= p;
    r = b[k] /= p;
    for (i = k1, ai = ak + nn; i < n; ++i, ai += nn) {
      q = ai[k];
      for (j = k1; j < n; ++j)
        ai[j] -= q * ak[j];
      b[i] -= q * r;
    }
  }
  x[n - 1] = b[n - 1] / ak[n - 1];
  for (k = n - 2, ak = a + nn * (n - 2); k >= 0; --k, ak -= nn) {
    k1 = k + 1;
    s = b[k];
    for (j = k1; j < n; ++j)
      s -= ak[j] * x[j];
    x[k] = s;
  }
}
