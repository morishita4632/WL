#include <math.h>
#include <chrono>
#include <iostream>
#include <random>
using namespace std;

chrono::system_clock::time_point chrono_start, chrono_end;
mt19937 mt(0);

static inline double* alloc_dvector(int n) {
  double* vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of int */
static inline int* alloc_ivector(int n) {
  int* vec;
  vec = (int*)calloc(n, sizeof(int));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_ivector\n");
    exit(1);
  }
  return vec;
}

/* allocate m x n row-major matrix of double */
static inline double** alloc_dmatrix(int m, int n) {
  int i;
  double** mat;
  mat = (double**)malloc((size_t)(m * sizeof(double*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  mat[0] = (double*)calloc(m * n, sizeof(double));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  for (i = 1; i < m; ++i)
    mat[i] = mat[i - 1] + n;
  return mat;
}

static inline double rand01() {
  return (double)(mt()) / (double)(mt.max());
}

static inline void vzero(double* v, int dim) {
  for (int i = 0; i < dim; i++) {
    v[i] = 0.0;
  }
}

// Measure time & Random initialize
static inline void START(int seed = 0) {
  mt.seed(seed);
  chrono_start = chrono::system_clock::now();
}

static inline void END() {
  chrono_end = chrono::system_clock::now();
  double time = static_cast<double>(
      chrono::duration_cast<chrono::microseconds>(chrono_end - chrono_start)
          .count() /
      1000000.0);
  printf("time %.2lf[s],  %.2lf[m]\n", time, time / 60.0);
}