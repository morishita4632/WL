#include "../utility.hpp"

static inline double f(const double* Js, double temperature) {
  return sinh(2.0 * Js[0] / temperature) * sinh(2.0 * Js[1] / temperature) -
         1.0;
}

static inline void normalize(double* Js) {
  double sum = 0.0;
  for (int i = 0; i < 2; i++)
    sum += Js[i];
  for (int i = 0; i < 2; i++)
    Js[i] /= sum * 2.0;
}

static inline double Js_to_Tc(const double* Js, double EPS) {
  double l = 0.01, r = 0.567296328554, c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? l : r) = c;
  }
  return (l + r) / 2.0;
}

static inline int Tc_to_ind(double Tc, double Tc_min, double bin_width,
                            double bins) {
  int ind = (int)((Tc - Tc_min) / bin_width);
  if (ind < 0 || ind >= bins)
    ind = -1;
  return ind;
}

static inline bool is_flat(const int* hist, int bins, double flat_coeff) {
  double ave = 0.0;
  for (int i = 0; i < bins; i++)
    ave += (double)hist[i] / (double)bins;
  bool flat = true;
  for (int i = 0; i < bins; i++)
    flat &= (hist[i] > ave * flat_coeff);
  return flat;
}

static inline double exp1(double cS, double nS) {
  return (cS - nS > 0) ? 1.0 : exp(cS - nS);
}