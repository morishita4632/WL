#include "utility.hpp"

double EPS = 1e-12;

double f(const double* Js, const double temperature) {
  return sinh(2.0 * Js[0] / temperature) * sinh(2.0 * Js[1] / temperature) -
         1.0;
}

void normalize(double* Js) {
  double sum = 0.0;
  for (int i = 0; i < 3; i++)
    sum += Js[i];
  for (int i = 0; i < 3; i++)
    Js[i] /= sum * 2.0;
}

double calc_Tc(const double* Js) {
  double l = 0.1, r = 0.567296328553, c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? l : r) = c;
  }
  return (l + r) / 2.0;
}

int main() {
  int bins = 30;
  double Tc_min = 0.1, Tc_max = 0.567296328553;
}