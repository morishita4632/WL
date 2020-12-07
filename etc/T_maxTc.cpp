#include "../utility.h"

double f(double* Js, double temperature) {
  return exp(-2.0 * Js[0] / temperature) * exp(-2.0 * Js[1] / temperature) +
         exp(-2.0 * Js[1] / temperature) * exp(-2.0 * Js[2] / temperature) +
         exp(-2.0 * Js[2] / temperature) * exp(-2.0 * Js[0] / temperature) -
         1.0;
}

// 0.606826151085

int main() {
  double Js[3] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
  double EPS = 1e-12;
  double l = 0.1, r = 0.7;

  double c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? r : l) = c;
  }
  printf("%.12f\n", c);
}