#include "../utility.hpp"

double f(double* Js, double temperature) {
  return sinh(2.0 * Js[0] / temperature) * sinh(2.0 * Js[1] / temperature) -
         1.0;
}

// 0.567296328553

int main() {
  START();

  double Js[2] = {1.0 / 4.0, 1.0 / 4.0};
  double EPS = 1e-12;
  double l = 0.1, r = 0.6;

  double c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? l : r) = c;
  }
  printf("%.12f\n", c);

  END();
}