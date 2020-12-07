#include "../utility.hpp"

double f(double* Js, double temperature) {
  return tanh(Js[0] / temperature) * tanh(Js[1] / temperature) +
         tanh(Js[1] / temperature) * tanh(Js[2] / temperature) +
         tanh(Js[2] / temperature) * tanh(Js[0] / temperature) - 1.0;
}

// 0.506217145000

int main() {
  double Js[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  double EPS = 1e-12;
  double l = 0.1, r = 0.7;

  double c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? l : r) = c;
  }
  printf("%.12f\n", c);
}