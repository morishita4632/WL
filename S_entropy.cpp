#include <vector>
#include "utility.hpp"

double EPS = 1e-12;

double f(const double* Js, const double temperature) {
  return sinh(2.0 * Js[0] / temperature) * sinh(2.0 * Js[1] / temperature) -
         1.0;
}

void normalize(double* Js) {
  double sum = 0.0;
  for (int i = 0; i < 2; i++)
    sum += Js[i];
  for (int i = 0; i < 2; i++)
    Js[i] /= sum * 2.0;
}

void k_to_Js(double k, double* Js) {
  Js[0] = 1.0, Js[1] = exp(k);
  normalize(Js);
}

double Js_to_Tc(const double* Js) {
  double l = 0.01, r = 0.567296328553, c;
  while (r - l > EPS) {
    c = (r + l) / 2.0;
    (f(Js, c) > 0 ? l : r) = c;
  }
  return (l + r) / 2.0;
}

int Tc_to_ind(double Tc, double Tc_min, double bin_width, double bins) {
  int ind = (int)((Tc - Tc_min) / bin_width + 0.5);
  if (ind < 0 || ind >= bins)
    ind = -1;
  return ind;
}

bool is_flat(const int* hist, int bins, double flat_coeff) {
  double ave = 0.0;
  for (int i = 0; i < bins; i++)
    ave += (double)hist[i] / (double)bins;
  bool flat = true;
  for (int i = 0; i < bins; i++)
    flat &= (hist[i] > ave * flat_coeff);
  return flat;
}

double exp1(double cS, double nS) {
  return (cS - nS > 0) ? 1.0 : exp(cS - nS);
}

int main() {
  START(1);

  int bins = 100;
  double Tc_min = 0.1, Tc_max = 0.567296328553;
  double f_min = 1e-9;
  double k_init = 1.0;
  double flat_coeff = 0.8;
  double RW_step = 1.0;

  int* hist = alloc_ivector(bins);
  double* S_s = alloc_dvector(bins);
  vector<double> Tc_s, k_s;
  vector<vector<double>> Js_s(0, vector<double>(2));

  double bin_width = (Tc_max - Tc_min) / (double)bins;
  double ck = k_init, dk, nk;
  double cJs[2], nJs[2];
  k_to_Js(ck, cJs);
  double cTc = Js_to_Tc(cJs), nTc;
  int cind = Tc_to_ind(cTc, Tc_min, bin_width, bins), nind;
  double f = 1.0;
  for (int i = 0; i < bins; i++) {
    hist[i] = 0, S_s[i] = 0.0;
  }

  int CNT = 0;
  int MAX_CNT = 100000000;
  while (f > f_min) {
    do {
      CNT++;
      dk = (rand01() * 2 - 1) * RW_step;
      nk = ck + dk;
      k_to_Js(nk, nJs);
      nTc = Js_to_Tc(nJs);
      nind = Tc_to_ind(nTc, Tc_min, bin_width, bins);
      if (nind != -1 && rand01() <= exp1(S_s[cind], S_s[nind])) {
        ck = nk, cTc = nTc, cind = nind;
        for (int i = 0; i < 2; i++)
          cJs[i] = nJs[i];
      }
      Tc_s.push_back(cTc), Js_s.push_back({cJs[0], cJs[1]});
      k_s.push_back(ck);
      hist[cind]++, S_s[cind] += f;
      if (CNT >= MAX_CNT)
        break;
    } while (!is_flat(hist, bins, flat_coeff));
    if (CNT >= MAX_CNT)
      break;
    for (int i = 0; i < bins; i++)
      hist[i] = 0;
    f /= 2.0;
    printf("%.12f\n", f);
  }


  double x;
  int y;
  FILE* gp;

  if (f > f_min) {
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set style fill solid border lc rgb \"black\"\n");
    fprintf(gp, "set xrange [%f:%f]\n", Tc_min, Tc_max);
    fprintf(gp, "set boxwidth %f\n", bin_width);
    fprintf(gp, "plot '-' with boxes title \"hist\"\n");
    for (int i = 0; i < bins; i++) {
      x = Tc_min + (i + 0.5) * bin_width;
      y = hist[i];
      fprintf(gp, "%.12f %d \n", x, y);
    }
    fprintf(gp, "e\n");
    pclose(gp);
  }

  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "plot '-' w l title \"Tc\" \n");
  for (int i = 0; i < Tc_s.size(); i += 10) {
    fprintf(gp, "%d %.12f \n", i + 1, Tc_s[i]);
  }
  fprintf(gp, "e\n");
  pclose(gp);


  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "plot '-' w l title \"k\" \n");
  for (int i = 0; i < k_s.size(); i += 10) {
    fprintf(gp, "%d %.12f \n", i + 1, k_s[i]);
  }
  fprintf(gp, "e\n");
  pclose(gp);


  END();
}