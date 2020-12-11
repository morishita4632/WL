#include <vector>
#include "S.hpp"

int main() {
  START(5);

  int bins = 100;
  double Tc_min = 0.1, Tc_max = 0.567296328554;
  double EPS = 1e-12;
  double RW_step = 4.0;
  double k_init = 1.0;

  double flat_coeff = 0.9;
  double f_min = 1e-9;

  int* hist = alloc_ivector(bins);
  double* S_s = alloc_dvector(bins);
  for (int i = 0; i < bins; i++) {
    hist[i] = 0, S_s[i] = 0.0;
  }

  vector<double> Tc_s, k_s;
  vector<vector<double>> Js_s(0, vector<double>(2));

  double bin_width = (Tc_max - Tc_min) / (double)bins;
  double ck = k_init, dk, nk;
  double cJs[2], nJs[2];
  k_to_Js(ck, cJs);
  double cTc = Js_to_Tc(cJs, EPS), nTc;
  int cind = Tc_to_ind(cTc, Tc_min, bin_width, bins), nind;
  double f = 1.0;

  int CNT = 0;
  int MAX_CNT = 100000000;
  while (f > f_min) {
    do {
      CNT++;
      dk = (rand01() * 2 - 1) * RW_step;
      nk = ck + dk;
      k_to_Js(nk, nJs);
      nTc = Js_to_Tc(nJs, EPS);
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
  cout << CNT << endl;

  // Write to file.
  FILE* fp = fopen("S_entropy.dat", "w");
  fprintf(fp, "%d\n", bins);
  fprintf(fp, "%.12f %.12f\n", Tc_min, Tc_max);
  fprintf(fp, "%.12f\n", EPS);
  fprintf(fp, "%.12f\n", RW_step);
  fprintf(fp, "%.12f\n", k_init);
  for (int i = 0; i < bins; i++)
    fprintf(fp, "%.12f\n", S_s[i]);
  fclose(fp);

  // Visualize
  double x, y, ymax = 0.0;
  for (int i = 0; i < bins; i++)
    ymax = max(ymax, (double)hist[i]);
  ymax *= 1.1;

  FILE* gp;
  if (f > f_min) {
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set style fill solid border lc rgb \"black\"\n");
    fprintf(gp, "set xrange [%f:%f]\n", Tc_min, Tc_max);
    fprintf(gp, "set yrange [0:%f]\n", ymax);
    fprintf(gp, "set boxwidth %f\n", bin_width);
    fprintf(gp, "plot '-' with boxes title \"hist\"\n");
    for (int i = 0; i < bins; i++) {
      x = Tc_min + (i + 0.5) * bin_width;
      y = hist[i];
      fprintf(gp, "%.12f %.1f \n", x, y);
    }
    fprintf(gp, "e\n");
    pclose(gp);
  }

  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "plot '-' w l title \"Tc\" \n");
  for (int i = 0; i < Tc_s.size(); i += 1000) {
    fprintf(gp, "%d %.12f \n", i + 1, Tc_s[i]);
  }
  fprintf(gp, "e\n");
  pclose(gp);


  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "plot '-' w l title \"k\" \n");
  for (int i = 0; i < k_s.size(); i += 1000) {
    fprintf(gp, "%d %.12f \n", i + 1, k_s[i]);
  }
  fprintf(gp, "e\n");
  pclose(gp);


  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "plot '-' w l title \"S\" \n");
  for (int i = 0; i < bins; i++) {
    x = Tc_min + (i + 0.5) * bin_width;
    y = S_s[i];
    fprintf(gp, "%.12f %.12f \n", x, y);
  }
  fprintf(gp, "e\n");
  pclose(gp);

  END();
}