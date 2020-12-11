#include <vector>
#include "T.hpp"

int main() {
  START(1);

  int samples = 10000;

  double S_Tc_max = 0.567296328553;
  double Tc_low_max = S_Tc_max + 0.003;
  double Tc_high_min = S_Tc_max + 0.002;
  double Tc_boundary = (Tc_low_max + Tc_high_min) / 2.0;

  double Tc_min = 0.1, Tc_max = 0.606826151086;
  int bins = 100;

  vector<double> Tc_s;
  vector<vector<double>> Js_s(0, vector<double>(3));
  FILE *fp_Tc, *fp_Js;
  int ret1, ret2;
  double Tc_tmp, Js_tmp[3];

  fp_Tc = fopen("../data/T_Tc_high.dat", "r");
  fp_Js = fopen("../data/T_Js_high.dat", "r");
  while (true) {
    ret1 = fscanf(fp_Tc, "%lf", &Tc_tmp);
    ret2 = fscanf(fp_Js, "%lf\t%lf\t%lf", &Js_tmp[0], &Js_tmp[1], &Js_tmp[2]);
    if (ret1 == -1 || ret2 == -1)
      break;
    if (Tc_tmp <= Tc_boundary)
      continue;
    Tc_s.push_back(Tc_tmp);
    Js_s.push_back({Js_tmp[0], Js_tmp[1], Js_tmp[2]});
  }

  fp_Tc = fopen("../data/T_Tc_low.dat", "r");
  fp_Js = fopen("../data/T_Js_low.dat", "r");
  while (true) {
    ret1 = fscanf(fp_Tc, "%lf", &Tc_tmp);
    ret2 = fscanf(fp_Js, "%lf\t%lf\t%lf", &Js_tmp[0], &Js_tmp[1], &Js_tmp[2]);
    if (ret1 == -1 || ret2 == -1)
      break;
    if (Tc_tmp >= Tc_boundary)
      continue;
    Tc_s.push_back(Tc_tmp);
    Js_s.push_back({Js_tmp[0], Js_tmp[1], Js_tmp[2]});
  }


  int exceed = Tc_s.size() - samples;
  int discard;
  for (int _ = 0; _ < exceed; _++) {
    discard = mt() % (Tc_s).size();
    Tc_s.erase(Tc_s.begin() + discard);
    Js_s.erase(Js_s.begin() + discard);
  }

  // Write to file
  FILE* fp = fopen("../data/T_Js.dat", "w");
  for (int i = 0; i < samples; i++) {
    fprintf(fp, "%.12f\t%.12f\t%.12f\n", Js_s[i][0], Js_s[i][1], Js_s[i][2]);
  }

  fp = fopen("../data/T_Tc.dat", "w");
  for (int i = 0; i < samples; i++) {
    fprintf(fp, "%.12f\n", Tc_s[i]);
  }

  // Visualize
  double bin_width = (Tc_max - Tc_min) / (double)bins;
  int* hist = alloc_ivector(bins);
  for (int i = 0; i < bins; i++)
    hist[i] = 0;
  for (int i = 0; i < samples; i++)
    hist[Tc_to_ind(Tc_s[i], Tc_min, bin_width, bins)]++;

  double x, y, ymax = 0.0;
  for (int i = 0; i < bins; i++)
    ymax = max(ymax, (double)hist[i]);
  ymax *= 1.1;

  FILE* gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set terminal pdfcairo color enhanced size 4in, 3in\n");
  fprintf(gp, "set output '../data/T_hist.pdf'\n");
  fprintf(gp, "set style fill solid border lc rgb \"black\"\n");
  fprintf(gp, "set xrange [%f:%f]\n", Tc_min - 0.02, Tc_max + 0.02);
  fprintf(gp, "set yrange [0:%f]\n", ymax);
  fprintf(gp, "set xlabel \"T_c\"\n");
  fprintf(gp, "set ylabel \"Frequency\"\n");
  fprintf(gp, "set boxwidth %f\n", bin_width);
  fprintf(gp, "plot '-' with boxes notitle \"hist\"\n");
  for (int i = 0; i < bins; i++) {
    x = Tc_min + (i + 0.5) * bin_width;
    y = hist[i];
    fprintf(gp, "%.12f %.1f \n", x, y);
  }
  fprintf(gp, "e\n");
  pclose(gp);

  END();
}