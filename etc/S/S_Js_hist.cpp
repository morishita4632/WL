#include <vector>
#include "S.hpp"

int main() {
  START(1);

  int samples = 10000;

  int bins = 100;
  double J_max = 0.5, J_min = 0.0;

  // Initialize
  int* hist = alloc_ivector(bins);
  for (int i = 0; i < bins; i++)
    hist[i] = 0;
  double** Js_s = alloc_dmatrix(samples, 2);

  double bin_width = (J_max - J_min) / (double)bins;

  // Read file
  FILE* fp = fopen("../data/S_Js.dat", "r");
  int ind;
  for (int i = 0; i < samples; i++) {
    fscanf(fp, "%.12lf\t%.12lf\n", Js_s[i][0], Js_s[i][1]);
    ind = Js_s[i][0] / bin_width;
    hist[ind]++;
  }

  // Visualize
  double x, y, ymax = 0.0;
  for (int i = 0; i < bins; i++)
    ymax = max(ymax, (double)hist[i]);
  ymax *= 1.1;

  FILE* gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set terminal pdfcairo color enhanced size 4in, 3in\n");
  fprintf(gp, "set output '../fig/S_Js.pdf'\n");
  fprintf(gp, "set style fill solid border lc rgb \"black\"\n");
  fprintf(gp, "set xrange [%f:%f]\n", J_min - 0.02, J_max + 0.02);
  fprintf(gp, "set yrange [0:%f]\n", ymax);
  fprintf(gp, "set xlabel \"T_c\"\n");
  fprintf(gp, "set ylabel \"Frequency\"\n");
  fprintf(gp, "set boxwidth %f\n", bin_width);
  fprintf(gp, "plot '-' with boxes notitle \"hist\"\n");
  for (int i = 0; i < bins; i++) {
    x = J_min + (i + 0.5) * bin_width;
    y = hist[i];
    fprintf(gp, "%.12f %.1f \n", x, y);
  }
  fprintf(gp, "e\n");
  pclose(gp);

  END();
}