#include <vector>
#include "T.hpp"

int main() {
  START(1);

  int bins;
  double Tc_min, Tc_max;
  double EPS, RW_step;
  double J_init[2];

  // Read file
  int dummy;
  FILE* fp = fopen("T_entropy.dat", "r");
  dummy = fscanf(fp, "%d", &bins);
  dummy = fscanf(fp, "%lf %lf", &Tc_min, &Tc_max);
  dummy = fscanf(fp, "%lf", &EPS);
  dummy = fscanf(fp, "%lf", &RW_step);
  dummy = fscanf(fp, "%lf %lf", &J_init[0], &J_init[1]);

  double* S_s = alloc_dvector(bins);
  for (int i = 0; i < bins; i++)
    dummy = fscanf(fp, "%lf", &S_s[i]);
  fclose(fp);

  double bin_width = (Tc_max - Tc_min) / (double)bins;

  double x, y;
  FILE* gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set terminal pdfcairo color enhanced size 4in, 3in\n");
  fprintf(gp, "set output '../fig/T_entropy.pdf'\n");
  fprintf(gp, "set xrange [%f:%f]\n", Tc_min - 0.02, Tc_max + 0.02);
  fprintf(gp, "set xlabel \"T_c\"\n");
  fprintf(gp, "set ylabel \"Entropy\"\n");
  fprintf(gp, "plot '-' with linespoint ps 0.3 notitle\n");
  for (int i = 0; i < bins; i++) {
    x = Tc_min + (i + 0.5) * bin_width;
    y = S_s[i];
    fprintf(gp, "%.12f %.12f \n", x, y);
  }
  fprintf(gp, "e\n");
  pclose(gp);

  END();
}