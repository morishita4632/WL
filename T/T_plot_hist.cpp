#include "T.hpp"
#include <vector>

int main() {
  START(1);

  int samples = 10000;

  int bins;
  double Tc_min, Tc_max;

  // Read file
  int dummy;
  FILE* fp = fopen("T_entropy.dat", "r");
  dummy = fscanf(fp, "%d", &bins);
  dummy = fscanf(fp, "%lf %lf", &Tc_min, &Tc_max);
  fclose(fp);

  double bin_width = (Tc_max - Tc_min) / (double)bins;

  fp = fopen("../data/T_Tc.dat", "r");
  int* hist = alloc_ivector(bins);
  double Tc;
  for (int i = 0; i < samples; i++) {
    dummy = fscanf(fp, "%lf", &Tc);
    hist[Tc_to_ind(Tc, Tc_min, bin_width, bins)]++;
  }
  fclose(fp);

  // Visualize
  double x, y, ymax = 0.0;
  for (int i = 0; i < bins; i++)
    ymax = max(ymax, (double)hist[i]);
  ymax *= 1.1;

  FILE* gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set terminal pdfcairo color enhanced size 4in, 3in\n");
  fprintf(gp, "set output '../fig/T_hist.pdf'\n");
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