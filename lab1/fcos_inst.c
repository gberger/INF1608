#include <stdio.h>
#include <math.h>

double fcos(double x);

double fcos_error(x) {
  return fabs(fabs(fcos(x)) - fabs(cos(x)));
}

double max_fcos_error() {
  double x;
  double max_error_x = 0;
  double max_error = fcos_error(0);
  double error;

  for (x = 0; x <= M_PI; x += M_PI/100000) {
    error = fcos_error(x);
    if (error > max_error) {
      max_error_x = x;
      max_error = error;
    }
  }

  printf("Max error is %f.\nfcos(%f) = %f\n cos(%f) = %f\n",
    max_error, max_error_x, fcos(max_error_x), max_error_x, cos(max_error_x));

  return max_error;
}

int main(void) {
  max_fcos_error();

  return 0;
}