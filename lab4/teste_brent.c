#include <stdio.h>
#include <math.h>
#include "brent.h"

double f1(double x) {
  return 4 * cos(x) - exp(2 * x);
}
double f2(double x) {
  return x/2.0 - tan(x);
}
double f3(double x) {
  return 1 - x * log(x + 20);
}
double f4(double x) {
  return pow(2, x) - 3 * x;
}
double f5(double x) {
  return pow(x, 3) + x - 1;
}

int main(void) {
  printf("f1 = 0 em x = %f\n\n", brent(-1, 1, 6, f1));
  printf("f2 = 0 em x = %f\n\n", brent(-1, 1, 6, f2));
  printf("f3 = 0 em x = %f\n\n", brent(-1, 1, 6, f3));
  printf("f4 = 0 em x = %f\n\n", brent(-1, 1, 6, f4));
  printf("f5 = 0 em x = %f\n\n", brent(-1, 1, 6, f5));

  return 0;
}