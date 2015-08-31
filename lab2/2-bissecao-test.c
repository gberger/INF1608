#include <stdio.h>
#include <math.h>
#include "1-bissecao.h"

double myfunc(double x) {
  return sin(x) - pow(x, 3);
}

int main(void) {
  printf("%f\n", bissecao(0.5, 1.5, 6, myfunc));
  return 0;
}