#include <stdio.h>
#include "4-pontofixo.h"

double myfunc(double x) {
  return -(x * x) + 1.8 * x + 2.5;
}

int main(void) {
  printf("%f\n", pontofixo(5, 0.0005, myfunc));
  return 0;
}