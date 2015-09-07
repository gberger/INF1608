#include <stdio.h>
#include <math.h>
#include "raiz2.h"

double f (double x) {
  return sin(x) - pow(x, 3);
}

double fl (double x) {
  return cos(x) - 3 * pow(x, 2);
}

int main (void) {
  double a = 0.5;
  double b = 1.5;
  double p = 6;

  double r;
  int i;

  i = (int) ceil( log( (b-a)/(0.5 * pow(10, -p)) ) / log(2) ) - 1;
  printf("BISSCAO\n");
  printf("iteracoes: %d\n\n", i);

  i = falsaposicao(a, b, p, f, &r);
  printf("FALSA POSICAO\n");
  printf("iteracoes: %d\nraiz: %f\n\n", i,  r);

  i = newtonraphson(b, p, f, fl, &r);
  printf("NEWTON RAPHSON\n");
  printf("iteracoes: %d\nraiz: %f\n", i,  r);

  return 0;
}