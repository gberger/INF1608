#include <stdio.h>
#include <math.h>
#include "1-bissecao.h"

#define G 9.8
#define C 15.0
#define V 35.0
#define T 9.0


double paraquedista(double m) {
  return G * m * (1 - exp(- C * T / m)) / C - V;
}

int main(void) {
  printf("%f\n", bissecao(50, 200, 6, paraquedista));
  return 0;
}