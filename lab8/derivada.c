#include <math.h>
#include <stdio.h>
#include <limits.h>

#define EPSILON 0.0000000000000002220446

double f(double x) {
  return 1 / (x + 1);
}

// Obtida pelo WolframAlpha
double df(double x) {
  return -1 / (x*(x+2) + 1);
}

// Obtida pelo WolframAlpha
double d3f(double x) {
  return -6 / (x*(x*(x*(x+4)+6)+4)+1);
}


double derivada_numerica(double x, double h, double (*f) (double x)) {
  return (f(x+h) - f(x-h)) / (2*h);
}

double get_erro_teorico(double x, double h, double (*d3f) (double x)) {
  return fabs(  h * h * d3f(x) / 6   +   EPSILON / h  );
}

double get_erro(double y_numerico, double y_analitico) {
  return fabs(y_numerico - y_analitico);
}


////////////////////////////////////////////////
int main(void) {
  int i;
  double h = 1;
  double x = 1;

  double y_numerico, y_analitico, erro, erro_teorico;
  int i_min_erro, i_min_erro_teorico;
  double min_erro = INT_MAX, min_erro_teorico = INT_MAX;


  y_analitico = df(x);

  printf("h\ty_numerico \t\ty_analitico \t\terro \t\t\terro_teorico\n");
  for (i = 1; i <= 12; i++) {
    h = h / 10;

    y_numerico = derivada_numerica(x, h, f);
    erro = get_erro(y_numerico, y_analitico);
    erro_teorico = get_erro_teorico(x, h, d3f);

    printf("10^-%d\t%.15f\t%.15f\t%.15f\t%.15f\n", i, y_numerico, y_analitico, erro, erro_teorico);

    if (erro < min_erro) {
      min_erro = erro;
      i_min_erro = i;
    }

    if (erro_teorico < min_erro_teorico) {
      min_erro_teorico = erro_teorico;
      i_min_erro_teorico = i;
    }
  }

  printf("\n");
  printf("h que minimiza o erro numérico: 10^-%d\n", i_min_erro);
  printf("h que minimiza o erro teórico:  10^-%d\n", i_min_erro_teorico);

  return 0;
}