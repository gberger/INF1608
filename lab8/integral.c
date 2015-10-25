#include <math.h>
#include <stdio.h>

double Simpson(double a, double b, double m, double (*f) (double x)) {
  int i;
  double h = (b - a)/(2 * m);
  double xi = a;

  double soma_pares = 0;
  double soma_impares = 0;
  for (i = 1; i < 2*m; i++) {
    xi += h;

    if(i % 2 == 0){
      soma_pares += f(xi);
    } else {
      soma_impares += f(xi);
    }
  }

  return h/3 * (f(a) + f(b) + 4 * soma_impares + 2 * soma_pares);
}

double PontoMedio(double a, double b, double m, double (*f) (double x)) {
  int i;
  double h = (b - a)/m;

  double w;
  double soma_fw = 0;
  double xnext;
  double xi = a;

  for (i = 1; i <= m; i++) {
    xnext = xi + h;
    w = (xnext + xi)/2;
    soma_fw += f(w);

    xi = xnext;
  }

  return h * soma_fw;
}

//////////////////////////////////
// 0, 4
double f1(double x) {
  return x / sqrt(x * x + 9);
}

// 1, 3
double f2(double x) {
  return x * x * log(x);
}

// 0, PI
double f3(double x) {
  return x * x * sin(x);
}

/////////////////////////////////
int main(void) {
  printf("f1\n");
  printf("Valor real: 2.0000000000000\n");
  printf("S[m=16]:    %.13f\n", Simpson(0, 4, 16, f1));
  printf("S[m=32]:    %.13f\n", Simpson(0, 4, 32, f1));
  printf("M[m=16]:    %.13f\n", PontoMedio(0, 4, 16, f1));
  printf("M[m=32]:    %.13f\n", PontoMedio(0, 4, 32, f1));

  printf("\nf2\n");
  printf("Valor real: 6.9986217091241\n");
  printf("S[m=16]:    %.13f\n", Simpson(1, 3, 16, f2));
  printf("S[m=32]:    %.13f\n", Simpson(1, 3, 32, f2));
  printf("M[m=16]:    %.13f\n", PontoMedio(1, 3, 16, f2));
  printf("M[m=32]:    %.13f\n", PontoMedio(1, 3, 32, f2));

  printf("\nf3\n");
  printf("Valor real: 5.8696044010894\n");
  printf("S[m=16]:    %.13f\n", Simpson(0, M_PI, 16, f3));
  printf("S[m=32]:    %.13f\n", Simpson(0, M_PI, 32, f3));
  printf("M[m=16]:    %.13f\n", PontoMedio(0, M_PI, 16, f3));
  printf("M[m=32]:    %.13f\n", PontoMedio(0, M_PI, 32, f3));


  return 0;
}