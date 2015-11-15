#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 9.81
#define D 0.5
#define E 2.0

double LinearInterpolation(double y0, double y1, double x0, double x1, double x) {
  return y0 + (y1-y0) * (x-x0) / (x1-x0);
}

double RungeKutta(double t0, double t1, double h, double y0, double (*f) (double t, double y)) {
  double s1, s2, s3, s4;
  double t, y, yn, tn;

  tn = t0;
  yn = y0;

  while (tn < t1) {
    t = tn;
    y = yn;

    s1 = h * f(t,       y);
    s2 = h * f(t + h/2, y + s1/2);
    s3 = h * f(t + h/2, y + s2/2);
    s4 = h * f(t + h,   y + s3);

    yn = y + 1.0/6 * (s1 + 2 * s2 + 2 * s3 + s4);
    tn = t + h;
  }

  return yn;
}


double TimeToY1(double t0, double y0, double y1, double h, double (*f) (double t, double y)) {
  double s1, s2, s3, s4;
  double t, y, yn, tn;

  tn = t0;
  yn = y0;

  int crescente = yn < y1;

  while ((crescente && yn < y1) || (!crescente && yn > y1)) {
    t = tn;
    y = yn;

    s1 = h * f(t,       y);
    s2 = h * f(t + h/2, y + s1/2);
    s3 = h * f(t + h/2, y + s2/2);
    s4 = h * f(t + h,   y + s3);

    yn = y + 1.0/6 * (s1 + 2 * s2 + 2 * s3 + s4);
    tn = t + h;
    // printf("%f, %f\n", tn, yn);
  }

  if (yn == y1) {
    return tn;
  } else {
    return LinearInterpolation(t, tn, y, yn, y1);
  }
}


/*************
 *** TESTE ***
 *************/

double yprime(double t, double y) {
  return t * y + pow(t, 3);
}

double yt(double t) {
  double t2 = pow(t, 2);
  return exp(t2/2) - t2 - 2;
}

double A(double y) {
  int i;
  double a[8] = {0, 18, 32, 45, 67, 97, 117, 137};

  for(i = 0; i <= 7; i++) {
    if(y == i) {
      return a[i];
    } else if(y < i) {
      return LinearInterpolation(a[i-1], a[i], i-1, i, y);
    }
  }

  return -1;
}

double reservatorio_prime(double t, double y) {
  return - (M_PI * pow(D, 2)) / (4 * A(y)) * sqrt(2 * G * (y + E));
}

int main(void) {
  double expected = yt(2.4);
  double actual, err;

  printf("RungeKutta\n");
  printf("esperado:  %f\n", expected);

  actual = RungeKutta(0, 2.4, 0.1, -1, yprime);
  err = fabs((actual-expected)/expected)*100;
  printf("h = 0.1:   %f, erro: %.3f%%\n", actual, err);

  actual = RungeKutta(0, 2.4, 0.01, -1, yprime);
  err = fabs((actual-expected)/expected)*100;
  printf("h = 0.01:  %f, erro: %.3f%%\n", actual, err);

  actual = RungeKutta(0, 2.4, 0.001, -1, yprime);
  err = fabs((actual-expected)/expected)*100;
  printf("h = 0.001: %f, erro: %.3f%%\n", actual, err);



  printf("\n\nTimeToY1 (a)\n");
  printf("esperado:  2.4\n");
  printf("h = 0.1:   %f\n", TimeToY1(0, -1, expected, 0.1, yprime));
  printf("h = 0.01:  %f\n", TimeToY1(0, -1, expected, 0.01, yprime));
  printf("h = 0.001: %f\n", TimeToY1(0, -1, expected, 0.001, yprime));


  printf("\n\nTimeToY1 (b) reservatorio\n");
  printf("%f horas\n", TimeToY1(0, 6, 0, 0.01, reservatorio_prime));

  return 0;
}
