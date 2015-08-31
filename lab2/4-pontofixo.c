#include <math.h>
#include <stdio.h>

double pontofixo(double x0, double epsilon, double (*g) (double x)) {
  double xn, xn1 = x0;

  do {
    xn = xn1;
    xn1 = sqrt(1.8 * xn + 2.5);;
  } while (fabs(xn - xn1) >= epsilon);

  return xn;
}
