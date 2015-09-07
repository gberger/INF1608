#include <math.h>

int falsaposicao (double a, double b, int p, double (*f) (double x), double* r) {
  int i = 0;
  double c, fc;
  double maxerr = 0.5 * pow(10.0, -p);
  double fb = f(b);

  do {
    c = a - ( (f(a)*(b-a)) / (f(b) - f(a)) );
    fc = f(c);
    if (f(c) * fb < 0) {
      a = c;
    } else {
      fb = fc;
      b = c;
    }
    i++;
  } while (fc >= maxerr);

  *r = c;
  return i;
}


int newtonraphson (double x0, int p, double (*f) (double x), double (*fl) (double x), double* r) {
  int i= 0;
  double xi;
  double maxerr = 0.5 * pow(10.0, -p);

  do {
    xi = x0 - f(x0) / fl(x0);
    x0 = xi;
    i++;
  } while (fabs(f(xi)) >= maxerr);

  *r = xi;
  return i;
}