#include <math.h>

double bissecao (double a, double b, int p, double (*f) (double x)) {
  int i;
  int n = ceil( log( (b-a)/(0.5 * pow(10, -p)) ) / log(2) ) - 1;

  double fa = f(a);
  double fc, c;

  for(i = 0; i < n; i++) {
    c = (a+b)/2;
    fc = f(c);
    if (fc == 0) {
      return c;
    } else if (fa * fc > 0) {
      a = c;
      fa = fc;
    } else {
      b = c;
    }
  }

  return c;
}