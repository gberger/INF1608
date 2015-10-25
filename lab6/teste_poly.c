#include <stdio.h>
#include <math.h>
#include "poly.h"

int main(void) {
  Sample* cheb = Chebyshev(6, 0, HALF_PI, cos);
  double* poly_lagrange = LagrangeCompute(cheb);
  double* poly_newton = NewtonCompute(cheb);

}