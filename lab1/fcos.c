#define CUBE(x) (x * x * x)

static const double PI = 3.141502653589793;
static const double PI_2 = PI/2;

// Aproximacao por serie de Taylor centrada em x0 = pi/2 da funcao cosseno
double fcos(double x) {
  // x0 = pi/2
  // cos(x0) = 0
  // sin(x0) = 1  
  double x_minus_x0 = x - PI_2;

  return (-1 * x_minus_x0) 
       + (1 * CUBE(x_minus_x0) / 6);
}