#include <math.h>

struct sample {
  int n;      // numero de amostras
  double* x;  // valores x das amostras
  double* y;  // valores y das amostras
};


Sample* Chebyshev(int n, double a, double b, double (*f) (double x)){
  int beta;
  for (beta = 1; beta < 2*n; beta += 2) {
    
  }
}

double* LagrangeCompute(Sample* s){
  
}

double LagrangeEval(Sample* s, double* den, double x){
  
}

double* NewtonCompute(Sample* s){
  
}

double NewtonEval(Sample* s, double* coef, double x){
  
}

