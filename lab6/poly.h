typedef struct sample Sample;

Sample* Chebyshev(int n, double a, double b, double (*f) (double x));

double* LagrangeCompute(Sample* s);

double LagrangeEval(Sample* s, double* den, double x);

double* NewtonCompute(Sample* s);

double NewtonEval(Sample* s, double* coef, double x);

