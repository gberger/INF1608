#include "sparse.h"

// Metodo do Gradiente Conjugado
// Retorna: numero de iteracoes
// Modifica: vetor x
int ConjugateGradient (int n, Sparse** A, double* b, double* x);


// Metodo do Gradiente Conjugado com uso de Pre-Condicionador
// Retorna: numero de iteracoes
// Modifica: vetor x
int ConjugateGradientPC (int n, Sparse** A, double* b, double* x, Sparse** M);
