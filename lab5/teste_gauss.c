#include <stdio.h>
#include "gauss.h"
#include "matriz.h"

double** mat_from_vector(double *v, int m, int n) {
  int i, j;
  double** A = mat_cria(m, n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A[i][j] = v[i*n + j];
    }
  }

  return A;
}

void mat_exibe(int n, double **M) {
  int i, j;
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      printf("%f ", M[i][j]);
    }
    printf("\n");
  }
  printf("\n"); 
}

void testa(int n, double* a, double* b) {
  int i;
  double** A = mat_from_vector(a, n, n);

  printf("A = \n");
  mat_exibe(n, A);

  double** P = fatoracao(n, A);
  double** LU = A;

  printf("P = \n");
  mat_exibe(n, P);

  printf("L | U = \n");
  mat_exibe(n, LU);

  printf("x = \n");
  double *x = substituicao(n, A, P, b);
  for (i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }

  printf("\n\n");
}

int main(void) {
  int n = 3;
  double a[9] = {1, 2, -1, 
                 2, 1, -2,
                 -3, 1, 1};
  double b[3] = {3, 3, -6};

  testa(n, a, b);


  int nn = 6;
  double aa[36] = {3, -1, 0, 0, 0, 0.5,
                   -1, 3, -1, 0, 0.5, 0,
                   0, -1, 3, -1, 0, 0,
                   0, 0, -1, 3, -1, 0,
                   0, 0.5, 0, -1, 3, -1,
                   0.5, 0, 0, 0, -1, 3};
  double bb[6] = {2.5, 1.5, 1, 1, 1.5, 2.5};

  testa(nn, aa, bb);

  return 0;
}