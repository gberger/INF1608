#include <stdlib.h>

double** mat_cria(int m, int n){
  int i;
  double** A = malloc(sizeof(double*) * m);

  for (i = 0; i < m; i++) {
    A[i] = calloc(n, sizeof(double));
  }

  return A;
}

void mat_libera(int m, double** A) {
  int i;
  for (i = 0; i < m; i++) {
    free(A[i]);
    A[i] = NULL;
  }
  free(A);
}

void mat_transposta(int m, int n, double** A, double** T){
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      T[j][i] = A[i][j];
    }
  }
}

void mat_multv(int m, int n, double** A, double* v, double* w){
  int i, j;

  for (i = 0; i < m; i++) {
    w[i] = 0.0;
    for (j = 0; j < n; j++) {
      w[i] += A[i][j] * v[j];
    }
  }
}

void mat_multm(int m, int n, int q, double** A, double** B, double** C){
  int i, j, k;

  for (i = 0; i < m; i++) {
    for (k = 0; k < q; k++) {
      C[i][k] = 0.0;
      for (j = 0; j < n; j++) {
        C[i][k] += A[i][j] * B[j][k];
      }
    }
  }
}
