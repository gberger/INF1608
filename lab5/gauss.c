#include <math.h>
#include <stdlib.h>
#include "matriz.h"

static double **pivot(int n, double** A) {
  int i, j, k, jmax, tmp;
  double** P = mat_cria(n, n);

  // Matriz identidade
  for (i = 0; i < n; i++) {
    P[i][i] = 1; 
  }

  for (i = 0; i < n; i++)  {
    jmax = i;
    for (j =  i; j <  n; j++) {
      if (fabs(A[j][i]) > fabs(A[jmax][i])) {
        jmax = j;
      }
    }

    if (jmax != i){
      for (k = 0; k < n; k++) { 
        tmp = P[i][k]; 
        P[i][k] = P[jmax][k]; 
        P[jmax][k] = tmp;
      }
    }
  }

  return P;
}

double **fatoracao (int n, double** A) {
  int i, j, k;
  double s;

  double **P = pivot(n, A);
  double **PA = mat_cria(n, n);
  mat_multm(n, n, n, P, A, PA);

  double **L = mat_cria(n, n);
  double **U = mat_cria(n, n);

  // A diagonal de L é 1
  for (i = 0; i < n; i++)  { 
    L[i][i] = 1; 
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j <= i) {
        s = 0; 
        for (k = 0; k <  j; k++) {
          s+=  L[j][k] * U[k][i];
        }
        U[j][i] = PA[j][i] - s;
      }
      if (j >= i) {
        s = 0;
        for (k = 0; k <  i; k++) {
          s+=  L[j][k] * U[k][i];
        }
        L[j][i] = (PA[j][i] - s) / U[i][i];
      }
    }
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i <= j) {
        A[i][j] = U[i][j];
      } else {
        A[i][j] = L[i][j];
      }
    }
  }

  // Liberar
  mat_libera(n, PA);
  mat_libera(n, L);
  mat_libera(n, U);

  return P;
}

double *substituicao (int n, double** A, double** P, double* b) {
  int i, j;
 
  // Reconstruir L e U
  double **L = mat_cria(n, n);
  double **U = mat_cria(n, n);

  for (i = 0; i < n; i++)  { 
    L[i][i] = 1; 
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i <= j) {
        U[i][j] = A[i][j];
      } else {
        L[i][j] = A[i][j];
      }
    }
  }

  // Achar y, top-down
  // Ly = Pb
  double *Pb = calloc(n, sizeof(double));
  double *x = calloc(n, sizeof(double));
  double *y = calloc(n, sizeof(double));
  mat_multv(n, n, P, b, Pb);

  for (i = 0; i < n; i++) {
    double soma = 0;
    for (j = 0; j < i; j++) {
      soma += L[i][j] * y[j];
    }
    y[i] = (Pb[i] - soma) / L[i][i];
  }
  
  // Achar x, bottom-up
  // Ux = y;
  for (i = n-1; i >= 0; i--) {
    double soma = 0;
    for (j = i+1; j < n; j++) {
      soma += U[i][j] * x[j];
    }
    x[i] = (y[i] - soma) / U[i][i];
  }

  // Liberar
  mat_libera(n, L);
  mat_libera(n, U);
  free(y);
  free(Pb);


  return x;
}
