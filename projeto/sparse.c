#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparse.h"

Sparse** sparse_cria(int n) {
  int i;
  Sparse** A = malloc(sizeof(Sparse*) * n);

  for (i = 0; i < n; i++) {
    A[i] = malloc(sizeof(Sparse));
    A[i][0].col = -1;
  }

  return A;
}

void sparse_libera(int n, Sparse** M) {
  int i;
  for (i = 0; i < n; i++) {
    
  }
  free(M);
}

void sparse_multv(int n, Sparse** M, double* v, double* w) {
  int i, j, col;
  for (i = 0; i < n; i++) {
    w[i] = 0.0;

    for (j = 0; j < n; j++) {
      col = M[i][j].col;
      if (col == -1) {
        break;
      }

      w[i] += M[i][j].val * v[col];
    }
  }
}

double sparse_get(int i, int j, Sparse** A) {
  int k = 0;

  while(A[i][k].col != -1) {
    if(A[i][k].col == j) {
      return A[i][k].val;
    }
    k++;
  }
  return 0;
}

Sparse** sparse_multm(int n, Sparse** A, Sparse** B) {
  int i, j, k, m;
  double num;
  Sparse** C = sparse_cria(n);

  for (i = 0; i < n; i++) {
    m = 0;
    C[i] = malloc(sizeof(Sparse) * (n+1));
    for (k = 0; k < n; k++) {
      num = 0.0;

      for (j = 0; j < n; j++) {
        if(A[i][j].col == -1) break;
        num += A[i][j].val * sparse_get(A[i][j].col, k, B);
      }

      if(num != 0) {
        C[i][m].col = k;
        C[i][m].val = num;
        m++;
      }
    }
    C[i][m].col = -1;
  }

  return C;
}


double* vet_cria(int n) {
  return (double*) calloc(n, sizeof(double));
}

void vet_soma(int n, double* a, double* b, double* c) {
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] + b[i];
  }
}

void vet_subtrai(int n, double* a, double* b, double* c) {
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] - b[i];
  }
}

void vet_copia(int n, double* src, double* dst) {
  int i;
  for (i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}

double vet_norma2 (int n, double* v) {
  int i;
  double soma = 0;
  for (i = 0; i < n; i++) {
    soma += v[i] * v[i];
  }
  return sqrt(soma);
}


void vet_mults(int n, double* v, double s, double* w) {
  int i;
  for (i = 0; i < n; i++) {
    w[i] = s * v[i];
  }
}

double vet_dot(int n, double* v, double* w) {
  int i;
  double s = 0;
  for (i = 0; i < n; i++) {
    s += v[i] * w[i];
  }
  return s;
}
