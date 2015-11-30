#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"


/******************
 *** AUXILIARES ***
 ******************/

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

void mat_print (int ii, int jj, double** A) {
  int i, j;
  for (i = 0; i < ii; i++) {
    printf("[ ");
    for (j = 0; j < jj; j++) {
      printf("%f ", A[i][j]);  
    }
    printf("]\n");
  }

  printf("\n");
}

void vet_print (int n, double* v) {
  int i;
  printf("[ ");
  for (i = 0; i < n; i++) {
    printf("%f ", v[i]);
  }
  printf("]\n\n");
}

double* vet_cria(int n) {
  return (double*) calloc(n, sizeof(double));
}

// M = v vT
void vet_vvT (int n, double* v, double** M) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      M[i][j] = v[i] * v[j];
    }
  }
}

// a + b = c
void vet_soma(int n, double* a, double* b, double* c) {
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] + b[i];
  }
}

// a - b = c
void vet_subtrai(int n, double* a, double* b, double* c) {
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] - b[i];
  }
}

// dst <- src
void vet_copia(int n, double* src, double* dst) {
  int i;
  for (i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}

// Norma 2
double norma2 (int n, double* v) {
  int i;
  double soma = 0;
  for (i = 0; i < n; i++) {
    soma += v[i] * v[i];
  }
  return sqrt(soma);
}

int vet_iszero(int n, double* v) {
  return norma2(n, v) == 0;
}

// w = sv
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



/***************
 *** METODOS ***
 ***************/

void Cholesky (int n, double** A) {
  int k, i, j, m;
  double** uuT;

  for (k = 0; k < n; k++) {
    A[k][k] = sqrt(A[k][k]);

    for (j = 0; j < k; j++) {
      A[k][j] = 0;
    }

    for (j = k+1; j < n; j++) {
      A[k][j] = A[k][j]/A[k][k];
    }

    m = n - k - 1;
    if(m > 0){
      uuT = mat_cria(m, m);
      vet_vvT(m, A[k] + 1 + k, uuT);

      for (i = k+1; i < n; i++) {
        for (j = k+1; j < n; j++) {
          A[i][j] = A[i][j] - uuT[i-k-1][j-k-1];
        }
      }
      mat_libera(m, uuT);
    }
  }
}


void ConjugateGradient (int n, double** A, double* b, double* x) {
  int k;
  double* ax = vet_cria(n);
  double* r = vet_cria(n);
  double* d = vet_cria(n);

  double alpha;
  double* Ad = vet_cria(n);

  double* xn = vet_cria(n);
  double* ad = vet_cria(n);
  
  double* rn = vet_cria(n);
  double* aAd = vet_cria(n);

  double beta;

  double* dn = vet_cria(n);
  double* bd = vet_cria(n);


  mat_multv(n, n, A, x, ax);
  vet_subtrai(n, b, ax, r);
  vet_copia(n, r, d);
  

  for (k = 0; k < n; k++) {

    if (vet_iszero(n, r)) {
      break;
    }
    
    // Alpha_k = rTk . rK / dTk . A dk

    mat_multv(n, n, A, d, Ad);
    alpha = vet_dot(n, r, r) / vet_dot(n, d, Ad);

    // x_{k+1} = xk + Alpha_k dk;
    vet_mults(n, d, alpha, ad);
    vet_soma(n, x, ad, xn);

    // r_{k+1} = rk - Alpha_k A dk;
    vet_mults(n, Ad, alpha, aAd);
    vet_subtrai(n, r, aAd, rn);

    // Beta_k = rTn rn / 
    beta = vet_dot(n, rn, rn) / vet_dot(n, r, r);

    // d_{k+1} = r_{k+1} + Beta_k dk;
    vet_mults(n, d, beta, bd);
    vet_soma(n, rn, bd, dn);


    vet_copia(n, xn, x);
    vet_copia(n, rn, r);
    vet_copia(n, dn, d);

  }
}


/**************
 *** TESTES ***
 **************/

int main (void) {
  int n = 3;
  double a1[9] = {1.0, -1.0, 0.0,
                 -1.0,  2.0, 1.0,
                  0.0,  1.0, 2.0};
  double **A1 = mat_from_vector(a1, n, n);
  double b1[3] = {0.0, 2.0, 3.0};
  double *x1 = vet_cria(n);

  double a2[9] = {1.0, -1.0, 0.0,
                 -1.0,  2.0, 1.0,
                  0.0,  1.0, 5.0};
  double **A2 = mat_from_vector(a2, n, n);
  double b2[3] = {3.0, -3.0, 4.0};
  double *x2 = vet_cria(n);


  Cholesky(n, A1);
  mat_print(n, n, A1); 

  A1 = mat_from_vector(a1, n, n);
  ConjugateGradient(n, A1, b1, x1);
  vet_print(n, x1);


  Cholesky(n, A2);
  mat_print(n, n, A2); 

  A2 = mat_from_vector(a2, n, n);
  ConjugateGradient(n, A2, b2, x2);
  vet_print(n, x2);

  return 0;
}
