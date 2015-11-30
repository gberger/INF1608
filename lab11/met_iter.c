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

double* vet_cria(int n) {
  return (double*) calloc(n, sizeof(double));
}

// a - b = c
void vet_subtrai(int n, double* a, double* b, double* c) {
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] - b[i];
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

// dst <- src
void vet_copia(int n, double* src, double* dst) {
  int i;
  for (i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}


/***************
 *** METODOS ***
 ***************/

int Jacobi (int n, double** A, double* b, double* x, double tol) {
  int iter = 0;
  int i, j;

  double** LU = mat_cria(n, n);
  double** Dinv = mat_cria(n, n);
  double* lux = vet_cria(n);
  double* blux = vet_cria(n);
  double* xn = vet_cria(n); 
  double* xdif = vet_cria(n);
  double nor;

  // obter inverso de D, e L+U
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) {
        Dinv[i][j] = 1.0/A[i][j];
      } else {
        LU[i][j] = A[i][j];
      }
    }
  }

  do {
    iter++;

    mat_multv(n, n, LU, x, lux);
    vet_subtrai(n, b, lux, blux);
    mat_multv(n, n, Dinv, blux, xn);

    vet_subtrai(n, x, xn, xdif);
    nor = norma2(n, xdif);

    vet_copia(n, xn, x);
  } while (nor > tol);

  mat_libera(n, LU);
  mat_libera(n, Dinv);
  free(lux);
  free(blux);
  free(xn);
  free(xdif);

  return iter;
}


int GaussSeidel (int n, double** A, double* b, double* x, double tol) {
  // Como esse metodo eh equivalente a SOR com w=1, essa funcao poderia
  // ser escrita como:
  // return SOR(n, A, b, x, tol, 1);

  int iter = 0;
  int i, j;

  double *xn = vet_cria(n);
  double *xdif = vet_cria(n);
  double sumu, suml, nor;

  do {
    iter++;

    // sumu eh U * x[k]
    // suml eh L * x[k+1]
    for (i = 0; i < n; i++) {
      sumu = 0;
      suml = 0;
      for (j = 0; j < n; j++) {
        if (j < i) {
          suml += A[i][j] * xn[j];
        } else if (j > i) {
          sumu += A[i][j] * x[j];
        }
      }

      // 1.0/A[i][i] é o inverso da diagonal
      xn[i] = 1.0/A[i][i] * (b[i] - sumu - suml);
    }

    vet_subtrai(n, x, xn, xdif);
    nor = norma2(n, xdif);

    vet_copia(n, xn, x);
  } while (nor > tol);

  free(xn);
  free(xdif);

  return iter;
}


int SOR (int n, double** A, double* b, double* x, double tol, double w) {
  int iter = 0;
  int i, j;

  double *xn = vet_cria(n);
  double *xdif = vet_cria(n);
  double sumu, suml, nor;

  do {
    iter++;

    // sumu eh U * x[k]
    // suml eh L * x[k+1]
    for (i = 0; i < n; i++) {
      sumu = 0;
      suml = 0;
      for (j = 0; j < n; j++) {
        if (j < i) {
          suml += A[i][j] * xn[j];
        } else if (j > i) {
          sumu += A[i][j] * x[j];
        }
      }

      xn[i] = (1.0 - w) * x[i] + w/A[i][i] * (b[i] - sumu - suml);
    }

    vet_subtrai(n, x, xn, xdif);
    nor = norma2(n, xdif);

    vet_copia(n, xn, x);
  } while (nor > tol);

  free(xn);
  free(xdif);

  return iter;
}


/**************
 *** TESTES ***
 **************/

int main (void) {
  double tol = pow(10, -7);
  double w = 1.07;
  int i;

  int n1 = 2;
  double a1[4] = {3.0, 1.0,
                 1.0, 2.0};
  double** A1 = mat_from_vector(a1, n1, n1);
  double b1[2] = {5.0, 5.0};
  double *x1;
  double r1[2] = {1.0, 2.0};


  int n2 = 6;
  double a2[36] = {3.0, -1.0, 0.0, 0.0, 0.0, 0.5,
                   -1.0, 3.0, -1.0, 0.0, 0.5, 0.0,
                   0.0, -1.0, 3.0, -1.0, 0.0, 0.0,
                   0.0, 0.0, -1.0, 3.0, -1.0, 0.0,
                   0.0, 0.5, 0.0, -1.0, 3.0, -1.0,
                   0.5, 0.0, 0.0, 0.0, -1.0, 3.0};
  double **A2 = mat_from_vector(a2, n2, n2);
  double b2[6] = {2.5, 1.5, 1.0, 1.0, 1.5, 2.5};
  double *x2;
  double r2[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};



  printf("*** JACOBI ***\n");
  x1 = vet_cria(n1);
  i = Jacobi(n1, A1, b1, x1, tol);
  printf("sistema 1: %d iter., [%.1f %.1f]\n", i, x1[0], x1[1]);
  x2 = vet_cria(n2);
  i = Jacobi(n2, A2, b2, x2, tol);
  printf("sistema 2: %d iter., [%.1f %.1f %.1f %.1f %.1f %.1f]\n", i, x2[0], x2[1], x2[2], x2[3], x2[4], x2[5]);


  printf("\n*** GAUSS-SEIDEL ***\n");
  x1 = vet_cria(n1);
  i = GaussSeidel(n1, A1, b1, x1, tol);
  printf("sistema 1: %d iter., [%.1f %.1f]\n", i, x1[0], x1[1]);
  x2 = vet_cria(n2);
  i = GaussSeidel(n2, A2, b2, x2, tol);
  printf("sistema 2: %d iter., [%.1f %.1f %.1f %.1f %.1f %.1f]\n", i, x2[0], x2[1], x2[2], x2[3], x2[4], x2[5]);

  
  printf("\n*** SOR ***\n");
  x1 = vet_cria(n1);
  i = SOR(n1, A1, b1, x1, tol, w);
  printf("sistema 1: %d iter., [%.1f %.1f]\n", i, x1[0], x1[1]);
  x2 = vet_cria(n2);
  i = SOR(n2, A2, b2, x2, tol, w);
  printf("sistema 2: %d iter., [%.1f %.1f %.1f %.1f %.1f %.1f]\n", i, x2[0], x2[1], x2[2], x2[3], x2[4], x2[5]);
  return 0;
}