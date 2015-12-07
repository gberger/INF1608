#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sparse.h"
#include "gradconj.h"

int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

Sparse** cria_primeiro_sistema(int n) {
  int i;
  Sparse** A = sparse_cria(n);

  for (i = 0; i < n; i++) {
    A[i] = malloc(sizeof(Sparse) * 6);

    A[i][0].col = i;
    A[i][0].val = i+1;

    A[i][1].col = mod(i+1, n);
    A[i][1].val = 0.5;

    A[i][2].col = mod(i+2, n);
    A[i][2].val = 0.5;

    A[i][3].col = mod(i-1, n);
    A[i][3].val = 0.5;

    A[i][4].col = mod(i-2, n);
    A[i][4].val = 0.5;

    A[i][5].col = -1;
  }
  return A;
}

Sparse** cria_segundo_sistema(int n) {
  int i, k, col;
  Sparse** A = sparse_cria(n);

  for (i = 0; i < n; i++) {
    A[i] = malloc(sizeof(Sparse) * 8);

    A[i][0].col = i;
    A[i][0].val = i+1;

    A[i][1].col = mod(i+1, n);
    A[i][1].val = 0.5;

    A[i][2].col = mod(i+2, n);
    A[i][2].val = 0.5;

    A[i][3].col = mod(i-1, n);
    A[i][3].val = 0.5;

    A[i][4].col = mod(i-2, n);
    A[i][4].val = 0.5;

    if(i+1 < n/2) {
      k = 5;

      col = (i+1)/2;
      if(col != i && col != mod(i+1, n)&& col != mod(i+2, n)&& col != mod(i-1, n)&& col != mod(i-2, n)){
        A[i][k].col = col;
        A[i][k].val = 0.5;
        k++;        
      }

      col = (i+1)*2;
      if(i * 2 < n && col != i && col != mod(i+1, n)&& col != mod(i+2, n)&& col != mod(i-1, n)&& col != mod(i-2, n)){
        A[i][k].col = col;
        A[i][k].val = 0.5;
        k++;
      }

      A[i][k].col = -1;
    } else {
      A[i][5].col = -1;
    }
  }
  return A;
}

Sparse** precond_jacobi(int n, Sparse** A, double* b) {
  int i;

  double* b_ = vet_cria(n);
  Sparse** Minv = sparse_cria(n);
  Sparse** C;

  for(i = 0; i < n; i++) {
    Minv[i] = malloc(sizeof(Sparse) * 2);

    Minv[i][0].col = i;
    Minv[i][0].val = 1.0/sparse_get(i, i, A);

    Minv[i][1].col = -1;
  }


  C = sparse_multm(n, Minv, A);

  sparse_multv(n, Minv, b, b_);
  vet_copia(n, b_, b);

  return C;
}

Sparse** precond_ssor(int n, Sparse** A, double* b, double w) {
  int i, j, k;

  double* b_ = vet_cria(n);

  Sparse** C;
  Sparse** M1 = sparse_cria(n);
  Sparse** M2 = sparse_cria(n);
  Sparse** M3 = sparse_cria(n);
  Sparse** Ma;
  Sparse** Mb;
  Sparse** Minv;

  for(i = 0; i < n; i++) {
    M1[i] = malloc(sizeof(Sparse) * n);

    M1[i][0].col = i;
    M1[i][0].val = sparse_get(i, i, A);

    k = 1;
    for(j = 0; j < n; j++) {
      if(A[i][j].col == -1) break;
      if(A[i][j].col < i) {
        M1[i][k].col = A[i][j].col;
        M1[i][k].val = A[i][j].val * w;
        k++;
      }
    }

    M1[i][k].col = -1;
  }

  for(i = 0; i < n; i++) {
    M2[i] = malloc(sizeof(Sparse) * 2);

    M2[i][0].col = i;
    M2[i][0].val = 1.0/sparse_get(i, i, A);

    M2[i][1].col = -1;
  }

  for(i = 0; i < n; i++) {
    M3[i] = malloc(sizeof(Sparse) * n);

    M3[i][0].col = i;
    M3[i][0].val = sparse_get(i, i, A);

    k = 1;
    for(j = 0; j < n; j++) {
      if(A[i][j].col == -1) break;
      if(A[i][j].col > i) {
        M3[i][k].col = A[i][j].col;
        M3[i][k].val = A[i][j].val * w;
        k++;
      }
    }

    M3[i][k].col = -1;
  }

  Ma = sparse_multm(n, M1, M2);
  Mb = sparse_multm(n, Ma, M3);
  // M = sparse_inv(n, Mb);
  Minv = Mb;

  C = sparse_multm(n, Minv, A);

  sparse_multv(n, Minv, b, b_);
  vet_copia(n, b_, b);

  return C;
}

double* cria_x(int n) {
  int i;
  double* x = vet_cria(n);
  for(i = 0; i < n; i++) {
    x[i] = 1;
  }
  return x;
}

void testa(int n, Sparse** A, double* b, double* xbarra, double* xsol) {
  int i, iter;
  double dif;

  iter = ConjugateGradient(n, A, b, xbarra);

  for(i = 0; i < n; i++) {
    dif += fabs((xbarra[i] - xsol[i])/xsol[i]);
  }
  dif = dif / n;
  printf("Erro: %f%%\n", dif*100);
  printf("Iteracoes: %d\n", iter);
  printf("\n");
}

int main(void) {
  int i;
  int n = 1000;
  Sparse** A;
  double* x;
  double* b;
  double* xbarra;


  printf("**** N = %d\n\n", n);

  printf("Sem pre-cond\n");
  A = cria_segundo_sistema(n);
  x = cria_x(n);  
  b = vet_cria(n);
  xbarra = calloc(n, sizeof(double));
  sparse_multv(n, A, x, b);
  testa(n, A, b, xbarra, x);


  printf("Jacobi\n");
  A = cria_segundo_sistema(n);
  x = cria_x(n);  
  b = vet_cria(n);
  xbarra = calloc(n, sizeof(double));
  sparse_multv(n, A, x, b);
  A = precond_jacobi(n, A, b);
  testa(n, A, b, xbarra, x);


  printf("Gauss-Seidel\n");
  A = cria_segundo_sistema(n);
  x = cria_x(n);  
  b = vet_cria(n);
  xbarra = calloc(n, sizeof(double));
  sparse_multv(n, A, x, b);
  A = precond_ssor(n, A, b, 1.0);
  testa(n, A, b, xbarra, x);
  

  printf("SSOR, w=1.1\n");
  A = cria_segundo_sistema(n);
  x = cria_x(n);  
  b = vet_cria(n);
  xbarra = calloc(n, sizeof(double));
  sparse_multv(n, A, x, b);
  A = precond_ssor(n, A, b, 1.1);
  testa(n, A, b, xbarra, x);

  return 0;
}
