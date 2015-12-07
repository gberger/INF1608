#ifndef SPARSE_H
#define SPARSE_H

struct sparse {
  int col;
  double val;
};
typedef struct sparse Sparse;

Sparse** sparse_cria(int n);

void sparse_libera(int n, Sparse** M);

// A[i][j]
double sparse_get(int i, int j, Sparse** A);

// w = M * v
void sparse_multv(int n, Sparse** M, double* v, double* w);

// C = A * B
Sparse** sparse_multm(int n, Sparse** A, Sparse** B);



double* vet_cria(int n);

// a + b = c
void vet_soma(int n, double* a, double* b, double* c);

// a - b = c
void vet_subtrai(int n, double* a, double* b, double* c);

// dst <- src
void vet_copia(int n, double* src, double* dst);

// w = sv
void vet_mults(int n, double* v, double s, double* w);

// v . w
double vet_dot(int n, double* v, double* w);

double vet_norma2 (int n, double* v);

#endif