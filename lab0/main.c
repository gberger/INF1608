#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
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

void test_mat_transposta() {
  double a[9] = {1.0, 2.0, 3.0, 
                 4.0, 5.0, 6.0,
                 7.0, 8.0, 9.0};
  double** A = mat_from_vector(a, 3, 3);

  double** TA = mat_cria(3, 3);
  mat_transposta(3, 3, A, TA);

  assert(TA[0][0] == 1.0);
  assert(TA[0][1] == 4.0);
  assert(TA[0][2] == 7.0);
  assert(TA[1][0] == 2.0);
  assert(TA[1][1] == 5.0);
  assert(TA[1][2] == 8.0);
  assert(TA[2][0] == 3.0);
  assert(TA[2][1] == 6.0);
  assert(TA[2][2] == 9.0);
}

void test_mat_libera() {
  int i;
  double** A = mat_cria(100, 200);
  mat_libera(100, A);
}

void test_mat_multv() {
  double a[6] = {1.0, 2.0, 3.0, 
                 4.0, 5.0, 6.0};
  double** A = mat_from_vector(a, 2, 3);

  double v[3] = {10.0, 20.0, 30.0};
  double w[2];

  mat_multv(2, 3, A, v, w);

  assert(w[0] == 140.0);
  assert(w[1] == 320.0);
}

void test_mat_multm() {
  double a[6] = {1.0, 2.0, 3.0, 
                 4.0, 5.0, 6.0};
  double** A = mat_from_vector(a, 2, 3);

  double b[6] = {10.0, 20.0, 
                 30.0, 40.0, 
                 50.0, 60.0};
  double** B = mat_from_vector(a, 3, 2);

  double** C = mat_cria(2, 2);

  mat_multm(2, 3, 2, A, B, C);

  assert(C[0][0] == 22.0);
  assert(C[0][1] == 28.0);
  assert(C[1][0] == 49.0);
  assert(C[1][1] == 64.0);
}


int main(void) {
  test_mat_transposta();
  test_mat_libera();
  test_mat_multv();
  test_mat_multm();
  return 0;
}
