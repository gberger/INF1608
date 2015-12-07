#include <stdio.h>
#include "gradconj.h"

int ConjugateGradient (int n, Sparse** A, double* b, double* x) {
  int k, i;
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


  sparse_multv(n, A, x, ax);
  vet_subtrai(n, b, ax, r);
  vet_copia(n, r, d);
  

  for (k = 0; k < n; k++) {
    if (vet_norma2(n, r) < 0.00000005) {
      break;
    }
    
    // Alpha_k = rTk . rK / dTk . A dk
    sparse_multv(n, A, d, Ad);
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

  return k;
}