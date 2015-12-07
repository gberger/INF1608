#include <stdio.h>
#include <math.h>
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


int ConjugateGradientPC (int n, Sparse** A, double* b, double* x, Sparse** M) {
  int k, i;
  double* ax = vet_cria(n);
  double* r = vet_cria(n);
  double* d = vet_cria(n);
  double* z = vet_cria(n);
  double p, pn;

  double* w = vet_cria(n);

  double alpha;
  double *ad = vet_cria(n);
  double *xn = vet_cria(n);
  double *rn = vet_cria(n);

  double *aw = vet_cria(n);
  double *pnpd = vet_cria(n);


  sparse_multv(n, A, x, ax);
  vet_subtrai(n, b, ax, r);
  vet_copia(n, r, d);


  z = vet_cria(n);
  ConjugateGradient(n, M, r, z);

  p = vet_dot(n, r, z);

  
  for (k = 0; k < n; k++) {
    if (vet_norma2(n, r) < 0.00000005 || sqrt(p) < 0.00000005) {
      break;
    }

    // w = A d
    sparse_multv(n, A, d, w);

    // alpha = p / wTd
    alpha = p / vet_dot(n, w, d);

    // x = x + alpha d
    vet_mults(n, d, alpha, ad);
    vet_soma(n, x, ad, xn);
    vet_copia(n, xn, x);

    // r = r - alpha w
    vet_mults(n, w, alpha, aw);
    vet_subtrai(n, r, aw, rn);
    vet_copia(n, rn, r);

    // z = Minv r
    // M z = r
    // usar ConjGrad para achar z
    z = vet_cria(n);
    ConjugateGradient(n, M, r, z);

    // pn = zT r
    pn = vet_dot(n, z, r);

    // d = z + p/pn d
    vet_mults(n, d, pn / p, pnpd);
    vet_soma(n, z, pnpd, d);

    // next
    p = pn;
  }

  return k;
}