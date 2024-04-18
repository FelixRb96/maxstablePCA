/*
 *
 */

#include <R.h>
#include <Rinternals.h>

#include "maxmatmul.h"

// macro to index array like matrix
#define IJ(arr, i, j, rowlen) arr[(i * rowlen) + j]

/*
 * R callable function for the max matrix multiplication
 *
 * @params
 * A: SEXP containing array of length n*l representing matrix A read rowwise.
 * B: SEXP containing array of length l*k representing matrix B read rowwise.
 * n: SEXP containing non-negative integer n.
 * l: SEXP containing non-negative integer l.
 * k: SEXP containing non-negative integer k.
 *
 * @returns 
 * SEXP with array containing the entries of the resultung product rowwise. 
 */
SEXP maxmatmulRC(SEXP A, SEXP B, SEXP n, SEXP l, SEXP k) 
{

  double *dA = REAL(A);
  double *dB = REAL(B); 

  int in = asInteger(n);
  int il = asInteger(l);
  int ik = asInteger(k);

  double *res = malloc(sizeof(double) * in * ik);
  maxmatmulC(res, dA, dB, in, il, ik);

  SEXP result = PROTECT(allocVector(REALSXP, in * ik));
  for(int i = 0; i < in * ik; i++) {
    REAL(result)[i] = res[i];
  }

  UNPROTECT(1);
  free(res);

  return result;
}


/*
 * C implementation of the max matrix multiplication. 
 *
 * @params
 *
 * res: double pointer containing array of length n*k representing matrix A read rowwise, where result will be stored to. 
 * A: double pointer containing array of length n*l representing matrix A read rowwise.
 * B: double pointer containing array of length l*k representing matrix B read rowwise.
 * n: int containing non-negative integer n.
 * l: int containing non-negative integer l.
 * k: int containing non-negative integer k.
 *
 * @returns
 *
 * void 
 */
void maxmatmulC(double *res, double *A, double *B, int n, int l, int k) 
{
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < k; j++) {

      // inner loop for determining the maximum
      double tmpres = 0;
      for(int m = 0; m < l; m++) {
        double compareval = IJ(A, i, m, l) * IJ(B, m, j, k);
        tmpres = (tmpres > compareval) ? tmpres : compareval;
      }

      IJ(res, i, j, k) = tmpres;
    }
  }
}
