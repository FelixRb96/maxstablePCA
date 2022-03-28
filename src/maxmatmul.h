/*
 *
 */

#ifndef maxmatmul_H
#define maxmatmul_H

SEXP maxmatmulRC(SEXP A, SEXP B, SEXP n, SEXP l, SEXP k);

void maxmatmulC(double *res, double *A, double *B, int n, int l, int k);

#endif
