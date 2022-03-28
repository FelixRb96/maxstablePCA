/*
 *
 */

#include <R.h>
#include <Rinternals.h>

#include "stabletaildep.h"

// macro to index array like matrix
#define IJ(arr, i, j, rowlen) arr[(i * rowlen) + j]


/*
 *
 */
SEXP stabletaildepRC(SEXP input, SEXP data, SEXP a, SEXP nrowdat, SEXP ncoldat) 
{

  double *dinput = REAL(input);
  double *ddata = REAL(data); 
  double aa = asReal(a);

  int inrowdat = asInteger(nrowdat);
  int incoldat = asInteger(ncoldat);

  double res = stabletaildep(dinput, ddata, aa, inrowdat, incoldat);

  SEXP result = PROTECT(allocVector(REALSXP, 1));
  REAL(result)[0] = res;

  UNPROTECT(1);
  return result;

}

/*
 *
 */
SEXP stabletaildep_vecRC(SEXP input, SEXP data, SEXP a, SEXP nrowdat, SEXP ncoldat, SEXP nevals) 
{

  double *dinput = REAL(input);
  double *ddata = REAL(data); 
  double aa = asReal(a);

  int inrowdat = asInteger(nrowdat);
  int incoldat = asInteger(ncoldat);
  int inevals = asInteger(nevals);

  double *res = malloc(sizeof(double) * inevals); 
  stabletaildep_vec(res, dinput, ddata, aa, inrowdat, incoldat, inevals);

  SEXP result = PROTECT(allocVector(REALSXP, inevals));
  for(int i = 0; i < inevals; i++) {
    REAL(result)[i] = res[i];
  }

  UNPROTECT(1);
  free(res); 

  return result;

}

/*
 *
 */
double stabletaildep(double *input, double *data, double a, int nrowdat, int ncoldat) 
{

  double res = 0;

  for(int i = 0; i < nrowdat; i++) {
    double rowsum = 0;
    double maxprod = 0; 
    for(int j = 0; j < ncoldat; j++) {
      rowsum += IJ(data, i, j, ncoldat);

      double tmp = IJ(data, i, j, ncoldat) * IJ(input, 0, j, ncoldat); 
      maxprod = (maxprod >= tmp)? maxprod : tmp; 

    }
    res += (rowsum > nrowdat / a)? maxprod / rowsum : 0;
  }
  return res / a;
}

/*
 * C implementation of the empirical stable tail dependence for a dataset with unit Fréchet margins. 
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
void stabletaildep_vec(double *res, double *input, double *data, double a, int nrowdat, int ncoldat, int nevals)
{

  double *tmp = malloc(sizeof(double) * ncoldat);

  for(int i = 0; i < nevals; i++) {
    for(int j = 0; j < ncoldat; j++) {
      tmp[j] = IJ(input, i, j, ncoldat);
    }
    res[i] = stabletaildep(tmp, data, a, nrowdat, ncoldat);
  }
  free(tmp);

  /*
  for(int k = 0; k < nevals; k++) {
    res[k] = 0;
  }

  for(int i = 0; i < nrowdat; i++) {
    double rowsum = 0;
    for(int k = 0; k < nevals; k++) {
      double maxprod = 0; 

      for(int j = 0; j < ncoldat; j++) {
        rowsum += IJ(data, i, j, ncoldat);

        double tmp = IJ(data, i, j, ncoldat) * IJ(input, k, j, ncoldat); 
        maxprod = (maxprod >= tmp)? maxprod : tmp; 

      }
      res[k] += (rowsum > nrowdat / a)? maxprod / rowsum : 0;
    }
  }
  for(int k = 0; k < nevals; k++) {
    res[k] = res[k] / a; 
  }
  */
}
