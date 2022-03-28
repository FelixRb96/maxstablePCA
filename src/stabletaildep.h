/*
 *
 */

#ifndef stabletaildep_H
#define stabletaildep_H

SEXP stabletaildepRC(SEXP input, SEXP data, SEXP a, SEXP nrowdat, SEXP ncoldat);
SEXP stabletaildep_vecRC(SEXP input, SEXP data, SEXP a, SEXP nrowdat, SEXP ncoldat, SEXP nevals); 


double stabletaildep(double *input, double *data, double a, int nrowdat, int ncoldat);
void stabletaildep_vec(double *res, double *input, double *data, double a, int nrowdat, int ncoldat, int nevals); 


#endif
