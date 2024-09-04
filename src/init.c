/*
 *
 */

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "maxmatmul.h"

static const R_CallMethodDef callMethods[]  = {
  {"maxmatmulRC", (DL_FUNC) &maxmatmulRC, -1},
  {NULL, NULL, 0}
};

void R_init_maxstablePCA(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

