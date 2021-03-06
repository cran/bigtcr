#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void getkendalltau(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"getkendalltau", (DL_FUNC) &getkendalltau, 3},
  {NULL, NULL, 0}
};

void R_init_bigtcr(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
