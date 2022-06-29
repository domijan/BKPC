#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .C calls */
extern void getKPCs(void *, void *, void *, void *, void *, void *);
extern void mainbkpc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rkern(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"getKPCs",  (DL_FUNC) &getKPCs,   6},
  {"mainbkpc", (DL_FUNC) &mainbkpc, 19},
  {"rkern",    (DL_FUNC) &rkern,     6},
  {NULL, NULL, 0}
};

void R_init_BKPC(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}