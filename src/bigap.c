
#include <R.h>
#include <Rmath.h>

//kendall's tau
void getkendalltau(double *uv, int *nrow, double *rst) {
  long int ncon=0, ndis=0;
  int i, j;

  for (i = 0; i < (*nrow-1); i++) {
    for (j = (i+1); j < *nrow; j++) {
      if ((uv[(i*2)] <= uv[(j*2)] &&
	  uv[(i*2)+1] <= uv[(j*2)+1]) ||
	(uv[(i*2)] >= uv[(j*2)] &&
	 uv[(i*2)+1] >= uv[(j*2)+1])) {
	ncon++;
      } else {
	ndis++;
      }
    }
  }
  //printf("%i,%i\n", ncon, ndis);
  *rst = (double)(ncon-ndis)/(ncon+ndis);
}
