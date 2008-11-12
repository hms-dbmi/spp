#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"

#undef DEBUG 1

// dout is npos-length output array.
// n - number of positions in pos (and length of tc cont array)
// spos - starting position
void cdensum(int *n, double *pos, double *tc, double *spos, int *bw,int *dw, int *npos, double *dout)
{
  int i,j;
 
  double epos= *spos + ((double) *npos);
  double dbw=(double) *bw;
  for(i = 0; i< *n; i++) {
    // size of the window to which the contributions should be added
    int in=(int) (pos[i]- *spos);
    int ic=tc[i];
    int whs=(*dw)*(*bw)*ic;
    int ws=in-whs;
    int we=in+whs;
    if(ws<0) { ws=0; } 
    if(we>= *npos) { we= *npos -1; }
    
    for(j=ws;j<we;j++) {
      double beta=((double)(j-in))/dbw;
      dout[j]+=((double)ic)*exp(-0.5*beta*beta);
    }
  }
}
