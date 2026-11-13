#include <vector>
#include <string.h>
#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include <cstdlib>
#include <algorithm>

extern "C" {
//#include "R.h"
//#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;

namespace {
  // sorted sparse histogram: positions pos[0..n) strictly increasing, counts cnt[]
  struct SparseHist {
    const int* pos;
    const int* cnt;
    int n;
    int ptr;
    SparseHist(const int* p, const int* c, int n_) : pos(p), cnt(c), n(n_), ptr(0) {}
    // value at index j; j must not decrease between calls (shared cursor)
    inline int get(int j) {
      while(ptr<n && pos[ptr]<j) { ptr++; }
      return (ptr<n && pos[ptr]==j) ? cnt[ptr] : 0;
    }
    // first position >= j after advancing the cursor, or -1 when none
    inline long next_pos(int j) {
      while(ptr<n && pos[ptr]<j) { ptr++; }
      return ptr<n ? (long) pos[ptr] : -1L;
    }
    // sum of counts over indices lo..hi inclusive (independent of the cursor)
    int range_sum(long lo, long hi) const {
      if(lo<0) { lo=0; }
      if(lo>hi) { return 0; }
      int a=(int)(lower_bound(pos,pos+n,(int)lo)-pos);
      int s=0;
      for(int k=a;k<n && pos[k]<=hi;k++) { s+=cnt[k]; }
      return s;
    }
  };
}

//#define DEBUG 1

extern "C" {

  /************************************************************************/
  /*
   * lwcc - calculate local window cross-correlation
   */

  SEXP spp_lwcc(SEXP x_R, // positive strand hist 
	    SEXP y_R, // negative strand hist of the same length
	    SEXP osize_R,       // outer boundary distance
	    SEXP isize_R,        // inner boundary distance
	    SEXP return_peaks_R, // whether all correlation values, or just peaks should be returned
	    SEXP min_peak_dist_R, // distance between closest peaks
	    SEXP min_peak_val_R, // min peak threshold
	    SEXP tag_weight_R,  // tag weight
	    SEXP bg_subtract_R, // a flag whether do background subtractio
	    SEXP bgp_R, // optional background hist for positive strand
	    SEXP bgn_R, // optional background hist for negative strand
	    SEXP bg_wsize_R, // window size for the background counts
	    SEXP bg_weight_R, // optional weighting for the background tags, must compensate for window size difference (including is cutout)
	    SEXP round_up_R // whether to round up fractional signal tag counts
	    )
  {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    
    int is=INTEGER(isize_R)[0];
    int os=INTEGER(osize_R)[0];
    double rs=((double)(2*os+1));
    int* x=INTEGER(x_R);
    int* y=INTEGER(y_R);
    int n_x=LENGTH(x_R);

    // background-related 
    int* bgp=INTEGER(bgp_R);
    int* bgn=INTEGER(bgn_R);
    int bg_whs=INTEGER(bg_wsize_R)[0];

    int return_peaks=*(INTEGER(return_peaks_R));
    double min_peak_val=*(REAL(min_peak_val_R));
    int min_peak_dist=*(INTEGER(min_peak_dist_R));
    double tag_weight=*(REAL(tag_weight_R));

    const int round_up=*(INTEGER(round_up_R));
    const int bg_subtract=*(INTEGER(bg_subtract_R));
    const double bg_weight=*(REAL(bg_weight_R));

    int i; // point at which the value is being calculated
    int start=os;
    int end=n_x-os-1;

    // bg tag counts within bg window
    int bg_pn1=0;
    int bg_nn1=0;
    int bg_pn2=0;
    int bg_nn2=0;


  
    // illustration for counting:
    //
    // 012345678901234567890123456789012
    // ==========------|------==========
    //
    //  osize=16; isize=6; 


    SEXP nv;
    PROTECT(nv=allocVector(REALSXP,n_x)); //putting declaration also here to avoid[-Wmaybe-uninitialized]

    double *d_nv;
    vector<int> ppos;
    vector<double> pval;
    if(!return_peaks) {
      PROTECT(nv=allocVector(REALSXP,n_x)); 
      d_nv=REAL(nv);
      for(int i=0;i<n_x;i++) {
	d_nv[i]=0;
      }
    }

#ifdef DEBUG  
    Rprintf("start=%d end=%d tag_weight=%f\n", start,end,tag_weight);
    Rprintf("x[1]=%d x[2]=%d y[1]=%d y[2]=%d\n",x[1],x[2],y[1],y[2]);
#endif

    int lpp=-1; // last peak position
    double lpv=-1e3; // last peak value
    
    double ppv=-1e3; // last value
    double pppv=-11e-3; // value before last

    int pn1,pn2,nn1,nn2;

    
    if(bg_subtract) {
      // pre-initialize background tag counts, 
      for(int i=0;i<bg_whs;i++) {
	if(i<n_x) {
	  bg_pn2+=bgp[i];
	  bg_nn2+=bgn[i];
	}
      }
    }


    for(i=0;i<end;i++) {
#ifdef DEBUG  
      //Rprintf("i=%d ", i);
#endif
      
      if(bg_subtract) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

      if(i >= start) {
	// update counts, taking into account masked out regions
	pn1=pn2=nn1=nn2=0;
	
	for(int k=0;k<=(os-is);k++) {
	  int xp1=x[i-os+k];
	  int xp2=x[i+os-k];
	  int xn1=y[i+os-k];
	  int xn2=y[i-os+k];

	  if(xp1!=-1 && xn1!=-1) {
	    pn1+=xp1;
	    nn1+=xn1;
	  }
	  if(xp2!=-1 && xn2!=-1) {
	    pn2+=xp2;
	    nn2+=xn2;
	  }
	}
      
	// calculate the means
	double mp=((double)(pn1+pn2))/rs;
	double mn=((double)(nn1+nn2))/rs;
#ifdef DEBUG  
	Rprintf("mp=%f mn=%f\n",mp,mn);
#endif
	// calculate correlation
	double varp=0;
	double varn=0;
	double num=0;
	double val=-1e3;
	if((mp>0) & (mn>0)) {
	  for(int k=0;k<=(os-is);k++) {
	    int xp1=x[i-os+k];
	    int xp2=x[i+os-k];
	    int xn1=y[i+os-k];
	    int xn2=y[i-os+k];

	    
	    if(xp1!=-1 && xn1!=-1) {  
	      double nnp1=((double) xp1)-mp;
	      double nnn1=((double) xn1)-mn;
	      num+=nnp1*nnn1;
	      varp+=nnp1*nnp1;
	      varn+=nnn1*nnn1;
	    }
	    
	    if(xp2!=-1 && xn2!=-1) {
	      double nnp2=((double) xp2)-mp;
	      double nnn2=((double) xn2)-mn;
	      num+=nnp2*nnn2;
	      varp+=nnp2*nnp2;
	      varn+=nnn2*nnn2;
	    }

	  }
	  double tagw;
	  double spn1=((double)pn1)*tag_weight;
	  double snn1=((double)nn1)*tag_weight;
	  double spn2=((double)pn2)*tag_weight;
	  double snn2=((double)nn2)*tag_weight;
	  if(round_up) {
	    if(pn1>0 && spn1<1) { spn1=1.0; }
	    //if(pn2>0 && spn2<1) { spn2=1.0; }
	    if(nn1>0 && snn1<1) { snn1=1.0; }
	    //if(nn2>0 && snn2<1) { snn2=1.0; }
	  }

	  if(bg_subtract) {
	    spn1-=((double)bg_pn1)*bg_weight;
	    snn1-=((double)bg_nn2)*bg_weight;
	    spn2-=((double)bg_pn2)*bg_weight;
	    snn2-=((double)bg_nn1)*bg_weight;

	    if(spn2<0) spn2=0;
	    if(snn2<0) snn2=0;
	    
	    if(spn1>0 && snn1>0) {
	      tagw=(2.0*sqrt(spn1*snn1)-(spn2+snn2+1.0));
	    } else {
	      tagw=-(spn2+snn2+1.0);
	    }
	    //cout<<"bg_pn1="<<bg_pn1<<"; bg_pn2="<<bg_pn2<<"; bg_nn1="<<bg_nn1<<"; bg_nn2="<<bg_nn2<<endl;
	  } else {
	    tagw=2.0*sqrt(spn1*snn1)-(spn2+snn2);
	  }

	  if(tagw<0) {
	    val=0.0; 
	  } else {
	    if(num==0.0) {
	      val=0;
	    } else {
	      val=num/(sqrt(varp*varn));
	    }
	    val=val*sqrt(tagw) + tagw;

	  }
	  //cout<<"val="<<val<<endl;

#ifdef DEBUG  
        Rprintf("pn1=%d pn2=%d nn1=%d nn2=%d tag.weight=%f tagw=%f\n",pn1,pn2,nn1,nn2,tag_weight,tagw);
	Rprintf("tagw=%f varp=%f varn=%f num=%f cor=%f val=%f\n",tagw,varp,varn,num,num/sqrt(varp*varn),val);
#endif
	}


	
	if(return_peaks) {
	  // determine if previous position was a peak
	  if(ppv>min_peak_val && ppv>val && ppv>pppv) {
	    if(lpp>0 && (i-lpp+1)>min_peak_dist) {
	      // record previous peak position
	      ppos.push_back(lpp);
	      pval.push_back(lpv);
#ifdef DEBUG  
	      Rprintf("recording peak x=%d y=%f d=%d\n",lpp,lpv,(i-lpp));
#endif	    
	      lpp=i-1; lpv=ppv;
#ifdef DEBUG  
	      Rprintf("updated peak to x=%d y=%f\n",lpp,lpv);
#endif	    
	    } else {
	      if(ppv>lpv) {
		// update last peak positions
#ifdef DEBUG  
		Rprintf("skipping peak x=%d y=%f d=%d in favor of x=%d y=%f\n",lpp,lpv,(i-lpp),i-1,ppv);
#endif
		lpp=i-1; lpv=ppv;
	      }
	    }
	  }

	  // update previous values
	  if(val!=ppv) {
	    pppv=ppv; ppv=val;
	  }
	} else {
	  d_nv[i]=val;
	}
      }
    }

    if(return_peaks) {
      // record last position
      if(lpp>0) {
#ifdef DEBUG  
	Rprintf("recording last peak x=%d y=%f\n",lpp,lpv);
#endif
	ppos.push_back(lpp);
	pval.push_back(lpv);
      }

      SEXP rpp_R,rpv_R;
      PROTECT(rpp_R=allocVector(INTSXP,ppos.size())); 
      PROTECT(rpv_R=allocVector(REALSXP,ppos.size())); 
      int* rpp=INTEGER(rpp_R);
      double* rpv=REAL(rpv_R);

	  int iloop=ppos.size(); //creating dummy variable for the loop below to avoid [-Wsign-compare]
      for(int i=0;i<iloop;i++) {
      //for(int i=0;i<ppos.size();i++) {
	rpp[i]=ppos[i];
	rpv[i]=pval[i];
      }
    
      SEXP ans_R, names_R;
      PROTECT(names_R = allocVector(STRSXP, 2));
      SET_STRING_ELT(names_R, 0, mkChar("x"));
      SET_STRING_ELT(names_R, 1, mkChar("v"));
    
      PROTECT(ans_R = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(ans_R, 0, rpp_R);
      SET_VECTOR_ELT(ans_R, 1, rpv_R);
      setAttrib(ans_R, R_NamesSymbol, names_R);
  
      UNPROTECT(4);
      return(ans_R);
    } else {
      UNPROTECT(1);
      return(nv);
    }

  }



  /************************************************************************/
  /*
   * wtd - window tag difference implementation
   */


  /*
   * wtd on sparse histograms, ignore_masking path with return_peaks=1 only.
   * Each histogram is given as sorted positions (0-based indices into the dense
   * histogram of length n_x) and counts. The window counts only change at steps
   * where a nonzero position enters or leaves a window, so val is constant
   * between those steps and the peak state is replayed only at the steps that
   * can change it.
   */
  SEXP spp_wtd_sparse(SEXP xpos_R, SEXP xcnt_R, SEXP ypos_R, SEXP ycnt_R,
	   SEXP n_x_R,
	   SEXP wsize_R,
	   SEXP min_peak_dist_R,
	   SEXP min_peak_val_R,
	   SEXP direct_count_R,
	   SEXP tag_weight_R,
	   SEXP bg_subtract_R,
	   SEXP bgxpos_R, SEXP bgxcnt_R, SEXP bgypos_R, SEXP bgycnt_R,
	   SEXP bg_wsize_R,
	   SEXP bg_weight_R,
	   SEXP round_up_R)
  {
    const int whs=INTEGER(wsize_R)[0];
    const int n_x=INTEGER(n_x_R)[0];
    const int bg_whs=INTEGER(bg_wsize_R)[0];
    const int direct_count=*(INTEGER(direct_count_R));
    const double min_peak_val=*(REAL(min_peak_val_R));
    const int min_peak_dist=*(INTEGER(min_peak_dist_R));
    const double tag_weight=*(REAL(tag_weight_R));
    const int round_up=*(INTEGER(round_up_R));
    const int bg_subtract=*(INTEGER(bg_subtract_R));
    const double bg_weight=*(REAL(bg_weight_R));

    SparseHist hx(INTEGER(xpos_R),INTEGER(xcnt_R),LENGTH(xpos_R));
    SparseHist hy(INTEGER(ypos_R),INTEGER(ycnt_R),LENGTH(ypos_R));
    SparseHist hbp(INTEGER(bgxpos_R),INTEGER(bgxcnt_R),bg_subtract ? LENGTH(bgxpos_R) : 0);
    SparseHist hbn(INTEGER(bgypos_R),INTEGER(bgycnt_R),bg_subtract ? LENGTH(bgypos_R) : 0);

    const int start=whs+1;
    const int end=n_x-whs-1;

    vector<int> ppos;
    vector<double> pval;

    if(start<end) {
      int pn1=hx.range_sum(0,whs-1);
      int pn2=hx.range_sum(whs+1,2L*whs);
      int nn1=hy.range_sum(0,whs-1);
      int nn2=hy.range_sum(whs+1,2L*whs);

      int bg_pn1=0, bg_pn2=0, bg_nn1=0, bg_nn2=0;
      if(bg_subtract) {
	// state of the dense loop after its pre-roll over i=0..start-1
	long lo1=(long)start-bg_whs-1;
	bg_pn1=hbp.range_sum(lo1,start-1);
	bg_nn1=hbn.range_sum(lo1,start-1);
	long hi0=(long)(bg_whs<n_x ? bg_whs : n_x)-1;
	long hi3=(long)start-1+bg_whs;
	if(hi3>n_x-1) { hi3=n_x-1; }
	bg_pn2=hbp.range_sum(0,hi0)-hbp.range_sum(0,(long)start-2)+hbp.range_sum(bg_whs,hi3);
	bg_nn2=hbn.range_sum(0,hi0)-hbn.range_sum(0,(long)start-2)+hbn.range_sum(bg_whs,hi3);
      }

      // one cursor per (histogram, offset): the index read at step i is i+c
      struct Acc { SparseHist* h; int c; };
      Acc acc[14];
      int nacc=0;
      SparseHist hxs[3]={hx,hx,hx}, hys[3]={hy,hy,hy};
      const int cmain[3]={-1,-whs-1,whs};
      for(int k=0;k<3;k++) { acc[nacc].h=&hxs[k]; acc[nacc].c=cmain[k]; nacc++; }
      for(int k=0;k<3;k++) { acc[nacc].h=&hys[k]; acc[nacc].c=cmain[k]; nacc++; }
      SparseHist hbps[4]={hbp,hbp,hbp,hbp}, hbns[4]={hbn,hbn,hbn,hbn};
      const int cbg[4]={-bg_whs-1,0,-1,bg_whs};
      if(bg_subtract) {
	for(int k=0;k<4;k++) { acc[nacc].h=&hbps[k]; acc[nacc].c=cbg[k]; nacc++; }
	for(int k=0;k<4;k++) { acc[nacc].h=&hbns[k]; acc[nacc].c=cbg[k]; nacc++; }
      }

      int lpp=-1;
      double lpv=-1000;
      double ppv=-1000;
      int ppl=-1;
      double pppv=-1000;

      int i=start;
      while(i<end) {
	if(bg_subtract) {
	  bg_pn1-=hbps[0].get(i-bg_whs-1);
	  bg_nn1-=hbns[0].get(i-bg_whs-1);
	  bg_pn1+=hbps[1].get(i);
	  bg_nn1+=hbns[1].get(i);
	  bg_pn2-=hbps[2].get(i-1);
	  bg_nn2-=hbns[2].get(i-1);
	  bg_pn2+=hbps[3].get(i+bg_whs);
	  bg_nn2+=hbns[3].get(i+bg_whs);
	}

	pn1+=hxs[0].get(i-1)-hxs[1].get(i-whs-1);
	pn2+=hxs[2].get(i+whs)-hx.get(i-1);
	nn1+=hys[0].get(i-1)-hys[1].get(i-whs-1);
	nn2+=hys[2].get(i+whs)-hy.get(i-1);

	double val;
	double spn1=((double)pn1)*tag_weight;
	double snn1=((double)nn1)*tag_weight;
	double spn2=((double)pn2)*tag_weight;
	double snn2=((double)nn2)*tag_weight;
	if(round_up) {
	  if(pn1>0 && spn1<1) { spn1=1.0; }
	  if(nn2>0 && snn2<1) { snn2=1.0; }
	}

	if(direct_count) {
	  val=spn1+snn2;
	  if(round_up && val<1) {
	    val=1.0;
	  }
	  if(bg_subtract) {
	    val-=((double) (bg_pn1+bg_nn2))*bg_weight;
	  }
	} else {
	  if(bg_subtract) {
	    spn1-=((double)bg_pn1)*bg_weight;
	    snn1-=((double)bg_nn1)*bg_weight;
	    spn2-=((double)bg_pn2)*bg_weight;
	    snn2-=((double)bg_nn2)*bg_weight;

	    if(spn2<0) spn2=0;
	    if(snn1<0) snn1=0;

	    if(spn1>0 && snn2>0) {
	      val=(2.0*sqrt(spn1*snn2)-(spn2+snn1+1.0));
	    } else {
	      val=-(spn2+snn1+1.0);
	    }
	  } else {
	    val=2.0*sqrt(spn1*snn2)-(spn2+snn1+tag_weight);
	  }
	}

	if(ppv>min_peak_val && ppv>val && ppv>pppv) {
	  if(lpp>0 && (i-lpp+1)>min_peak_dist) {
	    ppos.push_back(lpp);
	    pval.push_back(lpv);
	    if(ppl!=-1 && ppl!=i-1) {
	      lpp=(int) round((ppl+i-1)/2);
	    } else {
	      lpp=i-1;
	    }
	    lpv=ppv;
	  } else {
	    if(ppv>lpv) {
	      if(ppl!=-1 && ppl!=i-1) {
		lpp=(int) round((ppl+i-1)/2);
	      } else {
		lpp=i-1;
	      }
	      lpv=ppv;
	    }
	  }
	}

	if(val!=ppv) {
	  pppv=ppv; ppv=val; ppl=i;
	}

	// next step at which any window count can change; in between val is
	// constant, and with ppv==val the peak test and the state update are no-ops
	long nxt=-1;
	for(int a=0;a<nacc;a++) {
	  long p=acc[a].h->next_pos(i+1+acc[a].c);
	  if(p>=0) {
	    long cand=p-acc[a].c;
	    if(nxt<0 || cand<nxt) { nxt=cand; }
	  }
	}
	if(nxt<0 || nxt>=end) { break; }
	i=(int) nxt;
      }

      if(lpp>0) {
	ppos.push_back(lpp);
	pval.push_back(lpv);
      }
    }

    SEXP rpp_R,rpv_R;
    PROTECT(rpp_R=allocVector(INTSXP,ppos.size()));
    PROTECT(rpv_R=allocVector(REALSXP,ppos.size()));
    int* rpp=INTEGER(rpp_R);
    double* rpv=REAL(rpv_R);
    int iloop=ppos.size();
    for(int k=0;k<iloop;k++) {
      rpp[k]=ppos[k];
      rpv[k]=pval[k];
    }

    SEXP ans_R, names_R;
    PROTECT(names_R = allocVector(STRSXP, 2));
    SET_STRING_ELT(names_R, 0, mkChar("x"));
    SET_STRING_ELT(names_R, 1, mkChar("v"));
    PROTECT(ans_R = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans_R, 0, rpp_R);
    SET_VECTOR_ELT(ans_R, 1, rpv_R);
    setAttrib(ans_R, R_NamesSymbol, names_R);
    UNPROTECT(4);
    return(ans_R);
  }

  SEXP spp_wtd(SEXP x_R, // positive strand hist 
	   SEXP y_R, // negative strand hist of the same length
	   SEXP wsize_R,       // outer boundary distance
	   SEXP return_peaks_R, // whether all correlation values, or just peaks should be returned
	   SEXP min_peak_dist_R, // distance between closest peaks
	   SEXP min_peak_val_R, // min peak threshold
	   SEXP direct_count_R, // whether tag weighting should not be done
	   SEXP tag_weight_R,  // tag weight
	   SEXP ignore_masking_R,  // whether to ignore masked regions
	   SEXP bg_subtract_R, // a flag whether do background subtractio
	   SEXP bgp_R, // optional background hist for positive strand
	   SEXP bgn_R, // optional background hist for negative strand
	   SEXP bg_wsize_R, // window size for the background counts
	   SEXP bg_weight_R, // optional weighting for the background tags, must compensate for window size difference
	   SEXP round_up_R // whether to round up fractional signal tag counts
	   )
  {

#ifdef DEBUG  
    Rprintf("start\n");
#endif
    
    int whs=INTEGER(wsize_R)[0];
    int* x=INTEGER(x_R);
    int* y=INTEGER(y_R);
    int n_x=LENGTH(x_R);
    
    // background-related 
    int* bgp=INTEGER(bgp_R);
    int* bgn=INTEGER(bgn_R);
    int bg_whs=INTEGER(bg_wsize_R)[0];
    

    const int return_peaks=*(INTEGER(return_peaks_R));
    const int direct_count=*(INTEGER(direct_count_R));
    const int ignore_masking=*(INTEGER(ignore_masking_R));
    const double min_peak_val=*(REAL(min_peak_val_R));
    const int min_peak_dist=*(INTEGER(min_peak_dist_R));
    const double tag_weight=*(REAL(tag_weight_R));
    
    const int round_up=*(INTEGER(round_up_R));
    const int bg_subtract=*(INTEGER(bg_subtract_R));
    const double bg_weight=*(REAL(bg_weight_R));
    
    int i; // point at which the value is being calculated
    int start=whs+1;
    int end=n_x-whs-1;

    // tag counts to calculate the means
    int pn1=0;
    int pn2=0;
    int nn1=0;
    int nn2=0;

    // bg tag counts within bg window
    int bg_pn1=0;
    int bg_pn2=0;
    int bg_nn1=0;
    int bg_nn2=0;
    
    SEXP nv;
    PROTECT(nv=allocVector(REALSXP,n_x)); //putting declaration also here to avoid[-Wmaybe-uninitialized]
    double *d_nv;
    vector<int> ppos;
    vector<double> pval;
    if(!return_peaks) {
      PROTECT(nv=allocVector(REALSXP,n_x)); 
      d_nv=REAL(nv);
      for(int i=0;i<n_x;i++) {
	d_nv[i]=0;
      }
    }

#ifdef DEBUG  
    Rprintf("whs=%d start=%d end=%d tag_weight=%f ignore_masing=%d\n", whs, start,end,tag_weight,ignore_masking);
    Rprintf("x[1]=%d x[2]=%d y[1]=%d y[2]=%d\n",x[1],x[2],y[1],y[2]);
#endif

    int lpp=-1; // last peak position
    double lpv=-1000; // last peak value
    
    double ppv=-1000; // last value
    int ppl=-1; // position of the last value
    double pppv=-1000; // value before last


    if(ignore_masking==1) {
      for(int i=0;i<whs;i++) {
	pn1+=x[i];
	pn2+=x[i+whs+1];
	nn1+=y[i];
	nn2+=y[i+whs+1];

      }
    }

    if(bg_subtract) {
      // pre-initialize background tag counts, 
      for(int i=0;i<bg_whs;i++) {
	if(i<n_x) {
	  bg_pn2+=bgp[i];
	  bg_nn2+=bgn[i];
	}
      }
      // increment center of background count window to the start position
      for(int i=0;i<start;i++) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

    }

    
#ifdef DEBUG  
    Rprintf("initialization: i=%d pn1=%d, pn2=%d, nn1=%d, nn2=%d", i,pn1,pn2,nn1,nn2);
#endif

    for(i=start;i<end;i++) {
      if(bg_subtract) {
	// update background counts
	int nl=i-bg_whs-1;

	if(nl>=0) {
	  bg_pn1-=bgp[nl];
	  bg_nn1-=bgn[nl];
	}
	bg_pn1+=bgp[i];
	bg_nn1+=bgn[i];

	if(i>0) {
	  bg_pn2-=bgp[i-1];
	  bg_nn2-=bgn[i-1];
	}
	int nr=i+bg_whs;
	if(nr<n_x) {
	  bg_pn2+=bgp[nr];
	  bg_nn2+=bgn[nr];
	}
      }

      // update counts
      if(ignore_masking==1) {
	pn1+=x[i-1]-x[i-whs-1];
	pn2+=x[i+whs]-x[i-1];
	nn1+=y[i-1]-y[i-whs-1];
	nn2+=y[i+whs]-y[i-1];

      } else {

	pn1=pn2=nn1=nn2=0;
	
	for(int k=0;k<whs;k++) {
	  int xp1=x[i-k-1];
	  int xp2=x[i+k];
	  int xn1=y[i-k-1];
	  int xn2=y[i+k];

	  // omit masked positions
	  if(xp1!=-1 && xn1!=-1 && xp2!=-1 && xn2!=-1) {
	    pn1+=xp1;
	    nn1+=xn1;
	    pn2+=xp2;
	    nn2+=xn2;
	  }
	}
      }

      double val;
      double spn1=((double)pn1)*tag_weight;
      double snn1=((double)nn1)*tag_weight;
      double spn2=((double)pn2)*tag_weight;
      double snn2=((double)nn2)*tag_weight;
      if(round_up) {
	if(pn1>0 && spn1<1) { spn1=1.0; }
	//if(pn2>0 && spn2<1) { spn2=1.0; }
	//if(nn1>0 && snn1<1) { snn1=1.0; }
	if(nn2>0 && snn2<1) { snn2=1.0; }
      }

      if(direct_count) {
	val=spn1+snn2;
	if(round_up && val<1) {
	  val=1.0;
	}
	if(bg_subtract) {
	  val-=((double) (bg_pn1+bg_nn2))*bg_weight;
	}
      } else {
	if(bg_subtract) {
	  spn1-=((double)bg_pn1)*bg_weight;
	  snn1-=((double)bg_nn1)*bg_weight;
	  spn2-=((double)bg_pn2)*bg_weight;
	  snn2-=((double)bg_nn2)*bg_weight;

	  if(spn2<0) spn2=0;
	  if(snn1<0) snn1=0;

	  if(spn1>0 && snn2>0) {
	    val=(2.0*sqrt(spn1*snn2)-(spn2+snn1+1.0));
	  } else {
	    val=-(spn2+snn1+1.0);
	  }
	} else {
	  val=2.0*sqrt(spn1*snn2)-(spn2+snn1+tag_weight);
	}
      }	
      //double val=sqrt(pn1*nn2);
      //if(pn2>nn1) { val-=pn2; } else { val-=pn1; }
#ifdef DEBUG  
      Rprintf("update: i=%d pn1=%d pn2=%d nn1=%d nn2=%d val=%f\n",i,pn1,pn2,nn1,nn2,val);
#endif
      
      if(return_peaks) {
	// determine if previous position was a peak
	if(ppv>min_peak_val && ppv>val && ppv>pppv) {
	  if(lpp>0 && (i-lpp+1)>min_peak_dist) {
	    // record previous peak position
	    ppos.push_back(lpp);
	    pval.push_back(lpv);
#ifdef DEBUG  
	    Rprintf("recording peak x=%d y=%f d=%d\n",lpp,lpv,(i-lpp));
#endif	    
	    if(ppl!=-1 && ppl!=i-1) {
	      lpp=(int) round((ppl+i-1)/2);
	    } else {
	      lpp=i-1;
	    }
	    lpv=ppv;
#ifdef DEBUG  
	    Rprintf("updated peak to x=%d y=%f\n",lpp,lpv);
#endif	    
	  } else {
	    if(ppv>lpv) {
	      // update last peak positions
#ifdef DEBUG  
	      Rprintf("skipping peak x=%d y=%f d=%d in favor of x=%d y=%f\n",lpp,lpv,(i-lpp),i-1,ppv);
#endif
	      if(ppl!=-1 && ppl!=i-1) {
		lpp=(int) round((ppl+i-1)/2);
	      } else {
		lpp=i-1;
	      }
	      lpv=ppv;
	    }
	  }
	}
	
	// update previous values
	if(val!=ppv) {
	  pppv=ppv; ppv=val; ppl=i;
	}
      } else {
	d_nv[i]=val;
      }
    }

    if(return_peaks) {
      // record last position
      if(lpp>0) {
#ifdef DEBUG  
	Rprintf("recording last peak x=%d y=%f\n",lpp,lpv);
#endif
	ppos.push_back(lpp);
	pval.push_back(lpv);
      }

      SEXP rpp_R,rpv_R;
      PROTECT(rpp_R=allocVector(INTSXP,ppos.size())); 
      PROTECT(rpv_R=allocVector(REALSXP,ppos.size())); 
      int* rpp=INTEGER(rpp_R);
      double* rpv=REAL(rpv_R);

      int iloop=ppos.size(); //creating dummy variable for the loop below to avoid [-Wsign-compare]
	  for(int i=0;i<iloop;i++) {
      //for(int i=0;i<ppos.size();i++) {
	rpp[i]=ppos[i];
	rpv[i]=pval[i];
      }
    
      SEXP ans_R, names_R;
      PROTECT(names_R = allocVector(STRSXP, 2));
      SET_STRING_ELT(names_R, 0, mkChar("x"));
      SET_STRING_ELT(names_R, 1, mkChar("v"));
    
      PROTECT(ans_R = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(ans_R, 0, rpp_R);
      SET_VECTOR_ELT(ans_R, 1, rpv_R);
      setAttrib(ans_R, R_NamesSymbol, names_R);
  
      UNPROTECT(4);
      return(ans_R);
    } else {
      UNPROTECT(1);
      return(nv);
    }

  }


}




// per-shift sums for tag.scc: ptc ascending, ntc descending (both strictly), values ptv/ntv.
// Returns a 4 x nshift matrix: matched sum of ptv*ntv, unmatched ptv sum, unmatched ntv sum,
// matched count. Sums accumulate in long double in R's sum() order.
extern "C" SEXP tag_scc_shifts(SEXP ptc_, SEXP ntc_, SEXP ptv_, SEXP ntv_, SEXP shifts_) {
  const int np = LENGTH(ptc_), nn = LENGTH(ntc_), ns = LENGTH(shifts_);
  const int *ptc = INTEGER(ptc_), *ntc = INTEGER(ntc_);
  const double *ptv = REAL(ptv_), *ntv = REAL(ntv_), *sh = REAL(shifts_);
  SEXP res = PROTECT(allocMatrix(REALSXP, 4, ns));
  double *out = REAL(res);
  unsigned char *mark = (unsigned char *) R_alloc(nn > 0 ? nn : 1, 1);
  for (int k = 0; k < ns; k++) {
    const long long s = (long long) sh[k];
    long double a = 0, b = 0, c = 0;
    double nm = 0;
    for (int j = 0; j < nn; j++) mark[j] = 0;
    int j = nn - 1;
    for (int i = 0; i < np; i++) {
      const long long t = (long long) ptc[i] + s;
      while (j >= 0 && (long long) ntc[j] < t) j--;
      if (j >= 0 && (long long) ntc[j] == t) {
        volatile double p = ptv[i] * ntv[j];
        a += (double) p;
        mark[j] = 1;
        nm += 1;
      } else {
        b += ptv[i];
      }
    }
    if (nm > 0) {
      for (int q = 0; q < nn; q++) if (!mark[q]) c += ntv[q];
    }
    out[4*k] = (double) a; out[4*k+1] = (double) b; out[4*k+2] = (double) c; out[4*k+3] = nm;
  }
  UNPROTECT(1);
  return res;
}
