#include "pc.h"
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <algorithm>
#include <string>
#include <functional>
#include <utility>
#include <ext/hash_map>
#include <boost/tokenizer.hpp>

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 


class lessAbsoluteValue {
public:
  bool operator()(int a, int b) const {
    return abs(a) < abs(b);
  }
};





/**
 * Read in .bed data into a list chromosome of vectors representing 5' positions, with sign
 * corresponding to the strand.
 */

//#define DEBUG 1

extern "C" {
SEXP read_bed_ends(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");


  ifstream bed_file(fname);

#ifdef DEBUG  
  Rprintf("opened %s\n",fname);
#endif

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
    
  int fcount=0;
  while(getline(bed_file,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++; //chr=chr.substr(3,strlen(chr.c_str()));
      string str_start=*sit++;
      int fstart=atoi(str_start.c_str());
      string str_end=*sit++;
      int fend=atoi(str_end.c_str());
      int fpos=fstart;
      if(sit!=tok.end()) {
         string u0=*sit++;
         string nfield=*sit++;
         string strand=*sit++;
         if(strand=="-") { 
	   fpos=-1*fend;
         }
      }

      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d\n",chr.c_str(),cind,fpos);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  bed_file.close();
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
    sort(csi->begin(), csi->end(), lessAbsoluteValue());
  }

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    SEXP nv;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_nv=INTEGER(nv);
    int i=0;
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_nv[i++]=*pi;
    }
    SET_VECTOR_ELT(ans, csi-pos.begin(), nv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



SEXP read_meland_old(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");


  ifstream bed_file(fname);

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
    
  int fcount=0;
  while(getline(bed_file,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      sit++; sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string str_len=*sit++;
      int len=atoi(str_len.c_str());
      string chr=*sit++; chr=chr.substr(3,strlen(chr.c_str()));
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  bed_file.close();
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    
    
    
    SEXP tv,nv,lv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 3));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  int get_a_line(FILE *f,string& line) {
    line="";
    char cline[1024];
    if(fgets(cline,1024,f)) {
      line+=cline;
      return(1);
    } else {
      return(0);
    }
  }


  SEXP read_meland(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<int> > poslen; // length
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  
  Rprintf("opened %s\n",fname);


  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string str_len=*sit++;
      int len=atoi(str_len.c_str());
      string chr=*sit++; chr=chr.substr(3,strlen(chr.c_str()));
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	poslen.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      (poslen[cind]).push_back(len);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());
    lsi=poslen.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    SET_STRING_ELT(dnames_R, 2, mkChar("l"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 3, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,lv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 3+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 3, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



// reads regular eland files, recording mismatch positions
SEXP read_eland_mismatches(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > mm1; // position of the first mismatch (or 0 for none)
  vector< vector<int> > mm2; // position of the second mismatch

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }

  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      sit++; 
      string seq=*sit++; 
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string chr=*sit++; 
      // extract chromosome name from this
      int chrp=chr.find("chr");
      int pp=chr.find('.');
      chr=chr.substr(chrp+3,pp-chrp-3);
      
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());


      string strand=*sit++;
      int nstrand=0;
      if(strand=="R") { 
	fpos=-1*(fpos+seq.size()-1);
	nstrand=1;
      }

      sit++;
      
      int nm1=0; int nm2=0;
      if(sit!=tok.end()) {
	string nms=*sit++;
	nm1=atoi(nms.substr(0,nms.size()-1).c_str());
	if(nstrand) { nm1=seq.size()-nm1+1; }
      }
      if(sit!=tok.end()) {
	string nms=*sit++;
	nm2=atoi(nms.substr(0,nms.size()-1).c_str());
	if(nstrand) { nm2=seq.size()-nm2+1; }
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	mm1.push_back(vector<int>());
	mm2.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (mm1[cind]).push_back(nm1);
      (mm2[cind]).push_back(nm2);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm1=%d, nm2=%d\n",chr.c_str(),cind,fpos,nm1,nm2);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
    
    
#ifdef DEBUG  
  Rprintf("done. read %d fragments\n",fcount);
#endif

  Rprintf("done. read %d fragments\n",fcount);

    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi,lsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=mm1.begin()+(csi-pos.begin());
    lsi=mm2.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 3)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("f"));
    SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    
    
    
    SEXP tv,nv,lv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(lv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    int* i_lv=INTEGER(lv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    vector<int>::const_iterator ili=lsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i_lv[i]=*ili++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 3));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    SET_VECTOR_ELT(dv, 2, lv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in regular eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_eland(SEXP filename,SEXP read_tag_names_R,SEXP eland_tag_length_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
  int eland_tag_length=*(INTEGER(eland_tag_length_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string sequence=*sit++;
      int len=sequence.size();
      // adjust probe length if eland length limit was specified
      if(eland_tag_length>0 && len>eland_tag_length) {
	len=eland_tag_length;
      }
      string str_nm=*sit++;
      int nm=0;
      if(str_nm[0]=='U') {
	nm=atoi((str_nm.c_str()+1));
      } else {
	continue;
      }
      sit++; sit++; sit++;
      string chr=*sit++; 
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;

      if(str_strand[0]=='R') {
	fpos=-1*(fpos+len-1);
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



  // read in extended eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_eland_extended(SEXP filename,SEXP read_tag_names_R,SEXP eland_tag_length_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
  int eland_tag_length=*(INTEGER(eland_tag_length_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string machinename=*sit++;
      string runnumber=*sit++;
      string lanenumber=*sit++;
      *sit++;
      
      string str_x=*sit++;
      string str_y=*sit++;

      string tagname=machinename+"."+runnumber+"."+lanenumber+"."+str_x+"."+str_y;

      

      *sit++;
      *sit++;

      
      string sequence=*sit++;
      *sit++;
      
      string chr=*sit++; 
      string contig=*sit++; 
      chr=chr+contig;
      
      int len=sequence.size();
      // adjust probe length if eland length limit was specified
      if(eland_tag_length>0 && len>eland_tag_length) {
	len=eland_tag_length;
      }


      
      string str_pos=*sit++;
      if(str_pos.size()<1) { continue; }
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;

      if(str_strand[0]=='R') {
	fpos=-1*(fpos+len-1);
      }

      string str_nm=*sit++;
      // count non-digit characters
      int nm=0;
      for(int i=0;i<str_nm.size();i++) {
	if(!isdigit(str_nm[i])) { nm++; }
      }
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}


  // read in regular eland files, adjusting the negative strand coordinate by sequence length
  SEXP read_bowtie(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);
  boost::char_separator<char> sep2(",");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string str_strand=*sit++;
      string chr=*sit++; 

      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());

      string sequence=*sit++;
      sit++; sit++;
      string mm=*sit++;

      int len=sequence.size();
      if(str_strand[0]=='-') {
	fpos=-1*(fpos+len-1);
      }
      // determine number of mismatches
      int nm=0;
      if(mm.size()>0) {
	nm++;
	string::size_type tp(0);
	while(tp!=string::npos) {
	  tp = mm.find(",",tp);
	  if(tp!=string::npos) {
	    tp++;
	    ++nm;
	  }
	}
      }


      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



  // read in text version of maq map
  SEXP read_maqmap(SEXP filename,SEXP read_tag_names_R) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
  int read_names=*(INTEGER(read_tag_names_R));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches
  vector< vector<string> > tagnames;

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep("\t","",boost::keep_empty_tokens);

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string tagname=*sit++;
      string chr=*sit++;
      string str_pos=*sit++;
      int fpos=atoi(str_pos.c_str());
      string str_strand=*sit++;
      sit++; sit++; sit++; sit++; sit++; 
      string str_nm=*sit++;
      sit++; sit++; sit++; 
      string str_len=*sit++;
      int nm=atoi(str_nm.c_str());
      int len=atoi(str_len.c_str());

      if(str_strand[0]=='-') {
	fpos=-1*(fpos+len-1);
      }

      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
	if(read_names) {
	  tagnames.push_back(vector<string>());
	}
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
      if(read_names) {
	(tagnames[cind]).push_back(tagname);
      }
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d, nm=%d, len=%d\n",chr.c_str(),cind,fpos,nm,len);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  vector<vector<string> >::const_iterator ssi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2+read_names)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    if(read_names) {
      SET_STRING_ELT(dnames_R, 2, mkChar("s"));
    }
    
    
    
    SEXP tv,nv,sv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    if(read_names) {
      PROTECT(sv=allocVector(STRSXP,csi->size()));   np++;
    }
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    if(read_names) {
      int i=0;
      ssi=tagnames.begin()+(csi-pos.begin());
      for(vector<string>::const_iterator si=ssi->begin();si!=ssi->end();++si) {
	SET_STRING_ELT(sv,i,mkChar(si->c_str()));
	i++;
      }
    }
    PROTECT(dv = allocVector(VECSXP, 2+read_names));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    if(read_names) {
      SET_VECTOR_ELT(dv, 2, sv);
    }
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



  // read in tagalign file
  SEXP read_tagalign(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++;
      string str_spos=*sit++;
      string str_epos=*sit++;
      sit++; 
      string str_qual=*sit++;
      string str_strand=*sit;

      int fpos;
      if(str_strand[0]=='+') {
	fpos=atoi(str_spos.c_str());
      } else {
	fpos=-1*atoi(str_epos.c_str());
      }
      int nm=atoi(str_qual.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d nm=%d\n",chr.c_str(),cind,fpos,nm);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    
    
    SEXP tv,nv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 2));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}




  // arachne madness
  SEXP read_arachne(SEXP filename) {

#ifdef DEBUG  
  Rprintf("start\n");
#endif
  const char* fname=CHAR(asChar(filename));
#ifdef DEBUG  
  Rprintf("fname=%s\n",fname);
#endif

  // main data vector
  // chr - pos
  vector< vector<int> > pos;
  vector< vector<int> > posnm; // number of mismatches

  // chromosome map
  hash_map<string, int, hash<string>,equal_to<string> > cind_map;
  vector<string> cnames;
  

  typedef boost::tokenizer<boost::char_separator<char> >  tokType;
  boost::char_separator<char> sep(" \t");

  
  FILE *f=fopen(fname,"rb");
  if (!f)  { cout<<"can't open input file \""<<fname<<"\"\n"; }
  else {
  Rprintf("opened %s\n",fname);

  // read in bed line
  string line;
  int fcount=0;
  while(get_a_line(f,line)) {

#ifdef DEBUG  
    Rprintf("line: %s\n",line.c_str());
#endif


    tokType tok(line, sep);
    tokType::iterator sit=tok.begin();
    if(sit!=tok.end()) {
      string chr=*sit++;
      string str_spos=*sit++;
      string str_mm=*sit;
      
      int fpos=atoi(str_spos.c_str());;
      int nm=atoi(str_mm.c_str());
      
      // determine the chromosome index
      hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
      int cind=-1;
      if(li==cind_map.end()) {
	// register new chromosome
	cind=cnames.size();
	cnames.push_back(chr);
	cind_map[chr]=cind;
	// allocate new pos vector
	pos.push_back(vector<int>());
	posnm.push_back(vector<int>());
#ifdef DEBUG  
	Rprintf("registered new chromosome %s with cind=%d, pos.size=%d\n",chr.c_str(),cind,pos.size());
#endif
      } else {
	cind=li->second;
      }
      fcount++;
      (pos[cind]).push_back(fpos);
      (posnm[cind]).push_back(nm);
#ifdef DEBUG  
      Rprintf("read in position chr=%s cind=%d fpos=%d nm=%d\n",chr.c_str(),cind,fpos,nm);
      if(fcount>30) {
	break;
      }
#endif
      
    }
  }
  fclose(f);
     
  Rprintf("done. read %d fragments\n",fcount);
  }
    // construct output structures
  SEXP chnames;
  int np=0; // number of protections
  PROTECT(chnames = allocVector(STRSXP, cnames.size()));
  for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
    SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
  }
  np++;

  // sort
  //for(vector<vector<int> >::iterator csi=pos.begin();csi!=pos.end();++csi) {
  //  sort(csi->begin(), csi->end(), lessAbsoluteValue());
  //}

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
  vector<vector<int> >::const_iterator nsi;
  for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
    nsi=posnm.begin()+(csi-pos.begin());

    SEXP dv,dnames_R;
    PROTECT(dnames_R = allocVector(STRSXP, 2)); np++;
    SET_STRING_ELT(dnames_R, 0, mkChar("t"));
    SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    
    
    SEXP tv,nv;
    PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
    PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
    int* i_tv=INTEGER(tv);
    int* i_nv=INTEGER(nv);
    
    int i=0;
    vector<int>::const_iterator ini=nsi->begin();
    for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
      i_tv[i]=*pi;
      i_nv[i]=*ini++;
      i++;
    }
    PROTECT(dv = allocVector(VECSXP, 2));   np++;
    SET_VECTOR_ELT(dv, 0, tv);
    SET_VECTOR_ELT(dv, 1, nv);
    setAttrib(dv, R_NamesSymbol, dnames_R);
    
    SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
  }

  setAttrib(ans,R_NamesSymbol,chnames);

#ifdef DEBUG  
  Rprintf("unprotecting %d elements\n",np);
#endif
  
  UNPROTECT(np);
  return(ans);
}



}
