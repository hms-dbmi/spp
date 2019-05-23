#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    callfunction

  Most likely possible values need to be added below.
*/

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
/* .C calls */
extern void cdensum(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP ccdensum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cwindow_n_tags(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cwindow_n_tags_around(SEXP, SEXP, SEXP, SEXP);
extern SEXP find_poisson_enrichment_clusters(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_relative_coordinates(SEXP, SEXP, SEXP);
extern SEXP spp_lwcc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP points_withinC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_arachne(SEXP);
extern SEXP read_arachne_long(SEXP);
extern SEXP read_binmaqmap(SEXP, SEXP);
extern SEXP read_bowtie(SEXP, SEXP);
extern SEXP read_helicostabf(SEXP, SEXP);
extern SEXP read_maqmap(SEXP, SEXP);
extern SEXP read_meland(SEXP, SEXP);
extern SEXP read_tagalign(SEXP);
extern SEXP region_intersection(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spp_wtd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"cdensum", (DL_FUNC) &cdensum, 9},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"ccdensum",                         (DL_FUNC) &ccdensum,                          7},
    {"cwindow_n_tags",                   (DL_FUNC) &cwindow_n_tags,                    5},
    {"cwindow_n_tags_around",            (DL_FUNC) &cwindow_n_tags_around,             4},
    {"find_poisson_enrichment_clusters", (DL_FUNC) &find_poisson_enrichment_clusters,  8},
    {"get_relative_coordinates",         (DL_FUNC) &get_relative_coordinates,          3},
    {"spp_lwcc",                         (DL_FUNC) &spp_lwcc,                             14},
    {"points_withinC",                    (DL_FUNC) &points_withinC,                     6},
    {"read_arachne",                     (DL_FUNC) &read_arachne,                      1},
    {"read_arachne_long",                (DL_FUNC) &read_arachne_long,                 1},
    {"read_binmaqmap",                   (DL_FUNC) &read_binmaqmap,                    2},
    {"read_bowtie",                      (DL_FUNC) &read_bowtie,                       2},
    {"read_helicostabf",                 (DL_FUNC) &read_helicostabf,                  2},
    {"read_maqmap",                      (DL_FUNC) &read_maqmap,                       2},
    {"read_meland",                      (DL_FUNC) &read_meland,                       2},
    {"read_tagalign",                    (DL_FUNC) &read_tagalign,                     1},
    {"region_intersection",              (DL_FUNC) &region_intersection,               6},
    {"spp_wtd",                          (DL_FUNC) &spp_wtd,                              15},
    {NULL, NULL, 0}
};

void R_init_spp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}