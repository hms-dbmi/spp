AC_INIT([SPP], [1.7], [Peter_Kharchenko@hms.harvard.edu])
AC_CHECK_LIB(bz2, BZ2_bzDecompressInit)
AC_SUBST(HAVE_LIBBZ2)
BOOST_REQUIRE([1.41.0])
BOOST_CPPFLAGS="${BOOST_CPPFLAGS}"
AC_SUBST(BOOST_CPPFLAGS)


# find R and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`;

CXX11=`"${R_HOME}/bin/R" CMD config CXX11`
CXX11STD=`"${R_HOME}/bin/R" CMD config CXX11STD`
CXX="${CXX11} ${CXX11STD}"
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX11FLAGS`
AC_LANG(C++)

CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"
CXXFLAGS="${CXXFLAGS}  ${PKG_CXXFLAGS}"
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CFLAGS="${CFLAGS} ${PKG_CFLAGS}"

RINC=`${R_HOME}/bin/R CMD config --cppflags`
RCPP_CXX=`${R_HOME}/bin/Rscript -e "Rcpp:::CxxFlags()"`
CPPFLAGS="${CPPFLAGS} ${RCPP_CXX}"
RCPP_LIBS=`${R_HOME}/bin/Rscript -e "Rcpp:::LdFlags()"`
LIBS="${LIBS} ${PKG_LIBS} ${RCPP_LIBS}" 

AC_SUBST(R_HOME)
AC_SUBST(RINC)
AC_SUBST(RLD)

AC_ARG_VAR([PKG_CPPFLAGS],[additional pre-processor flags])
AC_ARG_VAR([PKG_LIBS],[additional linker library flags])

AC_CONFIG_FILES([src/Makevars])
cp confdefs.h src/config.h
AC_OUTPUT
