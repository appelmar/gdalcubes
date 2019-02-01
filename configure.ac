# Process this file with autoconf to produce a configure script.
# This file has been largely adapted from the R packages Cairo (c) Simon Urbanek and sf (c) Edzer Pebesma (see https://github.com/cran/Cairo/blob/master/configure.ac and https://github.com/r-spatial/sf/blob/master/configure.ac respectively).


# This script simply uses pkg-config if available, or tries to guess correct include / library paths.

AC_INIT(gdalcubes, 0.0.1, marius.appel@uni-muenster.de)
AC_LANG(C++)

# find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`;
CPP=`"${RBIN}" CMD config CPP`;
AC_SUBST(R_HOME)


# TODO: make sure that C++11 support is available


AC_PATH_PROG([PKGCONF],[pkg-config],[],[$PATH:/usr/local/bin:ext/bin:ext:/sw/bin:/opt/bin])
AC_PATH_PROG([FIND],[find],[],[$PATH:/usr/local/bin])

AC_MSG_CHECKING(pkg-config exists)
if test -r "${PKGCONF}"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_MSG_WARN([pkg-config not found - guessing library locations])
fi

# GDAL
if test -r "${PKGCONF}"; then
  AC_MSG_CHECKING([pkg-config knows gdal])
  if "${PKGCONF}" --exists gdal; then
    AC_MSG_RESULT([yes])
    GDAL_CFLAGS=`"${PKGCONF}" --cflags gdal`
    GDAL_LIBS=`"${PKGCONF}" --libs gdal`
  else
	  AC_MSG_RESULT([no])
  fi
fi

# if pkg-config is not available or pkg-config could not find GDAL
if test -z "${GDAL_LIBS}"; then
  AC_MSG_NOTICE([GDAL_LIBS is unset, attempting to guess it.])
    cpre=''
    for pre in /usr /usr/local /opt; do
      if test -e "${pre}/include/gdal/gdal.h"; then
        cpre=${pre}; break
      fi
    done
    if test -n "${cpre}"; then
      GDAL_CFLAGS="-I${cpre}/include/gdal"
      if test "${cpre}" = /usr; then
         GDAL_LIBS="-lgdal"
      else
         GDAL_LIBS="-L${cpre}/lib -lgdal"
      fi
   fi
fi

AC_MSG_NOTICE([GDAL_CFLAGS=${GDAL_CFLAGS}])
AC_MSG_NOTICE([GDAL_LIBS=${GDAL_LIBS}])







# proj
if test -r "${PKGCONF}"; then
  AC_MSG_CHECKING([pkg-config knows proj])
  if "${PKGCONF}" --exists proj; then
    AC_MSG_RESULT([yes])
    PROJ_CFLAGS=`"${PKGCONF}" --cflags proj`
    PROJ_LIBS=`"${PKGCONF}" --libs proj`
  else
	  AC_MSG_RESULT([no])
  fi
fi

# if pkg-config is not available or pkg-config could not find proj
if test -z "${PROJ_LIBS}"; then
  AC_MSG_NOTICE([PROJ_LIBS is unset, attempting to guess it.])
    cpre=''
    for pre in /usr /usr/local /opt; do
      if test -e "${pre}/include/proj_api.h"; then
        cpre=${pre}; break
      fi
    done
    if test -n "${cpre}"; then
      PROJ_CFLAGS="-I${cpre}/include"
      if test "${cpre}" = /usr; then
         PROJ_LIBS="-lproj"
      else
         PROJ_LIBS="-L${cpre}/lib -lproj"
      fi
   fi
fi
AC_MSG_NOTICE([PROJ_CFLAGS=${PROJ_CFLAGS}])
AC_MSG_NOTICE([PROJ_LIBS=${PROJ_LIBS}])





# NetCDF
if test -r "${PKGCONF}"; then
  AC_MSG_CHECKING([pkg-config knows netcdf])
  if "${PKGCONF}" --exists netcdf; then
    AC_MSG_RESULT([yes])
    NC_CFLAGS=`"${PKGCONF}" --cflags netcdf`
    NC_LIBS=`"${PKGCONF}" --libs netcdf`
  else
	  AC_MSG_RESULT([no])
  fi
fi

# if pkg-config is not available or pkg-config could not find netcdf
if test -z "${NC_LIBS}"; then
  AC_MSG_NOTICE([NC_LIBS is unset, attempting to guess it.])
    cpre=''
    for pre in /usr /usr/local /opt; do
      if test -e "${pre}/include/netcdf.h"; then
        cpre=${pre}; break
      fi
    done
    if test -n "${cpre}"; then
      NC_CFLAGS="-I${cpre}/include"
      if test "${cpre}" = /usr; then
         NC_LIBS="-lnetcdf"
      else
         NC_LIBS="-L${cpre}/lib -lnetcdf"
      fi
   fi
fi

AC_MSG_NOTICE([NC_CFLAGS=${NC_CFLAGS}])
AC_MSG_NOTICE([NC_LIBS=${NC_LIBS}])





# curl
if test -r "${PKGCONF}"; then
  AC_MSG_CHECKING([pkg-config knows libcurl])
  if "${PKGCONF}" --exists libcurl; then
    AC_MSG_RESULT([yes])
    CURL_CFLAGS=`"${PKGCONF}" --cflags libcurl`
    CURL_INCDIR=`"${PKGCONF}" --variable=includedir libcurl`
    CURL_CFLAGS="${CURL_CFLAGS} -I${CURL_INCDIR}/curl"
    CURL_LIBS=`"${PKGCONF}" --libs libcurl`
  else
	  AC_MSG_RESULT([no])
  fi
fi

# if pkg-config is not available or pkg-config could not find libcurl
if test -z "${CURL_LIBS}"; then
  AC_MSG_NOTICE([CURL_LIBS is unset, attempting to guess it.])
    cpre=''
    for pre in /usr /usr/local /opt; do
      if test -e "${pre}/include/curl/curl.h"; then
        cpre=${pre}; break
      fi
    done
    if test -n "${cpre}"; then
      CURL_CFLAGS="-I${cpre}/include/curl"
      if test "${cpre}" = /usr; then
         CURL_LIBS="-lcurl"
      else
         CURL_LIBS="-L${cpre}/lib -lcurl"
      fi
   fi
fi

AC_MSG_NOTICE([CURL_CFLAGS=${CURL_CFLAGS}])
AC_MSG_NOTICE([CURL_LIBS=${CURL_LIBS}])



# SQLite seems to be automatically available from R



# Bringing everything together

INLIBS="${LIBS}"
INCPPFLAGS="${CPPFLAGS}"
INPKG_CPPFLAGS="${PKG_CPPFLAGS}"
INPKG_LIBS="${PKG_LIBS}"

CPPFLAGS="$CPPFLAGS ${PKG_CPPFLAGS} ${GDAL_CFLAGS} ${NC_CFLAGS} ${PROJ_CFLAGS} ${CURL_CFLAGS} -DR_PACKAGE"
LIBS="$LIBS ${GDAL_LIBS} ${NC_LIBS} ${PROJ_LIBS} ${CURL_LIBS}"

AC_MSG_NOTICE([Package CPP flags: ${CPPFLAGS}])
AC_MSG_NOTICE([Package LIBS: ${LIBS}])


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT





