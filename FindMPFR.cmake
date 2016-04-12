# Try to find the MPFR librairies
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
		# Already in cache, be silent
		set(MPFR_FIND_QUIETLY TRUE)
endif (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)

find_path(MPFR_INCLUDE_DIR NAMES mpfr.h HINTS "./local/include")
find_library(MPFR_LIBRARIES NAMES mpfr libmpfr HINTS "./local/lib")
MESSAGE(STATUS "mpfr libs: " ${MPFR_LIBRARIES})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MPFR DEFAULT_MSG MPFR_INCLUDE_DIR MPFR_LIBRARIES)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)
