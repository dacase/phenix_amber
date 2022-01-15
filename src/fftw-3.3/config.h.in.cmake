/* config.h.in.  Used by CMake to generate config.h.  */

/* Defined by CMake to signal that Fortran is disabled. */
#cmakedefine DISABLE_FORTRAN 1

/* Wrapper around the CMake fortran mangling conversion functions that uses the automake name F77_FUNC */
#ifdef DISABLE_FORTRAN 
	#define F77_FUNC(namelcase, NAMEUCASE)
	#define F77_FUNC_(namelcase, NAMEUCASE) 
#else	
	#include "fortran-mangling.h"
	#define F77_FUNC(namelcase, NAMEUCASE) FORTRAN_MANGLEGLOBAL(namelcase, NAMEUCASE)
	#define F77_FUNC_(namelcase, NAMEUCASE) FORTRAN_MANGLEGLOBAL_(namelcase, NAMEUCASE)
#endif

/* Define to compile in long-double precision. */
#cmakedefine BENCHFFT_LDOUBLE 1

/* Define to compile in quad precision. */
#cmakedefine BENCHFFT_QUAD 1

/* Define to compile in single precision. */
#cmakedefine BENCHFFT_SINGLE 1

/* Define to disable Fortran wrappers. */
#cmakedefine DISABLE_FORTRAN 1

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. 
   TODO: add a check for this*/
/* #undef F77_DUMMY_MAIN */

/* Define if F77_FUNC and F77_FUNC_ are equivalent. */
#cmakedefine F77_FUNC_EQUIV 1

/* Define if F77 and FC dummy `main' functions are identical. */
#cmakedefine FC_DUMMY_MAIN_EQ_F77 1

/* C compiler name and flags */
#define FFTW_CC "${FFTW_CC}"

/* Version of library */
#define VERSION "${PACKAGE_VERSION}"
#define PACKAGE_VERSION "${PACKAGE_VERSION}"

/* Define to enable extra FFTW debugging code. */
#cmakedefine FFTW_DEBUG 1

/* Define to enable alignment debugging hacks. */
#cmakedefine FFTW_DEBUG_ALIGNMENT 1

/* Define to enable debugging malloc. */
#cmakedefine FFTW_DEBUG_MALLOC 1

/* Define to enable the use of alloca(). */
#cmakedefine FFTW_ENABLE_ALLOCA 1

/* Define to compile in long-double precision. */
#cmakedefine FFTW_LDOUBLE 1

/* Define to compile in quad precision. */
#cmakedefine FFTW_QUAD 1

/* Define to enable pseudorandom estimate planning for debugging. */
#cmakedefine FFTW_RANDOM_ESTIMATOR 1

/* Define to compile in single precision. */
#cmakedefine FFTW_SINGLE 1

/* Define to 1 if you have the `abort' function. */
#cmakedefine HAVE_ABORT 1

/* Define to 1 if you have `alloca', as a function or macro. */
#cmakedefine HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#cmakedefine HAVE_ALLOCA_H 1

/* Define to enable Altivec optimizations. */
#cmakedefine HAVE_ALTIVEC 1

/* Define to 1 if you have the <altivec.h> header file. */
#cmakedefine HAVE_ALTIVEC_H 1

/* Define to enable AVX optimizations. */
#cmakedefine HAVE_AVX 1

/* Define to 1 if you have the `BSDgettimeofday' function. */
#cmakedefine HAVE_BSDGETTIMEOFDAY 1

/* Define to 1 if you have the `clock_gettime' function. */
#cmakedefine HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the `cosl' function. */
#cmakedefine HAVE_COSL 1

/* Define to 1 if you have the <c_asm.h> header file. */
#cmakedefine HAVE_C_ASM_H 1

/* Define to 1 if you have the declaration of `cosl', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_COSL 1

/* Define to 1 if you have the declaration of `cosq', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_COSQ 1

/* Define to 1 if you have the declaration of `drand48', and to 0 if you
   don't. */
#cmakedefine HAVE_DECL_DRAND48 1

/* Define to 1 if you have the declaration of `memalign', and to 0 if you
   don't. */
#cmakedefine HAVE_DECL_MEMALIGN 1

/* Define to 1 if you have the declaration of `posix_memalign', and to 0 if
   you don't. */
#cmakedefine HAVE_DECL_POSIX_MEMALIGN 1

/* Define to 1 if you have the declaration of `sinl', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_SINL 1

/* Define to 1 if you have the declaration of `sinq', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_SINQ 1

/* Define to 1 if you have the declaration of `srand48', and to 0 if you
   don't. */
#cmakedefine HAVE_DECL_SRAND48 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
#cmakedefine HAVE_DOPRNT 1

/* Define to 1 if you have the `drand48' function. */
#cmakedefine HAVE_DRAND48 1

/* Define if you have a machine with fused multiply-add */
#cmakedefine HAVE_FMA 1

/* Define to 1 if you have the `gethrtime' function. */
#cmakedefine HAVE_GETHRTIME 1

/* Define to 1 if you have the `gettimeofday' function. */
#cmakedefine HAVE_GETTIMEOFDAY 1

/* Define to 1 if hrtime_t is defined in <sys/time.h> */
#cmakedefine HAVE_HRTIME_T 1

/* Define to 1 if you have the <intrinsics.h> header file. */
#cmakedefine HAVE_INTRINSICS_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H 1

/* Define if the isnan() function/macro is available. */
#cmakedefine HAVE_ISNAN 1

/* Define to 1 if you have the <libintl.h> header file. */
#cmakedefine HAVE_LIBINTL_H 1

/* Define to 1 if you have the `m' library (-lm). */
#cmakedefine HAVE_LIBM 1

/* Define to 1 if you have the `quadmath' library (-lquadmath). */
#cmakedefine HAVE_LIBQUADMATH 1

/* Define to 1 if you have the <limits.h> header file. */
#cmakedefine HAVE_LIMITS_H 1

/* Define to 1 if the compiler supports `long double' */
#cmakedefine HAVE_LONG_DOUBLE 1

/* Define to 1 if you have the `mach_absolute_time' function. */
#cmakedefine HAVE_MACH_ABSOLUTE_TIME 1

/* Define to 1 if you have the <mach/mach_time.h> header file. */
#cmakedefine HAVE_MACH_MACH_TIME_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H 1

/* Define to 1 if you have the `memalign' function. */
#cmakedefine HAVE_MEMALIGN 1

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#cmakedefine HAVE_MEMSET 1

/* Define to enable MIPS paired-single optimizations. */
#cmakedefine HAVE_MIPS_PS 1

/* Define to enable use of MIPS ZBus cycle-counter. */
#cmakedefine HAVE_MIPS_ZBUS_TIMER 1

/* Define if you have the MPI library. */
#cmakedefine HAVE_MPI 1

/* Define if OpenMP is enabled */
#cmakedefine HAVE_OPENMP 1

/* Define to 1 if you have the `posix_memalign' function. */
#cmakedefine HAVE_POSIX_MEMALIGN 1

/* Define if you have POSIX threads libraries and header files. */
#cmakedefine HAVE_PTHREAD 1

/* Define to 1 if you have the `read_real_time' function. */
#cmakedefine HAVE_READ_REAL_TIME 1

/* Define to 1 if you have the `sinl' function. */
#cmakedefine HAVE_SINL 1

/* Define to 1 if you have the `snprintf' function. */
#cmakedefine HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
#cmakedefine HAVE_SQRT 1

/* Define to enable SSE/SSE2 optimizations. */
#cmakedefine HAVE_SSE2 1

/* Define to 1 if you have the <stddef.h> header file. */
#cmakedefine HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H 1

/* Define to 1 if you have the `sysctl' function. */
#cmakedefine HAVE_SYSCTL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/sysctl.h> header file. */
#cmakedefine HAVE_SYS_SYSCTL_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#cmakedefine HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `tanl' function. */
#cmakedefine HAVE_TANL 1

/* Define if we have a threads library. */
#cmakedefine HAVE_THREADS 1

/* Define to 1 if you have the `time_base_to_time' function. */
#cmakedefine HAVE_TIME_BASE_TO_TIME 1

/* Define to 1 if the system has the type `uintptr_t'. */
#cmakedefine HAVE_UINTPTR_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#cmakedefine HAVE_VPRINTF 1

/* Define to 1 if you have the `_mm_free' function. */
#cmakedefine HAVE__MM_FREE 1

/* Define to 1 if you have the `_mm_malloc' function. */
#cmakedefine HAVE__MM_MALLOC 1

/* Define if you have the UNICOS _rtc() intrinsic. */
#cmakedefine HAVE__RTC 1

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
#cmakedefine PTHREAD_CREATE_JOINABLE ${PTHREAD_CREATE_JOINABLE}

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE ${SIZEOF_DOUBLE}

/* The size of `fftw_r2r_kind', as computed by sizeof. */
#define SIZEOF_FFTW_R2R_KIND ${SIZEOF_FFTW_R2R_KIND}

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT ${SIZEOF_FLOAT}

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT ${SIZEOF_INT}

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG ${SIZEOF_LONG}

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG ${SIZEOF_LONG_LONG}

/* The size of `MPI_Fint', as computed by sizeof. */
#cmakedefine SIZEOF_MPI_FINT ${SIZEOF_MPI_FINT}

/* The size of `ptrdiff_t', as computed by sizeof. */
#define SIZEOF_PTRDIFF_T ${SIZEOF_PTRDIFF_T}

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T ${SIZEOF_SIZE_T}

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT ${SIZEOF_UNSIGNED_INT}

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG ${SIZEOF_UNSIGNED_LONG}

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG ${SIZEOF_UNSIGNED_LONG_LONG}

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P}

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#cmakedefine TIME_WITH_SYS_TIME 1

/* Define if we have and are using POSIX threads. */
#cmakedefine USING_POSIX_THREADS 1

/* Use common Windows Fortran mangling styles for the Fortran interfaces. */
#cmakedefine WINDOWS_F77_MANGLING 1

/* Include g77-compatible wrappers in addition to any other Fortran wrappers.
   */
#cmakedefine WITH_G77_WRAPPERS 1

/* Use our own aligned malloc routine; mainly helpful for Windows systems
   lacking aligned allocation system-library routines. */
#cmakedefine WITH_OUR_MALLOC 1

/* Use low-precision timers, making planner very slow */
#cmakedefine WITH_SLOW_TIMER 1

/* Define to empty if `const' does not conform to ANSI C. */
#cmakedefine const

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name. 
#ifndef __cplusplus
#undef inline
#endif */