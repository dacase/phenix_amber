#CMake buildfile for libbench, a profiling library
set(LIBBENCH_SOURCES after-ccopy-from.c after-ccopy-to.c
	after-hccopy-from.c after-hccopy-to.c after-rcopy-from.c	
	after-rcopy-to.c aligned-main.c allocate.c aset.c			
	bench-cost-postprocess.c bench-exit.c bench-main.c can-do.c caset.c	
	dotens2.c info.c main.c mflops.c mp.c ovtpvt.c pow2.c problem.c		
	report.c speed.c tensor.c timer.c useropt.c util.c verify-dft.c		
	verify-lib.c verify-r2r.c verify-rdft2.c verify.c zero.c bench-user.h	
	bench.h verify.h my-getopt.c my-getopt.h)

# must be static because it contains undefined references
add_library(libfftw_bench STATIC ${LIBBENCH_SOURCES})

if(${LIBM})
	target_link_libraries(libfftw_bench ${LIBM})
endif()

make_pic_if_needed(libfftw_bench)
remove_prefix(libfftw_bench)