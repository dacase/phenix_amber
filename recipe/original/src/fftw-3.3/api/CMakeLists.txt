
add_compile_options(${STACK_ALIGN_CFLAGS})

# --------------------------------------------------------------------

#add the r2r_kind size into fftw.f03
set(FFTW_F03_HEADER ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03)
set(FFTW_F03_HEADER_TMP ${CMAKE_CURRENT_BINARY_DIR}/fftw3.f03.tmp)

math(EXPR R2R_KIND_NUM_BITS "${SIZEOF_FFTW_R2R_KIND} * 8")
set(FFTW_R2R_KIND "C_INT${R2R_KIND_NUM_BITS}_T")

# append the generated bit to the start of the file
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/fftw3.f03.in FFTW_F03_HEADER_IN_CONTENTS)

file(WRITE ${FFTW_F03_HEADER_TMP} "! Generated automatically.  DO NOT EDIT!

  integer, parameter :: C_FFTW_R2R_KIND = ${FFTW_R2R_KIND}
  ${FFTW_F03_HEADER_IN_CONTENTS}")
  
# replace the old generated header only if the old one has changed
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FFTW_F03_HEADER_TMP} ${FFTW_F03_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${FFTW_F03_HEADER_TMP})

install(FILES fftw3.h fftw3.f ${FFTW_F03_HEADER} DESTINATION ${INCDIR})

# --------------------------------------------------------------------

set(API_SOURCES apiplan.c configure.c execute-dft-c2r.c		
	execute-dft-r2c.c execute-dft.c execute-r2r.c execute-split-dft-c2r.c	
	execute-split-dft-r2c.c execute-split-dft.c execute.c			
	export-wisdom-to-file.c export-wisdom-to-string.c export-wisdom.c	
	f77api.c flops.c forget-wisdom.c import-system-wisdom.c			
	import-wisdom-from-file.c import-wisdom-from-string.c import-wisdom.c	
	malloc.c map-r2r-kind.c mapflags.c mkprinter-file.c mktensor-iodims.c	
	mktensor-rowmajor.c plan-dft-1d.c plan-dft-2d.c plan-dft-3d.c		
	plan-dft-c2r-1d.c plan-dft-c2r-2d.c plan-dft-c2r-3d.c plan-dft-c2r.c	
	plan-dft-r2c-1d.c plan-dft-r2c-2d.c plan-dft-r2c-3d.c plan-dft-r2c.c	
	plan-dft.c plan-guru-dft-c2r.c plan-guru-dft-r2c.c plan-guru-dft.c	
	plan-guru-r2r.c plan-guru-split-dft-c2r.c plan-guru-split-dft-r2c.c	
	plan-guru-split-dft.c plan-many-dft-c2r.c plan-many-dft-r2c.c		
	plan-many-dft.c plan-many-r2r.c plan-r2r-1d.c plan-r2r-2d.c		
	plan-r2r-3d.c plan-r2r.c print-plan.c rdft2-pad.c the-planner.c		
	version.c api.h f77funcs.h fftw3.h x77.h guru.h guru64.h		
	mktensor-iodims.h plan-guru-dft-c2r.h plan-guru-dft-r2c.h		
	plan-guru-dft.h plan-guru-r2r.h plan-guru-split-dft-c2r.h		
	plan-guru-split-dft-r2c.h plan-guru-split-dft.h plan-guru64-dft-c2r.c	
	plan-guru64-dft-r2c.c plan-guru64-dft.c plan-guru64-r2r.c		
	plan-guru64-split-dft-c2r.c plan-guru64-split-dft-r2c.c			
	plan-guru64-split-dft.c mktensor-iodims64.c)
	
add_library(fftw_api OBJECT ${API_SOURCES})
make_pic_if_needed(fftw_api)
