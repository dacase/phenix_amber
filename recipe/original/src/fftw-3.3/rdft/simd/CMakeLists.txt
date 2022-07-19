set(HC2CFDFTV hc2cfdftv_2.c hc2cfdftv_4.c hc2cfdftv_6.c hc2cfdftv_8.c	
	hc2cfdftv_10.c hc2cfdftv_12.c hc2cfdftv_16.c hc2cfdftv_32.c		
	hc2cfdftv_20.c)

set(HC2CBDFTV hc2cbdftv_2.c hc2cbdftv_4.c hc2cbdftv_6.c hc2cbdftv_8.c	
	hc2cbdftv_10.c hc2cbdftv_12.c hc2cbdftv_16.c hc2cbdftv_32.c		
	hc2cbdftv_20.c)

###########################################################################
set(SIMD_CODELETS ${HC2CFDFTV} ${HC2CBDFTV} codlist.c genus.c)

#build the codelet subdirs
if(HAVE_ALTIVEC)
	generate_and_build_fft_subdir(rdft_altivec_codelets common altivec simd-altivec.h ${SIMD_CODELETS})
	set_property(TARGET rdft_altivec_codelets PROPERTY COMPILE_OPTIONS ${ALTIVEC_CFLAGS})
endif()

if(HAVE_AVX)
	generate_and_build_fft_subdir(rdft_avx_codelets common avx simd-avx.h ${SIMD_CODELETS})
	set_property(TARGET rdft_avx_codelets PROPERTY COMPILE_OPTIONS ${AVX_CFLAGS})
endif()

if(HAVE_SSE2)
	generate_and_build_fft_subdir(rdft_sse2_codelets common sse2 simd-sse2.h ${SIMD_CODELETS})
	set_property(TARGET rdft_sse2_codelets PROPERTY COMPILE_OPTIONS ${SSE2_CFLAGS})
endif()