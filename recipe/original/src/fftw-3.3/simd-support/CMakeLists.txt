set(SIMD_SUPPORT_SOURCES taint.c sse2.c avx.c altivec.c sse2-nonportable.c)

# sse2-nonportable.c needs SSE2_CFLAGS
# CMake DOES support per-source-file flags, so none of this seperate-library business from the automake code.
# This will fail if SSE2_CFLAGS ever contains more than one flag.  However, there's currently no way this can happen.
set_property(SOURCE sse2-nonportable.c PROPERTY COMPILE_FLAGS ${SSE2_CFLAGS}) 

add_library(simd_support OBJECT ${SIMD_SUPPORT_SOURCES})
make_pic_if_needed(simd_support)