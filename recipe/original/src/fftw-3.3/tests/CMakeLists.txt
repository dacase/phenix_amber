# Sources used by fftw_benchmark and fftw_wisdom
set(FFTW_BENCHMARK_COMMON_SOURCES fftw-bench.c bench.c)

set(FFTW_BENCHMARK_SOURCES hook.c)

# --------------------------------------------------------------------

add_definitions(${DLL_IMPORT_DEFINITIONS})

add_library(fftw_benchmark_common_obj OBJECT ${FFTW_BENCHMARK_COMMON_SOURCES})

add_executable(fftw_benchmark $<TARGET_OBJECTS:fftw_benchmark_common_obj> ${FFTW_BENCHMARK_SOURCES})
target_link_libraries(fftw_benchmark libfftw_bench ${SUPPORT_PROGRAM_THREAD_LIB} fftw)

if(HAVE_OPENMP)
	target_link_libraries(fftw_benchmark openmp_c)
endif()