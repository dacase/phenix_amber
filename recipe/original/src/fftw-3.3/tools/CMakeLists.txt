add_definitions(${DLL_IMPORT_DEFINITIONS})

add_executable(fftw_wisdom $<TARGET_OBJECTS:fftw_benchmark_common_obj> fftw-wisdom.c)
target_link_libraries(fftw_wisdom libfftw_bench ${SUPPORT_PROGRAM_THREAD_LIB} fftw)

if(HAVE_OPENMP)
	target_link_libraries(fftw_wisdom openmp_c)
endif()

install(TARGETS fftw_wisdom DESTINATION ${BINDIR})

# --------------------------------------------------------------------

# generate fftw-wisdom-to-conf
set(CONFIGURED_WTC_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/fftw-wisdom-to-conf)

set(PACKAGE fftw)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fftw-wisdom-to-conf.in ${CONFIGURED_WTC_SCRIPT})

install(PROGRAMS ${CONFIGURED_WTC_SCRIPT} DESTINATION ${BINDIR})