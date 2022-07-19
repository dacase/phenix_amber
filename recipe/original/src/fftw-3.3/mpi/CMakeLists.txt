# CMake buildfile for fftw3's MPI support library

# --------------------------------------------------------------------
# Source lists

set(TRANSPOSE_SRC transpose-alltoall.c transpose-pairwise.c transpose-recurse.c transpose-problem.c transpose-solve.c)
set(DFT_SRC dft-serial.c dft-rank-geq2.c dft-rank-geq2-transposed.c dft-rank1.c dft-rank1-bigvec.c dft-problem.c dft-solve.c)
set(RDFT_SRC rdft-serial.c rdft-rank-geq2.c rdft-rank-geq2-transposed.c rdft-rank1-bigvec.c rdft-problem.c rdft-solve.c)
set(RDFT2_SRC rdft2-serial.c rdft2-rank-geq2.c rdft2-rank-geq2-transposed.c rdft2-problem.c rdft2-solve.c)
set(SRC any-true.c api.c block.c choose-radix.c conf.c dtensor.c rearrange.c wisdom-api.c f03-wrap.c)

set(FFTW_MPI_SOURCES ${TRANSPOSE_SRC} ${DFT_SRC} ${RDFT_SRC} ${RDFT2_SRC} ${SRC})

set(MPI_BENCH_SOURCES mpi-bench.c ../tests/hook.c ../tests/fftw-bench.c)

# --------------------------------------------------------------------
# Configure header

math(EXPR MPI_FINT_NUM_BITS "${SIZEOF_MPI_FINT} * 8")

set(FFTW_MPI_F03_HEADER ${CMAKE_CURRENT_BINARY_DIR}/fftw3-mpi.f03)

#substitute C_MPI_FINT into the intermediate include file
configuretime_file_replace(${CMAKE_CURRENT_SOURCE_DIR}/fftw3-mpi.f03.in ${FFTW_MPI_F03_HEADER} TO_REPLACE "C_MPI_FINT" REPLACEMENT "C_INT${MPI_FINT_NUM_BITS}_T")

install(FILES ${FFTW_MPI_F03_HEADER} DESTINATION ${INCDIR})

# --------------------------------------------------------------------
# Build library and program

include_directories(${MPI_C_INCLUDE_PATH})
add_compile_options(${MPI_C_COMPILE_FLAGS})
add_definitions(${DLL_IMPORT_DEFINITIONS}) # we USE libfftw

add_library(fftw_mpi ${FFTW_MPI_SOURCES})
target_include_directories(fftw_mpi PUBLIC . ${CMAKE_CURRENT_BINARY_DIR}) 
target_link_libraries(fftw_mpi ${SUPPORT_PROGRAM_THREAD_LIB} fftw mpi_c)

add_executable(fftw_mpi_benchmark ${MPI_BENCH_SOURCES})
target_link_libraries(fftw_mpi_benchmark libfftw_bench fftw_mpi)

if(HAVE_OPENMP)
	target_link_libraries(fftw_mpi openmp_c)
endif()

install_libraries(fftw_mpi)


