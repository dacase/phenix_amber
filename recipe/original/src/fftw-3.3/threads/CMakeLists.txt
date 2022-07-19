
# sources common to the OpenMP and non-OpenMP libraries
set(COMMON_SOURCES api.c conf.c hc2hc.c
	dft-vrank-geq1.c ct.c rdft-vrank-geq1.c 
	vrank-geq1-rdft2.c f77api.c)

set(LIBTHREADS_SOURCES threads.c ${COMMON_SOURCES})
	
set(OPENMP_SOURCES openmp.c ${COMMON_SOURCES})

# Note: these libraries depend on functions in FFTW. So, they can't be linked as shared libraries if fftw has not yet been built.
# So, we force them to be static libraries.
# Fixing this would require making two seperate CMakeLists in two seperate directories and including them before or after fftw has been linked depending on COMBINE_THREADS

if(HAVE_OPENMP)
	if(COMBINE_OPENMP_THREADS)
		#make an object library for combining into the main lib
		add_library(fftw3_omp_obj OBJECT ${OPENMP_SOURCES})
		make_pic_if_needed(fftw3_omp_obj)
		
		openmp_object_library(fftw3_omp_obj C)
	else()
		#make a standalone library
		add_library(fftw3_omp STATIC ${OPENMP_SOURCES})
		target_compile_definitions(fftw3_omp PRIVATE ${DLL_IMPORT_DEFINITIONS})
		
		install_libraries(fftw3_omp)
		
		target_link_libraries(fftw3_omp openmp_c)
	endif()
	
endif()

if(ENABLE_THREADS)
	if(COMBINE_REGULAR_THREADS)
		#make an object library for combining into the main lib
		add_library(fftw3_threads_obj OBJECT ${LIBTHREADS_SOURCES})
		make_pic_if_needed(fftw3_threads_obj)
	else()
		#make a standalone library
		add_library(fftw3_threads STATIC ${LIBTHREADS_SOURCES})
		target_compile_definitions(fftw3_threads PRIVATE ${DLL_IMPORT_DEFINITIONS})
		
		target_link_libraries(fftw3_threads Threads::Threads)
		install_libraries(fftw3_threads)
	endif()
endif()