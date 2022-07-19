#CMake buildfile for pysander
#this one is tricky because it builds a native library, and it uses the same source for two libraries
project(pysander)

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/pysander-build)

set(STAMP_FILE ${BUILD_DIR}/pysander-build.stamp)

set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${BUILD_DIR})

#find all python scripts and C files as dependencies
file(GLOB_RECURSE PYSANDER_SOURCES "*.py" "*.c")

set(RPATH_ARG "")

if(TARGET_OSX)
	# OS X: set RPATH to path up from AMBERHOME/lib/python2.7/site-packages/sander/pysander.so to AMBERHOME/lib
	set(RPATH_ARG LDFLAGS="-Wl,-rpath,@loader_path/../../..")
endif()

#The install script wants to build two python modules: one with LES defined, and one without.
#So, because apparantly there's no easy way to do this, it copies the sander dir to sanderles and builds it as a separate module.
# to convince distutils that there are actually two extensions.


#NOTE: the first -I argument to setup.py must be AMBERTOOLS_INC_DIR

add_custom_command(OUTPUT ${STAMP_FILE}
	COMMAND ${CMAKE_COMMAND} -E env --unset=AMBERHOME ${PYTHON_CXX_ENVVAR_ARG} ${RPATH_ARG}
	${PYTHON_EXECUTABLE} setup.py 
	"-I${AMBERTOOLS_INC_DIR}" 
	"-I${CMAKE_SOURCE_DIR}/AmberTools/src/sander"  
	"-L$<TARGET_FILE_DIR:libsander>"  "-L$<TARGET_FILE_DIR:libsanderles>" build -b ${BUILD_DIR} ${WIN64_DEFINE_ARG} ${PYTHON_COMPILER_ARG}
	COMMAND ${CMAKE_COMMAND} -E touch ${STAMP_FILE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${PYSANDER_SOURCES} 
	VERBATIM
	COMMENT "Building pysander native library")

#We want to build the python library during the build step so as to catch any build errors, and because custom commands support generator expressions
add_custom_target(pysander ALL DEPENDS ${STAMP_FILE})
	
install(CODE "
	${FIX_BACKSLASHES_CMD}
	execute_process(
    	COMMAND \"${CMAKE_COMMAND}\" -E env --unset=AMBERHOME \"${PYTHON_EXECUTABLE}\" 
    	\"${CMAKE_CURRENT_SOURCE_DIR}/setup.py\" build -b ${BUILD_DIR}
    	 install -f ${PYTHON_PREFIX_ARG} 
    	\"--install-scripts=\${CMAKE_INSTALL_PREFIX_BS}bin\"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})" COMPONENT Python)