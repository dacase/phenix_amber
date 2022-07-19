project(addles Fortran)

set(ADDLES_SOURCES lesmain.F addspace.F readprm.F
		writprm.F readcrd.F writcrd.F
		pick.F  	 rline.F   nxt.F
		intgr.F   find.F    of.F
		geti.F    unit.F    getc.F
		alert.F   echo.F    get4c.F
		getd.F    wlesprm.F lesprm.F random.F90
		les2prm.F checksz.F addspace_atm1st.F
		LINE.h TOP.h TTOP.h 
		MISC.h UNITS.h DEBUG.h 
		T3TOP.h SIZE.h TEMP.h)

apply_optimization_declarations(${ADDLES_SOURCES})
add_executable(addles ${ADDLES_SOURCES})
		
target_link_libraries(addles amber_common)
		
install(TARGETS addles DESTINATION ${BINDIR})
