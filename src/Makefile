#   top-level Makefile for msander

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd msander && $(MAKE) $(LIBSANDER)
	@cd msander && $(MAKE) $(LIBSANDERLES)
	@cd pysander && $(MAKE) $(PYSANDER)
	@cd parmed && $(MAKE) $(PARMED)
	@cd addles && $(MAKE) addles
	@cd leap && $(MAKE) tleap
	@cd antechamber && $(MAKE) ac
	@cd ambpdb && $(MAKE) ambpdb
	@cd sqm && $(MAKE) sqm
	@cp pdb4amber/pdb4amber $(BINDIR)

parallel:
	@cd msander && $(MAKE) msander.MPI
#	@cd msander && $(MAKE) msander.LES

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd addles && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd ambpdb && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )

uninstall: clean
	-(cd $(BINDIR) && /bin/rm -f msander msander.MPI \
        teLeap tleap addles pdb4amber)
	-(cd $(LIBDIR) && /bin/rm -f liblapack.a libblas.a \
        libarpack.a libsqm.a librism.a /libsys.a \
        libsander.so libsanderles.so)
	/bin/rm -f $(INCDIR)/* ../config.h
	/bin/rm -rf $(BINDIR)/x $(LIBDIR)/python2.7

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)
	-(/bin/rm -f netcdf*.log )
	-(cd netcdf-c-4.7.2 && make uninstall)
	-(cd netcdf-f-4.5.2 && make uninstall)
	-(cd fftw-3.3 && make uninstall)

