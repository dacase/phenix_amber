SHELL=/bin/sh
VPATH=
PREC=DPREC
include ../../config.h
.SUFFIXES: .F .o

OBJ=  nxtsec.o 

#  default rules for Fortran and C compilation:

.F.o:   $<
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

libsys:  $(LIBDIR)/libsys.a

$(LIBDIR)/libsys.a:  $(OBJ)
	@echo "[LIB]  AR $@"
	$(VB)$(AR) $@ $(OBJ)  > /dev/null
	$(VB)$(RANLIB) $@

clean:
	/bin/rm -f *.o  *.a  *.mod

uninstall:
	/bin/rm -f $(LIBDIR)/libsys.a

