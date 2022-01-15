.SUFFIXES : .F90 .o

AOBJ=   ambpdb_.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
        forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o nxtsec.o \
        amopen.o mexit.o

include ../../config.h

ambpdb: $(BINDIR)/ambpdb

$(BINDIR)/ambpdb: $(AOBJ)
	@echo "[AMBPDB]  FC $@"
	$(VB)$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(AOBJ)

clean:
	/bin/rm -f *.o *.mod

#  default rules for Fortran and C compilation:

.F90.o: $<
	@echo "[AMBPDB]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

