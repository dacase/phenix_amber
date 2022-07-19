include ../../config.h

COMMON_FILES = LINE.h TOP.h TTOP.h MISC.h UNITS.h DEBUG.h T3TOP.h SIZE.h TEMP.h

LESOBJ = lesmain.o addspace.o readprm.o \
        writprm.o readcrd.o writcrd.o  \
        pick.o  	 rline.o   nxt.o \
        intgr.o   find.o    of.o \
        geti.o    unit.o    getc.o \
        alert.o   echo.o    get4c.o \
        getd.o    wlesprm.o lesprm.o random.o \
        les2prm.o checksz.o addspace_atm1st.o \
        mexit.o  nxtsec.o

addles:  $(BINDIR)/addles$(SFX)

$(BINDIR)/addles$(SFX) : $(LESOBJ) 
	@echo "[ADDLES]  FC $@"
	$(VB)$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(LESOBJ) 

#  default rules for Fortran compilation:

.SUFFIXES:
.SUFFIXES: .F .F90 .o

.F.o:   $<
	@echo "[ADDLES]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

.F90.o:   $<
	@echo "[ADDLES]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

#-----------LIBS

clean:
	-rm -f *.o

uninstall:
	-rm -f $(BINDIR)/addles$(SFX)

$(LESOBJ) : $(COMMON_FILES)

nxtsec.o:  ../msander/nxtsec.F90
	@echo "[ADDLES]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

random.o: random.F90 random.h
	@echo "[ADDLES]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ random.F90
