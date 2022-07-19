include ../../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.f.o:
	@echo "[BLAS]  FC $<"
	$(VB)$(FC) -c $(FOPTFLAGS) $(FFLAGS) -o $@ $<

OBJS=    dasum.o daxpy.o dcabs1.o dcopy.o ddot.o \
         dgbmv.o dgemm.o dgemv.o dger.o dnrm2.o drot.o \
         drotg.o dsbmv.o dscal.o dspmv.o dspr.o dspr2.o \
         dswap.o dsymm.o dsymv.o dsyr.o dsyr2.o dsyr2k.o \
         dsyrk.o dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o \
         dtrmv.o dtrsm.o dtrsv.o dummy.o dzasum.o dznrm2.o \
         idamax.o isamax.o lsame.o xerbla.o 

libblas: $(LIBDIR)/libblas.a

$(LIBDIR)/libblas.a: $(OBJS)
	@echo "[BLAS]  AR $@"
	$(VB)$(AR) libblas.a $(OBJS) > /dev/null
	$(VB)$(RANLIB) libblas.a
	$(VB)mv libblas.a $(LIBDIR)

skip:
	@echo "skipping compilation of BLAS"

clean:
	-rm -f *.o 

uninstall:
	-rm -f $(LIBDIR)/libblas.a

