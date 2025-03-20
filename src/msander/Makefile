
include ../../config.h

SHELL=/bin/sh
VPATH=

# Add search paths for .mod files
LOCALFLAGS = $(CUDAFLAGS) -I../sqm

LESAPIOBJ= memory_module.o stack.o file_io_dat.o les.o \
    constants.o state.o wallclock.o veclib.o \
    rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o relax_mat.o nmr.LES.o multitmd.o \
    cshf.o nmrcal.o mtmdcall.o pearsn.o cpparm.LES.o \
    printe.o nonbond_list.LESAPI.o force.LES.o rdparm.LESAPI.o \
    mdread.LES.o locmem.LES.o runmd.LES.o quench.o degcnt.LES.o decnvh.o \
    fastwt.API.o parallel.LES.o shake.o rattle.o ene.LES.o \
    dynlib.LESAPI.o nmlsrc.o \
    ew_force.LES.o ew_setup.LESAPI.o ew_bspline.o\
    ew_fft.o short_ene.LES.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    amopen.o sander_lib.o random.o debug.API.o \
    mexit.o new_time.API.o extra_pts.LESAPI.o thermo_int.o matinv.o assert.o \
    erfcfun.o md_scheme.LES.o \
    ips.o sglds.o iwrap2.o emap.o \
    parms.o softcore.o mbar.o \
    linear_response.o molecule.o dssp.o \
    charmm.API.o crg_reloc.o barostats.o prmtop_type.o \
    bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o \
    mdwrit.LES.o lmod.LES.o trajene.o ew_box.LESAPI.o minrit.LES.o \
    set.LESAPI.o qm_mm.LESAPI.o amd.o scaledMD.o interface.LESAPI.o \
    music_module.o dsarpack.o nxtsec.o

NFE_OBJECTS = \
    nfe-abmd-hooks.o nfe-bbmd-ctxt.o nfe-bbmd-hooks.o  \
    nfe-colvar.o nfe-colvar-math.o nfe-colvar-type.o nfe-colvar-utils.o \
    nfe-constants.o nfe-cv-ANGLE.o nfe-cv-COM_ANGLE.o \
    nfe-cv-COM_DISTANCE.o nfe-cv-COM_TORSION.o nfe-cv-COS_OF_DIHEDRAL.o \
    nfe-cv-DISTANCE.o nfe-cv-HANDEDNESS.o nfe-cv-LCOD.o \
    nfe-cv-MULTI_RMSD.o nfe-cv-N_OF_BONDS.o nfe-cv-N_OF_STRUCTURES.o \
    nfe-cv-R_OF_GYRATION.o nfe-cv-TORSION.o  \
    nfe-pmd-hooks.o nfe-rmsd.o nfe-sander-hooks.o nfe-sander-proxy.o \
    nfe-stsm-hooks.o nfe-cv-PCA.o nfe-read-pca.o nfe-cv-SIN_OF_DIHEDRAL.o \
    nfe-cv-PAIR_DIHEDRAL.o nfe-cv-PATTERN_DIHEDRAL.o nfe-cv-DF_COM_DISTANCE.o \
    nfe-smd-hooks.o nfe-umbrella.o nfe-utils.o nfe-value.o mt19937.o \
    nfe-cv-ORIENTATION_ANGLE.o nfe-cv-ORIENTATION_PROJ.o \
	nfe-cv-QUATERNION0.o nfe-cv-QUATERNION1.o nfe-cv-QUATERNION2.o \
	nfe-cv-QUATERNION3.o nfe-cv-SPINANGLE.o nfe-cv-TILT.o

QMOBJ= \
    qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o \
    qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_genmpi_module.o tcpb_module.o \
    qm2_extern_quick_module.o  qm2_extern_reaxff_puremd_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o

APIOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    rmsgrd.o lmod.o rand2.o \
    lmodC.o xminC.o decomp.o icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o rdparm.API.o cpparm.o sander_lib.o \
    mdread.o locmem.o runmd.o quench.o getcor.o degcnt.o decnvh.o \
    fastwt.API.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
    set.API.o dynlib.API.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.API.o ew_box.API.o ew_bspline.o\
    ew_fft.o nonbond_list.API.o short_ene.o ew_recip.o \
    pcshift.o align.o csa.o rfree.o rgroup.o amopen.o \
    debug.API.o mexit.o new_time.API.o extra_pts.API.o thermo_int.o les.o \
    matinv.o assert.o erfcfun.o md_scheme.o \
    ips.o sglds.o iwrap2.o emap.o\
    AmberNetcdf.o bintraj.o binrestart.o \
    parms.o softcore.o mbar.o linear_response.o \
    molecule.o dssp.o \
    charmm.API.o crg_reloc.o amd.o scaledMD.o barostats.o interface.o \
    random.o qm_mm.API.o \
    prmtop_type.o music_module.o wallclock.o veclib.o dsarpack.o nxtsec.o

QMAPIOBJ = \
    qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o qm2_extern_nw_module.o \
    qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_genmpi_module.o tcpb_module.o \
    qm2_extern_quick_module.o  qm2_extern_reaxff_puremd_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o

MMOBJ = \
    align.o AmberNetcdf.o amd.o amopen.o assert.o \
    barostats.o binrestart.o bintraj.o charmm.o \
    constants.o crg_reloc.o csa.o cshf.o \
    debug.o decnvh.o degcnt.o dssp.o \
    dynlib.o egb.o emap.o ene.o erfcfun.o \
    ew_box.o ew_bspline.o ew_fft.o \
    ew_force.o ew_recip.o ew_setup.o extra_pts.o \
    fastwt.o file_io_dat.o findmask.o force.o getcor.o \
    icosasurf.o ips.o iwrap2.o les.o linear_response.o \
    lmodC.o lmod.o locmem.o matinv.o mbar.o \
    mdfil.o mdread.o md_scheme.o mdwrit.o memory_module.o \
    mexit.o minrit.o molecule.o mtmdcall.o multisander.o \
    multitmd.o music_module.o new_time.o nmlsrc.o nmrcal.o \
    nmr.o nonbond_list.o parallel.o parms.o pcshift.o \
    pearsn.o printe.o qm_mm.o quench.o rand2.o \
    random.o rattle.o rdparm.o relax_mat.o \
    remd.o rfree.o rgroup.o rmsgrd.o runmd.o \
    sander_lib.o sander.o scaledMD.o set.o \
    sglds.o shake.o $(SHORT_ENE) softcore.o \
    stack.o state.o thermo_int.o trajene.o \
    wallclock.o xminC.o veclib.o dsarpack.o nxtsec.o

LESOBJ= \
    align.o AmberNetcdf.o amd.o amopen.o \
    assert.o barostats.o binrestart.LES.o bintraj.LES.o charmm.o \
    constants.o crg_reloc.o csa.o \
    cshf.o debug.o decnvh.o degcnt.LES.o \
    dssp.o dynlib.LES.o egb.LES.o emap.o ene.LES.o \
    erfcfun.o ew_box.LES.o ew_bspline.o \
    ew_fft.o ew_force.LES.o ew_recip.LES.o ew_setup.LES.o \
    extra_pts.LES.o fastwt.o file_io_dat.o findmask.o force.LES.o \
    getcor.LES.o icosasurf.o ips.o iwrap2.o \
    les.o linear_response.o lmodC.o lmod.LES.o locmem.LES.o \
    matinv.o mbar.o mdfil.LES.o mdread.LES.o md_scheme.LES.o \
    mdwrit.LES.o memory_module.o mexit.o minrit.LES.o molecule.o \
    mtmdcall.o multisander.LES.o multitmd.o music_module.o new_time.o \
    nmlsrc.o nmrcal.o nmr.LES.o nonbond_list.LES.o parallel.LES.o \
    parms.o pcshift.o pearsn.o printe.o \
    qm_mm.LES.o quench.o rand2.o random.o \
    rattle.o rdparm.LES.o relax_mat.o remd.LES.o \
    rfree.o rgroup.o rmsgrd.o runmd.LES.o \
    sander.LES.o  sander_lib.o scaledMD.o set.LES.o sglds.o \
    shake.o short_ene.LES.o softcore.o \
    stack.o state.o thermo_int.o trajene.o wallclock.o \
    xminC.o veclib.o dsarpack.o nxtsec.o

.PHONY: dummy

msander$(SFX): $(BINDIR)/msander$(SFX)

msander.LES$(SFX): $(BINDIR)/msander.LES$(SFX)

msander.cuda$(SFX): $(BINDIR)/msander.cuda$(SFX)

libsander: $(LIBDIR)/libsander$(SHARED_SUFFIX)

libsanderles: $(LIBDIR)/libsanderles$(SHARED_SUFFIX)

libsander_skip::
	@echo "Skipping libsander"

libsanderles_skip::
	@echo "Skipping libsander"

dummy:

#---------------------------------------------------------------------------
$(BINDIR)/msander$(SFX): libs $(MMOBJ) $(QMOBJ) \
         $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) \
         dummy
	@echo "[MSANDER]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@  $(MMOBJ) $(QMOBJ) \
	     $(NFE_OBJECTS) \
	     -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/msander.cuda$(SFX): libs lib_cuda $(MMOBJ) $(QMOBJ) \
         $(NFE_OBJECTS) \
         $(PLUMED_DEPENDENCIES) dummy
	@echo "[MSANDER]  FC $@"
	$(VB)$(FC) -o $@  $(MMOBJ) $(QMOBJ) \
	     $(NFE_OBJECTS)  \
	     -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD) \
         -L$(CUDA_HOME)/lib64 -L$(CUDA_HOME)/lib \
         -lcurand -lcufft -lcudart -lcudadevrt -lgfortran -lstdc++

#---------------------------------------------------------------------------
$(BINDIR)/msander.LES$(SFX): libs $(LESOBJ) $(QMOBJ) \
         $(NFE_OBJECTS)  \
         $(PLUMED_DEPENDENCIES) dummy
	@echo "[MSANDER_LES]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LESOBJ) $(QMOBJ) \
	     $(NFE_OBJECTS) -L$(LIBDIR)  \
	     $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(LIBDIR)/libsander$(SHARED_SUFFIX): libs $(APIOBJ) \
          $(QMAPIOBJ) $(INCDIR)/sander_api.mod \
          force.o $(PLUMED_DEPENDENCIES) $(NFE_OBJECTS) 
	@echo "[MSANDER API]  FC $@"
	$(VB)$(FC) $(MAKE_SHARED) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  \
        $(APIOBJ) $(QMAPIOBJ) force.o \
        -L$(LIBDIR) $(NFE_OBJECTS) \
        $(LDFLAGS) $(FLIBSF) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(LIBDIR)/libsanderles$(SHARED_SUFFIX): libs $(LESAPIOBJ) \
        $(QMAPIOBJ) $(INCDIR)/sanderles_api.mod \
        $(NFE_OBJECTS) \
        $(PLUMED_DEPENDENCIES) 
	@echo "[MSANDER_LES API]  FC $@"
	$(VB)$(FC) $(MAKE_SHARED) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ \
        $(LESAPIOBJ) $(QMAPIOBJ) \
        -L$(LIBDIR) $(PLUMED_LOAD) \
        $(NFE_OBJECTS) $(LDFLAGS) $(FLIBSF) 

$(INCDIR)/sander_api.mod: interface.o
	$(VB)/bin/cp sander_api.mod $(INCDIR)

$(INCDIR)/sanderles_api.mod: interface.o
	$(VB)/bin/cp sanderles_api.mod $(INCDIR)

$(INCDIR)/sander.h: sander.h
	$(VB)/bin/cp sander.h $(INCDIR)

#---------------------------------------------------------------------------
libs:
#	$(VB)cd ../blas && $(MAKE) $(BLAS)
#	$(VB)cd ../lapack && $(MAKE) $(LAPACK)
#	$(VB)cd ../arpack && $(MAKE) $(ARPACK)
	$(VB)cd ../sqm && $(MAKE) libsqm
    
#---------------------------------------------------------------------------
Plumed.o: Plumed.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) \
	        -o Plumed.o Plumed.c

rand2.o:  rand2.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) -o $@ $<

lmodC.o:  lmodC.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

xminC.o:  xminC.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

ArrayList.o:  ArrayList.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

KFReader.o:  KFReader.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

qm2_read_adf_results.o:  qm2_read_adf_results.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

#---------------------------------------------------------------------------
clean:
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *__genmod.f90  # intel compiler-generated interface modules

#---------------------------------------------------------------------------
depend::
	./makedepend > depend

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
