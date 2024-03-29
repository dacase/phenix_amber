! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

module sgld

! MODULE: sgld
! ================== SELF-GUIDED MOLECULAR/LANGEVIN DYNAMICS ====================
! Xiongwu Wu, 2011

implicit none

!     Head file for the self-guided Langevin Dynamics simulation  
!
!      variables for SGLD simulation
!


! ... integers:
!
!
!    SGMD/SGLD applying range
!      ISGSTA      Begining atom index applying SGLD
!      ISGEND      Ending atom index applying SGLD
!
!
      integer,save:: ISGSTA,ISGEND,sgrndf,FIXCOM
!
! ... floats:
!
!
!    SGMD/SGLD VARIABLES
!     SGFT    !  Guiding factor 
!     TSGAVG  !  Local average time, ps
!     TEMPSG  !  Guiding temperature, K
!     TREFLF  !  Reference low frequency temperature, K
!     TSGAVP  !  Convergence time, ps
!
!
!     SGAVG0  !  Local average remains
!     SGAVG1  !  Local average factor, SGAVG1=1-SGAVG0
!     SGAVP0  !  Convergence average remains
!     SGAVP1  !  Convergency average factor, SGAVP1=1-SGAVP0
!     GAMMAS  !  friction coefficient
!     TEMPSG  !  Target guiding temperature 
!     TREFLF  !  Reference low frequency temperature
!

      _REAL_  SGFT,SGFF,SGFD,TSGAVG,TEMPSG,TSGSET,TREFLF,TSGAVP,GAMMAS,SGMASS, &
          SGAVG0,SGAVG1,SGAVP0,SGAVP1,SGFTI,SGFFI,SGFDI,TEMPSGI,EPOTLF,EPOTHF, &
          SGEFLF,SGEFHF,SGCFLF,AVGDD,AVGFF,AVGGD,AVGGF,AVGPP,AVGGP, &
          TEMPLFI,TEMPHFI,TREFLFI,TREFHFI,FRCLF,FRCHF,VIRSG,TEMPLLI,SGFTI2, &
          AVGELF,AVGTLF,AVGEFLF,AVGEFHF,AVGCFLF,AVGCFHF,SGWT, &
          TEMPRXLF,MYSCALSG,SGSCALE,FSGLDG

!  common block for parallel broadcast
      common/sgldr/SGFT,SGFF,SGFD,TSGAVG,TEMPSG,TSGSET,TREFLF,&
          SGFTI,SGFFI,SGFDI,TEMPSGI,SGMASS, &
          AVGDD,AVGFF,AVGGD,AVGGF,AVGPP,AVGGP, &
          TEMPLFI,TEMPHFI,TREFLFI,TREFHFI,TEMPLLI,SGFTI2, &
          AVGTLF,AVGEFLF,AVGEFHF,AVGCFLF,AVGCFHF, &
          TEMPRXLF,MYSCALSG,FSGLDG
!  Number of broadcasting variables
integer, parameter :: nsgld_real=32
       
!
! ... flags:
!
!
!     isgld   ! input control; positive values activate SGLD
!     TSGLDFP ! Perform SGLDfp to mantain canonical ensemble distribution
!     TSGLDG ! Perform SGLDg to mantain canonical ensemble distribution
!
!
      integer, save :: isgld
      LOGICAL, save, private :: TSGLDFP,TSGLDG
      LOGICAL, save :: TRXSGLD
      
!
! ...allocatable arrays:
!
!     VSG     ! local averages of velocity
!     DSG     ! interaction forces
!     FSG     ! local averages of force
!     GSG     ! guiding force
!     HSG     ! local averages of guiding force
!     SHKSG   ! constraint forces
!
      _REAL_, dimension(:), allocatable, private,save :: VSG,DSG,FSG,GSG,HSG,SHKSG


contains


    SUBROUTINE PSGLD(NATOM,AMASS,V, rem)
!-----------------------------------------------------------------------
!     This routine performs initiation for the Self-Guided        
!       Langevin Dynamcs (SGLD) simulaiton                 
!
      implicit none
#include "../include/md.h"
#ifdef MPI
#  include "parallel.h"
#endif
      INTEGER NATOM
      _REAL_ AMASS(*),V(*)
      integer rem
      INTEGER I,I3,M,ierror
      _REAL_ AMASSI,VI3,GAMM
      logical is_langevin  ! Is this a Langevin dynamics simulation
!
      TSGLDFP = (isgld == 2)
      tsgldg = (isgld == 3)
      trxsgld = rem > 0 .and. (isgld==1.or.isgld==3)
!  Check for invalid sgld setting
      IF(ISGSTA < 1)ISGSTA=1
      IF(ISGEND > NATOM .OR. ISGEND < 1)ISGEND=NATOM
      IF(TSGAVG.LT.DT)TSGAVG=DT
      IF(TSGAVP.LT.DT)TSGAVP=10.0D0*TSGAVG
      SGAVG1=DT/TSGAVG
      SGAVG0=1.0D0-SGAVG1
      SGAVP1=DT/TSGAVP
      SGAVP0=1.0D0-SGAVP1
      is_langevin = gamma_ln > 0.0d0
      SGFTI=SGFT
      SGFFI=SGFF
      SGFDI=SGFD
      IF(is_langevin)THEN
        IF(SGFT<-1.0D1)THEN
          SGFTI=1.0D0
          IF(TSGLDG)SGFTI=0.0D0
        ENDIF
        IF(FIXCOM<0)FIXCOM=0
      ELSE
        IF(SGFT<-1.0D1)THEN
          SGFTI=0.2D0
          IF(TSGLDG)SGFTI=1.0D0
        ENDIF
        IF(FIXCOM<0)FIXCOM=1
      ENDIF
      IF(SGFF<=-1.0D0)SGFFI=0.0D0
      IF(SGFD<=-1.0D0)SGFDI=0.0D0
      IF(TSGLDG)THEN
        IF(SGFTI<-0.0D0)SGFTI=1.0D0
        IF(TEMPSG<1.0D0)TEMPSG=temp0
      ENDIF
#ifdef MPI
      if(mytaskid.eq.0)THEN
#endif
        write(6,910)isgsta,isgend,tsgavg
        if(tempsg > 1.0d0)then
          write(6,920)tempsg
          if(tempsg<100.0d0)write(6,921)
          if(tsgldfp)then
            write(6,922)
          endif
        endif
        if(sgft > -1.0d2)then
          write(6,925)sgft
        endif
        if(sgff > -1.0d2)then
          write(6,926)sgff
        endif
        if(treflf < 1.0d-6)then
            if(tsgldfp)write(6,941)
        else
            if(tsgldfp)write(6,942)treflf
        endif
        if(is_langevin)then
          if(tsgldfp)then
            write(6,928)"  SGLDfp"
          else if(tsgldg)then
            write(6,929)"  SGLDg"
          else
            write(6,929)"  SGLD"
          endif
          write(6,930)gamma_ln
        else
          if(tsgldfp)then
            write(6,928)"  SGMDfp"
            if(treflf < 1.0d-6)then
              write(6,951)
              call mexit(6, 1)
            endif
          else if(tsgldg)then
            write(6,929)" SGMDg"
          else
            write(6,929)"  SGMD"
          endif
        endif
        if(fixcom>0)write(6,931)
        if(sgff>-1.0d1)then
            write(6,961)sgffi
        endif
        if(sgfd>-1.0d2)then
            write(6,962)sgfdi
        endif
        write(6,935)
#ifdef MPI
      ENDIF
#endif
      IF(GAMMA_LN > 0.0D0)THEN
        GAMMAS=GAMMA_LN/20.455d0
      ELSE
        GAMMAS=1.0D0/20.455d0
      ENDIF
      TSGSET=TEMP0
      IF(TSGSET<1.0D-6)TSGSET=MAX(TEMPSG,300.0D0)
      GAMM=SQRT(DT/TSGAVG)
      IF(TREFLF>1.0D-6)GAMM=SQRT(TREFLF/TSGSET)
!     Allocate working arrays
      allocate( VSG(3*natom ),DSG(3*natom ),FSG(3*natom ),&
                GSG(3*natom ),HSG(3*natom ),SHKSG(3*natom ),&
                stat=ierror)
      REQUIRE( ierror == 0 )
!    Initialize arrays
      SGMASS=0.0D0
      DO I=1,NATOM
        AMASSI = AMASS(I)
        SGMASS=SGMASS+AMASSI
        IF((I>=ISGSTA).AND.(I<=ISGEND))THEN
        ENDIF
        DO M=1,3
          I3=3*I-3+M
          VI3=V(I3)
          VSG(I3)=GAMM*VI3
          DSG(I3)=VI3
          FSG(I3)=0.0D0
          GSG(I3)=0.0D0
          HSG(I3)=0.0D0
          SHKSG(I3)=0.0D0
        END DO
        
      END DO
      temprxlf=0.0d0   !For RXSGLD
      TEMPSGI=TEMPSG
      IF(TEMPSG<1.0D-6)TEMPSGI=TSGSET
      TREFLFI=TREFLF
      IF(TREFLF<1.0D-6)TREFLFI=GAMM*TSGSET
      TEMPLFI=TREFLFI
      TREFHFI=TSGSET-TREFLFI
      TEMPHFI=TSGSET-TEMPLFI
      IF(TSGLDG)THEN
        TEMPLLI=GAMM*TSGSET
        IF(GAMMA_LN > 0.0D0)THEN
          SGFDI=(1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT((1.0D0-SGFTI)*(1.0D0+SGFFI))-1.0D0
          FSGLDG=((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT(1.0D0+SGFFI)-1.0D0)*SQRT(1.0D0-SGFTI)
        ELSE
          SGFTI2=0.0D0
          IF(SGFTI>0.0D0)SGFTI2=SQRT(SGFTI)
          FSGLDG=SGFTI2*((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT(1.0D0+SGFFI)-1.0D0)
          SGFDI=SGFTI2*((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0)) &
             *SQRT(1.0D0+SGFFI))
        ENDIF
      ENDIF
      EPOTLF=2.0D10
      SGWT=0.0D0
      VIRSG=0.0D0
      FRCLF=1.0D0
      FRCHF=1.0D0
      AVGELF=0.0D0
      AVGTLF=TREFLFI
      AVGEFLF=1.0D0
      AVGEFHF=1.0D0
      AVGCFLF=1.0D0
      AVGCFHF=1.0D0
      AVGFF=1.0D-6
      AVGDD=1.0D-6
      AVGGF=0.0D0
      AVGGD=0.0D0
      AVGPP=1.0D-6
      AVGGP=0.0D0
      SGEFLF=1.0D0
      SGEFHF=1.0D0
      SGCFLF=1.0D0
910   format("  _________________ SGLD parameters _________________"/  &
      "  Parameters for self-guided Langevin dynamics (SGLD) simulation"//  &
          "  Guiding range from ",i5,"  to ",i5 /  &
          "  Local averaging time: ",f10.4," ps ")
920   format("  Guiding temperature:",f8.2, " K" )
921   format("  *** WARNING: Guiding temperature is redefined since AMBER12 ***" / &
    "  Set tempsg > temp0 to enhance conformational search!"/ &
    "  tempsg defines a seaching ability comparable to a simulation at T=tempsg"/)
922   format("  *** WARNING: sgft, instead of tempsg, should be set for isgld == 2! ")
925   format("  Momentum guiding factor: ",f8.4)
926   format("  Force guiding factor: ",f8.4)
!927   format(" SGMDg is used for enhanced conformational search. ")
928   format(a8,"  method is used to mantain a canonical distribution. ")
929   format(a8,"  method is used to enhance conformational search. ")
930   format("  Collision frequency:",f8.2," /ps" )
931   format("  Translation of COM is freezed!" )
!933   format("  Guiding temperature is not defined.  Set tempsg=:",f8.2, " K" )
935   format("  SGMD/SGLD output properties:"    /  &
             "  SGLF=  SGFT   TEMPSG   TEMPLF   TREFLF   FRCLF   EPOTLF    SGWT" /  &
             "  SGHF=  SGFF   SGFD     TEMPHF   TREFHF   FRCHF   EPOTHF   VIRSG" /  &
             "         SGMD/SGLD weighting factor=exp(SGWT)"/  &
              " _______________________________________________________"/)
941   format("  WARNING: treflf is not defined and will be estimated from simulation. ")
942   format("  treflf= ",f8.2," K is input for guiding temperature calculation.")
951   format("  WARNING: treflf must be defined for SGMDfp. ")
961   format("  sgff is fixed at: ",f8.4)
962   format("  sgfd is fixed at: ",f8.4)
      RETURN
      END SUBROUTINE PSGLD

    subroutine sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
!-----------------------------------------------------------------------
!     Correct the total number of degrees of freedom for a translatable COM,
!       and compute the number of degrees of freedom in the SGLD part.
!       The latter is mostly done by the caller to avoid passing the long
!       argument list needed by routine degcnt which also differs between
!       sander and pmemd.
!
      implicit none
      _REAL_, intent(in)    :: sgsta_rndfp, sgend_rndfp
      integer, intent(in)   :: ndfmin
      _REAL_, intent(inout) :: rndf

      sgrndf = int(sgend_rndfp - sgsta_rndfp)
      if (fixcom > 0 .and. ndfmin == 0) then
         rndf = rndf - 3
         sgrndf = sgrndf - 3
      end if
      return
      end subroutine sg_fix_degree_count

    SUBROUTINE SGFSHAKE(ISTART,IEND,DT,AMASS,X,QCALC)
!-----------------------------------------------------------------------
!     This routine calculate SHAKE constraint forces needed for SGLD reweighting                 
!
      implicit none
      INTEGER ISTART,IEND,I,I3,M
      LOGICAL QCALC
      _REAL_ DT,AMASS(*),X(*)
      _REAL_ DT2,FACT
      IF(QCALC)THEN
!  Calculate shake forces
        DT2=1.0D0/(DT*DT)
        DO I=ISTART,IEND
          FACT=DT2*AMASS(I)
          I3=(I-1)*3
          DO M=1,3
            I3=I3+1
            SHKSG(I3)=FACT*(X(I3)-DSG(I3))
            !FSG(I3)=FSG(I3)+SGAVG1*SHKSG(I3)
          ENDDO
        ENDDO
      ELSE
!  Save old positions
        DO I=ISTART,IEND
          I3=(I-1)*3
          DO M=1,3
            I3=I3+1
            DSG(I3)=X(I3)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SGFSHAKE


      SUBROUTINE SGLDW(NATOM,ISTART,IEND,NTP, &
             DTX,TEMP0,ENER,AMASS,WINV,X,F,V)
!-----------------------------------------------------------------------
!     This routine perform SGLD integration        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
   include 'mpif.h'
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
#endif
      INTEGER NATOM,ISTART,IEND,NTP
      _REAL_ DTX,TEMP0
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),X(*),F(*),V(*)
!
      INTEGER I,M,I3,JSTA,JEND
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC,DRAGI,DRAGJ,PV1,PV2,RSD,FLN
      _REAL_ PTOT(3),FTOT(3),EKIN,EKINSG,EKINSGG
      _REAL_ SUMDD,SUMFF,SUMGD,SUMGF,SUMPP,SUMGP
      _REAL_ VIT,VI3,FI3,GI3,HI3,VSGI,DSGI,FSGI,PSGI,FRICI
      PARAMETER (BOLTZ = 1.987192d-3)
!
        JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
!   Estimate guiding factor or guiding temperature
        DO M=1,3
          PTOT(M)=0.0D0
          FTOT(M)=0.0D0
        ENDDO
        EKIN=0.0D0
        EKINSG=0.0D0
        EKINSGG=0.0D0
        PV1=0.0d0 
        PV2=0.0d0
        DO  I = 1,NATOM 
          AMASSI = AMASS(I)
          WFAC =  DTX*0.5D0*WINV(I)
          RSD = SQRT(2.D0*GAMMAS*BOLTZ*TEMP0*AMASSI/DTX)
          I3 = 3*(I-1)
          DO  M = 1,3
!   Keep random number series the same as that in a single cpu simulation
            CALL GAUSS( 0.D0, RSD, FLN )
            I3 = 3*(I-1)+M
            VI3=V(I3)
            IF(I>=JSTA.AND.I<=JEND)THEN
              VSGI=SGAVG0*VSG(I3)+SGAVG1*VI3
              VSG(I3)=VSGI
              PSGI=GAMMAS*AMASSI*VSGI
              FI3=F(I3)
              IF(TSGLDG)THEN
                GI3=SGAVG0*GSG(I3)+SGAVG1*VSGI
                GSG(I3)=GI3
                DSGI=SGAVG0*HSG(I3)+SGAVG1*FLN
                HSG(I3)=DSGI
                DRAGJ=FI3+SHKSG(I3)
                FSGI=SGAVG0*FSG(I3)+SGAVG1*DRAGJ
                FSG(I3)=FSGI
                EKINSGG=EKINSGG+AMASSI*GI3*GI3
                DRAGI=SGFTI*PSGI+SGFFI*FSGI+SGFDI*DSGI
                DRAGJ=SGFFI*FSGI+FSGLDG*DSGI
              ELSE
                DSGI=FI3+SHKSG(I3)
                FSGI=SGAVG0*FSG(I3)+SGAVG1*DSGI
                FSG(I3)=FSGI
                DSGI=DSGI-FSGI
                SHKSG(I3)=DSGI
                DRAGI=SGFTI*PSGI+SGFFI*FSGI+SGFDI*DSGI
                DRAGJ=DRAGI
              ENDIF
              DSG(I3)=DRAGJ
              FI3=FI3+fln+DRAGI
              F(I3)=FI3
              VIT = VI3 + FI3*wfac
              PV1=PV1+VIT*DRAGJ
              PV2=PV2+AMASSI*VIT*VIT
              EKIN=EKIN+AMASSI*VI3*VI3
              EKINSG=EKINSG+AMASSI*VSGI*VSGI
              PTOT(M)=PTOT(M)+AMASSI*VI3
              FTOT(M)=FTOT(M)+FI3
            ELSE IF(I>=ISTART.AND.I<=IEND)THEN
              FI3=F(I3)+fln
              F(I3)=FI3
              PTOT(M)=PTOT(M)+AMASSI*VI3
              FTOT(M)=FTOT(M)+FI3
            ENDIF
          END DO
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=PV1
          TEMP1(2)=PV2
          TEMP1(3)=PTOT(1)
          TEMP1(4)=PTOT(2)
          TEMP1(5)=PTOT(3)
          TEMP1(6)=FTOT(1)
          TEMP1(7)=FTOT(2)
          TEMP1(8)=FTOT(3)
          TEMP1(9)=EKIN
          TEMP1(10)=EKINSG
          TEMP1(11)=EKINSGG
          call mpi_allreduce(MPI_IN_PLACE,temp1,11,&
             MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          PV1=TEMP1(1)
          PV2=TEMP1(2)
          PTOT(1)=TEMP1(3)
          PTOT(2)=TEMP1(4)
          PTOT(3)=TEMP1(5)
          FTOT(1)=TEMP1(6)
          FTOT(2)=TEMP1(7)
          FTOT(3)=TEMP1(8)
          EKIN=TEMP1(9)
          EKINSG=TEMP1(10)
          EKINSGG=TEMP1(11)
        ENDIF
#endif
        IF(FIXCOM>0)THEN
          DO M=1,3
            PTOT(M)=PTOT(M)/SGMASS
            FTOT(M)=FTOT(M)/SGMASS
          ENDDO
        ELSE
          DO M=1,3
            PTOT(M)=0.0D0
            FTOT(M)=0.0D0
          ENDDO
        ENDIF
        sgscale=(1.0d0+0.5d0*gammas*dtx)*pv1/(pv2-0.5d0*dtx*pv1)
        SUMDD=0.0d0
        SUMFF=0.0d0
        SUMGD=0.0d0
        SUMGF=0.0d0
        SUMPP=0.0d0
        SUMGP=0.0d0
        DO  I = ISTART,IEND 
          WFAC = WINV(I)*DTX
          AMASSI=AMASS(I)
          IF(I<JSTA.OR.I>JEND)THEN
            FACT=0.5D0*DTX*GAMMAS
          ELSE
            FACT=0.5D0*DTX*(GAMMAS+SGSCALE)
          ENDIF
          DO  M = 1,3
            I3 = 3*(I-1)+M
            VI3=V(I3)-PTOT(M)
            FI3=F(I3)-AMASSI*FTOT(M)
            F(I3)=FI3
            VIT=((1.0D0-FACT)*VI3+FI3*WFAC)/(1.0D0+FACT)
            V(I3)=VIT
            VIT=0.5D0*(VI3+VIT)
            VSGI=VSG(I3)
            DSGI=SHKSG(I3)
            FSGI=FSG(I3)
            FRICI=AMASSI*SGSCALE*VIT
            PSGI=AMASSI*GAMMAS*VSGI
            FI3=SGFTI*PSGI-FRICI
            DRAGI=DSG(I3)-FRICI
            DSG(I3)=DRAGI
            IF(TSGLDG)CONTINUE
            GI3=SGAVG0*GSG(I3)+SGAVG1*FI3
            HI3=SGAVG0*HSG(I3)+SGAVG1*DRAGI
            GSG(I3)=GI3
            HSG(I3)=HI3
            SUMDD=SUMDD+DSGI*DSGI
            SUMFF=SUMFF+FSGI*FSGI
            SUMGD=SUMGD+(FI3-GI3)*DSGI
            SUMGF=SUMGF+(GI3)*FSGI
            SUMPP=SUMPP+PSGI*PSGI
            SUMGP=SUMGP+(HI3)*PSGI
          END DO
        END DO
    ! Estimate the viral due to the guiding forces
        IF(NTP>0)THEN
           VIRSG=0.0D0
           DO  I = ISTART,IEND 
             DO  M = 1,3
               I3 = 3*(I-1)+M
               VIRSG=VIRSG+X(I3)*DSG(I3)
             ENDDO
           ENDDO
           VIRSG=VIRSG/3.0D0
        ENDIF
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=SUMDD
          TEMP1(2)=SUMFF
          TEMP1(3)=SUMGD
          TEMP1(4)=SUMGF
          TEMP1(5)=SUMPP
          TEMP1(6)=SUMGP
          TEMP1(7)=VIRSG
          call mpi_allreduce(MPI_IN_PLACE,temp1,7,&
             MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          SUMDD=TEMP1(1)
          SUMFF=TEMP1(2)
          SUMGD=TEMP1(3)
          SUMGF=TEMP1(4)
          SUMPP=TEMP1(5)
          SUMGP=TEMP1(6)
          VIRSG=TEMP1(7)
        ENDIF
#endif
    ! Estimate low frequency temperatures
        TEMPI=EKIN/sgrndf/BOLTZ
        TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPHFI=TSGSET-TEMPLFI
        IF(TSGLDG)THEN
          TEMPLLI=SGAVP0*TEMPLLI+SGAVP1*EKINSGG/sgrndf/BOLTZ
          SGFDI=(1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT((1.0D0-SGFTI)*(1.0D0+SGFFI))-1.0D0
          FSGLDG=((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT(1.0D0+SGFFI)-1.0D0)*SQRT(1.0D0-SGFTI)
          FRCLF=1.0D0+SGFFI
          FRCHF=1.0D0
          SGCFLF=TSGSET/TEMPSG
          TREFLFI=TEMPLFI*SGCFLF
          TREFHFI=TSGSET-TREFLFI
          SGSCALE=SGFDI
          ! update accumulators
          CALL SGENERGY(ENER)
          RETURN
        ENDIF
    ! Estimate momentum guiding factor and guiding temperature
        AVGFF=SGAVP0*AVGFF+SGAVP1*SUMFF
        AVGDD=SGAVP0*AVGDD+SGAVP1*SUMDD
        AVGGF=SGAVP0*AVGGF+SGAVP1*SUMGF
        AVGGD=SGAVP0*AVGGD+SGAVP1*SUMGD
        AVGPP=SGAVP0*AVGPP+SGAVP1*SUMPP
        AVGGP=SGAVP0*AVGGP+SGAVP1*SUMGP
    ! Estimate SGLD factors
        SGEFLF=1.0D0+AVGGF/AVGFF
        SGEFHF=1.0D0+AVGGD/AVGDD
        SGCFLF=1.0D0-AVGGP/AVGPP
    ! Estimate reference temperatures
        TREFLFI=TREFLF
        IF(TREFLF<1.0D-5)THEN
          TREFLFI=MAX(TEMPLFI*SGCFLF,TEMPLFI/1.0D1)
          IF(TRXSGLD.AND.TEMPRXLF>1.0D-5)TREFLFI=TEMPRXLF
        ENDIF
        TREFHFI=TSGSET-TREFLFI
    ! Estimate weighting factor
        frclf=sgeflf+sgffi
        frchf=sgefhf+sgfdi
    ! Estimate guiding temperatures
        TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    ! Adjust guiding factors if wanted
        IF(TEMPSG>1.0D-5)THEN
           SGFTI=SGFTI+SGAVP1*(TEMPSG-TEMPSGI)*(TEMPSG-TSGSET) &
                                  /TSGSET/(TEMPSG+TEMPSGI)
           SGFTI=MIN(SGFTI,1.0D1)
           SGFTI=MAX(SGFTI,-1.0D1)
        ENDIF
        IF(TSGLDFP)THEN
    ! Estimate force guiding factor
          IF(SGFD<-1.0D2 .AND. SGFF<-1.0D2 )THEN
            SGFDI=SGAVP0*SGFDI+SGAVP1*(TEMPHFI/TREFHFI-SGEFHF)
            SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
          ELSE IF(SGFD<-1.0D2)THEN
            SGFDI=SGAVG0*SGFDI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF-SGFF)
          ELSE IF(SGFF<-1.0D2)THEN
            SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
          ENDIF
          SGSCALE=SGFDI
        ELSE
          ! Convert SGSCALE to unitless
          SGSCALE=SGSCALE/GAMMAS
        ENDIF
        ! update accumulators
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGLDW

        SUBROUTINE SGMDW(NATOM,ISTART,IEND,NTP, &
             DTX,ENER,AMASS,WINV,X,F,V)
!-----------------------------------------------------------------------
!     This routine calculate guiding force using SGLD method 
!     for MD simulation        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
   include 'mpif.h'
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
#endif
      INTEGER NATOM,ISTART,IEND,NTP
      _REAL_ DTX
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),X(*),F(*),V(*)
!
      INTEGER JSTA,JEND,I,M,I3
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC,FRICI,DRAGI,DRAGJ,PV1,PV2
      _REAL_ PTOT(3),FTOT(3),EKIN,EKINSG,EKINSGG
      _REAL_ SUMDD,SUMFF,SUMGD,SUMGF
      _REAL_ VI3,VIT,FI3,GI3,HI3,VSGI,DSGI,FSGI,PSGI,RSD,FLN
      PARAMETER (BOLTZ = 1.987192d-3)
!
        JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
!    Estimate guiding factor or guiding temperature
        DO M=1,3
          PTOT(M)=0.0D0
          FTOT(M)=0.0D0
        ENDDO
        EKIN=0.0D0
        EKINSG=0.0D0
        EKINSGG=0.0D0
        PV1=0.0d0 
        PV2=0.0d0
        DO  I = 1,NATOM 
          AMASSI = AMASS(I)
          WFAC =  DTX*0.5D0*WINV(I)
          RSD = SQRT(2.D0*GAMMAS*BOLTZ*TSGSET*AMASSI/DTX)
          DO  M = 1,3
            I3 = 3*(I-1)+M
            VI3=V(I3)
            FI3=F(I3)
!   Keep random number series the same as that in a single cpu simulation
            CALL GAUSS( 0.D0, RSD, FLN )
            IF(I>=JSTA.AND.I<=JEND)THEN
              VSGI=SGAVG0*VSG(I3)+SGAVG1*VI3
              VSG(I3)=VSGI
              PSGI=GAMMAS*AMASSI*VSGI
              IF(TSGLDG)THEN
                GI3=SGAVG0*GSG(I3)+SGAVG1*VSGI
                GSG(I3)=GI3
                DSGI=SGAVG0*HSG(I3)+SGAVG1*FLN
                HSG(I3)=DSGI
                DRAGJ=FI3+SHKSG(I3)
                FSGI=SGAVG0*FSG(I3)+SGAVG1*DRAGJ
                FSG(I3)=FSGI
                EKINSGG=EKINSGG+AMASSI*GI3*GI3
                DRAGI=-SGFTI*PSGI+SGFFI*FSGI+SGFDI*DSGI
                DRAGJ=SGFFI*FSGI+FSGLDG*DSGI
              ELSE
                DSGI=FI3+SHKSG(I3)
                FSGI=SGAVG0*FSG(I3)+SGAVG1*DSGI
                FSG(I3)=FSGI
                DSGI=DSGI-FSGI
                SHKSG(I3)=DSGI
                DRAGI=SGFTI*PSGI+SGFFI*FSGI+SGFDI*DSGI
                DRAGJ=DRAGI
              ENDIF
              DSG(I3)=DRAGJ
              FI3=FI3+DRAGI
              F(I3)=FI3
              VIT = VI3 + FI3*WFAC
              PV1=PV1+VIT*DRAGJ
              PV2=PV2+AMASSI*VIT*VIT
              EKIN=EKIN+AMASSI*VI3*VI3
              EKINSG=EKINSG+AMASSI*VSGI*VSGI
              PTOT(M)=PTOT(M)+AMASSI*VI3
              FTOT(M)=FTOT(M)+FI3
            ELSE IF(I>=ISTART.AND.I<=IEND)THEN
              PTOT(M)=PTOT(M)+AMASSI*VI3
              FTOT(M)=FTOT(M)+FI3
            ENDIF
          END DO
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=PV1
          TEMP1(2)=PV2
          TEMP1(3)=PTOT(1)
          TEMP1(4)=PTOT(2)
          TEMP1(5)=PTOT(3)
          TEMP1(6)=FTOT(1)
          TEMP1(7)=FTOT(2)
          TEMP1(8)=FTOT(3)
          TEMP1(9)=EKIN
          TEMP1(10)=EKINSG
          TEMP1(11)=EKINSGG
          call mpi_allreduce(MPI_IN_PLACE,temp1,11, &
            MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          PV1=TEMP1(1)
          PV2=TEMP1(2)
          PTOT(1)=TEMP1(3)
          PTOT(2)=TEMP1(4)
          PTOT(3)=TEMP1(5)
          FTOT(1)=TEMP1(6)
          FTOT(2)=TEMP1(7)
          FTOT(3)=TEMP1(8)
          EKIN=TEMP1(9)
          EKINSG=TEMP1(10)
          EKINSGG=TEMP1(11)
        ENDIF
#endif
        IF(FIXCOM>0)THEN
          DO M=1,3
            PTOT(M)=PTOT(M)/SGMASS
            FTOT(M)=FTOT(M)/SGMASS
          ENDDO
        ELSE
          DO M=1,3
            PTOT(M)=0.0D0
            FTOT(M)=0.0D0
          ENDDO
        ENDIF
        SGSCALE=PV1/(PV2-0.5D0*DTX*PV1)
        FACT=SGSCALE/(1.0D0+0.5D0*SGSCALE*DTX)
        SUMDD=0.0d0
        SUMFF=0.0d0
        SUMGD=0.0d0
        SUMGF=0.0d0
        DO  I = ISTART,IEND 
          WFAC = 0.5D0*WINV(I)*DTX
          AMASSI=AMASS(I)
          DO  M = 1,3
            I3 = 3*(I-1)+M
            VI3 = V(I3)-PTOT(M)
            V(I3)=VI3
            FI3=F(I3)-AMASSI*FTOT(M)
            IF(I<JSTA.OR.I>JEND)THEN
              F(I3)=FI3
              CYCLE
            ENDIF
            VIT = VI3 + FI3*WFAC
            FRICI=FACT*AMASSI*VIT
            F(I3)=FI3-FRICI
            DRAGI=DSG(I3)-FRICI
            DSG(I3)=DRAGI
            IF(TSGLDG)CONTINUE
            VSGI=VSG(I3)
            DSGI=SHKSG(I3)
            FSGI=FSG(I3)
            PSGI=AMASSI*GAMMAS*VSGI
            FI3=SGFTI*PSGI-FRICI
            GI3=SGAVG0*GSG(I3)+SGAVG1*FI3
            HI3=SGAVG0*HSG(I3)+SGAVG1*DRAGI
            GSG(I3)=GI3
            HSG(I3)=HI3
            SUMDD=SUMDD+DSGI*DSGI
            SUMFF=SUMFF+FSGI*FSGI
            SUMGD=SUMGD+(FI3-GI3)*DSGI
            SUMGF=SUMGF+(GI3)*FSGI
          END DO
        END DO
    ! Estimate the viral due to the guiding forces
        IF(NTP>0)THEN
           VIRSG=0.0D0
           DO  I = JSTA,JEND 
             DO  M = 1,3
               I3 = 3*(I-1)+M
               VIRSG=VIRSG+X(I3)*DSG(I3)
             ENDDO
           ENDDO
           VIRSG=VIRSG/3.0D0
        ENDIF
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=SUMDD
          TEMP1(2)=SUMFF
          TEMP1(3)=SUMGD
          TEMP1(4)=SUMGF
          TEMP1(5)=VIRSG
          call mpi_allreduce(MPI_IN_PLACE,temp1,5, &
            MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          SUMDD=TEMP1(1)
          SUMFF=TEMP1(2)
          SUMGD=TEMP1(3)
          SUMGF=TEMP1(4)
          VIRSG=TEMP1(5)
        ENDIF
#endif
    ! Estimate low frequency temperatures
        TEMPI=EKIN/sgrndf/BOLTZ
        TEMPLFI=SGAVP0*TEMPLFI+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPHFI=TSGSET-TEMPLFI
        IF(TSGLDG)THEN
          TEMPLLI=SGAVP0*TEMPLLI+SGAVP1*EKINSGG/sgrndf/BOLTZ
          FSGLDG=SGFTI2*((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0))  &
             *SQRT(1.0D0+SGFFI)-1.0D0)
          SGFDI=SGFTI2*((1.0D0+SQRT(TEMPLFI/TEMPLLI)*(SQRT(TEMPSG/TSGSET)-1.0D0)) &
             *SQRT(1.0D0+SGFFI))
          FRCLF=1.0D0+SGFFI
          FRCHF=1.0D0
          SGCFLF=TSGSET/TEMPSG
          TREFLFI=TEMPLFI*SGCFLF
          TREFHFI=TSGSET-TREFLFI
          SGSCALE=SGFDI
          ! update accumulators
          CALL SGENERGY(ENER)
          RETURN
        ENDIF
    ! Estimate momentum guiding factor and guiding temperature
        AVGFF=SGAVP0*AVGFF+SGAVP1*SUMFF
        AVGDD=SGAVP0*AVGDD+SGAVP1*SUMDD
        AVGGF=SGAVP0*AVGGF+SGAVP1*SUMGF
        AVGGD=SGAVP0*AVGGD+SGAVP1*SUMGD
    ! Estimate SGLD factors
        SGEFLF=SGAVP0*SGEFLF+SGAVP1*(AVGGF/AVGFF+1.0D0)
        SGEFHF=SGAVP0*SGEFHF+SGAVP1*(AVGGD/AVGDD+1.0D0)
    ! Estimate reference temperatures
        TREFLFI=TREFLF
        IF(TREFLF<1.0D-5)THEN
          TREFLFI=TEMPLFI
          IF(TRXSGLD.AND.TEMPRXLF>1.0D-5)TREFLFI=TEMPRXLF
        ENDIF
        TREFHFI=TSGSET-TREFLFI
    ! Estimate guiding temperatures
        TEMPSGI=SGAVP0*TEMPSGI+SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
    ! Estimate weighting factor
        FRCLF=SGEFLF+SGFFI
        FRCHF=SGEFHF+SGFDI
    ! Adjust guiding factors if wanted
        IF(TEMPSG>1.0D-5)THEN
           SGFTI=SGFTI+SGAVP1*(TEMPSG-TEMPSGI)*(TEMPSG-TSGSET) &
                                  /TSGSET/(TEMPSG+TEMPSGI)
           SGFTI=MIN(SGFTI,1.0D1)
           SGFTI=MAX(SGFTI,-1.0D1)
        ENDIF
        IF(TSGLDFP)THEN
    ! Estimate force guiding factor
          IF(SGFD<-1.0D2 .AND. SGFF<-1.0D2 )THEN
            SGFDI=SGAVP0*SGFDI+SGAVP1*(TEMPHFI/TREFHFI-SGEFHF)
            SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
          ELSE IF(SGFD<-1.0D2)THEN
            SGFDI=SGAVG0*SGFDI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF-SGFF)
          ELSE IF(SGFF<-1.0D2)THEN
            SGFFI=SGAVP0*SGFFI+SGAVP1*(TEMPLFI/TREFLFI-SGEFLF)
          ENDIF
          SGSCALE=SGFDI
        ELSE
          ! Convert SGSCALE to unitless
          SGSCALE=SGSCALE/GAMMAS
        ENDIF
        ! update accumulators
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGMDW
     
        SUBROUTINE SGENERGY(ENER)
!-----------------------------------------------------------------------
!     This routine set the ener fields of SGLD variables
!
      use state
      implicit none
      type(state_rec) :: ener
      _REAL_ BOLTZ
      PARAMETER (BOLTZ = 1.987192d-3)
      _REAL_ EPOTI
    ! Weighting accumulators
        EPOTI=ENER%POT%TOT
        IF(EPOTLF>1.0D10)THEN
          EPOTLF=EPOTI
          AVGEFLF=FRCLF
          AVGEFHF=FRCHF
          AVGCFLF=TREFLFI/TEMPLFI
          AVGCFHF=TREFHFI/(TSGSET-TEMPLFI)
          AVGELF=EPOTLF
          AVGTLF=TEMPLFI
        ELSE
          EPOTLF=SGAVG0*EPOTLF+SGAVG1*EPOTI
          AVGEFLF=SGAVP0*AVGEFLF+SGAVP1*FRCLF
          AVGEFHF=SGAVP0*AVGEFHF+SGAVP1*FRCHF
          AVGCFLF=SGAVP0*AVGCFLF+SGAVP1*TREFLFI/TEMPLFI
          AVGCFHF=SGAVP0*AVGCFHF+SGAVP1*TREFHFI/(TSGSET-TEMPLFI)
          AVGELF=SGAVP0*AVGELF+SGAVP1*EPOTLF
          AVGTLF=SGAVP0*AVGTLF+SGAVP1*TEMPLFI
        ENDIF
        EPOTHF=EPOTI-EPOTLF
        SGWT=((AVGEFLF*AVGCFLF-1.0D0)*(EPOTLF-AVGELF)+  &
         (AVGEFHF*AVGCFHF-1.0D0)*EPOTHF+VIRSG)/(BOLTZ*TSGSET)
    ! Update ENER structure
        ENER%SGLD%SGFT=SGFTI
        ENER%SGLD%SGFF=SGFFI
        ENER%SGLD%SGSCAL=SGSCALE
        ENER%SGLD%TEMPSG=TEMPSGI
        ENER%SGLD%TEMPLF=TEMPLFI
        ENER%SGLD%TEMPHF=TEMPHFI
        ENER%SGLD%TREFLF=TREFLFI
        ENER%SGLD%TREFHF=TREFHFI
        ENER%SGLD%FRCLF=FRCLF
        ENER%SGLD%FRCHF=FRCHF
        ENER%SGLD%EPOTLF=EPOTLF
        ENER%SGLD%EPOTHF=EPOTHF
        ENER%SGLD%VIRSG=VIRSG
        ENER%SGLD%SGWT=SGWT
       RETURN
       END SUBROUTINE SGENERGY

#ifdef MPI

!*********************************************************************
!               SUBROUTINE SGLD_EXCHG
!*********************************************************************
!  exchange all data in the SGLDR common block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sgld_exchg(irep)

   implicit none
   include 'mpif.h'
#  include "parallel.h"

   integer, intent(in) :: irep
   integer  ierror, istatus(mpi_status_size)
   call mpi_sendrecv_replace(sgft,nsgld_real, mpi_double_precision, &
                   irep, 511, irep, 511, commmaster, istatus, ierror)
    !call mpi_barrier(commmaster, ierror)
   return
   end subroutine sgld_exchg

!*********************************************************************
!               SUBROUTINE REMD_SCALE_VELO
!*********************************************************************
! Scale velocities based on new temps after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rxsgld_scale(stagid,nr,myscaling,amass,v)

   implicit none
   include 'mpif.h'
#  include "parallel.h"

   integer, intent(in) :: stagid,nr
   _REAL_, intent(in) :: myscaling
   _REAL_, dimension (*), intent(inout) :: v
   _REAL_, dimension (*), intent(in)    :: amass
   integer ierror
   integer i,j,jsta,jend
   _REAL_ ek,elf,ehf,elh,amassi,vi,vlf,vhf,chk,hfscale
   _REAL_ temp1(10),temp2(10)
!--------------------
         !if (sanderrank==0) then
         !   write (6,'(a,i4,2x,f8.3,a,2f8.3)') &
         !      "RXSGLD: stagid, scalsg  ",stagid,myscalsg,&
         !      " to match sgft,tempsg,: ",sgft,tempsg
         !endif
      call mpi_bcast(stagid,1,mpi_integer,0,commsander,ierror)

         ! All processes scale velocities.
         ! DAN ROE: This could potentially be divided up as in runmd
         !  since when there are mutiple threads per group each thread 
         !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
         if (sanderrank==0) then
            write (6,'(a,f8.3,a,f8.3)') &
               "| RXSGLD: scaling guiding properties by ",myscalsg,&
               " to match a new guiding factor ",sgfti
         endif
#endif
! ---=== Broadcast RXSGLD guiding effect ===---
      IF(SANDERSIZE > 1)call mpi_bcast(sgft,nsgld_real,mpi_double_precision,&
                                              0,commsander,ierror)
      if (myscalsg > 0.0d0) then
        jsta = iparpt(mytaskid) + 1
        jend = iparpt(mytaskid+1)
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
         ek=0.0d0
         ehf=0.0d0
         elf=0.0d0
         elh=0.0d0
         do i = jsta,jend
            amassi=amass(i)
            do j=3*i-2,3*i
               vi=v(j)
               vlf=vsg(j)
               vhf=vi-vlf
               ek=ek+amassi*vi*vi
               ehf=ehf+amassi*vhf*vhf
               elf=elf+amassi*vlf*vlf
               elh=elh+amassi*vhf*vlf
            enddo
         enddo
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=ek
          TEMP1(2)=ehf
          TEMP1(3)=elf
          TEMP1(4)=elh
          call mpi_allreduce(MPI_IN_PLACE,temp1,4,&
             MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierror)
          ek=TEMP1(1)
          ehf=TEMP1(2)
          elf=TEMP1(3)
          elh=TEMP1(4)
        ENDIF
         ! solve high frequency scaling factor
         chk=elh*elh+ehf*(myscaling*myscaling*ek-myscalsg*myscalsg*elf)
         hfscale=1.0d0
         if(chk>0.0d0)hfscale=(sqrt(chk)-elh)/ehf
         hfscale=max(hfscale,0.5d0)
         hfscale=min(hfscale,2.0d0)
         !write(6,*)"scales:",stagid,sanderrank, &
         !        myscaling,myscalsg,hfscale,ek,elf,ehf
         do i = 1,nr
           do j=3*i-2,3*i
             if(i<jsta.or.i>jend)then
               v(j)=myscaling*v(j)
             else
               vsg(j) = vsg(j) * myscalsg
               gsg(j) = gsg(j) * myscalsg
               hsg(j) = hsg(j) * myscalsg
               v(j)=hfscale*v(j)+(myscalsg-hfscale)*vsg(j)
             endif
           enddo
         enddo
      endif
   return

end subroutine rxsgld_scale

!*********************************************************************
!               FUNCTION TEMPSGLOOKUP
!*********************************************************************
! lookup temp in templist and return its index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function tempsglookup(numreps,temp,tsg,sgft,temps,tsgs,sgfts)

   implicit none
#  include "parallel.h"

   integer numreps
   _REAL_, intent(in) :: temp,tsg,sgft
   _REAL_, dimension(numreps), intent(in) :: temps,tsgs,sgfts

   integer i
   
   tempsglookup=0
   do i = 1, numreps
      if(abs(temp-temps(i)) < 1.0d-6 &
      .and. abs(tsg-tsgs(i)) < 1.0d-6 &
      .and. abs(sgft-sgfts(i)) < 1.0d-6) then
         if(tempsglookup>0)then
            write (6,*) "================================"
            write (6,*) "Two replicas are the same: ",tempsglookup,i
            write (6,*) "================================"
            call mexit(6,1)
         endif
         tempsglookup = i
      end if
   end do
   return
end function tempsglookup

!*********************************************************************
!               FUNCTION STAGIDLOOKUP
!*********************************************************************
! lookup stagid in stagidlist and return its neighboring index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function stagidlookup(numreps,id,idtable)

   implicit none
#  include "parallel.h"

   integer numreps
   INTEGER, dimension(numreps), intent(in) :: idtable
   INTEGER, intent(in) :: id

   integer i
   
   stagidlookup=-1
   do i = 1, numreps
      if(id==idtable(i)) stagidlookup=i
   end do
   return
end function stagidlookup

!*********************************************************************
!               SUBROUTINE SORTTEMPSG
!*********************************************************************
! sort temp ascendingly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sorttempsg(numreps,temp,tsgs,sgfts)

   implicit none

#  include "parallel.h"

   integer numreps
   _REAL_, dimension(numreps), intent(inout) :: temp,tsgs,sgfts
   _REAL_, dimension(numreps) :: tmp
   INTEGER, dimension(numreps) :: tmpid

   _REAL_ tempt
   integer i, j, ii

   do i = 1, numreps
     tmp(i)=1000000*temp(i)+100*tsgs(i)+sgfts(i)
     tmpid(i)=i
   enddo
   do i = 1, numreps
      do j = i + 1, numreps
         if(tmp(j) < tmp(i)) then
            tempt = tmp(i)
            tmp(i) = tmp(j)
            tmp(j) = tempt
            ii=tmpid(i)
            tmpid(i)=tmpid(j)
            tmpid(j)=ii
         end if
      end do
   end do
   do i = 1, numreps
     ii=tmpid(i)
     tmp(i)=temp(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     temp(i)=tmp(i)
     tmp(i)=tsgs(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     tsgs(i)=tmp(i)
     tmp(i)=sgfts(ii)
   enddo
   do i = 1, numreps
     sgfts(i)=tmp(i)
   enddo
   return
end subroutine sorttempsg

#endif /* MPI */

end module sgld

