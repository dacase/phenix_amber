!<compile=optimized>
#include "nfe-utils.h"
#include "nfe-config.h" 

! By Sishi Tang and Lin Fu
! June 17, 2011 
! This modules contains functions that handles I/O for PCA 
! In order to calculate PCA, we need: 
! cv_value = Transpose(A) * (x(fitted) - xavg) 
!          = Transpose(A) * (R(x-x(cm)) - xavg)  

module nfe_read_pca

implicit none 

private 

!===================================================================

public :: read_evec 
public :: read_avgcrd
public :: read_refcrd
public :: read_index  

!===================================================================

contains 

!===================================================================

subroutine read_evec(cv, evec_file) 

  use nfe_colvar_type
  use nfe_constants 

  implicit none 
 
  type(colvar_t), intent(inout)  :: cv  
 
  character(len = *), intent(in)  :: evec_file 
 
!  integer, intent(in) :: first, last 
  integer :: i 
  character :: dummy 

  OPEN(EVEC_UNIT1, FILE = evec_file) 
  ! skip the first two lines 
  read(EVEC_UNIT1, '(A1)') dummy 
  read(EVEC_UNIT1, '(A1)') dummy 
  
  read(EVEC_UNIT1, '(7F11.5)') (cv%evec(i),i=1, 3*cv%i(3)) 

  close(EVEC_UNIT1)

end subroutine read_evec 

subroutine read_avgcrd(cv, avgcrd_file)

  use nfe_colvar_type 
  use nfe_constants

  implicit none

  type(colvar_t), intent(inout)  :: cv

  character(len = *), intent(in)  :: avgcrd_file

!  integer, intent(in) :: first, last

  integer :: i
  character :: dummy 

  open(CRD_UNIT1, FILE = avgcrd_file)
  ! skip the first two lines of CRD files
  read(CRD_UNIT1, '(A1)') dummy
  read(CRD_UNIT1, '(A1)') dummy

!  read(CRD_UNIT1, '(7F11.5)') (cv%avgcrd(i),i=first*3-2,(last-first+1)*3)
  read(CRD_UNIT1, '(7F11.5)') (cv%avgcrd(i), i=1, 3*cv%i(3))


  close(CRD_UNIT1)
end subroutine read_avgcrd

subroutine read_refcrd(cv, refcrd)

  use nfe_colvar_type
  use nfe_constants

  implicit none

  type(colvar_t), intent(inout)  :: cv

  character(len = *), intent(in)  :: refcrd

!  integer, intent(in) :: first, last

  integer :: i
  character :: dummy

  open(REF_UNIT1, FILE = refcrd)
  ! skip the first two lines of the CRD files 
  read(REF_UNIT1, '(A1)') dummy
  read(REF_UNIT1, '(A1)') dummy

!  read(REF_UNIT1, '(6F12.7)') (cv%r(i), i=first*3-2, (last-first+1)*3)
  read(REF_UNIT1, '(6F12.7)') (cv%r(i), i=1, 3*cv%i(1))
  
  close(REF_UNIT1)

end subroutine read_refcrd

subroutine read_index(cv, index_file)

  use nfe_colvar_type
  use nfe_constants
  use nfe_utils

  implicit none

  type(colvar_t), intent(inout)  :: cv

!  character(len = *), intent(in)  :: crd_file

  character(len = *), intent(in)  :: index_file

!  integer, intent(in) :: first, last

!  integer, intent (in) :: nsolut

  integer :: i, tempi, nref, npca

!  character :: dummy

!  open(REF_UNIT1, FILE = crd_file)
  open(IDX_UNIT1, FILE = index_file)
  ! skip the first two lines of the CRD files 
!  read(REF_UNIT1, '(A1)') dummy
!  read(REF_UNIT1, '(A1)') dummy

!  read(REF_UNIT1, '(6F12.7)') (cv%r(i), i=first*3-2, (last-first+1)*3)

   read(IDX_UNIT1, *) (tempi, cv%state_ref(i), cv%state_pca(i), i=1, cv%i(1))
   
   nref = 0 
   npca = 0
   do i = 1, cv%i(1)
    
      if(cv%state_ref(i) == 1) then
         nref = nref + 1
      end if
      
      if(cv%state_pca(i) == 1) then
         npca = npca + 1
         cv%ipca_to_i(npca) = i
      end if
      
   end do
  
  ! write(*,*) "nref=", nref, " npca=", npca, " cv%i(1)=", cv%i(1), " cv%i(2)=", cv%i(2), " cv%i(3)=", cv%i(3) 
  
   nfe_assert( nref == cv%i(2) )
   nfe_assert( npca == cv%i(3) )
   nfe_assert( size(cv%ipca_to_i) == cv%i(3) )
  
  close(IDX_UNIT1)

end subroutine read_index

end module nfe_read_pca



