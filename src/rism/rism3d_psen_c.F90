!<compile=optimized>

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Partial series expansion of order-n (PSE-n) closure class for 3D-RISM.
!!!Stefan M. Kast and Thomas Kloss. J. Chem. Phys. 129, 236101 (2008)
!!!
!!!Interpolates between the Kovalenko-Hirata (KH) and hypernetted chain equation
!!!closures (HNC).  Order-1 gives KH and order-infinity gives HNC.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_psen_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the PSEN type
    type rism3d_psen
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
       !order    : order of the series expansion.  >= 1
       !order1   : order+1
       integer :: order=0, order1
       !order1fac: (order+1)!
       _REAL_ :: order1fac
    end type rism3d_psen

    public rism3d_psen_new, rism3d_psen_destroy!, rism3d_psen_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_new(this,pot, order)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      integer, intent(in) :: order
      integer :: i, orderfac
      this%pot => pot
      this%grid => this%pot%grid
      if(order <1)then
         call rism_report_error('(a,i4)',"PSE-n closure: order must be > 0:",order)
      end if
      this%order = order
      this%order1 = order+1
      this%order1fac = 1
      orderfac=1

      !calculate the factorial of (order+1) for future reference.  Check that we 
      !are not losing precision or overflowing the precision
      do i=2,this%order1
         orderfac = this%order1fac
         this%order1fac = this%order1fac*i
         !keep 16 significant digits
         if( abs((this%order1fac/i)/orderfac -1d0) >1d-16 )then
            call rism_report_error('(a,i4)',&
                 "PSE-n: numerical precision exceeded.  Use a smaller order. Try ",i-1)
         end if
      end do
    end subroutine rism3d_psen_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the PSEN closure
!!!IN:
!!!   this : the PSEN closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_psen_guv(this,guv, huv, cuv)
    use constants_rism, only: omp_num_threads
    implicit none
    type(rism3d_psen), intent(in) :: this
    _REAL_, intent(out) :: guv(:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer :: i, iv, ir, ix, iy, iz, ig
    _REAL_ :: tuv, orderfac

!$omp parallel do private(iv,ix,iy,iz,ig,i,orderfac,tuv)  &
!$omp&        num_threads(omp_num_threads)
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
#ifdef MPI
             ig = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                 (iz-1)*(this%grid%localDimsR(1)+2)*this%grid%localDimsR(2)
#else
             ig = ix + (iy-1)*this%grid%localDimsR(1) + &
                 (iz-1)*this%grid%localDimsR(1)*this%grid%localDimsR(2)
#endif
             do iv = 1,this%pot%solvent%numAtomTypes
                tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                if(tuv >= 0d0)then
                   guv(ig,iv) = 1d0
                   orderfac = 1d0
                   do i=1,this%order
                      orderfac = orderfac*i
                      guv(ig,iv) = guv(ig,iv) + (tuv**i)/orderfac
                   end do
                else
                   guv(ig,iv) = exp(tuv)
                endif
             end do
          end do
       end do
    end do
!$omp end parallel do
  end subroutine rism3d_psen_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_psen_excessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
    implicit none
    type(rism3d_psen), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%pot%solvent%numAtomTypes)
    _REAL_ :: tuv, tsuv, hk0
    integer :: ix, iy, iz, iv, igk
    excessChemicalPotential = 0.d0
!$omp parallel do private (iv,iz,iy,ix,igk,tuv,hk0,tsuv) &
!$omp&   num_threads(this%pot%solvent%numAtomTypes)
    do iv=1,this%pot%solvent%numAtomTypes
       hk0 = 1.d0 + this%pot%huvk0(1,iv)/2.d0
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                   (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
                igk = ix + (iy-1)*this%grid%localDimsR(1) + &
                   (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif
                tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                excessChemicalPotential(iv) = excessChemicalPotential(iv) &
                  + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)*hk0
                if (huv(igk,iv) > 0d0)  then
                   tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
                   excessChemicalPotential(iv) = excessChemicalPotential(iv) &
                       - tsuv**(this%order1)/this%order1fac
                endif
             end do
          end do
       end do
       excessChemicalPotential(iv) =  this%pot%solvent%density(iv)*&
            excessChemicalPotential(iv)*this%grid%voxelVolume
    enddo
!$omp end parallel do
  end function rism3d_psen_excessChemicalPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_destroy(this)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_psen_destroy
  end module rism3d_psen_c
