! <compile=optimized>

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Hypernetted chain (HNC) equation closure class for 3D-RISM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_hnc_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the HNC type
    type rism3d_hnc
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
    end type rism3d_hnc

    public rism3d_hnc_new, rism3d_hnc_destroy!, rism3d_hnc_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the HNC closure
!!!IN:
!!!   this : HNC object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_new(this,pot)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      this%pot => pot
      this%grid => this%pot%grid
    end subroutine rism3d_hnc_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the HNC closure
!!!IN:
!!!   this : the HNC closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_guv(this,guv, huv, cuv)
      use constants_rism, only: omp_num_threads
      implicit none
      type(rism3d_hnc), intent(in) :: this
      _REAL_, intent(out) :: guv(:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer :: iv, ir, ix, iy, iz, ig
      _REAL_ :: exponent

!$omp parallel do private(iv,ix,iy,iz,ig,exponent)  &
!$omp&        num_threads(omp_num_threads)
      do iz = 1, this%grid%localDimsR(3)
         do iy = 1, this%grid%localDimsR(2)
            do ix = 1, this%grid%localDimsR(1)
#if defined(MPI)
               ig = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                  (iz-1)*(this%grid%localDimsR(1)+2)*this%grid%localDimsR(2)
#else
               ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                    (iz - 1) * this%grid%localDimsR(1) * this%grid%localDimsR(2)
#endif /*defined(MPI)*/
               do iv = 1,this%pot%solvent%numAtomTypes
                  exponent = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                  guv(ig,iv) = exp(exponent)
               end do
            end do
         end do
      end do
!$omp end parallel do
    end subroutine rism3d_hnc_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_hnc_excessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
    implicit none
    type(rism3d_hnc), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%pot%solvent%numAtomTypes)
    _REAL_ :: tuv, hk0
    integer :: ix, iy, iz, iv, ig, igk
    excessChemicalPotential = 0.d0
!$omp parallel do private (iv,iz,iy,ix,ig,igk,tuv,hk0) &
!$omp&   num_threads(this%pot%solvent%numAtomTypes)
    do iv=1,this%pot%solvent%numAtomTypes
       hk0 = 1.d0 + this%pot%huvk0(1,iv)/2.d0
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                     (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                   (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) + &
                    (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                excessChemicalPotential(iv) = excessChemicalPotential(iv) + &
                     0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)*hk0
             end do
          end do
       end do
       excessChemicalPotential(iv) =  this%pot%solvent%density(iv)&
            *excessChemicalPotential(iv)*this%grid%voxelVolume
    enddo
!$omp end parallel do
  end function rism3d_hnc_excessChemicalPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the HNC closure
!!!IN:
!!!   this : HNC object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_destroy(this)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_hnc_destroy
  end module rism3d_hnc_c
