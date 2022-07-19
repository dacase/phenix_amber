
! epilogue: 12-10 LJ terms.


!  --- this code allows 10-12 terms; in many (most?) (all?) cases, the
!       only "nominal" 10-12 terms are on waters, where the asol and bsol
!       parameters are always zero; hence we can skip this part.

do im_new = 1,icount
   j = cache_bckptr(im_new)

   df =   cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)

   dfx = delx*df
   dfy = dely*df
   dfz = delz*df
   dumx = dumx + dfx
   dumy = dumy + dfy
   dumz = dumz + dfz
   force(1,j) = force(1,j) + dfx
   force(2,j) = force(2,j) + dfy
   force(3,j) = force(3,j) + dfz
end do  !  im_new = 1,icount
