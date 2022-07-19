      subroutine dsarpack(n_dim,n_eig_in,n_eig_out,ncv_in,itr_in,  &
                          eigval_tol,eigvals,eigvecs,spectrum,  &
                          need_eigvecs,ierr,debug_arpack,  &
                          v,workl,workd,d,resid,ax,select,  &
                          xyz,grad,return_flag,label)
!
      implicit none
!
!     %-----------------%
!     | Dummy Arguments |
!     %-----------------%
!
      integer n_dim,n_eig_in,n_eig_out,ncv_in,itr_in,spectrum, &
              need_eigvecs,ierr,debug_arpack,return_flag,label
      Double precision eigval_tol
      Double precision eigvals(n_eig_in),eigvecs(n_dim * n_eig_in)
      Double precision v(n_dim,ncv_in), &
                       workl(ncv_in*(ncv_in+8)),workd(3*n_dim), &
                       d(ncv_in,2),resid(n_dim),ax(n_dim), &
                       xyz(n_dim),grad(n_dim)
      logical select(ncv_in)
!
      save
!
!     %---------------%
!     | Include Files |
!     %---------------%
!
!
!     %---------------------------------%
!     | See debug.doc for documentation |
!     %---------------------------------%
      integer  logfil, ndigit, mgetv0,  &
               msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
               mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
               mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ logfil, ndigit, mgetv0, &
               msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
               mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
               mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     This code shows how to use ARPACK to find a few eigenvalues 
!     (lambda) and corresponding eigenvectors (x) for the standard 
!     eigenvalue problem:
!          
!                        A*x = lambda*x
! 
!     where A is an n by n real symmetric matrix.
!
!     The main points illustrated here are 
!
!        1) How to declare sufficient memory to find NEV 
!           eigenvalues of largest magnitude.  Other options
!           are available.
!
!        2) Illustration of the reverse communication interface 
!           needed to utilize the top level ARPACK routine DSAUPD 
!           that computes the quantities needed to construct
!           the desired eigenvalues and eigenvectors(if requested).
!
!        3) How to extract the desired eigenvalues and eigenvectors
!           using the ARPACK routine DSEUPD.
!
!     The only thing that must be supplied in order to use this
!     routine on your problem is to change the array dimensions 
!     appropriately, to specify WHICH eigenvalues you want to compute 
!     and to supply a matrix-vector product
!
!                         w <-  Av
!
!     in place of the call to AV( ) below.
!
!     Once usage of this routine is understood, you may wish to explore
!     the other available options to improve convergence, to solve generalized
!     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
!     This codes implements  
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is derived from the central difference discretization
!         of the 2-dimensional Laplacian on the unit square with
!         zero Dirichlet boundary condition.
!     ... OP = A  and  B = I.
!     ... Assume "call av (n,x,y)" computes y = A*x
!     ... Use mode 1 of DSAUPD.
!
!\BeginLib
!
!\Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: %Z%
! FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R%
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!     %-------------------------------------------------------%
!     | Storage Declarations:                                 |
!     |                                                       |
!     | The maximum dimensions for all arrays are             |
!     | set here to accommodate a problem size of             |
!     | N .le. MAXN                                           |
!     |                                                       |
!     | NEV is the number of eigenvalues requested.           |
!     |     See specifications for ARPACK usage below.        |
!     |                                                       |
!     | NCV is the largest number of basis vectors that will  |
!     |     be used in the Implicitly Restarted Arnoldi       |
!     |     Process.  Work per major iteration is             |
!     |     proportional to N*NCV*NCV.                        |
!     |                                                       |
!     | You must set:                                         |
!     |                                                       |
!     | MAXN:   Maximum dimension of the A allowed. (dynamic) |
!     | MAXNEV: Maximum NEV allowed. (dynamic)                |
!     | MAXNCV: Maximum NCV allowed. (dynamic)                |
!     %-------------------------------------------------------%
!
!     %--------------------------------------%
!     | F90 Allocatable Arrays (on the heap) |
!     %--------------------------------------%
!
!     Double precision,allocatable,save :: v(:,:)
!     integer,save :: v_row_allocated = 0, v_col_allocated = 0
!
!     %----------------------------------------------%
!     | Originally, as F77 parameters, the following |
!     | integers were used to dimension work arrays. |
!     | They are replaced by dummy arguments used to |
!     | dimension the work arrays as F90 automatic   |
!     | arrays, but the integers are still used for  |
!     | passing the dimensions to lower level ARPACK |
!     | routines dsaupd, dseupd and dmout.           |
!     %----------------------------------------------%
!
      integer          maxn, maxnev, maxncv, ldv
!
!     %-------------------------------------------%
!     | Local F90 Automatic Arrays (on the stack) |
!     %-------------------------------------------%
!
      Double precision cg_dstat(4)
      integer          iparam(11), ipntr(11), cg_istat(4)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, &
                       i, j, nx, ishfts, maxitr, mode1, nconv
      integer          L12, L18, ARPACK_ERROR, status_flag
      data             L12, L18, ARPACK_ERROR /1, 2, -2/
!     integer          v_row_needed, v_col_needed
      logical          rvec
      Double precision tol, sigma
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision zero
      parameter        (zero = 0.0D+0)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision dnrm2
      external         dnrm2, daxpy, hessvec
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic        abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if ( label.eq.0 ) go to 1
      go to (12,18) label
  1   continue
!
!     %------------------------------------------------%
!     | Values used to calculate work array dimensions |
!     %------------------------------------------------%
!
      maxn = n_dim
      maxnev = n_eig_in
      maxncv = ncv_in
      ldv = maxn
!
!     %---------------------------------------------------%
!     | The include debug.h statement above and           |
!     | assignments here initiate trace output from the   |
!     | internal actions of ARPACK.  See debug.doc in the |
!     | DOCUMENTS directory for usage.  Initially, the    |
!     | most useful information will be a breakdown of    |
!     | time spent in the various stages of computation   |
!     | given by setting msaupd = 1.                      |
!     %---------------------------------------------------%
!
      ndigit = -5
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      if ( debug_arpack.eq.1 ) then
        msaupd = 1
      else
        msaupd = 0
      endif
      msaup2 = 0
      mseigt = 0
      mseupd = 0
!     
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      n = n_dim
!
!     %----------------------------------------------%
!     |                                              | 
!     | Specifications for ARPACK usage are set      | 
!     | below:                                       |
!     |                                              |
!     |    1) NEV = N_EIG_IN  asks for N_EIG_IN      |  
!     |       eigenvalues to be computed.            | 
!     |                                              |
!     |    2) NCV = NCV_IN sets the length of the    |
!     |       Arnoldi factorization                  |
!     |                                              |
!     |    3) This is a standard problem             |
!     |         (indicated by bmat  = 'I')           |
!     |                                              |
!     |    4) Ask for the NEV eigenvalues of         |
!     |       smallest magnitude                     |
!     |         (indicated by which = 'SM')          |
!     |       See documentation in SSAUPD for the    |
!     |       other options SA, LA, LM, BE.          | 
!     |                                              |
!     | Note: NEV and NCV must satisfy the following |
!     | conditions:                                  |
!     |              NEV <= MAXNEV                   |
!     |          NEV + 1 <= NCV <= MAXNCV            |
!     %----------------------------------------------%
!
      nev   = n_eig_in
      ncv   = ncv_in 
      bmat  = 'I'
      if ( spectrum .eq. 1 ) then
         which = 'SM'
      else if ( spectrum .eq. 2 ) then
         which = 'SA'
      else if ( spectrum .eq. 3 ) then
         which = 'LM'
      else if ( spectrum .eq. 4 ) then
         which = 'LA'
      else if ( spectrum .eq. 5 ) then
         which = 'BE'
      else
          print *, ' ERROR with _SSIMP: Spectrum .NE. (SM|SA|LA|LM|BE)'
         go to 9000
      end if
!
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DSAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from DSAUPD. (See usage below.)                |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DSAUPD.                                | 
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     | 
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID.)  | 
!     |                                                     |
!     | The work array WORKL is used in DSAUPD as           | 
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = eigval_tol 
      info = 0
      ido = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting PARAM(1) = 1).              |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DSAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = itr_in 
      mode1 = 1
!
      iparam(1) = ishfts
!                
      iparam(3) = maxitr
!                  
      iparam(7) = mode1
!
!     %------------------------------------------------%
!     | M A I N   L O O P (Reverse communication loop) |
!     %------------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                       ncv, v, ldv, iparam, ipntr, workd, workl,&
                       lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            status_flag = 0
 11         continue
               call hessvec ( n, workd(ipntr(1)), workd(ipntr(2)), &
                              xyz, grad, return_flag, status_flag )
               if ( status_flag.eq.0 ) go to 13
               if ( status_flag.lt.0 ) go to 9000
               label = L12
               return
 12         go to 11
 13         continue
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if 
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
         go to 9000
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        |                                           |
!        | The routine DSEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!           
         if ( need_eigvecs .eq. 1 ) then
            rvec = .true.
         else
            rvec = .false.
         end if
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,  &
              bmat, n, which, nev, tol, resid, ncv, v, ldv, &
              iparam, ipntr, workd, workl, lworkl, ierr )
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NCONV (=IPARAM(5)) columns of the  |
!        | two dimensional array V if requested.        |
!        | Otherwise, an orthogonal basis for the       |
!        | invariant subspace corresponding to the      |
!        | eigenvalues in D is returned in V.           |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DSEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' '
            go to 9000
!
         else if ( debug_arpack.eq.1 ) then
!
            nconv =  iparam(5)
            n_eig_out = nconv
            if ( nconv .le. 0 ) then
               print *, ' '
               print *, ' ARPACK: Not a single mode converged.'
               print *, ' '
               go to 9000
            endif
!
!           %--------------------------------------------%
!           | "UnDO" DO 20 j=1,nconv loop, because it is |
!           | illegal to jump in and out from a DO loop. |
!           %--------------------------------------------%
!
            j = 1
 16         continue
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               status_flag = 0
 17            continue
                  call hessvec ( n, v(1,j), ax, xyz, grad,  &
                                 return_flag, status_flag )
                  if ( status_flag.eq.0 ) go to 19
                  if ( status_flag.lt.0 ) go to 9000
                  label = L18
                  return
 18            go to 17
 19            continue
!
               call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               d(j,2) = dnrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
!
               j = j + 1
               if ( j .gt. nconv ) go to 20
!
               go to 16
!
 20         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call dmout(6, nconv, 2, d, maxncv, -6,  &
                 'Ritz values and relative residuals')
!
!           %-------------------------------------------%
!           | Print additional convergence information. |
!           %-------------------------------------------%
!
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit', &
                        ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
!
            print *, ' '
            print *, ' _SSIMP '
            print *, ' ====== '
            print *, ' '
            print *, ' Size of the matrix is ', n
            print *, ' The number of Ritz values requested is ', nev
            print *, ' The number of Arnoldi vectors generated',  &
                     ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
            print *, ' The number of converged Ritz values is ',  &
                       nconv
            print *, ' The number of Implicit Arnoldi update',  &
                     ' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol
            print *, ' '
         end if
!
!        %----------------------------%
!        | Return eigvals and eigvecs |
!        %----------------------------%
!
         nconv =  iparam(5)
         n_eig_out = nconv
         if ( nconv .le. 0 ) then
            print *, ' '
            print *, ' ARPACK: Not a single mode converged.'
            print *, ' '
            go to 9000
         endif
!
         do 40 j=1, nconv
             eigvals(j) = d(j,1)
!
             do 30 i=1, n
                eigvecs((j-1)*n+i) = v(i,j)
 30          continue
 40      continue
!
      end if
!
!     %--------------------------------%
!     | Done with subroutine dsarpack. |
!     %--------------------------------%
!
      label = 0
      return
!
 9000 continue !!! Error
!
      if( status_flag.eq.0 ) status_flag = ARPACK_ERROR
!
      label = status_flag
      return
!
      end
! 
! ------------------------------------------------------------------
