$FIXEDFORMLINESIZE:132
c #    #    ##     #####  #    #
c ##  ##   #  #      #    #    #
c # ## #  #    #     #    ######
c #    #  ######     #    #    #
c #    #  #    #     #    #    #
c #    #  #    #     #    #    #
c
c
c --------------------------------------------------------------------
      real function ran2 (idum)
c --------------------------------------------------------------------
c
      include 'implicit.h'
c
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
c
      real am, eps, rnmx
c
      parameter (im1=2147483563, 
     >           im2=2147483399,
     >           am=1./im1,
     >           imm1=im1-1,
     >           ia1=40014,
     >           ia2=40692,
     >           iq1=53668,
     >           iq2=52774, 
     >           ir1=12211,
     >           ir2=3791,
     >           ntab=32,
     >           ndiv=1+imm1/ntab,
     >           eps=1.0e-7,
     >           rnmx=1.-eps)
c
      integer idum2,j,k,iv(ntab),iy
c
      save iv,iy,idum2
c
      data idum2/123456789/, iv/ntab*0/, iy/0/
c
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 10 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
 10     continue
        iy=iv(1)
      endif
c
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
c
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
c
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
c
      ran2=min(am*iy,rnmx)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine linsolve (arr, rhs, vec, mplate, nd)
c ---------------------------------------------------------------------
c
c solves multiple linear systems of equations, vectorizing
c over the number of systems. basic gaussian elimination is 
c used, with no pivoting (relies on all diagonal elements
c being and staying significantly non-zero)
c
c a template array mplate is used to detect when an operation 
c is not necessary (element already zero or would add zeros),
c assuming that every system has the same pattern of zero
c elements
c
c this template is first copied to mplatex since it 
c must be updated during the procedure in case an original-zero
c pattern location becomes non-zero
c
c the first subscript in arr, rhs, vec is over the multiple
c systems, and the others are the usual row, column subscripts
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (input-output)
c
      integer nd                  ! number of equations (supplied)
c
      integer  mplate(nd,nd)      ! pattern of zero elements of arr (supplied)
c
      real arr(nd,nd),       ! equation coefficients (supplied, overwritten)
     >     rhs(nd),          ! equation right-hand sides (supplied, overwritten) 
     >     vec(nd)           ! solution (returned)
c 
c local variables
c
      integer ndx,                ! Max number of equations
     >        j, i, id         ! loop indices
c
      parameter (ndx=9)
c
      integer mplatex(ndx,ndx)
c
      real f

cc	write(30,*)"linsolve 1"
c
      if (nd.gt.ndx) then
         write(*,900) nd, ndx
  900    format(/' *** fatal error ***'/
     >          /' number of linsolve eqns',i4,' exceeds limit',i4)
         call endrun
      endif

cc       write(30,*)"linsolve 2"
c
c copy the zero template so it can be changed below
c
      do 6 j=1,nd
        do 5 i=1,nd
          mplatex(i,j) = mplate(i,j)
    5   continue
    6 continue

cc	write(30,*)"linsolve 3"
c
c zero all array elements below the diagonal, proceeding from
c the first row to the last. note that mplatex is set non-zero
c for changed (i,j) locations, in loop 20
c
      do 10 id=1, nd-1
         do 12 i=id+1,nd
c
            if (mplatex(i,id).ne.0) then
cc               do 14 m=1,npoi
                  f = arr(i,id) / arr(id,id)
cc   14          continue
c
               do 20 j=id,nd
                  if (mplatex(id,j).ne.0) then
cc                     do 22 m=1,npoi
                        arr(i,j) = arr(i,j) - f*arr(id,j)
cc   22                continue
                     mplatex(i,j) = 1
                  endif
   20          continue
c
cc               do 30 m=1,npoi
                  rhs(i) = rhs(i) - f*rhs(id)
cc   30          continue
            endif
c
   12    continue
   10 continue

cc	write(30,*)"linsolve 4"
c
c all array elements below the diagonal are zero, so can
c immediately solve the equations in reverse order
c
      do 50 id=nd,1,-1
c
cc         call const (f, 0.0)
	   f = 0
         if (id.lt.nd) then
            do 52 j=id+1,nd
               if (mplatex(id,j).ne.0) then
cc                  do 54 m=1,npoi
                     f = f + arr(id,j)*vec(j)
cc   54             continue
               endif
   52       continue
         endif
c
cc         do 56 m=1,npoi
            vec(id) = (rhs(id) - f) / arr(id,id)
cc   56    continue
c
   50 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine tridia (ne, a, b, c, y, x, alpha, gamma)
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
c
c     purpose:
c     to compute the solution of many tridiagonal linear systems.
c
c      arguments:
c
c      ns ..... the number of systems to be solved.
c
c      nd ..... first dimension of arrays (ge ns).
c
c      ne ..... the number of unknowns in each system.
c               this must be > 2. second dimension of arrays.
c
c      a ...... the subdiagonals of the matrices are stored
c               in locations a(j,2) through a(j,ne).
c
c      b ...... the main diagonals of the matrices are stored
c               in locations b(j,1) through b(j,ne).
c
c      c ...... the super-diagonals of the matrices are stored in
c               locations c(j,1) through c(j,ne-1).
c
c      y ...... the right hand side of the equations is stored in
c               y(j,1) through y(j,ne).
c
c      x ...... the solutions of the systems are returned in
c               locations x(j,1) through x(j,ne).
c
c      alpha .. work array dimensioned alpha(nd,ne)
c
c      gamma .. work array dimensioned gamma(nd,ne)
c
c       history:  based on a streamlined version of the old ncar
c                 ulib subr trdi used in the phoenix climate
c                 model of schneider and thompson (j.g.r., 1981).
c                 revised by starley thompson to solve multiple
c                 systems and vectorize well on the cray-1.
c                 later revised to include a parameter statement
c                 to define loop limits and thus enable cray short
c                 vector loops.
c
c       algorithm:  lu decomposition followed by solution.
c                   note: this subr executes satisfactorily
c                   if the input matrix is diagonally dominant
c                   and non-singular.  the diagonal elements are
c                   used to pivot, and no tests are made to determine
c                   singularity. if a singular or numerically singular
c                   matrix is used as input a divide by zero or
c                   floating point overflow will result.
c
c       last revision date:      4 february 1988
c
c
c Arguments
c
      integer ns,     ! number of systems to be solved.
cc     >        nd,     ! first dimension of arrays (ge ns)
     >        ne      ! number of unknowns in each system. (>2)
      
      real 
     >  a(1,ne),     ! subdiagonals of matrices stored in a(j,2)...a(j,ne).
     >  b(1,ne),     ! main diagonals of matrices stored in b(j,1)...b(j,ne).
     >  c(1,ne),     ! super-diagonals of matrices stored in c(j,1)...c(j,ne-1).
     >  y(1,ne),     ! right hand side of equations stored in y(j,1)...y(j,ne).
     >  x(1,ne),     ! solutions of the systems returned in x(j,1)...x(j,ne).
     >  alpha(1,ne), ! work array 
     >  gamma(1,ne)  ! work array
c
c local variables
c
      integer nm1,    !
     >  j, i, ib      ! loop indices

c
      nm1 = ne-1
c
c obtain the lu decompositions
c
      ns = 1

      do 10 j=1,ns
         alpha(j,1) = 1./b(j,1)
         gamma(j,1) = c(j,1)*alpha(j,1)
   10 continue
      do 11 i=2,nm1
         do 12 j=1,ns
            alpha(j,i) = 1./(b(j,i)-a(j,i)*gamma(j,i-1))
            gamma(j,i) = c(j,i)*alpha(j,i)
   12    continue
   11 continue
c
c solve
c
      do 20 j=1,ns
         x(j,1) = y(j,1)*alpha(j,1)
   20 continue
      do 21 i=2,nm1
         do 22 j=1,ns
            x(j,i) = (y(j,i)-a(j,i)*x(j,i-1))*alpha(j,i)
   22    continue
   21 continue
      do 23 j=1,ns
         x(j,ne) = (y(j,ne)-a(j,ne)*x(j,nm1))/
     >             (b(j,ne)-a(j,ne)*gamma(j,nm1))
   23 continue
      do 24 i=1,nm1
         ib = ne-i
         do 25 j=1,ns
            x(j,ib) = x(j,ib)-gamma(j,ib)*x(j,ib+1)
   25    continue
   24 continue
c
      return
      end
c
	REAL FUNCTION RAN(IDUM)
!  UNIFORM PSEUDORANDOM NUMBER GENERATOR
!   A LINEAR CONGRUENTIAL GENERATOR X(N+1) = MOD( A*X(N), 2**31 - 1 )
!  CODING FOLLOWS ALGORITHM 1 OF HORMANN AND DERFLINGER, WITH
!   PERSONAL MODIFICATIONS TO AVOID OVERFLOWS
!
!  W. HORMANN & G. DERFLINGER (1993) "A PORTABLE RANDOM NUMBER GENERATOR
!   WELL SUITED FOR THE REJECTION METHOD," ACM TOMS VOL 19, PP.489-495.
!
!  ARGUMENT 
!   IDUM    INTEGER  FIRST CALL SETS SEED, IGNORED IN SUBSEQUENT CALLS
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDUM
      INTEGER, PARAMETER :: AHI = 12121          ! PART OF MULTIPLIER
      INTEGER, PARAMETER :: ALOW = 23166          ! PART OF MULTIPLIER
      INTEGER, PARAMETER :: ALOW2 = 46332         ! PART OF MULTIPLIER
      INTEGER, PARAMETER :: B15 = 32768           ! 2**15
      INTEGER, PARAMETER :: B16 = 65536           ! 2**16
      INTEGER, PARAMETER :: P = 2147483647        ! MODULUS 2**31 - 1
      INTEGER, SAVE :: XHI = 0
      INTEGER, SAVE :: XLOW = 0
      INTEGER MID1,MID2,MID,X,K
!
!                   IF NOT FIRST CALL, THEN SKIP SETTING SEED
      IF( (XHI .EQ. 0) .AND. (XLOW .EQ. 0) ) THEN
      XHI = IDUM / 65536
      XLOW = IDUM - XHI*65536
                                             END IF ! ( FIRST CALL )
!                   MULTIPLIER IS A = 397204094 = AHI*(2**15) + ALOW
      MID1 = AHI*XLOW 
      MID2 = ALOW2*XHI                                                     
!                   TEST FOR OVERFLOW
      IF( MID1-1 .GT. P-MID2 ) THEN
!                   HERE MID IS > 2**31, SO WRITE AS MID - 2**31
      MID = (MID1-1) - (P-MID2)
      K = 1
                               ELSE
!                   HERE MID IS < 2**31, SO WRITE AS POSITIVE
      MID = MID1 + MID2
      K = 0
                               END IF ! ( MID1-1 .GT. P-MID2 )
!                   NOW GET TO MAIN STEP
!                    SUBTRACT P = 2**31 - 1 IN ADVANCE
      MID2 = MID / B16
      X = (ALOW*XLOW - P) + AHI*XHI + MID2 + K*B15
      IF( X .LT. 0 ) X = X + P
      X = X + ( ( MID - MID2*B16 )*B15 - P )
      IF( X .LT. 0 ) X = X + P
!                   WILL NEED TWO PARTS OF X FOR NEXT CALL
      XHI = X / B16
      XLOW = X - XHI*B16
!                                          2**16      2**15
      RAN = ( REAL(XHI) + ( REAL(XLOW) / 65536. ) )/32768.
	 
      RETURN
      END FUNCTION RAN	