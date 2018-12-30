$FIXEDFORMLINESIZE:132
c #    #   #####     #    #          #     #####     #    ######   ####
c #    #     #       #    #          #       #       #    #       #
c #    #     #       #    #          #       #       #    #####    ####
c #    #     #       #    #          #       #       #    #            #
c #    #     #       #    #          #       #       #    #       #    #
c  ####      #       #    ######     #       #       #    ######   ####
c
c ---------------------------------------------------------------------
      subroutine scopy (nt, arr, brr)
cc      subroutine scopy (arr, brr)
c ---------------------------------------------------------------------
c
c copies array arr to brr,for 1st nt words of arr
c
      include 'implicit.h'
c
c Arguments
c
      integer nt     
      real arr(1),    ! input
     >     brr(1)     ! output
c
c Local variables
c
      integer ia

cc	nt = 1   !Yuan
c
      do 100 ia = 1, nt
        brr(ia) = arr(ia)
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine const (arr, nar, value)
c ---------------------------------------------------------------------
c
c sets all elements of real vector arr to value
c
      include 'implicit.h'
c
c Arguments
c
      integer nar
c     
      real value
      real arr(nar)
c
c Local variables
c
      integer j
c
      do 100 j = 1, nar
        arr(j) = value
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine endrun
c ---------------------------------------------------------------------
c
c stops gracefully
c
      stop
      end
c
c 
c ---------------------------------------------------------------------
      real function cvmgt (x,y,l)
c ---------------------------------------------------------------------
c
c chooses between two things.  Used in canopy.f
c
      include 'implicit.h'
c
      logical l
      real x, y
c
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
c lenchr - find index of last non-blank, non-null
c ---------------------------------------------------------------------
c
      integer function lenchr (ch)
c
c returns position of last non-blank,null character in ch,
c or 1 if ch is all blanks
c
      include 'implicit.h'
c
c Arguments
c
      character*(*) ch
c
      integer i
c
      do i = len(ch), 1, -1
        if (ch(i:i).ne.' '.and.ch(i:i).ne.char(0)) then
           lenchr = i
           return
        endif
      enddo
c
      lenchr = 1
c
      return
      end
c
