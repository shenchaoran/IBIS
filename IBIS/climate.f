
c  ####   #          #    #    #    ##     #####  ######
c #    #  #          #    ##  ##   #  #      #    #
c #       #          #    # ## #  #    #     #    #####
c #       #          #    #    #  ######     #    #
c #    #  #          #    #    #  #    #     #    #
c  ####   ######     #    #    #  #    #     #    ######
c
c ---------------------------------------------------------------------
      subroutine climanl
c ---------------------------------------------------------------------
c
c this subsroutine is only used to initialize growing degree days,
c coldest temp, and warmest temp at very beginning - provides a
c climate 'history' based on monthly mean values
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
c
c Local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >     k,lda      ! loop indices
c
      real rwork,        ! work variable (1/ndaypm)
     >     dt,           ! used for interpolation
     >     dtemp         ! interpolated temperature
c
c initialize values
c
      gdd0 = 0.0
      gdd5 = 0.0
c
cc      do 100 i = 1, npoi
c
c coldest monthly temperature (year 0) in deg c
c
        tc = min (xint(1),  xint(2),  xint(3),
     >               xint(4),  xint(5),  xint(6),
     >               xint(7),  xint(8),  xint(9),
     >               xint(10), xint(11), xint(12))
c
c warmest monthly temperature (year 0) in deg c
c
        tw = max (xint(1),  xint(2),  xint(3),
     >               xint(4),  xint(5),  xint(6),
     >               xint(7),  xint(8),  xint(9),
     >               xint(10), xint(11), xint(12))
c
        tcmin = tc + deltat
c
cc 100  continue 
c
c interpolating climatological monthly input values to daily
c
cc      do 200 i = 1, npoi 
c
        do 210 k = 1, 12
c
          rwork = 1. / float(ndaypm(k))
c
          do 220 lda = 1, ndaypm(k)
c
            if (float(lda).lt.float(ndaypm(k)+1)*0.5) then
              it1w = k - 1
              it2w = k
              dt   = (float(lda) - 0.5) * rwork + 0.5
            else
              it1w = k
              it2w = k + 1
              dt   = (float(lda) - 0.5) * rwork - 0.5
            end if
c
            if (it1w.lt. 1) it1w = 12
            if (it2w.gt.12) it2w = 1
c
            dtemp = xint(it1w) + dt * (xint(it2w) - xint(it1w))
c
c growing degree days, using deg c
c
            gdd0 = gdd0 + max(0.0, dtemp)
            gdd5 = gdd5 + max(0.0, (dtemp - 5.0))
c
 220      continue
 210    continue
c
cc 200  continue
c
c call routine to determine pft existence arrays
c
      call existence
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine climanl2
c ---------------------------------------------------------------------
c
c this subroutine updates the growing degree days, coldest temp, and
c warmest temp if monthly anomalies or daily values are used
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
c
c local variables
c
cc      integer i             ! loop indice
c
      real zweigc,          ! 30-year e-folding time-avarage
     >     zweigw,          ! 30-year e-folding time-avarage
     >     rworkc,          ! 30-year e-folding time-avarage
     >     rworkw 
c
c calculate a 30-year e-folding time-avarage
c
      zweigc = exp(-1./30.)
      zweigw = exp(-1./30.)
c
      rworkc = 1. - zweigc
      rworkw = 1. - zweigw
c
c update critical climatic parameters with running average
c
cc      do 100 i = 1, npoi
c
        tc = zweigc * tc + rworkc * tcthis
        tw = zweigw * tw + rworkw * twthis
c
        tcmin = tc + deltat
c
        gdd0 = zweigc * gdd0 +
     >            rworkc * gdd0this
c
        gdd5 = zweigc * gdd5 +
     >            rworkc * gdd5this
c
cc 100  continue
c
      call existence
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine existence
c ---------------------------------------------------------------------
c
c this routine determines which plant functional types (pft's) are allowed
c to exist in each gridcell, based on a simple set of climatic criteria
c
c the logic here is based on the biome3 model of haxeltine and prentice
c
c plant functional types:
c
c 1)  tropical broadleaf evergreen trees
c 2)  tropical broadleaf drought-deciduous trees
c 3)  warm-temperate broadleaf evergreen trees
c 4)  temperate conifer evergreen trees
c 5)  temperate broadleaf cold-deciduous trees
c 6)  boreal conifer evergreen trees
c 7)  boreal broadleaf cold-deciduous trees
c 8)  boreal conifer cold-deciduous trees
c 9)  evergreen shrubs
c 10) deciduous shrubs
c 11) warm (c4) grasses
c 12) cool (c3) grasses
c
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'compft.h'
c
c Local variables
c
cc      integer i      ! loop indice
c
c ---------------------------------------------------------------------
c
cc      do 100 i = 1, npoi
c
c determine which plant types can exist in a given gridcell
c
        exist(1)  = 0.
        exist(2)  = 0.
        exist(3)  = 0.
        exist(4)  = 0.
        exist(5)  = 0.
        exist(6)  = 0.
        exist(7)  = 0.
        exist(8)  = 0.
        exist(9)  = 0.
        exist(10) = 0.
        exist(11) = 0.
        exist(12) = 0.
c
c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
*        if (tcmin.gt.0.0)           exist(1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
*        if (tcmin.gt.0.0)           exist(2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
*        if ((tcmin.lt.0.0).and.
*     >      (tcmin.gt.-10.0))       exist(3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin.lt.0.0).and.
*     >      (tcmin.gt.-45.0).and.
*     >      (gdd5.gt.1200.0))       exist(4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin.lt.0.0).and.
*     >      (tcmin.gt.-45.0).and.
*     >      (gdd5.gt.1200.0))       exist(5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*        if (((tcmin.lt.-45.0).or.(gdd5.lt.1200.0)).and.
*     >       (tcmin.gt.-57.5).and.
*     >       (gdd5.gt.350.0))       exist(6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*        if (((tcmin.lt.-45.0).or.(gdd5.lt.1200.0)).and.
*     >       (tcmin.gt.-57.5).and.
*     >       (gdd5.gt.350.0))       exist(7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
*        if (((tcmin.lt.-45.0).or.(gdd5.lt.1200.0)).and.
*     >       (gdd5.gt.350.0))       exist(8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0.gt.100.0)          exist(9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0.gt.100.0)          exist(10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
*        if ((tw.gt.22.0).and.
*     >      (gdd0.gt.100.0))        exist(11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
*        if (gdd0.gt.100.0)          exist(12) = 1.0
c
c
**** DTP 2001/06/07: Modified version of above code reads in PFT
*    existence criteria from external parameter file "params.veg"
*    These are copied here for reference.... 
*------------------------------------------------------------------
*  TminL    TminU    Twarm    GDD    PFT
*------------------------------------------------------------------
*    0.0   9999.0   9999.0   9999  !   1
*    0.0   9999.0   9999.0   9999  !   2
*  -10.0      0.0   9999.0   9999  !   3
*  -45.0      0.0   9999.0   1200  !   4
*  -45.0      0.0   9999.0   1200  !   5
*  -57.5    -45.0   9999.0    350  !   6
*  -57.5    -45.0   9999.0    350  !   7
* 9999.0    -45.0   9999.0    350  !   8
* 9999.0   9999.0   9999.0    100  !   9
* 9999.0   9999.0   9999.0    100  !  10
* 9999.0   9999.0     22.0    100  !  11
* 9999.0   9999.0   9999.0    100  !  12
*------------------------------------------------------------------

c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
        if (tcmin.gt.TminL(1))      exist(1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
        if (tcmin.gt.TminL(2))      exist(2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
        if ((tcmin.lt.TminU(3)).and.
     >      (tcmin.gt.TminL(3)))    exist(3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin.lt.TminU(4)).and.
     >      (tcmin.gt.TminL(4)).and.
     >      (gdd5.gt.GDD(4)))       exist(4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin.lt.TminU(5)).and.
     >      (tcmin.gt.TminL(5)).and.
     >      (gdd5.gt.GDD(5)))       exist(5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin.lt.TminU(6)).or.
     >      (gdd5.lt.GDD(4))).and.
     >      (tcmin.gt.TminL(6)).and.
     >      (gdd5.gt.GDD(6)))       exist(6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin.lt.TminU(7)).or.
     >      (gdd5.lt.GDD(5))).and.
     >      (tcmin.gt.TminL(7)).and.
     >      (gdd5.gt.GDD(7)))       exist(7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
        if (((tcmin.lt.TminU(8)).or.
     >      (gdd5.lt.TminL(4))).and.
     >      (gdd5.gt.GDD(8)))       exist(8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0.gt.GDD(9))         exist(9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0.gt.GDD(10))        exist(10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
        if ((tw.gt.Twarm(11)).and.
     >      (gdd0.gt.GDD(11)))      exist(11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
        if (gdd0.gt.GDD(12))        exist(12) = 1.0

cc 100  continue
c
      return
      end