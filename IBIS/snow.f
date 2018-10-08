$FIXEDFORMLINESIZE:132
c  ####   #    #   ####   #    #
c #       ##   #  #    #  #    #
c  ####   # #  #  #    #  #    #
c      #  #  # #  #    #  # ## #
c #    #  #   ##  #    #  ##  ##
c  ####   #    #   ####   #    #
c
c ---------------------------------------------------------------------
      subroutine snow
c ---------------------------------------------------------------------
c
c steps snow model through one timestep
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables:
c
      integer k,             ! loop indices
     >     npn                  ! index indsno, npcounter for pts with snow

      real rwork, rwork2,       ! working vaiable
     >     finew,               ! storing variable for fi
     >     zhh,                 ! 0.5*hsnomin
     >     zdh                  ! max height of snow above hsnomin (?)     

      integer indsno      ! index of points with snow in current 1d strip
c
      real hinit(nsnolay),      ! initial layer thicknesses when snow first forms
     >     hsnoruf,       ! heigth of snow forced to cover lower canopy (?)
     >     fiold,         ! old fi at start of this timestep
     >     fhtop,         ! heat flux into upper snow surface
     >     sflo(nsnolay+2),! heat flux across snow and buried-lower-veg layer bdries
     >     zmelt,         ! liquid mass flux increments to soil, at temperature 
     >                          ! tmelt, due to processes occuring during this step
     >     zheat,         ! heat flux to soil, due to processes occuring this step
     >     dfi,           ! change in fi
     >     xl,            ! lower veg density
     >     xh,            ! temporary arrays
     >     xm,            ! "
     >     ht,            ! "
     >     x1,            ! "
     >     x2,            ! "
     >     x3             ! "
c
cc      do 10 i = 1, npoi
        hsnoruf =  min (0.70, max (hsnomin+.05, fl*ztop(1)))
        xl = fl * 2.0 * (lai(1) + sai(1))
        x1 = tlsub
cc 10   continue
c
      hinit(1) = hsnotop
c
      do 15 k = 2, nsnolay
        hinit(k) = (hsnomin - hsnotop) / (nsnolay-1)
   15 continue
c
cc      do 20 i = 1, npoi
        fiold = fi
cc   20 continue
c
c zero out arrays
c
      call const (sflo, (nsnolay+2), 0.0)
cc      call const (zmelt, 0.0)
cc      call const (zheat, 0.0)
     
	 zmelt = 0
	 zheat = 0

c
c set up index indsno, npn for pts with snow - indsno is used
c only by vadapt - elsewhere below, just test on npn > 0
c
      npn = 0                                        
c
cc      do 30 i = 1, npoi 
        if (fi.gt.0.) then
          npn = npn + 1
          indsno = 1
        endif 
cc   30 continue
c
c set surface heat flux fhtop and increment top layer thickness
c due to rainfall, snowfall and sublimation on existing snow
c
      if (npn.gt.0) then
c
        rwork = dtime / rhos
c
cc        do 40 i = 1, npoi
c
          fhtop = heati +
     >               rainl * (ch2o * (trainl - tmelt)     + hfus +
     >                           cice * (tmelt     - tsno(1)))       +
     >               snowl *  cice * (tsnowl - tsno(1))
c
          if (fi.gt.0.) hsno(1) = hsno(1) +
     >                     (rainl + snowl - fvapi) * rwork
c
cc   40   continue
c
      endif
c
c step temperatures due to heat conduction, including buried
c lower-veg temperature tlsub
c
      if (npn.gt.0) then
c
cc        call scopy (1,tlsub, x1)
	   x1 = tlsub
        call snowheat (tlsub, fhtop, sflo, xl, chl)
c
      endif
c
c put snowfall from 1-fi snow-free area onto side of existing
c snow, or create new snow if current fi = 0. also reset index.
c (assumes total depth of newly created snow = hsnomin.)
c (fi will not become gt 1 here if one timestep's snowfall
c <= hsnomin, but protect against this anyway.)
c
c if no adjacent snowfall or fi = 1, dfi = 0, so no effect
c
cc      call const (ht, 0.0)

	ht = 0

      do 190 k=1,nsnolay
cc        do 192 i=1,npoi
          ht = ht + hsno(k)
cc  192   continue
  190 continue
c
cc      do 195 i=1,npoi
        if (ht.eq.0.) ht = hsnomin
cc  195 continue
c
      rwork = dtime / rhos
cc      do 200 i=1,npoi
        dfi = (1.-fi)*rwork*snowg / ht
        dfi = min (dfi, 1.-fi)
cc  200 continue
c
      do 210 k=1,nsnolay
cc        do 212 i=1,npoi
          if (fi+dfi.gt.0.)
     >      tsno(k) = (tsno(k)*fi + tsnowg*dfi) / (fi+dfi)
c
c set initial thicknesses for newly created snow
c
          if (fi.eq.0. .and. dfi.gt.0.) hsno(k) = hinit(k)
cc  212   continue
  210 continue
c
      npn = 0
cc      do 220 i=1,npoi
        fi = fi + dfi
        if (fi.gt.0.) then 
          npn = npn + 1
          indsno = 1
        endif
cc  220 continue
c
c melt from any layer (due to implicit heat conduction, any
c layer can exceed tmelt, not just the top layer), and reduce
c thicknesses (even to zero, and give extra heat to soil)
c
c ok to do it for non-snow points, for which xh = xm = 0
c
      if (npn.gt.0) then
c
        rwork = 1. / rhos
        do 300 k=1,nsnolay
cc          do 302 i=1,npoi
            xh = rhos*hsno(k)*cice * max(tsno(k)-tmelt, 0.)
            xm = min (rhos*hsno(k), xh/hfus)
            hsno(k) = hsno(k) - xm*rwork
            tsno(k) = min (tsno(k),tmelt)
            zmelt = zmelt + fi*xm
            zheat = zheat + fi*(xh-hfus*xm)
cc  302     continue
  300   continue
c
c adjust fi and thicknesses for coverage-vs-volume relation
c ie, total thickness = hsnomin for fi < fimax, and fi <= fimax.
c (ok to do it for no-snow points, for which ht=fi=finew=0.)
c
cc        call const (ht,  0.0)

	  ht = 0

        do 400 k=1,nsnolay
cc          do 402 i=1,npoi
            ht = ht + hsno(k)
cc  402     continue
  400   continue
c
c linear variation  for 0 < fi < 1
c
        zhh = 0.5*hsnomin
cc        do 404 i=1,npoi
          zdh = hsnoruf-hsnomin
          finew = ( -zhh + sqrt(zhh**2 + zdh*fi*ht) ) / zdh

          finew = max (0., min (fimax, finew))
          x1 =  fi / max (finew, epsilon)
          fi =  finew
cc  404   continue
c
        do 406 k=1,nsnolay
cc          do 408 i=1,npoi
            hsno(k) = hsno(k) * x1
cc  408     continue
  406   continue
c
      endif
c
c re-adapt snow thickness profile, so top thickness = hsnotop
c and other thicknesses are equal
c
c adjust temperature to conserve sensible heat
c
      call vadapt (hsno, tsno, hsnotop, indsno, npn, nsnolay)
c
c if fi is below fimin, melt all snow and adjust soil fluxes
c
      if (npn.gt.0) then
cc        call scopy (1,fi, x1)
	  x1 = fi

        do 500 k=1,nsnolay
cc          do 502 i=1,npoi
            if (x1.lt.fimin) then
              xm = x1 * rhos * hsno(k)
              zmelt = zmelt + xm
              zheat = zheat - xm*(cice*(tmelt-tsno(k))+hfus)
              hsno(k) = 0.
              tsno(k) = tmelt
              fi = 0.
            endif
cc  502     continue
  500   continue
      endif
c
c adjust buried lower veg for fi changes. if fi has increased,
c incorporate newly buried intercepted h2o into bottom-layer 
c snow, giving associated heat increment to soil, and mix the
c specific heat of newly buried veg (at tl) into tlsub. if fi
c has decreased, change temp of newly exhumed veg to tl, giving
c assoc heat increment to soil, and smear out intercepted h2o
c
      if (npn.gt.0) then
cc        do 600 i=1,npoi
          dfi = fi - fiold
c
          if (dfi.gt.0.) then
c
c factor of xl*chl has been divided out of next line
c
            tlsub= (tlsub*fiold+ tl*dfi) / fi
            zheat = zheat + dfi*xl
     >               * ( wliql * (ch2o*(tl-tmelt) + hfus
     >                              +cice*(tmelt-tsno(nsnolay)))
     >                   + wsnol *  cice*(tl-tsno(nsnolay)) )
c
            hsno(nsnolay) = hsno(nsnolay)
     >                      + dfi*xl*(wliql+wsnol)
     >                        / (rhos*fi)
          endif
c
          if (dfi.lt.0.) then
            zheat = zheat - dfi*xl*chl*(tlsub-tl)
            rwork = (1.-fiold) / (1.-fi)
            wliql = wliql * rwork
            wsnol = wsnol * rwork
          endif
c
cc  600   continue
      endif
c
c areally average fluxes to be used by soil model. (don't use
c index due to mix call, but only need at all if npn > 0)
c
      if (npn.gt.0) then
c
        rwork = 1. / dtime
cc        do 700 i=1,npoi
          rwork2 = 1. - fiold
          heatg = rwork2*heatg
     >             + fiold*sflo(nsnolay+2)
     >             + zheat*rwork
          solg  = rwork2 * solg
          fvapg = rwork2 * fvapg
          x1    = rwork2 * raing
          x2    = zmelt*rwork
          x3    = tmelt
cc  700   continue
c
        call mix (raing,traing, x1,traing, x2,x3, vzero,vzero)
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine snowheat (tlsub, fhtop, sflo, xl, chl)
c ---------------------------------------------------------------------
c
c sets up call to tridia to solve implicit snow heat conduction,
c using snow temperatures in tsno (in comsno). adds an extra
c buried-lower-veg layer to the bottom of the snow with 
c conduction coefficient conbur/xl and heat capacity chl*xl
c

      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
      include 'comsoi.h'
c
c Arguments
c
      real chl                   ! specific heat of lower veg per l/s area (supplied)
      
      real tlsub,          ! temperature of buried lower veg (supplied, returned)
     >     fhtop,          ! heat flux into top snow layer from atmos (supplied)
     >     sflo(nsnolay+2), ! downward heat flow across layer boundaries (returned)
     >     xl              ! (lai(i,1)+sai(i,1))*fl(i), lower-veg density(supplied)
c
c Local variables
c
      integer k,               ! loop indices
     >        km1,               ! used to avoid layer 0
     >        kp1                ! used to avoid layer nsnolay+2
c
      real rimp,                 ! implicit fraction of the calculation (0 to 1)
     >     conbur,               ! conduction coeff of buried lower veg layer 
     >                           ! for unit density xl=(lai+sai)*fl, in w m-2 k-1
     >     hfake,                ! arbitrary small thickness to allow processing 
     >                           ! for zero snow. (doesn't use index since tridia
     >                           ! not set up for index.)
     >     rwork,                ! to compute matrix diagonals and right-hand side
     >     dt,                   ! '
     >     dti                   ! '
c
      real con(nsnolay+2),  ! conduction coefficents between layers
     >     temp(nsnolay+1), ! combined snow and buried-veg temperatures
     >     d1(nsnolay+1),   ! diagonals of tridiagonal systems of equations 
     >     d2(nsnolay+1),   ! '
     >     d3(nsnolay+1),   ! '
     >     rhs(nsnolay+1),  ! right-hand sides of systems of equations
     >     w1(nsnolay+1),   ! work array needed by tridia
     >     w2(nsnolay+1)    ! '
c
c conbur (for xl=1) is chosen to be equiv to 10 cm of snow
c
      data rimp, conbur, hfake /1.0, 2.0, .01/
c
c copy snow and buried-lower-veg temperatures into combined
c array temp
c
      call scopy (nsnolay, tsno,  temp             )
cc      call scopy (1, tlsub, temp(nsnolay+1))
       temp(nsnolay + 1) = tlsub
c
c set conduction coefficients between layers
c
      do 100 k=1,nsnolay+2
        if (k.eq.1) then
          con(k) = 0.0
c
        else if (k.le.nsnolay) then
          rwork = 0.5 / consno
cc          do 102 i=1,npoi
            con(k) = 1. / (   max(hsno(k-1),hfake)*rwork
     >                        + max(hsno(k)  ,hfake)*rwork )
cc  102     continue
c
        else if (k.eq.nsnolay+1) then
          rwork = 0.5 / consno
cc          do 104 i=1,npoi
            con(k) = 1. / (   max(hsno(k-1),hfake)*rwork
     >                        + 0.5*xl/conbur )
cc  104     continue
c
        else if (k.eq.nsnolay+2) then
          rwork = 0.5 / conbur
cc          do 106 i=1,npoi
            con(k) = 1. / (   xl*rwork
     >                        + 0.5*hsoi(1) / consoi(1) )
cc  106     continue
        endif
  100 continue
c
c set matrix diagonals and right-hand side. for layer nsnolay+1
c (buried-lower-veg layer), use explicit contact with soil, and
c multiply eqn through by xl*chl/dtime to allow zero xl.
c
      do 200 k=1,nsnolay+1
        km1 = max (k-1,1)
        kp1 = min (k+1,nsnolay+1)
c
        if (k.le.nsnolay) then
          rwork = dtime /(rhos*cice)
cc          do 202 i=1,npoi
            dt = rwork / (max(hsno(k),hfake))
            d1(k) =    - dt*rimp* con(k)
            d2(k) = 1. + dt*rimp*(con(k)+con(k+1))
            d3(k) =    - dt*rimp* con(k+1)
c
            rhs(k) = temp(k) + dt
     >               * ( (1.-rimp)*con(k)  *(temp(km1)-temp(k))
     >                 + (1.-rimp)*con(k+1)*(temp(kp1)-temp(k)) )
cc  202     continue
c
          if (k.eq.1) then 
            rwork = dtime /(rhos*cice)
cc            do 204 i=1,npoi
              dt = rwork / (max(hsno(k),hfake))
              rhs(k) = rhs(k) + dt*fhtop
cc  204       continue
          endif
c
        else if (k.eq.nsnolay+1) then
c
          rwork = chl / dtime
cc          do 206 i=1,npoi
            dti = xl*rwork
            d1(k) =     -  rimp* con(k)
            d2(k) = dti +  rimp*(con(k)+con(k+1))
            d3(k) = 0.
            rhs(k) = dti*temp(k)
     >               + ( (1.-rimp)*con(k)*(temp(km1)-temp(k))
     >                 + con(k+1)*(tsoi(1)-(1.-rimp)*temp(k)) )
cc  206     continue
        endif
  200 continue
c
c solve the tridiagonal systems
c
      call tridia (nsnolay+1, d1,d2,d3, rhs, temp, w1,w2)
c
c deduce downward heat fluxes between layers
c
cc      call scopy (1,fhtop, sflo(1))
      sflo(1) = fhtop
c
      do 400 k=1,nsnolay+1
        if (k.le.nsnolay) then
          rwork = rhos*cice/dtime
cc          do 402 i=1,npoi
            sflo(k+1) = sflo(k) - rwork*hsno(k)
     >                                *(temp(k)-tsno(k))
cc  402     continue
c
        else
          rwork = chl/dtime
cc          do 404 i=1,npoi
            sflo(k+1) = sflo(k)
     >                  - xl*rwork*(temp(nsnolay+1)-tlsub)
cc  404     continue
        endif
  400 continue
c
c copy temperature solution to tsno and tlsub, but not for
c points with no snow
c
      do 500 k=1,nsnolay
cc        do 502 i=1,npoi
          if (fi.gt.0.) tsno(k) = temp(k) 
cc  502   continue
  500 continue
c
cc      do 510 i=1,npoi
        if (fi.gt.0.) tlsub = temp(nsnolay+1)
cc  510 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vadapt (hcur, tcur, htop, indp, np, nlay) 
c ---------------------------------------------------------------------
c
c re-adapt snow layer thicknesses, so top thickness
c equals hsnotop and other thicknesses are equal
c
c also adjusts profile of tracer field tcur so its vertical
c integral is conserved (eg, temperature)
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments
c
      integer np,            ! number of snow pts in current strip (supplied)
     >        nlay           ! # of layer
c
      integer indp     ! index of snow pts in current strip (supplied)
c
      real htop              ! prescribed top layer thickness (supplied)
c
      real hcur(nlay),  ! layer thicknesses (supplied and returned)     
     >     tcur(nlay)   ! tracer field (supplied and returned)
c
c local variables
c
      integer i, j, k, ko    ! loop indices
c
      real dz, rwork
c
      real ht,         ! storing variable for zold        
     >     h1,         ! to compute new layer thickness
     >     za,         ! 
     >     zb,         ! 
     >     zheat       !
c
      real hnew(nsnolay),    ! new layer thickness
     >     tnew(nsnolay),    ! new temperatures of layers
     >     zold(nsnolay+1)   ! distances from surface to old layer boundaries
c
c if no snow or seaice points in current 1d strip, return. note
c that the index is not used below (for cray vec and efficiency)
c except in the final loop setting the returned values
c
      if (np.eq.0) return
c
c set distances zold from surface to old layer boundaries
c
cc      call const (zold(1), 0.0)
      zold(1) = 0
c
      do 300 k=1,nlay
cc        do 302 i=1,npoi
          zold(k+1) = zold(k) + hcur(k)
cc  302   continue
  300 continue
c
c set new layer thicknesses hnew (tot thickness is unchanged).
c if total thickness is less than nlay*htop (which should be
c le hsnomin), make all new layers equal including
c top one, so other layers aren't so thin. use epsilon to 
c handle zero (snow) points
c
cc      call scopy (1,zold(nlay+1), ht)
	ht = zold(nlay + 1)
c
      rwork = nlay*htop
cc      do 304 i=1,npoi
        if (ht.ge.rwork) then
          h1 = (ht-htop)/(nlay-1)
        else
          h1 = max (ht/nlay, epsilon)
        endif
cc  304 continue
c
      do 306 k=1,nlay
cc        do 308 i=1,npoi
          hnew(k) = h1
cc  308   continue
  306 continue
c
      rwork = nlay*htop
cc      do 310 i=1,npoi
        if (ht.ge.rwork) hnew(1) = htop
cc  310 continue
c
c integrate old temperature profile (loop 410) over each
c new layer (loop 400), to get new field tnew
c
cc      call const (zb, 0.0)
	zb = 0
c
      do 400 k=1,nlay
c
cc        do 402 i=1,npoi
          za = zb
          zb = za + hnew(k)
cc  402   continue
cc        call const (zheat, 0.0)
	zheat = 0 
c
        do 410 ko=1,nlay
cc          do 412 i=1,npoi
            if (za.lt.zold(ko+1) .and. zb.gt.zold(ko)) then
              dz = min(zold(ko+1),zb) - max(zold(ko),za)
              zheat = zheat + tcur(ko)*dz
            endif
cc  412     continue
  410   continue
c
cc        do 420 i=1,npoi
          tnew(k) = zheat / hnew(k)
cc  420   continue
c
  400 continue
c
c use index for final copy to seaice or snow arrays, to avoid
c changing soil values (when called for seaice) and to avoid
c changing nominal snow values for no-snow points (when called
c for snow)
c
      do 500 k=1,nlay
        do 502 j=1,np
cc          i = indp(j)
          i = indp

          hcur(k) = hnew(k)
          tcur(k) = tnew(k)
  502   continue
  500 continue
c
      return
      end
