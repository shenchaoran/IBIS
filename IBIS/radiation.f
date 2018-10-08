$FIXEDFORMLINESIZE:132
c #####     ##    #####      #      ##     #####     #     ####   #    #
c #    #   #  #   #    #     #     #  #      #       #    #    #  ##   #
c #    #  #    #  #    #     #    #    #     #       #    #    #  # #  #
c #####   ######  #    #     #    ######     #       #    #    #  #  # #
c #   #   #    #  #    #     #    #    #     #       #    #    #  #   ##
c #    #  #    #  #####      #    #    #     #       #     ####   #    #
c
c ---------------------------------------------------------------------
      subroutine solset
c ---------------------------------------------------------------------
c
c zeros albedos and internal absorbed solar fluxes, and sets
c index for other solar routines. the index indsol, with number
c of points nsol, points to current 1d strip arrays whose coszen 
c values are gt 0 (indsol, nsol are in com1d)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c
c zero albedos returned just as a niceity
c
      call const (asurd, nband, 0.0)
      call const (asuri, nband, 0.0)
c
c zeros absorbed solar fluxes sol[u,s,l,g,i]1 since only points
c with +ve coszen will be set in solarf, and since
c sol[u,l,s,g,i]1 are summed over wavebands in solarf
c
c similarly zero par-related arrays set in solarf for turvap
c
      solu = 0.0
      sols = 0.0
      soll = 0.0
	solg = 0.0
      soli = 0.0
c
      topparu = 0.0
      topparl = 0.0
c
c set canopy scaling coefficients for night-time conditions
c
      call const (scalcoefl, 4, 0.0)
      call const (scalcoefu, 4, 0.0)
c
c set index of points with positive coszen
c
      nsol = 0
c
cc      do 300 i = 1, npoi
        if (coszen.gt.0.) then
          nsol = nsol + 1
          indsol = 1
        endif
cc  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solsur (ib)
c ---------------------------------------------------------------------
c
c sets surface albedos for soil and snow, prior to other
c solar calculations
c
c ib = waveband number
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'com1d.h'
c
c input variable
c
      integer ib    ! waveband number. 1 = visible, 2 = near IR
c
c local variables
c     
cc      integer j    ! loop indice on number of points with >0 coszen
cc     >        i     ! indice of point in (1, npoi) array. 
c
      real a7svlo,  ! snow albedo at low threshold temp., visible
     >     a7snlo,  !                                   , near IR
     >     a7svhi,  !                 high              , visible
     >     a7snhi,  !                                   , near-IR
     >     t7shi,   ! high threshold temperature for snow albed
     >     t7slo,   ! low  threshold temperature for snow albedo
     >     dinc,    ! albedo correction du to soil moisture
     >     zw       ! liquid moisture content

      real x, zfac
c
c set the "standard" snow values:
c
      data a7svlo, a7svhi /0.90, 0.70/
      data a7snlo, a7snhi /0.60, 0.40/
c
c     t7shi ... high threshold temperature for snow albedo
c     t7slo ... low  threshold temperature for snow albedo
c
      t7shi = tmelt
      t7slo = tmelt - 15.0
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) then
        return
      endif
c
      if (ib.eq.1) then
c
c soil albedos (visible waveband)
c
cc        do 100 j = 1, nsol
c
cc          i = indsol(j)
c
c change the soil albedo as a function of soil moisture
c
          zw = wsoi(1) * (1.-wisoi(1))
c
          dinc = 1.0 + 1.0 * min (1., max (0.0, 1. - (zw /.50) ))
c
          albsod = min (albsav * dinc, .80)
          albsoi = albsod
c
cc  100   continue
c
c snow albedos (visible waveband)
c
cc        do 110 j = 1, nsol
c
cc          i = indsol(j)
c
          x = (a7svhi*(tsno(1)-t7slo) + a7svlo*(t7shi-tsno(1)))
     >           / (t7shi-t7slo)
c
          x = min (a7svlo, max (a7svhi, x))
c
          zfac   = max ( 0., 1.5 / (1.0 + 4.*coszen) - 0.5 )
          albsnd = min (0.99, x + (1.-x)*zfac)
          albsni = min (1., x)
c
cc  110   continue
c
      else
c
c soil albedos (near-ir waveband)
c
cc        do 200 j = 1, nsol
cc          i = indsol(j)
c
c lsx.2 formulation (different from lsx.1)
c
          zw = wsoi(1) * (1. - wisoi(1))
c
          dinc = 1.0 + 1.0 * min (1., max (0.0, 1.0 - (zw / .50)  ))
c
          albsod = min (albsan * dinc, .80)
          albsoi = albsod
c
cc  200   continue
c
c snow albedos (near-ir waveband)
c
cc        do 210 j = 1, nsol
c
cc          i = indsol(j)
c
          x = (a7snhi*(tsno(1)-t7slo) + a7snlo*(t7shi-tsno(1)))
     >           / (t7shi-t7slo)
          x = min (a7snlo, max (a7snhi, x))
c
          zfac = max ( 0., 1.5/(1.+4.*coszen) - 0.5 )
c
          albsnd = min (0.99, x + (1.-x)*zfac)
          albsni = min (1., x)
c
cc  210   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solalb (ib)
c ---------------------------------------------------------------------
c
c calculates effective albedos of the surface system,
c separately for unit incoming direct and diffuse flux -- the 
c incoming direct zenith angles are supplied in comatm array 
c coszen, and the effective albedos are returned in comatm
c arrays asurd, asuri -- also detailed absorbed and reflected flux
c info is stored in com1d arrays, for later use by solarf
c
c the procedure is first to calculate the grass+soil albedos,
c then the tree + (grass+soil+snow) albedos. the labels
c (a) to (d) correspond to those in the description doc
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c 
c Arguments
c 
      integer ib     ! waveband number (1= visible, 2= near-IR)
c
c local variables
c     
cc      integer j    ! loop indice on number of points with >0 coszen
cc     >        i     ! indice of point in (1, npoi) array. 
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c (a) obtain albedos, etc, for two-stream lower veg + soil
c     system, for direct and diffuse incoming unit flux
c
cc      do 100 j = 1, nsol
c
cc        i = indsol(j)
c
        asurd(ib) = albsod
        asuri(ib) = albsoi
c
cc  100 continue
c
      call twostr (ablod, abloi,  relod, reloi,  flodd,  dummy,
     >             flodi, floii,  asurd,  asuri,    1,   coszen, ib)
c
c (b) areally average surface albedos (lower veg, soil, snow)
c
cc      do 200 j = 1, nsol
c
cc        i = indsol(j)
c
        asurd(ib) = fl*(1.-fi)*relod
     >              + (1.-fl)*(1.-fi)*albsod
     >              + fi*albsnd   
c
        asuri(ib) = fl*(1.-fi)*reloi
     >              + (1.-fl)*(1.-fi)*albsoi
     >              + fi*albsni    
c
cc  200 continue
c
c (c) obtain albedos, etc, for two-stream upper veg + surface
c     system, for direct and diffuse incoming unit flux
c
      call twostr (abupd, abupi,  reupd, reupi,  fupdd,  dummy,
     >             fupdi, fupii,  asurd,  asuri,    2,   coszen, ib)
c
c (d) calculate average overall albedos 
c
cc      do 300 j = 1, nsol
c
cc        i = indsol(j)
c
        asurd(ib) = fu*reupd + (1.-fu)*asurd(ib)
c
        asuri(ib) = fu*reupi + (1.-fu)*asuri(ib)
c
cc  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine solarf (ib)
c ---------------------------------------------------------------------
c
c calculates solar fluxes absorbed by upper and lower stories,
c soil and snow
c
c zenith angles are in comatm array coszen, and must be the same
c as supplied earlier to solalb
c
c solarf uses the results obtained earlier by solalb and 
c stored in com1d arrays. the absorbed fluxes are returned in
c com1d arrays sol[u,s,l,g,i]
c
c the procedure is first to calculate the upper-story absorbed
c fluxes and fluxes below the upper story, then the lower-story
c absorbed fluxes and fluxes below the lower story, then fluxes
c absorbed by the soil and snow
c
c ib = waveband number
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c 
c Arguments
c 
      integer ib     ! waveband number (1= visible, 2= near-IR)
c
c local variables
c     
cc      integer j    ! loop indice on number of points with >0 coszen
cc     >        i     ! indice of point in (1, npoi) array. 
c
      real x, y, xd, xi, 
     >     xaiu,    ! total single-sided lai+sai, upper
     >     xail     ! total single-sided lai+sai, lower
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c (f) calculate fluxes absorbed by upper leaves and stems,
c     and downward fluxes below upper veg, using unit-flux
c     results of solalb(c) (apportion absorbed flux between
c     leaves and stems in proportion to their lai and sai)
c
cc     do 600 j=1,nsol
c
cc        i = indsol(j)
        x = solad(ib)*abupd + solai(ib)*abupi
        y = lai(2) / max (lai(2)+sai(2), epsilon)
        solu = solu + x * y
        sols = sols + x * (1.-y)
        sol2d = solad(ib)*fupdd
        sol2i = solad(ib)*fupdi + solai(ib)*fupii
c
cc  600 continue
c
c (g) areally average fluxes to lower veg, soil, snow
c
cc      do 700 j=1,nsol
c
cc        i = indsol(j)
        sol3d = fu*sol2d + (1.-fu)*solad(ib)
        sol3i = fu*sol2i + (1.-fu)*solai(ib)
c
cc  700 continue
c
c (h,i) calculate fluxes absorbed by lower veg, snow-free soil
c       and snow, using results of (g) and unit-flux results
c       of solalb(a)
c
cc      do 800 j=1,nsol
c
cc        i = indsol(j)
        soll = soll + sol3d*ablod + sol3i*abloi
c
        xd = (fl*flodd + 1.-fl) * sol3d
c
        xi = fl*(sol3d*flodi + sol3i*floii) + (1.-fl) * sol3i
c
        solg = solg + (1.-albsod)*xd + (1.-albsoi)*xi
c
        soli = soli + (1.-albsnd)*sol3d + (1.-albsni)*sol3i
c
cc  800 continue
c
c estimate absorbed pars at top of canopy, toppar[u,l] and
c some canopy scaling parameters
c
c this neglects complications due to differing values of dead vs 
c live elements, averaged into rhoveg, tauveg in vegdat, and 
c modifications of omega due to intercepted snow in twoset
c
c do only for visible band (ib=1)
c
      if (ib.eq.1) then
c
cc        do 900 j = 1, nsol
c
cc          i = indsol(j)
c
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c
c some of the required terms (i.e. term[u,l] are calculated in the subroutine 'twostr'.
c in the equations below, 
c
c   A = scalcoefu(i,1) = term[u,l](i,1) * ipardir(0)
c   B = scalcoefu(i,2) = term[u,l](i,2) * ipardir(0) + term[u,l](i,3) * ipardif(0)
c   C = scalcoefu(i,3) = term[u,l](i,4) * ipardir(0) + term[u,l](i,5) * ipardif(0)
c   A + B + C = scalcoefu(i,4) = also absorbed par at canopy of canopy by leaves & stems
c
c upper canopy:
c
c total single-sided lai+sai
c
          xaiu = max (lai(2)+sai(2), epsilon)
c
c some terms required for use in canopy scaling:
c
          scalcoefu(1) = termu(1) * solad(ib)
c
          scalcoefu(2) = termu(2) * solad(ib) + termu(3) * solai(ib)
c
          scalcoefu(3) = termu(4) * solad(ib) + termu(5) * solai(ib)
c
          scalcoefu(4) = scalcoefu(1) +  scalcoefu(2) +  scalcoefu(3)
c
c apar of the "top" leaves of the canopy
c
          topparu = scalcoefu(4) * lai(2) / xaiu
c
c lower canopy:
c
c total single-sided lai+sai
c
          xail = max (lai(1)+sai(1), epsilon)
c
c some terms required for use in canopy scaling:
c
          scalcoefl(1) = terml(1) * sol3d
c
          scalcoefl(2) = terml(2) * sol3d + terml(3) * sol3i
c
          scalcoefl(3) = terml(4) * sol3d + terml(5) * sol3i
c
          scalcoefl(4) = scalcoefl(1) +
     >                     scalcoefl(2) +
     >                     scalcoefl(3)
c
c apar of the "top" leaves of the canopy
c
          topparl = scalcoefl(4) * lai(1) / xail
c
cc  900   continue
c
      endif
c
      return
      end
c
c
c ------------------------------------------------------------------------
      subroutine twostr (abvegd, abvegi, refld, refli, fbeldd, fbeldi,
     >                   fbelid, fbelii, asurd, asuri, iv, coszen, ib)
c ------------------------------------------------------------------------
c
c solves canonical radiative transfer problem of two-stream veg
c layer + underlying surface of known albedo, for unit incoming
c direct or diffuse flux. returns flux absorbed within layer,
c reflected flux, and downward fluxes below layer. note that all
c direct fluxes are per unit horizontal zrea, ie, already 
c including a factor cos (zenith angle)
c
c the solutions for the twostream approximation follow Sellers (1985),
c and Bonan (1996) (the latter being the LSM documentation)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments
c
      integer ib,             ! waveband number (1= visible, 2= near-IR)
     >        iv              ! 1 for lower, 2 for upper story params (supplied)
c
      real abvegd,      ! direct flux absorbed by two-stream layer (returned)
     >     abvegi,      ! diffuse flux absorbed by two-stream layer (returned)
     >     refld,       ! direct flux reflected above two-stream layer (returned)
     >     refli,       ! diffuse flux reflected above two-stream layer (returned)
     >     fbeldd,      ! downward direct  flux below two-stream layer(returned)
     >     fbeldi,      ! downward direct  flux below two-stream layer(returned)
     >     fbelid,      ! downward diffuse flux below two-stream layer(returned)
     >     fbelii,      ! downward diffuse flux below two-stream layer(returned)
     >     asurd(nband), ! direct  albedo of underlying surface (supplied)
     >     asuri(nband), ! diffuse albedo of underlying surface (supplied)
     >     coszen       ! cosine of direct zenith angle (supplied, must be gt 0)
c
c local variables
c
cc      integer j    ! loop indice on number of points with >0 coszen
cc     >        i     ! indice of point in (1, npoi) array. 
c
      real b, c, c0, d, f, h, k, q, p, sigma
c
      real ud1, ui1, ud2, ui2, ud3, xai, s1, s2, p1, p2, p3, p4, 
     >     rwork, dd1, di1, dd2, di2, h1, h2, h3, h4, h5, h6, h7, h8, 
     >     h9, h10, absurd, absuri
c
c [d,i] => per unit incoming direct, diffuse (indirect) flux
c
      real omega,       !
     >     betad,       !
     >     betai,       !
     >     avmu,        !
     >     gdir,        !
     >     tmp0         !
c
c do nothing if all points in current strip have coszen le 0
c
      if (nsol.eq.0) return
c
c calculate two-stream parameters omega, betad, betai, avmu, gdir
c
      call twoset (omega, betad, betai, avmu, gdir, coszen, iv, ib)
c
cc      do 100 j=1,nsol
c
cc        i = indsol(j)
c
c the notations used here are taken from page 21 of Bonan's LSM documentation:
c Bonan, 1996: A Land Surface Model (LSM version 1.0) for ecological, hydrological,
c and atmospheric studies: Technical description and user's guide. NCAR Technical
c Note. NCAR/TN-417+STR, January 1996.
c
c some temporary variables are also introduced, which are from the original
c lsx model.
c
        b = 1. - omega * (1.-betai)
        c = omega * betai
c
        tmp0 = b*b-c*c
c
        q = sqrt ( max(0.0, tmp0) )
        k = gdir / max(coszen, 0.01)
        p = avmu * k
c
c next line perturbs p if p = q
c
        if ( abs(p-q) .lt. .001*p )
     >  p = (1.+sign(.001,p-q)) * p
c
        c0 = omega * p
        d = c0 * betad
        f = c0 * (1.-betad)
        h = q / avmu
c
        sigma = p*p - tmp0
c
c direct & diffuse parameters are separately calculated
c
        ud1 = b - c/asurd(ib)
        ui1 = b - c/asuri(ib)
        ud2 = b - c*asurd(ib)
        ui2 = b - c*asuri(ib)
        ud3 = f + c*asurd(ib)
c
        xai = max (lai(iv) + sai(iv), epsilon)
c
        s1 = exp(-1.*h*xai)
        s2 = exp(-1.*k*xai)
c
        p1 = b + q
        p2 = b - q
        p3 = b + p
        p4 = b - p
        rwork = 1./s1
c
c direct & diffuse parameters are separately calculated
c
        dd1 = p1*(ud1-q)*rwork - p2*(ud1+q)*s1
        di1 = p1*(ui1-q)*rwork - p2*(ui1+q)*s1
        dd2 = (ud2+q)*rwork - (ud2-q)*s1
        di2 = (ui2+q)*rwork - (ui2-q)*s1
        h1 = -1.*d*p4 - c*f
        rwork = s2*(d-c-h1*(ud1+p)/sigma)
        h2 = 1./dd1*( (d-h1*p3/sigma)*(ud1-q)/s1 - 
     >       p2*rwork )
        h3 = -1./dd1*( (d-h1*p3/sigma)*(ud1+q)*s1 - 
     >       p1*rwork )
        h4 = -1.*f*p3 - c*d
        rwork = s2*(ud3-h4*(ud2-p)/sigma)
        h5 = -1./dd2*( h4*(ud2+q)/(sigma*s1) +
     >       rwork )
        h6 = 1./dd2*( h4*s1*(ud2-q)/sigma +
     >       rwork )
        h7 = c*(ui1-q)/(di1*s1)
        h8 = -1.*c*s1*(ui1+q)/di1
        h9 = (ui2+q)/(di2*s1)
        h10= -1.*s1*(ui2-q)/di2
c
c save downward direct, diffuse fluxes below two-stream layer
c
        fbeldd = s2
        fbeldi = 0.
        fbelid = h4/sigma*s2 + h5*s1 + h6/s1
        fbelii = h9*s1 + h10/s1
c
c save reflected flux, and flux absorbed by two-stream layer
c
        refld = h1/sigma + h2 + h3
        refli = h7 + h8
        absurd = (1.-asurd(ib)) * fbeldd
     >         + (1.-asuri(ib)) * fbelid
        absuri = (1.-asuri(ib)) * fbelii
c
        abvegd = max (0., 1. - refld - absurd)
        abvegi = max (0., 1. - refli - absuri)
c
c if no veg, make sure abveg (flux absorbed by veg) is exactly zero
c if this is not done, roundoff error causes small (+/-)
c sols, soll values in solarf and subsequent problems in turvap
c via stomata
c
        if (xai.lt.epsilon) abvegd = 0.0
        if (xai.lt.epsilon) abvegi = 0.0
c
c some terms needed in canopy scaling
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c
c in the equations below, 
c
c   k = term[u,l](6)
c   h = term[u,l](7)
c
c   A = term[u,l](1) * ipardir(0)
c   B = term[u,l](2) * ipardir(0) + term[u,l](3) * ipardif(0)
c   C = term[u,l](4) * ipardir(0) + term[u,l](5) * ipardif(0)
c
c calculations performed only for visible (ib=1)
c
      if (ib.eq.1) then
c
        if (iv.eq.1) then
          terml(1) = k * (1. + (h4-h1) / sigma)
          terml(2) = h * (h5 - h2)
          terml(3) = h * (h9 - h7)
          terml(4) = h * (h3 - h6)
          terml(5) = h * (h8 - h10)
          terml(6) = k
          terml(7) = h
        else
          termu(1) = k * (1. + (h4-h1) / sigma)
          termu(2) = h * (h5 - h2)
          termu(3) = h * (h9 - h7)
          termu(4) = h * (h3 - h6)
          termu(5) = h * (h8 - h10)
          termu(6) = k
          termu(7) = h
        endif
c
      end if
c
cc  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine twoset (omega, betad, betai, avmu, gdir,
     >                   coszen, iv, ib)  
c ---------------------------------------------------------------------
c
c sets two-stream parameters, given single-element transmittance
c and reflectance, leaf orientation weights, and cosine of the
c zenith angle, then adjusts for amounts of intercepted snow
c
c the two-stream parameters omega,betad,betai are weighted 
c combinations of the "exact" values for the 3 orientations:
c all vertical, all horizontal, or all random (ie, spherical)
c
c the vertical, horizontal weights are in oriev,orieh (comveg)
c
c the "exact" expressions are as derived in my notes(8/6/91,p.6).
c note that values for omega*betad and omega*betai are calculated
c and then divided by the new omega, since those products are 
c actually used in twostr. also those depend *linearly* on the
c single-element transmittances and reflectances tauveg, rhoveg,
c which are themselves linear weights of leaf and stem values 
c
c for random orientation, omega*betad depends on coszen according
c to the function in array tablemu
c
c the procedure is approximate since omega*beta[d,i] and gdir
c should depend non-linearly on the complete leaf-angle
c distribution. then we should also treat leaf and stem angle
c distributions separately, and allow for the cylindrical
c shape of stems (norman and jarvis, app.b; the expressions 
c below are appropriate for flat leaves)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (all quantities are returned unless otherwise note)
c
      integer ib,             ! waveband number (1= visible, 2= near-IR)
     >        iv              ! 1 for lower, 2 for upper story params (supplied)
c
      real omega,       ! fraction of intercepted radiation that is scattered
     >     betad,       ! fraction of scattered *direct* radiation that is
     >                        !  scattered into upwards hemisphere
     >     betai,       ! fraction of scattered downward *diffuse* radiation
     >                        ! that is scattered into upwards hemisphere (or fraction
     >                        ! of scattered upward diffuse rad. into downwards hemis)
     >     avmu,        ! average diffuse optical depth
     >     gdir,        ! average projected leaf area into solar direction
     >     coszen       ! cosine of solar zenith angle (supplied)
c
c local variables
c
      integer     ! loop indice on number of points with >0 coszen
     >        ntmu, !
     >        itab
cc     >        i,    ! indice of point in (1, npoi) array. 

c
      real zrho, ztau, orand, ztab, rwork, y, o, x, betadsno, betaisno
c
      real otmp
c
      parameter (ntmu=100)
      real tablemu(ntmu+1), omegasno(nband)
      save tablemu, omegasno, betadsno, betaisno
c
      data tablemu /
     >   0.5000, 0.4967, 0.4933, 0.4900, 0.4867, 0.4833, 0.4800, 0.4767,
     >   0.4733, 0.4700, 0.4667, 0.4633, 0.4600, 0.4567, 0.4533, 0.4500,
     >   0.4467, 0.4433, 0.4400, 0.4367, 0.4333, 0.4300, 0.4267, 0.4233,
     >   0.4200, 0.4167, 0.4133, 0.4100, 0.4067, 0.4033, 0.4000, 0.3967,
     >   0.3933, 0.3900, 0.3867, 0.3833, 0.3800, 0.3767, 0.3733, 0.3700,
     >   0.3667, 0.3633, 0.3600, 0.3567, 0.3533, 0.3500, 0.3467, 0.3433,
     >   0.3400, 0.3367, 0.3333, 0.3300, 0.3267, 0.3233, 0.3200, 0.3167,
     >   0.3133, 0.3100, 0.3067, 0.3033, 0.3000, 0.2967, 0.2933, 0.2900,
     >   0.2867, 0.2833, 0.2800, 0.2767, 0.2733, 0.2700, 0.2667, 0.2633,
     >   0.2600, 0.2567, 0.2533, 0.2500, 0.2467, 0.2433, 0.2400, 0.2367,
     >   0.2333, 0.2300, 0.2267, 0.2233, 0.2200, 0.2167, 0.2133, 0.2100,
     >   0.2067, 0.2033, 0.2000, 0.1967, 0.1933, 0.1900, 0.1867, 0.1833,
     >   0.1800, 0.1767, 0.1733, 0.1700, 0.1667 /
c
      data omegasno /0.9, 0.7/
      data betadsno, betaisno /0.5, 0.5/
c
c set two-stream parameters omega, betad, betai, gdir and avmu
c as weights of those for 100% vert,horiz,random orientations
c
cc      do 100 j=1,nsol
cc        i = indsol(j)
c
        zrho = rhoveg(ib,iv)
        ztau = tauveg(ib,iv)
c
c weight for random orientation is 1 - those for vert and horiz
c
        orand = 1. - oriev(iv) - orieh(iv)
c
        omega = zrho + ztau
c
c ztab is transmittance coeff - for random-orientation omega*betad,
c given by tablemu as a function of coszen
c
        itab = nint (coszen*ntmu + 1)
        ztab = tablemu(itab)
        rwork = 1./omega
c
        betad = (  oriev(iv) * 0.5*(zrho + ztau)
     >              + orieh(iv) * zrho
     >              + orand       * ((1.-ztab)*zrho + ztab*ztau) )
     >             * rwork
c
        betai = (  oriev(iv) * 0.5*(zrho + ztau)
     >              + orieh(iv) * zrho
     >              + orand       * ((2./3.)*zrho + (1./3.)*ztau) )
     >             * rwork
c
        gdir  = oriev(iv) * (2./pi) *
     >             sqrt ( max (0., 1.-coszen*coszen) )
     >           + orieh(iv) * coszen
     >           + orand       * 0.5
c
        avmu = 1.
c
cc  100 continue
c
c adjust omega, betad and betai for amounts of intercepted snow
c (omegasno decreases to .6 of cold values within 1 deg of tmelt)
c
      if (iv.eq.1) then
c
c lower story
c
cc        do 210 j=1,nsol
cc          i = indsol(j)
          y = fwetl*(1.-rliql)
          o = omegasno(ib)*(.6 + .4*max(0.,min(1.,(tmelt-tl)/1.0)))
          otmp  = omega
          rwork = y * o
          omega =  (1-y)*otmp          + rwork
          betad = ((1-y)*otmp*betad + rwork*betadsno) /
     >               omega  
          betai = ((1-y)*otmp*betai + rwork*betaisno) /
     >               omega  
cc  210   continue
c
      else
c
c upper story
c
cc        do 220 j=1,nsol
cc          i = indsol(j)
          x = lai(iv) / max (lai(iv)+sai(iv), epsilon)
          y = x * fwetu*(1.-rliqu) + (1-x) *fwets*(1.-rliqs)
          o = (     x  * min (1., max (.6, (tmelt-tu)/0.1))
     >         + (1-x) * min (1., max (.6, (tmelt-ts)/0.1)) )
     >      *  omegasno(ib) 
c
          otmp  = omega
          rwork = y * o
          omega =  (1-y)*otmp          + rwork
          betad = ((1-y)*otmp*betad + rwork*betadsno) /
     >               omega
          betai = ((1-y)*otmp*betai + rwork*betaisno) /
     >               omega
c
cc  220   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine irrad
c ---------------------------------------------------------------------
c
c calculates overall emitted ir flux, and net absorbed minus
c emitted ir fluxes for upper leaves, upper stems, lower story,
c soil and snow. assumes upper leaves, upper stems and lower
c story each form a semi-transparent plane, with the upper-leaf
c plane just above the upper-stem plane. the soil and snow 
c surfaces have emissivities of 0.95.
c
c the incoming flux is supplied in comatm array fira
c
c the emitted ir flux by overall surface system is returned in
c com1d array firb - the ir fluxes absorbed by upper leaves,
c upper stems, lower veg, soil and snow are returned in com1d 
c arrays firu, firs, firl, firg and firi
c 
c other com1d arrays used are:
c
c emu, ems, eml  = emissivities of the vegetation planes
c fup, fdown     = upward and downward fluxes below tree level
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Local arrays:
c
cc      integer i           ! loop indice
c
      real emisoil,       ! soil emissivity
     >     emisnow,       ! snow emissivity
     >     avmuir         ! average diffuse optical depth
c
      real emu,     ! ir emissivity of upper-leaves veg plane
     >     ems,     ! ir emissivity of upper-stems veg plane
     >     eml,     ! ir emissivity of lower-story veg plane
     >     emg,     ! ir emissivity (gray) of soil surface
     >     emi,     ! ir emissivity (gray) of snow surface
     >     fdown,   ! downward ir flux below tree level per overall area
     >     fdowng,  ! upward   ir flux below tree level per overall area
     >     fup,     ! downward ir flux below lower-story veg
     >     fupg,    ! upward   ir flux below lower-story veg
     >     fupgb,   ! upward   ir flux above bare soil surface
     >     fupi     ! upward   ir flux above snow surface
c
c set emissivities of soil and snow
c
      data emisoil, emisnow
     >    /0.95, 0.95/
c
c use uniform value 1.0 for average diffuse optical depth
c (although an array for solar, all values are set to 1 in twoset).
c
      save avmuir
      data avmuir /1./
c
cc      do 100 i=1,npoi
c
        emu = 1. - exp ( -lai(2) / avmuir )
        ems = 1. - exp ( -sai(2) / avmuir )
        eml = 1. - exp ( -(lai(1)+sai(1)) / avmuir )
c
        emg = emisoil
        emi = emisnow
c
        ts = min(max(ts, 253.16), 323.16)
        tg = min(max(tg, 253.16), 323.16)
        ti = min(max(ti, 253.16), 323.16)
        tl = min(max(tl, 253.16), 323.16)
        tu = min(max(tu, 253.16), 323.16)



        fdown =  (1.-fu) * fira
     >            + fu * ( (1.-emu)*(1.-ems)*fira
     >                       +    emu* (1.-ems)*stef*(tu**4)
     >                       +    ems*stef*(ts**4) )
c
        fdowng = (1.-eml)*fdown  + eml*stef*(tl**4)
c
        fupg   = (1.-emg)*fdowng + emg*stef*(tg**4)
c
        fupgb  = (1.-emg)*fdown  + emg*stef*(tg**4)
c
        fupi   = (1.-emi)*fdown  + emi*stef*(ti**4)
c
        fup = (1.-fi)*(      fl*(       eml *stef*(tl**4)
     >                                     + (1.-eml)*fupg )
     >                        +(1.-fl)*fupgb
     >                      )
     >         +     fi * fupi
c
        firb =   (1.-fu) * fup
     >            + fu  * ( (1.-emu)*(1.-ems)*fup
     >                        +    emu*stef*(tu**4)
     >                        +    ems*(1.-emu)*stef*(ts**4) )
c
        firu =   emu*ems*stef*(ts**4)
     >            + emu*(1.-ems)*fup
     >            + emu*fira
     >            - 2*emu*stef*(tu**4)
c
        firs =   ems*emu*stef*(tu**4)
     >            + ems*fup
     >            + ems*(1.-emu)*fira
     >            - 2*ems*stef*(ts**4)
c
        firl =   eml*fdown
     >            + eml*fupg
     >            - 2*eml*stef*(tl**4)
c
        firg =       fl  * (fdowng - fupg)
     >            + (1.-fl) * (fdown  - fupgb)
c
        firi =   fdown - fupi
c
cc  100 continue
c
      return
      end
