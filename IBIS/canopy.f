$FIXEDFORMLINESIZE:132
c  ####     ##    #    #   ####   #####    #   #
c #    #   #  #   ##   #  #    #  #    #    # #
c #       #    #  # #  #  #    #  #    #     #
c #       ######  #  # #  #    #  #####      #
c #    #  #    #  #   ##  #    #  #          #
c  ####   #    #  #    #   ####   #          #
c
c ---------------------------------------------------------------------
      subroutine canopy(iyear, imonth, iday, pft)
c ---------------------------------------------------------------------
c
c calculates sensible heat and moisture flux coefficients,
c and steps canopy temperatures through one timestep
c
c atmospheric conditions at za are supplied in comatm
c arrays ta, qa, psurf and scalar siga (p/ps)
c
c downward sensible heat and moisture fluxes at za
c are returned in com1d arrays fsena, fvapa
c
c sensible heat and moisture fluxes from solid objects to air
c are stored (for other models and budget) in com1d arrays
c fsen[u,s,l,g,i], fvap[u,s,l,g,i]
c
c the procedure is first to compute wind speeds and aerodynamic
c transfer coefficients in turcof, then call turvap to solve an
c implicit linear system for temperatures and specific
c humidities and the corresponding fluxes - this is iterated
c niter times for non-linearities due to stratification,
c implicit/explicit (h2o phase), dew, vpd and max soil
c moisture uptake - t12 and q12 are changed each iteration,
c and tu, ts, tl, tg, ti can be adjusted too
c
c initialize aerodynamic quantities
c
       include 'implicit.h'
c
c Local variables
c
      integer niter,       ! total number of ierations
     >        iter,         ! number of iteration
     >        pft

	integer iyear, imonth, iday

      real scaleu,scalel

      call canini
c
c estimate soil moisture stress parameters
c
      call drystress (iyear, imonth, iday)
c
c iterate the whole canopy physics solution niter times:
c
cc      niter = 3

	niter = 1
c
      do 100 iter = 1, niter
c
c calculate wind speeds and aerodynamic transfer coeffs
c
        call turcof (iter)
cc	  write(100,*)"turcof over", pft
c
c calculate canopy photosynthesis rates and conductance
c
        call scaler(scaleu, scalel) ! Yuan added this

cc        write(100,*)"scaler over", scalel, scaleu, pft

	  if ((pft.le.3).or.(pft.eq.5).or.(pft.eq.7)) then            
	      call stomata_B (iyear,imonth,iday,scaleu)
        endif

        if ((pft.eq.4).or.(pft.eq.6).or.(pft.eq.8)) then           
	     call stomata_N (iyear, scaleu) 
        endif

	  if((pft.eq.9).or.(pft.eq.10)) then
	      call stomata_S (scalel)
	  end if

	  if(pft.eq.12) then
	      call stomata_C3 (scalel,iyear)
	  end if

	  if(pft.eq.11) then
	      call stomata_C4 (scalel)
	  end if

cc	  write(100,*)"stomata over"
c
c solve implicit system of heat and water balance equations
c
        call turvap (iter, niter)
cc	  write(100,*)"turvap over"
c
  100 continue

      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine canini
c ---------------------------------------------------------------------
c
c initializes aerodynamic quantities that remain constant 
c through one timestep
c
c note that some quantities actually are
c constant as long as the vegetation amounts and fractional
c coverage remain unchanged, so could re-arrange code for
c efficiency - currently all arrays initialized here are in
c com1d which can be overwritten elsewhere
c
c rwork is used throughout as a scratch variable to reduce number of
c computations
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
c Local variables
c
c
      real siga,      ! sigma level of atmospheric data
     >     pa,        ! pressure at level of atmospheric data
     >     x,         ! density of vegetation (without distinction between
     >                ! lai,sai)
     >     x1,        ! density of vegetation (different max)
     >     rwork,     ! difference between top and bottom of canopy 
     >     cvegl,     !
     >     dvegl,     ! diffusion coefficient for lower canopy
     >     bvegl,     ! e-folding depth in canopy for lower canopy
     >     cvegu,     !
     >     dvegu,     ! diffusion coefficient for upper canopy
     >     bvegu      ! e-folding depth in canopy for upper canopy
c
c define sigma level of atmospheric data
c
c currently, the value of siga is set to 0.999. This is roughly 10 meters
c above ground, which is the typical height for the CRU05 input wind speed data
c
      siga = 0.999
c
      tfac = 1.0 / (siga**cappa)
c
c atmospheric conditions at za
c za is variable, although siga = p/ps is constant
c
cc      do 100 i = 1, npoi
c
        pa = psurf * siga
c
        rhoa = pa / ( rair * ta * 
     >            (1.0 + (rvap / rair - 1.0) * qa) )
c
        cp = cair * (1.0 + (cvap / cair - 1.0) * qa)
c
        za = (psurf - pa) / (rhoa * grav)
c
c make sure that atmospheric level is higher than canopy top
c
        za = max (za, ztop(2) + 1.0)
c
cc 100  continue 
c
c aerodynamic coefficients for the lower story
c
c cvegl (drag coeff for momentum) is proportional, and dvegl
c (diffusion coeff for momentum) inversely proportional,
c to x = density of vegetation (without distinction between
c lai,sai and fl*(1-fi)) - x is not allowed to be exactly
c zero to avoid divide-by-zeros, and for x>1 dvegl is 
c proportional to 1/x**2 so that roughness length tends to
c zero as x tends to infinity
c
c also the top, bottom and displacement heights z3,z4,
c displ tend to particular values as the density tends to
c zero, to give same results as equations for no veg at all.
c
cc      do 200 i = 1, npoi
c
        x = fl * (1.0 - fi) * 2.0 * (lai(1) + sai(1)) / alaiml
c
        x  = min (x, 3.0)
        x1 = min (x, 1.0)
c
        rwork = max(ztop(1)-zbot(1),0.01)
        cvegl = (0.4 / rwork) *
     >           max(1.e-5, x)
c
        dvegl = (0.1 * rwork) / 
     >           max(1.e-5, x, x**2)
c
c e-folding depth in canopy
c
        bvegl = sqrt (2.0 * cvegl / dvegl )
c
c [(tau/rho)/u**2] for inf canopy
c
        bdl = 0.5 * bvegl * dvegl
c
c 1 / diffusion coefficient
c
        dil = 1. / dvegl
c
        rwork = (1.0 - x1) * (max (z0soi,z0sno) + 0.01) 
c
        z3 = x1 * ztop(1) + rwork
c
        z4 = x1 * zbot(1) + rwork
c
        z34 = 0.5 * (z3 + z4)
c
        exphl = exp (0.5 * bvegl * (z3-z4))
        expl  = exphl**2
c
        displ = x1 * 0.7 * z3
c
cc 200  continue 
c
c aerodynamic coefficients for the upper story
c same comments as for lower story
c
cc      do 300 i = 1, npoi
c
        x = fu * 2.0 * (lai(2)+sai(2)) / alaimu
c
        x  = min (x, 3.0)
        x1 = min (x, 1.0)
c
        rwork = max(ztop(2)-zbot(2),.01)
        cvegu = (0.4 / rwork) * 
     >           max(1.e-5,x)
c
        dvegu = (0.1 * rwork) / 
     >           max(1.e-5,x,x**2)
c
        rwork = 1. / dvegu
        bvegu  = sqrt (2.0 * cvegu * rwork)
        bdu = 0.5 * bvegu * dvegu
        diu = rwork
c
        rwork = (1.0 - x1) * (z3 + 0.01)
        z1 = x1 * ztop(2) + rwork
        z2 = x1 * zbot(2) + rwork
c
        z12 = 0.5 * (z1 + z2)
c
        exphu = exp (0.5 * bvegu * (z1 - z2))
        expu  = exphu**2
c
        dispu = x1 * 0.7 * z1 + (1.0 - x1) * displ
c
cc 300  continue 
c
c mixing-length logarithms
c
cc      do 400 i = 1, npoi
c
        alogg  = alog (z0soi)
        alogi  = alog (z0sno)
        alogav = (1.0 - fi) * alogg + fi * alogi
c
c alog4 must be > z0soi, z0sno to avoid possible problems later 
c
        alog4 = alog ( max (z4, 1.1*z0soi, 1.1*z0sno) )
        alog3 = alog (z3-displ)
        alog2 = alog (z2-displ)
        alog1 = alog (z1-dispu)
        aloga = alog (za-dispu)
c
c initialize u2, alogu, alogl for first iteration's fstrat
c
        u2    = ua/exphu
        alogu = alog (max(.01, .1*(z1-z2)))
        alogl = alog (max(.01, .1*(z3-z4)))
c
cc  400 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine turcof (iter)
c ---------------------------------------------------------------------
c
c solves for wind speeds at various levels
c
c also computes upper and lower-region air-air transfer coefficients
c and saves them in com1d arrays cu and cl for use by turvap,
c and similarly for the solid-air transfer coefficients
c su, ss, sl, sg and si
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (input)
c
      integer iter          !current iteration number
c
c
c Local variables
c
cc      integer i             ! loop indice
c
      real xfac,            !
     >     x,               !
     >     rwork,           ! working variable
     >     cdmax,           ! max value for cd
     >     tauu,            !
     >     a,b,c,d,         !
     >     taul,            !
     >     ca,              ! to compute inverse air-air transfer coeffs
     >     cai, cbi, cci,   !
     >     cdi, cei, cfi,   !
     >     sg0,             ! to compute air-solid transfer coeff for soil
     >     si0              ! to compute air-solid transfer coeff for ice
c
      real yu, yl
c
c set stratification factors for lower and upper regions
c using values from the previous iteration
c
      xfac = 1.0
c
      call fstrat (t34, t12, xfac, q34, q12, z3, z2, 
     >             alogl, alogl, alog2, u2, richl, straml, strahl, iter)
c
      call fstrat (t12, ta,  tfac, q12, qa,  z1, za, 
     >             alogu, alogu, aloga, ua, richu, stramu, strahu, iter)
c
c eliminate c/d from eq (28), tau_l/rho from (26),(27), to get
c lower-story roughness alogl. yl/bdl is (tau_l/rho)/(c+d)
c
c equation numbers correspond to lsx description section 4.e
c
cc      do 100 i = 1, npoi
c
        x = ((alog4-alogav)/vonk)**2 * bdl
c
        rwork = 1. / expl
        yl = ((x+1)*expl + (x-1)*rwork)
     >        / ((x+1)*expl - (x-1)*rwork)
c
        alogl = alog3 - vonk * sqrt(yl/bdl)
c
cc 100  continue 
c
c eliminate tau_l/rho from (24),(25), tau_u/rho and a/b from
c (22),(23), to get upper-story roughness alogu
c 
c yu/bdu is (tau_u/rho)/(a+b)
c
cc      do 110 i = 1, npoi
c          
        x = ((alog2-alogl)/vonk)**2 * bdu / straml
c
        rwork = 1. / expu
        yu = ((x+1)*expu + (x-1)*rwork)
     >        / ((x+1)*expu - (x-1)*rwork)
c
        alogu = alog1 - vonk * sqrt(yu/bdu)
c
cc 110  continue
c
c define the maximum value of cd
c
      cdmax = 300.0 / (2.0 * dtime)
c
c get tauu (=tau_u/rho) from (21), a and b from (22),(23),
c taul (=tau_u/rho) from (25), c and d from (26),(27)
c
c changed the following to eliminate small errors associated with
c moving this code to single precision - affected c and d,
c which made u_ become undefined, as well as affecting some
c other variables
c
cc      do 200 i = 1, npoi
c
        tauu = (ua * vonk/(aloga-alogu))**2 * stramu
c
        a = 0.5 * tauu * (yu+1)/bdu
        b = 0.5 * tauu * (yu-1)/bdu
c
        taul = bdu * (a/expu - b*expu)
c
        c = 0.5 * taul * (yl+1)/bdl
        d = 0.5 * taul * (yl-1)/bdl
c
c evaluate wind speeds at various levels, keeping a minimum 
c wind speed of 0.01 m/s at all levels
c   
        u1  = max (0.01, sqrt (max (0.0, (a+b))))
        u12 = max (0.01, sqrt (max (0.0, (a/exphu+b*exphu))))
        u2  = max (0.01, sqrt (max (0.0, (a/expu +b*expu))))
        u3  = max (0.01, sqrt (max (0.0, (c+d))))
        u34 = max (0.01, sqrt (max (0.0, (c/exphl+d*exphl))))
        u4  = max (0.01, sqrt (max (0.0, (c/expl +d*expl))))
c
cc 200  continue
c
c compute inverse air-air transfer coeffs
c
c use of inverse individual coeffs cai, cbi, cci, cdi, cei, cfi avoids
c divide-by-zero as vegetation vanishes - combine into
c upper-region coeff cu from za to z12, and lower-region coeff
c cl from z34 to z12, and also coeffs
c
cc      do 300 i = 1, npoi
c
        ca = ua*strahu*vonk**2  /
     >       ((aloga-alogu) * (aloga-alog1))
c
        ca = min (cdmax, ca / (1. + ca * 1.0e-20))
c
        cai = 1.0 / (rhoa*ca)
c
        cbi = diu * (z1-z12) / (rhoa * 0.5*(u1+u12))
        cci = diu * (z12-z2) / (rhoa * 0.5*(u12+u2))
c
        cdi = (alog2-alogl) * (alog2-alog3) /
     >        (rhoa*u2*strahl*vonk**2)
c
        cei = dil * (z3-z34) / (rhoa * 0.5*(u3+u34))
        cfi = dil * (z34-z4) / (rhoa * 0.5*(u34+u4))
c
        cu = 1.0 / (cai + cbi)
        cl = 1.0 / (cci + cdi + cei)
c
c compute air-solid transfer coeffs for upper leaves, upper
c stems, lower story (su,ss,sl)
c
        su = rhoa * cleaf  * sqrt (u12 / dleaf(2))
        ss = rhoa * cstem  * sqrt (u12 / dstem(2))
        sl = rhoa * cgrass * sqrt (u34 / dleaf(1))
c
c compute air-solid transfer coeffs for soil and snow (sg,si)
c
c old technique
c
c       sg0 = rhoa * u4 * (vonk/(alog4-alogg))**2
c       si0 = rhoa * u4 * (vonk/(alog4-alogi))**2
c
c replace above formulations which depend on the log-wind profile
c (which may not work well below a canopy), with empirical formulation
c of Norman's. In the original LSX, turcof.f solves for the winds at
c the various levels from the momentum equations. This gives the transfer
c coefficients for heat and moisture. Heat and moisture eqns are then solved 
c in subroutine turvap. Using the empirical formulation of John Norman is 
c not consistent with the earlier solution for u4 (based on a logarithmic 
c profile just above the ground. However, this is used here because it 
c improved a lot simulations of the sensible heat flux over the 
c HAPEX-MOBILHY and FIFE sites
c
        sg0 = rhoa * (0.004 + 0.012 * u4)
        si0 = rhoa * (0.003 + 0.010 * u4)
c
c modify the cofficient to deal with cfi (see above)
c
        sg = 1.0 / (cfi + 1.0 / sg0)
        si = 1.0 / (cfi + 1.0 / si0)
c
cc  300 continue
c
c JAF:  not necessary 
c
c if no veg, recalculate coefficients appropriately for a
c single logarithmic profile, and 2 fictitious levels just
c above soil/snow surface. these levels are arbitrary but are
c taken as z2 and z4, preset in vegdat to a few cm height
c for bare ground and ice. use strahu from above, which used
c t12 and alogu (ok after first iteration)
c
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine turvap (iter, niter)
c ---------------------------------------------------------------------
c
c solves canopy system with linearized implicit sensible heat and
c moisture fluxes
c
c first, assembles matrix arr of coeffs in linearized equations
c for tu,ts,tl,t12,t34,q12,q34,tg,ti and assembles the right hand
c sides in the rhs vector
c
c then calls linsolve to solve this system, passing template mplate of
c zeros of arr 
c 
c finally calculates the implied fluxes and stores them 
c for the agcm, soil, snow models and budget calcs
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comhyd.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
c
c Arguments (input)
c
      integer niter,      ! total # of iteration
     >        iter        ! # of iteration
c
c local variables
c
      integer k
c
      real rwork, zwtot, rwork2, tgav, tiav, tuav, 
     >     tsav, tlav, quav, qsav, qlav, qgav, qiav, zwpud, zwsoi,
     >     psig, hfac, hfac2, zwopt, zwdry, betaw, emisoil, e, qs1,
     >     dqs1, xnumer, xdenom, betafac, betas
c
      real
     >  xu,
     >  xs,
     >  xl,
     >  chux, 
     >  chsx, 
     >  chlx,
     >  chgx,
     >  wlgx,
     >  wigx,
     >  fradu, 
     >  frads, 
     >  fradl,
     >  wu,
     >  ws,
     >  wl,      
     >  wg,
     >  wi, 
     >  qu,
     >  qs,  
     >  ql,  
     >  qg,  
     >  qi,
     >  dqu,
     >  dqs, 
     >  dql, 
     >  dqg, 
     >  dqi,
     >  tuold,
     >  tsold,
     >  tlold,
     >  tgold, 
     >  tiold,
     >  tupre,
     >  tspre,
     >  tlpre,
     >  tgpre,
     >  tipre,
     >  suw,
     >  ssw,
     >  slw,
     >  sut,
     >  slt,
     >  slt0,
     >  suh,
     >  ssh, 
     >  slh,
     >  cog,
     >  coi,
     >  zirg,
     >  ziri,
     >  qgfac,
     >  qgfac0
c
      save 
     >  xu,
     >  xs,
     >  xl,
     >  chux, 
     >  chsx, 
     >  chlx,
     >  chgx,
     >  wlgx,
     >  wigx,
     >  cog,
     >  coi,
     >  zirg,
     >  ziri,
     >  wu,  
     >  ws,    
     >  wl,    
     >  wg, 
     >  wi,
     >  tuold, 
     >  tsold, 
     >  tlold, 
     >  tgold, 
     >  tiold
c
      integer nqn
c
      parameter (nqn=9)
c
      real arr(nqn,nqn),      !    
     >     rhs(nqn),          ! right hand side
     >     vec(nqn)           ! 
c
      integer  mplate(nqn,nqn)
cc
c                  tu  ts  tl t12 t34 q12 q34  tg  ti 
c                  ----------------------------------
      data mplate / 1,  0,  0,  1,  0,  1,  0,  0,  0, !tu
     >              0,  1,  0,  1,  0,  1,  0,  0,  0, !ts
     >              0,  0,  1,  0,  1,  0,  1,  0,  0, !tl
     >              1,  1,  0,  1,  1,  0,  0,  0,  0, !t12
     >              0,  0,  1,  1,  1,  0,  0,  1,  1, !t34
     >              1,  1,  0,  0,  0,  1,  1,  0,  0, !q12
     >              0,  0,  1,  0,  0,  1,  1,  1,  1, !q34
     >              0,  0,  0,  0,  1,  0,  1,  1,  0, !tg
     >              0,  0,  0,  0,  1,  0,  1,  0,  1  !ti
     >            /
c
      include 'comsat.h'
c
c if first iteration, save original canopy temps in t*old
c (can use tsoi,tsno for original soil/snow skin temps), for
c rhs heat capacity terms in matrix soln, and for adjustment
c of canopy temps after each iteration
c
c also initialize soil/snow skin temps tg, ti to top-layer temps
c
c the variables t12, t34, q12, q34, for the first iteration
c are saved via global arrays from the previous gcm timestep,
c this is worth doing only if the agcm forcing is
c smoothly varying from timestep to timestep
c
      if (iter.eq.1) then
c
c weights for canopy coverages
c
cc      do 10 i = 1, npoi
          xu = 2.0 * lai(2) * fu
          xs = 2.0 * sai(2) * fu
          xl = 2.0 * (lai(1) + sai(1)) * fl * (1.0 - fi)
cc0     continue 
c
c specific heats per leaf/stem area
c
cc      do 20 i = 1, npoi
          chux = chu + ch2o * wliqu + cice * wsnou
          chsx = chs + ch2o * wliqs + cice * wsnos
          chlx = chl + ch2o * wliql + cice * wsnol
cc0     continue
c
cc      do 30 i = 1, npoi 
c
          rwork = poros(1) * rhow
c
          chgx = ch2o * wpud + cice * wipud
     >              + ((1.-poros(1))*csoi(1)*rhosoi(1)
     >              + rwork*(1.-wisoi(1))*wsoi(1)*ch2o
     >              + rwork*wisoi(1)*cice
     >                ) * hsoi(1)
c
          wlgx = wpud +
     >              rwork * (1. - wisoi(1)) *
     >              wsoi(1) * hsoi(1)
c
          wigx = wipud + rwork * wisoi(1) * hsoi(1)
c
cc0     continue 
c
c conductivity coeffs between ground skin and first layer
c
cc      do 40 i = 1, npoi
          cog = consoi(1) / (0.5 * hsoi(1))
          coi = consno      / (0.5 * max (hsno(1), hsnotop))
cc0     continue
c
c d(ir emitted) / dt for soil
c
        rwork = 4. * 0.95 * stef
c
cc      do 50 i = 1, npoi
          zirg = rwork * (tg**3)
          ziri = rwork * (ti**3)
cc0     continue
c
c updated temperature memory
c
cc      do 60 i = 1, npoi
          tuold = tu
          tsold = ts
          tlold = tl
          tgold = tg
          tiold = ti
cc0     continue
c
      endif
c
c set implicit/explicit factors w* (0 to 1) for this iteration
c w* is 1 for fully implicit, 0 for fully explicit
c for first iteration, impexp and impexp2 set w* to 1
c
      call impexp (wu, tu, chux, wliqu, wsnou, iter)
      call impexp (ws, ts, chsx, wliqs, wsnos, iter)
      call impexp (wl, tl, chlx, wliql, wsnol, iter)
      call impexp (wg, tg, chgx, wlgx,  wigx,  iter)
c
c call impexp2 for snow model
c
      call impexp2 (wi, ti, tiold, iter)

cc      write(100,*)"turvap 1"
c
c adjust t* for this iteration 
c
c in this routine we are free to choose them, 
c since they are just the central values about which the 
c equations are linearized - heat is conserved in the matrix
c solution because t*old are used for the rhs heat capacities
c
c here, let t* represent the previous soln if it was fully
c implicit, but weight towards t*old depending on the amount
c (1-w*) the previous soln was explicit
c
c this weighting is necessary for melting/freezing surfaces, for which t*
c is kept at t*old, presumably at or near tmelt
c
cc    do 80 i = 1, npoi
        tu = wu * tu + (1.0 - wu) * tuold
        ts = ws * ts + (1.0 - ws) * tsold
        tl = wl * tl + (1.0 - wl) * tlold
        tg = wg * tg + (1.0 - wg) * tgold
        ti = wi * ti + (1.0 - wi) * tiold
cc0   continue 
c
c save current "central" values for final flux calculations
c
cc      do 90 i = 1, npoi
        tupre = tu
        tspre = ts
        tlpre = tl
        tgpre = tg
        tipre = ti
cc 90   continue
c
c calculate various terms occurring in the linearized eqns,
c using values of t12, t34, q12, q34 from
c the previous iteration
c
c specific humidities for canopy and ground, and derivs wrt t
c for canopy
c
c limit derivs to avoid -ve implicit q's below,
c as long as d(temp)s in one iteration are le 10 deg k
c
cc      do 100 i = 1, npoi
c
        e      = esat(tu)
        qu  = qsat (e, psurf)
        dqu = dqsat (tu, qu)
        dqu = min (dqu, qu * 0.1)

cc      write(100,*)"turvap 1.4"

c
        e      = esat(ts)
        qs  = qsat (e, psurf)
        dqs = dqsat (ts, qs)
        dqs = min (dqs, qs * 0.1)
c
        e      = esat(tl)
        ql  = qsat (e, psurf)
        dql = dqsat (tl, ql)
        dql = min (dql, ql * 0.1)
c
        e      = esat(tg)
        qg  = qsat (e, psurf)
        dqg = dqsat (tg, qg)
        dqg = min (dqg, qg * 0.1)
c
        e      = esat(ti)
        qi  = qsat (e, psurf)
        dqi = dqsat (ti, qi)
        dqi = min (dqi, qi * 0.1)

cc      write(100,*)"turvap 1.5"
c
cc 100  continue
c
c set qgfac0, factor by which soil surface specific humidity
c is less than saturation
c
c it is important to note that the qgfac expression should
c satisfy timestep cfl criterion for upper-layer soil moisture
c for small wsoi(1)
c
c for each iteration, qgfac is set to qgfac0, or to 1 if
c condensation onto soil is anticipated (loop 110 in canopy.f)
c
c Evaporation from bare soil is calculated using the "beta method"
c (e.g., eqns 5 & 7 of Mahfouf and Noilhan 1991, JAM 30 1354-1365),
c but converted to the "alpha method" (eqns 2 & 3 of M&N), to match
c the structure in IBIS. The conversion from the beta to alpha
c method is through the relationship:
c   alpha * qgs - q34 = beta * (hfac * qgs - q34),
c from which one solves for alpha (which is equal to qgfac0):
c   qgfac0 = alpha = (beta * hfac) + (1 - beta)*(q34/qgs)
c
cc        do 105 i = 1, npoi
c
c first calculate the total saturated fraction at the soil surface
c (including puddles ... see soil.f)
c
          zwpud = max (0.0, min (0.5, 0.5*(wpud+wipud)/wpudmax) )
          zwsoi = (1.0 - wisoi(1)) * wsoi(1) + wisoi(1)
          zwtot = zwpud + (1. - zwpud) * zwsoi
c
c next calculate the matric potential (from eqn 9.3 of Campbell and
c Norman), multiply by gravitational acceleration to get in units
c of J/kg, and calculate the relative humidity at the soil water
c surface (i.e., within the soil matrix), based on thermodynamic
c theory (eqn 4.13 of C&N)
c
	    zwtot = max(zwtot, 0.01)

          psig = -grav * suction(1) * (zwtot ** (-bex(1)))

cc         write(100,*)"turvap 1.6", bex, zwtot, psig, rvap, tg

          hfac = exp(psig/(rvap*tg))

cc      write(100,*)"turvap 1.7"

c
c then calculate the relative humidity of the air (relative to
c saturation at the soil temperature). Note that if hfac2 > 1
c (which would imply condensation), then qgfac is set to 1
c later in the code (to allow condensation to proceed at the
c "potential rate")
c
          hfac2 = q34/qg
c
c set the "beta" factor and then calculate "alpha" (i.e., qgfac0)
c as the beta-weighted average of the soil water RH and the "air RH"
c First calculate beta_w:
c
          zwopt = 1.0
          zwdry = swilt(1)
          betaw = max(0.0, min(1., (zwtot - zwdry)/(zwopt - zwdry)) )
cc      write(100,*)"turvap 2"
c
c Next convert beta_w to beta_s (see Milly 1992, JClim 5 209-226):
c
          emisoil = 0.95
          e      = esat(t34)
          qs1    = qsat (e, psurf)
          dqs1   = dqsat (t34, qs1)
          xnumer = hvap * dqs1
          xdenom = cp + (4.0 * emisoil * stef * (t34)**3) / sg
          betafac = xnumer / xdenom
          betas = betaw / (1.0 + betafac * (1.0 - betaw))
c
c Combine hfac and hfac2 into qgfac0 ("alpha") using beta_s
c
          qgfac0 = betas * hfac + (1. - betas) * hfac2
cc  105   continue
c
c set fractions covered by intercepted h2o to 1 if dew forms
c
c these fwet*x are used only in turvap, and are distinct from
c the real fractions fwet* that are set in fwetcal
c
c they must be exactly 1 if q12 > qu or q34 > ql, to zero transpiration
c by the factor 1-fwet[u,l]x below, so preventing "-ve" transp
c
c similarly, set qgfac, allowing for anticipated dew formation
c to avoid excessive dew formation (which then infiltrates) onto
c dry soils
c
cc      do 110 i = 1, npoi
c
        fwetux = fwetu
        if (q12.gt.qu) fwetux = 1.0
c
        fwetsx = fwets
        if (q12.gt.qs) fwetsx = 1.0
c
        fwetlx = fwetl
        if (q34.gt.ql) fwetlx = 1.0
c
        qgfac = qgfac0
        if (q34.gt.qg) qgfac = 1.0
c
c set net absorbed radiative fluxes for canopy components
c
        fradu = 0.0
c
        if (lai(2).gt.epsilon)
     >     fradu = (solu + firu) / (2.0 * lai(2))
c
        frads = 0.0
c
        if (sai(2).gt.epsilon)
     >     frads = (sols + firs) / (2.0 * sai(2))
c
        fradl = 0.0
c
        if ((lai(1)+sai(1)).gt.epsilon)
     >     fradl = (soll + firl) /
     >                (2.0 * (lai(1) + sai(1)))
cc      write(100,*)"turvap 3"
c
cc 110  continue
c
c calculate canopy-air moisture transfer coeffs for wetted
c leaf/stem areas, and for dry (transpiring) leaf areas
c
c the wetted-area coeffs suw,ssw,slw are constrained to be less
c than what would evaporate 0.8 * the intercepted h2o mass in 
c this timestep (using previous iteration's q* values)
c
c this should virtually eliminate evaporation-overshoots and the need
c for the "negative intercepted h2o"  correction in steph2o2
c        
cc      do 200 i = 1, npoi
c
c coefficient for evaporation from wet surfaces in the upper canopy:
c
        suw = min ( fwetux * su, 
     >                 0.8 * (wliqu + wsnou) /
     >                 max (dtime * (qu - q12), epsilon))
c
c coefficient for transpiration from average upper canopy leaves:
c
        sut = (1.0 - fwetux) * 0.5 *
     >           ( totcondub * frac(1) +
     >             totcondub * frac(2) +
     >             totcondub * frac(3) +
     >             totconduc * frac(4) +
     >             totcondub * frac(5) +
     >             totconduc * frac(6) +
     >             totcondub * frac(7) +
     >             totconduc * frac(8) )	 !xjz 	totcondub=>totconduc
c
        sut = max (0.0, sut)
c
c coefficient for sensible heat flux from upper canopy:
c
        suh = suw * (rliqu  * hvapf(tu,ta)  +
     >                 (1.-rliqu) * hsubf(tu,ta)) +
     >           sut *              hvapf(tu,ta)
c
c coefficient for evaporation from wet surfaces on the stems:
c
        ssw = min (fwetsx * ss, 
     >                0.8 * (wliqs + wsnos)
     >                / max (dtime * (qs - q12), epsilon))
c
c coefficient for sensible heat flux from stems:
c
        ssh = ssw * (rliqs  * hvapf(ts,ta) +
     >                 (1.-rliqs) * hsubf(ts,ta))
c
c coefficient for evaporation from wet surfaces in the lower canopy:
c
        slw = min (fwetlx * sl, 
     >                0.8 * (wliql + wsnol)
     >                / max (dtime * (ql - q34), epsilon))
c
c coefficient for transpiration from average lower canopy leaves:
c
        slt0 = (1. - fwetlx) * 0.5 *
     >            ( totcondls * frac(9)  +
     >              totcondls * frac(10) +
     >              totcondl4 * frac(11) +
     >              totcondl3 * frac(12) )
c
        slt0 = max (0., slt0)
c
c averaged over stems and lower canopy leaves:
c 
        slt = slt0 * lai(1) / max (lai(1)+sai(1), epsilon)
c
c coefficient for sensible heat flux from lower canopy:
c
        slh = slw * (  rliql  * hvapf(tl,ta)  +
     >                   (1.-rliql) * hsubf(tl,ta)) +
     >           slt *                hvapf(tl,ta)
c
cc 200  continue
c
c set the matrix of coefficients and the right-hand sides
c of the linearized equations
c
      call const(arr, nqn*nqn, 0.0)
      call const(rhs, nqn, 0.0)
c

cc      write(100,*)"turvap 4"
      rwork = 1. / dtime
c
c upper leaf temperature tu
c
cc      do 300 i = 1, npoi
c
        rwork2 = su*cp
        arr(1,1) = chux*rwork
     >             + wu*rwork2
     >             + wu*suh*dqu
        arr(1,4) = -rwork2
        arr(1,6) = -suh
        rhs(1) = tuold*chux*rwork
     >           - (1.-wu)*rwork2*tu
     >           - suh * (qu-wu*dqu*tu)
     >           + fradu - pfluxu
c 
cc 300  continue
c
c upper stem temperature ts
c
cc      do 310 i = 1, npoi
c
        rwork2 = ss*cp
        arr(2,2) = chsx*rwork
     >             + ws*rwork2
     >             + ws*ssh*dqs
        arr(2,4) = -rwork2
        arr(2,6) = -ssh
        rhs(2) = tsold*chsx*rwork
     >           - (1.-ws)*rwork2*ts
     >           - ssh * (qs-ws*dqs*ts)
     >           + frads - pfluxs
c
cc 310  continue
c
c lower veg temperature tl
c
cc      do 320 i = 1, npoi
c
        rwork2 = sl*cp
        arr(3,3) = chlx*rwork
     >             + wl*rwork2
     >             + wl*slh*dql
        arr(3,5) = -rwork2
        arr(3,7) = -slh
        rhs(3) = tlold*chlx*rwork
     >           - (1.-wl)*rwork2*tl
     >           - slh * (ql-wl*dql*tl)
     >           + fradl - pfluxl
c
cc 320  continue
c
c upper air temperature t12
c
cc      do 330 i = 1, npoi
c
        rwork = xu*su
        rwork2 = xs*ss
        arr(4,1) = -wu*rwork
        arr(4,2) = -ws*rwork2
        arr(4,4) = cu + cl + rwork + rwork2
        arr(4,5) = -cl
        rhs(4) = cu*ta*tfac
     >           + (1.-wu)*rwork*tu
     >           + (1.-ws)*rwork2*ts
cc      write(100,*)"turvap 5"
c
cc 330  continue
c
c lower air temperature t34
c
cc      do 340 i = 1, npoi
c
        rwork = xl*sl
        rwork2 = fi*si
        arr(5,3) = -wl*rwork
        arr(5,4) = -cl
        arr(5,5) = cl + rwork
     >             + (1.-fi)*sg + rwork2
        arr(5,8) = -wg*(1.-fi)*sg
        arr(5,9) = -wi*rwork2
        rhs(5) = (1.-wl)*rwork           *tl
     >           + (1.-wg)*(1.-fi)*sg*tg
     >           + (1.-wi)*rwork2          *ti
c
cc 340  continue
c
c upper air specific humidity q12
c
cc      do 350 i = 1, npoi
c
        rwork = xu*(suw+sut)
        rwork2 = xs*ssw
        arr(6,1) = -wu*rwork *dqu
        arr(6,2) = -ws*rwork2*dqs
        arr(6,6) = cu + cl
     >             + rwork + rwork2
        arr(6,7) = -cl
        rhs(6) = cu*qa
     >           + rwork  * (qu-wu*dqu*tu)
     >           + rwork2 * (qs-ws*dqs*ts)
c
cc  350 continue
c
c lower air specific humidity q34
c
cc      do 360 i = 1, npoi
c
        rwork  = xl*(slw+slt)
        rwork2 = (1.-fi)*sg
        arr(7,3) = -wl*rwork*dql
        arr(7,6) = -cl
        arr(7,7) = cl + rwork
     >             + rwork2 +fi*si
        arr(7,8) = -wg*rwork2*qgfac*dqg
        arr(7,9) = -wi*fi*si*dqi
        rhs(7)= rwork           *(ql-wl*dql*tl)
     >          + rwork2*qgfac *(qg-wg*dqg*tg)
     >          + fi *si    *(qi-wi*dqi*ti)

cc      write(100,*)"turvap 6"
c
cc  360 continue
c
c soil skin temperature
c
c (there is no wg in this eqn since it solves for a fully
c implicit tg. wg can be thought of as the fractional soil
c area using a fully implicit soln, and 1-wg as that using a
c fully explicit soln. the combined soil temperature is felt
c by the lower air, so wg occurs in the t34,q34 eqns above.)
c
cc      do 370 i = 1, npoi
c
        rwork  = sg*cp
        rwork2 = sg*hvasug
        arr(8,5) = -rwork
        arr(8,7) = -rwork2
        arr(8,8) = rwork + rwork2*qgfac*dqg
     >             + cog + zirg
        rhs(8) = -rwork2*qgfac*(qg-dqg*tg)
     >           + cog*tsoi(1)
     >           + solg + firg + zirg * tgold
c
cc  370 continue
c
c snow skin temperature
c
c (there is no wi here, for the same reason as for wg above.)
c
cc      do 380 i = 1, npoi
c
	rwork  = si*cp
        rwork2 = si*hvasui
        arr(9,5) = -rwork
        arr(9,7) = -rwork2
        arr(9,9) = rwork + rwork2*dqi
     >             + coi + ziri
        rhs(9) = -rwork2*(qi-dqi*ti)
     >           + coi*tsno(1)
     >           + soli + firi + ziri * tiold
c
cc  380 continue
c
c solve the systems of equations
c
      call linsolve (arr, rhs, vec, mplate, nqn)
c
c copy this iteration's solution to t*, q12, q34
c
cc      do 400 i = 1, npoi
c
        tu  = vec(1)
        ts  = vec(2)
        tl  = vec(3)
        t12 = vec(4)
        t34 = vec(5)
        tg  = vec(8)
        ti  = vec(9)
c
        q12 = vec(6)
        q34 = vec(7)

cc      write(100,*)"turvap 7"
c
cc  400 continue
c
c all done except for final flux calculations,
c so loop back for the next iteration (except the last)
c
      if (iter.lt.niter) return
c
c evaluate sensible heat and moisture fluxes (per unit
c leaf/stem/snow-free/snow-covered area as appropriate)
c
c *******************************
c diagnostic sensible heat fluxes
c *******************************
c
cc      do 500 i = 1, npoi
c
        fsena = cp * cu * (ta*tfac - t12)
c
        tgav = wg*tg + (1.-wg)*tgpre
        fseng = cp * sg * (tgav - t34)
c
        tiav = wi*ti + (1.-wi)*tipre
        fseni = cp * si * (tiav - t34)

        tuav = wu*tu + (1. - wu)*tupre
        fsenu = cp * su * (tuav - t12)
c
        tsav = ws*ts + (1. - ws)*tspre
        fsens = cp * ss * (tsav - t12)
c
        tlav = wl*tl + (1. - wl)*tlpre
        fsenl = cp * sl * (tlav - t12)
c
cc 500  continue
c
c *************************
c calculate moisture fluxes
c *************************
c
cc      do 510 i = 1, npoi
c
c total evapotranspiration from the entire column
c
        fvapa  = cu * (qa-q12)
c
c evaporation from wet surfaces in the upper canopy
c and transpiration per unit leaf area - upper canopy
c
        quav = qu + wu*dqu*(tu-tupre)
        fvapuw = suw * (quav-q12)
        fvaput = max (0.0, sut * (quav-q12))
c
c evaporation from wet surfaces on stems
c
        qsav = qs + ws*dqs*(ts-tspre)
        fvaps = ssw * (qsav-q12)

cc      write(100,*)"turvap 8"
c
c evaporation from wet surfaces in the lower canopy
c and transpiration per unit leaf area - lower canopy
c
        qlav = ql + wl*dql*(tl-tlpre)
        fvaplw = slw  * (qlav-q34)
        fvaplt = max (0.0, slt0 * (qlav-q34))
c
c evaporation from the ground
c
        qgav = qg + wg*dqg*(tg-tgpre)
        fvapg = sg * (qgfac*qgav - q34)
c
c evaporation from the snow
c
        qiav = qi + wi*dqi*(ti-tipre)
        fvapi = si * (qiav-q34)
c
cc 510  continue
c 
c adjust ir fluxes
c
cc      do 520 i = 1, npoi
c
        firg = firg - wg*zirg*(tg - tgold)
        firi = firi - wi*ziri*(ti - tiold)
        firb = firb + (1.-fi)*wg*zirg*(tg-tgold)
     >                    +     fi *wi*ziri*(ti-tiold)
c
c impose constraint on skin temperature
c
        ti = min (ti, tmelt)
c
cc 520  continue
c
c set upsoi[u,l], the actual soil water uptake rates from each
c soil layer due to transpiration in the upper and lower stories,
c for the soil model 
c
      do 600 k = 1, nsoilay
cc        do 610 i = 1, npoi
c
          upsoiu(k) = fvaput * 2.0 * lai(2) * fu *
     >                  stressu(k) / max (stresstu, epsilon)
c
          upsoil(k) = fvaplt * 2.0 * lai(1) * fl *
     >                  (1. - fi) *
     >                  stressl(k) / max (stresstl, epsilon)
c
cc 610    continue
 600  continue
c
c set net evaporation from intercepted water, net evaporation
c from the surface, and net transpiration rates
c
cc      do 700 i = 1, npoi
c
c evaporation from intercepted water
c
        ginvap = fvapuw * 2.0 * lai(2) * fu +
     >              fvaps  * 2.0 * sai(2) * fu +
     >              fvaplw * 2.0 * (lai(1) + sai(1)) * 
     >                                 fl * (1. - fi)
c
c evaporation from soil and snow surfaces
c
        gsuvap = fvapg  * (1. - fi) + fvapi  * fi
c
c transpiration
c
        gtrans = fvaput * 2.0 * lai(2) * fu +
     >              fvaplt * 2.0 * lai(1) * fl * (1.-fi)
c
        gtransu = fvaput * 2.0 * lai(2) * fu
        gtransl = fvaplt * 2.0 * lai(1) * fl * (1.-fi)
c
cc 700  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine fstrat (tb, tt, ttfac, qb, qt, zb, zt, 
     >                   albm, albh, alt, u, rich, stram, strah, iter)
c ---------------------------------------------------------------------
c
c computes mixing-length stratification correction factors
c for momentum and heat/vapor, for current 1d strip, using
c parameterizations in louis (1979),blm,17,187. first computes
c richardson numbers. sets an upper limit to richardson numbers
c so lower-veg winds don't become vanishingly small in very
c stable conditions (cf, carson and richards,1978,blm,14,68)
c
c system  is as in louis(1979). system (vi) is improved as
c described in louis(1982), ecmwf workshop on planetary boundary
c layer parameterizations,november 1981,59-79 (qc880.4 b65w619)
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
c
c input variables
c
      integer iter   ! current iteration number
c
      real ttfac     ! pot. temp factor for ttop (relative to bottom,supplied)
c
      real
     > tb,     ! bottom temperature (supplied)
     > tt,     ! top temperature (supplied)
     > qb,     ! bottom specific humidity (supplied)
     > qt,     ! top specific humidity (supplied)
     > zb,     ! height of bottom (supplied)
     > zt,     ! height of top (supplied)
     > albm,   ! log (bottom roughness length) for momentum (supplied)
     > albh,   ! log (bottom roughness length) for heat/h2o (supplied)
     > alt,    ! log (z at top) (supplied)
     > u,      ! wind speed at top (supplied)
     > rich,   ! richardson number (returned)
     > stram,  ! stratification factor for momentum (returned)
     > strah,  ! stratification factor for heat/vap (returned)
     > stramx, !
     > strahx  !
c
c local variables
c
      integer 
     > indp,   !
     > indq    !
c
      integer np, nq
c
      real zht, zhb, xm, xh, rwork, ym, yh, z, w
c ---------------------------------------------------------------------
      np = 0
      nq = 0
c
c do for all points
c
cc      do 100 i = 1, npoi
c
c calculate richardson numbers
c
        zht = tt*ttfac*(1.+.622*qt)
        zhb = tb*      (1.+.622*qb)
c
        rich = grav * max (zt-zb, 0.)
     >            * (zht-zhb) / (0.5*(zht+zhb) * u**2)
c
c bound richardson number between -2.0 (unstable) to 1.0 (stable)
c
        rich = max (-2.0, min (rich, 1.0))
c
cc 100  continue
c
c set up indices for points with negative or positive ri
c
cc      do 110 i = 1, npoi
c
        if (rich.le.0.) then
          np = np + 1
          indp = 1
        else
          nq = nq + 1
          indq = 1
        endif
c
cc  110 continue
c
c calculate momentum and heat/vapor factors for negative ri
c
      if (np.gt.0) then
c
cc        do 200 j = 1, np
c
cc          i = indp(j)
c
          xm = max (alt-albm, .5)
          xh = max (alt-albh, .5)
c
          rwork = sqrt(-rich)
c
          ym = (vonk/xm)**2 * exp (0.5*xm) * rwork
          yh = (vonk/xh)**2 * exp (0.5*xh) * rwork
c
c system (vi)
c
          stramx =   1.0 - 2*5*rich / (1.0 + 75*ym)
          strahx =   1.0 - 3*5*rich / (1.0 + 75*yh)
c
cc  200   continue
c
      endif
c
c calculate momentum and heat/vapor factors for positive ri
c
      if (nq.gt.0) then
c
cc        do 300 j=1,nq
c
cc          i = indq(j)
c
c system (vi)
c
          z = sqrt(1.0 + 5 * rich)
c
          stramx = 1.0 / (1.0 + 2*5*rich / z)
          strahx = 1.0 / (1.0 + 3*5*rich * z)
c
cc  300   continue
c
      endif
c
c except for the first iteration, weight results with the
c previous iteration's values. this improves convergence by
c avoiding flip-flop between stable/unstable stratif, eg,
c with cold upper air and the lower surface being heated by
c solar radiation
c
      if (iter.eq.1) then
c
cc        do 400 i = 1, npoi
c
          stram = stramx
          strah = strahx
c
cc  400   continue
c
      else
c
        w = 0.5
c
cc        do 410 i = 1, npoi
c
          stram = w * stramx + (1.0 - w) * stram
          strah = w * strahx + (1.0 - w) * strah
c
cc  410   continue
c
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine impexp (wimp, tveg, ch, wliq, wsno, iter)
c ---------------------------------------------------------------------
c
c sets the implicit vs explicit fraction in turvap calcs for
c upper leaves, upper stems or lower veg. this is to account for
c temperatures of freezing/melting intercepted h2o constrained
c at the melt point. if a purely implicit calc is used for such
c a surface, the predicted temperature would be nearly the atmos
c equil temp with little sensible heat input, so the amount of
c freezing or melting is underestimated. however, if a purely
c explicit calc is used with only a small amount of intercepted
c h2o, the heat exchange can melt/freeze all the h2o and cause
c an unrealistic huge change in the veg temp. the algorithm
c below attempts to avoid both pitfalls
c
c common blocks
c
      include 'implicit.h'
      include 'compar.h'
c
c input/output variables
c
      integer iter  ! current iteration number (supplied)
c
      real      
     >  wimp, ! implicit/explicit fraction (0 to 1) (returned)
     >  tveg, ! temperature of veg (previous iteration's soln) (supp)
     >  ch,   ! heat capacity of veg (supplied)
     >  wliq, ! veg intercepted liquid (supplied)
     >  wsno  ! veg intercepted snow (supplied)
c
c local variables
c
cc      integer i
c
      real h, z, winew
c
c for first iteration, set wimp to fully implicit, and return
c
      if (iter.eq.1) then
        wimp = 1.0

        return
      endif
c
c for second and subsequent iterations, estimate wimp based on
c the previous iterations's wimp and its resulting tveg.
c
c calculate h, the "overshoot" heat available to melt any snow
c or freeze any liquid. then the explicit fraction is taken to
c be the ratio of h to the existing h2o's latent heat (ie, 100%
c explicit calculation if not all of the h2o would be melted or
c frozen). so winew, the implicit amount, is 1 - that ratio.
c but since we are using the previous iteration's t* results
c for the next iteration, to ensure convergence we need to damp
c the returned estimate wimp by averaging winew with the 
c previous estimate. this works reasonably well even with a
c small number of iterations (3), since for instance with large
c amounts of h2o so that wimp should be 0., a good amount of 
c h2o is melted or frozen with wimp = .25
c
cc      do 100 i = 1, npoi
c
        h = ch * (tveg - tmelt)
        z = max (abs(h), epsilon)
c
        winew = 1.0
c
        if (h.gt.epsilon)  winew = 1. - min (1., hfus * wsno / z)
        if (h.lt.-epsilon) winew = 1. - min (1., hfus * wliq / z)
c
        wimp = 0.5 * (wimp + winew)
c
cc  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine impexp2 (wimp, t, told, iter)
c ---------------------------------------------------------------------
c
c sets the implicit vs explicit fraction in turvap calcs for
c seaice or snow skin temperatures, to account for temperatures
c of freezing/melting surfaces being constrained at the melt
c point
c
c unlike impexp, don't have to allow for all h2o 
c vanishing within the timestep
c
c wimp   = implicit fraction (0 to 1) (returned)
c
      include 'implicit.h'
      include 'compar.h'
c
c input variables
c
      integer iter
      real 
     >  wimp, t, told
c
c local variables
c
cc      integer i    ! loop indice
c
c for first iteration, set wimp to fully implicit, and return
c
      if (iter.eq.1) then
        wimp = 1.0
        return
      endif
c
cc      do 100 i = 1, npoi
c
        if ((t-told).gt.epsilon) wimp = (tmelt - told) / 
     >                                           (t  - told)
        wimp = max (0.0, min (1.0, wimp))
c
cc 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine fwetcal
c ---------------------------------------------------------------------
c
c calculates fwet[u,s,l], the fractional areas wetted by 
c intercepted h2o (liquid and snow combined) -  the maximum value
c fmax (<1) allows some transpiration even in soaked conditions
c
c use a linear relation between fwet* and wliq*,wsno* (at least
c for small values), so that the implied "thickness" is constant
c (equal to wliq*max, wsno*max as below) and the typical amount
c evaporated in one timestep in steph2o will not make wliq*,wsno*
c negative and thus cause a spurious unrecoverable h2o loss
c
c (the max(w*max,.01) below numericaly allows w*max = 0 without
c blowup.) in fact evaporation in one timestep *does* sometimes
c exceed wliq*max (currently 1 kg/m2), so there is an additional
c safeguard in turvap that limits the wetted-area aerodynamic
c coefficients suw,ssw,slw -- if that too fails, there is an 
c ad-hoc adjustment in steph2o2 to reset negative wliq*,wsno*
c amounts to zero by taking some water vapor from the atmosphere.
c
c also sets rliq[u,s,l], the proportion of fwet[u,s,l] due to
c liquid alone. fwet,rliq are used in turvap, rliq in steph2o. 
c (so rliq*fwet, (1-rliq)*fwet are the fractional areas wetted
c by liquid and snow individually.) if fwet is 0, choose rliq
c = 1 if t[u,s,l] ge tmelt or 0 otherwize, for use by turvap and
c steph2o in case of initial dew formation on dry surface.
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
cc      integer i           ! loop indice
c
      real fmax,          ! maximum water cover on two-sided leaf
     >     xliq,          ! fraction of wetted leaf (liquid only)
     >     xtot           ! fraction of wetted leaf (liquid and snow)
c
c maximum water cover on two-sided leaf
c
      parameter (fmax = 0.25)
c
c upper leaves
c
cc      do 100 i = 1, npoi
c
        xliq = wliqu / max (wliqumax, 0.01)
        xtot = xliq + wsnou / max (wsnoumax, 0.01)
c
        fwetu = min (fmax, xtot)
        rliqu = xliq / max (xtot, epsilon)
c
        if (fwetu.eq.0.0) then
          rliqu = 1.0
          if (tu.lt.tmelt) rliqu = 0.0
        endif
c
cc  100 continue
c
c upper stems
c
cc      do 200 i = 1, npoi
c
        xliq = wliqs / max (wliqsmax, 0.01)
        xtot = xliq + wsnos / max (wsnosmax, 0.01)
c
        fwets = min (fmax, xtot)
        rliqs = xliq / max (xtot, epsilon)
c
        if (fwets.eq.0.0) then
          rliqs = 1.0
          if (ts.lt.tmelt) rliqs = 0.0
        endif
c
cc  200 continue
c
c lower veg
c
cc      do 300 i = 1, npoi
c
        xliq = wliql / max (wliqlmax, 0.01)
        xtot = xliq + wsnol / max (wsnolmax, 0.01)
c
        fwetl = min (fmax, xtot)
        rliql = xliq / max (xtot, epsilon)
c
        if (fwetl.eq.0.) then
          rliql = 1.0
          if (tl.lt.tmelt) rliql = 0.0
        endif
c
cc  300 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine cascade
c ---------------------------------------------------------------------
c
c steps intercepted h2o due to drip, precip, and min/max limits
c
c calls steph2o for upper leaves, upper stems and lower veg in
c iurn, adjusting precips at each level
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
cc      integer i            ! loop indice
c
      real twet3           ! Function: wet bulb temperature (K)
      real twetbulb        ! wet bulb temperature (K)
c    
      real
     >  xai,         !lai and/or sai for veg component
                           ! (allows steph2o to work on any veg component)
     >  rain,        !rainfall at appropriate level (modified by steph2o)
     >  train,       !temperature of rain (modified by steph2o)  
     >  snow,        !snowfall at appropriate level (modified by steph2o)
     >  tsnow,       !temperature of snow (modified by steph2o)
     >  x1,        ! 
     >  x2,        ! 
     >  x3,        ! 
     >  x4         ! 
c
c adjust rainfall and snowfall rates at above-tree level
c
c set wliqmin, wsnomin -- unlike wliq*max, wsno*max, these are
c part of the lsx numerical method and not from the vegetation
c database, and they are the same for all veg components
c
c the value 0.0010 should be small compared to typical precip rates
c times dtime to allow any intercepted h2o to be initiated, but
c not too small to allow evap rates to reduce wliq*, wsno* to
c that value in a reasonable number of time steps
c
      wliqmin = 0.0010 * (dtime/3600.) * (wliqumax / 0.2)
      wsnomin = 0.0010 * (dtime/3600.) * (wsnoumax / 2.0)
c
cc      do 50 i=1,npoi
        rainu = raina
c
c set rain temperature to the wet bulb temperature
c
        if (ta .gt. tmelt) then
           twetbulb = twet3( ta, qa, psurf )
        else
           twetbulb = tmelt
        endif
        trainu = max (twetbulb, tmelt)
        x1 = 0.0
        x2 = max (t12, tmelt)
cc   50 continue
c
      call mix (rainu,trainu, rainu,trainu, x1,x2, vzero,vzero)
c
cc      do 52 i=1,npoi
        snowu = snowa
        tsnowu = min (ta, tmelt)
        x1 = 0.0
        x2 = min (t12, tmelt)
cc   52 continue
c
      call mix (snowu,tsnowu, snowu,tsnowu, x1,x2, vzero,vzero)
c
c set up for upper leaves
c
cc      do 100 i = 1, npoi
        xai   = 2.0 * lai(2)
        rain  = rainu
        train = trainu
        snow  = snowu
        tsnow = tsnowu
cc  100 continue
c
c step upper leaves
c
      call steph2o
     >  (tu,  wliqu,  wsnou,  xai,  pfluxu,  rain, train, snow, tsnow,
     >   tdripu, tblowu, wliqumax, wsnoumax, wliqmin, wsnomin)
c
c set up for upper stems
c the upper stems get precip as modified by the upper leaves
c
cc      do 200 i=1,npoi
        xai = 2.0 * sai(2)
cc  200 continue
c
c step upper stems
c
      call steph2o
     >  (ts,  wliqs,  wsnos,  xai,  pfluxs,  rain, train, snow, tsnow,
     >   tdrips, tblows, wliqsmax, wsnosmax, wliqmin, wsnomin)
c
c adjust rainfall and snowfall rates at below-tree level
c allowing for upper-veg interception/drip/belowoff
c
cc      do 300 i=1,npoi
        x1 = fu*rain
        x2 = (1.-fu)*rainu
        x3 = 0.0
        x4 = max (t34, tmelt)
cc  300 continue
c
      call mix (rainl,trainl, x1,train, x2,trainu, x3,x4)
c
cc      do 310 i=1,npoi
        x1 = fu*snow
        x2 = (1.-fu)*snowu
        x3 = 0.0
        x4 = min (t34, tmelt)
cc  310 continue
c
      call mix (snowl,tsnowl, x1,tsnow, x2,tsnowu, x3,x4)
c
c set up for lower veg
c
cc      do 400 i = 1, npoi
        xai   = 2.0 * (lai(1) + sai(1))
        rain  = rainl
        train = trainl
        snow  = snowl
        tsnow = tsnowl
cc  400 continue
c
c step lower veg
c
      call steph2o
     >  (tl,  wliql,  wsnol,  xai,  pfluxl,  rain, train, snow, tsnow,
     >   tdripl, tblowl, wliqlmax, wsnolmax, wliqmin, wsnomin)
c
c adjust rainfall and  snowfall rates at soil level,
c allowing for lower-veg interception/drip/blowoff
c
cc      do 500 i=1,npoi
        x1 = fl * rain
        x2 = (1.-fl) * rainl
cc  500 continue
c
      call mix (raing,traing, x1,train, x2,trainl, vzero,vzero)
c
cc      do 510 i=1,npoi
        x1 = fl * snow
        x2 = (1.-fl) * snowl
cc  510 continue
c
      call mix (snowg,tsnowg, x1,tsnow, x2,tsnowl, vzero,vzero)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine steph2o
     >  (tveg,  wliq,  wsno,  xai,  pflux,  rain, train, snow, tsnow,
     >   tdrip, tblow, wliqmax, wsnomax, wliqmin, wsnomin)
c ---------------------------------------------------------------------
c
c steps intercepted h2o for one canopy component (upper leaves, 
c upper stems, or lower veg) through one lsx time step, adjusting
c for h2o sensible heat and phase changes. also modifies precip
c due to interception and drip,blowoff
c
c 
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (all arguments are supplied (unchanged) unless otherwise noted
c
      real tdrip,       ! e-folding time of liquid drip  tdrip[u,s,l]
     >     tblow,       ! e-folding time of snow blowoff tblow[u,s,l]
     >     wliqmax,     ! max amount of intercepted liquid wliq[u,s,l]max
     >     wsnomax,     ! max amount of intercepted snow   wsno[u,s,l]max
     >     wliqmin,     ! min amount of intercepted liquid (same name for u,s,l)
     >     wsnomin      ! min amount of intercepted snow (same name for u,s,l)
c
      real       
     >  tveg,     ! temperature of veg component t[u,s,l]
     >  wliq,     ! intercepted liquid amount wliq[u,s,l] (returned)
     >  wsno,     ! intercepted snow amount wsno[u,s,l] (returned)
     >  xai,      ! lai, sai, lai+sai for upper leaves/stems,lower veg
     >  pflux,    ! ht flux due to adjust of intercep precip (returned)
     >  rain,     ! rainfall rate. Input: above veg, Output: below veg
     >  train,    ! temperature of rain. (returned)
     >  snow,     ! snowfall rate. Input: above veg, output: below veg
     >  tsnow     ! temperature of snow (returned)
c
c local variables:
c
cc      integer i         ! loop indice
c
      real rwork,       ! 1/dtime
     >     x,           ! work variable
     >     rwork2,      ! work variable: ch2o - cice
     >     dw           ! correction: freezing liguid or melting snow
c
      real fint,  ! precip fraction intercepted by unit leaf/stem area
     >     drip,  ! rate of liquid drip
     >     blow   ! rate of snow blowoff
c
c ---------------------------------------------------------------------
c
c calculate fint, the intercepted precip fraction per unit
c leaf/stem area -- note 0.5 * lai or sai (similar to irrad)
c 
cc      do 50 i = 1, npoi
c
        if (xai.ge.epsilon) then
          fint = ( 1.-exp(-0.5*xai) )/ xai
        else
          fint = 0.5
        endif
c
cc   50 continue
c
c step intercepted liquid and snow amounts due to drip/blow,
c intercepted rainfall/snowfall, and min/max limits. also 
c adjust temperature of intercepted precip to current veg temp,
c storing the heat needed to do this in pflux for use in turvap
c 
c without these pfluxes, the implicit turvap calcs could not
c account for the heat flux associated with precip adjustments,
c especially changes of phase (see below), and so could not
c handle equilibrium situations such as intercepted snowfall
c being continuously melted by warm atmos fluxes, with the veg 
c temp somewhat lower than the equil atmos temp to supply heat
c that melts the incoming snow; (turvap would just change veg 
c temp to atmos equil, with little sensible heat storage...then
c final phase adjustment would return veg temp to melt point)
c
c the use of the current (ie, previous timestep's) veg temp 
c gives the best estimate of what this timestep's final temp
c will be, at least for steady conditions
c
      rwork = 1. / dtime
c
cc      do 100 i=1,npoi
c    
c liquid
c
        drip = xai*wliq/tdrip
        wliq = wliq * (1.-dtime/tdrip)
c
        wliq = wliq + dtime*rain*fint
        pflux = rain*fint * (tveg-train)*ch2o
        rain = rain*(1.-xai*fint)
c
        x = wliq
        wliq = min (wliq, wliqmax)
        if (wliq.lt.wliqmin) wliq = 0.
        drip = drip + xai*(x-wliq)*rwork
c
c snow
c
        blow = xai*wsno/tblow
        wsno = wsno * (1.-dtime/tblow)
c
        wsno = wsno + dtime*snow*fint
        pflux = pflux + snow*fint * (tveg-tsnow)*cice
        snow = snow*(1.-xai*fint)
c
        x = wsno
        wsno = min (wsno, wsnomax)
        if (wsno.lt.wsnomin) wsno = 0. 
        blow = blow + xai*(x-wsno)*rwork
c
cc  100 continue
c
c change phase of liquid/snow below/above melt point, and add
c required heat to pflux (see comments above). this will only
c affect the precip intercepted in this timestep, since original
c wliq, wsno must have been ge/le melt point (ensured in later
c call to cascad2/steph2o2)
c
      rwork2 = ch2o - cice
c
cc      do 300 i=1,npoi
c
c liquid below freezing
c
        dw = 0.
        if (tveg.lt.tmelt)  dw = wliq
c
        pflux = pflux
     >           + dw * (rwork2*(tmelt-tveg) - hfus) * rwork
        wliq = wliq - dw
        wsno = wsno + dw
c
c snow above freezing
c
        dw = 0.
        if (tveg.gt.tmelt)  dw = wsno
c
        pflux = pflux
     >           + dw * (rwork2*(tveg-tmelt) + hfus) * rwork
        wsno = wsno - dw
        wliq = wliq + dw
c
cc  300 continue
c
c adjust rainfall, snowfall below veg for interception 
c and drip, blowoff
c
      call mix (rain,train, rain,train, drip,tveg, vzero,vzero)
      call mix (snow,tsnow, snow,tsnow, blow,tveg, vzero,vzero)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine cascad2
c ---------------------------------------------------------------------
c
c at end of timestep, removes evaporation from intercepted h2o,
c and does final heat-conserving adjustment for any liquid/snow 
c below/above melt point. calls steph2o2 for upper leaves, 
c upper stems and lower veg in turn.
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
      include 'comveg.h'
      include 'com1d.h'
c
c local variables
c
cc      integer i           ! loop indice
c
      real fveg,    ! fractional areal coverage of veg component
     >     xai      ! lai and/or sai for veg component
c
c ---------------------------------------------------------------------
c
c set up for upper leaves
c
cc      do 100 i=1,npoi
        fveg = fu
        xai = 2.0 * lai(2)
cc  100 continue
c
c step upper leaves
c
      call steph2o2 (tu,wliqu,wsnou,fveg,xai,rliqu,fvapuw,chu)
c
c set up for upper stems
c
cc      do 200 i=1,npoi
        fveg = fu
        xai = 2.0 * sai(2)
cc  200 continue
c
c step upper stems
c
      call steph2o2 (ts,wliqs,wsnos,fveg,xai,rliqs,fvaps,chs)
c
c set up for lower veg
c
cc      do 400 i=1,npoi
        fveg = (1.-fi)*fl
        xai = 2.0 * (lai(1) + sai(1))
cc  400 continue
c
c step lower veg
c
      call steph2o2 (tl,wliql,wsnol,fveg,xai,rliql,fvaplw,chl)
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine steph2o2 (tveg,wliq,wsno,fveg,xai,rliq,fvapw,cveg)
c ---------------------------------------------------------------------
c
c removes evaporation from intercepted h2o, and does final
c heat-conserving adjustment for any liquid/snow below/above
c melt point, for one veg component
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
c
c Arguments (all arguments are supplied unless otherwise noted)
c
      real cveg        ! specific heat of veg component ch[u,s,l] 
c
      real      
     >  tveg,    ! temperature of veg component t[u,s,l] (returned)
     >  wliq,    ! intercepted liquid amount wliq[u,s,l] (returned)
     >  wsno,    ! intercepted snow amount wsno[u,s,l] (returned)
     >  fveg,    ! fractional areal coverage, fu or (1-fi)*fl
     >  xai,     ! lai, sai, lai+sai for upper leaves/stems,lower veg
     >  rliq,    ! ratio of area wetted by liquid to total wetted area
     >  fvapw    ! wetted evap h2o flx per leaf/stem area fvap[uw,s,lw]
c
c local variables
c
cc      integer i        ! loopi indice
c
      real zm,         ! to compute corrective fluxes
     >     rwork,      ! 1/specific heat of fusion 
     >     chav        ! average specific heat for veg, liw and snow
c
      real dh,   ! correct heat flux for liquid below melt point and opposite
     >     dw    ! correct water flux for liquid below melt point and opposite
c
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------
c
c step intercepted h2o due to evaporation/sublimation.
c (fvapw already has been multiplied by fwet factor in turvap,
c so it is per unit leaf/stem area.)
c
c due to linear fwet factors (see comments in fwetcal) and
c the cap on suw,ssw,slw in turvap, evaporation in one timestep
c should hardly ever make wliq or wsno negative -- but if this
c happens, compensate by increasing vapor flux from atmosphere, 
c and decreasing sensib heat flux from atmos (the former is
c dangerous since it could suck moisture out of a dry atmos,
c and both are unphysical but do fix the budget) tveg in hvapf
c and hsubf should be pre-turvap-timestep values, but are not
c
cc      do 100 i = 1, npoi
c
        wliq = wliq - dtime *     rliq  * fvapw
        wsno = wsno - dtime * (1.-rliq) * fvapw
c
c check to see if predicted wliq or wsno are less than zero
c
        if ((wliq.lt.0. or. wsno.lt.0.)
     >      .and. fveg*xai.gt.0. )  then
c
c         write (*,9999) i, wliq, wsno
c9999     format(' ***warning: wliq<0 or wsno<0 -- steph2o2 9999',
c    >           ' i, wliq, wsno:',i4, 2f12.6)
c
c calculate corrective fluxes
c
          zm = max (-wliq, 0.) * fveg * xai / dtime
          fvapa = fvapa + zm
          fsena = fsena - zm*hvapf(tveg,ta)
          wliq = max (wliq, 0.)
c
          zm = max (-wsno, 0.) * fveg * xai / dtime
          fvapa = fvapa + zm
          fsena = fsena - zm*hsubf(tveg,ta)
          wsno = max (wsno, 0.)
c
        endif
c
cc  100 continue
c
c final heat-conserving correction for liquid/snow below/above
c melting point
c
      rwork = 1. / hfus
c
cc      do 200 i=1,npoi
c
        chav = cveg + ch2o*wliq + cice*wsno
c
c correct for liquid below melt point
c
c (nb: if tveg > tmelt or wliq = 0, nothing changes.)
c
        if (tveg.lt.tmelt .and. wliq.gt.0.0) then
          dh = chav*(tmelt - tveg)
          dw = min (wliq, max (0., dh*rwork))
          wliq = wliq - dw
          wsno = wsno + dw 
          chav = cveg + ch2o*wliq + cice*wsno
          tveg = tmelt - (dh-hfus*dw)/chav
        endif
c
c correct for snow above melt point
c
c (nb: if tveg < tmelt or wsno = 0, nothing changes.)
c
        if (tveg.gt.tmelt .and. wsno.gt.0.0) then
          dh = chav*(tveg - tmelt)
          dw = min (wsno, max (0., dh*rwork))
          wsno = wsno - dw
          wliq = wliq + dw
          chav = cveg + ch2o*wliq + cice*wsno
          tveg = tmelt + (dh-hfus*dw)/chav
        endif
c
cc  200 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine mix (xm,tm, x1,t1, x2,t2, x3,t3)
c ---------------------------------------------------------------------
c
c calorimetrically mixes masses x1,x2,x3 with temperatures
c t1,t2,t3 into combined mass xm with temperature tm
c
c xm,tm may be returned into same location as one of x1,t1,..,
c so hold result temporarily in xtmp,ttmp below
c
c will work if some of x1,x2,x3 have opposite signs, but may 
c give unphysical tm's
c
      include 'implicit.h'
c
      include 'compar.h'
c
c Arguments (input except for xm, tm)
c
      real xm,     ! resulting mass  
     >     tm,     ! resulting temp
     >     x1,     ! mass 1
     >     t1,     ! temp 1
     >     x2,     ! mass 2
     >     t2,     ! temp 2
     >     x3,     ! mass 3
     >     t3      ! temp 3
c
c local variables
c
cc      integer i          ! loop indice
c
      real xtmp,         ! resulting mass (storing variable)
     >     ytmp,         !  "
     >     ttmp          ! resulting temp
c
c ---------------------------------------------------------------------
c
cc      do 100 i=1,npoi
c
        xtmp = x1 + x2 + x3
c
        ytmp = sign (max (abs(xtmp), epsilon), xtmp)
c
        if (abs(xtmp).ge.epsilon) then
          ttmp = (t1*x1 + t2*x2 + t3*x3) / ytmp
        else
          ttmp = 0.
          xtmp = 0.
        endif
c
        xm = xtmp
        tm = ttmp
c
cc  100 continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine noveg
c ---------------------------------------------------------------------
c
c if no veg surfaces exist, set prog vars to nominal values
c
c (sensible fluxes fsen[u,s,l], latent fluxes fvap[u,s,l]*, 
c temperature t[u,s,l], and intercepted liquid, snow amounts 
c wliq[u,s,l], wsno[u,s,l] have been calculated for a unit 
c leaf/stem surface, whether or not one exists.)
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
c
c local variables
c
cc      integer i   ! loop indice
c
      real tav,   ! average temp for soil and snow 
     >     x,     ! total lai + sai
     >     y      ! fraction of lower canopy not snow covered 
c
cc      do 100 i = 1, npoi
c
        tav = (1.-fi)*tg + fi*ti
c
        if (lai(2).eq.0. .or. fu.eq.0.) then
          tu = tav
          wliqu = 0.
          wsnou = 0.
        endif
c
        if (sai(2).eq.0. .or. fu.eq.0.) then
          ts = tav
          wliqs = 0.
          wsnos = 0.
        endif 
c
        x = 2.0 * (lai(1) + sai(1))
        y = fl*(1.-fi)
c
        if (x .eq.0. .or. y.eq.0.) then
          tl = tav 
          wliql = 0.
          wsnol = 0.
        endif
c
cc  100 continue
c
      return
      end
c
c ------------------------------------------------------------------------
      real function twet3(tak, q, p)
c ------------------------------------------------------------------------
c
c twet3.f last update 8/30/2000 C Molling
c
c This function calculates the wet bulb temperature given
c air temp, specific humidity, and air pressure.  It needs the function esat
c in order to work (in comsat.h).  The function is an approximation to
c the actual wet bulb temperature relationship.  It agrees well with the
c formula in the Smithsonian Met. Tables for moderate humidities, but differs
c by as much as 1 K in extremely dry or moist environments.
c
c INPUT
c     tak - air temp in K
c     q - specific humidity in kg/kg
c     p - air pressure in Pa (Pa = 100 * mb)
c
c OUTPUT
c     twet3 - wet bulb temp in K, accuracy?
c
      include 'implicit.h'
      include 'compar.h'
c
cc      integer i
c
      real tak, q, p, ta, twk, twold, diff
c
      include 'comsat.h'
c
c temperatures in twet3 equation must be in C
c pressure in qsat function must be in Pa
c temperatures in esat,hvapf functions must be in K
c
c     Air temp in C
c     -------------
      ta = tak - 273.16
c
c     First guess for wet bulb temp in C, K
c     -------------------------------------
      twet3 = ta * q / qsat(esat(tak),p)
      twk = twet3 + 273.16
c
c     Iterate to converge
c     -------------------
cc      do 100 i = 1, 20
         twold = twk - 273.16
         twet3 = ta - (hvapf(twk,tak)/cair) * ( qsat( esat(twk),p )-q )
         diff = twet3 - twold
c
c below, the 0.2 is the relaxation parameter that works up to 40C (at least)
c
         twk = twold + 0.2 * diff + 273.16
         if (abs(twk-273.16-twold) .lt. 0.02) goto 999
cc 100  continue
c
cyuan      print *, 'Warning, twet3 failed to converge after 20 iterations!'
cyuan      print *, 'twet3, twetold: ', twk, twold+273.16
cyuan      print *, 'twetbulb is being set to the air temperature'
c
      twet3 = tak
c
c     Return wet bulb temperature in K
c     --------------------------------
 999  twet3 = twk
c
      return
      end
c

