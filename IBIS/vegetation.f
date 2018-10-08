$FIXEDFORMLINESIZE:132
c
c #    #  ######   ####   ######   #####    ##     #####     #     ####   #    #
c #    #  #       #    #  #          #     #  #      #       #    #    #  ##   #
c #    #  #####   #       #####      #    #    #     #       #    #    #  # #  #
c #    #  #       #  ###  #          #    ######     #       #    #    #  #  # #
c  #  #   #       #    #  #          #    #    #     #       #    #    #  #   ##
c   ##    ######   ####   ######     #    #    #     #       #     ####   #    #
c
c ---------------------------------------------------------------------
      subroutine pheno(iday,laisum,jday)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c local variables
c
      integer iday,jday,laimaxjd,flag0
cc     >  i,          !
cc     >  imonth     !
cc     >  iday        !
c
      real
     >  ddays,        !
     >  ddfac,        !
     >  tthreshold,   ! temperature threshold for budburst and senescence
     >  gthreshold,   ! temperature threshold for budburst and senescence
     >  avglaiu,      ! average lai of upper canopy 
     >  avglail,      ! average lai of lower canopy 
**** DTP 2000/06/28 Modified this following discussion with Navin. We now
*    retain fu derived from dynaveg and constrain it to a local value
*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.
     >  fu_phys       ! Local value of fu constrained to range 0.25 to 0.975
                      ! to keep physics calculations stable.
	 real laid,	laisum
c
c define 'drop days' -- number of days to affect phenology change
c
      ddays = 15.0
      ddfac = 1.0 / ddays
c
c ---------------------------------------------------------------------
c * * * upper canopy winter phenology * * *
c ---------------------------------------------------------------------
c
c temperature threshold for budburst and senescence
c
c temperature threshold is assumed to be 0 degrees C 
c or 5 degrees warmer than the coldest monthly temperature
c
        tthreshold = max (0.0 + 273.16, tc + 5.0 + 273.16)
c
c gdd threshold temperature for leaf budburst
c with a growing degree threshold of 100 units
c
        gthreshold = 0.0 + 273.16
c 
c determine if growing degree days are initiated
c
        if (a10td.lt.gthreshold) then
          agddu  = 0.0
        else
          agddu = agddu + td - gthreshold
        endif
c
c determine leaf display
c
        if (a10td.lt.tthreshold) then
          tempu  = max (0.0, tempu - ddfac)  !pheno factor
        else
          tempu = min (1., max (0.0, agddu - 100.0) / 50.0)
        endif
c
c ---------------------------------------------------------------------
c * * * lower canopy winter phenology * * *
c ---------------------------------------------------------------------
c
c temperature threshold for budburst and senescence
c
c temperature threshold is assumed to be 0 degrees C 
c
        tthreshold = 0.0 + 273.16
c
c gdd threshold temperature for leaf budburst
c with a growing degree threshold of 150 units
c
        gthreshold = -5.0 + 273.16
c 
c determine if growing degree days are initiated
c
        if (a10td.lt.gthreshold) then
          agddl  = 0.0
        else
          agddl = agddl + td - gthreshold
        endif
c
c determine leaf display
c
        if (a10td.lt.tthreshold) then
          templ  = max (0.0, templ - ddfac)
        else
          templ = min (1., max (0.0, agddl - 150.0) / 50.0)
        endif
c
c ---------------------------------------------------------------------
c * * * drought canopy winter phenology * * *
c ---------------------------------------------------------------------
c
        if (a10ancub.lt.0.0) dropu = max (0.1, dropu - ddfac)
        if (a10ancub.ge.0.0) dropu = min (1.0, dropu + ddfac)
c
        if (a10ancls.lt.0.0) dropls = max (0.1, dropls - ddfac)
        if (a10ancls.ge.0.0) dropls = min (1.0, dropls + ddfac)
c
        if (a10ancl4.lt.0.0) dropl4 = max (0.1, dropl4 - ddfac)
        if (a10ancl4.ge.0.0) dropl4 = min (1.0, dropl4 + ddfac)
c
        if (a10ancl3.lt.0.0) dropl3 = max (0.1, dropl3 - ddfac)
        if (a10ancl3.ge.0.0) dropl3 = min (1.0, dropl3 + ddfac)
c
c ---------------------------------------------------------------------
c * * * update lai and canopy fractions * * *
c ---------------------------------------------------------------------
c
c upper canopy single sided leaf area index (area-weighted)
c
        avglaiu = plai(1)             +
     >            plai(2) * dropu  +
     >            plai(3)             +
     >            plai(4)             +
     >            plai(5) * tempu  +
     >            plai(6)             +
     >            plai(7) * tempu  +
     >            plai(8) * tempu

cc        write(100,"(10f10.2)")plai(1),plai(2),plai(3),plai(4),plai(5),plai(6),plai(7),plai(8),
cc     >              dropu, tempu	   
c
c upper canopy fractions
c
        frac(1) = plai(1)            / max (avglaiu, epsilon)
        frac(2) = plai(2) * dropu / max (avglaiu, epsilon)
        frac(3) = plai(3)            / max (avglaiu, epsilon)
        frac(4) = plai(4)            / max (avglaiu, epsilon)
        frac(5) = plai(5) * tempu / max (avglaiu, epsilon)
        frac(6) = plai(6)            / max (avglaiu, epsilon)
        frac(7) = plai(7) * tempu / max (avglaiu, epsilon)
        frac(8) = plai(8) * tempu / max (avglaiu, epsilon)
c
c lower canopy single sided leaf area index (area-weighted)
c
        avglail = plai(9)                              +
     >            plai(10) * min (templ, dropls) +
     >            plai(11) * min (templ, dropl4) +
     >            plai(12) * min (templ, dropl3)
c
c lower canopy fractions
c
        frac(9)  = plai(9)                              /
     >               max (avglail, epsilon)
c
        frac(10) = plai(10) * min (templ, dropls) /
     >               max (avglail, epsilon)
c
        frac(11) = plai(11) * min (templ, dropl4) /
     >               max (avglail, epsilon)
c
        frac(12) = plai(12) * min (templ, dropl3) /
     >               max (avglail, epsilon)
c
c calculate the canopy leaf area index using the fractional vegetation cover
c
        lai(1) = avglail / fl

**** DTP 2000/06/28 Modified this following discussion with Navin. We now
*    retain fu derived from dynaveg and constrain it to a local value
*    "fu_phys" in the range 0.25 to 0.975 used in the canopy physics calcs.

        fu_phys = max (0.25, min (0.975, fu))
        lai(2) = avglaiu / fu_phys
        lai(2) = avglaiu / fu
        
cc	  write(100,*)iday,fu,fu_phys,avglaiu,lai(2)

c
c put a fix on canopy lais to avoid problems in physics
c
        lai(1) = min (lai(1), 6.0)  !grassland
        lai(2) = min (lai(2), 8.0)  !forestc      phenology and carbon allocation

c
c ---------------------------------------------------------------------
c * * * update canopy height parameters * * *
c ---------------------------------------------------------------------
c
c update lower canopy height parameters
c
c note that they are based on vegetation fraction and not
c averaged over the entire gridcell
c
        zbot(1)   =  0.05
        ztop(1)   =  max (0.25, lai(1) * 0.25)
c        
c constrain ztop to be at least 0.5 meter lower than 
c zbot for upper canopy
c
        ztop(1) = min (ztop(1), zbot(2) - 0.5) 

c
c end of loop
c
cc 100  continue
c
c return to main program
c 
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine dynaveg (iyear, isimfire) ! , isim_ac, year)	  
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'compft.h'
      include 'comage.h'
c
c Arguments
c
      integer isimfire  ! fire switch
!     >        isim_ac,   ! age-class dynamics switch
!     >        year      ! year of simulation

      real    pfire     ! probability of fire -- should be determined externally.
cc      real        pdist      ! probability of other disturbance types....
       
      PARAMETER (pfire = 1.0) ! for now we just assume it occurs all the time
      
c
c local variables
c
      integer
     >  j,           ! gridcell counter
     >  iyear
c
      real
     >  sapspeed,      ! in mm/day
     >  trans,         ! (2.5 mm/day) 
     >  saparea,       ! in m**2
     >  sapvolume,     ! in m**3
     >  denswood,      ! kg/m**3
     >  wood,          ! total amount of woody biomass in gridcell
     >  taufin         !
*     >  xminlai        !
c
*      real
*     >  aleaf(npft),   ! allocation fraction to leaves
*     >  aroot(npft),   ! allocation fraction to fine roots
*     >  awood(npft),   ! allocation fraction to wood
*     >  tauleaf(npft), ! turnover time of carbon in leaves (years)
*     >  tauroot(npft), ! turnover time of carbon in fine roots (years)
*     >  tauwood(npft)   ! turnover time of carbon in wood (years)
*     >  tauwood0(npft) ! normal (unstressed) turnover time
c
c ibis uses a small number of plant functional types:
c
c  1: tropical broadleaf evergreen tree
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen tree
c  4: temperate conifer evergreen tree
c  5: temperate broadleaf cold-deciduous tree
c  6: boreal conifer evergreen tree
c  7: boreal broadleaf cold-deciduous tree
c  8: boreal conifer cold-deciduous tree
c  9: evergreen shrub
c 10: deciduous shrub
c 11: warm (c4) grass
c 12: cool (c3) grass
c
c ---------------------------------------------------------------------
c * * * specify biomass turnover parameters (years) * * *
c ---------------------------------------------------------------------
c
      data tauleaf / 1.00,   ! tropical broadleaf evergreen trees
     >               1.00,   ! tropical broadleaf drought-deciduous trees
     >               1.00,   ! warm-temperate broadleaf evergreen trees
     >               2.00,   ! temperate conifer evergreen trees
     >               1.00,   ! temperate broadleaf cold-deciduous trees
     >               2.50,   ! boreal conifer evergreen trees
     >               1.00,   ! boreal broadleaf cold-deciduous trees
     >               1.00,   ! boreal conifer cold-deciduous trees
     >               1.50,   ! evergreen shrubs
     >               1.00,   ! deciduous shrubs
     >               1.25,   ! warm (c4) grasses
     >               1.50 /  ! cool (c3) grasses
c
      data tauwood0 / 25.0,  ! tropical broadleaf evergreen trees
     >                25.0,  ! tropical broadleaf drought-deciduous trees
     >                25.0,  ! warm-temperate broadleaf evergreen trees
     >                50.0,  ! temperate conifer evergreen trees
     >                50.0,  ! temperate broadleaf cold-deciduous trees
     >               100.0,  ! boreal conifer evergreen trees
     >               100.0,  ! boreal broadleaf cold-deciduous trees
     >               100.0,  ! boreal conifer cold-deciduous trees
     >                 5.0,  ! evergreen shrubs
     >                 5.0,  ! deciduous shrubs
     >               999.0,  ! warm (c4) grasses
     >               999.0 / ! cool (c3) grasses
c
cc 999  do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * initialize vegetation dynamics pools * * *
c ---------------------------------------------------------------------
c
c zero out litter fall fields
c
        falll = 0.0
        fallr = 0.0
        fallw = 0.0
c
c zero out carbon lost due to disturbance
c 
        cdisturb = 0.0
c
        wood = 0.001
c
c ---------------------------------------------------------------------
c * * * update npp, and pool losses  * * *
c ---------------------------------------------------------------------
c
c go through all the pfts
c
        do 110 j = 1, npft
c
c apply this year's existence arrays to npp
c
          aynpp(j)  = exist(j) * aynpp(j)
c
c determine above-ground npp for each plant type
c
          ayanpp(j) = (aleaf(j) + awood(j)) * aynpp(j)

cc          write(100,*)"dynaveg 0",i,j,tauwood0(j),taufin,exist(j)
c
c determine turnover rates for woody biomass:
c
c if pft can exist,    then tauwood = tauwood0 (normal turnover),
c if pft cannot exist, then tauwood = taufin years (to kill off trees)
c
c          taufin     = 5.0
          taufin     = tauwood0(j)/2.0
          tauwood(j) = tauwood0(j) - (tauwood0(j) - taufin) * (1.0 - exist(j))
c
c assume a constant fine root turnover time
c
          tauroot(j) = 5.0
c
c determine litter fall rates
cc          write(100,*)"dynaveg 1",tauleaf(j),tauroot(j),tauwood(j)
c
          falll = falll + cbiol(j) / tauleaf(j)
          fallr = fallr + cbior(j) / tauroot(j)
          fallw = fallw + cbiow(j) / tauwood(j)
c
c ---------------------------------------------------------------------
c * * * update biomass pools  * * *
c ---------------------------------------------------------------------
c
c update carbon reservoirs using an analytical solution
c to the original carbon balance differential equation
c
          cbiol(j) = cbiol(j) * exp(-1./tauleaf(j))  +
     >                 aleaf(j) * tauleaf(j) * max (0., aynpp(j)) *
     >                 (1. - exp(-1./tauleaf(j)))
c
          cbiow(j) = cbiow(j) * exp(-1./tauwood(j))  +
     >                 awood(j) * tauwood(j) * max (0., aynpp(j)) *
     >                 (1. - exp(-1./tauwood(j)))
c
          cbior(j) = cbior(j) * exp(-1./tauroot(j))  +
     >                 aroot(j) * tauroot(j) * max (0., aynpp(j)) *
     >                 (1. - exp(-1./tauroot(j)))
c
          if (j.le.8) wood = wood + max (0.0, cbiow(j))

 110    continue

cc	write(100,"(6f10.2)")cbiol(5),cbiol(12),cbiow(5),cbiow(12),cbior(5),cbior(12)

c
c ---------------------------------------------------------------------
c * * * apply disturbances * * *
c ---------------------------------------------------------------------
c
c set fixed disturbance regime
c
        disturbf = 0.005
        disturbo = 0.005
        
**** DTP 2000/08/10. One can do a decent test of ACME by setting
*    these disturbance rates to zero. With these values, the area
*    disturbed each year will be zero so the distribution of biomass
*    and PFTs across the domain should be identical to those 
*    resulting from a run of standard IBIS with zero disturbance

*        disturbf = 0.0  ! Test with zero disturbance rate 
*        disturbo = 0.0  ! (This should equal reference sim).

c
c call fire disturbance routine
**** DTP 2001/03/06: In general isimfire should be set to zero if isim_ac
*    is set to 1 (but what does isimfire really do?)
c
        if (isimfire.eq.1) call fire
 
          do 116 j = 1, npft 
c
c calculate biomass (vegetations) carbon lost to atmosphere   
c used to balance net ecosystem exchange  
c
c ---------------------------------------------------------------------
**** DTP 2000/04/22 QUESTION: 
c ---------------------------------------------------------------------
* Shouldn't a portion of the destroyed material be added to litter fall?
c
            cdisturb = cdisturb + 
     >                    cbiol(j) * (disturbf + disturbo) +
     >                    cbiow(j) * (disturbf + disturbo) +
     >                    cbior(j) * (disturbf + disturbo)                  
c          
c adjust biomass pools due to disturbances
c
            cbiol(j) = cbiol(j) * (1. - disturbf - disturbo)
            cbiow(j) = cbiow(j) * (1. - disturbf - disturbo)
            cbior(j) = cbior(j) * (1. - disturbf - disturbo)
c
c constrain biomass fields to be positive
c
            cbiol(j) = max (0.0, cbiol(j))
            cbiow(j) = max (0.0, cbiow(j))
            cbior(j) = max (0.0, cbior(j))

 116      continue
c
c ---------------------------------------------------------------------
c * * * check and update biomass pools following disturbance * * *
c ---------------------------------------------------------------------
c
        do 120 j = 1, npft
c
c maintain minimum value of leaf carbon in areas that plants exist
c
*          xminlai = 0.010
c
          cbiol(j) = max (exist(j) * xminlai / specla(j),cbiol(j))
c
c update vegetation's physical characteristics
c
          plai(j)    = cbiol(j) * specla(j)
          biomass(j) = cbiol(j) + cbiow(j) + cbior(j)
c
 120    continue
c
c ---------------------------------------------------------------------
c * * * update annual npp, lai, and biomass * * *
c ---------------------------------------------------------------------
c
c adjust annual net ecosystem exchange (calculated in stats.f) 
c by loss of carbon to atmosphere due to biomass burning (fire)
c
        ayneetot = ayneetot - cdisturb
c
c determine total ecosystem above-ground npp
c
        ayanpptot = ayanpp(1)  + ayanpp(2) + ayanpp(3)  + ayanpp(4) +
     >              ayanpp(5)  + ayanpp(6) + ayanpp(7)  + ayanpp(8) +
     >              ayanpp(9)  + ayanpp(10) + ayanpp(11) + ayanpp(12)
c

c update total canopy leaf area
c
        totlaiu = plai(1)  + plai(2) + plai(3)  + plai(4) +
     >               plai(5)  + plai(6) + plai(7)  + plai(8)
c
        totlail = plai(9)  + plai(10) + plai(11) + plai(12)
c
c update total biomass
        totbiou = biomass(1) + biomass(2) + biomass(3) + biomass(4) +
     >            biomass(5) + biomass(6) + biomass(7) + biomass(8)
c
        totbiol = biomass(9)  + biomass(10) + biomass(11) + biomass(12)
c
c ---------------------------------------------------------------------
c * * * update fractional cover and vegetation height parameters * * *
c ---------------------------------------------------------------------
c
**** Added these in temporarily for comparison with original code.
**** Delete these from production version....
        fu = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
        fu = fu * (1. - disturbf - disturbo)

c constrain the fractional cover (upper canopy)
c
        fu = max (0.25, min (0.975, fu))
c
c update fractional cover of herbaceous (lower) canopy:
c 
        fl = totlail / 1.0
c
c apply disturbances to fractional cover (lower canopy)
c
        fl = fl * (1. - disturbf - disturbo)
c
c constrain the fractional cover (lower canopy)
c
        fl = max (0.25, min (0.975, fl))
c
c
c annual update upper canopy height parameters
c should be calculated based on vegetative fraction and not the
c average over the entire grid cell
c
        zbot(2) = 3.0
        ztop(2) = max(zbot(2) + 1.00, 2.50 * totbiou / fu * 0.75)
c
c ---------------------------------------------------------------------
c * * * update stem area index and sapwood fraction * * *
c ---------------------------------------------------------------------
c
c estimate stem area index (sai) as a fraction of the lai
c
        sai(1) = 0.050 * totlail
        sai(2) = 0.250 * totlaiu
c
c estimate sapwood fraction of woody biomass
c
        sapspeed  = 25.0                        ! (m/day)
        trans     = 0.0025                      ! (2.5 mm/day) 
        saparea   = (trans / sapspeed)          ! m**2
c
        sapvolume = saparea * ztop(2) * 0.75  ! m**3
c
        denswood  = 400.0                       ! kg/m**3
c
        sapfrac = min (0.50, max (0.05, sapvolume * denswood / wood))
c
cc 100  continue
c
c ---------------------------------------------------------------------
c * * * map out vegetation classes for this year * * *
c ---------------------------------------------------------------------
c
      call vegmap
c
c
c return to the main program
c
      return

      end  ! DYNAVEG
c
c
c ---------------------------------------------------------------------
      subroutine fire
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
c
c local variables
c
cc      integer i
c
      real burn
c
c begin global grid
c
cc     do 100 i = 1, npoi
c
        burn = firefac * min (1.0, totlit / 0.200)
c
        disturbf = 1.0 - exp(-0.5 * burn)
c
        disturbf = max (0.0, min (1.0, disturbf))
c
cc 100  continue
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine vegmap
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comveg.h'
c
c local variables
c
      integer 
     >     j,            ! loop indice
     >     domtree          ! dominant tree
c
      real maxlai,          ! maximum lai
     >     totlai,          ! total ecosystem lai
cc     >     grassfrac,       ! fraction of total lai in grasses
cc     >     treefrac,        ! fraction of total lai in trees
     >     treelai,         ! lai of trees
     >     shrublai,        ! lai of shrubs
     >     grasslai,        ! lai of grass
     >     ratio
c
c classify vegetation cover into standard ibis vegetation classes 
c
c ---------------------------------------------------
c  1: tropical evergreen forest / woodland
c  2: tropical deciduous forest / woodland
c  3: temperate evergreen broadleaf forest / woodland
c  4: temperate evergreen conifer forest / woodland
c  5: temperate deciduous forest / woodland
c  6: boreal evergreen forest / woodland
c  7: boreal deciduous forest / woodland
c  8: mixed forest / woodland
c  9: savanna
c 10: grassland / steppe 
c 11: dense shrubland
c 12: open shrubland
c 13: tundra
c 14: desert 
c 15: polar desert / rock / ice
c ---------------------------------------------------
c
c begin global grid
c
cc      do 100 i = 1, npoi
c
c determine total lai and tree, shrub, and grass fractions
c
        treelai   = totlaiu 
        shrublai  = plai(9)  + plai(10)
        grasslai  = plai(11) + plai(12)
c
        totlai    = max (0.01, totlail + totlaiu)
c
c determine dominant tree type by lai dominance
c
        domtree = 0
        maxlai = 0.0
c
        do 110 j = 1, 8
          if (plai(j).gt.maxlai) then
            domtree = j
            maxlai = plai(j)
          endif
 110    continue
c
c assign initial vegetation type
c
        vegtype0 = -999.99
c
c dominant type:  tropical broadleaf evergreen tree
c
        if (domtree.eq.1) then
          if (treelai.gt.2.5)         vegtype0 =  1.0  ! tropical evergreen forest / woodland
          if (treelai.le.2.5)         vegtype0 =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  tropical broadleaf drought-deciduous tree
c
        if (domtree.eq.2) then
          if (treelai.gt.2.5)         vegtype0 =  2.0  ! tropical deciduous forest / woodland
          if (treelai.le.2.5)         vegtype0 =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  warm-temperate broadleaf evergreen tree
c
        if (domtree.eq.3) then
          if (treelai.gt.2.5)         vegtype0 =  3.0  ! temperate evergreen broadleaf forest / woodland
          if (treelai.le.2.5)         vegtype0 =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate conifer evergreen tree
c
        if (domtree.eq.4) then
          if (treelai.gt.1.5)         vegtype0 =  4.0  ! temperate evergreen conifer forest / woodland
          if (treelai.le.1.5)         vegtype0 =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate broadleaf deciduous tree
c
        if (domtree.eq.5) then
          if (treelai.gt.1.5)         vegtype0 =  5.0  ! temperate deciduous forest / woodland
          if (treelai.le.1.5)         vegtype0 =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  boreal conifer evergreen tree
c
        if (domtree.eq.6)             vegtype0 =  6.0  ! boreal evergreen forest / woodland
c
c       if (domtree.eq.6) then
c         if (treelai.gt.1.0)         vegtype0 =  6.0  ! boreal evergreen forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
c         endif
c       endif
c
c dominant type:  boreal broadleaf cold-deciduous tree
c
        if (domtree.eq.7)             vegtype0 =  7.0  ! boreal deciduous forest / woodland
c
c       if (domtree.eq.7) then
c         if (treelai.gt.1.0)         vegtype0 =  7.0  ! boreal deciduous forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
c         endif
c       endif
c
c dominant type:  boreal conifer cold-deciduous tree
c
        if (domtree.eq.8)             vegtype0 =  7.0  ! boreal deciduous forest / woodland
c
c       if (domtree.eq.8) then
c         if (treelai.gt.1.0)         vegtype0 =  7.0  ! boreal deciduous forest / woodland
c         if (treelai.le.1.0) then
c           if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
c           if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
c         endif
c       endif
c
c temperate/boreal forest mixtures
c
        if ((domtree.ge.4).and.(domtree.le.8)) then
          ratio = (plai(5) + plai(7) + plai(8)) / 
     >            (plai(4) + plai(5) + plai(6) + 
     >             plai(7) + plai(8))
          if (treelai.gt.1.0) then
            if ((ratio.gt.0.45).and.(ratio.lt.0.55)) vegtype0 = 8.
          endif
          if ((domtree.le.5).and.(treelai.le.1.0)) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c no tree is dominant
c
        if (domtree.eq.0) then
          if (treelai.gt.1.0)         vegtype0 =  9.0  ! savanna
          if (treelai.le.1.0) then
            if (grasslai.ge.shrublai) vegtype0 = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0 = 11.0  ! closed shrubland
          endif
        endif
c
c overriding vegtation classifications
c
        if (totlai.lt.1.0)            vegtype0 = 12.0  ! open shrubland
        if (totlai.le.0.4)            vegtype0 = 14.0  ! desert
c
c overriding climatic rules
c
        if (gdd5.lt.350.0) then
          if (totlai.ge.0.4)          vegtype0 = 13.0  ! tundra
          if (totlai.lt.0.4)          vegtype0 = 15.0  ! polar desert
        endif
c
        if (gdd0.lt.100.0)         vegtype0 = 15.0  ! polar desert
c
cc 100  continue
c
c return to the main program
c
      return
      end
c
