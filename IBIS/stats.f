$FIXEDFORMLINESIZE:132
c  ####    #####    ##     #####   ####
c #          #     #  #      #    #
c  ####      #    #    #     #     ####
c      #     #    ######     #         #
c #    #     #    #    #     #    #    #
c  ####      #    #    #     #     ####
c
c
c ---------------------------------------------------------------------
       subroutine sumnow(iyear, imonth, iday, istep)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'comsum.h'
c
c local variables
c
      integer k         ! loop indices
      integer iyear, imonth, iday, istep
c
      real rwood,    ! maintenance respiration coefficient for wood (/s)
     >     rroot,    ! maintenance respiration coefficient for root (/s)
     >     rgrowth,  ! growth respiration coefficient (fraction)
     >     stemtemp, ! stem temperature
     >     roottemp, ! average root temperature for all roots
     >     funca,    ! temperature function for aboveground biomass (stems)
     >     funcb,    ! temperature function for belowground biomass (roots)
     >     zweight,  ! 10-day time averaging factor
     >     smask     ! 1 - fi
c
c ---------------------------------------------------------------------
c * * * define working variables * * *
c ---------------------------------------------------------------------
c
c maintenance respiration coefficients (per second)
c
c initially, we pick values for respiration coefficients that
c defined in units of  / year
c
c   rwood ~ 0.0125 
c   rroot ~ 1.2500
c
c however, we convert the unitsconvert to have resulting respiration
c fluxes in units of mol-C / m**2 / second
c
c this requires we convert the time unit to seconds and add an additional
c factor to convert biomass units from kilograms to moles
c
      rwood   = 0.0125 / (ndaypy * 86400.0) * (1000.0 / 12.0)
      rroot   = 1.2500 / (ndaypy * 86400.0) * (1000.0 / 12.0)
c
c growth respiration coefficient (fraction)
c
      rgrowth = 0.30
c
c 10-day time averaging factor
c
      zweight = exp(-1. / (10.0 * 86400.0 / dtime)) 
c
c begin global grid
c
cc      do 100 i = 1, npoi
c
c calculate instantaneous carbon flux parameters, including
c npp (net primary production) and nee (net ecosystem exchange)
c
c in this routine, all of the fluxes are calculated in the units
c of mol-C / m**2 / sec
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous GPP * * *
c ---------------------------------------------------------------------
c
c snow masking for lower canopy vegetation
c
        smask = 1.0 - fi
c
c note that the following plants types follow different physiological paths
c
c   - broadleaf trees   :  types 1, 2, 3, 5, 7, 8 
c   - conifer   trees   :  types 4, 6
c   - shrubs            :  types 9, 10
c   - c4 grasses        :  type 11
c   - c3 grasses        :  type 12
c
c note that plant type 8 is actually a deciduous conifer (e.g., Larix), but
c we are assuming that it's physiological behavior is like a broadleaf tree
c
c nppdummy is canopy npp before accounting for stem & root respirtation
c Navin Sept 02
c
        nppdummy(1)  = frac(1)  * ancub * lai(2) * fu
        nppdummy(2)  = frac(2)  * ancub * lai(2) * fu
        nppdummy(3)  = frac(3)  * ancub * lai(2) * fu
        nppdummy(4)  = frac(4)  * ancuc * lai(2) * fu
        nppdummy(5)  = frac(5)  * ancub * lai(2) * fu
        nppdummy(6)  = frac(6)  * ancuc * lai(2) * fu
        nppdummy(7)  = frac(7)  * ancub * lai(2) * fu
        nppdummy(8)  = frac(8)  * ancuc * lai(2) * fu
        nppdummy(9)  = frac(9)  * ancls * lai(1) * fl * smask 
        nppdummy(10) = frac(10) * ancls * lai(1) * fl * smask
        nppdummy(11) = frac(11) * ancl4 * lai(1) * fl * smask
        nppdummy(12) = frac(12) * ancl3 * lai(1) * fl * smask
c
c Navin's correction to compute npp using tgpp via agXXX
c agXXX should be used 
c
        tgpp(1)  = frac(1)  * agcub * lai(2) * fu
        tgpp(2)  = frac(2)  * agcub * lai(2) * fu
        tgpp(3)  = frac(3)  * agcub * lai(2) * fu
        tgpp(4)  = frac(4)  * agcuc * lai(2) * fu
        tgpp(5)  = frac(5)  * agcub * lai(2) * fu
        tgpp(6)  = frac(6)  * agcuc * lai(2) * fu
        tgpp(7)  = frac(7)  * agcub * lai(2) * fu
        tgpp(8)  = frac(8)  * agcuc * lai(2) * fu
        tgpp(9)  = frac(9)  * agcls * lai(1) * fl * smask 
        tgpp(10) = frac(10) * agcls * lai(1) * fl * smask
        tgpp(11) = frac(11) * agcl4 * lai(1) * fl * smask
        tgpp(12) = frac(12) * agcl3 * lai(1) * fl * smask

cc        if(iyear.eq.2006)	then
cc           write(110,"(4I6,5f10.6)")1,iyear,imonth,iday,tgpp(2),agcub, lai(2), fu, frac(2)
cc	  end if
c
c calculate total gridcell gpp
c
        tgpptot = 0.0
c
        do 110 k = 1, npft
          tgpptot = tgpptot + tgpp(k)
 110    continue
c
c ---------------------------------------------------------------------
c * * * calculate temperature functions for respiration * * *
c ---------------------------------------------------------------------
c
c calculate the stem temperature
c
        stemtemp = ts
c
c calculate average root temperature (average of all roots)
c
        roottemp = 0.0
c
        do 120 k = 1, nsoilay
          roottemp = roottemp + tsoi(k) * 0.5 *
     >               (froot(k,1) + froot(k,2))
 120    continue

cc        write(100,*)"sumnow 2", stemtemp, roottemp

c
c calculate respiration terms on a 15 degree base
c following respiration parameterization of Lloyd and Taylor
c
        stemtemp = max(stemtemp, 253.0)
        roottemp = max(roottemp, 253.0)

        funca = exp(3500.0 * (1. / 288.16 - 1. / stemtemp))
        funcb = exp(3500.0 * (1. / 288.16 - 1. / roottemp))

cc        write(100,*)"sumnow 2.1"
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous NPP * * *
c ---------------------------------------------------------------------
c
c the basic equation for npp is
c
c   npp = (1 - growth respiration term) * (gpp - maintenance respiration terms)
c
c here the respiration terms are simulated as
c
c   growth respiration = rgrowth * (gpp - maintenance respiration terms)
c
c where
c
c   rgrowth is the construction cost of new tissues
c
c and
c
c   root respiration = rroot * cbior(k) * funcb
c   wood respiration = rwood * cbiow(k) * funca * sapwood fraction
c
c where
c 
c   funca = temperature function for aboveground biomass (stems)
c   funcb = temperature function for belowground biomass (roots)
c
c note that we assume the sapwood fraction for shrubs is 1.0
c
c also note that we apply growth respiration, (1 - rgrowth), 
c throughout the year; this may cause problems when comparing
c these npp values with flux tower measurements
c
c also note that we need to convert the mass units of wood and
c root biomass from kilograms of carbon to moles of carbon
c to maintain consistent units (done in rwood, rroot)
c
c finally, note that growth respiration is only applied to 
c positive carbon gains (i.e., when gpp-rmaint is positive)
c
c Navin fix Sept 02 using nppdummy
        tnpp(1)  = nppdummy(1)                               -
     >               rwood * cbiow(1) * sapfrac * funca -
     >               rroot * cbior(1)              * funcb
c
        tnpp(2)  = nppdummy(2)                               -
     >               rwood * cbiow(2) * sapfrac * funca -
     >               rroot * cbior(2)              * funcb
c
        tnpp(3)  = nppdummy(3)                               -
     >               rwood * cbiow(3) * sapfrac * funca -
     >               rroot * cbior(3)              * funcb
c
        tnpp(4)  = nppdummy(4)                               -
     >               rwood * cbiow(4) * sapfrac * funca -
     >               rroot * cbior(4)              * funcb
c
        tnpp(5)  = nppdummy(5)                               -
     >               rwood * cbiow(5) * sapfrac * funca -
     >               rroot * cbior(5)              * funcb

cc        write(100,*)"sumnow 2.2"

c
        tnpp(6)  = nppdummy(6)                               -
     >               rwood * cbiow(6) * sapfrac * funca -
     >               rroot * cbior(6)              * funcb
c
        tnpp(7)  = nppdummy(7)                               -
     >               rwood * cbiow(7) * sapfrac * funca -
     >               rroot * cbior(7)              * funcb
c
        tnpp(8)  = nppdummy(8)                               -
     >               rwood * cbiow(8) * sapfrac * funca -
     >               rroot * cbior(8)              * funcb
c
        tnpp(9)  = nppdummy(9)                               -
     >               rwood * cbiow(9)              * funca -
     >               rroot * cbior(9)              * funcb

cc        write(100,*)"sumnow 2.5"
c
        tnpp(10) = nppdummy(10)                              -
     >               rwood * cbiow(10)             * funca -
     >               rroot * cbior(10)             * funcb
c
cc        write(100,*)"sumnow 2.6"

        tnpp(11) = nppdummy(11)                              -
     >               rroot * cbior(11)            * funcb
c
        tnpp(12) = nppdummy(12)                              -
     >               rroot * cbior(12)            * funcb
c
c apply growth respiration and calculate total gridcell npp
c
        tnpptot = 0.0
c
        do 130 k = 1, npft
          if (tnpp(k).gt.0.0) tnpp(k) = tnpp(k)  * (1.0 - rgrowth)
          tnpptot = tnpptot + tnpp(k)
 130    continue

cc        write(100,*)"sumnow 3"

c
c ---------------------------------------------------------------------
c * * * calculate total fine root respiration * * *
c ---------------------------------------------------------------------
c
        tco2root = 0.0
c
        do 140 k = 1, npft
          tco2root = tco2root + rroot * cbior(k) * funcb
 140    continue
c
c ---------------------------------------------------------------------
c * * * calculate instantaneous NEE * * *
c ---------------------------------------------------------------------
c
c microbial respiration is calculated in biogeochem.f
c
        tneetot = tnpptot - tco2mic
c
c ---------------------------------------------------------------------
c * * * update 10-day running-mean parameters * * *
c ---------------------------------------------------------------------
c
c 10-day daily air temperature
c
        a10td    = zweight * a10td    + (1. - zweight) * td
c
c 10-day canopy photosynthesis rates
c
        a10ancub = zweight * a10ancub + (1. - zweight) * ancub
        a10ancuc = zweight * a10ancuc + (1. - zweight) * ancuc
        a10ancls = zweight * a10ancls + (1. - zweight) * ancls
        a10ancl3 = zweight * a10ancl3 + (1. - zweight) * ancl3
        a10ancl4 = zweight * a10ancl4 + (1. - zweight) * ancl4

cc        write(100,*)"sumnow 4"

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
      subroutine sumday (iyear,imonth,iday,istep)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments
c
      integer
     >  istep   ! daily timestep number (passed in)
c
c local variables
c
      integer k      ! loop indices
	integer iyear, imonth, iday
c
      real 
     >  rwork,      !working time variable
     >  rwork2,     ! "
     >  rwork3,     ! " 
     >  rwork4,     ! "
     >  tconst,     ! constant for Lloyd and Taylor (1994) function
     >  bconst,     ! base temperature used for carbon decomposition
     >  btemp,      ! maximum value of decomposition factor
     >  rdepth,     ! total depth of the 4 1st soil layers
     >  rdepth2,    ! total depth of the 2 1st soil layers
     >  snodpth,    ! total snow depth
     >  soiltemp,   ! average soil temp for 2 1st layers
     >  soilmois,   ! average soil moisture (fraction of porosity) for 2 1st layers
     >  soilice,    ! average soil ice for 2 1st layers
     >  soitempc,   ! average soil temp over 6 layers
     >  soimoisc,   ! average soil moisture over 6 layers
     >  factor,     ! temperature decomposition factor for ltter/soil carbon
     >  wfps,       ! water filled pore space
     >  moist,      ! moisture effect on decomposition
     >  precipfac
c
cc      real plens
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c
c reset sumday if the first timestep of the day 
c
      if (istep .eq. 1) ndtimes = 0
c
c accumulate daily output (at this point for soil decomposition)
c
      ndtimes = ndtimes + 1
c
c working variables
c
      rwork  = 1. / float(ndtimes)
      rwork2 = 86400.
      rwork3 = 86400. * 12.e-3
      rwork4 = 86400. * 14.e-3
c
c constants used in temperature function for c decomposition
c (arrhenius function constant) 
c
cc      tconst  = 344.00  ! constant for Lloyd and Taylor (1994) function
      tconst  = 200.00  ! constant for Lloyd and Taylor (1994) function

      btemp   = 288.16  ! base temperature used for carbon decomposition
c
      bconst  = 10.0    ! maximum value of decomposition factor

cc	write(100,*)"in sumday1"

c
c soil weighting factors
c
      rdepth  = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
      rdepth2 = 1. / (hsoi(1) + hsoi(2))
c
c soil moisture stress
c
cc      adstresstu = rwork * ((nytimes-1) * adstresstu + stresstu)
cc      adstresstl = rwork * ((nytimes-1) * adstresstl + stresstl)
c
c ---------------------------------------------------------------------
c * * * daily water budget terms * * *
c ---------------------------------------------------------------------
c
cc        adrain    = ((ndtimes-1) * adrain + raina * 86400.) * rwork
cc        adsnow    = ((ndtimes-1) * adsnow + snowa * 86400.) * rwork
c
c Verification of the weather generator algorithm
c
cc        niter = int (86400./dtime)
c
cc        if (istep .eq. niter) then
c
cc          precipfac = precip - adrain - adsnow
c
cc          if ((precipfac .lt. -0.1) .OR. (precipfac .gt. 0.1)) print *,
cc     >         'ERROR in sumday:', adrain + adsnow, precip
cc        endif
c
c End of verification
c
c
c ---------------------------------------------------------------------
c * * * determine daily gpp * * *
c ---------------------------------------------------------------------
c
c gross primary production of each plant type
c
      adgpp(1)  = ((ndtimes-1) * adgpp(1)  + tgpp(1)  * rwork3) * rwork
      adgpp(2)  = ((ndtimes-1) * adgpp(2)  + tgpp(2)  * rwork3) * rwork
      adgpp(3)  = ((ndtimes-1) * adgpp(3)  + tgpp(3)  * rwork3) * rwork
      adgpp(4)  = ((ndtimes-1) * adgpp(4)  + tgpp(4)  * rwork3) * rwork
      adgpp(5)  = ((ndtimes-1) * adgpp(5)  + tgpp(5)  * rwork3) * rwork
      adgpp(6)  = ((ndtimes-1) * adgpp(6)  + tgpp(6)  * rwork3) * rwork
      adgpp(7)  = ((ndtimes-1) * adgpp(7)  + tgpp(7)  * rwork3) * rwork
      adgpp(8)  = ((ndtimes-1) * adgpp(8)  + tgpp(8)  * rwork3) * rwork
      adgpp(9)  = ((ndtimes-1) * adgpp(9)  + tgpp(9)  * rwork3) * rwork
      adgpp(10) = ((ndtimes-1) * adgpp(10) + tgpp(10) * rwork3) * rwork
      adgpp(11) = ((ndtimes-1) * adgpp(11) + tgpp(11) * rwork3) * rwork
      adgpp(12) = ((ndtimes-1) * adgpp(12) + tgpp(12) * rwork3) * rwork

cc	write(100,*)"in sumday2", tgpp(8), adgpp(8)

c
c gross primary production of the entire gridcell
c
      adgpptot = adgpp(1) + adgpp(2) + adgpp(3) + adgpp(4) + adgpp(5) + adgpp(6) + 
     >           adgpp(7) + adgpp(8) + adgpp(9) + adgpp(10) + adgpp(11) + adgpp(12)

cc      if(iyear.eq.2006)	then
cc           write(110,"(4I6,3f10.6)")2,iyear,imonth,iday,tgpp(2),adgpp(2), adgpptot
cc      end if

c
c net primary production of each plant type
c
        adnpp(1)  = ((ndtimes-1) * adnpp(1)  +
     >                  tnpp(1)  * rwork3) * rwork
        adnpp(2)  = ((ndtimes-1) * adnpp(2)  +
     >                  tnpp(2)  * rwork3) * rwork
        adnpp(3)  = ((ndtimes-1) * adnpp(3)  +
     >                  tnpp(3)  * rwork3) * rwork
        adnpp(4)  = ((ndtimes-1) * adnpp(4)  +
     >                  tnpp(4)  * rwork3) * rwork
        adnpp(5)  = ((ndtimes-1) * adnpp(5)  +
     >                  tnpp(5)  * rwork3) * rwork
        adnpp(6)  = ((ndtimes-1) * adnpp(6)  +
     >                  tnpp(6)  * rwork3) * rwork
        adnpp(7)  = ((ndtimes-1) * adnpp(7)  +
     >                  tnpp(7)  * rwork3) * rwork
        adnpp(8)  = ((ndtimes-1) * adnpp(8)  +
     >                  tnpp(8)  * rwork3) * rwork
        adnpp(9)  = ((ndtimes-1) * adnpp(9)  +
     >                  tnpp(9)  * rwork3) * rwork
        adnpp(10) = ((ndtimes-1) * adnpp(10) +
     >                  tnpp(10) * rwork3) * rwork
        adnpp(11) = ((ndtimes-1) * adnpp(11) +
     >                  tnpp(11) * rwork3) * rwork
        adnpp(12) = ((ndtimes-1) * adnpp(12) +
     >                  tnpp(12) * rwork3) * rwork
c
c net primary production of the entire gridcell
c
        adnpptot = adnpp(1)  + adnpp(2)  + adnpp(3)  + 
     >                adnpp(4)  + adnpp(5)  + adnpp(6)  + 
     >                adnpp(7)  + adnpp(8)  + adnpp(9)  +
     >                adnpp(10) + adnpp(11) + adnpp(12)

        adaet     = ((ndtimes-1) * adaet  - fvapa * 86400.) * rwork

c	  
        adtrunoff  = ((ndtimes-1) * adtrunoff  + (grunof + gdrain) * 86400.) * rwork

        adsrunoff  = ((ndtimes-1) * adsrunoff + grunof * 86400.) * rwork
        addrainage = ((ndtimes-1) * addrainage + gdrain * 86400.) * rwork

	  adtrans = ((ndtimes-1) * adtrans + gtrans * 86400.) * rwork  !Yuan Add transpiration
	  adinvap = ((ndtimes-1) * adinvap + ginvap * 86400.) * rwork  !Yuan Add intercepted evaporation
	  adsuvap = ((ndtimes-1) * adsuvap + gsuvap * 86400.) * rwork  !Yuan Add soil/snow surface evaporation	

c
c
c ---------------------------------------------------------------------
c * * * daily snow parameters * * *
c ---------------------------------------------------------------------
c
        snodpth = hsno(1) + hsno(2) + hsno(3)
c
        adsnod = ((ndtimes-1) * adsnod + snodpth) * rwork
        adsnof = ((ndtimes-1) * adsnof + fi)   * rwork
c
c ---------------------------------------------------------------------
c * * * soil parameters * * *
c ---------------------------------------------------------------------
c
c initialize average soil parameters
c
        soiltemp = 0.0
        soilmois = 0.0
        soilice  = 0.0
c
 	  soitempc = 0.0
	  soimoisc = 0.0
c
c averages for first 2 layers of soil
c
        do 110 k = 1, 2
          soiltemp =  soiltemp + tsoi(k)  * hsoi(k)
          soilmois =  soilmois + wsoi(k)  * hsoi(k)
          soilice  =  soilice  + wisoi(k) * hsoi(k)
 110    continue
c
c weighting on just thickness of each layer
c
        soilmois = soilmois * rdepth2
        soilice  = soilice  * rdepth2
        soiltemp = soiltemp * rdepth2
c
c calculate average root temperature, soil temperature and moisture and 
c ice content based on rooting profiles (weighted) from jackson et al
c 1996
c
c these soil moisture and temperatures are used in biogeochem.f 
c we assume that the rooting profiles approximate
c where carbon resides in the soil
c
        do 120 k = 1, nsoilay
          soitempc = soitempc + tsoi(k)  * 0.5 * (froot(k,1) + froot(k,2))
          soimoisc = soimoisc + wsoi(k)  * 0.5 * (froot(k,1) + froot(k,2))
 120    continue

cc 	write(100,*)"in sumday4"

c
c calculate daily average soil moisture and soil ice
c using thickness of each layer as weighting function
c
        adwsoi  = ((ndtimes-1) * adwsoi  + soilmois) * rwork
        adtsoi  = ((ndtimes-1) * adtsoi  + soiltemp) * rwork
        adwisoi = ((ndtimes-1) * adwisoi + soilice)  * rwork
c
c calculate daily average for soil temp/moisture of top layer
c
        adtlaysoi = ((ndtimes-1) * adtlaysoi + tsoi(1)) * rwork
        adwlaysoi = ((ndtimes-1) * adwlaysoi + wsoi(1)) * rwork
c
c calculate separate variables to keep track of weighting using 
c rooting profile information
c
c note that these variables are only used for diagnostic purposes
c and that they are not needed in the biogeochemistry code
c
        adwsoic  = ((ndtimes-1) * adwsoic + soimoisc) * rwork
        adtsoic  = ((ndtimes-1) * adtsoic + soitempc) * rwork
c
c ---------------------------------------------------------------------
c * * * calculate daily soil co2 fluxes * * *
c ---------------------------------------------------------------------
c
c increment daily total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
        adco2mic = ((ndtimes-1) * adco2mic + tco2mic * rwork3) * rwork
c
c increment daily total co2 respiration from fine roots
c tco2root is instantaneous value of co2 flux calculated in stats.f
c
        adco2root = ((ndtimes-1) * adco2root +  tco2root * rwork3) * rwork
c 
c calculate daily total co2 respiration from soil
c
        adco2soi  = adco2root + adco2mic
	  adneetot = adnpptot - adco2mic
c
c calculate daily ratio of total root to total co2 respiration
c
        if (adco2soi.gt.0.0) then
          adco2ratio = adco2root / adco2soi
        else
          adco2ratio = -999.99
        endif
c
c ---------------------------------------------------------------------
c * * * calculate daily litter decomposition parameters * * *
c ---------------------------------------------------------------------
c
c calculate litter carbon decomposition factors
c using soil temp, moisture and ice for top soil layer
c
c calculation of soil biogeochemistry decomposition factors 
c based on moisture and temperature affects on microbial
c biomass dynamics
c
c moisture function based on water-filled pore space (wfps)  
c williams et al., 1992 and friend et al., 1997 used in the
c hybrid 4.0 model; this is based on linn and doran, 1984
c
c temperature functions are derived from arrhenius function
c found in lloyd and taylor, 1994 with a 15 c base 
c
c calculate temperature decomposition factor
c CD impose lower limit to avoid division by zero at tsoi=227.13
c
        if (tsoi(1) .gt. 237.13) then
           factor = min(exp(tconst * ((1./(btemp - 227.13)) - (1./(tsoi(1)-227.13)))), bconst)
        else
           factor = exp(tconst * ((1./(btemp - 227.13)) - (1./(237.13-227.13))))
        end if


c
c calculate water-filled pore space (in percent)
c
c wsoi is relative to pore space not occupied by ice and water
c thus must include the ice fraction in the calculation
c	
        wfps = (1.0 - wisoi(1)) * wsoi(1) * 100.0	
c
c calculate moisture decomposition factor
c
        if (wfps .ge. 60.0) then
          moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
        else
          moist = exp((wfps - 60.0)**2 / (-800.0))	
        endif

cc		write(100,*)"decompl",moist,factor

c
c calculate combined temperature / moisture decomposition factor
c
        factor = max (0.001, min (bconst, factor * moist))
c
c calculate daily average litter decomposition factor
c
        decompl = ((ndtimes-1) * decompl + factor) * rwork

cc	 	if(iyear.eq.2006) then
cc            write(120,"(3I6,8f10.2)")iyear,imonth,iday, wfps, moist, factor, tsoi(1), decompl, 
cc     >  		rwork, wisoi(1), wsoi(1)
cc	     end if
c
c ---------------------------------------------------------------------
c * * * calculate daily soil carbon decomposition parameters * * *
c ---------------------------------------------------------------------
c
c calculate soil carbon decomposition factors
c using soil temp, moisture and ice weighted by rooting profile scheme 
c
c calculation of soil biogeochemistry decomposition factors 
c based on moisture and temperature affects on microbial
c biomass dynamics
c
c moisture function based on water-filled pore space (wfps)  
c williams et al., 1992 and friend et al., 1997 used in the
c hybrid 4.0 model; this is based on linn and doran, 1984
c
c temperature functions are derived from arrhenius function
c found in lloyd and taylor, 1994 with a 15 c base 
c
c calculate temperature decomposition factor
c CD: impose lower limit to avoid division by zero at tsoi=227.13
c
        if (soiltemp .gt. 237.13) then
           factor = min (exp(tconst * ((1./(btemp - 227.13)) - (1./(soiltemp - 227.13)))), bconst)
        else
           factor = exp(tconst * ((1./ (btemp - 227.13)) - (1./(237.13-227.13))))
        end if
c
c calculate water-filled pore space (in percent)
c
c wsoi is relative to pore space not occupied by ice and water
c thus must include the ice fraction in the calculation
c	
        wfps = (1. - soilice) * soilmois * 100.0	
c
c calculate moisture decomposition factor
c
      if (wfps .ge. 60.0) then
          moist = 0.000371 * (wfps**2) - (0.0748 * wfps) + 4.13
        else
          moist = exp((wfps - 60.0)**2 / (-800.0))
	end if
c
c calculate combined temperature / moisture decomposition factor
c
        factor = max (0.001, min (bconst, factor * moist))
c
c calculate daily average soil decomposition factor
c
        decomps = ((ndtimes-1) * decomps + factor) * rwork
c
c ---------------------------------------------------------------------
c * * * calculate other daily biogeochemical parameters * * *
c ---------------------------------------------------------------------
c
c increment daily total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
        adnmintot = ((ndtimes-1) * adnmintot + tnmin * rwork4) * rwork
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
      subroutine summonth (iday, imonth)
c ---------------------------------------------------------------------
c
c first convert to units that make sense for output
c
c   - convert all temperatures to deg c
c   - convert all liquid or vapor fluxes to mm/day
c   - redefine upwd directed heat fluxes as positive
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'comsum.h'
c
c Arguments (input)
c
      integer
     >  istep,     ! daily timestep number (passed in)
     >  iday,      ! day number  (passed in)
     >  imonth     ! month number (passed in)
c
c local variables
c
      integer
     >  k       ! loop indices
c
      real 
     >  rwork,     ! time work variable
     >  rwork2,    !
     >  rwork3,    !
     >  rwork4,    !
     >  rdepth    ! 1/total soil depth over 4 1st layers
cc     >  solartot,  ! total incoming radiation (direct + diffuse, visible + nearIR)
cc     >  soiltemp,  ! average soil temp for 4 1st layers
cc     >  soilmois,  ! average soil moisture for 4 1st layers 
cc     >  soilice,   ! average soil ice for 4 1st layers 
cc     >  vwc,       ! total liquid + ice content of 4 1st layers
cc     >  awc,       ! total available water (+ ice) content of 4 1st layer
cc     >  snodpth    ! total snow depth
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c 
c if the first timestep of the month then reset averages
c
      if (iday.eq.1) nmtimes = 0
cc      if ((istep.eq.1).and.(iday.eq.1)) nmtimes = 0
c
c accumulate terms
c
      nmtimes = nmtimes + 1
c
c working variables
c
c rwork4 for conversion of nitrogen mineralization (moles)
c
      rwork  = 1. / float(nmtimes)
      rwork2 = float(ndaypm(imonth)) 
      rwork3 = float(ndaypm(imonth)) * 12.e-3
      rwork4 = float(ndaypm(imonth)) * 14.e-3
c
cc      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
c
c soil moisture stress
c
cc      amstresstu = rwork * ((nytimes-1) * amstresstu + adstresstu)
cc      amstresstl = rwork * ((nytimes-1) * amstresstl + adstresstl)
c
c ---------------------------------------------------------------------
c * * * monthly water budget terms * * *
c ---------------------------------------------------------------------
c
cc        amrain    = ((nmtimes-1) * amrain + adrain) * rwork
cc        amsnow    = ((nmtimes-1) * amsnow + adsnow) * rwork
        amaet     = ((nmtimes-1) * amaet  - adaet) * rwork
cc        amtrunoff  = ((nmtimes-1) * amtrunoff  + adtrunoff) * rwork
cc        amsrunoff  = ((nmtimes-1) * amsrunoff  + adsrunoff) * rwork
cc        amdrainage = ((nmtimes-1) * amdrainage + addrainage) * rwork
c
c ---------------------------------------------------------------------
c * * * monthly atmospheric terms * * *
c ---------------------------------------------------------------------
c
cc        amtemp  = ((nmtimes-1) * amtemp  + ta - 273.16)  * rwork
cc        amcloud = ((nmtimes-1) * amcloud + cloud * 100.) * rwork
cc        amqa    = ((nmtimes-1) * amqa    + qa)           * rwork
cc        amrh    = ((nmtimes-1) * amrh    + rh)           * rwork
c
c ---------------------------------------------------------------------
c * * * energy budget terms * * *
c ---------------------------------------------------------------------
c
cc        solartot = solad(1) + solad(2) + solai(1) + solai(2)
cc        amsolar  = amsolar + solartot
cc        amirup   = amirup + firb
cc        amirdown = amirdown + fira
cc        amsens   = amsens -  fsena
cc        amlatent = amlatent - fvapa * hvap
c
c ---------------------------------------------------------------------
c * * * monthly vegetation parameters * * *
c ---------------------------------------------------------------------
c
cc        amlaiu = ((nmtimes-1) * amlaiu + fu * lai(2)) * rwork
cc        amlail = ((nmtimes-1) * amlail + fl * lai(1)) * rwork
c
c ---------------------------------------------------------------------
c * * * monthly soil parameters * * *
c ---------------------------------------------------------------------
c
cc        amtsoi  = ((nmtimes-1) * amtsoi  + adtsoi) * rwork
cc        amwsoi  = ((nmtimes-1) * amwsoi  + adwsoi) * rwork
cc        amwisoi = ((nmtimes-1) * amwisoi + adwisoi) * rwork
cc        amvwc   = ((nmtimes-1) * amvwc   + advwc) * rwork
cc        amawc   = ((nmtimes-1) * amawc   + adawc) * rwork
c
c ---------------------------------------------------------------------
c * * * snow parameters * * *
c ---------------------------------------------------------------------
c
cc        snodpth = hsno(1) + hsno(2) + hsno(3)		 
cc        amsnod = ((nmtimes-1) * amsnod + adsnod) * rwork
cc        amsnof = ((nmtimes-1) * amsnof + adsnof) * rwork
c
c ---------------------------------------------------------------------
c * * * determine monthly npp * * *
c ---------------------------------------------------------------------
c
cc        amnpp(1)  = amnpp(1)  + adnpp(1)
cc        amnpp(2)  = amnpp(2)  + adnpp(2)
cc        amnpp(3)  = amnpp(3)  + adnpp(3)
cc        amnpp(4)  = amnpp(4)  + adnpp(4)
cc        amnpp(5)  = amnpp(5)  + adnpp(5)
cc        amnpp(6)  = amnpp(6)  + adnpp(6)
cc        amnpp(7)  = amnpp(7)  + adnpp(7)
cc        amnpp(8)  = amnpp(8)  + adnpp(8)
cc        amnpp(9)  = amnpp(9)  + adnpp(9)
cc        amnpp(10) = amnpp(10) + adnpp(10)
cc        amnpp(11) = amnpp(11) + adnpp(11)
cc        amnpp(12) = amnpp(12) + adnpp(12)
c
cc        amnpptot = amnpp(1) + amnpp(2) + amnpp(3) + amnpp(4) + amnpp(5) + amnpp(6) + 
cc     >             amnpp(7) + amnpp(8) + amnpp(9) + amnpp(10) + amnpp(11) + amnpp(12)

        amnpptot = amnpptot + adnpptot

cc        amgpp(1)  = amgpp(1)  + adgpp(1)
cc        amgpp(2)  = amgpp(2)  + adgpp(2)
cc        amgpp(3)  = amgpp(3)  + adgpp(3)
cc        amgpp(4)  = amgpp(4)  + adgpp(4)
cc        amgpp(5)  = amgpp(5)  + adgpp(5)
cc        amgpp(6)  = amgpp(6)  + adgpp(6)
cc        amgpp(7)  = amgpp(7)  + adgpp(7)
cc        amgpp(8)  = amgpp(8)  + adgpp(8)
cc        amgpp(9)  = amgpp(9)  + adgpp(9)
cc        amgpp(10) = amgpp(10) + adgpp(10)
cc        amgpp(11) = amgpp(11) + adgpp(11)
cc        amgpp(12) = amgpp(12) + adgpp(12)
c
cc        amgpptot = amgpp(1) + amgpp(2) + amgpp(3) + amgpp(4) + amgpp(5) + amgpp(6) + 
cc     >             amgpp(7) + amgpp(8) + amgpp(9) + amgpp(10) + amgpp(11) + amgpp(12)

        amgpptot = amgpptot + adgpptot

c ---------------------------------------------------------------------
c * * * monthly biogeochemistry parameters * * *
c ---------------------------------------------------------------------
c
c increment monthly total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
cc        amco2mic = amco2mic + adco2mic
c
c increment monthly total co2 respiration from roots
c tco2root is instantaneous value of co2 flux calculated in stats.f
c
cc        amco2root = amco2root + adco2root
c
c calculate average total co2 respiration from soil
c
        amco2soi  = amco2root + amco2mic

c******* Yuan added more monthly variables ****/
        amclitlm = ((nmtimes-1) * amclitlm + clitlm) * rwork     
        amclitls = ((nmtimes-1) * amclitls + clitls) * rwork     	  
        amclitll = ((nmtimes-1) * amclitll + clitll) * rwork 
        amclitrm = ((nmtimes-1) * amclitrm + clitrm) * rwork 
        amclitrs = ((nmtimes-1) * amclitrs + clitrs) * rwork 
        amclitrl = ((nmtimes-1) * amclitrl + clitrl) * rwork
        amclitwm = ((nmtimes-1) * amclitwm + clitwm) * rwork
        amclitws = ((nmtimes-1) * amclitws + clitws) * rwork 
        amclitwl = ((nmtimes-1) * amclitwl + clitwl) * rwork
        amtlit = ((nmtimes-1) * amtlit + totlit) * rwork
        amalit = ((nmtimes-1) * amalit + totalit) * rwork
        amblit = ((nmtimes-1) * amblit + totrlit) * rwork 
        amtotcsoi = ((nmtimes-1) * amtotcsoi + totcsoi) * rwork 
        amtotfall = ((nmtimes-1) * amtotfall + totfall) * rwork 
c  
c  calculate ratio of root to total co2 respiration
c
cc        if (amco2soi.gt.0.0) then
cc          amco2ratio = amco2root / amco2soi
cc        else
cc          amco2ratio = -999.99
cc        endif
c 
c  monthly net ecosystem co2 flux -- npp total minus microbial respiration 
c  the npp total includes losses from root respiration
c
        amneetot  = amnpptot - amco2mic 
c
c increment monthly total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
cc        amnmintot = amnmintot + adnmintot

c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine sumyear (imonth, iday)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'com1d.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments (input)
c
      integer
     >  imonth,     ! month number (passed in)
     >  iday
c
c local variables
c
      integer
     >  k        ! loop indices
c
      real 
     >  rwork,     !
     >  rwork2,    !
     >  rwork3,    !
     >  rwork4,    !
     >  rdepth,    ! 1/total soil depth over 4 1st layers
cc     >  solartot,  ! total incoming radiation (direct + diffuse, visible + nearIR)
cc     >  soiltemp,  ! average soil temp for 4 1st layers
cc     >  soilmois,  ! average soil moisture for 4 1st layers 
cc     >  soilice,   ! average soil ice for 4 1st layers 
cc     >  vwc,       ! total liquid + ice content of 4 1st layers
cc     >  awc,       ! total available water (+ ice) content of 4 1st layer
cc     >  water,     ! fire factor: total water content of 1st layer (liquid+ice)
cc     >  waterfrac, ! fire factor: available water content of 1st layer
cc     >  fueldry,   ! fire factor
     >  allroots  ! annual average root biomass
cc     >  wtotp      ! total water stored in soil+vegetation+snow at previous timestep
c
c ---------------------------------------------------------------------
c * * * update counters and working variables * * *
c ---------------------------------------------------------------------
c
cc      write(100,*)"in sumyear 1"
c reset sumyear if the first timestep of the year
c
cc      if ((istep.eq.1).and.(iday.eq.1).and.(imonth.eq.1)) nytimes = 0
      if ((imonth.eq.1).and.(iday.eq.1)) nytimes = 0
c
c accumulate yearly output
c
      nytimes = nytimes + 1
c
c working variables
c
c rwork4 is for nitrogen mineralization conversion
c
      rwork  = 1. / float(nytimes)
      rwork2 = 12
      rwork3 = 12 * 12.e-3
      rwork4 = 12 * 14.e-3

cc      rwork2 = float(ndaypy) * 86400.
cc      rwork3 = float(ndaypy) * 86400. * 12.e-3
cc      rwork4 = float(ndaypy) * 86400. * 14.e-3

cc      rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
c
c begin global grid
c
cc      do 100 i = 1, npoi
c
c ---------------------------------------------------------------------
c * * * annual energy budget terms * * *
c ---------------------------------------------------------------------
c
cc        solartot = solad(1) + solad(2) + solai(1) + solai(2)
cc        aysolar  = aysolar  + amsolar
cc        ayirup   = ayirup   + amirup
cc        ayirdown = ayirdown + amirdown
cc        aysens   = aysens   - amsens
cc        aylatent = aylatent - amlatent
c
c ---------------------------------------------------------------------
c * * * annual water budget terms * * *
c ---------------------------------------------------------------------
c
cc        ayprcp     = ayprcp  + amprcp
cc        ayaet      = ayaet   - amaet
cc        aytrans    = aytrans + amtransc
cc        aytrunoff  = aytrunoff  + amtrunoff
cc        aysrunoff  = aysrunoff  + amsrunoff
cc        aydrainage = aydrainage + amdrainage

        ayaet      = ayaet   + adaet
        aytrans    = aytrans + adtrans
	  ayinvap = ayinvap + adinvap
	  aysuvap = aysuvap + adsuvap
        aytrunoff  = aytrunoff  + adtrunoff
        aysrunoff  = aysrunoff  + adsrunoff
        aydrainage = aydrainage + addrainage

	

c
c---------------------------------------------------------------------
c CD
c estimate the change in soil-vegetation water content. Used to check 
c mass conservation
c---------------------------------------------------------------------
c
cc        wtotp = wtot
cc        wtot = (wliqu+wsnou) * fu * 2.0 * lai(2) +
cc     >            (wliqs+wsnos) * fu * 2.0 * sai(2) +
cc     >            (wliql+wsnol) * fl * 2.0 *
cc     >            (lai(1) + sai(1)) * (1. - fi)
c
cc        wtot = wtot + wpud + wipud
c
cc        do 10 k = 1, nsoilay
cc          wtot = wtot + poros(k)*wsoi(k)*(1.-wisoi(k))*hsoi(k)*rhow +
cc     >           poros(k)*wisoi(k)*hsoi(k)*rhow
cc 10     continue
c
cc        do 20 k = 1, nsnolay
cc          wtot = wtot + fi*rhos*hsno(k)
cc 20     continue
c
cc        aydwtot = ((nytimes-1) * aydwtot + wtot - wtotp) * rwork
c
c ---------------------------------------------------------------------
c * * * annual soil parameters * * *
c ---------------------------------------------------------------------
c annual average soil moisture and soil ice
c
cc        aywsoi  = ((nytimes-1) * aywsoi  +  amwsoi) * rwork
cc        aywisoi = ((nytimes-1) * aywisoi + amwisoi) * rwork
cc        aytsoi  = ((nytimes-1) * aytsoi  +  amtsoi) * rwork
cc        ayvwc   = ((nytimes-1) * ayvwc   +   amvwc) * rwork
cc        ayawc   = ((nytimes-1) * ayawc   +   amawc) * rwork
c
c soil moisture stress
c
cc        aystresstu = rwork * ((nytimes-1) * aystresstu + amstresstu)
cc        aystresstl = rwork * ((nytimes-1) * aystresstl + amstresstl)

        aywsoi  = ((nytimes-1) * aywsoi  +  adwsoi) * rwork
        aywisoi = ((nytimes-1) * aywisoi + adwisoi) * rwork
        aytsoi  = ((nytimes-1) * aytsoi  +  adtsoi) * rwork
cc        ayvwc   = ((nytimes-1) * ayvwc   +   advwc) * rwork
cc        ayawc   = ((nytimes-1) * ayawc   +   adawc) * rwork
c
c ---------------------------------------------------------------------
c * * * determine annual gpp * * *
c ---------------------------------------------------------------------
c
c gross primary production of each plant type
c
        aygpp(1)  = aygpp(1)  + adgpp(1)
        aygpp(2)  = aygpp(2)  + adgpp(2)
        aygpp(3)  = aygpp(3)  + adgpp(3)
        aygpp(4)  = aygpp(4)  + adgpp(4)
        aygpp(5)  = aygpp(5)  + adgpp(5)
        aygpp(6)  = aygpp(6)  + adgpp(6)
        aygpp(7)  = aygpp(7)  + adgpp(7)
        aygpp(8)  = aygpp(8)  + adgpp(8)
        aygpp(9)  = aygpp(9)  + adgpp(9)
        aygpp(10) = aygpp(10)  + adgpp(10)
        aygpp(11) = aygpp(11)  + adgpp(11)
        aygpp(12) = aygpp(12)  + adgpp(12)
c
c gross primary production of the entire gridcell
c
        aygpptot = aygpp(1) + aygpp(2) + aygpp(3) + aygpp(4) + aygpp(5)  + aygpp(6)  + 
     >             aygpp(7) + aygpp(8) + aygpp(9) + aygpp(10) + aygpp(11) + aygpp(12)

cc	  write(100,*)imonth,aygpptot,amgpptot
c
c ---------------------------------------------------------------------
c * * * determine annual npp * * *
c ---------------------------------------------------------------------
c
c net primary production of each plant type
c
        aynpp(1)  = aynpp(1)  + adnpp(1)
        aynpp(2)  = aynpp(2)  + adnpp(2)
        aynpp(3)  = aynpp(3)  + adnpp(3)
        aynpp(4)  = aynpp(4)  + adnpp(4)
        aynpp(5)  = aynpp(5)  + adnpp(5)
        aynpp(6)  = aynpp(6)  + adnpp(6)
        aynpp(7)  = aynpp(7)  + adnpp(7)
        aynpp(8)  = aynpp(8)  + adnpp(8)
        aynpp(9)  = aynpp(9)  + adnpp(9)
        aynpp(10) = aynpp(10)  + adnpp(10)
        aynpp(11) = aynpp(11)  + adnpp(11)
        aynpp(12) = aynpp(12)  + adnpp(12)
c
c gross primary production of the entire gridcell
c
        aynpptot = aynpp(1) + aynpp(2) + aynpp(3) + aynpp(4) + aynpp(5)  + aynpp(6)  + 
     >             aynpp(7) + aynpp(8) + aynpp(9) + aynpp(10) + aynpp(11) + aynpp(12)

	  aylai1  = ((nytimes-1) * aylai1  + lai(1))  * rwork
	  aylai2  = ((nytimes-1) * aylai2  + lai(2))  * rwork

        aylail  = ((nytimes-1) * aylail  + fl * lai(1))  * rwork
	  aylaiu  = ((nytimes-1) * aylaiu  + fu * lai(2))  * rwork

c
c ---------------------------------------------------------------------
c * * * annual carbon budget terms * * *
c ---------------------------------------------------------------------
c
c fire factor used in vegetation dynamics calculations
c
cc        water     = wisoi(1) + (1. - wisoi(1)) * wsoi(1)
cc        waterfrac = (water - swilt(1)) / (1. - swilt(1))
cc        fueldry = max (0.0, min (1.0, -2.0 * (waterfrac - 0.5)))
cc        firefac = ((nytimes-1) * firefac + fueldry) * rwork
c
c increment annual total co2 respiration from microbes
c tco2mic is instantaneous value of co2 flux calculated in biogeochem.f
c
        ayco2mic = ayco2mic + adco2mic
c
c increment annual total co2 respiration from roots
c
        ayco2root = ayco2root + adco2root
c
c calculate annual total co2 respiration from soil
c
        ayco2soi  = ayco2root + ayco2mic
c  
c annual net ecosystem co2 flux -- npp total minus microbial respiration 
c the npp total includes losses from root respiration
c
        ayneetot  = aynpptot - ayco2mic 

cc	write(100,*)"in sumyear 4"
c
c annual average root biomass
c
        allroots = cbior(1)  + cbior(2)  + cbior(3)  +
     >             cbior(4)  + cbior(5)  + cbior(6)  +
     >             cbior(7)  + cbior(8)  + cbior(9)  +
     >             cbior(10) + cbior(11) + cbior(12)
c
        ayrootbio = ((nytimes-1) * ayrootbio + allroots) * rwork
c
c ---------------------------------------------------------------------
c * * * annual biogeochemistry terms * * *
c ---------------------------------------------------------------------
c
c increment annual total of net nitrogen mineralization
c value for tnmin is calculated in biogeochem.f
c
cc        aynmintot = aynmintot +amnmintot
c
c other biogeochemistry variables
c
        ayalit  = ayalit  + totalit
        ayblit  = ayblit  + totrlit
	  aytlit  = aytlit  + totlit	   !Yuan added this

        aycsoi  = ((nytimes-1) * aycsoi  + totcsoi)  * rwork
        ayanlit = ((nytimes-1) * ayanlit + totanlit) * rwork
        aybnlit = ((nytimes-1) * aybnlit + totrnlit) * rwork
        aynsoi  = ((nytimes-1) * aynsoi  + totnsoi)  * rwork

c******* Yuan added more monthly variables ****
        ayclitlm = ((nytimes-1) * ayclitlm + clitlm) * rwork     
        ayclitls = ((nytimes-1) * ayclitls + clitls) * rwork     	  
        ayclitll = ((nytimes-1) * ayclitll + clitll) * rwork
        ayclitrm = ((nytimes-1) * ayclitrm + clitrm) * rwork	
        ayclitrs = ((nytimes-1) * ayclitrs + clitrs) * rwork	
        ayclitrl = ((nytimes-1) * ayclitrl + clitrl) * rwork
        ayclitwm = ((nytimes-1) * ayclitwm + clitwm) * rwork
        ayclitws = ((nytimes-1) * ayclitws + clitws) * rwork	
        ayclitwl = ((nytimes-1) * ayclitwl + clitwl) * rwork
        ayfalll = ((nytimes-1) * ayfalll + falll) * rwork
        ayfallr = ((nytimes-1) * ayfallr + fallr) * rwork
        ayfallw = ((nytimes-1) * ayfallw + fallw) * rwork
        aycsoipas = ((nytimes-1) * aycsoipas + csoipas) * rwork	
        aycsoislop = ((nytimes-1) * aycsoislop + csoislop) * rwork
        aycsoislon = ((nytimes-1) * aycsoislon + csoislon) * rwork 
        aycmic  = ((nytimes-1) * aycmic  + totcmic)  * rwork

c******* end of Yuan added more monthly variables ****/

c
cc 100  continue
c
      return
      end	  