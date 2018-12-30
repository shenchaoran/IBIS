
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata_B (iyear,imonth,iday,scaleu)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
      include 'comsum.h'
      include 'compft.h'
      include 'comsat.h'
c
c local variables
c
      integer iyear,imonth,iday
c
      real rwork,  ! 3.47e-03 - 1. / tu
     >     tempvm

      real rh12,   ! relative humidity in upper canopy air 
     >     gbco2u ! bound. lay. conductance for CO2 in upper canopy
      real vmax, vmaxub, rdarkub, agub, anub
      real duma, dumb, dumc, dume, dumq, dump
      real cscub, scaleu, gscub 
	real tau, tleaf, esat12, qsat12 
c
      real kc,     ! co2 kinetic parameter (mol/mol)
     >     ko,     ! o2  kinetic parameter (mol/mol)
     >     je,     ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jc,     ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
     >     gamstar ! gamma*, the co2 compensation points for c3 plants
c
c model parameters
c
c intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
c
cc      data alpha3 /0.060/
      alpha3 = 0.060

c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
      tau15 = 4500.0     
c
c o2/co2 kinetic parameters (mol/mol)
c
      kc15 = 1.5e-04 
      ko15 = 2.5e-01
c
c leaf respiration coefficients
c
c      real gammaub, gammauc, gammals, gammal3, gammal4
c
      gammaub = 0.0150   ! broadleaf trees
c
c 'm' coefficients for stomatal conductance relationship
c
cc      real coefmub, coefmuc, coefmls, coefml3, coefml4
c
      coefmub = 10.0     ! broadleaf trees
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
      coefbub =0.010     ! broadleaf trees
c
c absolute minimum stomatal conductances
c
      gsubmin = 0.00001  ! broadleaf trees
c
c
      cimax = 2000.e-06  ! maximum values for ci
c
	  theta3 = 0.970     ! c3 photosynthesis
	   
	  tu = ta
        rwork = 3.47e-03 - 1. / tu
c
        tau = tau15 * exp(-5000.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1400.0 * rwork)
c
        tleaf = tu - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))

cc          tempvm = (tleaf - 0)*(tleaf - 40)/((tleaf - 0)*(tleaf - 40) - (tleaf - 20)*(tleaf - 20))
c
c upper canopy gamma-star values (mol/mol)
c

cc        write(100,*)tau, kc15, ko15
        gamstar = o2conc / (2. * tau)
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2u = min (10.0, max (0.1, su * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat12 = esat (t12)
        qsat12 = qsat (esat12, psurf)
        rh12   = max (0.30, q12 / qsat12)

c
c constrain ci values to acceptable bounds -- to help ensure numerical stability
c
        ciub = max (1.05 * gamstar, min (cimax, ciub))
c
c ---------------------------------------------------------------------
c broadleaf (evergreen & deciduous) tree physiology 
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c tropical broadleaf trees          60.0 e-06 mol/m**2/sec
c warm-temperate broadleaf trees    40.0 e-06 mol/m**2/sec
c temperate broadleaf trees         25.0 e-06 mol/m**2/sec
c boreal broadleaf trees            25.0 e-06 mol/m**2/sec
c
c
c vmax and dark respiration for current conditions
c
	  vmaxub = 30.0e-06
cc        vmaxub = vmaxub * 2  !yuan

        vmax  = vmaxub * tempvm * stresstu

cc	  write(110,"(5f10.2)")tleaf, tu, rwork, tempvm, stresstu

        rdarkub = gammaub * vmaxub * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu * 4.59e-06 * alpha3 * (ciub - gamstar) / 
     >       (ciub + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciub - gamstar) / 
     >       (ciub + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agub = min (dumq/duma, dumc/dumq)
        anub = agub - rdarkub
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csub = 0.5 * (csub + co2conc - anub / gbco2u)
        csub = max (1.05 * gamstar, csub)
c
c calculate new value of gs using implicit scheme
c
        gsub = 0.5 * (gsub  +  (coefmub * anub * rh12 / csub + 
     >                                coefbub * stresstu))
c
        gsub = max (gsubmin, coefbub * stresstu, gsub)
c
c calculate new value of ci using implicit scheme
c
        ciub = 0.5 * (ciub + csub - 1.6 * anub / gsub)
        ciub = max (1.05 * gamstar, min (cimax, ciub))

c ---------------------------------------------------------------------
c upper canopy scaling
c ---------------------------------------------------------------------
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcub = agub * scaleu
        ancub = anub * scaleu

cc        agcub = agub
cc        ancub = anub
c
c calculate diagnostic canopy average surface co2 concentration 
c (big leaf approach)
c
        cscub = max (1.05 * gamstar, co2conc - ancub / gbco2u)
c
c calculate diagnostic canopy average stomatal conductance (big leaf approach)
c
        gscub = coefmub * ancub * rh12 / cscub + coefbub * stresstu
        gscub = max (gsubmin, coefbub * stresstu, gscub)
c
c calculate total canopy and boundary-layer total conductance for 
c water vapor diffusion
c
        rwork = 1. / su
        dump  = 1. / 0.029
        totcondub = 1. / (rwork + dump / gscub)

c
c multiply canopy photosynthesis by wet fraction - this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1 - fwetu

cc	  if(iyear.eq.2006) then
cc	     write(120,"(3I6,10f10.6)")iyear,imonth,iday,agcub,agub,scaleu,vmax,je,jc,tempvm, stresstu,
cc     >                              dumq/duma, dumc/dumq
cc	  end if

        agcub = rwork * agcub
        ancub = rwork * ancub


      return
      end


c
c ---------------------------------------------------------------------
      subroutine drystress (iyear, imonth, iday)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comveg.h'
c
c local variables
c
      integer k    ! loop indices
	integer iyear, imonth, iday
c
      real stressfac, ! to calculate moisture stress factor 
     >     awc,       ! available water content (fraction)
     >     znorm,     ! normalizing factor
     >     zwilt      ! function of awc, =1 if awc = 1 (no stress)
c
c stressfac determines the 'strength' of the soil moisture
c stress on physiological processes
c
c strictly speaking, stresst* is multiplied to the vmax
c parameters used in the photosynthesis calculations
c
c stressfac determines the shape of the soil moisture response
c
      stressfac = -5.0
c
      znorm = 1.0 - exp(stressfac)
c
cc      do 100 i = 1, npoi
c
c initialize stress parameter
c
        stresstl = 0.0
        stresstu = 0.0
c
c fraction of soil water uptake in each layer
c
        do 110 k = 1, nsoilay
c
c plant available water content (fraction)
c
          awc = min (1.0, max (0.0,
     >              (wsoi(k)*(1 - wisoi(k))   - swilt(k)) /
     >              (sfield(k) - swilt(k))
     >              )         )
c
          zwilt = (1. - exp(stressfac * awc)) / znorm
cc          zwilt = (1. - exp(stressfac * wsoi(k))) / znorm

c
c update for each layer
c
          stressl(k) = froot(k,1) * max (0.0, min (1.0, zwilt))
          stressu(k) = froot(k,2) * max (0.0, min (1.0, zwilt))
c
c integral over rooting profile
c
          stresstl = stresstl + stressl(k)
          stresstu = stresstu + stressu(k)
		
          stresstl = stresstl + froot(k,1)*max(min(1.0, wsoi(k)), swilt(k))
          stresstu = stresstu + froot(k,2)*max(min(1.0, wsoi(k)), swilt(k))
c
 110    continue

cc	    if((iyear.eq.2003).or.(iyear.eq.2008)) then
cc	       write(110,"(3I6, 5f10.4)")iyear, imonth, iday, stresstl, wsoi(1), wisoi(1), swilt(1), sfield(1)
cc	    end if

cc          stresstl = 1
cc          stresstu = 1
c
cc 100  continue
c
c return to main program
c
      return
      end
c
c
