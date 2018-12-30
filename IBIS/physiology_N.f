
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata_N (iyear, scaleu)
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

c
c include water vapor functions
c
      include 'comsat.h'
c
c local variables
c
      integer iyear
c
      real rwork,  ! 3.47e-03 - 1. / tu
     >     tempvm !

      real rh12,   ! relative humidity in upper canopy air 
     >     gbco2u, ! bound. lay. conductance for CO2 in upper canopy
     >     gscuc  !
      real vmax, vmaxuc, rdarkuc, aguc, anuc
      real duma, dumb, dumc, dume, dumq, dump
      real cscuc
      real scaleu
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
cc      real alpha3, alpha4
c

c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
      tau15 = 4500.0     
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
c leaf respiration coefficients
c
      gammauc = 0.0150   ! conifer trees
c
c 'm' coefficients for stomatal conductance relationship
c
      coefmuc = 6.0     ! conifer trees
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
      coefbuc = 0.010    ! conifer trees
c
c absolute minimum stomatal conductances
c
      gsucmin = 0.00001  ! conifer trees
c
c maximum values for ci (for model stability)
c
cc      real cimax
c
      cimax = 2000.e-06  ! maximum values for ci

c
c o2/co2 kinetic parameters (mol/mol)
c
cc       data kc15 /1.5e-04/ 
cc       data ko15 /2.5e-01/

	 kc15 = 1.5e-04
	 ko15 = 2.5e-01
	 alpha3 = 0.060

cc       write(100,*)"stomata_N 5", tau15, gammauc, coefmuc, coefbuc, gsucmin, cimax, kc15, ko15, alpha3
c
c ---------------------------------------------------------------------
c * * * upper canopy physiology calculations * * *
c ---------------------------------------------------------------------
cc        write(100,*)"stomata_N 1", tau15, kc15, ko15, tu, tau, tleaf, ta

        theta3 = 0.970      ! c3 photosynthesis

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
c
c upper canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)

cc        write(100,*)"stomata_N 2"
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
        ciuc = max (1.05 * gamstar, min (cimax, ciuc))

cc        write(100,*)"stomata_N 3", gamstar, ko, ciuc
c
c ---------------------------------------------------------------------
c conifer tree physiology 
c ---------------------------------------------------------------------
c
c vmax and dark respiration for current conditions
c
        vmaxuc = 30.0e-06
cc        vmaxuc = vmaxuc * 2  !yuan

        vmax  = vmaxuc * tempvm * stresstu
        rdarkuc = gammauc * vmaxuc * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu * 4.59e-06 * alpha3 * (ciuc - gamstar) / 
     >       (ciuc + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciuc - gamstar) / 
     >       (ciuc + kc * (1. + o2conc / ko))

cc        write(100,*)"stomata_N 4"
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
        aguc = min (dumq/duma, dumc/dumq) 
        anuc = aguc - rdarkuc

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
        csuc = 0.5 * (csuc + co2conc - anuc / gbco2u)
        csuc = max (1.05 * gamstar, csuc)
c
c calculate new value of gs using implicit scheme
c
        gsuc = 0.5 * (gsuc  +  (coefmuc * anuc * rh12 / csuc + 
     >                                coefbuc * stresstu))
c
        gsuc = max (gsucmin, coefbuc * stresstu, gsuc)
c
c calculate new value of ci using implicit scheme
c
        ciuc = 0.5 * (ciuc + csuc - 1.6 * anuc / gsuc)
        ciuc = max (1.05 * gamstar, min (cimax, ciuc))

cc        write(100,*)"stomata_N 6"
c
c ---------------------------------------------------------------------
c upper canopy scaling
c ---------------------------------------------------------------------
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcuc = aguc * scaleu
        ancuc = anuc * scaleu

c
c calculate diagnostic canopy average surface co2 concentration 
c (big leaf approach)
c
        cscuc = max (1.05 * gamstar, co2conc - ancuc / gbco2u)
c
c calculate diagnostic canopy average stomatal conductance (big leaf approach)
c
        gscuc = coefmuc * ancuc * rh12 / cscuc + coefbuc * stresstu
        gscuc = max (gsucmin, coefbuc * stresstu, gscuc)

cc        write(100,*)"stomata_N 7"
c
c calculate total canopy and boundary-layer total conductance for 
c water vapor diffusion
c
        rwork = 1. / su
        dump  = 1. / 0.029
        totconduc = 1. / (rwork + dump / gscuc)
c
c multiply canopy photosynthesis by wet fraction - this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1 - fwetu
        agcuc = rwork * agcuc
        ancuc = rwork * ancuc

cc	  write(100,*)"stomata_N 8", agcuc, fwetu, aguc, scaleu 

c
      return
      end
c
