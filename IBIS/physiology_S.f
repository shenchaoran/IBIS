$FIXEDFORMLINESIZE:132
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata_S (scalel)
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
cc      integer i
c
      real rwork,  ! 3.47e-03 - 1. / tu
     >     tempvm !

      real rh34,   ! relative humidity in lower canopy air 
     >     gbco2l, ! bound. lay. conductance for CO2 in lower canopy
     >     gscls  !
      real vmax, vmaxls, rdarkls, agls, anls
      real duma, dumb, dumc, dume, dumq, dump
      real cscls, scalel 
	real tau, tleaf, esat34, qsat34 
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
      gammals = 0.0150   ! shrubs
c
c 'm' coefficients for stomatal conductance relationship
c
      coefmls = 9.0     ! shrubs
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
      coefbls = 0.010    ! shrubs
c
c absolute minimum stomatal conductances
c
      gslsmin = 0.00001  ! shrubs
c
c photosynthesis coupling coefficients (dimensionless)
c
      cimax = 2000.e-06  ! maximum values for ci
c
c include water vapor functions
c
	  
	  write(100,*)"phy 0"

        theta3 = 0.970     ! c3 photosynthesis
c
        cils = max (1.05 * gamstar, min (cimax, cils))
c ---------------------------------------------------------------------
c shrub physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c ---------------------------------------------------------------------
c * * * lower canopy physiology calculations * * *
c ---------------------------------------------------------------------
c calculate physiological parameter values which are a function of temperature
c
        write(100,*)"phy 0.5", tl, tleaf

	  tl = min(max(tl, 253.16), 323.16)

	  tleaf = min(max(tleaf, 253.16), 323.16)

        rwork = 3.47e-03 - 1. / tl
        tau = tau15 * exp(-5000.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1400.0 * rwork)
c
        tleaf = tl - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c lower canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2l = min (10.0, max (0.1, sl * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat34 = esat (t34)
        qsat34 = qsat (esat34, psurf)
        rh34   = max (0.30, q34 / qsat34)
c 
c vmax and dark respiration for current conditions
c

	  vmaxls = 27.5e-06	
cc       vmaxls = vmaxls * 2  !yuan

        vmax  = vmaxls * tempvm * stresstl
        rdarkls = gammals * vmaxls * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl * 4.59e-06 * alpha3 * (cils - gamstar) / 
     >       (cils + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cils - gamstar) / (cils + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
 
        write(100,*)"phy 1"
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agls = min (dumq/duma, dumc/dumq)
        anls = agls - rdarkls
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
        csls = 0.5 * (csls + co2conc - anls / gbco2l)
        csls = max (1.05 * gamstar, csls)
c
c calculate new value of gs using implicit scheme
c
        gsls = 0.5 * (gsls + coefmls * anls * rh34 / csls + coefbls * stresstl)
        gsls = max (gslsmin, coefbls * stresstl, gsls)
c
c calculate new value of ci using implicit scheme
c
        cils = 0.5 * (cils + csls - 1.6 * anls / gsls)
        cils = max (1.05 * gamstar, min (cimax, cils))
c
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcls = agls * scalel
        ancls = anls * scalel

        write(100,*)"phy 2"

c
c calculate canopy average surface co2 concentration
c CD: For numerical stability (to avoid division by zero in gscl4),
c cscl4 is limited to 1e-8 mol_co2/mol_air.
c
        cscls = max (1.05 * gamstar, co2conc - ancls / gbco2l)
c
c calculate canopy average stomatal conductance
c
        gscls = coefmls * ancls * rh34 / cscls + coefbls * stresstl
        gscls = max (gslsmin, coefbls * stresstl, gscls)
c
c calculate canopy and boundary-layer total conductance for water vapor diffusion
c
        rwork = 1. / sl
        dump =  1. / 0.029
        totcondls = 1. / (rwork + dump / gscls)
c
c multiply canopy photosynthesis by wet fraction -- this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1. - fwetl
c
        agcls = rwork * agcls
        ancls = rwork * ancls

        write(100,*)"phy 3"

c
c return to main program
c
      return
      end
