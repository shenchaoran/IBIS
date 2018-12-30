
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata_C4 (scalel)
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
     >     gscl4   !
      real vmax, vmaxl4, rdarkl4, agl4, anl4
      real duma, dumb, dumc, dume, dumq, dump
      real cscl4, scalel 
	real tleaf, gamstar, tau, esat34, qsat34

      real kco2,   ! initial c4 co2 efficiency (mol-co2/m**2/s)
     >     je,     ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jc,     ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jp,
     >     ji
c
c model parameters
c
c intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
c
c
c leaf respiration coefficients
c
      gammal4 = 0.0300   ! c4 grasses

      alpha4 = 0.050
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
      tau15 = 4500.0     
c
c 'm' coefficients for stomatal conductance relationship
c
      coefml4 = 4.0     ! c4 grasses
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
      coefbl4 = 0.040    ! c4 grasses
c
c absolute minimum stomatal conductances
c
      gsl4min = 0.00001  ! c4 grasses
c
c maximum values for ci (for model stability)
c
      cimax = 2000.e-06  ! maximum values for ci
c
c include water vapor functions
c

	  theta4 = 0.970     ! c4 photosynthesis
        beta4 = 0.800     ! c4 photosynthesis

        cil4 = max (0.0 , min (cimax, cil4))

c ---------------------------------------------------------------------
c * * * lower canopy physiology calculations * * *
c ---------------------------------------------------------------------
c calculate physiological parameter values which are a function of temperature
c
            rwork = 3.47e-03 - 1. / tl
            tau = tau15 * exp(-5000.0 * rwork)
            tleaf = tl - 273.16
            tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * ( 10.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
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
c ---------------------------------------------------------------------
c c4 grass physiology
c ---------------------------------------------------------------------
c
c vmax and dark respiration for current conditions
cc	  write(100,*)"C4 1"
 	  vmaxl4 = 15.0e-06
cc        vmaxl4 = vmaxl4 * 2  !yuan

        vmax  = vmaxl4 * tempvm * stresstl
        rdarkl4 = gammal4 * vmaxl4 * tempvm
c
c initial c4 co2 efficiency (mol/m**2/s)
c
        kco2 = 18.0e+03 * vmax
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
cc        je = topparl * 4.59e-06 * alpha4
        je = topparl * 4.59e-06 * 0.05
c
c 'rubisco limited' rate of photosynthesis
c
        jc = vmax

c
c solve for intermediate photosynthesis rate
c
        duma = theta4
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
        jp = min (dumq/duma, dumc/dumq)

cc	  write(100,*)"C4 2"
c
c 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
c
        ji = kco2 * cil4
c
c solution to quadratic equation
c
        duma = beta4
        dumb = jp + ji
        dumc = jp * ji
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl4 = min (dumq/duma, dumc/dumq)
        anl4 = agl4 - rdarkl4
cc	  write(100,*)"C4 3", gbco2l,csl4
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c CD: For numerical stability (to avoid division by zero in gsl4), 
c csl4 is limited to 1e-8 mol_co2/mol_air.
c  
        csl4 = 0.5 * (csl4 + co2conc - anl4 / gbco2l)
        csl4 = max (1.e-8, csl4)
c
c calculate new value of gs using implicit scheme
c
        gsl4 = 0.5 * (gsl4 + coefml4 * anl4 * rh34 / csl4 +
     >                   coefbl4 * stresstl)
c
        gsl4 = max (gsl4min, coefbl4 * stresstl, gsl4)
cc	  write(100,*)"C4 4"
c
c calculate new value of ci using implicit scheme
c
        cil4 = 0.5 * (cil4 + csl4 - 1.6 * anl4 / gsl4)
        cil4 = max (0.0, min (cimax, cil4))
c
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c perform scaling on all carbon fluxes from upper canopy
c
        
        agcl4 = agl4 * scalel
        ancl4 = anl4 * scalel
c
c calculate canopy average surface co2 concentration
c CD: For numerical stability (to avoid division by zero in gscl4),
c cscl4 is limited to 1e-8 mol_co2/mol_air.
c
        cscl4 = max (1.e-8         , co2conc - ancl4 / gbco2l)
cc	  write(100,*)"C4 5"
c
c calculate canopy average stomatal conductance
c
        gscl4 = coefml4 * ancl4 * rh34 / cscl4 + coefbl4 * stresstl
        gscl4 = max (gsl4min, coefbl4 * stresstl, gscl4)
c
c calculate canopy and boundary-layer total conductance for water vapor diffusion
c
        rwork = 1. / sl
        dump =  1. / 0.029
        totcondl4 = 1. / (rwork + dump / gscl4)
c
c multiply canopy photosynthesis by wet fraction -- this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1. - fwetl
c
        agcl4 = rwork * agcl4
        ancl4 = rwork * ancl4
c
c return to main program
c
      return
      end