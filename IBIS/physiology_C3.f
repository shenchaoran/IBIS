$FIXEDFORMLINESIZE:132
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata_C3 (scalel, iyear)
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
      real rwork,  ! 3.47e-03 - 1. / tu
     >     tempvm !

      integer iyear

      real rh34,   ! relative humidity in lower canopy air 
     >     gbco2l, ! bound. lay. conductance for CO2 in lower canopy
     >     gscl3  !
      real vmax, vmaxl3, rdarkl3, agl3, anl3
      real duma, dumb, dumc, dume, dumq, dump
      real cscl3, scalel
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
c leaf respiration coefficients
c
c      real gammaub, gammauc, gammals, gammal3, gammal4
c
      gammal3 =0.0150   ! c3 grasses

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
c 'm' coefficients for stomatal conductance relationship
c
      coefml3 = 9.0     ! c3 grasses
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
      coefbl3 = 0.010    ! c3 grasses
c
c absolute minimum stomatal conductances
c
      gsl3min = 0.00001  ! c3 grasses
c
c maximum values for ci (for model stability)
c
      cimax = 2000.e-06  ! maximum values for ci
c
c include water vapor functions
c

	  theta3 = 0.970     ! c3 photosynthesis

        cil3 = max (1.05 * gamstar, min (cimax, cil3))

cc	  write(100,*)"start physiology"
c
c ---------------------------------------------------------------------
c c3 grass physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c ---------------------------------------------------------------------
c * * * lower canopy physiology calculations * * *
c ---------------------------------------------------------------------
c calculate physiological parameter values which are a function of temperature
c
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

cc	  write(100,*)"physiology 1"

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
	  vmaxl3 = 25.0e-06	

cc        vmaxl3 = vmaxl3 * 2  !yuan

        vmax  = vmaxl3 * tempvm * stresstl
        rdarkl3 = gammal3 * vmaxl3 * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl * 4.59e-06 * alpha3 * (cil3 - gamstar) / 
     >       (cil3 + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cil3 - gamstar) / (cil3 + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15

cc	  write(100,*)"physiology 2"

c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl3 = min (dumq/duma, dumc/dumq)
        anl3 = agl3 - rdarkl3
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
        csl3 = 0.5 * (csl3 + co2conc - anl3 / gbco2l)
        csl3 = max (1.05 * gamstar, csl3)
c
c calculate new value of gs using implicit scheme
c
        gsl3 = 0.5 * (gsl3 + coefml3 * anl3 * rh34 / csl3 +
     >                   coefbl3 * stresstl)
c
        gsl3 = max (gsl3min, coefbl3 * stresstl, gsl3)
c
c calculate new value of ci using implicit scheme
c
        cil3 = 0.5 * (cil3 + csl3 - 1.6 * anl3 / gsl3)
        cil3 = max (1.05 * gamstar, min (cimax, cil3))
c
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcl3 = agl3 * scalel
        ancl3 = anl3 * scalel

cc        if((iyear.eq.2003).or.(iyear.eq.2008)) then
cc       	  write(110,*)iyear, agcl3, scalel, stresstl,vmax 
cc	  end if

c
c calculate canopy average surface co2 concentration
c CD: For numerical stability (to avoid division by zero in gscl4),
c cscl4 is limited to 1e-8 mol_co2/mol_air.
c
        cscl3 = max (1.05 * gamstar, co2conc - ancl3 / gbco2l)
c
c calculate canopy average stomatal conductance
c
        gscl3 = coefml3 * ancl3 * rh34 / cscl3 +
     >          coefbl3 * stresstl
c
        gscl3 = max (gsl3min, coefbl3 * stresstl, gscl3)
c
c calculate canopy and boundary-layer total conductance for water vapor diffusion
c
        rwork = 1. / sl
        dump =  1. / 0.029
c
        totcondl3 = 1. / (rwork + dump / gscl3)
c
c multiply canopy photosynthesis by wet fraction -- this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1. - fwetl
c
        agcl3 = rwork * agcl3
        ancl3 = rwork * ancl3
c
      return
      end
c
c
