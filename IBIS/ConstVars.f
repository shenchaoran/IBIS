$FIXEDFORMLINESIZE:132
c
c ---------------------------------------------------------------------
      subroutine constvars
c ---------------------------------------------------------------------
c
      include 'implicit.h'    
      include 'compar.h'
      include 'comatm.h'
      include 'comhyd.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'soilbgc.h'
      include 'compft.h'
      include 'comtex.h'
c
c set physical constants (mks)
c
      stef  = 5.67e-8 
      vonk  = 0.4
      grav  = 9.80616
      tmelt = 273.16
      hvap  = 2.5104e+6
      hfus  = 0.3336e+6
      hsub  = hvap + hfus
      ch2o  = 4.218e+3
      cice  = 2.106e+3
      cair  = 1.00464e+3
      cvap  = 1.81e+3
      rair  = 287.04
      rvap  = 461.0
      cappa = rair / cair
      rhow  = 1.0e+3
	vzero = 0
c
c specify the epsilon value for the model
c
      epsilon = 1.0e-7
c
      rhos = 0.25 * rhow
c
c consno is thermal conductivity of snow
c
      consno = 0.20
c
c hsnotop is "adaptive-grid" thickness of top snow layer
c
      hsnotop = 0.05
c
c hsnomin is minimum total snow thickness. total thickness
c is constrained to hsnomin for less than 100% cover. (hsnomin
c should be ge nsnolay*hsnotop for vadapt to work properly.)
c
      hsnomin = max (0.15, nsnolay * hsnotop)
c
c fimin and fimax are minimum and maximum snowcover fractions
c
      fimin = 0.00002 * (dtime / 1800.) * (0.15 / hsnomin)
      fimax = 1.000
c
c z0sno is roughness lenth of snow cover
c
      z0sno = 0.0005
c
c
cc      real forganic          ! fraction of organic matter in soil
cc      real 
cc     >     xdat(ndat),        ! % of sand in each textural class
cc     >     ydat(ndat),        ! % of silt in each textural class
cc     >     zdat(ndat)         ! % of clay in each textural class
c
c
c Rawls et al. (1992) soil properties data
c
c      ------------------
c       sand  silt  clay
c      ------------------

      data texdat /
     >  0.92, 0.05, 0.03,  ! sand
     >  0.81, 0.12, 0.07,  ! loamy sand
     >  0.65, 0.25, 0.10,  ! sandy loam
     >  0.42, 0.40, 0.18,  ! loam
     >  0.20, 0.65, 0.15,  ! silt loam
     >  0.60, 0.13, 0.27,  ! sandy clay loam
     >  0.32, 0.34, 0.34,  ! clay loam
     >  0.09, 0.58, 0.33,  ! silty clay loam
     >  0.53, 0.07, 0.40,  ! sandy clay
     >  0.10, 0.45, 0.45,  ! silty clay
     >  0.20, 0.20, 0.60   ! clay
     >  /
c
c porosity (fraction)
c
      data porosdat /
     >  0.437,             ! sand
     >  0.437,             ! loamy sand
     >  0.453,             ! sandy loam
     >  0.463,             ! loam
     >  0.501,             ! silt loam
     >  0.398,             ! sandy clay loam
     >  0.464,             ! clay loam
     >  0.471,             ! silty clay loam
     >  0.430,             ! sandy clay
     >  0.479,             ! silty clay
     >  0.475              ! clay
     >  /
c
c field capacity (fraction)
c
      data sfielddat /
     >  0.091,             ! sand
     >  0.125,             ! loamy sand
     >  0.207,             ! sandy loam
     >  0.270,             ! loam
     >  0.330,             ! silt loam
     >  0.255,             ! sandy clay loam
     >  0.318,             ! clay loam
     >  0.366,             ! silty clay loam
     >  0.339,             ! sandy clay
     >  0.387,             ! silty clay
     >  0.396              ! clay
     >  /
c
c wilting point (fraction)
c
      data swiltdat /
     >  0.033,             ! sand
     >  0.055,             ! loamy sand
     >  0.095,             ! sandy loam
     >  0.117,             ! loam
     >  0.133,             ! silt loam
     >  0.148,             ! sandy clay loam
     >  0.197,             ! clay loam
     >  0.208,             ! silty clay loam
     >  0.239,             ! sandy clay
     >  0.250,             ! silty clay
     >  0.272              ! clay
     >  /
c
c "b" exponent for the Campbell moisture-release equation
c
      data bexdat /
     >  1.7,               ! sand
     >  2.1,               ! loamy sand
     >  3.1,               ! sandy loam
     >  4.5,               ! loam
     >  4.7,               ! silt loam
     >  4.0,               ! sandy clay loam
     >  5.2,               ! clay loam
     >  6.6,               ! silty clay loam
     >  6.0,               ! sandy clay
     >  7.9,               ! silty clay
     >  7.6                ! clay
     >  /
c
c saturated (air entry) potential (m-h2o)
c
      data suctiondat /
     >  0.070,             ! sand
     >  0.090,             ! loamy sand
     >  0.150,             ! sandy loam
     >  0.110,             ! loam
     >  0.210,             ! silt loam
     >  0.280,             ! sandy clay loam
     >  0.260,             ! clay loam
     >  0.330,             ! silty clay loam
     >  0.290,             ! sandy clay
     >  0.340,             ! silty clay
     >  0.370              ! clay
     >  /
c
c saturated hydraulic conductivity (m s-1)
c
      data hydrauldat /
     >  5.8330e-05,        ! sand
     >  1.6972e-05,        ! loamy sand
     >  7.1944e-06,        ! sandy loam
     >  3.6667e-06,        ! loam
     >  1.8889e-06,        ! silt loam
     >  1.1944e-06,        ! sandy clay loam
     >  6.3889e-07,        ! clay loam
     >  4.1667e-07,        ! silty clay loam
     >  3.3333e-07,        ! sandy clay
     >  2.5000e-07,        ! silty clay
     >  1.6667e-07         ! clay
     >  /
c
c set sand/silt/clay vectors (xdat,ydat,zdat) for 11 data points
c
cc      do 100 l = 1, ndat
cc        xdat(l) = texdat(1,l)
cc        ydat(l) = texdat(2,l)
cc        zdat(l) = texdat(3,l)
cc 100  continue
c
c initialization and normalization constant for puddle model (kg m-2)
c
      wpudmax = 4.5
c
cc      if (irestart .eq. 0) then
        wpud =  0.0
        wipud = 0.0
cc      end if
c
c set prescribed soil layer thicknesses
c
      hsoi(1)  = 0.10
      hsoi(2)  = 0.15
      hsoi(3)  = 0.25
      hsoi(4)  = 0.50
      hsoi(5)  = 1.00
      hsoi(6)  = 2.00

c set physical parameters of soil
c
      z0soi = 0.005
c
c leaf litter constants
c
      klm = 0.15 		!dpm leaf --> microbial biomass
      kls = 0.01 		!spm leaf --> microbial biomass
cc      kls = 0.002 		!spm leaf --> microbial biomass
      kll = 0.01		!rpm leaf --> non or protected om
c
c root litter constants
c
      krm = 0.10		!dpm root --> microbial biomass
      krs = 0.005 		!spm root --> microbial biomass
cc      krs = 0.001 		!spm root --> microbial biomass
      krl = 0.005		!rpm root --> non or protected om 
c
c woody litter constants
c
      kwm = 0.001		!dpm wood --> microbial biomass
cc      kwm = 0.0002		!dpm wood --> microbial biomass
      kws = 0.001	 	!spm wood --> microbial biomass
      kwl = 0.001		!rpm wood --> non or protected om 
c
c biomass constants
c
      kbn = 0.045		!biomass --> non protected organic matter 
      kbp = 0.005		!biomass --> protected organic matter
c
c slow and passive c pools
c
      knb = 0.001		!non protected om --> biomass
      kns = 0.000001		!non protected om --> stablized om
      kpb = 0.0001 		!protected om     --> biomass
      kps = 0.000001		!protected om     --> stablized om
      ksb = 8.0e-07		!stablized om     --> biomass
c
c ---------------------------------------------------------------------
c  yield (efficiency) with which microbes gain biomass from c source
c  the rest is driven off as co2 respiration (microbial respiration)
c  all of the respiration produced by microbes is assumed to leave
c  the soil profile over the course of a year
c  taken primarily from the models of Verberne and CENTURY
c ---------------------------------------------------------------------
c
      ylm = 0.4       ! metabolic material efficiencies
      yrm = 0.4
      ywm = 0.4
      yls = 0.3       ! structural efficiencies
      yrs = 0.3
      yws = 0.3
c
      yll = 1.0       ! resistant fraction
      yrl = 1.0 
      ywl = 1.0 
      ybn = 1.0       ! biomass       --> non-protected pool
      ybp = 1.0       ! biomass       --> protected pool
      yps = 1.0       ! protected     --> passive
      yns = 1.0       ! non-protected --> passive

cc      yll = 0.5       ! resistant fraction
cc      yrl = 0.5 
cc      ywl = 0.5 
cc     ybn = 0.5       ! biomass       --> non-protected pool
cc      ybp = 0.5       ! biomass       --> protected pool
cc     yps = 0.5
cc	yns = 0.5
cc
cc      ysb = 0.20       ! passive pool  --> biomass
cc      ypb = 0.20       ! protected     --> biomass
cc      ynb = 0.25       ! non-protected --> biomass
c
c -------------------------------------------------------------------
c split of lignified litter material between protected/non-protected
c slow OM pools
c -------------------------------------------------------------------
c
      lig_frac = 0.50 
c
c -------------------------------------------------------------------
c protected biomass as a fraction of total soil organic carbon
c from Verberne et al., 1990
c -------------------------------------------------------------------
c
      fbsom = 0.017
c
c ---------------------------------------------------------------------
c (effac) --> efficiency of microbial biomass reincorporated
c into biomass pool.(from NCSOIL parameterizations; Molina et al., 1983)
c ---------------------------------------------------------------------
c
       effac = 0.40 
c
c ---------------------------------------------------------------------
c define C:N ratios of substrate pools and biomass
c metabolic, structural, and lignin are for Leaves and roots
c values from Parton et al., 1987 and Whitmore and Parry, 1988
c index: 1 - biomass, 2 - passive pool, 3- slow protected c,
c 4 - slow carbon, non-protected, 5 - resistant, 6 - structural plant
c leaf and root litter, 7 - metabolic plant and root litter, 
c 8- woody biomass
c ---------------------------------------------------------------------
c
       cnr(1)  = 8.0       !c:n ratio of microbial biomass
       cnr(2)  = 15.0      !c:n ratio of passive soil carbon
       cnr(3)  = 10.0      !c:n ratio of protected slow soil carbon
       cnr(4)  = 15.0      !c:n ratio of non-protected slow soil C
       cnr(5)  = 100.0     !c:n ratio of resistant litter lignin
       cnr(6)  = 150.0     !c:n ratio of structural plant litter
       cnr(7)  = 6.0	   !c:n ratio of metabolic plant litter
       cnr(8)  = 250.0     !c:n Ratio of woody components
c
       fmax   = 0.45
       rconst = 1200.0
       cnleaf = 40.0      ! average c:n ratio for leaf litterfall
       cnroot = 60.0      ! average c:n ratio for root turnover
       cnwood = 200.0     ! average c:n ratio for woody debris
c
      data specla  / 25.0,  ! tropical broadleaf evergreen trees
     >               25.0,  ! tropical broadleaf drought-deciduous trees
     >               25.0,  ! warm-temperate broadleaf evergreen trees
     >               12.5,  !12.5,  ! temperate conifer evergreen trees   (xjz)	
     >               25.0,  ! temperate broadleaf cold-deciduous trees
     >               12.5,  !12.5,  ! boreal conifer evergreen trees
     >               25.0,  ! boreal broadleaf cold-deciduous trees  
     >               25.0,  ! boreal conifer cold-deciduous trees
     >               12.5,  !12.5,  ! evergreen shrubs 
     >               25.0,  ! deciduous shrubs 
     >               20.0,  ! warm (c4) grasses
     >               20.0 / ! cool (c3) grasses

      woodnorm = 7.5
c
cc     	real leaf0(12), stem0(12), root0(12)

c set c allocation coefficients for natural vegetation

      aleaf(1)  = 0.30
      aroot(1)  = 0.20
      awood(1)  = 1. - aleaf(1) - aroot(1)
c
      aleaf(2)  = 0.30
      aroot(2)  = 0.20
      awood(2)  = 1. - aleaf(2) - aroot(2)
c
      aleaf(3)  = 0.30
      aroot(3)  = 0.20
      awood(3)  = 1. - aleaf(3) - aroot(3)
c
      aleaf(4)  = 0.30
      aroot(4)  = 0.40
      awood(4)  = 1. - aleaf(4) - aroot(4)

c
      aleaf(5)  = 0.30
      aroot(5)  = 0.20
      awood(5)  = 1. - aleaf(5) - aroot(5)

c
      aleaf(6)  = 0.30
      aroot(6)  = 0.40
      awood(6)  = 1. - aleaf(6) - aroot(6)
c
      aleaf(7)  = 0.30
      aroot(7)  = 0.20
      awood(7)  = 1. - aleaf(7) - aroot(7)

c
      aleaf(8)  = 0.30
      aroot(8)  = 0.20
      awood(8)  = 1. - aleaf(8) - aroot(8)

c
c allocation coefficients for shrubs
c
      aleaf(9)  = 0.45
      aroot(9)  = 0.40
      awood(9)  = 1. - aleaf(9) - aroot(9)

c
      aleaf(10) = 0.45
      aroot(10) = 0.35
      awood(10) = 1. - aleaf(10) - aroot(10)
c
c allocation coefficients for grasses
c
      aleaf(11) = 0.45
      aroot(11) = 0.55
      awood(11) = 0.00

      aleaf(12) = 0.45
      aroot(12) = 0.55
      awood(12) = 0.00

c
c ************************************************************************
c assign some physical properties of vegetation
c ************************************************************************
c
c leaf optical properties were taken from Sellers et al., 1996
c and Bonan, 1995
c
      rhoveg(1,1) = 0.10     ! vis leaf reflectance, lower story
      rhoveg(1,2) = 0.10     ! vis leaf reflectance, upper story 

      rhoveg(2,1) = 0.60     ! nir leaf reflectance, lower story
      rhoveg(2,2) = 0.40     ! nir leaf reflectance, upper story

      tauveg(1,1) = 0.07     ! vis leaf transmittance, lower story
      tauveg(1,2) = 0.05     ! vis leaf transmittance, upper story

      tauveg(2,1) = 0.25     ! nir leaf transmittance, lower story
      tauveg(2,2) = 0.20     ! nir leaf transmittance, upper story

      chiflz = -0.5          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
      chifuz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)

      oriev(1) = max (-chiflz, 0.)
      oriev(2) = max (-chifuz, 0.)

      orieh(1) = max ( chiflz, 0.)
      orieh(2) = max ( chifuz, 0.)

      dleaf(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
      dstem(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization

      dleaf(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
      dstem(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization

      chu = ch2o *  2.0      ! heat capacity of upper leaves
      chl = ch2o *  2.0      ! heat capacity of lower leaves
      chs = ch2o * 50.0      ! heat capacity of stems

      alaimu = 8.0           ! normalization constant for upper canopy aerodynamics
      alaiml = 8.0           ! normalization constant for lower canopy aerodynamics

      cleaf  = 0.01          ! constant in leaf-air aero transfer parameterization
      cgrass = 0.01          ! constant in leaf-air aero transfer parameterization
      cstem  = 0.01          ! constant in leaf-air aero transfer parameterization

      wliqumax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
      wliqsmax = 0.40        ! intercepted water capacity (mm h2o per unit leaf area)
      wliqlmax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)

      wsnoumax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
      wsnosmax = 4.00        ! intercepted snow capacity (mm h2o per unit leaf area)
      wsnolmax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)

      tdripu =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
      tdrips =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
      tdripl =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)

      tblowu = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
      tblows = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
      tblowl = 12.0 * 3600.0 ! decay time for snow blowoff (sec)

      beta1 = 0.950  ! for lower layer herbaceous plants
      beta2 = 0.975  ! for upper layer trees

      return
      end