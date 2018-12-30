c
c -----
c soilbgc
c -----
c
c
c leaf litter constants
c
      real klm, 		!dpm leaf --> microbial biomass
     >kls,	!spm leaf --> microbial biomass
     >kll 		!rpm leaf --> non or protected om

      common /soilbgc1/ klm,kls,kll
c
c root litter constants
c
      real krm, 	!dpm root --> microbial biomass
     > krs, 		!spm root --> microbial biomass
     > krl		!rpm root --> non or protected om
	 
      common /soilbgc2/ krm,krs,krl
c
c woody litter constants
c
      real kwm,		!dpm wood --> microbial biomass
     > kws,	 	    !spm wood --> microbial biomass
     > kwl	        !rpm wood --> non or protected om 

      common /soilbgc3/ kwm,kws,kwl
c
c biomass constants
c
      real kbn,		!biomass --> non protected organic matter 
     > kbp,	!biomass --> protected organic matter
     > knb,	!non protected om --> biomass
     > kns,		!non protected om --> stablized om
     >  kpb, 		!protected om     --> biomass
     > kps,		!protected om     --> stablized om
     > ksb 		!stablized om     --> biomass

      common /soilbgc4/ kbn,kbp,knb,kns,kpb,kps,ksb
c
c ---------------------------------------------------------------------
c  yield (efficiency) with which microbes gain biomass from c source
c  the rest is driven off as co2 respiration (microbial respiration)
c  all of the respiration produced by microbes is assumed to leave
c  the soil profile over the course of a year
c  taken primarily from the models of Verberne and CENTURY
c ---------------------------------------------------------------------
c
      real ylm,      ! metabolic material efficiencies
     >  yrm,
     > ywm,
     > yls,       ! structural efficiencies
     > yrs,
     > yws,
     > yll,     ! resistant fraction
     > yrl,
     > ywl, 
     > ybn,      ! biomass       --> non-protected pool
     > ybp,      ! biomass       --> protected pool
     > yps,       ! protected     --> passive
     > yns,       ! non-protected --> passive
     > ysb,       ! passive pool  --> biomass
     > ypb,      ! protected     --> biomass
     > ynb      ! non-protected --> biomass

      common /soilbgc5/ ylm,yrm,ywm,yls,yrs,yws,
     >	     yll,yrl,ywl,ybn,ybp,yps,yns,ysb,ypb,ynb
c
c -------------------------------------------------------------------
c split of lignified litter material between protected/non-protected
c slow OM pools
c -------------------------------------------------------------------
c
      real lig_frac
      common /soilbgc6/ lig_frac
c
c -------------------------------------------------------------------
c protected biomass as a fraction of total soil organic carbon
c from Verberne et al., 1990
c -------------------------------------------------------------------
c
       real fbsom
       common /soilbgc7/ fbsom

c
c ---------------------------------------------------------------------
c (effac) --> efficiency of microbial biomass reincorporated
c into biomass pool.(from NCSOIL parameterizations; Molina et al., 1983)
c ---------------------------------------------------------------------
c
        real effac
       common /soilbgc8/ effac

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
        real cnr(10)
       common /soilbgc9/ cnr

c
c ---------------------------------------------------------------------
c calculate the fraction of wood, roots and leaves that are structural,
c decomposable, and resistant based on equations presented in Verberne
c model discussion (Geoderma, December 1997 special issue).  fmax is the
c maximum fraction allowed in resistant fraction, rconst is a constant
c defined as 1200.  The cnratio of each plant part has to be less than
c the value of structural defined above (i.e. 150) otherwise the equations
c are unstable...thus the wood litter pool value for cnr(6) is substituted
c with a value higher than that for cnwood (i.e. 250).  this is 
c insignificant for wood since 97% is structural anyways.
c
c ** NOTE ******** 
c Would like to incorporate different C:N ratios of residue/roots for
c different biome types based on literature search
c average c:n ratio would be based on litter inputs from each pft
c ****************
c ---------------------------------------------------------------------
c
c equations were changed on 1-26-99 for erratum in literature (Whitmore
c et al. 1997) which had an error in equations to split litterfall into
c the correct three fractions
c 
c
        real fmax,
     > rconst,
     > cnleaf,      ! average c:n ratio for leaf litterfall
     > cnroot,      ! average c:n ratio for root turnover
     > cnwood    ! average c:n ratio for woody debris
      common /soilbgc11/ fmax,rconst,cnleaf,cnroot,cnwood

