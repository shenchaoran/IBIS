c
c ------
c comveg 
c ------
c
      real 
     >  alaiml,             ! lower canopy leaf & stem maximum area (2 sided) for normalization of drag coefficient (m2 m-2)
     >  alaimu,             ! upper canopy leaf & stem area (2 sided) for normalization of drag coefficient (m2 m-2)
     >  cgrass,             ! empirical constant in lower canopy-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  chl,                ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
     >  chs,                ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
     >  chu,                ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
     >  cleaf,              ! empirical constant in upper canopy leaf-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  cstem,              ! empirical constant in upper canopy stem-air aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
     >  tblowl,             ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
     >  tblows,             ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
     >  tblowu,             ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
     >  tdripl,             ! decay time for dripoff of liquid intercepted by lower canopy leaves & stem (sec)
     >  tdrips,             ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
     >  tdripu,             ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
     >  wliqmin,            ! minimum intercepted water on unit vegetated area (kg m-2)
     >  wliqlmax,           ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
     >  wliqsmax,           ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
     >  wliqumax,           ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
     >  wsnomin,            ! minimum intercepted snow on unit vegetated area (kg m-2)
     >  wsnolmax,           ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
     >  wsnosmax,           ! intercepted snow capacity for upper canopy stems (kg m-2)
     >  wsnoumax,           ! intercepted snow capacity for upper canopy leaves (kg m-2)
     >  woodnorm	    ! value of woody biomass for upper canopy closure (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
c
      common /comveg1/ alaiml, alaimu, cgrass, chl, chs, chu, cleaf, cstem, tblowl,
     >     tblows, tblowu, tdripl, tdrips, tdripu, wliqmin, wliqlmax,
     >     wliqsmax, wliqumax, wsnomin, wsnolmax, wsnosmax, wsnoumax, 
     >     woodnorm
c
      real
     >  q12,          ! specific humidity of air at z12
     >  q34,          ! specific humidity of air at z34
     >  sl,           ! air-vegetation transfer coefficients (*rhoa) for lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
     >  ss,           ! air-vegetation transfer coefficients (*rhoa) for upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
     >  su,           ! air-vegetation transfer coefficients (*rhoa) for upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
     >  topparl,      ! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
     >  topparu,      ! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
     >  tl,           ! temperature of lower canopy leaves & stems(K)
     >  ts,           ! temperature of upper canopy stems (K)
     >  tu,           ! temperature of upper canopy leaves (K)
     >  tlsub,        ! temperature of lower canopy vegetation buried by snow (K)
     >  t12,          ! air temperature at z12 (K)
     >  t34,          ! air temperature at z34 (K)
     >  wliql,        ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
     >  wliqs,        ! intercepted liquid h2o on upper canopy stem area (kg m-2)
     >  wliqu,        ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
     >  wsnol,        ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
     >  wsnos,        ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
     >  wsnou         ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
c
      common /comveg2/ q12, q34, sl, ss, su, topparl, topparu, tl, ts, tu, tlsub, 
     >     t12, t34, wliql, wliqs, wliqu, wsnol, wsnos, wsnou
c
c all photosynthesis rates are per unit leaf area
c
      real
     >  agcub,        ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
     >  agcuc,        ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
     >  agcls,        ! canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
     >  agcl3,        ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
     >  agcl4,        ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
     >  ancub,        ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
     >  ancuc,        ! canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
     >  ancls,        ! canopy average net photosynthesis rate - shrubs       (mol_co2 m-2 s-1)
     >  ancl3,        ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
     >  ancl4,        ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
     >  totcondub,     ! 
     >  totconduc,     !
     >  totcondls,     ! 
     >  totcondl3,     !
     >  totcondl4      !
c
      common /comveg3/ agcub, agcuc, agcls, agcl3, agcl4, 
     >     ancub, ancuc, ancls, ancl3, ancl4, 
     >     totcondub, totconduc, totcondls, totcondl3, totcondl4
c
      real
     >  ciub,         ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
     >  ciuc,         ! intercellular co2 concentration - conifer   (mol_co2/mol_air)
     >  cils,         ! intercellular co2 concentration - shrubs    (mol_co2/mol_air)
     >  cil3,         ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
     >  cil4,         ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
     >  csub,         ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
     >  csuc,         ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
     >  csls,         ! leaf boundary layer co2 concentration - shrubs    (mol_co2/mol_air)
     >  csl3,         ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
     >  csl4,         ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
     >  gsub,         ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
     >  gsuc,         ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
     >  gsls,         ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
     >  gsl3,         ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
     >  gsl4          ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
c
      common /comveg4/ ciub, ciuc, cils, cil3, cil4, 
     >     csub, csuc, csls, csl3, csl4,
     >     gsub, gsuc, gsls, gsl3, gsl4
c
      real 
     >  agddl,        ! annual accumulated growing degree days for bud burst, lower canopy (day-degrees)
     >  agddu,        ! annual accumulated growing degree days for bud burst, upper canopy (day-degrees)
     >  fl,           ! fraction of snow-free area covered by lower  canopy
     >  fu,           ! fraction of overall area covered by upper canopy
     >  gdd0,         ! growing degree days > 0C 
     >  gdd0this,     ! annual total growing degree days for current year
     >  gdd5,         ! growing degree days > 5C
     >  gdd5this,     ! annual total growing degree days for current year
     >  sapfrac,      ! fraction of woody biomass that is in sapwood
     >  tc,           ! coldest monthly temperature (C)
     >  tcthis,       ! coldest monthly temperature of current year (C)
     >  tcmin,        ! coldest daily temperature of current year (C)
     >  totlail,      ! total leaf area index for the lower canopy
     >  totlaiu,      ! total leaf area index for the upper canopy
     >  totbiol,      ! total biomass in the lower canopy (kg_C m-2)
     >  totbiou,      ! total biomass in the upper canopy (kg_C m-2)
     >  tw,           ! warmest monthly temperature (C)
     >  twthis,       ! warmest monthly temperature of current year (C)
     >  disturbf,     ! annual fire disturbance regime (m2/m2/yr)
     >  disturbo,     ! fraction of biomass pool lost every year to disturbances other than fire
     >  firefac,      ! factor that respresents the annual average fuel dryness of a grid cell, and hence characterizes the readiness to burn
     >  tco2mic,      ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
     >  tco2root,     ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
     >  tneetot,      ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
     >  tnmin,        ! instantaneous nitrogen mineralization (kg_N m-2/timestep)
     >  tnpptot,      ! instantaneous npp (mol-CO2 / m-2 / second)
     >  tgpptot,      ! instantaneous gpp (mol-CO2 / m-2 / second)
     >  totalit,	    ! total standing aboveground litter (kg_C m-2)
     >  totanlit,	    ! total standing aboveground nitrogen in litter (kg_N m-2)
     >  totcmic,      ! total carbon residing in microbial pools (kg_C m-2)
     >  totcsoi,      ! total carbon in all soil pools (kg_C m-2)
     >  totfall,	    ! total litterfall and root turnover (kg_C m-2/year)
     >  totlit,       ! total carbon in all litter pools (kg_C m-2)
     >  totnlit,      ! total nitrogen in all litter pools (kg_N m-2)
     >  totnmic,      ! total nitrogen residing in microbial pool (kg_N m-2)
     >  totnsoi,      ! total nitrogen in soil (kg_N m-2)
     >  totrlit,      ! total root litter carbon belowground (kg_C m-2)
     >  totrnlit,     ! total root litter nitrogen belowground (kg_N m-2)
     >  tempu,        ! cold-phenology trigger for trees (non-dimensional)
     >  templ,        ! cold-phenology trigger for grasses/shrubs (non-dimensional)
     >  dropu,        ! drought-phenology trigger for trees (non-dimensional)
     >  dropls,       ! drought-phenology trigger for shrubs (non-dimensional)
     >  dropl4,       ! drought-phenology trigger for c4 grasses (non-dimensional)
     >  dropl3,       ! drought-phenology trigger for c3 grasses (non-dimensional)
     >  vegtype0      ! annual vegetation type - ibis classification
c
      common /comveg5/ agddl, agddu, fl, fu, 
     >     gdd0, gdd0this, gdd5, gdd5this, 
     >     sapfrac, tc, tcthis, tcmin, 
     >     totlail, totlaiu, totbiol, totbiou, 
     >     tw, twthis, 
     >     disturbf, disturbo, firefac, 
     >     tco2mic, tco2root, tneetot, tnmin, tnpptot, tgpptot, 
     >     totalit, totanlit, totcmic, totcsoi, totfall, totlit, 
     >     totnlit, totnmic, totnsoi, totrlit, totrnlit, 
     >     tempu, templ, dropu, dropls, dropl4, dropl3, vegtype0
c
      real 
     >  cbiol(npft),   ! carbon in leaf biomass pool (kg_C m-2)
     >  cbior(npft),   ! carbon in fine root biomass pool (kg_C m-2)
     >  cbiow(npft)    ! carbon in woody biomass pool (kg_C m-2)
c
      common /comveg6/ cbiol, cbior, cbiow
c
      real 
     >  clitll,       ! carbon in leaf litter pool - lignin          (kg_C m-2)
     >  clitlm,       ! carbon in leaf litter pool - metabolic       (kg_C m-2)
     >  clitls,       ! carbon in leaf litter pool - structural      (kg_C m-2)
     >  clitrl,       ! carbon in fine root litter pool - lignin     (kg_C m-2)
     >  clitrm,       ! carbon in fine root litter pool - metabolic  (kg_C m-2)
     >  clitrs,       ! carbon in fine root litter pool - structural (kg_C m-2)
     >  clitwl,       ! carbon in woody litter pool - lignin         (kg_C m-2)
     >  clitwm,       ! carbon in woody litter pool - metabolic      (kg_C m-2)
     >  clitws,       ! carbon in woody litter pool - structural     (kg_C m-2)
     >  csoipas,      ! carbon in soil - passive humus               (kg_C m-2)
     >  csoislo,      ! carbon in soil - slow humus                  (kg_C m-2)
     >  csoislon,     ! carbon in soil - slow nonprotected humus     (kg_C m-2)
     >  csoislop,     ! carbon in soil - slow protected humus        (kg_C m-2)
     >  decompl,      ! litter decomposition factor                  (dimensionless)
     >  decomps,      ! soil organic matter decomposition factor     (dimensionless)
     >  falll,        ! annual leaf litter fall                      (kg_C m-2/year)
     >  fallr,        ! annual root litter input                     (kg_C m-2/year)
     >  fallw,        ! annual wood litter fall                      (kg_C m-2/year)
     >  cdisturb      ! annual amount of vegetation carbon lost 
                            ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
c
      common /comveg7/ clitll, clitlm, clitls, clitrl, clitrm, clitrs, clitwl,
     >     clitwm, clitws, csoipas, csoislo, csoislon, csoislop, decompl
     >     , decomps, falll, fallr, fallw, cdisturb
c
      real 
     >  biomass(npft), ! total biomass of each plant functional type  (kg_C m-2)
     >  frac(npft),    ! fraction of canopy occupied by each plant functional type
     >  plai(npft),    ! total leaf area index of each plant functional type
     >  tnpp(npft),    ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
     >  nppdummy(npft), ! canopy NPP before accounting for stem and root respiration
     >  tgpp(npft)     ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
c
      common /comveg8/ biomass, frac, plai, tnpp, nppdummy, tgpp
c
c 1: lower canopy
c 2: upper canopy
c
      real
     >  lai(2),        ! canopy single-sided leaf area index (area leaf/area veg)
     >  sai(2),        ! current single-sided stem area index
     >  zbot(2),       ! height of lowest branches above ground (m)
     >  ztop(2),       ! height of plant top above ground (m)
     >  ztopmx(2)      ! maximum annual height of plant top above ground (m) 
c
      common /comveg9/ lai, sai, zbot, ztop, ztopmx
c
      real
     >  dleaf(2),           ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
     >  dstem(2),           ! typical linear stem dimension in aerodynamic transfer coefficient (m)
     >  orieh(2),           ! fraction of leaf/stems with horizontal orientation
     >  oriev(2)            ! fraction of leaf/stems with vertical
c
      common /comveg10/  dleaf, dstem, orieh, oriev
c
      real 
     >  rhoveg(nband,2),    ! reflectance of an average leaf/stem
     >  tauveg(nband,2)     ! transmittance of an average leaf/stem
c
      common /comveg11/ rhoveg, tauveg
c
      real 
     >  froot(nsoilay,2)    ! fraction of root in soil layer 
c
      common /comveg12/ froot
c
      real 
     >  exist(npft)    ! probability of existence of each plant functional 
c	                            type in a gridcell
c
      common /comveg13/ exist
c
      real 
     >  specla(npft),        ! specific leaf area (m**2/kg) 
     >  aleaf(npft),         ! carbon allocation fraction to leaves
     >  aroot(npft),         ! carbon allocation fraction to fine roots
     >  awood(npft)          ! carbon allocation fraction to wood 
c
      common /comveg14/ specla, aleaf, aroot, awood
c
c
	  real leaf0(npft), stem0(npft), root0(npft)
c
      common /comveg15/ leaf0, root0, stem0
c