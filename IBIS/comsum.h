c
c ------
c comsum
c ------
c
      integer 
     >  ndtimes,	    ! counter for daily average calculations
     >  nmtimes,            ! counter for monthly average calculations
     >  nytimes             ! counter for yearly average calculations
c
      common /comsum1/ ndtimes, nmtimes, nytimes
c
c daily average fields
c
      real 
     >  adrain,       ! daily average rainfall rate (mm/day)
     >  adsnow,       ! daily average snowfall rate (mm/day)
     >  adaet,        ! daily average aet (mm/day)
     >  adtrans,
     >	adinvap,
     >  adsuvap,
     >  adtrunoff,    ! daily average total runoff (mm/day)
     >  adsrunoff,    ! daily average surface runoff (mm/day)
     >  addrainage,   ! daily average drainage (mm/day)
     >  adrh,         ! daily average rh (percent)
     >  adsnod,       ! daily average snow depth (m)
     >  adsnof,       ! daily average snow fraction (fraction)
     >  adwsoi,       ! daily average soil moisture (fraction)
     >  adwisoi,      ! daily average soil ice (fraction)
     >  adtsoi,       ! daily average soil temperature (c)
     >  adwsoic,      ! daily average soil moisture using root profile weighting (fraction)
     >  adtsoic,      ! daily average soil temperature (c) using profile weighting
     >  adco2mic,     ! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
     >  adco2root,    ! daily accumulated co2 respiration from roots (kg_C m-2 /day)
     >  adco2soi,     ! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
     >  adco2ratio,   ! ratio of root to total co2 respiration
     >  adnmintot,    ! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
     >  adtlaysoi,    ! daily average soil temperature (c) of top layer
     >  adwlaysoi,    ! daily average soil moisture of top layer(fraction)
     >  adneetot,     ! daily accumulated net ecosystem exchange of co2 in ecosystem (kg-C/m**2/day)
     >  adgpptot,     ! daily average GPP (kg_C m-2 /day)
     >  adnpptot

c
      common /comsum2/ adrain, adsnow, adaet, adtrans, adinvap,adsuvap, adtrunoff, adsrunoff, addrainage,
     >     adrh, adsnod, adsnof, adwsoi, adwisoi, adtsoi, adwsoic, adtsoic,
     >     adco2mic, adco2root, adco2soi, adco2ratio, adnmintot,
     >     adtlaysoi, adwlaysoi, adneetot, adgpptot, adnpptot

      real
     >  adnpp(npft),
     >  adgpp(npft)

	  common /comsum3/ adnpp, adgpp

c
c monthly average fields
c
      real
     >  amtemp,       ! monthly average air temperature (C)
     >  amrain,       ! monthly average rainfall rate (mm/day)
     >  amsnow,       ! monthly average snowfall rate (mm/day)
     >  amcloud,      ! monthly average cloudiness (percent)
     >  amrh,         ! monthly average rh (percent)
     >  amqa,         ! monthly average specific humidity (kg-h2o/kg-air)
     >  amaet,        ! monthly average aet (mm/day)
     >  amtrunoff,    ! monthly average total runoff (mm/day)
     >  amsrunoff,    ! monthly average surface runoff (mm/day)
     >  amdrainage,   ! monthly average drainage (mm/day)
     >  amwsoi,       ! monthly average 1m soil moisture (fraction)
     >  amwisoi,      ! monthly average 1m soil ice (fraction)
     >  amvwc,        ! monthly average 1m volumetric water content (fraction)
     >  amawc,        ! monthly average 1m plant-available water content (fraction)
     >  amtsoi,       ! monthly average 1m soil temperature (C)
     >  amsnod,       ! monthly average snow depth (m)
     >  amsnof,       ! monthly average snow fraction (fraction)
     >  amlaiu,       ! monthly average lai for upper canopy (m**2/m**2)
     >  amlail,       ! monthly average lai for lower canopy (m**2/m**2)
     >  amsolar,      ! monthly average incident solar radiation (W/m**2)
     >  amalbedo,     ! monthly average solar albedo (fraction)
     >  amirdown,     ! monthly average downward ir radiation (W/m**2)
     >  amirup,       ! monthly average upward ir radiation (W/m**2)
     >  amsens,       ! monthly average sensible heat flux (W/m**2)
     >  amlatent,     ! monthly average latent heat flux (W/m**2)
     >  amnpptot,     ! monthly total npp for ecosystem (kg-C/m**2/month)
     >  amneetot,	    ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
     >  amco2mic,     ! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
     >  amco2root,    ! monthly total CO2 flux from soil due to root respiration (kg-C/m**2/month)
     >  amco2soi,     ! monthly total soil CO2 flux from microbial and root respiration (kg-C/m**2/month)
     >  amco2ratio,   ! monthly ratio of root to total co2 flux
     >  amnmintot,     ! monthly total N mineralization from microbes (kg-N/m**2/month)
     >  amclitlm,		 ! start of Yuan's 
     >  amclitls,
     >  amclitll,
     >  amclitrm,
     >  amclitrs,
     >  amclitrl,
     >  amclitwm,
     >  amclitws,
     >  amclitwl,
     >  amtotcsoi,
     >  amtotfall,		! end of Yuan's 
     >  amalit,
     >  amblit,
     >  amtlit,
     >  amgpptot
c
      common /comsum4/ amtemp, amrain, amsnow, amcloud, amrh, amqa, amaet,
     >     amtrunoff, amsrunoff, 
     >     amdrainage, amwsoi, amwisoi, amvwc, amawc, amtsoi, amsnod,
     >     amsnof, amlaiu, amlail, amsolar, amalbedo, amirdown, amirup, 
     >     amsens, amlatent, amnpptot, amneetot, amco2mic, amco2root, 
     >     amco2soi, amco2ratio, amnmintot,amclitlm,amclitls,amclitll,amclitrm,
     >     amclitrs,amclitrl,amclitwm,amclitws,amclitwl,
     >     amtotcsoi,amtotfall, amalit, amblit, amtlit, amgpptot
c
      real
     >  amnpp(npft),    ! monthly total npp for each plant type (kg-C/m**2/month)
     >  amgpp(npft)
c
      common /comsum5/ amnpp, amgpp
c
c annual average fields
c
      real
     >  ayprcp,       ! annual average precipitation (mm/yr)
     >  ayaet,        ! annual average aet (mm/yr)
     >  aytrans,      ! annual average transpiration (mm/yr)
     >	ayinvap,
     >  aysuvap,
     >  aytrunoff,    ! annual average total runoff (mm/yr)
     >  aysrunoff,    ! annual average surface runoff (mm/yr)
     >  aydrainage,   ! annual average drainage (mm/yr)
     >  aydwtot,      ! annual average soil+vegetation+snow water recharge (mm/yr or kg_h2o/m**2/yr)
     >  aywsoi,       ! annual average 1m soil moisture (fraction)
     >  aywisoi,      ! annual average 1m soil ice (fraction)
     >  ayvwc,        ! annual average 1m volumetric water content (fraction)
     >  ayawc,        ! annual average 1m plant-available water content (fraction)
     >  aytsoi,       ! annual average 1m soil temperature (C)
     >  ayrratio,     ! annual average runoff ratio (fraction)
     >  aytratio,     ! annual average transpiration ratio (fraction)
     >  aysolar,      ! annual average incident solar radiation (w/m**2)
     >  ayalbedo,     ! annual average solar albedo (fraction)
     >  ayirdown,     ! annual average downward ir radiation (w/m**2)
     >  ayirup,       ! annual average upward ir radiation (w/m**2)
     >  aysens,       ! annual average sensible heat flux (w/m**2)
     >  aylatent,     ! annual average latent heat flux (w/m**2)
     >  aystresstu,   ! annual average soil moisture stress parameter for upper canopy (dimensionless)
     >  aystresstl,   ! annual average soil moisture stress parameter for lower canopy (dimensionless)
     >  ayanpptot,    ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
     >  aynpptot,     ! annual total npp for ecosystem (kg-c/m**2/yr)
     >  aygpptot,     ! annual total gpp for ecosystem (kg-c/m**2/yr)
     >  ayalit,       ! aboveground litter (kg-c/m**2)
     >  ayblit,       ! belowground litter (kg-c/m**2)
     >  aycsoi,       ! total soil carbon (kg-c/m**2)
     >  aycmic,       ! total soil carbon in microbial biomass (kg-c/m**2)
     >  ayanlit,      ! aboveground litter nitrogen (kg-N/m**2)
     >  aybnlit,      ! belowground litter nitrogen (kg-N/m**2)
     >  aynsoi,       ! total soil nitrogen (kg-N/m**2)
     >  ynleach,
     >  ayneetot,     ! annual total NEE for ecosystem (kg-C/m**2/yr)
     >  ayco2mic,     ! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
     >  ayco2root,    ! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
     >  ayco2soi,     ! annual total soil CO2 flux from microbial and root respiration (kg-C/m**2/yr)
     >  aynmintot,    ! annual total nitrogen mineralization (kg-N/m**2/yr)
     >  ayrootbio,     ! annual average live root biomass (kg-C / m**2)
     >  ayclitlm,		 ! start of Yuan's 
     >  ayclitls,
     >  ayclitll,
     >  ayclitrm,
     >  ayclitrs,
     >  ayclitrl,
     >  ayclitwm,
     >  ayclitws,
     >  ayclitwl,
     >  aytlit,
     >  ayfalll,		 
     >  ayfallr,
     >  ayfallw,
     >  aycsoipas,
     >  aycsoislop,
     >  aycsoislon,
     >  aylai1,
     >  aylai2, aylail, aylaiu

c
      common /comsum6/ ayprcp, ayaet, aytrans, ayinvap,aysuvap,aytrunoff, aysrunoff, aydrainage, aydwtot,
     >     aywsoi, aywisoi, ayvwc, ayawc, aytsoi, ayrratio, aytratio,
     >     aysolar, ayalbedo, ayirdown, ayirup, aysens, aylatent,
     >     aystresstu, aystresstl, ayanpptot, aynpptot, aygpptot, 
     >     ayalit, ayblit, aycsoi, aycmic, ayanlit, aybnlit, aynsoi,ynleach,
     >     ayneetot, ayco2mic, ayco2root, ayco2soi, aynmintot, ayrootbio,
     >     ayclitlm, ayclitls, ayclitll, ayclitrm,
     >     ayclitrs,ayclitrl,ayclitwm,ayclitws,
     >     ayclitwl,aytlit,ayfalll,ayfallr,ayfallw,aycsoipas,aycsoislop,
     >	   aycsoislon, aylai1, aylai2, aylail, aylaiu
c
      real 
     >  ayanpp(npft),  ! annual above-ground npp for each plant type(kg-c/m**2/yr)
     >  aynpp(npft),   ! annual total npp for each plant type(kg-c/m**2/yr)
     >  aygpp(npft)    ! annual gross npp for each plant type(kg-c/m**2/yr)
c
      common /comsum7/ ayanpp, aynpp, aygpp
c
c other time average fields
c
      real 
     >  a10td,        ! 10-day average daily air temperature (K)
     >  a10ancub,     ! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
     >  a10ancuc,     ! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
     >  a10ancls,     ! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
     >  a10ancl3,     ! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
     >  a10ancl4,     ! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
     >  a10scalparamu,    ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
     >  a10scalparaml,    ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
     >  a10daylightu,    ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
     >  a10daylightl     ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
c
      common /comsum8/ a10td, a10ancub, a10ancuc, a10ancls, a10ancl3, a10ancl4, 
     >     a10scalparamu, a10scalparaml, a10daylightu, a10daylightl
c
c biogeochem summations
c
      real 
     >  storedn,      ! total storage of N in soil profile (kg_N m-2) 
     >  yrleach       ! annual total amount C leached from soil profile (kg_C m-2/yr)
c
      common /comsum9/ storedn, yrleach
c
