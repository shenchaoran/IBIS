c
c ------
c comatm
c ------
c      
      real 
     >  coszen,      ! cosine of solar zenith angle
     >  fira         ! incoming ir flux (W m-2)
c
      common /comatm1/ coszen, fira
c
      real 
     >  solad(nband), ! direct downward solar flux (W m-2)
     >  solai(nband), ! diffuse downward solar flux (W m-2)
     >  asurd(nband), ! direct albedo of surface system
     >  asuri(nband)  ! diffuse albedo of surface system 
c
      common /comatm2/ solad, solai, asurd, asuri
c
      real 
     >  ua,          ! wind speed (m s-1)
     >  ta,          ! air temperature (K)
     >  qa,          ! specific humidity (kg_h2o/kg_air)
     >  raina,       ! rainfall rate (mm/s or kg m-2 s-1)
     >  rh,          ! relative humidity(%)
     >  snowa        ! snowfall rate (mm/s or kg m-2 s-1 of water)
c
      common /comatm3/ ua, ta, qa, raina, rh, snowa
c
      real 
     >  psurf,       ! surface pressure (Pa)
     >  cloud,       ! cloud fraction
     >  td,          ! daily average temperature (K)
     >  tmax,        ! maximum daily temperature (K)
     >  tmin,        ! maximum daily temperature (K)
     >  qd,          ! daily average specific humidity (kg_h2o/kg_air)
     >  ud,          ! daily average wind speed (m/sec)
     >  precip,      ! daily precitation (mm/day)
     >  precipday(31),       
     >  precipdaysum       
c
      common /comatm4/ psurf, cloud, td, tmax, tmin, qd, ud, precip,
     >                   precipday, precipdaysum
c
      real 
     >  xstore(3)     ! weather generator 'memory' matrix
c
      common /comatm5/ xstore
c
      integer 
     >  iwet,        ! wet day / dry day flag
     >  iwetday(31), 
     >  iwetdaysum 
c
      common /comatm6/ iwet, iwetday, iwetdaysum
c
      real 
     >  co2conc,           ! co2 concentration (mol/mol)
     >  o2conc             ! o2 concentration (mol/mol)
c
      common /comatm7/ co2conc, o2conc
c
