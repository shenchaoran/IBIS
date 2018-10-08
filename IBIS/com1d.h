c
c -----
c com1d
c -----
c
c canopy scaling terms (Ramankutty, 1998 in progress)
c
      real    
     >  terml(7),     ! term needed in lower canopy scaling
     >  termu(7)      ! term needed in upper canopy scaling
c
      common /com1d1/ termu,   terml
c
      real 
     >  scalcoefl(4), ! term needed in lower canopy scaling
     >  scalcoefu(4)  ! term needed in upper canopy scaling
c
      common /com1d2/ scalcoefl,scalcoefu
c
c variables for internal fluxes: passed to soil and snow models
c
      real  
     >  firb ,       ! net upward ir radiation at reference atmospheric level za (W m-2)
     >  firg ,       ! ir radiation absorbed by soil/ice (W m-2)
     >  firi ,       ! ir radiation absorbed by snow (W m-2)
     >  firl ,       ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
     >  firs ,       ! ir radiation absorbed by upper canopy stems (W m-2)
     >  firu ,       ! ir raditaion absorbed by upper canopy leaves (W m-2)
     >  fsena ,      ! downward sensible heat flux between za & z12 at za (W m-2)
     >  fseng ,      ! upward sensible heat flux between soil surface & air at z34 (W m-2)
     >  fseni,       ! upward sensible heat flux between snow surface & air at z34 (W m-2)
     >  fsenl ,      ! sensible heat flux from lower canopy to air (W m-2)
     >  fsens ,      ! sensible heat flux from upper canopy stems to air (W m-2)
     >  fsenu ,      ! sensible heat flux from upper canopy leaves to air (W m-2)
     >  fvapa ,      ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
     >  fvapg ,      ! h2o vapor flux (evaporation) between soil & air at z34 (kg m-2 s-1/bare ground fraction)
     >  fvapi ,      ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
     >  fvaplt ,     ! h2o vapor flux (transpiration) between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
     >  fvaplw ,     ! h2o vapor flux (evaporation from wet surface) between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
     >  fvaps ,      ! h2o vapor flux (evaporation from wet surface) between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
     >  fvaput ,     ! h2o vapor flux (transpiration from dry parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
     >  fvapuw ,     ! h2o vapor flux (evaporation from wet parts) between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
     >  raing,       ! rainfall rate at soil level (kg m-2 s-1)
     >  rainl,       ! rainfall rate below upper canopy (kg m-2 s-1)
     >  rainu,       ! rainfall rate above upper canopy (kg m-2 s-1)
     >  snowg,       ! snowfall rate at soil level (kg h2o m-2 s-1)
     >  snowl,       ! snowfall rate below upper canopy (kg h2o m-2 s-1)
     >  snowu,       ! snowfall rate above upper canopy (kg h2o m-2 s-1)
     >  soli ,       ! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
     >  solg ,       ! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
     >  soll ,       ! solar flux (direct + diffuse) absorbed by lower canopy leaves and stems per unit canopy area (W m-2)
     >  sols ,       ! solar flux (direct + diffuse) absorbed by upper canopy stems per unit canopy area (W m-2)
     >  solu ,       ! solar flux (direct + diffuse) absorbed by upper canopy leaves per unit canopy area (W m-2)
     >  traing,      ! rainfall temperature at soil level (K)
     >  trainl,      ! rainfall temperature below upper canopy (K)
     >  trainu,      ! rainfall temperature above upper canopy (K)
     >  tsnowg,      ! snowfall temperature at soil level (K) 
     >  tsnowl,      ! snowfall temperature below upper canopy (K)
     >  tsnowu       ! snowfall temperature above upper canopy (K)
c
      common /com1d3/ firb,  firg,  firi,  firl,  firs, firu, fsena, fseng, fseni,
     >     fsenl, fsens, fsenu, fvapa, fvapg, fvapi, fvaplt, fvaplw,
     >     fvaps, fvaput, fvapuw, soli, solg, soll, sols, solu,
     >     raing, rainl, rainu, snowg, snowl, snowu, traing, trainl,
     >     trainu, tsnowg, tsnowl, tsnowu
c
c variables for solar calculations
c
c note that all direct fluxes are per unit horizontal area
c (i.e., already including a factor of cos(zen angle))
c
      integer nsol               ! number of points in indsol
c
      common /com1d4/ nsol
c
      real 
     >  abupd,       ! fraction of direct  radiation absorbed by upper canopy
     >  abupi,       ! fraction of diffuse radiation absorbed by upper canopy
     >  ablod,       ! fraction of direct  radiation absorbed by lower canopy
     >  abloi,       ! fraction of diffuse radiation absorbed by lower canopy
     >  albsnd,      ! direct  albedo for snow surface (visible or IR)
     >  albsni,      ! diffuse albedo for snow surface (visible or IR)
     >  albsod,      ! direct  albedo for soil surface (visible or IR)
     >  albsoi,      ! diffuse albedo for soil surface (visible or IR)
     >  dummy,       ! placeholder, always = 0: no direct flux produced for diffuse incident
     >  flodd,       ! downward direct radiation per unit incident direct radiation on lower canopy (W m-2)
     >  flodi,       ! downward diffuse radiation per unit incident direct radiation on lower canopy (W m-2)
     >  floii,       ! downward diffuse radiation per unit incident diffuse radiation on lower canopy
     >  fupdd,       ! downward direct radiation per unit incident direct beam on upper canopy (W m-2)
     >  fupdi,       ! downward diffuse radiation per unit icident direct radiation on upper canopy (W m-2)
     >  fupii,       ! downward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
     >  relod,       ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
     >  reloi,       ! upward diffuse radiation per unit incident diffuse radiation on lower canopy (W m-2)
     >  reupd,       ! upward direct radiation per unit incident direct radiation on upper canopy (W m-2)
     >  reupi,       ! upward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
     >  sol2d,       ! direct downward radiation  out of upper canopy per unit vegetated (upper) area (W m-2)
     >  sol2i,       ! diffuse downward radiation out of upper canopy per unit vegetated (upper) area(W m-2)
     >  sol3d,       ! direct downward radiation  out of upper canopy + gaps per unit grid cell area (W m-2)
     >  sol3i        ! diffuse downward radiation out of upper canopy + gaps per unit grid cell area (W m-2)
c
      common /com1d5/ abupd, abupi, ablod, abloi, albsnd, albsni, albsod, albsoi,
     >     dummy, flodd, flodi, floii, fupdd, fupdi, fupii, relod,  
     >     reloi, reupd, reupi, sol2d, sol2i, sol3d, sol3i
c
      integer indsol         ! index of current strip for points with positive coszen
c
      common /com1d6/ indsol
c
c variables for aerodynamic calculations
c 
      real      
     >  aloga,       ! log (za - dispu) 
     >  alogav,      ! average of alogi and alogg 
     >  alogg,       ! log of soil roughness
     >  alogi,       ! log of snow roughness
     >  alogl,       ! log (roughness length of lower canopy)
     >  alogu,       ! log (roughness length of upper canopy)
     >  alog1,       ! log (z1 - dispu) 
     >  alog2,       ! log (z2 - displ)
     >  alog3,       ! log (z3 - displ)
     >  alog4,       ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
     >  bdl,         ! aerodynamic coefficient ([(tau/rho)/u**2] for laower canopy (A31/A30 Pollard & Thompson 1995)
     >  bdu,         ! aerodynamic coefficient ([(tau/rho)/u**2] for upper canopy (A31/A30 Pollard & Thompson 1995)
     >  cl,          ! air transfer coefficient (*rhoa) (m s-1 kg m-3) between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
     >  cp,          ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
     >  cu,          ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
     >  dil,         ! inverse of momentum diffusion coefficient within lower canopy (m)
     >  displ,       ! zero-plane displacement height for lower canopy (m)
     >  dispu,       ! zero-plane displacement height for upper canopy (m)
     >  diu,         ! inverse of momentum diffusion coefficient within upper canopy (m)
     >  exphl,       ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
     >  exphu,       ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
     >  expl,        ! exphl**2
     >  expu,        ! exphu**2
     >  pfluxl,      ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
     >  pfluxs,      ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
     >  pfluxu,      ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
     >  rhoa,        ! air density at za (allowing for h2o vapor) (kg m-3)
     >  richl,       ! richardson number for air above upper canopy (z3 to z2)
     >  richu,       ! richardson number for air between upper & lower canopy (z1 to za)
     >  sg,          ! air-soil transfer coefficient
     >  si,          ! air-snow transfer coefficient
     >  strahl,      ! heat/vap correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
     >  strahu,      ! heat/vap correction factor for stratif above upper canopy (z1 to za) (louis et al.)
     >  straml,      ! momentum correction factor for stratif between upper & lower canopy (z3 to z2) (louis et al.)
     >  stramu,      ! momentum correction factor for stratif above upper canopy (z1 to za) (louis et al.)
     >  u1,          ! wind speed at level z1 (m s-1)
     >  u12,         ! wind speed at level z12 (m s-1)
     >  u2,          ! wind speed at level z2 (m s-1)
     >  u3,          ! wind speed at level z3 (m s-1)
     >  u34,         ! wind speed at level z34 (m s-1)
     >  u4,          ! wind speed at level z4 (m s-1)
     >  za,          ! height above the surface of atmospheric forcing (m)
     >  z1,          ! effective top of upper canopy (for momentum) (m)
     >  z12,         ! effective middle of the upper canopy (for momentum) (m)
     >  z2,          ! effective bottom of the upper canopy (for momentum) (m)
     >  z3,          ! effective top of the lower canopy (for momentum) (m)
     >  z34,         ! effective middle of the lower canopy (for momentum) (m)
     >  z4           ! effective bottom of the lower canopy (for momentum) (m)
c
      common /com1d8/ aloga, alogav, alogg, alogi, alogl, alogu, alog1, alog2,
     >     alog3, alog4, bdl, bdu, cl, cp, cu, dil, displ, dispu, diu,
     >     expl, expu,
     >     exphl, exphu, pfluxl, pfluxs, pfluxu, rhoa, richl, richu, sg
     >     ,si,strahl, strahu, straml, stramu, u1, u12, u2, u3, u34, u4
     >     ,za,z1, z12, z2, z3, z34, z4
c
      real tfac            ! (ps/p) ** (rair/cair) for atmospheric level  (const)
c
      common /com1d9/ tfac               
c
c variables for intercepted water
c
      real      
     >  fwetl ,      ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
     >  fwets ,      ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
     >  fwetu ,      ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
     >  fwetlx ,     ! fraction of lower canopy leaf and stem area wetted if dew forms
     >  fwetsx ,     ! fraction of upper canopy stem area wetted if dew forms
     >  fwetux ,     ! fraction of upper canopy leaf area wetted if dew forms
     >  rliql ,      ! proportion of fwetl due to liquid
     >  rliqs ,      ! proportion of fwets due to liquid
     >  rliqu        ! proportion of fwetu due to liquid
c
      common /com1d10/ fwetl, fwets, fwetu, fwetlx, fwetsx, fwetux, rliql, rliqs,
     >     rliqu

c
