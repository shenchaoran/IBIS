c
c ------
c combcs
c ------
c
      real
     >  xintopo,    ! topography (m)
     >  xinveg,     ! fixed vegetation map
     >  deltat      ! absolute minimum temperature -
     >                    ! temp on average of coldest month (C)
c
      common /combcs1/ xintopo, xinveg, deltat
c
cc      integer
cc     >  lmask(nlon,nlat)  ! landmask 0=water, 1=land
c
cc      common /combcs2/ lmask
c
      real
     >  xint(12),    ! climatological temp + anomaly (C)
     >  xinq(12),    ! climatological relative humidity + anomaly (%)
     >  xinprec(12), ! climatological precipition + anomaly (mm/day)
     >  xinwind(12), ! climatological wind speed + anomaly (m s-1)
     >  xincld(12),  ! climatological cloudiness + anomaly(%)
     >  xinwet(12),  ! climatological wet days + anomaly (days/month)
     >  xintrng(12)  ! climatological temp range + anomaly(C)
c
      common /combcs3/ xint, xinq, xinprec, xinwind, xincld, xinwet, xintrng
c
      real
     >  clmt(12),    ! climatological temperature (C)
     >  clmq(12),    ! climatological relative humidity (%)
     >  clmprec(12), ! climatological precipitation (mm/day)
     >  clmw(12),    ! climatological wind speed (m s-1)
     >  clmwet(12),  ! climatological wet days (days/month)
     >  clmcld(12),  ! climatological cloudiness (%)
     >  clmtrng(12)  ! climatological temp range (C)
c
      common /combcs4/ clmt, clmq, clmprec, clmw, clmwet, clmcld, clmtrng
c
      real
     >  xintd,      ! daily climatological temperature (C) +
     >                    !   anomaly (C) + daily anomaly (C)
     >  xinqd,      ! daily climatological relative humidity (%) +
     >                    !   anomaly (%) * daily anomaly (fraction)
     >  xinprecd,   ! daily climatological precipitation (mm/day) +
     >                    !   anomaly (mm/day) * daily anomaly (fraction)
     >  xinwindd,   ! daily climatological windspeed (m/s) +
     >                    !   anomaly (m/s) * daily anomaly (fraction)
     >  xincldd,    ! daily climatological cloud fraction (fraction) +
     >                    !   anomaly (fraction) + daily anomaly (fraction)
     >  xintrngd    ! daily climatological temp range (C) +
     >                    !   anomaly (C) * daily anomaly (C)
c
      common /combcs5/ xintd, xinqd, xinprecd, xinwindd, xincldd, xintrngd
