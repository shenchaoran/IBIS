
c    #    #    #     #     #####     #      ##    #
c    #    ##   #     #       #       #     #  #   #
c    #    # #  #     #       #       #    #    #  #
c    #    #  # #     #       #       #    ######  #
c    #    #   ##     #       #       #    #    #  #
c    #    #    #     #       #       #    #    #  ######
c
c ---------------------------------------------------------------------
      subroutine initial (isimveg,irestart,sand0,clay0)
c ---------------------------------------------------------------------
c
      include 'implicit.h'    
c
c Arguments (input)     
c
      integer isimveg,          ! 0 = static veg, 1 = dynamic veg, 
     >                          ! 2 = dynamic veg with cold start
     >        irestart         ! 0 = initial run, 1 = restart run
cc     >     iyrlast              ! last year of previous run (for restart)
      real sand0,
     >     clay0
c
cc      if (irestart .eq. 0) then
cc       write(*,*)"start initial"
        call coldstart
cc      else
cc        call restart (iyrlast)
cc      end if
c
c initialize physical consts, dimensions, unit numbers, lsx model
c
cc	write(30,*)"inisurf"
      call inisurf(irestart)
c
c initialize snow model
c
cc	write(30,*)"inisnow"
      call inisnow
c
c initialize soil model
c
cc	write(30,*)"inisoil"
      call inisoil (irestart,sand0,clay0)
c
c initialize vegetation parameters
c
cc      write(*,*)"inimveg"
      call iniveg (isimveg, irestart)
c
c initialize variables for time averaging
c
cc      write(*,*)"inisum"	
      call inisum
c
c return to main program
c
      return
      end

c
c ---------------------------------------------------------------------
      subroutine inisurf(irestart)
c ---------------------------------------------------------------------
c
c does initialization for model
c
      include 'implicit.h'    
c
      include 'compar.h'
      include 'comatm.h'
      include 'comhyd.h'
      include 'comsum.h'
      include 'comveg.h'
c
c Arguments (input)     
c
      integer irestart          ! 0 = initial run, 1 = restart run
c
c local variables
c
      integer j       ! loop indice
c
c initialize integer variables (can't use const for this)
c
c wet day / dry day flag initialized to dry day (0)
c
cc      do 100 i = 1, npoi 
        iwet = 0
        do 105 j = 1,31
          iwetday(j) = 0
          precipday(j) = 0
 105    continue 
cc 100  continue
c
c zero flux arrays, and global diagnostic arrays
c
      call const (asurd, nband, 0.0)
      call const (asuri, nband, 0.0)
c
      
	totcondub = 0.0
      totconduc = 0.0
      totcondls = 0.0
      totcondl3 = 0.0
      totcondl4 = 0.0   

      ginvap = 0.0
      gsuvap = 0.0
      gtrans = 0.0
      grunof = 0.0
      gdrain = 0.0
c
c initialize vegetation prognostic variables
c
c initialize all temperature fields to 10 degrees C
c

      tu = 283.16
      ts = 283.16
      tl = 283.16

c
c initialize weather generator 'memory'
c
      call const (xstore, 3, 0.0)

c
c initialize temperature of lower canopy buried by
c snow to 0 degrees C
c
      tlsub = 273.16
c
c initialize canopy air conditions (used in turvap)
c
      t12 = 283.16
      t34 = 283.16
      q12 = 0.0
      q34 = 0.0
c
c initialize all co2 concentrations (mol/mol)
c
      ciub = 350.0e-06
      ciuc = 350.0e-06
      cils = 350.0e-06
      cil3 = 350.0e-06
      cil4 = 350.0e-06
      csub = 350.0e-06
      csuc = 350.0e-06
      csls = 350.0e-06
      csl3 = 350.0e-06
	csl4 = 350.0e-06
c
c initialize stomatal conductance (mol-h2o/m**2/sec)
c
      gsub = 0.5
      gsuc = 0.5
      gsls = 0.5
      gsl3 = 0.5 
      gsl4 = 0.5
c
c initialize soil biogeochemistry variables
c
      totlit = 0.0
      totfall = 0.0
      totalit = 0.0
      totrlit = 0.0
      totcsoi =  0.0

      totnlit = 0.0
      totnmic = 0.0
      totanlit = 0.0
      totrnlit = 0.0
      tco2mic =  0.0
      tnpptot =  0.0
      tneetot =  0.0
      tnmin =    0.0
c
c initialize carbon lost to atmosphere due
c to biomass burning
c
	cdisturb = 0.0
c
c initialize phenology flags
c
      if (irestart .eq. 0) then

      tempu = 1.0
      templ = 1.0
      dropu = 1.0
      dropls = 1.0
      dropl4 = 1.0
      dropl3 = 1.0

      end if
c
c initialize water and snow interception fractions
c
      wliqu = 0.0
	wliqs = 0.0
      wliql = 0.0
c
      wsnou = 0.0
      wsnos = 0.0
      wsnol = 0.0
c
      su = 0.0
      ss = 0.0
      sl = 0.0
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisum
c ---------------------------------------------------------------------
c CD
c does initialization for time averaging
c
      include 'implicit.h'    
c
      include 'compar.h'
      include 'comhyd.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c local variables
c
      integer k        ! loop indices
c
c initialize total water content in soil+snow+vegetation (for mass conservation check)
c
cc      do 20 i = 1, npoi
      wtot = (wliqu+wsnou) * fu * 2.0 * lai(2) +
     >       (wliqs+wsnos) * fu * 2.0 * sai(2) +
     >       (wliql+wsnol) * fl * 2.0 *
     >       (lai(1) + sai(1)) * (1. - fi)
c
      wtot = wtot + wpud + wipud
c
        do 10 k = 1, nsoilay
          wtot = wtot +
     >              poros(k)*wsoi(k)*(1.-wisoi(k))*hsoi(k)*rhow +
     >              poros(k)*wisoi(k)*hsoi(k)*rhow
 10     continue
c
        do 20 k = 1, nsnolay
          wtot = wtot + fi*rhos*hsno(k)
 20     continue
c
c Daily means
c
    
      adrain = 0.0
      adsnow = 0.0
      adaet =  0.0
      adtrunoff =  0.0
      adsrunoff =  0.0 
      addrainage = 0.0
      adrh =  0.0
      adsnod = 0.0
      adsnof =  0.0
      adwsoi = 0.0
	adtsoi =  0.0
      adwisoi = 0.0
	adtlaysoi = 0.0
      adwlaysoi = 0.0
      adwsoic = 0.0
      adtsoic = 0.0
	adco2mic = 0.0
      adco2root = 0.0
      decompl = 0.0
      decomps = 0.0
      adnmintot = 0.0
c
c Annual mean quantities
c
      aysolar=   0.0 
      ayirup=    0.0 
      ayirdown=  0.0 
      aysens=    0.0 
      aylatent=  0.0 
      ayprcp=    0.0 
      ayaet=     0.0 
      aytrans=   0.0 
      aytrunoff= 0.0 
      aysrunoff= 0.0 
      aydrainage=0.0 
      aywsoi=    0.0 
      aywisoi=   0.0
      aytsoi=    0.0 
      ayvwc=     0.0
      ayawc=     0.0
      aystresstu=0.0
      aystresstl=0.0
      firefac=   0.0
      ayco2mic=  0.0 
      ayco2root= 0.0 
      ayrootbio= 0.0 
      aynmintot= 0.0 
      ayalit=    0.0
      ayblit=    0.0
      aycsoi=    0.0
      aycmic=    0.0
      ayanlit=   0.0
      aybnlit=   0.0
      aynsoi=    0.0
	aytlit=  0.0

      call const (aygpp,      npft, 0.0)     
      call const (aynpp,      npft, 0.0) 

      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisnow
c ---------------------------------------------------------------------
c
c does initialization for snow model
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsno.h'
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
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisoil (irestart,sand0,clay0)
c ---------------------------------------------------------------------
c
c does initialization for soil database
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
      include 'comtex.h'
	include 'soilbgc.h'
c
c Arguments (input)     
c
      integer irestart          ! 0 = initial run, 1 = restart run
c
c local variables
c
      integer k,
     >        l,
cc     >        ndat,           ! number of textural classes
     >        msand,          ! % of sand in grid point
     >        mclay,          ! % of clay in grid point
     >        lmin,           ! closest textural class from texture of 
     >                        ! grid point
     >        textcls         ! textural class assignment (1..11)
c
      real fsand,             ! fraction of sand in grid point
     >     fsilt,             ! fraction of silt in grid point
     >     fclay,             ! fraction of clay in grid point
     >     forganic          ! fraction of organic matter in soil
* M. El Maayar modified this....
cc     >     dmin               ! small number

cc      parameter (ndat=11)
c
      real 
     >     xdat(ndat),        ! % of sand in each textural class
     >     ydat(ndat),        ! % of silt in each textural class
     >     zdat(ndat)         ! % of clay in each textural class
	 real sand0,
     > 	  clay0
c
c set sand/silt/clay vectors (xdat,ydat,zdat) for 11 data points
c
      do 100 l = 1, ndat
        xdat(l) = texdat(1,l)
        ydat(l) = texdat(2,l)
        zdat(l) = texdat(3,l)
 100  continue
c
c initialization and normalization constant for puddle model (kg m-2)
c
      wpudmax = 4.5
c
      if (irestart .eq. 0) then
        wpud =  0.0
        wipud = 0.0
      end if
c
c initialize soil water and soil temperature fields
c
      if (irestart .eq. 0) then
c
        call const (wsoi,  nsoilay, 0.50  )
        call const (wisoi, nsoilay, 0.00  )
c
        call const (tsoi,  nsoilay, 278.13)
c
        tg = 278.13
        ti = 273.13
      else
cc         do 150 i = 1, npoi
            tg = tsoi(1)
            ti = tsno(1)
cc 150     continue
c
      end if
c
c set soil surface parameters for the global domain
c
cc      do 200 i = 1, npoi
c
c Convert input sand and clay percents to fractions
c
        msand = nint(sand0)
        mclay = nint(clay0) 
c
        fsand = 0.01 * msand
        fclay = 0.01 * mclay
        fsilt = 0.01 * (100 - msand - mclay)
* M. El Maayar modified this.
c
c soil surface albedo:
c
c from bats table 3.ii assuming albedo depends on texture
c
        albsav = fsand * 0.120 +
     >              fsilt * 0.085 +
     >              fclay * 0.050
c
c
* M. El Maayar modified this.
*      if (nint(forganic).eq.1) then
*        albsan = 1.0
*      else
        albsan = 2.0 * albsav
*      endif
c
cc 200  continue   
c
c create soil properties look-up table
c
c set soil parameters at each layer for the global domain
c soita.nc file is for layers only to 4 m; currently this
c is for a total of six layers. If there are any remaining  
c layers below that, set texture to be equal to that of the
c last layer (layer 6)
c analysis of the current WISE-IGBP soil textural dataset
c reveals very little information below 4 m.
c
      do 300 k = 1, nsoilay 
cc        do 310 i = 1, npoi 
c
c Convert input sand and clay percents to fractions
c
          if (k.le.6) then
	      sand(k) = sand0
	      clay(k) = clay0

            msand = nint(sand(k))
            mclay = nint(clay(k)) 
          else
            msand = nint(sand(6)) 
            mclay = nint(clay(6)) 
          endif
c
* M. El Maayar modified this
*       if ((msand.ge.99).AND.(mclay.ge.99)) then
*          fsand = 0.
*          fclay = 0.
*          fsilt = 0.
*          forganic = 1.
*       else
          fsand = 0.01 * msand
          fclay = 0.01 * mclay
          fsilt = 0.01 * (100 - msand - mclay)
c
c for now, we assume that all soils have a 1% organic content -- 
c this is just a place holder until we couple the soil carbon
c dynamics to the soil physical properties
c
          forganic = 0.010
*       endif
c
c density of soil material (without pores, not bulk) (kg m-3)
c from Campbell and Norman, 1998
c
          rhosoi(k) = 2650.0 * (1.0 - forganic) + 1300.0 * forganic 
c
c specific heat of soil material (j kg-1 k-1):
c from Campbell and Norman, 1998
c
          csoi(k) =  870.0 * (1.0 - forganic) + 1920.0 * forganic 
c
c cjk
c match textural fractions with soil textural class 
c calls two functions to match sand and clay fractions
c with proper soil textural class based on the usda
c classification system
c
          lmin = textcls (msand,mclay)
c
c porosity (fraction):
c 
          poros(k) = porosdat(lmin)
c
c field capacity (defined relative to the porosity):
c
          sfield(k) = 1.0 / poros(k) * sfielddat(lmin)
c
c wilting point (defined relative to the porosity):
c
          swilt(k)  = 1.0 / poros(k) * swiltdat(lmin)
c
c "b" exponent for the Campbell moisture-release equation:
c
          bex(k) = bexdat(lmin)
c
c nearest integer of "b" exponent (for computational efficiency):
c
          ibex(k) = nint(bex(k))
c
c saturated matric (air entry) potential (m-h2o):
c
          suction(k) = suctiondat(lmin)
c
c saturated hydraulic conductivity (m s-1):
c
          hydraul(k) = hydrauldat(lmin)
c
 300  continue

      return
      end

c-------------------------------------------------------------------------
      integer function textcls (msand,mclay)
c
c adapted for ibis by cjk 01/11/01
c-------------------------------------------------------------------------
c |
c |                         T R I A N G L E
c | Main program that calls WHAT_TEXTURE, a function that classifies soil
c | in the USDA textural triangle using sand and clay %
c +-----------------------------------------------------------------------
c | Created by: aris gerakis, apr. 98 with help from brian baer
c | Modified by: aris gerakis, july 99: now all borderline cases are valid
c | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
c |              main program
c +-----------------------------------------------------------------------
c | COMMENTS
c | o Supply a data file with two columns, in free format:  1st column sand,
c |   2nd column clay %, no header.  The output is a file with the classes.
c +-----------------------------------------------------------------------
c | You may use, distribute and modify this code provided you maintain
c ! this header and give appropriate credit.
c +-----------------------------------------------------------------------
c
c code adapted for IBIS by cjk 01-11-01
c
c     include 'compar.h'
c     include 'comsoi.h'
c
      integer  msand,
     >         mclay
c
      logical inpoly
c
      real    silty_loam(1:7,1:2),
     >        sandy(1:7,1:2),
     >        silty_clay_loam(1:7,1:2), 
     >        loam(1:7,1:2),
     >        clay_loam(1:7,1:2),
     >        sandy_loam(1:7,1:2),
     >        silty_clay(1:7,1:2),
     >        sandy_clay_loam(1:7,1:2), 
     >        loamy_sand(1:7,1:2),
     >        clayey(1:7,1:2),
c    >        silt(1:7,1:2), 
     >        sandy_clay(1:7,1:2)
c
c initalize polygon coordinates:
c each textural class reads in the sand coordinates (1,7) first, and
c then the corresponding clay coordinates (1,7)

c     data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
c
c because we do not have a separate silt category, have to redefine the
c polygon boundaries for the silt loam  
c
      data sandy        /85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
      data loamy_sand   /70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
      data sandy_loam   /50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 
     >	  15, 0/
      data loam       /43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
      data silty_loam     /0, 0, 23, 50, 0, 0, 0, 0, 27, 27, 0, 0, 0, 0/ 
c     data silt            /0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
      data sandy_clay_loam /52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 
     >	 20, 0, 0/
      data clay_loam  /20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data silty_clay_loam /0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 
     > 	0, 0/
      data sandy_clay   /45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/
      data silty_clay   /0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
      data clayey    /20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
c
c +-----------------------------------------------------------------------
c | figure out what texture grid cell and layer are part of  
c | classify a soil in the triangle based on sand and clay %
c +-----------------------------------------------------------------------
c | Created by: aris gerakis, apr. 98
c | Modified by: aris gerakis, june 99.  Now check all polygons instead of
c | stopping when a right solution is found.  This to cover all borderline 
c | cases.
c +-----------------------------------------------------------------------
c
c find polygon(s) where the point is.  
c
      textcls = 0 
c
      if (msand .gt. 0.0 .and. mclay .gt. 0.0) then
         if (inpoly(sandy, 3, msand, mclay)) then
            textcls = 1      ! sand
         endif
         if (inpoly(loamy_sand, 4, msand, mclay)) then
            textcls = 2      ! loamy sand
         endif
         if (inpoly(sandy_loam, 7, msand, mclay)) then
            textcls = 3      ! sandy loam
         endif
         if (inpoly(loam, 5, msand, mclay)) then
            textcls = 4      ! loam
         endif
         if (inpoly(silty_loam, 4, msand, mclay)) then
            textcls = 5      ! silt loam
         endif
         if (inpoly(sandy_clay_loam, 5, msand, mclay)) then
            textcls = 6      ! sandy clay loam
         endif
         if (inpoly(clay_loam, 4, msand, mclay)) then
            textcls = 7      ! clay loam
         endif
         if (inpoly(silty_clay_loam, 4, msand, mclay)) then
            textcls = 8      ! silty clay loam
         endif
         if (inpoly(sandy_clay, 3, msand, mclay)) then
            textcls = 9      ! sandy clay
         endif
         if (inpoly(silty_clay, 3, msand, mclay)) then
            textcls = 10     ! silty clay
         endif
         if (inpoly(clayey, 5, msand, mclay)) then
            textcls = 11     ! clay
         endif
      endif
c
      if (textcls .eq. 0) then
         textcls = 5         ! silt loam
c
c        write (*, 1000) msand, mclay
c 1000   format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
      endif
c
      return
      end
c
c---------------------------------------------------------------------------
      logical function inpoly (poly, npoints, xt, yt)
c
c adapted for ibis by cjk 01/11/01
c---------------------------------------------------------------------------
c
c                            INPOLY
c   Function to tell if a point is inside a polygon or not.
c--------------------------------------------------------------------------
c   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
c
c   Please feel free to use this source code for any purpose, commercial
c   or otherwise, as long as you don't restrict anyone else's use of
c   this source code.  Please give credit where credit is due.
c
c   Point-in-polygon algorithm, created especially for World-Wide Web
c   servers to process image maps with mouse-clickable regions.
c
c   Home for this file:  http://www.gcomm.com/develop/inpoly.c
c
c                                       6/19/95 - Bob Stein & Craig Yap
c                                       stein@gcomm.com
c                                       craig@cse.fau.edu
c--------------------------------------------------------------------------
c   Modified by:
c   Aris Gerakis, apr. 1998: 1.  translated to Fortran
c                            2.  made it work with real coordinates
c                            3.  now resolves the case where point falls
c                                on polygon border.
c   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
c   Aris Gerakis, july 1999: Now all borderline cases are valid
c--------------------------------------------------------------------------
c   Glossary:
c   function inpoly: true=inside, false=outside (is target point inside
c                    a 2D polygon?)
c   poly(*,2):  polygon points, [0]=x, [1]=y
c   npoints: number of points in polygon
c   xt: x (horizontal) of target point
c   yt: y (vertical) of target point
c--------------------------------------------------------------------------
c
c declare arguments  
c
      integer  npoints,
     >         xt,
     >         yt 
c
      real     poly(7, 2)
c
c local variables
c
      real     xnew,
     >         ynew,
     >         xold,
     >         yold,
     >         x1,
     >         y1,
     >         x2,
     >         y2
c
      integer  i
c
      logical inside,
     >        on_border

      inside = .false.
      on_border = .false.
c
      if (npoints .lt. 3)  then
        inpoly = .false.
        return
      end if
c
      xold = poly(npoints,1)
      yold = poly(npoints,2)

      do 300  i = 1 , npoints
        xnew = poly(i,1)
        ynew = poly(i,2)

        if (xnew .gt. xold)  then
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
        else
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
        end if

c the outer IF is the 'straddle' test and the 'vertical border' test.
c the inner IF is the 'non-vertical border' test and the 'north' test.  

c the first statement checks whether a north pointing vector crosses  
c (stradles) the straight segment.  There are two possibilities, depe-
c nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
c must be "open" at left, which is necessary to keep correct count when 
c vector 'licks' a vertix of a polygon.  

        if ((xnew .lt. xt .and. xt .le. xold) .or. (.not. xnew .lt.  
     >     xt .and. .not. xt .le. xold)) then
c
c the test point lies on a non-vertical border:
c
          if ((yt-y1)*(x2-x1) .eq. (y2-y1)*(xt-x1)) then
            on_border = .true. 
c
c check if segment is north of test point.  If yes, reverse the 
c value of INSIDE.  The +0.001 was necessary to avoid errors due   
c arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
c
          elseif ((yt-y1)*(x2-x1) .lt. (y2-y1)*(xt-x1) + 0.001) then
            inside = .not.inside ! cross a segment
          endif
c
c this is the rare case when test point falls on vertical border or  
c left edge of non-vertical border. The left x-coordinate must be  
c common.  The slope requirement must be met, but also point must be
c between the lower and upper y-coordinate of border segment.  There 
c are two possibilities,  depending on whether ynew < yold or ynew > 
c yold:
c
        elseif ((xnew .eq. xt .or. xold .eq. xt) .and. (yt-y1)*(x2-x1)  
     > .eq.(y2-y1)*(xt-x1) .and. ((ynew .le. yt .and. yt .le. yold) .or. 
     >    (.not. ynew .lt. yt .and. .not. yt .lt. yold))) then
          on_border = .true. 
        endif
c
        xold = xnew
        yold = ynew
c
 300    continue  
c
c If test point is not on a border, the function result is the last state 
c of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
c inside the polygon if it falls on any of its borders:
c
      if (.not. on_border) then
         inpoly = inside
      else
         inpoly = .true.
      endif
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine iniveg (isimveg, irestart)
c ---------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'compft.h'
c
c Arguments (input)
c
      integer irestart,     ! 0: not a restart run 1: restart run
     >        isimveg
c
c local variables
c
      integer ideci,        ! # deciduous plant functional types (pft)
     >        ievgr,        ! # evergreen pft 
     >        ishrub,       ! # shrub pft 
     >        igrass,       ! # herbaceous pft 
     >        ilower,       ! possible # pft for lower canopy
     >        iupper,       ! possible # pft for upper canopy
     >        j,k         ! loop indices
c
      real plaievgr,        ! potential lai of evergreen trees
     >     plaideci,        ! potential lai of deciduous trees
     >     plaishrub,       ! potential lai of shrubs
     >     plaigrass,       ! potential lai of grasses
     >     wood,            ! total wood biomas in grid cell
     >     totdepth,        ! total soil depth
     >     frootnorm1,      ! normalization factor for Jackson rooting profile,low
     >     frootnorm2       ! normalization factor for Jackson rooting profile, up
     
      real depth(nsoilay)   ! soil layer depth (cm)
c
cc      do 100 i = 1, npoi
c
c initialize a few climatic variables needed for vegetation
c
cc         if (irestart .eq. 0) then
cc           agddu = 1000.0
cc           agddl = 1000.0
cc         end if
c
c initialize the moisture stress factors
c
        stresstu = 1.0
        stresstl = 1.0
c
c initialize running-mean air temperature
c
cc        if (irestart .eq. 0) then
c
           a10td = 273.16
c
c initialize running-mean values of canopy photosynthesis rates
c
           a10ancub = 10.0e-06
           a10ancuc = 10.0e-06
           a10ancls = 10.0e-06
           a10ancl4 = 10.0e-06
           a10ancl3 = 10.0e-06
c 
c initialize running-mean values of the scaling parameter
c
           a10scalparamu = 0.5 * 5.
           a10scalparaml = 0.5 * 5.
           a10daylightu = 5.
           a10daylightl = 5.
c
c initialize litter fall
c
cc           falll = 0.1
cc           fallr = 0.1
cc           fallw = 0.1
c
cc        end if
      

c
c reset counters
c
        ievgr  = 0
        ideci  = 0
        ishrub = 0
        igrass = 0
        ilower = 0
        iupper = 0
c
c determine number of evergreen plant functional types
c
        if (nint(exist(1)).eq.1) ievgr = ievgr + 1
        if (nint(exist(3)).eq.1) ievgr = ievgr + 1
        if (nint(exist(4)).eq.1) ievgr = ievgr + 1
        if (nint(exist(6)).eq.1) ievgr = ievgr + 1
c
c determine number of deciduous plant functional types
c
        if (nint(exist(2)).eq.1) ideci = ideci + 1
        if (nint(exist(5)).eq.1) ideci = ideci + 1
        if (nint(exist(7)).eq.1) ideci = ideci + 1
        if (nint(exist(8)).eq.1) ideci = ideci + 1
c
c make sure counter is at least 1 (to avoid division by zero)
c
        ievgr = max (1, ievgr)
        ideci = max (1, ideci)
c
c determine number of shrub functional types
c
        if (nint(exist(9)).eq.1)  ishrub = ishrub + 1
        if (nint(exist(10)).eq.1) ishrub = ishrub + 1
c
c determine number of herbaceous plant functional types
c
        if (nint(exist(11)).eq.1) igrass = igrass + 1
        if (nint(exist(12)).eq.1) igrass = igrass + 1
c
c make sure counter is at least 1 (to avoid division by zero)
c
        ishrub = max (1, ishrub)
        igrass = max (1, igrass)
c
c total number of possible pfts for each canopy
c
        iupper = ievgr  + ideci
        ilower = ishrub + igrass 
c
c make sure counter is at least 1 (to avoid division by zero)
c
        iupper = max (1, iupper)
        ilower = max (1, ilower)
c
c for initialization purposes, set the predicted vegetation type
c to the initial vegetation type
c
c ---------------------------------------------------
c  1: tropical evergreen forest / woodland
c  2: tropical deciduous forest / woodland
c  3: temperate evergreen broadleaf forest / woodland
c  4: temperate evergreen conifer forest / woodland
c  5: temperate deciduous forest / woodland
c  6: boreal evergreen forest / woodland
c  7: boreal deciduous forest / woodland
c  8: mixed forest / woodland
c  9: savanna
c 10: grassland / steppe
c 11: dense shrubland
c 12: open shrubland
c 13: tundra
c 14: desert
c 15: polar desert / rock / ice
c ---------------------------------------------------
c
c these classes consist of some combination of 
c plant functional types:
c
c ---------------------------------------------------
c  1: tropical broadleaf evergreen trees
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen trees
c  4: temperate conifer evergreen trees
c  5: temperate broadleaf cold-deciduous trees
c  6: boreal conifer evergreen trees
c  7: boreal broadleaf cold-deciduous trees
c  8: boreal conifer cold-deciduous trees
c  9: evergreen shrubs
c 10: cold-deciduous shrubs
c 11: warm (c4) grasses
c 12: cool (c3) grasses
c ---------------------------------------------------

**** DTP 2001/05/25. The following code replaces the 450+
*    lines of stuff that follows it (hence the temporary goto
*    statement). Note that values of plai_init are read in as
*    parameters from params.veg. Note also that the declarations
*    of the four local variables plaievgr, plaideci, plaishrub 
*    and plaigrass can all be dropped.

cc            plai(1)  = exist(1)  / float(ievgr)  * plai_init(1,inveg)
cc            plai(2)  = exist(2)  / float(ideci)  * plai_init(2,inveg)
cc            plai(3)  = exist(3)  / float(ievgr)  * plai_init(1,inveg)
cc            plai(4)  = exist(4)  / float(ievgr)  * plai_init(1,inveg)
cc            plai(5)  = exist(5)  / float(ideci)  * plai_init(2,inveg)
cc            plai(6)  = exist(6)  / float(ievgr)  * plai_init(1,inveg)
cc            plai(7)  = exist(7)  / float(ideci)  * plai_init(2,inveg)
cc            plai(8)  = exist(8)  / float(ideci)  * plai_init(2,inveg)
cc            plai(9)  = exist(9)  / float(ishrub) * plai_init(3,inveg)
cc            plai(10) = exist(10) / float(ishrub) * plai_init(3,inveg)
*            plai(11) = exist(11) / float(igrass) * plai_init(4,inveg)
*            plai(12) = exist(12) / float(igrass) * plai_init(4,inveg)
cc            if ((inveg.eq.9).or.(inveg.eq.10)) then
cc              if (tw.gt.22.0) then
cc                plai(11) = exist(11) * 0.80 * plai_init(4,inveg)
cc                plai(12) = exist(12) * 0.20 * plai_init(4,inveg)
cc              else
cc                plai(11) = exist(11) * 0.00 * plai_init(4,inveg)
cc                plai(12) = exist(12) * 1.00 * plai_init(4,inveg)
cc              endif
cc            else
cc              plai(11) = exist(11) / float(igrass) * plai_init(4,inveg)
cc              plai(12) = exist(12) / float(igrass) * plai_init(4,inveg)
cc            endif
c
c set minimum lai for each existing plant type
c
          xminlai = 0.010
c
cc          plai(1)  = max (plai(1) , exist(1)  * xminlai)
cc          plai(2)  = max (plai(2) , exist(2)  * xminlai)
cc          plai(3)  = max (plai(3) , exist(3)  * xminlai)
cc          plai(4)  = max (plai(4) , exist(4)  * xminlai)
cc          plai(5)  = max (plai(5) , exist(5)  * xminlai)
cc          plai(6)  = max (plai(6) , exist(6)  * xminlai)
cc          plai(7)  = max (plai(7) , exist(7)  * xminlai)
cc          plai(8)  = max (plai(8) , exist(8)  * xminlai)
cc          plai(9)  = max (plai(9) , exist(9)  * xminlai)
cc          plai(10) = max (plai(10), exist(10) * xminlai)
cc          plai(11) = max (plai(11), exist(11) * xminlai)
cc          plai(12) = max (plai(12), exist(12) * xminlai)
c
c set sapwood fraction and biomass characteristics
c
          sapfrac = sapfrac_init ! 0.1 from params.veg
c
          wood = 0.0 
c
          do 120 j = 1, npft
            cbiol(j) = max (exist(j) * xminlai / specla(j),cbiol(j))
            plai(j)    = cbiol(j) * specla(j)

cc            cbiol(j) = plai(j) / specla(j)
cc            cbior(j) = 0.5 * cbiol(j)
c
cc            cbiow(j) = 0.0
c
cc            if (j.lt.9) cbiow(j) = plai(j) * 10.0 / 6.0


            biomass(j) = cbiol(j) + cbiow(j) + cbior(j) 
            wood = wood + cbiow(j)
 120      continue
c
c ************************************************************************
c determine basic vegetation structure characteristics
c ************************************************************************
c 
c total leaf area for upper and lower canopies
c
        totlaiu  =  plai(1) + plai(2) +
     >                 plai(3) + plai(4) +
     >                 plai(5) + plai(6) +
     >                 plai(7) + plai(8) 
c
        totlail  =  plai(9)  + plai(10) + plai(11) + plai(12) 
c
        totbiou  = biomass(1) + biomass(2) + biomass(3) + biomass(4) +
     >                biomass(5) + biomass(6) + biomass(7) + biomass(8) 
c
        totbiol  = biomass(9) + biomass(10) + biomass(11) + biomass(12) 
c
c initial single-sided sai for upper and lower canopies
c
        sai(1)  =  0.050 * totlail
        sai(2)  =  0.250 * totlaiu
c
c fractional cover
c
c       fu = wood / woodnorm
c
        fu = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
c
        fl = totlail / 1.0
c
        fu = max (0.25, min (0.975, fu))
        fl = max (0.25, min (0.975, fl))
c
c initial lai for canopy physics
c
        lai(1) = totlail / fl
        lai(2) = totlaiu / fu
c
c specify canopy height parameters
c calculated as a function of only the vegetative fraction
c of each grid cell
c
        zbot(1) =  0.05
c       ztop(1) =  max (0.25, totlail * 0.25)
        ztop(1) =  max (0.25, lai(1) * 0.25)
c
        zbot(2) =  ztop(1) + 1.0 
c       ztop(2) =  max (zbot(2) + 1.00, 2.50 * totbiou * 0.75)
        ztop(2) =  max (zbot(2) + 1.00, 2.50 * totbiou / fu * 0.75)
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
c
      rhoveg(2,1) = 0.60     ! nir leaf reflectance, lower story
      rhoveg(2,2) = 0.40     ! nir leaf reflectance, upper story
c
      tauveg(1,1) = 0.07     ! vis leaf transmittance, lower story
      tauveg(1,2) = 0.05     ! vis leaf transmittance, upper story
c
      tauveg(2,1) = 0.25     ! nir leaf transmittance, lower story
      tauveg(2,2) = 0.20     ! nir leaf transmittance, upper story
c
      chiflz = -0.5          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
      chifuz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
c
      oriev(1) = max (-chiflz, 0.)
      oriev(2) = max (-chifuz, 0.)
c
      orieh(1) = max ( chiflz, 0.)
      orieh(2) = max ( chifuz, 0.)
c
      dleaf(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
      dstem(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
c
      dleaf(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
      dstem(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
c
      chu = ch2o *  2.0      ! heat capacity of upper leaves
      chl = ch2o *  2.0      ! heat capacity of lower leaves
      chs = ch2o * 50.0      ! heat capacity of stems
c
      alaimu = 8.0           ! normalization constant for upper canopy aerodynamics
      alaiml = 8.0           ! normalization constant for lower canopy aerodynamics
c
      cleaf  = 0.01          ! constant in leaf-air aero transfer parameterization
      cgrass = 0.01          ! constant in leaf-air aero transfer parameterization
      cstem  = 0.01          ! constant in leaf-air aero transfer parameterization
c
      wliqumax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
      wliqsmax = 0.40        ! intercepted water capacity (mm h2o per unit leaf area)
      wliqlmax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
c
      wsnoumax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
      wsnosmax = 4.00        ! intercepted snow capacity (mm h2o per unit leaf area)
      wsnolmax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
c
      tdripu =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
      tdrips =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
      tdripl =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
c
      tblowu = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
      tblows = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
      tblowl = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
c
c ************************************************************************
c define rooting profiles
c ************************************************************************
c
c define rooting profiles based upon data published in:
c
c Jackson et al., 1996:  A global analysis of root distributions
c for terrestrial biomes, Oecologia, 108, 389-411.
c
c and
c
c Jackson et al., 1997:  A global budget for fine root biomass, 
c surface area, and nutrient contents, Proceedings of the National
c Academy of Sciences, 94, 7362-7366.
c
c rooting profiles are defined by the "beta" parameter
c
c beta1 is assigned to the lower vegetation layer (grasses and shrubs)
c beta2 is assigned to the upper vegetation layer (trees)
c
c according to Jackson et al. (1996, 1997), the values of beta
c typically fall in the following range
c
c note that the 1997 paper specifically discusses the distribution
c of *fine roots* (instead of total root biomass), which may be more
c important for water and nutrient uptake
c
c --------------                 ------------   ------------
c forest systems                 beta2 (1996)   beta2 (1997)
c --------------                 ------------   ------------
c tropical evergreen forest:        0.962          0.972
c tropical deciduous forest:        0.961          0.982
c temperate conifer forest:         0.976          0.980
c temperate broadleaf forest:       0.966          0.967
c all tropical/temperate forest:    0.970  
c boreal forest:                    0.943          0.943
c all trees:                                       0.976
c
c -------------------------      ------------   ------------
c grassland / shrub systems      beta1 (1996)   beta1 (1997)
c -------------------------      ------------   ------------
c tropical grassland / savanna:     0.972          0.972
c temperate grassland:              0.943          0.943
c all grasses:                      0.952          0.952
c schlerophyllous shrubs:           0.964          0.950
c all shrubs:                       0.978          0.975
c crops:                            0.961
c desert:                           0.975          0.970
c tundra:                           0.914
c
c --------------                 ------------
c all ecosystems                 beta  (1996)
c --------------                 ------------
c all ecosystems:                   0.966
c
c for global simulations, we typically assign the following
c values to the beta parameters
c
c beta1 = 0.950, which is typical for tropical/temperate grasslands
c beta2 = 0.970, which is typical for tropical/temperate forests
c
c however, these values could be (and should be) further refined
c when using the model for specific regions
c 
      beta1 = 0.950  ! for lower layer herbaceous plants
      beta2 = 0.975  ! for upper layer trees
c
c calculate total depth in centimeters
c
      totdepth = 0.0
c
      do 300 k = 1, nsoilay
        totdepth = totdepth + hsoi(k) * 100.0
 300  continue
c
c normalization factors
c
      frootnorm1 = 1. - beta1 ** totdepth
      frootnorm2 = 1. - beta2 ** totdepth
c
c calculate rooting profiles
c
      do 400 k = 1, nsoilay
c
        if (k.eq.1) then
c
          depth(k) = hsoi(k) * 100.0
c
          froot(k,1) = 1. - beta1 ** depth(k)
          froot(k,2) = 1. - beta2 ** depth(k)
c
        else
c
          depth(k) = depth(k-1) + hsoi(k) * 100.0
c
          froot(k,1) = (1. - beta1 ** depth(k)) - 
     >                 (1. - beta1 ** depth(k-1)) 
c
          froot(k,2) = (1. - beta2 ** depth(k)) - 
     >                 (1. - beta2 ** depth(k-1)) 
c
        endif
c
        froot(k,1) = froot(k,1) / frootnorm1
        froot(k,2) = froot(k,2) / frootnorm2
c
 400  continue
c
c return to main program
c
      return
      end