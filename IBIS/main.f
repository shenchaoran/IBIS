$FIXEDFORMLINESIZE:132
c #    #    ##       #    #    #
c ##  ##   #  #      #    ##   #
c # ## #  #    #     #    # #  #
c #    #  ######     #    #  # #
c #    #  #    #     #    #   ##
c #    #  #    #     #    #    #
c
c ---------------------------------------------------------------
	program main
c ---------------------------------------------------------------
        include 'implicit.h'
        include 'compar.h'
        include 'comatm.h'
        include 'comdiag.h'
        include 'comsum.h'
        include 'comveg.h'
        include 'compft.h'
        include 'comsoi.h'
        include 'comsat.h'
	    include 'soilbgc.h'

        integer argc, argcI, filePathStart
        character(len=10):: fileTag
        character(len=100):: tempStr
        character(len=100), dimension(3):: argvs
377     format("-----this is an error identification-----",/, a)
920     format("-----Progress:",I3,"%-----")

	    integer, parameter:: sitesum = 5595 !252	!格点数目总和
	    integer, parameter:: runsum = 100


        integer isite, iyear, imonth, j
	    integer ihr, imin, isec, i100th
	    integer yearnum, daysum1, irun,	daysum10
	    real year, jd, id
	    real rgd, gppv(3600), resv(3600), neev(3600), lev(3600), hv(3600) 
	    real  gppdobs(3600),resdobs(3600),needobs(3600),ledobs(3600),hdobs(3600),laiobs(3600)
	    real gppobs, reobs, neeobs, leobs, tsfactor, wsfactor

        integer istep,           ! timestep counter (per day)
     >        iday,              ! daily loop counter
     >        iy1,               ! first year for year loop
     >        iy2,iy20,          ! last year for year loop
     >        irestart,          ! 0: normal mode 1: restart mode
     >        isimveg,           ! 0: static veg 1: dynam veg initialized w/ fixed
     >                           ! 2: dynam veg initialized w/ cold start
     >        isimfire,          ! 0: fixed fire  1: dynamic fire
     >        jday,              ! julian day of the simulation
     >        niter,             ! total number of time iterations per day
     >        plen,              ! length of precipitation event in timesteps (see plens)
     >        plenmax,           ! upper bound on plen
     >        plenmin,           ! lower bound on plen
     >        seed,              ! value used to initialize the random # generator
     >        soilcspin          ! 0: no spinup procedure for soil c  1: acceleration procedure used
        real    plens,           ! length of precipitation event in seconds
     >        startp,            ! time to start precipitation event (seconds since midnight)
     >        endp,              ! time to end precipitation event (seconds since midnight)
     >        time               ! time of day since midnight (in seconds)

cc define middle variables
        real rhd, dran, lat, sand0, clay0
	    real soilc, soiln, lon, laid, sun, m, lod
	    real tav(19360), tmaxv(19360),tminv(19360), rhv(19360),laiv(3600)
	    real cloudv(19360), precv(19360), windv(19360), psurfv(3600), sunv(3600),rnwjv(19360)
	    integer row, col, pft, ifpt, idayy, dayy, k, dayy0, jday0,istep0,laimaxjd,flag0
	    integer	cadayid,castep


c/****************allocation********************
        real leafratio, stemratio, rootratio, nppleaf ,sumrootstem
	    real precsum0, cldsum0,etsum0,tasum0,cldmean,tamean,taidx
	    real wateridx,lgt,NitroTs,NitroWs,nitrogen0,water0,minwn
	    real stemnppsum,rootnppsum,leafnppsum ,nppsum
	    real rhsum0,wssum0,rhmean,wsmean,tameanK
        real pwdmid,pwd0,vpd0,psy0,ps0,es0,s0,epr0,rn0,epa0,ep0,pet0
	    real m0,rnwsum0,rnwmean,lod0,rnw0,laisum

c/****************allocation*********************

        real ran2                               ! Function : random number generator
        data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
	    character (len=80)::metFile, outfile, siteFile
	    character (len=6) :: site
        irestart = 0                            ! irestart   0: not a restart run  1: restart run
        soilcspin = 1                           ! soilcspin  0: no soil spinup, 1: acceleration procedure used 
        isimveg = 1                             ! isimveg    0: static veg, 1: dynamic veg, 2: dynamic veg-cold start
        isimfire = 0                            ! isimfire   0: fixed fire, 1: dynam fire
        dtime = 3600.                           ! dtime      time step in secondscc      read (lun,*) idiag
        co2conc = 0.000350	                    ! co2 concentration (mol/mol)
        o2conc  = 0.209000	                    ! o2 concentration (mol/mol)
        niter = int (86400.0 / dtime)           !24 hours

! region argv scr
        argc = iargc()
        if(argc /= 3) then
            print 377, 'invalid argc number, you must input 3 argc'
            stop
        end if
      
        do argcI=1, argc
            tempStr = argvs(argcI)
            call getarg(argcI, tempStr)
            filePathStart = index(tempStr, '=')
            if(filePathStart == 0) then
                print 377, 'invalid file path argvs, the file path must as a prefix of "--tag="!'
                stop
            end if
          
            fileTag = trim(tempStr(:filePathStart))
            if(fileTag == '--input=' .or. fileTag == '-i=') then
                metFile = tempStr(filePathStart+1:)
            elseif(fileTag == '--output=' .or. fileTag == '-o=') then
                outfile = tempStr(filePathStart+1:)
            elseif(fileTag == '--site=' .or. fileTag == '-s=') then
                siteFile = tempStr(filePathStart+1:)
            end if
            print *,tempStr
        end do
      
        if(len_trim(siteFile) == 0 .or. len_trim(metFile) == 0 .or. len_trim(outfile) == 0 ) then
            print 377, 'invalid input file!'
            stop
        end if
! endregion

cc      open(100,file="C:\IBIS\Data\test.txt",action="write")
cc      open(120,file="D:\South_Drought\input\Photo_test2.txt",action="write")
c
c ---------------------------------------------------------------------
c also take care of calculation of texfact, which is a leaching
c parameter based on the average sand fraction of the top 1 m of
c soil
c ---------------------------------------------------------------------

	    !do isite= 1, sitesum
		
	    open(21, file=siteFile,action="read")           
        read(21, *)  lat
        read(21, *)  sand0
        read(21, *)  clay0
c read(21,*)soilc   
c read(21,*)soiln   
        read(21, *)  iy1
        read(21, *)  yearnum
        read(21, *)  pft
        read(21, *)  daysum1
        daysum10 = daysum1
cc	    write(100, "(5f20.2)")   lat,sand0,clay0,soilc,soiln

        call constvars  !constant parameters in constvars.f

	    iy2 = iy1 + yearnum-1
	    iy20=iy2

		sand0 = max(min(95.0, sand0), 5.0)
		clay0 = max(min(95.0, clay0), 5.0)

	    if((pft.le.0).or.(pft.gt.12)) then 
	        pft = 12
	    end if	    

	    open(22, file=metFile,action="read")
	    open(23, file=outfile,action="write")

cc	    write(23,*)"run year day gppsim nppsim co2mic neesim lai gppobs reobs neeobs leobs"

	    !read(22,*)	! read the first line

        do idayy = 1, daysum1
	        read(22,*) site, lat, year, jd, id, tav(idayy), tmaxv(idayy), tminv(idayy), 
     >		         rhv(idayy), precv(idayy), windv(idayy), cloudv(idayy)
	    end do 

cc 必须的要素：tav, tmaxv, tminv, rhv, precv, windv, cloudv

	    clitlm =        0.0
        clitls =        0.0
        clitll =        0.0
        clitrm =        0.0
        clitrs =        0.0
        clitrl =        0.0
        clitwm =        0.0
        clitws =        0.0
        clitwl =        0.0
        totcmic =       0.0
        csoislop =      0.0
        csoislon =      0.0
        csoipas =       0.0
        falll =         0.0
        fallr =         0.0
        fallw =         0.0

	    do ifpt = 1, 12		        !initialize biomass parameters
	        exist(ifpt) = 0
		    cbiol(ifpt) = 0.0       !carbon in leaf biomass pool (kg_C m-2)
		    cbiow(ifpt) = 0.0
		    cbior(ifpt) = 0.0
	    end do

	    exist(pft) = 1	   ! fixed pft

cc      write(100,*) "main1"

	    call initial (isimveg,irestart,sand0,clay0)

cc	    write(100,*) "initial"
c
c initialize random number generator seed
c
        seed = -1
        dran = ran2 (seed)

	    do irun = 1, runsum	        !spin up
c
c initialize this year's values of gdd0, gdd5, tc, tw
c
            dayy = 0
            dayy0 =0

            do iyear = iy1, iy2                 ! start of yearly loop
                tcthis = 100
                twthis = -100
                gdd0this = 0.0
                gdd5this = 0.0

                do j = 1, npft
                    aygpp(j) = 0
                    aynpp(j) = 0
                end do

                aygpptot =          0.0
                aynpptot =          0.0
                ayneetot =          0.0
                ayco2mic =          0.0
                aycsoi =            0.0
                aylai1 =            0.0
                aylai2 =            0.0
                aylail =            0.0
                aylaiu =            0.0
                ayaet =             0.0
                aytrans =           0.0
                aytrunoff =         0.0
                aysrunoff =         0.0
                aydrainage =        0.0
                aywsoi =            0.0
                aywisoi =           0.0
                aysuvap =           0.0
                ayinvap	=           0.0
                ayclitlm =          0.0    
                ayclitls =          0.0   	  
                ayclitll =          0.0
                ayclitrm =          0.0	
                ayclitrs =          0.0	
                ayclitrl =          0.0
                ayclitwm =          0.0
                ayclitws =          0.0	
                ayclitwl =          0.0
                ayfalll =           0.0
                ayfallr =           0.0
                ayfallw =           0.0
                aycsoipas =         0.0	
                aycsoislop =        0.0
                aycsoislon =        0.0
                aycmic =            0.0 

	            write(*,920) irun		            !检查运行步骤

                do k = 1, nsoilay
                    wsoi(k) = swilt(k) + (sfield(k) - swilt(k))/2
                end do

                ndaypm(2) = 28
                ndaypy = 365
                if (mod(iyear,4).eq.0) then
                    if (mod(iyear,100).ne.0) then
                    ndaypm(2) = 29
                    ndaypy = 366
                    else if (mod(iyear/100,4).eq.0) then
                    ndaypm(2) = 29
                    ndaypy = 366
                    end if
                end if

                jday =                  0	                    ! reset julian date
                jday0 =                 0           
                laimaxjd =              200         
                laisum =                0.51            

                do imonth = 1, 12                               ! start of monthly loop
                    do iday = 1, ndaypm(imonth)                 ! start of daily loop
                        dayy = dayy + 1
                        if(dayy.gt.daysum1) then
                            dayy = 1
                        end if

                        td = tav(dayy) + 273.16
                        tmax = tmaxv(dayy) + 273.16
                        tmin = tminv(dayy) + 273.16
                        rhd = max(0.01,min(rhv(dayy),1.0))
                        precip = precv(dayy)
                        ud = windv(dayy)
                        cloud = cloudv(dayy) 
                        psurf = 102350
                    
                        qd = rhd * qsat(esat(td), psurf)
                        jday = jday + 1
c		 
c calculated temperature extremes -- for vegetation limits (deg c)
c for this purpose, use the 10-day running mean temperature
c
                        tcthis = min (tcthis, (a10td - 273.16))
                        twthis = max (twthis, (a10td - 273.16))
c
c update this year's growing degree days
c
                        gdd0this = gdd0this + max(0., (td - 273.16))
                        gdd5this = gdd5this + max(0., (td - 278.16))
c
c determine the daily vegetation cover characteristics
c
                        call pheno(iday,laisum,jday)

cc                      write(100,*) iyear, jday, laid

                        tsfactor = min(exp(300 * ((1./(287 - 227.13)) - (1./(adtsoic-227.13)))), 4.0)
	     		        wsfactor = (1-exp(-6.0*adwsoic))/(1-exp(-6.0))
                        ! adtsoic:  daily average soil temperature (c) using profile weighting
                        ! adwsoic: daily average soil moisture using root profile weighting (fraction)

                        decompl = tsfactor * wsfactor	            ! litter decomposition factor
                        decomps = tsfactor * wsfactor	            ! soil organic matter decomposition factor

                        call soilbgc (irun, iyear,iy1,imonth,iday,sand0, clay0)
cc	                    write(100,*) "soilbgc"
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
                        plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
                        plenmax = max (int (24.0 * 3600. / dtime), plenmin)
                        plen    = min (plenmax, int (plenmin + ran2(seed) * (plenmax-plenmin+1) ))
                        plens   = dtime * plen
                        startp  = dtime * min (niter-plen, int(ran2(seed)*(niter-plen+1)))
                        endp    = startp + plens
                        
                        do istep = 1, niter             ! start of hourly loop, niter=24
c
c calculate the time of day since midnight (in seconds)
c
                            time = (istep - 1) * dtime  !dtime = 3600. ,  time step in secondscc 
c
c determine climatic conditions for the given timestep
c                  
                            call diurnal (iyear,lat, time, jday, plens, startp, endp, seed)
cc		                    write(100,*)"diurnal"

                            call lsxmain(iyear, imonth, iday, pft)
cc		                    write(100,*)"lsxmain"
c
c accumulate some variables every timestep
c
                            call sumnow(iyear, imonth, iday, istep)
cc		                    write(100,*)"sumnow"

                            call sumday (iyear, imonth, iday, istep)
cc   	                    write(100,*)"sumday"
                        end do          ! end of the hourly loop
c
c write out daily output

cc****************************Carbon allocation****************************
      	                adnpptot = max(adnpptot, 0.0) 
cc	                    write(100,*) "adnpptot", adnpptot
c/**********************************************************/
  
                        call sumyear  (imonth, iday)
cc		                write(100,*) "sumyear"
                    end do              ! end of the daily loop
c
cc                  write(100,*)iyear,imonth,iday,"sumyear over"
                end do                  ! end of the monthly loop

c
c perform vegetation dynamics
c/**********************************************************/

cc	            write(100,"(5f20.2)") aleaf(pft), awood(pft), aroot(pft), nppsum, leafnppsum
cc              write(100,*)aleaf(pft), awood(pft), aroot(pft), nppsum
                leafnppsum = 0
                stemnppsum = 0
                rootnppsum =0
                nppsum = 0
                lai(1) = 0
                lai(2) = 0
c/**********************************************************/
                if (isimveg.ne.0) call dynaveg (iyear, isimfire)	 
                if(irun.eq.runsum) then           
		            write(23,"(I6,5f13.8)") iyear,falll,fallw,aylail,aylaiu,ayco2mic
  	            end if									  
            end do              ! end of year loop
        end do                  ! end of spin-up loop	
	  
	    close(21)
	    close(22)
        close(23)
        !end do                 ! end of site loop
c
c end of the simulation
c
cc      return
        !pause
        end
c
c
c ---------------------------------------------------------------
      subroutine lsxmain(iyear, imonth, iday, pft)
c ---------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
c
      include 'compar.h'
      include 'com1d.h'
      include 'comatm.h'
      include 'comsoi.h'
c
c Local variables
c
      integer ib,       ! waveband number (1= visible, 2= near-IR)
     >        pft       ! loop indice

	integer iyear, imonth, iday

c
c set physical soil quantities
c
      call setsoi
cc      write(100,*)"start lsx"
c
c calculate areal fractions wetted by intercepted h2o
c
      
      call fwetcal
cc	write(100,*)"fwetcal over"
c
c set up for solar calculations
c
      call solset
cc	write(100,*)"solset over"
c
c solar calculations for each waveband
c
      do ib = 1, nband
c
c solsur sets surface albedos for soil and snow
c solalb performs the albedo calculations
c solarf uses the unit-incident-flux results from solalb
c to obtain absorbed fluxes sol[u,s,l,g,i] and 
c incident pars sunp[u,l]
c
        call solsur (ib)
cc	  write(100,*)"solsur over"

        call solalb (ib)
cc	  write(100,*)"solalb over"

        call solarf (ib)
cc	  write(100,*)"solarf over"
c
      end do
c
c calculate ir fluxes
c
      call irrad
cc	write(100,*)"irrad over"
c
c step intercepted h2o
c
      call cascade
cc	write(100,*)"cascade over"
c
c re-calculate wetted fractions, changed by cascade
c
      call fwetcal
cc	write(100,*)"fwetcal over"
c
c step vegetation canopy temperatures implicitly
c and calculate sensible heat and moisture fluxes
c
cc       write(100,*)"canopy over", pft

      call canopy(iyear, imonth, iday, pft)
c
c step intercepted h2o due to evaporation
c
      call cascad2
cc	write(100,*)"cascad2 over"
c
c arbitrarily set veg temps & intercepted h2o for no-veg locations
c
      call noveg
cc	write(100,*)"noveg over"
c
c set net surface heat fluxes for soil and snow models
c
cc      do 110 i = 1, npoi
c
        heatg = solg + firg - fseng - hvasug*fvapg
c
        heati = soli + firi - fseni - hvasui*fvapi
c
cc 110  continue
c
c step snow model

c
      call snow
cc	write(100,*)"snow over"
c
c step soil model
c
      call soilctl
cc	write(100,*)"soilctl over"
c
c return to main program
c							   
      return
      end
c 
c ---------------------------------------------------------------------
      subroutine coldstart
c ---------------------------------------------------------------------
c  
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsno.h'
c
c initialize some model variables for cold start conditions
c
      fi = 0.0
c
      call const (hsno, nsnolay, 0.0)
      call const (tsno, nsnolay, 273.16)
c
      call const (tsoi, nsoilay, 278.16)
      call const (wsoi, nsoilay, 0.50)
      call const (wisoi, nsoilay, 0.00)
c
c return to main program
c
      return
      end
