$FIXEDFORMLINESIZE:132
c #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
c #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       ##  ##
c #####      #    #    #  #       #####   #    #  #       ######  #####   # ## #
c #    #     #    #    #  #  ###  #       #    #  #       #    #  #       #    #
c #    #     #    #    #  #    #  #       #    #  #    #  #    #  #       #    #
c #####      #     ####    ####   ######   ####    ####   #    #  ######  #    #
c
c
c --------------------------------------------------------------------------
      subroutine soilbgc (irun, iyear, iyear0, imonth, iday, sand0, clay0)
c --------------------------------------------------------------------------
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
      include 'comatm.h'
	include 'soilbgc.h'
c
c Arguments (input)
c
      integer   iday,           ! day in month
     >          iyear,          ! current year
     >          iyear0,         ! initial year
     >          imonth         ! current month
cc     >          nspinsoil,      ! year when soil carbon spinup stops
cc     >          spin,           ! # of times soilbgc has been called in the current day
cc     >          spinmax         ! total # of times soilbgc is called per day (spinup)
c 
c local variables
c
      integer irun            ! loop indice
c
      real totts,          ! 1/ndaypy
     >     fracll,         ! lignin fraction of leaves  木质素
     >     fracls,         ! structural fraction of leaves   结构
     >     fraclm,         ! metabolic fraction of leaves   新陈代谢
     >     fracrl,         ! lignin fraction of roots
     >     fracrs,         ! structural fraction of roots 
     >     fracrm,         ! metabolic fraction of roots
     >     fracwl,         ! lignin fraction of wood
     >     fracws,         ! structural fraction of wood
     >     fracwm          ! metabolic fraction of wood 

      real      outclm,   ! c leaving leaf metabolic pool 
     >          outcls,   ! c leaving leaf structural pool 
     >          outcll,   ! c leaving leaf lignin pool
     >          outcrm,   ! c leaving root metabolic pool 
     >          outcrs,   ! c leaving root structural pool 
     >          outcrl,   ! c leaving root lignin pool
     >          outcwm,   ! c leaving woody metabolic carbon pool
     >          outcws,   ! c leaving woody structural carbon pool
     >          outcwl,   ! c leaving woody lignin carbon pool
     >          outcsb,   ! flow of passive c to biomass
     >          outcps,   ! flow of protected om to passive pool 
     >          outcns,   ! flow of non-protected om to passive pool
     >          outcnb,   ! flow of non-protected om to biomass 
     >          outcpb,   ! flow of protected om to biomass
     >          outcbp,   ! c leaving protected biomass pool  
     >          outcbn,   ! c leaving non-protected biomass pool
     >          totc      ! total c in soil
c
      real      dbdt,     ! change of c in biomass pools with time 
     >          dcndt,    ! change of c in non-protected om with time
     >          dcpdt,    ! change of c in protected om with time
     >          dcsdt,    ! change of c in passive om with time
     >          totmin,   ! total nitrogen mineralization
     >          totimm,   ! total nitrogen immobilization 
     >          netmin,   ! net nitrogen mineralization
     >          nbiors,
     >          nbiols,
     >          nbiows,
     >          nbiowm,
     >          nbiolm,
     >          nbiorm,
     >          nbioslon,
     >          nbioslop,
     >          nbiopas,
     >          nminrs,
     >          nminls,
     >          nminws,
     >          nminwm,
     >          nminlm,
     >          nminrm,
     >          nminslon,
     >          nminslop,
     >          nminpas,
     >          nrelps,
     >          nrelns,
     >          nrelbn,
     >          nrelbp,
     >          nrelll,
     >          nrelrl,
     >          nrelwl,
     >          totnrel,
     >          ymintot,
     >          yminmic
c
c nitrogen in litter and soil pools
c
      real      nlitlm,
     >          nlitls,
     >          nlitll,
     >          nlitrm,
     >          nlitrs, 
     >          nlitrl,
     >          nlitwm,
     >          nlitws,
     >          nlitwl,
     >          nsoislop,
     >          nsoipas,
     >          nsoislon
cc     >          budgetn
c
c variables controlling constraints on microbial biomass 
c
      real      cmicn,
     >          cmicp,
     >          cmicmx
c
c variables controlling leaching, calculating co2 respiration and n deposition
c
      real      cleach,
     >          totcbegin,
     >          totcend,
     >          totcin,
     >          fixsoin,
     >          deposn
c
      real      fleach,
     >          h20
c
c decay constants for c pools
c
c 
c variables added to do daily time series of some values
c
cc      integer   gridpt,
      integer          kk
c
c variables dealing with soil texture and algorithms
c
      integer   msand,
     >          mclay
c
      real      fsand,
     >          fclay,
     >          cfrac,
     >          sand0,
     >          clay0,
     >          texfact,
     >          fbpom,
cc     >          fbsom,
     >          rdepth
cc     >          effac,
cc     >          lig_frac
c
cc       gridpt = npoi		   ! total number of gridpoints used
c
c total timesteps (daily) used to divide litterfall into daily fractions 
c
       totts=1./float(ndaypy)
c
c -------------------------------------------------------------------------------------
c specific maximum decay rate or growth constants; rates are per day
c constants are taken from Parton et al., 1987 and Verberne et al., 1990
c and special issue of Geoderma (comparison of 9 organic matter models) in Dec. 1997
c
c leaching parameterization was changed to agree with field data,
c this caused a changing of the below constants.  
c
c approximate factors for Verberne et al. model where efficiencies are 100%
c for some of the transformations: one problem was that their rate constants were
c based on 25C, and our modifying functions are based on 15 C...thus the rate constants
c are somewhat smaller compared to the Verberne et al. (1990) model parameters
c rates are based on a daily decomposition timestep (per day)
c -------------------------------------------------------------------------------------
c
c leaf litter 
c
      write(100,*)"bgc 1", rconst,cnr(7)

       fracll = fmax * (cnleaf**2)/(rconst + cnleaf**2)
       fracls = (1./cnleaf - fracll/cnr(5) - (1.-fracll)/cnr(7))/
     >          (1./cnr(6) - 1./cnr(7))
       fraclm = 1.0 - fracll - fracls

      write(100,*)"bgc 1.2"
c
c root litter
c
       fracrl = fmax * (cnroot**2)/(rconst + cnroot**2)
       fracrs = (1./cnroot - fracrl/cnr(5) - (1.-fracrl)/cnr(7))/
     >          (1./cnr(6) - 1./cnr(7))
       fracrm = 1.0 - fracrl - fracrs

      write(100,*)"bgc 1.3"
c
c wood litter
c
       fracwl = fmax * (cnwood**2)/(rconst + cnwood**2)
       fracws = (1./cnwood - fracwl/cnr(5) - (1.-fracwl)/cnr(7))/
     >          (1./cnr(8) - 1./cnr(7))
       fracwm = 1.0 - fracwl - fracws

c
      write(100,*)"bgc 1.5"
c
        rdepth   = 1./(hsoi(1) + hsoi(2) + hsoi(3) + hsoi(4))
        cfrac    = 0.0
        texfact  = 0.0 

        do 90 kk = 1, 4                    ! top 1 m of soil -- 4 layers
          msand    = nint(sand0) 
          mclay    = nint(clay0) 
          fclay    = 0.01 * mclay
          fsand    = 0.01 * msand 
          cfrac    = cfrac   + fclay * hsoi(kk)
          texfact  = texfact + fsand * hsoi(kk)
 90     continue
c
        cfrac   = cfrac   * rdepth
        texfact = texfact * rdepth
c
c ---------------------------------------------------------------------
c fraction of decomposing microbial biomass into protected organic
c matter; taken from the model of Verberne et al., 1990
c this is the proportion of decomposing dead microbial biomass that
c is transferred to a protected pool vs. a non-protected pool
c related to the clay content of the soil. in sandy soils, fbpom = 0.3,
c whereas in clay soils fbpom = 0.7.  created a linear function based
c on clay fraction of soil to adjust according to amount of clay in
c the top 1 m of profile (weighted average according to depth of each
c layer)
c
c
c if cfrac is greater than 0.4, set fbpom = 0.7, if cfrac is less
c than 0.17, set fbpom = 0.30 (sandy soil)
c
         fbpom = 0.50

      write(100,*)"bgc 1.6"

c
c ------------------------------------------------------------------------
c total soil carbon initialized to 0 at beginning of model run
c used in calculation of soil co2 respiration from microbial decomposition 
c ------------------------------------------------------------------------
c
       if (iday .eq. 1 .and. imonth .eq. 1 .and. iyear .eq. iyear0) then
         totcbegin = 0.0
         storedn   = 0.0 
       endif
c
c ------------------------------------------------------------------------
c initialize yearly summation of net mineralization and co2 respiration
c to 0 at beginning of each year; because these quantities are usually 
c reported on a yearly basis, we wish to do the same in the model so we
c can compare easily with the data.
c ------------------------------------------------------------------------
c
       if (iday .eq. 1 .and. imonth .eq. 1) then
         yrleach = 0.0
         cleach  = 0.0
         ynleach = 0.0
         ymintot = 0.0
         yminmic = 0.0
       endif
c       
c determine amount of substrate available to microbial growth
c
c calculate the total amount of litterfall entering soil(C)
c
       totcin =  falll*totts + fallr*totts + fallw*totts
c
c
c calculate the current total amount of carbon at each grid cell
c
       totc = clitlm + clitls + clitrm + clitrs +
     >        clitwm + clitws + csoislop + csoislon +
     >        csoipas + totcmic + clitll + clitrl + clitwl

       write(100,"(14f15.4)")totc, clitlm,clitls,clitrm,clitrs,
     >       clitwm, clitws, csoislop, csoislon,
     >       csoipas, totcmic, clitll, clitrl, clitwl

c
c beginning amount of soil C at each timestep (used for respiration
c calculation)
c
       totcbegin = totc
c
c ------------------------------------------------------------------------
c split current amount of total soil microbes
c maximum amount of biomass is a function of the total soil C
c from Verberne et al., 1990
c ------------------------------------------------------------------------
c
c      totcmic(i) = cmicp(i) + cmicn(i)
       cmicmx = fbsom * totc 
c
c calculate the amount of protected and unprotected biomass
c
       if (totcmic .ge. cmicmx) then
         cmicp = cmicmx
         cmicn = totcmic - cmicmx
       else
         cmicn = 0.0
         cmicp = totcmic
       endif
c
c ---------------------------------------------------------------
c litter pools 
c
c add in the amount of litterfall, and root turnover
c ---------------------------------------------------------------
c
      write(100,*)"bgc 2"

       clitlm = clitlm + (fraclm * falll*totts)  
       clitls = clitls + (fracls * falll*totts)  
       clitll = clitll + (fracll * falll*totts)  
       clitrm = clitrm + (fracrm * fallr*totts)  
       clitrs = clitrs + (fracrs * fallr*totts)  
       clitrl = clitrl + (fracrl * fallr*totts)  
       clitwm = clitwm + (fracwm * fallw*totts)  
       clitws = clitws + (fracws * fallw*totts)  
       clitwl = clitwl + (fracwl * fallw*totts)  
c
c ---------------------------------------------------------------
c calculate microbial growth rates based on available C sources
c to microbes (substrate : litter, C in slow, passive pools)
c the amount of biomass added cannot be larger than the amount of
c available carbon from substrates and other pools at this point.
c ---------------------------------------------------------------
c
c
       outcrs = min(decomps * krs * clitrs,clitrs)
       outcws = min(decompl * kws * clitws,clitws)
       outcls = min(decompl * kls * clitls,clitls)
       outclm = min(decompl * klm * clitlm,clitlm)
       outcrm = min(decomps * krm * clitrm,clitrm)
       outcwm = min(decompl * kwm * clitwm,clitwm)
       outcnb = min(decomps * knb * csoislon,csoislon) 
       outcpb = min(decomps * kpb * csoislop,csoislop)
       outcsb = min(decomps * ksb * csoipas,csoipas)
c
c ---------------------------------------------------------------
c calculate turnover of microbial biomass
c two disctinct pools: one with rapid turnover, and one with slow
c turnover rate
c ---------------------------------------------------------------
c
       outcbp = min(kbp * cmicp,cmicp)
       outcbn = min(kbn * cmicn,cmicn)
c
c ---------------------------------------------------------------------
c recycle microbes back to respective microbial pools based on effac as
c discussed in NCSOIL model from Molina et al., 1983
c ---------------------------------------------------------------------
c
       outcbp = outcbp *  effac
       outcbn = outcbn *  effac
c
c -------------------------------------------------------------------------
c have to adjust inputs into microbial pool for the slow
c and passive carbon amounts that are leaving their respective
c pools at an increased rate during the spinup procedure.
c these values should be decreased by the respective spinup factors
c because the microbial pools will otherwise become larger without
c scientific reason due to the spinup relationships used.
c 3 main pools: outcpb, outcnb, outcsb
c -------------------------------------------------------------------------
c
       dbdt =  outcrs * yrs + outcws * yws +
     >            outcls * yls + outclm * ylm +
     >            outcrm * yrm + outcwm * ywm +
     >            outcnb * ynb + 
     >            outcpb * ypb +
     >            outcsb * ysb - outcbp -
     >            outcbn
c 
c -------------------------------------------------------------------------
c change in non-protected organic matter from growth in microbial
c biomass, lignin input, and stablized organic matter pool
c the flow out of the pool from its decomposition is always less
c the yield--which is factored into the pool it is flowing into
c -------------------------------------------------------------------------
c
       outcll = min(decompl * kll * clitll,clitll)
       outcrl = min(decomps * krl * clitrl,clitrl)
       outcwl = min(decompl * kwl * clitwl,clitwl)
       outcns = min(decomps * kns * csoislon,
     >                 csoislon)

	
cc	  if(iyear.gt.2005) then
           write(110,"(3I6,2f10.2)")iyear,imonth,iday,decompl,decomps
cc	  end if
c
c ------------------------------------------------------------ 
c the lig_frac  factor only applies to lignin content...half goes to
c protected slow OM, and half goes to non protected slow OM
c ------------------------------------------------------------
c
       dcndt =  (lig_frac * (outcll * yll + outcrl * yrl +
     >             outcwl * ywl) +
     >             (1. - fbpom) * (ybn * outcbn +
     >             ybp * outcbp)) - outcnb - outcns

      write(100,*)"bgc 3"

c
c ------------------------------------------------------------
c change in protected organic matter from growth in microbial 
c biomass, lignin input, and stablized organic matter pool
c ------------------------------------------------------------
c
       outcps = min(decomps * kps * csoislop, csoislop)
c
c ------------------------------------------------------------
c the lig_frac factor only applies to lignin content...half goes to
c protected slow OM, and half goes to non protected slow OM
c ------------------------------------------------------------
c
       dcpdt = (lig_frac * (outcll*yll+outcrl*yrl +
     >            outcwl * ywl) +
     >            fbpom * (ybn * outcbn +
     >            ybp * outcbp)) - outcpb - outcps
c
c ----------------------------------------------------------------------
c change in stablized organic matter (passive pool) from growth
c in microbial biomass, and changes in protected and unprotected
c SOM
c
c add a loss of C due to leaching out of the profile, based
c on approximation of CENTURY model below 1 m in depth
c based on water in the profile, and texture of soil
c tuned to known outputs or leaching that has been measured in the field
c at Arlington-WI (Courtesy K. Brye, MS) and applied to the global scale
c on average, this calibration yields about 10-50 Kg C ha-1 yr-1 leaching
c depending on C in soil...will need to be tied to an amount of water
c flowing through the profile based upon precipitation eventually
c ----------------------------------------------------------------------
c
         h20    = 0.30e-03
c
c h20 is a constant relating to the flow of water through the top 1 m of the
c profile 
c use texfact -- the % sand -- or texture factor effect on leaching (see Parton
c et al. (1991) calculated from the average sand content of top 1 m of soil
c in the model
c
        fleach = h20/18.0 * (0.01 + 0.04 * texfact)
c
c --------------------------------------------------------------------
c change in passive organic carbon pool
c ---------------------------------------------------------------------
c
       dcsdt = ((yns * outcns) + (yps * outcps)) -
     >            outcsb -  (fleach * csoipas)
c
       cleach = fleach * csoipas + fleach * csoislop +
     >             fleach * csoislon

       write(100,"(5f10.2)")texfact, fleach, csoipas, csoislop, csoislon
c
       ynleach = ynleach + fleach * csoipas/cnr(2) +
     >              fleach * csoislop/cnr(3) +
     >              fleach * csoislon/cnr(4)
c
c update slow pools of carbon for leaching losses
c
       dcndt = dcndt - fleach * csoislon		
       dcpdt = dcpdt - fleach * csoislop		
c
cc      if (spin .eq. spinmax) then

         yrleach =  cleach + yrleach
c
cc      endif
c
c ---------------------------------------------------------------------
c calculate the amount of net N mineralization or immobilization
c ---------------------------------------------------------------------
c
c uptake of n by growth of microbial biomass
c
c immobilized n used for requirements of microbial growth
c is based on flow of carbon and the difference of C/N ratio of
c the microbes and their efficiency versus the C/N ratio of the
c material that is being decomposed 
c
c
c ------------------------------
c structural root decomposition 
c ------------------------------
c
       if (yrs/cnr(1) .gt. 1./cnr(6)) then
         nbiors = (1./cnr(6) - yrs/cnr(1)) * outcrs
         nminrs = 0.0
c
       else
         nminrs = (1./cnr(6) - yrs/cnr(1)) * outcrs
         nbiors = 0.0
       endif
c
c ------------------------------
c structural leaf decomposition
c ------------------------------
c
       if (yls/cnr(1) .gt. 1./cnr(6)) then
         nbiols = (1./cnr(6) - yls/cnr(1))
     >                * outcls
         nminls = 0.0
c
       else
         nminls = (1./cnr(6) - yls/cnr(1))
     >                * outcls
         nbiols = 0.0
       endif
c
c
c ------------------------------
c structural wood decomposition
c ------------------------------
c	
       if (yws/cnr(1) .gt. 1./cnr(8)) then
         nbiows = (1./cnr(8) - yws/cnr(1))
     >                * outcws
         nminws = 0.0
c
       else
         nminws = (1./cnr(8) - yws/cnr(1))
     >                * outcws
         nbiows = 0.0
       endif
c
c ------------------------------
c metabolic wood decomposition
c ------------------------------
c
       if (ywm/cnr(1) .gt. 1./cnr(8)) then
         nbiowm = (1./cnr(8) - ywm/cnr(1))
     >                * outcwm
         nminwm = 0.0
c
       else
         nminwm = (1./cnr(8) - ywm/cnr(1))
     >                * outcwm
         nbiowm = 0.0
       endif
c
c ------------------------------
c metabolic leaf decomposition
c ------------------------------
c
       if (ylm/cnr(1) .gt. 1./cnr(7)) then
         nbiolm = (1./cnr(7) - ylm/cnr(1))
     >                * outclm
         nminlm = 0.0
c
       else
         nminlm = (1./cnr(7) - ylm/cnr(1))
     >                * outclm
         nbiolm = 0.0
       endif
c
c
c ------------------------------
c metabolic root decomposition
c ------------------------------
c
       if (yrm/cnr(1) .gt. 1./cnr(7)) then
         nbiorm = (1./cnr(7) - yrm/cnr(1))
     >                * outcrm
         nminrm = 0.0
c
       else
         nminrm = (1./cnr(7) - yrm/cnr(1))
     >                * outcrm
         nbiorm = 0.0
       endif
c
      write(100,*)"bgc 4"
c
c ----------------------------------------------
c non-protected organic matter decomposition
c ----------------------------------------------
c
       if (ynb/cnr(1) .gt. 1./cnr(4)) then
         nbioslon = (1./cnr(4) - ynb/cnr(1))
     >                  * outcnb
         nminslon = 0.0
c
       else
         nminslon = (1./cnr(4) - ynb/cnr(1))
     >                  * outcnb
         nbioslon = 0.0
       endif
c
c
c ----------------------------------------------
c protected organic matter decomposition
c ----------------------------------------------
c
       if (ypb/cnr(1) .gt. 1./cnr(3)) then
         nbioslop = (1./cnr(3) - ypb/cnr(1))
     >                  * outcpb
         nminslop = 0.0
c
       else
         nminslop = (1./cnr(3) - ypb/cnr(1))
     >                  * outcpb
         nbioslop = 0.0
       endif
c
c
c ----------------------------------------------
c stablized organic matter decomposition
c ----------------------------------------------
c
       if (ysb/cnr(1) .gt. 1./cnr(2)) then
         nbiopas = (1./cnr(2) - ysb/cnr(1))
     >                 * outcsb
         nminpas = 0.0
c
       else
         nminpas = (1./cnr(2) - ysb/cnr(1))
     >                 * outcsb
         nbiopas = 0.0
       endif
c
c ----------------------------------------------
c total immobilized N used for biomass growth
c ----------------------------------------------
c
       totimm = nbiors + nbiols + nbiows + nbiowm
     >           + nbiolm + nbiorm + nbioslon + nbioslop
     >           + nbiopas
c
c -----------------------------------------------------------------------------
c gross amount of N mineralized by decomposition of C by microbial biomass
c assume that N is attached to the flow of C by the C/N ratio of the substrate
c also assume that the amount of N attached to CO2 that is respired is also
c mineralized (i.e. the amount of N mineralized is related to the total outflow
c of carbon, and not the efficiency or yield)..see Parton et al., 1987
c -----------------------------------------------------------------------------
c
       totmin = nminrs + nminls + nminws + nminwm
     >           + nminlm + nminrm + nminslon + nminslop
     >           + nminpas
c
c -----------------------------------------------------------------------------
c when carbon is transferred from one pool to another, each pool has a distinct
c C:N ratio.  In the case of pools where carbon is moving from the pool to 
c the microbial biomass (used for growth/assimilation), net mineralization
c takes place (N is released) after the requirements of building the biomass
c are met.  In the cases of other transformations of C, N is not conserved
c if it follows from one pool to another which has a different C:N ratio;
c either N is released or is needed to make the transformation and keep N
c conserved in the model. 
c
c other calculations of either N release or immobilization to keep track of
c the budget
c
        nrelps = outcps * (1./cnr(3) - 1./cnr(2))
        nrelns = outcns * (1./cnr(4) - 1./cnr(2))
        nrelbn = (1.-fbpom) * outcbn * (1./cnr(1) - 1./cnr(4)) +
     >              (1.-fbpom) * outcbp * (1./cnr(1) - 1./cnr(4))
        nrelbp = fbpom * outcbp * (1./cnr(1) - 1./cnr(3)) +
     >              fbpom * outcbn * (1./cnr(1) - 1./cnr(3))
        nrelll = lig_frac * outcll * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcll * (1./cnr(5) - 1./cnr(4))
        nrelrl = lig_frac * outcrl * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcrl * (1./cnr(5) - 1./cnr(4))
        nrelwl = lig_frac * outcwl * (1./cnr(5) - 1./cnr(3)) +
     >              lig_frac * outcwl * (1./cnr(5) - 1./cnr(4))
c
        totnrel = nrelps + nrelns + nrelbn +
     >               nrelbp + nrelll + nrelrl + nrelwl
c
c -----------------------------------------------------------------------------
c calculate whether net mineralization or immobilization occurs
c on a grid cell basis -- tnmin is an instantaneous value for each time step
c it is passed along to stats to calculate, daily, monthly and annual totals
c of nitrogen mineralization
c
c this is for mineralization/immobilization that is directly related to 
c microbial processes (oxidation of carbon)
c
c the value of totnrel would need to be added to complete the budget
c of N in the model. Because it can add/subtract a certain amount of N
c from the amount of net mineralization.  However, these transformations
c are not directly related to microbial decomposition, so do we add them
c into the value or not?
c -----------------------------------------------------------------------------
c
           netmin = totmin + totimm + totnrel 
           if (netmin .gt. 0.0) then
c
              tnmin = netmin
c	
           else
c	  
              tnmin = 0.0 
c
           endif
c
c convert value of tnmin of Kg-N/m2/dtime to mole-N/s
c based on N = .014 Kg/mole -- divide by the number of seconds in daily timestep
c
            tnmin = tnmin/(86400. * 0.014)

      write(100,*)"bgc 5"

c
c ---------------------------------------------------
c update soil c pools for transformations of c and n
c ---------------------------------------------------
c
           totcmic  = max(totcmic  + dbdt, 0.0)
           csoislon = max(csoislon + dcndt,0.0)
           csoislop = max(csoislop + dcpdt,0.0)
           csoipas  = max(csoipas  + dcsdt,0.0)
           clitlm   = max(clitlm  - outclm,0.0)
           clitls   = max(clitls  - outcls,0.0)
           clitll   = max(clitll  - outcll,0.0)
           clitrm   = max(clitrm  - outcrm,0.0)
           clitrs   = max(clitrs  - outcrs,0.0)
           clitrl   = max(clitrl  - outcrl,0.0)
           clitwm   = max(clitwm  - outcwm,0.0)
           clitws   = max(clitws  - outcws,0.0)
           clitwl   = max(clitwl  - outcwl,0.0)

c
c -----------------------------------------------------------
c update soil n pools based on c:n ratios of each pool
c this approach is assuming that the c:n ratios are remaining
c constant through the simulation. flow of nitrogen is attached
c to carbon 
c -----------------------------------------------------------
c
           totnmic  = totcmic /cnr(1)
           nsoislon = csoislon/cnr(4)
           nsoislop = csoislop/cnr(3)
           nsoipas  = csoipas /cnr(2)
           nlitlm   = clitlm  /cnr(7)
           nlitls   = clitls  /cnr(6)
           nlitll   = clitll  /cnr(5)
           nlitrm   = clitrm  /cnr(7)
           nlitrs   = clitrs  /cnr(6)
           nlitrl   = clitrl  /cnr(5)
           nlitwm   = clitwm  /cnr(8)
           nlitws   = clitws  /cnr(8)
           nlitwl   = clitwl  /cnr(8)
c
c total above and belowground litter
c
           totlit =  clitlm + clitls + clitll +
     >                  clitrm + clitrs + clitrl +
     >                  clitwm + clitws + clitwl
c
c sum total aboveground litter (leaves and wood)
c
           totalit = clitlm + clitls + clitwm +
     >                  clitll + clitws + clitwl
c
c sum total belowground litter (roots) 
c
           totrlit = clitrm + clitrs + clitrl
c
c determine total soil carbon amounts (densities are to 1 m depth; Kg/m-2)
c	
           totcsoi = csoipas + csoislop + totcmic + csoislon
c
c calculate total amount of litterfall occurring (total for year)
c
           totfall = falll + fallr + fallw

cc	     if(imonth.eq.40) then
	         write(100,"(3I6,8f10.2)")iyear, imonth, iday, totalit, totrlit, 
     >			 totfall, totcsoi, csoipas, csoislop, totcmic, csoislon 
cc	     end if
c
c nitrogen 
c
c total nitrogen in litter pools (above and belowground)
c
          totnlit =  nlitlm + nlitls + nlitrm + nlitrs +
     >                  nlitwm + nlitws + nlitll + nlitrl +
     >                  nlitwl
c
c sum total aboveground litter   (leaves and wood)
c
          totanlit = nlitlm + nlitls + nlitwm + nlitll + nlitws + nlitwl
c
c sum total belowground litter  (roots)
c
          totrnlit = nlitrm + nlitrs + nlitrl
c
c total soil nitrogen to 1 m depth (kg-N/m**2)
c
          totnsoi = nsoislop + nsoislon + nsoipas  + totnmic + totnlit
c
c --------------------------------------------------------------------------
c calculate running sum of yearly net mineralization, and nitrogen in pool
c available to plants for uptake--during spin up period, can only count one
c of the cycles for each timestep--otherwise false additions will result
c values of yearly mineralization are in Kg/m-2
c --------------------------------------------------------------------------
c
cc        if (spin .eq. spinmax) then
c
          storedn  = storedn + tnmin
c
cc        endif
c
c calculate total amount of carbon in soil at end of cycle
c this is used to help calculate the amount of carbon that is respired
c by decomposing microbial biomass
c
        totcend = totlit + totcsoi
c
c --------------------------------------------------------------------------
c the amount of co2resp is yearly value and is dependent on the amount
c of c input each year, the amount in each pool at beginning of the year,
c and the amount left in the pool at the end of the year
c along with the amount of root respiration contributing to the flux from
c calculations performed in stats.f
c --------------------------------------------------------------------------
c
cc        if (spin .eq. spinmax) then 
c
c --------------------------------------------------------------------------
c only count the last cycle in the spin-up for co2soi
c when the iyear is less than the nspinsoil value...otherwise
c an amount of CO2 respired will be about 10 times the actual
c value because this routine is called articially 10 extra times
c each time step to spin up the soil carbon
c
c add n-deposition due to rainfall once each day, and
c the amount of N fixed through N-fixers.  These equations
c are based on the annual precip input (cm) and are from
c the CENTURY model...Parton et al., 1987.
c The base equations are in units of (g) N m-2 so have to
c divide by 1000 to put in units of Kg.
c
c the values in the equation of 0.21 and -0.18 were adjusted to reflect
c average daily inputs when no precipitation was falling - the original
c constants are for the entire year 
c --------------------------------------------------------------------------
c
            deposn    = (0.0005753  + 0.0028 * (precip*0.1))*1.e-3
            fixsoin   = (-0.0004932 + 0.14   * (precip*0.1))*1.e-3
c
c --------------------------------------------------------------------------
c add to the daily total of co2 flux leaving the soil from microbial
c respiration -- instantaneous value for each timestep
c since this subroutine gets called daily...instantaneous fluxes
c the fluxes need to be put on a per second basis, which will be dependent
c on the timestep.  Furthermore, because the biogeochem subroutine does
c not get called each timestep...an approximation for a timestep average
c microbial flux and nmineralization rate will be applied
c --------------------------------------------------------------------------
c
c calculate daily co2 flux due to microbial decomposition
c
          tco2mic = totcbegin + totcin - totcend - cleach

          write(100,"(5f15.6)") tco2mic, totcbegin, totcin, totcend, cleach
c
c convert co2 flux from kg C/day  (seconds in a daily timestep) to mol-C/s
c based on .012 Kg C/mol
c
          tco2mic = tco2mic/(86400. * 0.012)

cc	    if(iyear.eq.2006) then
	       write(120,"(3I6,3f10.6)")iyear,imonth,iday,tco2mic,decomps,decompl
cc	    end if
	   
c
cc        endif
c
cc 100  continue
c
c return to main
c
        return
        end
c
