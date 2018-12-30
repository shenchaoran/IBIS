$FIXEDFORMLINESIZE:132
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine scaler(scaleu, scalel)
c ---------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsno.h'
      include 'comsoi.h'
      include 'comveg.h'
      include 'com1d.h'
      include 'comsum.h'
      include 'compft.h'
c
c local variables
      real zweight 

      real pxaiu, plaiu, pxail, plail
      real extpar, scaleu, scalel
c
c ---------------------------------------------------------------------
c upper canopy scaling
c ---------------------------------------------------------------------
c
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c an(x) is proportional to apar(x)
c
c therefore, an(x) = an(0) * apar(x) / apar(0)
c an(x) = an(0) * (A exp(-k x) + B exp(-h x) + C exp(h x)) / 
c                 (A + B + C)
c
c this equation is further simplified to
c an(x) = an(0) * exp (-extpar * x)
c
c an(0) is calculated for a sunlit leaf at the top of the canopy using
c the full-blown plant physiology model (Farquhar/Ball&Berry, Collatz).
c then the approximate par extinction coefficient (extpar) is calculated
c using parameters obtained from the two-stream radiation calculation.
c
c an,canopy avg.= integral (an(x), from 0 to xai) / lai
c               = an(0) * (1 - exp (-extpar * xai )) / (extpar * lai)
c
c the term '(1 - exp (-extpar * xai )) / lai)' scales photosynthesis from leaf
c to canopy level (canopy average) at day time. A 10-day running mean of this
c scaling parameter (weighted by light) is then used to scale the respiration
c during night time.
c
c once canopy average photosynthesis is calculated, then the canopy average
c stomatal conductance is calculated using the 'big leaf approach',i.e. 
c assuming that the canopy is a big leaf and applying the leaf-level stomatal
c conductance equations to the whole canopy.
c
c calculate the approximate par extinction coefficient:
c
c extpar = (k * A + h * B - h * C) / (A + B + C)
c
        extpar = (termu(6) * scalcoefu(1) +
     >            termu(7) * scalcoefu(2) -
     >            termu(7) * scalcoefu(3)) /
     >            max (scalcoefu(4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))
c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxaiu = extpar * (lai(2) + sai(2))
        plaiu = extpar *  lai(2)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
         zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c

        if (plaiu.gt.0.0) then
c
c day-time conditions, use current scaling coefficient

c ! total photosynthetically active raditaion absorbed(APAR) by top leaves of upper canopy (W m-2)

          if (topparu.gt.10.) then
c
            scaleu = (1. - exp(-pxaiu)) / plaiu
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparamu = zweight * a10scalparamu + 
     >                         (1. - zweight) * scaleu * topparu
c
            a10daylightu  = zweight * a10daylightu + 
     >                         (1. - zweight) * topparu
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else
c
            scaleu = a10scalparamu / a10daylightu
c
          endif
c
c if no lai present
c
        else
c
          scaleu = 0.0
c
        endif

c				  
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c calculate the approximate extinction coefficient
c
        extpar = (terml(6) * scalcoefl(1) + 
     >            terml(7) * scalcoefl(2) -
     >            terml(7) * scalcoefl(3)) /
     >            max (scalcoefl(4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))

c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxail = extpar * (lai(1) + sai(1))
        plail = extpar *  lai(1)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
        zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c
        if (plail.gt.0.0) then
c
c day-time conditions, use current scaling coefficient
c
          if (topparl.gt.10.) then

c
            scalel = (1. - exp(-pxail)) / plail
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparaml = zweight * a10scalparaml + 
     >                         (1. - zweight) * scalel * topparl
c
            a10daylightl  = zweight * a10daylightl + 
     >                         (1. - zweight) * topparl
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else

            scalel = a10scalparaml / a10daylightl
c
          endif
c
c if no lai present
c
        else
c
          scalel = 0.0
c
        endif
c
c perform scaling on all carbon fluxes from upper canopy
c

      return
      end
