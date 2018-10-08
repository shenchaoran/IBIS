c
c ------
c comhyd
c ------
c
      real 
     >  ginvap,  ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
     >  gsuvap,  ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
     >  gtrans,  ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
     >  gtransu, ! transpiration from upper canopy (kg_h2o m-2 s-1)
     >  gtransl, ! transpiration from lower canopy (kg_h2o m-2 s-1)
     >  grunof,  ! surface runoff rate (kg_h2o m-2 s-1)
     >  gdrain,  ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
     >  gadjust, ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
     >  wtot     ! total amount of water stored in snow, soil, puddels, and on vegetation (kg_h2o)
c
      common /comhyd1/ ginvap, gsuvap, gtrans, gtransu, gtransl, grunof, gdrain,
     >     gadjust, wtot


