$FIXEDFORMLINESIZE:132
c   ####     ##	 #		 #		  ####	  #### 	   ##	  #####   #    ####   #    #
c  #    #   #  # 	 #		 #		 #	  #	 #	  #	  #  #		#     #   #    #  ##   #
c  #       #    #  #		 #		 #	  #	 #		 #    #		#     #   #    #  # #  #
c  #       ######	 #		 #		 #	  #	 #		 ######		#     #   #    #  #  # #
c  #    #  #    #	 #		 #		 #	  #	 #	  #	 #	  #		#     #   #    #  #	  ##
c   ####   #    #	 ######  ######	  ####	  ####	 #	  #		#     #    ####   #    #
c  
c ------------------------------------------------------------------------------------
      subroutine Callocation()
c ------------------------------------------------------------------------------------
c
c common blocks
c
      include 'implicit.h'
c
      include 'compar.h'
      include 'comatm.h'
      include 'comsoi.h'
      include 'comsum.h'
      include 'comveg.h'
c
c local variables
c
      integer iday,jday,laimaxjd,flag0
cc     >  i,          !
cc     >  imonth     !
cc     >  iday        !
c

c Arguments
c
      integer isimfire  ! fire switch
!     >        isim_ac,   ! age-class dynamics switch
!     >        year      ! year of simulation

      real    pfire     ! probability of fire -- should be determined externally.
cc      real        pdist      ! probability of other disturbance types....
       
      PARAMETER (pfire = 1.0) ! for now we just assume it occurs all the time
      
c
c local variables
c


c ibis uses a small number of plant functional types:
c
c  1: tropical broadleaf evergreen tree
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen tree
c  4: temperate conifer evergreen tree
c  5: temperate broadleaf cold-deciduous tree
c  6: boreal conifer evergreen tree
c  7: boreal broadleaf cold-deciduous tree
c  8: boreal conifer cold-deciduous tree
c  9: evergreen shrub
c 10: deciduous shrub
c 11: warm (c4) grass
c 12: cool (c3) grass
c

c
c return to the main program
c
      return
      end
c
