INCLUDE "sophia_legacy.f"

c****************************************************************************
c
c   SOPHIAEVENT
c
c   interface between Sophia and CRPropa
c   simulate an interaction between p/n of given energy and the CMB
c
c   Eric Armengaud, 2005
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sophiaeventmod(nature,Ein,eps,OutPart,OutPartType,NbOut
     &     Part)

c**********************************
c nature, Ein = input nature and energy of the nucleon
c        nature = 0 -> p ; 1 -> n
c        Ein : in GeV (SOPHIA standard energy unit)
c        eps : in GeV (SOPHIA standard energy unit)
c OutPart,OutPartType,NbOutPart = output data:
c        P(2000,5) list of 4-momenta + masses of output particles
c        LList(2000) list of output particle IDs
c        NP nb of output particles
c**********************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE
      
      COMMON/input/ tbb,E0,alpha1,alpha2,
     &     epsm1,epsm2,epsb,L0

      COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
      COMMON /S_MASS1/ AM(49), AM2(49)
      COMMON /S_CHP/  S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
      COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)
      
      CHARACTER*6 NAMPRES
      COMMON /RES_PROP/ AMRES(9), SIG0(9),WIDTH(9), 
     +                    NAMPRES(0:9)

      CHARACTER*6 NAMPRESp
      COMMON /RES_PROPp/ AMRESp(9), BGAMMAp(9),WIDTHp(9),  
     +                    RATIOJp(9),NAMPRESp(0:9)

      CHARACTER*6 NAMPRESn
      COMMON /RES_PROPn/ AMRESn(9), BGAMMAn(9),WIDTHn(9),  
     +                    RATIOJn(9),NAMPRESn(0:9)


      external sample_s,eventgen,initial
      
      integer nature
      double precision Ein,Pp,eps
      double precision OutPart(2000,5)
      integer OutPartType(2000)
      integer NbOutPart

      DATA pi /3.141593D0/

      if (nature.eq.0) then 
         L0=13
      else if (nature.eq.1) then
         L0=14
      else
         print*,'sophiaevent: incoming particle incorrectly specified'
         stop
      endif

      print*,'sophiaeventmod version'

      call initial(L0)

      E0 = Ein
      pm = AM(L0)

      call sample_s(s,eps)

      Pp = sqrt(E0*E0-pm*pm)
      theta = ((pm*pm-s)/2.D0/eps+E0)/Pp
      if (theta.gt.1.D0) then
         print*,'sophiaevent: theta > 1.D0: ', theta
         theta = 0.D0
      else if (theta.lt.-1.D0) then
         print*,'sophiaevent: theta < -1.D0: ', theta
         theta = 180.D0
      else
          theta = acos(theta)*180.D0/pi 
      endif

      call eventgen(L0,E0,eps,theta,Imode)

      do i=1,2000
         do j=1,5
            OutPart(i,j)=P(i,j)
         end do
         OutPartType(i)=LLIST(i)
      end do
      NbOutPart=NP

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc