      module fluffycloud
        use phoenix_variables
        implicit none
C        include 'physconst.inc' ! Rgas is in erg K-1 mol-1 (convert to J by x10^-7)
C        include 'cloudvars_phoenix.inc'
C        include 'declaration_phoenix.inc'
        include 'cloudvars.inc'
        include 'declaration.inc'
C        include 'param.inc'
!
        real*8, dimension (:), allocatable :: teplota,telak ! temperature, pressure
        real*8, dimension (:), allocatable :: stredni,cislo ! mean radius, number density
!
C        common /snmod/ wltau,teff,rtau1,vtau1,vfold                              !ALLOC-OFF
C     &/snmod/ tk(layer),the(layer),pg(layer),pe(layer),radius(layer)    !ALLOC-OFF
C     &/snmod/ rhoexp                                                    !ALLOC-OFF
C     &/snmod/ identyp,vswfac                                            !ALLOC-OFF
        private
        public :: Pz, Tz, layers, zed, Qc, Qv, Qt, Qs
        public :: teplota,telak,stredni,cislo
        public :: FLUFFY
        contains
!----------------------------------------------------------------------!
!       Find temp at given alt by interpolation mark dos               !
!----------------------------------------------------------------------!
      function TatZ2(Pin) result(Tout)
      real*8, intent(in) :: Pin ! input
      real*8 :: prod,sum
      integer :: j,k,lo
      integer :: l ! only a location
      real*8 :: Tout ! output
      l = minloc(abs(telak-Pin),1) ! find pressure value closest to input
!
      sum = 0.d0
      if ( l .eq. (size(telak)) ) then
        lo = l-4
      elseif ( l .eq. (size(telak)-1) ) then
        lo = l-3
      elseif ( l .eq. 2) then
        lo = l-1
      elseif ( l .eq. 1) then
        lo = l
      else
        lo = l-2
      endif
!
      do k = 0, 4
        prod = 1.d0
        do j = 0, 4
          if ( j .eq. k ) then
            cycle
          endif
          prod = prod*((Pin-telak(lo+j))/(telak(lo+k)-telak(lo+j)))
        enddo
        prod = prod*teplota(lo+k)
        sum = sum + prod
      enddo
!
      Tout = sum
!
      end function TatZ2
!----------------------------------------------------------------------!
!       Find local temperature gradient                                !
!       This is probably redundant now                                 !
!----------------------------------------------------------------------!
      function GRADI(Pz) result(tau)
      real*8, intent(in) :: Pz ! input
      real*8 :: tau ! output
      integer :: l ! only a location
      l = minloc(abs(telak-Pz),1)
      if ( abs(telak(l+1)-Pz) <= abs(telak(l-1)-Pz) ) then
        tau = (teplota(l+1)-teplota(l))/(telak(l+1)-telak(l))
      elseif ( abs(telak(l+1)-Pz) > abs(telak(l-1)-Pz) ) then
        tau = (teplota(l)-teplota(l-1))/(telak(l)-telak(l-1))
      end if
      end function GRADI
!----------------------------------------------------------------------!
!    Find eddy diffusion coefficient at given alt                      !
!    This is Equation 5 from Ackerman & Marley 2001                    !
!----------------------------------------------------------------------!
      function EDDY(H, len, ro) 	  result(Kout)
      real*8, intent(in) :: H, len, ro
      real*8, parameter :: F = (stefco*1.D-3)*(Tf_jup)**4.        ! Flux [J/m2 s]
      real*8, parameter :: g_r = Rgas*1.D-7
      real*8 :: Kout
      Kout = H*1.D5/3.*((len/H)**(4./3.))*(1.D3*g_r*F/
     & (mu_jup*ro*cp_jup))**(1./3.)
      end function EDDY
!----------------------------------------------------------------------!
!       A function to calculate the saturation vapour mixing ratio for !
!       a given temperature and ambient pressure                       !
!       Currently this is specific to ammonia                          !
!----------------------------------------------------------------------!
      function SatVap(T,P)  result(qs)
        real*8, intent(in) :: T, P  ! input temperature, pressure    !
        real*8 :: qs    ! output is saturation vapour MR !
        qs = (exp(10.53 - 2161./T - 86596./T**2.))/P
      end function SatVap
!----------------------------------------------------------------------!
!       find sedimentation velocities                                  !
!----------------------------------------------------------------------!
      subroutine VSED(w, Tin, Pin, r, vs)
        real*8, intent(in) :: w  ! input mixing velocity
        real*8, dimension(1 : count3), intent(out) :: vs, r ! output sedimentation velocity
        real*8 :: x,y,Tin,Pin,dr,dl
        integer :: j
!       dl : Knudsen number
!       dr : difference in condensate and atm. densities
!
        dl = 1.d-6*Rgas*Tin/(sqrt(2.d0)*pi*d_h2**2*Navog*Pin) ! molec. MFP in [cm]
        dr = rho_amm-rho_a ! density contrast btwn condensate & atm [g/cm^3]
        rw = (-1.26*dl+sqrt((1.26*dl)**2+4*w*eta*9./(2.*g_jup*dr)))/2.
        do j=1,count3
          r(j) = rw/sig_AM + rw*(1.-1./sig_AM)*(j-1.)/(count3-1.)
          x = log(32.d2*rho_a*g_jup*dr*r(j)**3/(3.*eta**2))
          y = .8*x-.01*x**2
            if ( dexp(y) <= 1000 ) then
            vs(j) = dexp(y)*eta/(2.*r(j)*rho_a)
            else
            vs(j) = (1+1.26*dl/r(j))
     & *sqrt(8.d2*g_jup*r(j)*dr/(1.35*rho_a))
            endif
        enddo
      end subroutine VSED
!----------------------------------------------------------------------!
!     find the exponent of the exponential function
!     y = x^a
!----------------------------------------------------------------------!
      function ALPH(x, y) result(a)
         real*8, dimension(1 : count3), intent(in) :: x, y
         real*8 :: a
        a = size(x)*sum(log(x)*log(y))-sum(log(x))*sum(log(y))
        a = a/(size(x)*sum((log(x))**2)-(sum(log(x)))**2)
      end function ALPH
!----------------------------------------------------------------------!
!     This is a model to compute the vertical size distribution and number
!     density of a given condensate. It is based on the model of Ackerman &
!     Marley (2001).
!     A pressure-temperature profile is input.
!     Requires mixing ratio of condensible vapour

!     Written by Jessie Brown
!     Version control begins FEBRUARY 14, 2017
!     JUNE 20, 2017: Converting from program to subroutine
!     JULY 11, 2017: Integrating into PHOENIX
!----------------------------------------------------------------------!
!
      subroutine FLUFFY(meansize,numberden)
        implicit none
!
        real*8, dimension(:), allocatable :: meansize
        real*8, dimension(:), allocatable :: numberden
!       temperature and pressure received from phoenix via phoenix_variables
C        teplota = temp
C        telak = pg
!
!     time to allocate my fancy arrays
      allocate( zed(0) )
      allocate( dt(0) )
      allocate( Qc(0) )
      allocate( Qs(0) )
      allocate( Qv(0) )
      allocate( Qt(0) )
      allocate( Pz(0) )
      allocate( Tz(0) )
      allocate( meansize(0) )
      allocate( numberden(0) )
!
C      open (unit=20,file="cloudout.txt",action="write",status="replace") ! FOR TESTING
!
!     give some arrays some initial values (ideally from PHOENIX input...)
      Tz = [ Tz, teplota(size(teplota)) ]
      Pz = [ Pz, telak(size(telak)) ]
      Qv = [ Qv, qbelow ]
      Qc = [ Qc, 0.d0 ]
      Qt = [ Qt, (Qv(1) + Qc(1)) ]
      dt = [ dt, 0.d0 ]
      zed = [ zed, 0.d0 ]
!
      Hs = (Rgas*1.D-7)*Tz(1)/(mu_jup*g_jup)
!
!     NOW CALCULATE THE LOT
      i=1
      do while ( Pz(i) .GE. telak(2) ) !j = 1, 4 ! ( Tz(i) .GE. 162 )!
        zed = [ zed, i*Dz ]
!       ESTIMATING PRESSURE & DENSITY USING HYDROSTATIC EQUILIBRIUM
        Pz = [ Pz, Pz(size(Pz))*exp(-Dz/Hs) ]
        rho_a = rho_0*exp(-zed(size(zed))/Hs)
!       GETTING TEMPERATURE FROM MEASURED PROFILE BY INTERPOLATION
        Tz = [ Tz, TatZ2(Pz(size(Pz))) ]
C        write (20,*) i, Pz(i), TatZ2(Pz(i)), Tz(i) ! FOR TESTING
!       FIND LOCAL SCALE HEIGHT ... makes the match to A&M better
        Hs = (Rgas*1.D-7)*Tz(size(Tz))/(mu_jup*g_jup)
!       CALCULATE LOCAL TEMPERATURE LAPSE RATE dT/dz
!       USING
!       dT/dz = [dT/dP]*[dP/dz] = [dT/dP]*[-rho_a * g_jup]
C        lapse_rate = GRADI(Pz(size(Pz)))*(-rho_a*g_jup)*1.d1 ! I think the new method should be better..?
        lapse_rate = (Pz(size(Pz))-Pz(size(Pz)-1))
     & /(Tz(size(Tz))-Tz(size(Tz)-1))
        lapse_rate = lapse_rate*(-rho_a*g_jup)*1.d1
!       CALCULATE LOCAL MIXING LENGTH
        mix_L = Hs*max(A,lapse_rate/Tau_a)
!       CALCULATE EDDY DIFFUSION COEFFICIENT
        diff_K = EDDY(Hs,mix_L,rho_a)
         if ( diff_K < Kmin ) then
	        diff_K = Kmin ! have prescribed a minimum value for eddy diff. coeff.
       	 endif
!       FIND THE CONVECTIVE VELOCITY FROM MIXING LENGTH THEORY
        w = (diff_K/mix_L)*(1.D-7)
!
!       AT LAST, THE GOOD BIT
!       SATURATION MIXING RATIO
        Qs = [ Qs, SatVap(Tz(size(Tz)),Pz(size(Pz)))]
!       ACTUAL MIXING RATIO OF VAPOUR => EQ 2 OF ACKERMAN & MARLEY 2001
        Qv = [ Qv, (min( Qv(size(Qv)), (1.D0+Sc)*Qs(size(Qs)) )) ]
!       FIND AMOUNT OF "NEW" CONDENSATE BORN IN THIS LAYER
        cond = max( 0.d0, Qv(size(Qv))-(1.D0+Sc)*Qs(size(Qs)) )
!       FIND RATE OF CHANGE OF MIXING RATIO WITH HEIGHT
        Dqt = -(frain*w*(Qc(size(Qc))+cond)/(diff_K))*1.d4
!       UPDATE THE TOTAL MIXING RATIO
        Qt = [ Qt, (Qt(size(Qt)) + Dqt*Dz*1.d3) ]
!       FINALLY, SET THE NEW MIXING RATIO OF CONDENSATE
        Qc = [ Qc, (Qt(size(Qt))-Qv(size(Qv))) ]
!
!       NOW THE SECOND PART
        if ( Qc(size(Qc)) > 0.d0 ) then
          eta = ((kb*Tz(size(Tz))/eps_jup)**0.16)/(pi*(d_h2**2)*1.22)
          eta = (5./16.)*sqrt(pi*em*kb*Tz(size(Tz))/Navog)*eta
          call VSED(w,Tz(size(Tz)),Pz(size(Tz)),radyuss,vf)
          alf = ALPH(radyuss/rw,vf/(w*1.d2))
!       FIND GEOMETRIC MEAN RADIUS
          r_g = rw*frain**(1./alf)*exp(-(alf+6.)*(log(sig_AM))**2./2.)
          meansize = [ meansize, r_g ]
!       FIND TOTAL NUMBER CONCENTRATION OF PARTICLES
          EN = (3.*(mu_a/mu_jup)*rho_a*Qc(size(Qc))
     & /(4.*pi*rho_amm*r_g**3))
          EN = EN*exp(-9.*(log(sig_AM))**2./2.)
          numberden = [ numberden, EN ]
        else
          meansize = [ meansize, 0.d0 ]
          numberden = [ numberden, 0.d0 ]
        endif
!
C      write (20,*) i, Pz(i), Tz(i), meansize(i) ! FOR TESTING
      i = i+1
      enddo
C      do i = 1, size(temp)
C        print *, i, temp(i)! Pz(i), Tz(i), meansize(i)
C        write (6,*) i, temp(i)!
C      enddo
C      close (20) ! FOR TESTING
!
      end subroutine FLUFFY
!
      end module fluffycloud
