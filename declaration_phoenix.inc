!---------------------------------------------------------------------!
!     CALCULATION STUFF
!---------------------------------------------------------------------!
      integer :: i,j
      real*8 :: Dqt ! rate of change of total mixing ratio with height [1/km]
!---------------------------------------------------------------------!
!     REAL, ACTUAL THINGS IN THE WORLD
!---------------------------------------------------------------------!
      real*8, dimension(:), allocatable :: zed ! elevation [km]
      real*8, dimension(:), allocatable :: Pz, Tz ! interpolated pressure/temperature profile [bars, K]
      real*8, dimension(1 : count3) :: radyuss ! particle radius [cm]
      real*8, dimension(1 : count3) :: vf ! sedimentation velocity [cm/s]
      real*8 :: rho_a ! atmospheric density [g/cm3]
      real*8 :: Hs ! atmospheric scale height [km]
      real*8 :: lapse_rate ! local temperature lapse rate dT/dz ![K/km]
      real*8 :: mix_L ! mixing length [km]
      real*8 :: diff_K ! eddy diffusion coefficient [cm2/s]
      real*8 :: w ! convective velocity [m/s]
!---------------------------------------------------------------------!
!     THE MIXING RATIOS
!---------------------------------------------------------------------!
      real*8, dimension(:), allocatable :: Qs ! saturation vapour mixing ratio
      real*8, dimension(:), allocatable :: Qv ! vapour mixing ratio
      real*8, dimension(:), allocatable :: Qc ! condensate mixing ratio
      real*8, dimension(:), allocatable :: Qt ! total mixing ratio {Qt = Qv + Qc}
      real*8 :: cond
!---------------------------------------------------------------------!
!     NOT SURE WHERE THESE BELONG
!---------------------------------------------------------------------!
      real*8 :: eta ! dynamic viscosity of atm [g/cm s] or [Poise]
      real*8 :: rw ! value of r @ wch vf = w [cm]
      real*8 :: alf !
      real*8 :: r_g ! geometric mean radius [cm] (units may be optimistic...)
      real*8 :: EN ! total number concentration of particles [1/cm3]
      real*8, dimension(:), allocatable :: dt
      real*8 :: anum ! dN/dln(r) the size distribution   [cm^-3 cm^-1]
