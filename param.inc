c Time-stamp: <1998-05-11 15:00:55 baron>
c--
c-- number of layers
c--
      integer,parameter :: layer=64
      integer,parameter :: layerend=100
      integer, parameter ::  tma=layer
cjunk      parameter ( tma=247)

c--
c-- EOS parameters:
c--
c-- mist: min. ionization stage
c-- mxst: max. ionization stage
c-- z_max: max element Z value.
c--
      integer,parameter :: mist=-1,mxst=30,z_max=111
c
      integer,parameter :: nvemx=z_max+1,nionmx=30*z_max,nmolmx=1000
      integer, parameter :: npdmx=nvemx + nionmx + nmolmx 
      integer,parameter :: nelmx=nvemx + nionmx 
      integer,parameter :: nvemx1=nvemx - 1 
      integer,parameter :: nvionmx=nvemx + nionmx
      integer,parameter :: nopmax=150
      integer,parameter :: nopc=84 
      integer,parameter :: nrest=nopmax-nopc
      integer,parameter :: ntab=72+nrest 
      integer,parameter :: nkmax= 162 
      integer,parameter :: nrall=13 
      integer,parameter :: nsrt=npdmx,nvsrt=nvemx
clarge      integer,parameter :: nthdim=4096,npgdim=512
      integer,parameter :: nthdim=1000,npgdim=96
c nmolmx: total number of molecular + molec ions + condensate species available
c--
c-- number of elements for "useful" block:
c--
      integer,parameter :: maxel=103
c--
c-- JOLA parameter
c--
      integer,parameter :: nbmax=3000
      integer,parameter :: nbmol=10 
c nbmax : maximum number of molecular bands
c nbmol : maximum number of molecules with bands data
c--
c-- number of grains for which we have opacities:
c--
c-- ahw 3-24-2005 changes made to reflect addition of Akemi's data
c--   originals:  ndust=37; ndust_opac=46
      integer,parameter :: ndust=43
      integer,parameter :: npoints=40
      integer,parameter :: ndust_opacities=52
c-- jwf 19 Aug, 2004 added to make the number of n and k data scalable
      integer,parameter :: ndust_wave=6000
c--
c-- DOME stuff:
c--
c--
c-- parameter for the dome radiative transfer (mev.f):
c-----------------------------------------------------
c-- muemax: maximum number of angles per half-sphere
c-- maxpol: maximum degree of planck-polynomial (.le. 2)
c-- mueomx: maximum number of euler-frame angles
c-- layst: number of layer for j,h and k storage (debug)
c-- nlt   : total no. of species
c-- nlimax: max. no. of nlte lines explicitly handled
      integer,parameter :: maxpol=2
      integer,parameter :: muemax=4,lmax=1
      integer,parameter :: np=muemax*lmax,iw=np !,np2=2*np
c     -- max. number of wavelength points in cm & obs frames
c     -- max number of points for long wavelength approximation
      integer,parameter :: maxcm=5000000,maxobs=110000,maxapp=1000
c--
c-- junk
c--
      integer,parameter :: ntl=37,nlimax=100
