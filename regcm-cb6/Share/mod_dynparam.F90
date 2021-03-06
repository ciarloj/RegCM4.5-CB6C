!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
module mod_dynparam

  use mod_stdio
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
  use netcdf

  implicit none

  public
  !
  ! PARAMETER definitions
  !
  integer(ik4) , parameter :: ipunit = 255
  !
  !################### GRID DIMENSION ####################################
  !

  ! Point in Y (latitude) direction

  integer(ik4) :: iy

  ! Point in X (longitude) direction

  integer(ik4) :: jx

  ! Point in vertical

  integer(ik4) :: kz

  ! If not 14 , 18 or 23 (precalculated), hint for custom calculation

  real(rk8) :: dsmax , dsmin

  ! Sub grid decomposition

  integer(ik4) :: nsg

  ! Dynamical core

  integer(ik4) :: idynamic

  ! Projection
  !
  ! One in : 'LAMCON', Lambert conformal
  !          'POLSTR', Polar stereographic
  !          'NORMER', Normal  Mercator (ROTMER w/ plat = clat
  !          'ROTMER', Rotated Mercator
  !
  character(len=6) :: iproj

  ! Control flag for tropical band option.

  integer(ik4) :: i_band

  ! Control flag for creating bathymetry for lake model
  !    (Hostetler, etal. 1991, 1993a,b, 1995)

  logical :: lakedpth = .false.

  ! Control flag for crating teture dataset for aerosol dust

  logical :: ltexture = .false.

  ! Control flag for crating initial soil moisture dataset

  logical :: lsmoist = .false.

  ! Grid point horizontal resolution in km

  real(rk8) :: ds

  ! Pressure of model top in cbar

  real(rk8) :: ptop

  ! Central latitude  of model domain in degrees, north hem. is positive

  real(rk8) :: clat

  ! Central longitude of model domain in degrees, west is negative

  real(rk8) :: clon

  ! Pole latitude (only for rotated Mercator Proj, else set = clat)

  real(rk8) :: plat

  ! Pole longitude (only for rotated Mercator Proj, else set = clon)

  real(rk8) :: plon

  ! Lambert / Polar Cone factor

  real(rk8) :: xcone

  ! Lambert true latitude (low latitude side)

  real(rk8) :: truelatl

  ! Lambert true latitude (high latitude side)

  real(rk8) :: truelath

  ! Smoothness level

  integer(ik4) :: ismthlev

  !###################### DEBUG I/O control flag #########################

  ! Set amount of printout (still unused, sorry)

  integer(ik4) :: debug_level = 0
  integer(ik4) :: dbgfrq = 3600

  !###################### I/O control flag ###############################

  ! Buffer Zone Depth
  ! nspgx-1,nspgd-1 represent the number of cross/dot point slices
  ! on the boundary sponge or relaxation boundary conditions.

  integer(ik4) :: nspgx = 12
  integer(ik4) :: nspgd = 12

  ! Nudge control coefficients

  real(rk8) :: high_nudge   = 3.0D0
  real(rk8) :: medium_nudge = 2.0D0
  real(rk8) :: low_nudge    = 1.0D0

  ! Number od split exp modes

  integer(ik4) :: nsplit

  ! Type of global analysis datasets used in Pre processing
  ! One in: ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
  !         NRP2W,GFS11,FVGCM,FNEST,EH5OM

  character(len=5) :: dattyp

  !Type of Global chemistry boundary conditions
  !      MZ6HR is for MOZART 6 hourly boundary conditions
  !      MZCLM is for MOZART climatology

  character(len=5) :: chemtyp

  ! Type of Sea Surface Temperature used
  ! One in: GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
  !         EH5A2,EH5B1,EHA1B,ERSST,ERSKT

  character(len=5) :: ssttyp

  ! Land Surface Legend number

  integer(ik4) :: nveg

  ! Tracer parameters: number of tracers

  integer(ik4) :: ntr    = 0 ! Total number of chemical tracers
  integer(ik4) :: ngstr  = 0 ! Total number of gas tracers for sorgam
  integer(ik4) :: nsoatr = 0 ! Total number of soa tracers

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of configureation. Below this point things are
  !    calculated from above or should be considered as fixed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer(ik4) :: iym1
  integer(ik4) :: iym2
  integer(ik4) :: iym3
  integer(ik4) :: jxm1
  integer(ik4) :: jxm2
  integer(ik4) :: jxm3
  integer(ik4) :: kzm1
  integer(ik4) :: kzm2
  integer(ik4) :: kzp1
  integer(ik4) :: kzp2
  integer(ik4) :: kzp3
  integer(ik4) :: kzp4
  integer(ik4) :: iysg
  integer(ik4) :: jxsg
  integer(ik4) :: iym1sg
  integer(ik4) :: jxm1sg
  integer(ik4) :: iym2sg
  integer(ik4) :: jxm2sg
  integer(ik4) :: iym3sg
  integer(ik4) :: jxm3sg
  integer(ik4) :: nnsg

  integer(ik4) :: njcross , njdot , njout , njoutsg
  integer(ik4) :: nicross , nidot , niout , nioutsg

  integer(ik4) :: jcross1 , icross1
  integer(ik4) :: jcross2 , icross2
  integer(ik4) :: jdot1 , idot1
  integer(ik4) :: jdot2 , idot2
  integer(ik4) :: jout1 , iout1
  integer(ik4) :: jout2 , iout2
  integer(ik4) :: joutsg1 , ioutsg1
  integer(ik4) :: joutsg2 , ioutsg2

  ! D stands for DOT
  integer(ik4) :: ide1 , ide2 ! External i (included bdy) (latitude)
  integer(ik4) :: jde1 , jde2 ! External j (included bdy) (longitude)
  integer(ik4) :: idi1 , idi2 ! Internal (excluded first and last line) i
  integer(ik4) :: jdi1 , jdi2 ! Internal (excluded first and last column) j
  integer(ik4) :: idii1 , idii2 ! Internal (excluded 2 lines and cols) i
  integer(ik4) :: jdii1 , jdii2 ! Internal (excluded 2 lines and cols) j

  ! C stands for CROSS
  integer(ik4) :: ice1 , ice2 ! External (included bdy) i (latitude)
  integer(ik4) :: jce1 , jce2 ! External (included bdy) j (longitude)
  integer(ik4) :: ici1 , ici2 ! Internal (excluded first and last line) i
  integer(ik4) :: jci1 , jci2 ! Internal (excluded first and last column) j
  integer(ik4) :: icii1 , icii2 ! Internal (excluded 2 lines and cols) i
  integer(ik4) :: jcii1 , jcii2 ! Internal (excluded 2 lines and cols) j

  ! J index Dot points Full Domain  = jde1 : begin , jde2 : end
  ! I index Cross points Internal Domain = ici1 : begin , ici2 : end

  ! Global reference in global grid jx*iy of dot points
  ! The CROSS grid is contained within
  integer(ik4) :: global_dot_jstart
  integer(ik4) :: global_dot_jend
  integer(ik4) :: global_dot_istart
  integer(ik4) :: global_dot_iend
  integer(ik4) :: global_cross_jstart
  integer(ik4) :: global_cross_jend
  integer(ik4) :: global_cross_istart
  integer(ik4) :: global_cross_iend
  integer(ik4) :: global_out_jstart
  integer(ik4) :: global_out_jend
  integer(ik4) :: global_out_istart
  integer(ik4) :: global_out_iend

  !####################### MPI parameters ################################

  integer(ik4) :: mycomm
  integer(ik4) :: nproc
  integer(ik4) :: myid
  integer(ik4) :: njxcpus , niycpus
  integer(ik4) :: iyp , jxp
  integer(ik4) :: iypsg , jxpsg

  !####################### MPI parameters ################################

  ! Surface minimum H2O percent to be considered water

  real(rk8) :: h2opct

  ! Allow water pixels to have an elevation

  logical :: h2ohgt

  ! Smoothing Control flag
  !     true  -> Perform extra smoothing in boundaries

  logical :: smthbdy

  ! Fudging for landuse and texture for grid and subgrid

  logical :: fudge_lnd
  logical :: fudge_lnd_s
  logical :: fudge_tex
  logical :: fudge_tex_s
  logical :: fudge_lak
  logical :: fudge_lak_s

  ! Terrain output files

  character(len=64) :: domname

  ! Global Begin and End date for Input Pre processing

  type(rcm_time_and_date) , save :: globidate1 ! BEGIN
  type(rcm_time_and_date) , save :: globidate2 ! END

  ! Days per year and degrees per day

  character(len=12) :: calendar
  integer(ik4) :: ical
  real(rk8) :: dayspy
  real(rk8) :: half_dayspy
  real(rk8) :: sixteenth_dayspy
  real(rk8) :: dpd

  ! Fixed dimensions

  integer(ik4) , parameter :: mpy = 12         ! Months per Year

  ! Number of Soil texture categories, leave it to 17

  integer(ik4) , parameter :: ntex = 17
  integer(ik4) , parameter :: nats = 12 ! Should be ntex-5. Soil classes.

  ! Maximum number of depths in lake model

  integer(ik4) , parameter :: ndpmax = 200

  ! Number of bins in solar spectra

  integer(ik4) , parameter :: nspi = 19

#ifdef CLM45
  ! Soil layer thickness discretization (m)
  real(rk8) , parameter :: scalez = 0.025D0
  integer(ik4) , parameter :: num_soil_layers = 10
#else
  integer(ik4) , parameter :: num_soil_layers = 2
#endif

  ! Shall we use this to port?

  character(len=1), parameter :: pthsep = '/'

  ! Paths

  character(len=256) :: dirter , inpter
  character(len=256) :: dirglob , inpglob
  character(len=256) :: dirout
  character(len=256) :: moist_filename
  character(len=8)   :: tersrc
#ifdef NETCDF4_HDF5
  integer(ik4) :: iomode = ior(nf90_clobber, &
                               ior(nf90_netcdf4,nf90_classic_model))
  integer(ik4) :: deflate_level = 3
#else
  integer(ik4) :: iomode = ior(nf90_clobber, nf90_64bit_offset)
#endif

  ! Model output control parameters

  logical :: ifsave
  logical :: ifatm
  logical :: ifrad
  logical :: ifsrf
  logical :: ifsub
  logical :: ifsts
  logical :: iflak
  logical :: ifopt
  logical :: ifchem

  real(rk8) :: savfrq
  real(rk8) :: atmfrq
  real(rk8) :: radfrq
  real(rk8) :: lakfrq
  real(rk8) :: subfrq
  real(rk8) :: srffrq
  real(rk8) :: chemfrq

  integer(ik4) :: ibdyfrq

  logical :: ensemble_run

  logical :: lperturb_topo
  logical :: lperturb_ts
  logical :: lperturb_ps
  logical :: lperturb_t
  logical :: lperturb_q
  logical :: lperturb_u
  logical :: lperturb_v

  real(rk8) :: perturb_frac_topo
  real(rk8) :: perturb_frac_ts
  real(rk8) :: perturb_frac_ps
  real(rk8) :: perturb_frac_t
  real(rk8) :: perturb_frac_q
  real(rk8) :: perturb_frac_u
  real(rk8) :: perturb_frac_v

#ifdef CLM45
  logical :: enable_megan_emission = .false.
  logical :: enable_urban_landunit = .true.
  logical :: enable_more_crop_pft = .false.
#endif

  contains

  subroutine initparam(filename, ierr)
    implicit none
    character (len=*) , intent(in) :: filename
    integer(ik4) , intent(out) :: ierr
    integer(ik4) :: gdate1 , gdate2 , iresult

    namelist /dimparam/ iy , jx , kz , dsmax , dsmin , nsg , njxcpus , niycpus
    namelist /coreparam/ idynamic
    namelist /geoparam/ iproj , ds , ptop , clat , clon , plat ,    &
      plon , truelatl, truelath , i_band
    namelist /terrainparam/ domname , smthbdy , ltexture , lakedpth,  &
      lsmoist , fudge_lnd , fudge_lnd_s , fudge_tex , fudge_tex_s ,   &
      fudge_lak , fudge_lak_s , h2opct , h2ohgt , ismthlev , dirter , &
      inpter , moist_filename , tersrc
    namelist /debugparam/ debug_level , dbgfrq
    namelist /boundaryparam/ nspgx , nspgd , high_nudge , &
      medium_nudge , low_nudge
    namelist /globdatparam/ dattyp , chemtyp, ssttyp , gdate1 , gdate2 , &
      dirglob , inpglob , calendar , ibdyfrq , ensemble_run
    namelist /perturbparam/ lperturb_ts , perturb_frac_ts ,         &
      lperturb_topo , perturb_frac_topo ,         &
      lperturb_ps , perturb_frac_ps , lperturb_t , perturb_frac_t , &
      lperturb_q , perturb_frac_q , lperturb_u , perturb_frac_u ,   &
      lperturb_v , perturb_frac_v
#ifdef CLM45
    namelist /clm_regcm/ enable_megan_emission , enable_urban_landunit, &
      enable_more_crop_pft
#endif

    open(ipunit, file=filename, status='old', &
                 action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      ierr = 1
      return
    end if
!
    dsmax = 0.05D0
    dsmin = 0.0025D0
    njxcpus = -1
    niycpus = -1

    rewind(ipunit)
    read(ipunit, nml=dimparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading dimparam namelist in ',trim(filename)
      ierr = 2
      return
    end if
    if ( nsg < 1 ) then
      nsg = 1
    end if

    idynamic = 1
    rewind(ipunit)
    read(ipunit, nml=coreparam, iostat=iresult)

    i_band = 0
    rewind(ipunit)
    read(ipunit, nml=geoparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading geoparam namelist in ',trim(filename)
      ierr = 3
      return
    end if

!   Setup all convenience dimensions

    iym1 = iy - 1
    iym2 = iy - 2
    iym3 = iy - 3
    jxm1 = jx - 1
    jxm2 = jx - 2
    jxm3 = jx - 3
    kzm1 = kz - 1
    kzm2 = kz - 2
    kzp1 = kz + 1
    kzp2 = kz + 2
    kzp3 = kz + 3
    kzp4 = kz + 4
    iysg = iy * nsg
    jxsg = jx * nsg
    iym1sg = iym1 * nsg
    jxm1sg = jxm1 * nsg
    iym2sg = iym2 * nsg
    jxm2sg = jxm2 * nsg
    iym3sg = iym3 * nsg
    jxm3sg = jxm3 * nsg
    nnsg = nsg*nsg
    jdot1 = 1
    jdot2 = jx
    jcross1 = 1
    if ( i_band == 1 ) then
      jcross2 = jx
      jout1 = 1
      jout2 = jx
      joutsg1 = 1
      joutsg2 = jx*nsg
    else
      jcross2 = jxm1
      jout1 = 2
      jout2 = jxm2
      joutsg1 = nsg+1
      joutsg2 = jxm2*nsg
    end if
    idot1 = 1
    idot2 = iy
    icross1 = 1
    icross2 = iym1
    iout1 = 2
    iout2 = iym2
    ioutsg1 = nsg+1
    ioutsg2 = iym2*nsg
    njcross = jcross2-jcross1+1
    nicross = icross2-icross1+1
    njdot = jdot2-jdot1+1
    nidot = idot2-idot1+1
    njout = jout2-jout1+1
    niout = iout2-iout1+1
    njoutsg = joutsg2-joutsg1+1
    nioutsg = ioutsg2-ioutsg1+1

    nveg = 22

    if ( i_band.eq.1 ) then
      ds = (2.0D0*mathpi*erkm)/dble(jx)
      iproj = 'NORMER'
      clat  =   0.0D0
      clon  = 180.0D0
    end if

    ! Defaults to have SAME behaviour of V3 if not specified
    inpter  = '../DATA'
    inpglob = '../DATA'
    dirter  = '../../Input'
    dirglob = '../../Input'
    moist_filename = 'moist.nc'
    tersrc = 'GMTED'

    h2ohgt = .true.
    h2opct = 50.0D0
    ismthlev = 1
    rewind(ipunit)
    read(ipunit, nml=terrainparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading terrainparam namelist in ',trim(filename)
      ierr = 4
      return
    end if

    ! Set convenient defaults for debug I/O parameters
    rewind(ipunit)
    read(ipunit, nml=debugparam, iostat=iresult)
    if ( iresult /= 0 ) then
      ! We have defaults
      continue
    end if

    ! We assume 12 points is "OK" for a resolution of 50 km.
    ! Let us assume a "default" for the selected ds, not getting less
    ! than the OK number of points. If the domain is REALLY small,
    ! use 1/4 of the overall points, or AT LEAST 3 points...
    nspgx = max(min(max(int(dble(nspgx*50)/ds),nspgx),min(jx,iy)/4),3)
    nspgd = max(min(max(int(dble(nspgd*50)/ds),nspgd),min(jx,iy)/4),3)
    ! Anyway the user specify this...
    rewind(ipunit)
    read(ipunit, nml=boundaryparam, iostat=iresult)
    if ( iresult /= 0 ) then
      ! We have defaults
      continue
    end if
    ! Just double check ;)
    nspgx = max(nspgx,3)
    nspgd = max(nspgd,3)

    ! Removed modesparam. It must not be changed, so avoid cluttering.
    nsplit = 2

    ibdyfrq = 6 ! Convenient default
    calendar = 'gregorian'
    ensemble_run = .false.

    rewind(ipunit)
    read(ipunit, nml=globdatparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading globdatparam namelist in ',trim(filename)
      ierr = 6
      return
    end if
    if (calendar == 'gregorian') then
      dayspy = 365.2422D+00
      ical = gregorian
    else if (calendar == 'noleap' .or. calendar == '365_day') then
      dayspy = 365.0D+00
      ical = noleap
    else if (calendar == '360_day') then
      dayspy = 360.0D+00
      ical = y360
    else
      write(stderr,*) 'No calendar specified. Assuming gregorian'
      dayspy = 365.2422D+00
      ical = gregorian
    end if
    dpd = 360.0D0/dayspy
    half_dayspy = dayspy/2.0D0
    sixteenth_dayspy = dayspy/16.0D0
    globidate1 = gdate1
    globidate2 = gdate2
    call setcal(globidate1,ical)
    call setcal(globidate2,ical)
    if ( ensemble_run ) then
      lperturb_topo = .false.
      lperturb_ts = .false.
      lperturb_ps = .false.
      lperturb_t = .false.
      lperturb_q = .false.
      lperturb_u = .false.
      lperturb_v = .false.
      perturb_frac_topo = d_r1000
      perturb_frac_ts = d_r1000
      perturb_frac_ps = d_r1000
      perturb_frac_t = d_r1000
      perturb_frac_q = d_r1000
      perturb_frac_u = d_r1000
      perturb_frac_v = d_r1000
      rewind(ipunit)
      read(ipunit, nml=perturbparam, iostat=iresult)
      if ( iresult /= 0 ) then
        write (stderr,*) 'Error reading perturbparam namelist in ' , &
          trim(filename)
        ierr = 6
        return
      end if
    end if

#ifdef CLM45
    rewind(ipunit)
    read(ipunit, nml=clm_regcm, iostat=iresult)
    if ( iresult /= 0 ) then
      ! No error here. we have defualts. ;)
      continue
    end if
#endif

    close(ipunit)
    ierr = 0
    return
  end subroutine initparam

  subroutine init_fnestparam(filename,coarse_outdir,coarse_domname)
    implicit none
    character(len=*) , intent(in) :: filename
    character(len=256) , intent(out) :: coarse_outdir , coarse_domname
    integer(ik4) :: iresult
    namelist /fnestparam/ coarse_outdir , coarse_domname

    coarse_outdir = '        '
    coarse_domname = '        '
    open(ipunit, file=filename, status='old', &
                 action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      return
    end if
    read(ipunit, nml=fnestparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stdout,*) 'fnestparam not present: Using defaults'
    end if
    close(ipunit)
  end subroutine init_fnestparam

  subroutine init_globwindow(filename,lat0,lon0,lat1,lon1)
    implicit none
    character(len=*) , intent(in) :: filename
    real(rk8) , intent(out) :: lat0 , lat1 , lon0 , lon1
    integer(ik4) :: iresult
    namelist /globwindow/ lat0 , lat1 , lon0 , lon1

    lat0 = 0.0D0
    lon0 = 0.0D0
    lat1 = 0.0D0
    lon1 = 0.0D0

    open(ipunit, file=filename, status='old', &
                 action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      return
    end if
    read(ipunit, nml=globwindow, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stdout,*) 'Globwindow namelist not present: Assuming Global'
    end if
    close(ipunit)
  end subroutine init_globwindow

end module mod_dynparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
