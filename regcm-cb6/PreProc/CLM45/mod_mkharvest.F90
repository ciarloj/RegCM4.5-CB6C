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
module mod_mkharvest
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_getwindow
  use mod_bilinear
  use mod_nchelper
  use mod_memutil
  use netcdf

  implicit none

  private

  public :: mkharvest

  integer , parameter :: nvarc = 6

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter , dimension(nvarc):: varname = &
          (/ 'HARVEST_VH1' , 'HARVEST_VH2' , 'HARVEST_SH1' , &
             'HARVEST_SH2' , 'HARVEST_SH3' , 'GRAZING    '/)
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rk8) :: vmin = -1.0D0
  real(rk8) :: vmisdat = -9999.0D0

  contains

  subroutine mkharvest(harvest,iyear)
    implicit none
    real(rk8) , dimension(:,:,:) , intent(out) :: harvest
    integer(ik4) :: iyear
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarmask , illvar , ncid
    integer(ik4) , dimension(6) :: ivarid
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) :: istatus , i , j , li , lo
    real(rk8) , dimension(:,:,:) , allocatable :: rvar
    real(rk8) , dimension(:,:) , allocatable :: rmask
    real(rk8) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain
    character(len=32) :: p1 , p2
    character(len=4) :: cy

    character(len=256) :: inpfile

    p1 = 'dynamic'
    p2 = '.'
    write(cy,'(i0.4)') iyear

    if ( iyear > 2005 ) then
      select case (dattyp(4:5))
        case ('RF')
          continue
        case ('26')
          p2 = 'SCENARIO'//pthsep//'RCP2.6'
        case ('45')
          p2 = 'SCENARIO'//pthsep//'RCP4.5'
        case ('60')
          p2 = 'SCENARIO'//pthsep//'RCP6.0'
        case ('85')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          call die(__FILE__,'Dynamic landuse only supported for CMIP5',__LINE__)
      end select
    end if

    inpfile = trim(inpglob)//pthsep//'CLM45'//pthsep//'surface'// &
            pthsep//trim(p1)//pthsep//trim(p2)//pthsep//&
            'mksrf_landuse_'//cy//'.nc'
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    if ( size(harvest,3) < nvarc ) then
      call die(__FILE__,'Size too small for harvest in mkharvest',__LINE__)
    end if

    istatus = nf90_inq_dimid(ncid,latdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lat in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lat in file '//trim(inpfile))

    allocate(glat(nlat))

    istatus = nf90_inq_varid(ncid,latvar,illvar)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,illvar,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,illvar)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,illvar,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lon in file '//trim(inpfile))

    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0D0 )
      glon = glon - 360.0D0
    end where
    where ( glon < -180.0D0 )
      glon = glon + 360.0D0
    end where

    call get_window(glat,glon,domain)

    allocate(rvar(sum(domain%ni),domain%nj,nvarc))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    do j = 1 , nvarc
      istatus = nf90_inq_varid(ncid,varname(j),ivarid(j))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find variable '//trim(varname(j))//' in file '//trim(inpfile))
    end do
    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      lo = li+domain%ni(i)-1
      do j = 1 , nvarc
        istatus = nf90_get_var(ncid,ivarid(j),rvar(li:lo,:,j),istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot read variable '//trim(varname(j))//&
          ' from file '//trim(inpfile))
      end do
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar,rmask,rlon,rlat,harvest,xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkharvest
#endif

end module mod_mkharvest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
