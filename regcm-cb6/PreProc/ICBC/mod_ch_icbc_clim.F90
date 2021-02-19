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

module mod_ch_icbc_clim

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_ch_param
  use netcdf

  private
!
  integer(ik4) :: chilon , chjlat , chilev

  real(rk8) , pointer , dimension(:) :: cht42lon
  real(rk8) , pointer , dimension(:) :: cht42lat
  real(rk8) , pointer , dimension(:) :: cht42hyam , cht42hybm
!
! Oxidant climatology variables
!
  real(rk8) :: p0
  real(rk8) , pointer , dimension(:,:) :: pchem_3
  real(rk8) , pointer , dimension(:,:,:,:) :: chv3
  real(rk8) , pointer , dimension(:,:) :: xps
  real(rk8) , pointer , dimension(:,:,:,:) :: xinp
  real(rk8) , pointer , dimension(:,:,:,:) :: chv4_1
  real(rk8) , pointer , dimension(:,:,:,:) :: chv4_2
  real(rk8) , pointer , dimension(:,:,:,:) :: chv4_3

  real(rk8) :: prcm , pmpi , pmpj
  real(rk8) :: r4pt
  integer(ik4) :: ism
  type (rcm_time_and_date) , save :: iref1 , iref2

  public :: header_ch_icbc_clim , close_ch_icbc_clim
  public :: get_ch_icbc_clim , get_c6_icbc_clim

  contains

  subroutine header_ch_icbc_clim(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type (rcm_time_and_date) :: imonmidd
    integer(ik4) :: ivarid , idimid
    integer(ik4) :: nyear , month , nday , nhour
    character(len=256) :: chfilename
    integer(ik4) :: ncid , istatus
    integer(ik4) :: im1 , im2

    call split_idate(idate,nyear,month,nday,nhour)
    imonmidd = monmiddle(idate)
    im1 = month
    im2 = month
    if ( idate >= imonmidd ) then
      im2 = inextmon(im2)
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idate))
    else
      im1 = iprevmon(im1)
      iref1 = monmiddle(prevmon(idate))
      iref2 = imonmidd
    end if
    ism = im1

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_19990401.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
       'Error open file chemical '//trim(chfilename))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chjlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
    istatus = nf90_inq_dimid(ncid,'lev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lev')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilev)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lev')

    call getmem2d(pchem_3,1,jx,1,iy,'mod_ch_icbc:pchem_3_1')
    call getmem4d(chv3,1,jx,1,iy,1,chilev,1,nchsp,'mod_ch_icbc:chv3')
    call getmem4d(chv4_1,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_1')
    call getmem4d(chv4_2,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_2')
    call getmem4d(chv4_3,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_2')

    call getmem1d(cht42lon,1,chilon,'mod_ch_icbc:cht42lon')
    call getmem1d(cht42lat,1,chjlat,'mod_ch_icbc:cht42lat')
    call getmem1d(cht42hyam,1,chilev,'mod_ch_icbc:cht42hyam')
    call getmem1d(cht42hybm,1,chilev,'mod_ch_icbc:cht42hybm')
    call getmem2d(xps,1,chilon,1,chjlat,'mod_ch_icbc:xps1')
    call getmem4d(xinp,1,chilon,1,chjlat,1,chilev,1,nchsp,'mod_ch_icbc:xinp')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,cht42lon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,cht42lat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,cht42hyam)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,cht42hybm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var P0')
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')

    p0 = 1000.0
    r4pt = real(ptop)
    write(stdout,*) 'Static read OK.'

    call read2m(im1,im2)

  end subroutine header_ch_icbc_clim

  subroutine get_ch_icbc_clim(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: nyear , month , nday , nhour
    logical :: doread
    type (rcm_time_and_date) :: imonmidd
    type (rcm_time_interval) :: tdif
    real(rk8) :: xfac1 , xfac2 , odist
    integer(ik4) :: im1 , im2

    call split_idate(idate,nyear,month,nday,nhour)
    imonmidd = monmiddle(idate)
    im1 = month
    im2 = month
    if ( idate >= imonmidd ) then
      im2 = inextmon(im2)
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idate))
    else
      im1 = iprevmon(im1)
      iref1 = monmiddle(prevmon(idate))
      iref2 = imonmidd
    end if
    doread = .false.
    if ( ism /= im1 ) then
      ism = im1
      doread = .true.
    end if

    if ( doread ) then
      call read2m(im1,im2)
    end if

    tdif = idate-iref1
    xfac1 = tohours(tdif)
    tdif = iref2-idate
    xfac2 = tohours(tdif)
    odist = xfac1 + xfac2
    xfac1 = xfac1/odist
    xfac2 = d_one-xfac1
! rq: pppv(mozart) to pppm
    chv4_3 = (chv4_1*xfac2+chv4_2*xfac1)
chv4(:,:,:,cb_O3)    = chv4_3(:,:,:,mz_O3)*w_o3/amd
chv4(:,:,:,cb_NO)    = chv4_3(:,:,:,mz_NO)*w_no/amd
chv4(:,:,:,cb_NO2)   = chv4_3(:,:,:,mz_NO2)*w_no2/amd
chv4(:,:,:,cb_HNO3)  = chv4_3(:,:,:,mz_HNO3)*w_hno3/amd
chv4(:,:,:,cb_HNO4)  = chv4_3(:,:,:,mz_HO2NO2)*w_hno4/amd
chv4(:,:,:,cb_N2O5)  = chv4_3(:,:,:,mz_N2O5)*w_n2o5/amd
chv4(:,:,:,cb_H2O2)  = chv4_3(:,:,:,mz_H2O2)*w_h2o2/amd
chv4(:,:,:,cb_CH4)   = chv4_3(:,:,:,mz_CH4)*w_ch4/amd
chv4(:,:,:,cb_CO)    = chv4_3(:,:,:,mz_CO)*w_co/amd
chv4(:,:,:,cb_SO2)   = chv4_3(:,:,:,mz_SO2)*w_so2/amd
chv4(:,:,:,cb_H2SO4) = chv4_3(:,:,:,mz_SO4)*w_h2so4/amd
chv4(:,:,:,cb_DMS)   = chv4_3(:,:,:,mz_DMS)*w_dms/amd
chv4(:,:,:,cb_PAR)   = (3*chv4_3(:,:,:,mz_C3H8)+4*chv4_3(:,:,:,mz_BIGALK)+chv4_3(:,:,:,mz_C3H6)+chv4_3(:,:,:,mz_BIGENE))*w_par/amd
chv4(:,:,:,cb_C2H6)  = chv4_3(:,:,:,mz_C2H6)*w_c2h6/amd
chv4(:,:,:,cb_ETH)   = chv4_3(:,:,:,mz_C2H4)*w_eth/amd
chv4(:,:,:,cb_OLET)  = chv4_3(:,:,:,mz_BIGENE)*w_olet/amd
chv4(:,:,:,cb_OLEI)  = chv4_3(:,:,:,mz_BIGENE)*w_olei/amd
chv4(:,:,:,cb_TOL)   = chv4_3(:,:,:,mz_TOLUENE)*w_tol/amd
chv4(:,:,:,cb_XYL)   = chv4_3(:,:,:,mz_TOLUENE)*w_xyl/amd
chv4(:,:,:,cb_ISOP)  = chv4_3(:,:,:,mz_ISOP)*w_isop/amd
chv4(:,:,:,cb_CRES)  = chv4_3(:,:,:,mz_CRESOL)*w_cres/amd
chv4(:,:,:,cb_OPEN)  = chv4_3(:,:,:,mz_BIGALD)*w_open/amd
chv4(:,:,:,cb_ISOPN) = chv4_3(:,:,:,mz_ISOPNO3)*w_isopn/amd
chv4(:,:,:,cb_ISOPRD)= (chv4_3(:,:,:,mz_MVK)+chv4_3(:,:,:,mz_MACR)+chv4_3(:,:,:,mz_HYDRALD))*w_isoprd/amd
chv4(:,:,:,cb_ONIT)  = (chv4_3(:,:,:,mz_ONIT)+chv4_3(:,:,:,mz_ONITR))*w_onit/amd
chv4(:,:,:,cb_MGLY)  = chv4_3(:,:,:,mz_CH3COCHO)*w_mgly/amd
chv4(:,:,:,cb_AONE)  = (chv4_3(:,:,:,mz_CH3COCH3)+chv4_3(:,:,:,mz_HYAC)+chv4_3(:,:,:,mz_MEK))*w_aone/amd
chv4(:,:,:,cb_PAN)   = chv4_3(:,:,:,mz_PAN)*w_pan/amd
chv4(:,:,:,cb_CH3OOH)= chv4_3(:,:,:,mz_CH3OOH)*w_ch3ooh/amd
chv4(:,:,:,cb_ETHOOH)= chv4_3(:,:,:,mz_C2H5OOH)*w_ethooh/amd
chv4(:,:,:,cb_ALD2)  = (chv4_3(:,:,:,mz_CH3CHO)+chv4_3(:,:,:,mz_GLYALD))*w_ald2/amd
chv4(:,:,:,cb_HCHO)  = chv4_3(:,:,:,mz_CH2O)*w_hcho/amd
chv4(:,:,:,cb_CH3OH) = chv4_3(:,:,:,mz_CH3OH)*w_ch3oh/amd



    call write_ch_icbc(idate)

  end subroutine get_ch_icbc_clim

  subroutine get_c6_icbc_clim(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: nyear , month , nday , nhour
    logical :: doread
    type (rcm_time_and_date) :: imonmidd
    type (rcm_time_interval) :: tdif
    real(rk8) :: xfac1 , xfac2 , odist
    integer(ik4) :: im1 , im2

    call split_idate(idate,nyear,month,nday,nhour)
    imonmidd = monmiddle(idate)
    im1 = month
    im2 = month
    if ( idate >= imonmidd ) then
      im2 = inextmon(im2)
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idate))
    else
      im1 = iprevmon(im1)
      iref1 = monmiddle(prevmon(idate))
      iref2 = imonmidd
    end if
    doread = .false.
    if ( ism /= im1 ) then
      ism = im1
      doread = .true.
    end if

    if ( doread ) then
      call read2m(im1,im2)
    end if

    tdif = idate-iref1
    xfac1 = tohours(tdif)
    tdif = iref2-idate
    xfac2 = tohours(tdif)
    odist = xfac1 + xfac2
    xfac1 = xfac1/odist
    xfac2 = d_one-xfac1
    ! rq: pppv(mozart) to pppm
    chv4_3 = (chv4_1*xfac2+chv4_2*xfac1)
    !Lumping of Mozart species to CB6C species
!   chv4(:,:,:,cb6_OZN)   = (0.6D0*chv4_3(:,:,:,mz_O3))*w_ozn/amd
    chv4(:,:,:,cb6_OZN)   = chv4_3(:,:,:,mz_O3)*w_ozn/amd
    chv4(:,:,:,cb6_NMOX)  = chv4_3(:,:,:,mz_NO)*w_nmox/amd
!   chv4(:,:,:,cb6_NDOX)  = (0.8D0*chv4_3(:,:,:,mz_NO2))*w_ndox/amd
    chv4(:,:,:,cb6_NDOX)  = chv4_3(:,:,:,mz_NO2)*w_ndox/amd
    chv4(:,:,:,cb6_NTRC)  = chv4_3(:,:,:,mz_HNO3)*w_ntrc/amd
    chv4(:,:,:,cb6_PNA)   = chv4_3(:,:,:,mz_HO2NO2)*w_pna/amd
    chv4(:,:,:,cb6_DNPO)  = chv4_3(:,:,:,mz_N2O5)*w_dnpo/amd
    chv4(:,:,:,cb6_HPOX)  = chv4_3(:,:,:,mz_H2O2)*w_hpox/amd
    chv4(:,:,:,cb6_MTHA)  = chv4_3(:,:,:,mz_CH4)*w_mtha/amd
    chv4(:,:,:,cb6_CMON)  = chv4_3(:,:,:,mz_CO)*w_cmon/amd
    chv4(:,:,:,cb6_SDIO)  = chv4_3(:,:,:,mz_SO2)*w_sdio/amd
    chv4(:,:,:,cb6_SULF)  = chv4_3(:,:,:,mz_SO4)*w_sulf/amd
!   chv4(:,:,:,cb6_PRPA)  = (3*chv4_3(:,:,:,mz_C3H8)+4*chv4_3(:,:,:,mz_BIGALK)+chv4_3(:,:,:,mz_C3H6)+chv4_3(:,:,:,mz_BIGENE))*w_alka/amd
    chv4(:,:,:,cb6_PRPA)  = (3.0D0*chv4_3(:,:,:,mz_C3H8))*w_prpa/amd
!  chv4(:,:,:,cb6_ALKA)  = (3*chv4_3(:,:,:,mz_C3H8)+4*chv4_3(:,:,:,mz_BIGALK)+chv4_3(:,:,:,mz_C3H6)+chv4_3(:,:,:,mz_BIGENE))*w_alka/amd
!   chv4(:,:,:,cb6_ALKA)  = chv4_3(:,:,:,mz_BIGALK)*w_alka/amd ! CBMZ PAR mod
    chv4(:,:,:,cb6_ETHA)  = chv4_3(:,:,:,mz_C2H6)*w_etha/amd
    chv4(:,:,:,cb6_ETHE)  = chv4_3(:,:,:,mz_C2H4)*w_ethe/amd
    chv4(:,:,:,cb6_IOLE)  = chv4_3(:,:,:,mz_BIGENE)*w_iole/amd
!   chv4(:,:,:,cb6_OLE)   = chv4_3(:,:,:,mz_C3H6)*w_ole/amd
    chv4(:,:,:,cb6_OLE)   = chv4_3(:,:,:,mz_BIGENE)*w_ole/amd
    chv4(:,:,:,cb6_TOLN)  = chv4_3(:,:,:,mz_TOLUENE)*w_toln/amd
    chv4(:,:,:,cb6_XYLN)  = chv4_3(:,:,:,mz_TOLUENE)*w_xyln/amd
!   chv4(:,:,:,cb6_ISPR)  = chv4_3(:,:,:,mz_ISOP)*w_ispr/amd
    chv4(:,:,:,cb6_CRSL)  = chv4_3(:,:,:,mz_CRESOL)*w_crsl/amd
    chv4(:,:,:,cb6_ROPN)  = chv4_3(:,:,:,mz_BIGALD)*w_ropn/amd
!   chv4(:,:,:,cb6_NTR2)  = chv4_3(:,:,:,mz_ISOPNO3)*w_ntr2/amd ! CBMZ likely ISOPN
    chv4(:,:,:,cb6_NTR2)  = (chv4_3(:,:,:,mz_ONIT)+chv4_3(:,:,:,mz_ONITR))*w_cron/amd
    chv4(:,:,:,cb6_ISPD)  = (4.0D0*(chv4_3(:,:,:,mz_MVK)+chv4_3(:,:,:,mz_MACR)+chv4_3(:,:,:,mz_HYDRALD)))*w_ispd/amd ! CBMZ likely ISOPRD
!   chv4(:,:,:,cb6_CRON)  = (chv4_3(:,:,:,mz_ONIT)+chv4_3(:,:,:,mz_ONITR))*w_cron/amd ! CBMZ ONIT
    chv4(:,:,:,cb6_MEGY)  = chv4_3(:,:,:,mz_CH3COCHO)*w_megy/amd
    chv4(:,:,:,cb6_ACET)  = (chv4_3(:,:,:,mz_CH3COCH3)+chv4_3(:,:,:,mz_HYAC)+chv4_3(:,:,:,mz_MEK))*w_acet/amd ! CBMZ AONE original
!   chv4(:,:,:,cb6_ACET)  = chv4_1(:,:,:,mz_CH3COCH3)*w_acet/amd ! CBMZ AONE mod
    chv4(:,:,:,cb6_PACN)  = chv4_3(:,:,:,mz_PAN)*w_pacn/amd
!   chv4(:,:,:,cb6_AACD)  = chv4_3(:,:,:,mz_CH3COOH)*w_aacd/amd
    chv4(:,:,:,cb6_AALD)  = (chv4_3(:,:,:,mz_CH3CHO)+chv4_3(:,:,:,mz_GLYALD))*w_aald/amd ! CBMZ ALD2 original
!   chv4(:,:,:,cb6_AALD)  = chv4_1(:,:,:,mz_CH3CHO)*w_aald/amd ! CBMZ ALD2 mod
    chv4(:,:,:,cb6_FORM)  = chv4_3(:,:,:,mz_CH2O)*w_form/amd
    chv4(:,:,:,cb6_MEOH)  = chv4_3(:,:,:,mz_CH3OH)*w_meoh/amd

    chv4(:,:,:,cb6_GLYD)  = (4.0D0*chv4_3(:,:,:,mz_GLYALD))*w_glyd/amd
    chv4(:,:,:,cb6_MEPX)  = chv4_3(:,:,:,mz_CH3OOH)*w_mepx/amd
!!   chv4(:,:,:,cb6_NTR2)  = chv4_3(:,:,:,mz_ONIT)*w_ntr2/amd
!!   chv4(:,:,:,cb6_INTR)  = chv4_3(:,:,:,mz_ONIT)*w_intr/amd
    chv4(:,:,:,cb6_INTR)  = chv4_3(:,:,:,mz_ISOPNO3)*w_ntr2/amd 
!    chv4(:,:,:,cb6_KET)   = chv4_3(:,:,:,mz_MEK)*w_ket/amd
!    chv4(:,:,:,cb6_ISPD)  = chv4_3(:,:,:,mz_MACR)*w_ispd/amd
!!   chv4(:,:,:,cb6_HPLD)  = chv4_3(:,:,:,mz_XOOH)*w_hpld/amd
!    chv4(:,:,:,cb6_XOPN)  = chv4_3(:,:,:,mz_BIGALD)*w_xopn/amd
!    chv4(:,:,:,cb6_NTR1)  = chv4_3(:,:,:,mz_ONITR)*w_ntr1/amd
!    !chv4(:,:,:,cb6_O)     = chv4_3(:,:,:,mz_O)*w_o/amd
!    !chv4(:,:,:,cb6_OSNG)  = chv4_3(:,:,:,mz_O1D)*w_osng/amd
!!   chv4(:,:,:,cb6_NTOX)  = chv4_3(:,:,:,mz_NO3)*w_ntox/amd
!!   chv4(:,:,:,cb6_HOX)   = chv4_3(:,:,:,mz_OH)*w_hox/amd
!!   chv4(:,:,:,cb6_POX)   = chv4_3(:,:,:,mz_HO2)*w_pox/amd
!!   chv4(:,:,:,cb6_MEO2)  = chv4_3(:,:,:,mz_CH3O2)*w_meo2/amd
    chv4(:,:,:,cb6_ETOH)  = chv4_3(:,:,:,mz_C2H5OH)*w_etoh/amd
!!   chv4(:,:,:,cb6_ACOO)  = chv4_3(:,:,:,mz_CH3CO3)*w_acoo/amd
!!   chv4(:,:,:,cb6_ROO)   = chv4_3(:,:,:,mz_ENEO2)*w_roo/amd
!!   chv4(:,:,:,cb6_ISO2)  = chv4_3(:,:,:,mz_ISOPO2)*w_iso2/amd
!!   chv4(:,:,:,cb6_ISPX)  = chv4_3(:,:,:,mz_ISOPOOH)*w_ispx/amd
    chv4(:,:,:,cb6_TERP)  = chv4_3(:,:,:,mz_C10H16)*w_terp/amd
!!   chv4(:,:,:,cb6_TOLR)  = chv4_3(:,:,:,mz_TOLO2)*w_tolr/amd
!!   chv4(:,:,:,cb6_XYLR)  = chv4_3(:,:,:,mz_TOLO2)*w_xylr/amd
!!   chv4(:,:,:,cb6_GLY)   = chv4_3(:,:,:,mz_GLYOXAL)*w_gly/amd
    chv4(:,:,:,cb6_BENZ)  = (0.29D0*chv4_3(:,:,:,mz_TOLUENE))*w_benz/amd

    call write_c6_icbc(idate)
  end subroutine get_c6_icbc_clim


  subroutine read2m(im1,im2)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2
    integer(ik4) :: i , is , j , k , l , k0
    character(len=256) :: chfilename
    real(rk8) :: wt1 , wt2
    integer(ik4) :: ncid , istatus , ivarid

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_avg_1999-2009_',im1,'.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open file chemical')
    write(stdout, *) trim(chfilename)
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    xps = xps*0.01
    do is = 1 , nchsp
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_inst',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//trim(chspec(is)))
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//trim(chspec(is)))
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev)
    end do
    call bilinx2(pchem_3,xps,xlon,xlat,cht42lon,cht42lat, &
                 chilon,chjlat,iy,jx)
    do i = 1 , iy
      do j = 1 , jx
        do l = 1 , kz
          prcm=((pchem_3(j,i)*0.1-r4pt)*sigmah(l)+r4pt)*10.
          k0 = -1
          do k = chilev , 1 , -1
            pmpi = cht42hyam(k)*p0+pchem_3(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == chilev) then
            pmpj = cht42hyam(chilev-1)*p0+pchem_3(j,i)*cht42hybm(chilev-1)
            pmpi = cht42hyam(chilev  )*p0+pchem_3(j,i)*cht42hybm(chilev  )
            do is = 1 , nchsp
              chv4_1(j,i,l,is) = chv3(j,i,chilev,is) + &
                 (chv3(j,i,chilev-1,is) - chv3(j,i,chilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+pchem_3(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+pchem_3(j,i)*cht42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , nchsp
              chv4_1(j,i,l,is) = chv3(j,i,k0+1,is)*wt1 + chv3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')
    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_avg_1999-2009_',im2,'.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open file chemical')
    write(stdout, *) trim(chfilename)
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    xps = xps*0.01
    do is = 1 , nchsp
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_inst',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//trim(chspec(is)))
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//trim(chspec(is)))
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev)
    end do
    call bilinx2(pchem_3,xps,xlon,xlat,cht42lon,cht42lat, &
                 chilon,chjlat,iy,jx)
    do i = 1 , iy
      do j = 1 , jx
        do l = 1 , kz
          prcm=((pchem_3(j,i)*0.1-r4pt)*sigmah(l)+r4pt)*10.
          k0 = -1
          do k = chilev , 1 , -1
            pmpi = cht42hyam(k)*p0+pchem_3(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == chilev) then
            pmpj = cht42hyam(chilev-1)*p0+pchem_3(j,i)*cht42hybm(chilev-1)
            pmpi = cht42hyam(chilev  )*p0+pchem_3(j,i)*cht42hybm(chilev  )
            do is = 1 , nchsp
              chv4_2(j,i,l,is) = chv3(j,i,chilev,is) + &
                 (chv3(j,i,chilev-1,is) - chv3(j,i,chilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+pchem_3(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+pchem_3(j,i)*cht42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , nchsp
              chv4_2(j,i,l,is) = chv3(j,i,k0+1,is)*wt1 + chv3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')
  end subroutine read2m

  subroutine close_ch_icbc_clim
    implicit none
    return
  end subroutine close_ch_icbc_clim

  integer(ik4) function inextmon(im)
    implicit none
    integer(ik4) , intent(in) :: im
    inextmon = im+1
    if ( inextmon == 13 ) inextmon = 1
  end function inextmon

  integer(ik4) function iprevmon(im)
    implicit none
    integer(ik4) , intent(in) :: im
    iprevmon = im-1
    if ( iprevmon == 0 ) iprevmon = 12
  end function iprevmon

end module mod_ch_icbc_clim
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
