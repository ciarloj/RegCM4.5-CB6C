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

module mod_che_common

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_species
  use mod_che_indices
  use mod_cbmz_global , only : xr , xrin , xrout , c
  use mod_cbmz_parameters , only : nfix
  use mod_cb6_global , only : xrb , xrinb , xroutb , c_cb6 , &
                  chemtmpn , chemtmpa , chemtmpt , &
                  xrin_nu0 , xrin_ac0 , xrin_nu3 , xrin_ac3
  use mod_cb6_parameters , only : nfix_cb6

  implicit none

  public

  integer(ik4) , parameter :: sbin = 2
  integer(ik4) , parameter :: maxntr = 50

  ! tracer indices :
  type tracer
    integer(ik4) , pointer , dimension(:) :: index
    integer(ik4) , pointer , dimension(:) :: indcbmz
    integer(ik4) , pointer , dimension(:) :: indcb6
    integer(ik4) , pointer , dimension(:) :: indchbdy
    real(rk8)    , pointer , dimension(:) :: mw
    integer(ik4) , pointer , dimension(:) :: soa_category
  end type tracer
  type(tracer) trac

  ! tracer variables

  real(rk8) , pointer , dimension(:,:,:,:) :: chi 
  real(rk8) , pointer , dimension(:,:,:,:) :: chic , chiten , chiten0 , chemten
  real(rk8) , pointer , dimension(:,:,:) :: chemsrc, tmpsrc
  real(rk8) , pointer , dimension(:,:,:,:) :: chia , chib 
  real(rk8) , pointer , dimension(:,:,:) :: dtrace , wdwout , &
                                           wdrout , wxaq , wxsg , ddv_out
  real(rk8) , pointer , dimension(:,:,:) :: drydepv

  real(rk8) , pointer , dimension(:,:,:,:) :: chemall , jphoto
  real(rk8) , pointer , dimension(:,:,:,:) :: dgn , dga , dnsn , dnsa
  real(rk8) , pointer , dimension(:,:,:) :: ch_nu0 , ch_ac0 , ch_nu3 , ch_ac3
  integer(ik4) , pointer , dimension(:,:) :: kcumtop , kcumbot , cveg2d

  real(rk8) , pointer , dimension(:,:)   :: chtrsize
  real(rk8) , pointer , dimension(:)     :: chtrsol

  real(rk8), pointer , dimension(:,:,:)  :: chifxuw

  integer(ik4) , pointer , dimension(:) :: isslt , icarb , idust
  integer(ik4) , parameter :: nphoto = 56

  real(rk8) , pointer , dimension(:,:,:) :: convcldfra , cemtrac , remdrd

  ! diagnostic
  real(rk8) , pointer , dimension(:,:,:,:) :: washout , rainout , &
                                             rxsaq1 , rxsaq2 , rxsg
  real(rk8) , pointer , dimension(:,:,:,:) :: chemdiag , cadvhdiag , &
          cadvvdiag , cdifhdiag , cconvdiag , cbdydiag , ctbldiag ,  &
          cseddpdiag , cemisdiag


!*****************************************************************************
! INTERFACE VARIABLES  for chemistry / regcm
!   the pointer targets are defined in mod_che_interface
!*****************************************************************************

  real(rk8) , pointer , dimension(:,:,:,:) ::chib3d 
  real(rk8) , pointer , dimension(:,:,:,:) :: cqxb3d
  real(rk8) , pointer , dimension(:,:,:) :: bndt0 , bndt1 , bndq0 , bndq1
  real(rk8) , pointer , dimension(:,:) :: bndp0 , bndp1
  real(rk8) , pointer , dimension(:,:) :: sp0 , sp1
  real(rk8) , pointer , dimension(:,:,:) :: ctb3d , cubx3d , cvbx3d ,  &
         crhob3d , cpb3d , cpf3d , cfcc , cza , czq , cdzq , ccldfra , &
         crembc , cremrat ,  cconvpr , crhb3d , cdrydepflx , cwetdepflx , tvirt
  real(rk8) , pointer , dimension(:,:) :: cpsb , ctg , ctga , clndcat , cht , &
         cssw2da , cvegfrac , cxlai2d , csol2d , csdeltk2d , csdelqk2d , &
         cuvdrag , csfracv2d , csfracb2d , csfracs2d , cxlat , crainc ,  &
         cps2d , cps0 , cptrop
  real(rk8) , pointer , dimension(:,:) :: psbb0 , psbb1
  real(rk8) , pointer , dimension(:,:) :: czen
  real(rk8) , pointer , dimension(:,:,:,:) :: ctaucld
#if defined CLM45
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rk8) , pointer , dimension(:,:,:) :: cvoc_em_clm
  real(rk8) , pointer , dimension(:,:,:) :: cdustflx_clm
#endif
#if (defined CLM && defined VOC)
  ! Tracer mask that uses MEGAN indices
  integer(ik4) , pointer , dimension(:) :: bvoc_trmask
  real(rk8) , pointer , dimension(:,:) :: cvoc_em0
  real(rk8) , pointer , dimension(:,:) :: cvoc_em1
  real(rk8) , pointer , dimension(:,:) :: cvoc_em2
#endif
#if defined CLM
  real(rk8) , pointer , dimension(:,:,:) :: cdep_vels
#endif

  contains

  subroutine allocate_mod_che_common(isladvec)
    implicit none
    integer(ik4) , intent(in) :: isladvec

    if ( ichem == 1 ) then
      call getmem1d(trac%index,1,ntr,'mod_che_common:trac%index')
      if ( chemsimtype(1:4) == 'CB6C' ) then
        call getmem1d(trac%indcb6,1,ntr,'mod_che_common:trac%indcb6')
        call getmem1d(trac%soa_category,1,ntr,'mod_che_common:trac%soa_category')
      else
        call getmem1d(trac%indcbmz,1,ntr,'mod_che_common:trac%indcbmz')
      end if
      call getmem1d(trac%mw,1,ntr,'mod_che_common:trac%mw')
      call getmem1d(trac%indchbdy,1,ntr,'mod_che_common:trac%indchbdy')

      call getmem4d(chia,jce1-ma%jbl2,jce2+ma%jbr2, &
                    ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chia')
      if ( isladvec == 1 ) then
        call getmem4d(chib,jce1-ma%jbl4,jce2+ma%jbr4, &
                      ice1-ma%ibb4,ice2+ma%ibt4,1,kz,1,ntr,'che_common:chib')
        call getmem4d(chi,jce1-ma%jbl4,jce2+ma%jbr4, &
                      ice1-ma%ibb4,ice2+ma%ibt4,1,kz,1,ntr,'che_common:chi')
      else
        call getmem4d(chib,jce1-ma%jbl2,jce2+ma%jbr2, &
                      ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,ntr,'che_common:chib')
        call getmem4d(chi,jce1-ma%jbl1,jce2+ma%jbr1, &
                      ice1-ma%ibb1,ice2+ma%ibt1,1,kz,1,ntr,'che_common:chi')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(tvirt,jce1,jce2,ice1,ice2,1,kz,'che_common:tvirt')
        call getmem2d(sp0,jce1,jce2,ice1,ice2,'che_common:sp0')
        call getmem2d(sp1,jce1,jce2,ice1,ice2,'che_common:sp1')
      end if
      call getmem4d(chic,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chic')
      call getmem4d(chiten,jce1,jce2,ice1,ice2,1,kz,1,ntr,'che_common:chiten')
      call getmem4d(chemten,jce1,jce2, &
                    ice1,ice2,1,kz,1,ntr,'che_common:chemten')
      call getmem3d(chemsrc,jce1,jce2,ice1,ice2, &
                    1,ntr,'mod_che_common:chemsrc')

      call getmem3d(tmpsrc,jce1,jce2,ice1,ice2, &
                    1,ntr,'mod_che_common:tmpsrc')

      call getmem3d(chifxuw,jci1,jci2,ici1,ici2, &
                    1,ntr,'mod_che_common:chifxuw')

      call getmem3d(convcldfra,jci1,jci2,ici1,ici2, &
                    1,kz,'mod_che_common:convcldfra')

      call getmem4d(rxsg,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsg')
      call getmem4d(rxsaq1,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsaq1')
      call getmem4d(rxsaq2,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rxsaq2')
      call getmem4d(rainout,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:rainout')
      call getmem4d(washout,jce1,jce2,ice1,ice2,1,kz,1,ntr, &
                    'che_common:washout')
      call getmem3d(remdrd,jce1,jce2,ice1,ice2,1,ntr,'che_common:remdrd')

      call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
      call getmem1d(idust,1,nbin,'mod_che_common:idust')
      call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
      call getmem1d(icarb,1,7,'mod_che_common:icarb')
      call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')

      if ( igaschem == 1 .and. ichsolver > 0 ) then
        call getmem4d(jphoto,jci1,jci2,ici1,ici2,1,kz, &
                      1,nphoto,'che_common:jphoto')
        if ( chemsimtype(1:4) == 'CB6C' ) then
          if ( isorgam == 1 ) then
            call getmem4d(chemall,jci1,jci2,ici1,ici2, &
                          1,kz,1,3*76,'mod_che_common:chemall')
            call getmem1d(xrb,1,3*76,'che_common:xrb')
            call getmem1d(xrinb,1,3*76,'che_common:xrinb')
            call getmem1d(xroutb,1,3*76,'che_common:xroutb')
            call getmem1d(c_cb6,1,totch+nfix_cb6,'che_common:c_cb6')
            call getmem1d(chemtmpn,1,76,'che_common:chemtmpn')
            call getmem1d(chemtmpa,1,76,'che_common:chemtmpa')
            call getmem1d(chemtmpt,1,12,'che_common:chemtmpt')
            call getmem4d(dnsn,jci1,jci2,ici1,ici2,1,kz,1,3*76,'mod_che_common:dnsn')
            call getmem4d(dnsa,jci1,jci2,ici1,ici2,1,kz,1,3*76,'mod_che_common:dnsa')
            call getmem4d(dgn,jci1,jci2,ici1,ici2,1,kz,1,3*76,'mod_che_common:dgn')
            call getmem4d(dga,jci1,jci2,ici1,ici2,1,kz,1,3*76,'mod_che_common:dga')
            call getmem3d(ch_nu0,jci1,jci2,ici1,ici2,1,kz,'che_common:ch_nu0')
            call getmem3d(ch_ac0,jci1,jci2,ici1,ici2,1,kz,'che_common:ch_ac0')
            call getmem3d(ch_nu3,jci1,jci2,ici1,ici2,1,kz,'che_common:ch_nu3')
            call getmem3d(ch_ac3,jci1,jci2,ici1,ici2,1,kz,'che_common:ch_ac3')
            call getmem1d(xrin_nu0,1,1,'che_common:xrin_nu0')
            call getmem1d(xrin_ac0,1,1,'che_common:xrin_ac0')
            call getmem1d(xrin_nu3,1,1,'che_common:xrin_nu3')
            call getmem1d(xrin_ac3,1,1,'che_common:xrin_ac3')
          else
            call getmem4d(chemall,jci1,jci2,ici1,ici2, &
                          1,kz,1,totch,'mod_che_common:chemall')
            call getmem1d(xrb,1,totch,'che_common:xrb')
            call getmem1d(xrinb,1,totch,'che_common:xrinb')
            call getmem1d(xroutb,1,totch,'che_common:xroutb')
            call getmem1d(c_cb6,1,totch+nfix_cb6,'che_common:c_cb6')
          end if
        else
          call getmem4d(chemall,jci1,jci2,ici1,ici2, &
                        1,kz,1,totsp,'mod_che_common:chemall')
          call getmem1d(xr,1,totsp,'che_common:xr')
          call getmem1d(xrin,1,totsp,'che_common:xrin')
          call getmem1d(xrout,1,totsp,'che_common:xrout')
          call getmem1d(c,1,totsp+nfix,'che_common:c')
        end if
      end if

      call getmem3d(dtrace,jce1,jce2,ice1,ice2,1,ntr,'che_common:dtrace')
      call getmem3d(wdrout,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdrout')
      call getmem3d(wdwout,jce1,jce2,ice1,ice2,1,ntr,'che_common:wdwout')

      call getmem3d(wxsg,jce1,jce2,ice1,ice2,1,ntr,'che_common:wxsg')
      call getmem3d(wxaq,jce1,jce2,ice1,ice2,1,ntr,'che_common:wxaq')
      call getmem3d(cemtrac,jce1,jce2,ice1,ice2,1,ntr,'che_common:cemtrac')
      call getmem3d(drydepv,jce1,jce2,ice1,ice2,1,ntr,'che_common:drydepv')
      call getmem3d(ddv_out,jce1,jce2,ice1,ice2,1,ntr,'che_common:ddv_out')

      if ( ichdiag > 0 ) then
        call getmem4d(chiten0,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chiten0')
        call getmem4d(chemdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:chemdiag')
        call getmem4d(cadvhdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cadvhdiag')
        call getmem4d(cadvvdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cadvvdiag')
        call getmem4d(cdifhdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cdifhdiag')
        call getmem4d(cconvdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cconvdiag')
        call getmem4d(ctbldiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:ctbldiag')
        call getmem4d(cbdydiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cbdydiag')
        call getmem4d(cseddpdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cseddpdiag')
        call getmem4d(cemisdiag,jce1,jce2, &
                      ice1,ice2,1,kz,1,ntr,'che_common:cemisdiag')
      end if
#if defined CLM45 || (defined CLM && defined VOC)
      call getmem1d(bvoc_trmask,1,ntr,'mod_che_common:bvoc_trmask')
#endif
    end if
  end subroutine allocate_mod_che_common
!
  subroutine chem_config
    implicit none
    integer(ik4) :: i
    ! Define here the possible types of simulation and fix the dimension
    ! of relevant tracer dimension and parameters
    igaschem = 0
    iaerosol = 0
    iisoropia = 0
    if ( chemsimtype(1:4) == 'DUST' ) then
      nbin = 4
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST simulation'
    else if ( chemsimtype(1:4) == 'DU12' ) then
      nbin = 12
      ntr = nbin
      allocate(chtrname(nbin))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04',&
                               'DUST05','DUST06','DUST07','DUST08',&
                               'DUST09','DUST10','DUST11','DUST12' /)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUST 12 bins simulation'
    else if ( chemsimtype(1:4) == 'SSLT' ) then
      ntr = sbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'SSLT01','SSLT02'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'SSLT simulation'
    else if ( chemsimtype(1:4) == 'DUSS' ) then
      nbin = 4
      ntr = sbin + nbin
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'DUST01','DUST02','DUST03','DUST04', &
                               'SSLT01','SSLT02'/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'DUSS simulation'
    else if ( chemsimtype(1:4) == 'CARB' ) then
      ntr = 4
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB '/)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'CARB simulation'
    else if ( chemsimtype(1:4) == 'SULF' ) then
      ntr = 2
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'SO2   ','SO4   '/)
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SULF simulation'
    else if ( chemsimtype(1:4) == 'SUCA' ) then
      ntr = 6
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                               'SO2   ','SO4   '/)
      iaerosol = 1
      ioxclim  = 1
      if ( myid == italk ) write(stdout,*) 'SUCA simulation'
    else if ( chemsimtype(1:4) == 'AERO' ) then
      nbin = 4
      ntr = 12
      allocate(chtrname(ntr))
      iaerosol = 1
      ioxclim  = 1
      chtrname(1:ntr)(1:6) = (/'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                               'SO2   ','SO4   ','DUST01','DUST02', &
                               'DUST03','DUST04','SSLT01','SSLT02' /)
      if ( myid == italk ) write(stdout,*) 'AERO simulation'
    else if ( chemsimtype(1:4) == 'DCCB' ) then
      nbin = 4
      ntr = 50
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'NO    ','NO2   ','N2O5  ','HNO2  ',&
                               'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                               'CO    ','SO2   ','DMS   ','H2SO4 ',&
                               'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                               'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                               'OLET  ','OLEI  ','TOL   ','XYL   ',&
                               'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                               'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                               'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                               'CRES  ','NH3   ','DUST01','DUST02',&
                               'DUST03','DUST04','BC_HL ','BC_HB ',&
                               'OC_HL ','OC_HB ','SSLT01','SSLT02',&
                               'ANO3  ','ANH4  ' /)
      iaerosol = 1
      igaschem = 1
      iisoropia = 1
      if ( myid == italk ) write(stdout,*) 'DCCB simulation'
    else if ( chemsimtype(1:4) == 'CBMZ' ) then
      ! This does not include any aerosol(NH3) or monoterpens(APIN, LIMO)
      ntr = 37
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'NO    ','NO2   ','N2O5  ','HNO2  ',&
                               'HNO3  ','HNO4  ','O3    ','H2O2  ',&
                               'CO    ','SO2   ','DMS   ','H2SO4 ',&
                               'CH4   ','C2H6  ','PAR   ','CH3OH ',&
                               'HCHO  ','ALD2  ','AONE  ','ETH   ',&
                               'OLET  ','OLEI  ','TOL   ','XYL   ',&
                               'ISOP  ','ONIT  ','PAN   ','HCOOH ',&
                               'RCOOH ','CH3OOH','ETHOOH','ROOH  ',&
                               'MGLY  ','ISOPRD','ISOPN ','OPEN  ',&
                               'CRES  '/)
      igaschem = 1
      if ( myid == italk ) write(stdout,*) 'CBMZ simulation'
    else if ( chemsimtype(1:4) == 'CB6C' ) then
      if ( isorgam == 1 ) then
        ngstr  = 76
        nsoatr = 12              !aggregated aerosols
        ntr    = 3*ngstr+nsoatr  !(        1:ngstr  ) = gas species
                                 !(  ngstr+1:2*ngstr) = nuc aerosols
                                 !(2*ngstr+1:3*ngstr) = acc aerosols
                                 !(3+ngstr+1:ntr    ) = aggregates
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'OZN   ','CMON  ','HPOX  ','NMOX  ',&
                                 'NDOX  ','DNPO  ','HONO  ','NTRC  ',&
                                 'CDOX  ','NTR1  ','NTR2  ','INTR  ',&
                                 'SDIO  ','SULF  ','MEOH  ','MEPX  ',&
                                 'FACD  ','FORM  ','ETHE  ','ETHA  ',&
                                 'ETHY  ','ETOH  ','ACET  ','AACD  ',&
                                 'ALKA  ','PRPA  ','PACD  ','RPOX  ',&
                                 'BENZ  ','TOLN  ','XYLN  ','XOPN  ',&
                                 'EPOX  ','CAT1  ','CRON  ','CRSL  ',&
                                 'TERP  ','GLY   ','MEGY  ','GLYD  ',&
                                 'AALD  ','ALDX  ','KET   ','HPLD  ',&
                                 'PACN  ','OPAN  ','ROPN  ','PNA   ',&
                                 'ISPR  ','ISPD  ','ISPX  ','IOLE  ',&
                                 'OLE   ','NTOX  ','PANX  ','OSNG  ',&
                                 'MTHA  ','HCO3  ','CRER  ','BZO2  ',&
                                 'ROR   ','TOLR  ','XO2R  ','EPX2  ',&
                                 'XYLR  ','OPO3  ','XO2N  ','ISO2  ',&
                                 'XO2H  ','MEO2  ','ROO   ','HOX   ',&
                                 'POX   ','CXO3  ','ACOO  ','O     ',&
                                 'OZNnu ','CMONnu','HPOXnu','NMOXnu',&
                                 'NDOXnu','DNPOnu','HONOnu','NTRCnu',&
                                 'CDOXnu','NTR1nu','NTR2nu','INTRnu',&
                                 'SDIOnu','SULFnu','MEOHnu','MEPXnu',&
                                 'FACDnu','FORMnu','ETHEnu','ETHAnu',&
                                 'ETHYnu','ETOHnu','ACETnu','AACDnu',&
                                 'ALKAnu','PRPAnu','PACDnu','RPOXnu',&
                                 'BENZnu','TOLNnu','XYLNnu','XOPNnu',&
                                 'EPOXnu','CAT1nu','CRONnu','CRSLnu',&
                                 'TERPnu','GLYnu ','MEGYnu','GLYDnu',&
                                 'AALDnu','ALDXnu','KETnu ','HPLDnu',&
                                 'PACNnu','OPANnu','ROPNnu','PNAnu ',&
                                 'ISPRnu','ISPDnu','ISPXnu','IOLEnu',&
                                 'OLEnu ','NTOXnu','PANXnu','OSNGnu',&
                                 'MTHAnu','HCO3nu','CRERnu','BZO2nu',&
                                 'RORnu ','TOLRnu','XO2Rnu','EPX2nu',&
                                 'XYLRnu','OPO3nu','XO2Nnu','ISO2nu',&
                                 'XO2Hnu','MEO2nu','ROOnu ','HOXnu ',&
                                 'POXnu ','CXO3nu','ACOOnu','Onu   ',&
                                 'OZNac ','CMONac','HPOXac','NMOXac',&
                                 'NDOXac','DNPOac','HONOac','NTRCac',&
                                 'CDOXac','NTR1ac','NTR2ac','INTRac',&
                                 'SDIOac','SULFac','MEOHac','MEPXac',&
                                 'FACDac','FORMac','ETHEac','ETHAac',&
                                 'ETHYac','ETOHac','ACETac','AACDac',&
                                 'ALKAac','PRPAac','PACDac','RPOXac',&
                                 'BENZac','TOLNac','XYLNac','XOPNac',&
                                 'EPOXac','CAT1ac','CRONac','CRSLac',&
                                 'TERPac','GLYac ','MEGYac','GLYDac',&
                                 'AALDac','ALDXac','KETac ','HPLDac',&
                                 'PACNac','OPANac','ROPNac','PNAac ',&
                                 'ISPRac','ISPDac','ISPXac','IOLEac',&
                                 'OLEac ','NTOXac','PANXac','OSNGac',&
                                 'MTHAac','HCO3ac','CRERac','BZO2ac',&
                                 'RORac ','TOLRac','XO2Rac','EPX2ac',&
                                 'XYLRac','OPO3ac','XO2Nac','ISO2ac',&
                                 'XO2Hac','MEO2ac','ROOac ','HOXac ',&
                                 'POXac ','CXO3ac','ACOOac','Oac   ',&
                                 'CNYLnu','STHCnu','OLEFnu','AROMnu',&
                                 'NTROnu','MFNCnu','CNYLac','STHCac',&
                                 'OLEFac','AROMac','NTROac','MFNCac'/)
        iaerosol = 1
        igaschem = 1
        if ( myid == italk ) write(stdout,*) 'SORGAM simulation'
      else
        ntr = 56  ! = totch ! if changed check with max_incb6_tracers
       !ntr = 76  ! = totch ! if changed check with max_incb6_tracers
        allocate(chtrname(ntr))
        chtrname(1:ntr)(1:6) = (/'OZN   ','CMON  ','HPOX  ','NMOX  ',&
                                 'NDOX  ','DNPO  ','HONO  ','NTRC  ',&
                                 'CDOX  ','NTR1  ','NTR2  ','INTR  ',&
                                 'SDIO  ','SULF  ','MEOH  ','MEPX  ',&
                                 'FACD  ','FORM  ','ETHE  ','ETHA  ',&
                                 'ETHY  ','ETOH  ','ACET  ','AACD  ',&
                                 'ALKA  ','PRPA  ','PACD  ','RPOX  ',&
                                 'BENZ  ','TOLN  ','XYLN  ','XOPN  ',&
                                 'EPOX  ','CAT1  ','CRON  ','CRSL  ',&
                                 'TERP  ','GLY   ','MEGY  ','GLYD  ',&
                                 'AALD  ','ALDX  ','KET   ','HPLD  ',&
                                 'PACN  ','OPAN  ','ROPN  ','PNA   ',&
                                 'ISPR  ','ISPD  ','ISPX  ','IOLE  ',&
                                 'OLE   ','NTOX  ','PANX  ','MTHA  '/)
       !chtrname(1:ntr)(1:6) = (/'OZN   ','CMON  ','HPOX  ','NMOX  ',&
       !                         'NDOX  ','DNPO  ','HONO  ','NTRC  ',&
       !                         'CDOX  ','NTR1  ','NTR2  ','INTR  ',&
       !                         'SDIO  ','SULF  ','MEOH  ','MEPX  ',&
       !                         'FACD  ','FORM  ','ETHE  ','ETHA  ',&
       !                         'ETHY  ','ETOH  ','ACET  ','AACD  ',&
       !                         'ALKA  ','PRPA  ','PACD  ','RPOX  ',&
       !                         'BENZ  ','TOLN  ','XYLN  ','XOPN  ',&
       !                         'EPOX  ','CAT1  ','CRON  ','CRSL  ',&
       !                         'TERP  ','GLY   ','MEGY  ','GLYD  ',&
       !                         'AALD  ','ALDX  ','KET   ','HPLD  ',&
       !                         'PACN  ','OPAN  ','ROPN  ','PNA   ',&
       !                         'ISPR  ','ISPD  ','ISPX  ','IOLE  ',&
       !                         'OLE   ','NTOX  ','PANX  ','OSNG  ',&
       !                         'MTHA  ','HCO3  ','CRER  ','BZO2  ',&
       !                         'ROR   ','TOLR  ','XO2R  ','EPX2  ',&
       !                         'XYLR  ','OPO3  ','XO2N  ','ISO2  ',&
       !                         'XO2H  ','MEO2  ','ROO   ','HOX   ',&
       !                         'POX   ','CXO3  ','ACOO  ','O     '/)
        igaschem = 1
        if ( myid == italk ) write(stdout,*) 'CB6C simulation'
      end if
    else if ( chemsimtype(1:6) == 'POLLEN' ) then
      ntr = 1
      allocate(chtrname(ntr))
      chtrname(1:ntr)(1:6) = (/'POLLEN' /)
      iaerosol = 1
      if ( myid == italk ) write(stdout,*) 'POLLEN simulation'
    else
      if ( myid == italk ) then
        write (stderr,*) 'Not a valid chemtype simulation : STOP !'
        write (stderr,*) 'Valid simulations are : ' , &
           'DUST DU12 SSLT DUSS CARB SULF SUCA AERO CBMZ CB6C DCCB POLLEN'
      end if
      call fatal(__FILE__,__LINE__,'INVALID CHEM CONFIGURATION')
    end if

  end subroutine chem_config

end module mod_che_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
