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

module mod_che_emission
  !
  ! Chemical and aerosol surface emission
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_service
  use mod_dynparam
  use mod_che_common
  use mod_che_param
  use mod_che_mppio
  use mod_che_dust
  use mod_che_ncio
  use mod_che_indices
  use netcdf
  !
  implicit none

  private
  !
  public :: chem_emission , emis_tend
  !
  contains
  !
  ! SURFACE EMIOSSION flux
  !
  subroutine chem_emission(lyear,lmonth,lday,lhour)
    implicit none
    integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
    integer(ik4) , save :: curry = -1
    integer(ik4) , save :: currm = -1
    integer(ik4) , save :: currd = -1
    integer(ik4) , save :: ifreq = -1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_emission'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! read the aerosol emission files
    !
    if ( chemsimtype(1:2) == 'DU' .or. &
         chemsimtype(1:4) == 'SSLT' .or. &
         chemsimtype(1:4) == 'DUSS' ) return
    if ( ifreq == ifrqmon ) then
      if ( curry == lyear .and. currm == lmonth ) then
        if ( myid == italk ) then
          write(stdout,*) &
            'EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,' ready', &
            ' from ',curry*1000000+currm*10000+100*currd
        end if
        return
      end if
    else if ( ifreq == ifrqday ) then
      if ( curry == lyear .and. currm == lmonth .and. currd == lday ) return
    end if
    curry = lyear
    currm = lmonth
    currd = lday
    if ( myid == italk ) then
      write(stdout,*)'READ CHEM EMISSION for ',lyear*1000000+lmonth*10000+lday
    end if
    ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
    ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES
    call read_emission(ifreq,lyear,lmonth,lday,lhour,chemsrc)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_emission
  !
  !
  ! Calculation of emission tendency
  !
  !FAB: no dirunal evol for now:
  ! subroutine emis_tend(ktau,j,lmonth,xlat,coszrs,declin,dsigma)
  subroutine emis_tend(ktau,j,lmonth,declin)
    implicit none

    integer(ik4) , intent(in) :: j , lmonth
    integer(ik8) , intent(in) :: ktau

    real(rk8) , intent(in) ::declin
    integer(ik4)  :: i , itr
    real(rk8) :: daylen , fact , maxelev , amp,dayhr
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'emis_tend'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! calculate the tendency linked to emissions from emission fluxes
    ! In the future split these calculations in corresponding module  ??
    ! Modify chemsrc for species that need a dirunal cycle
    if ( iisop > 0 .or. iispr > 0  ) then
#ifdef CLM45
      !overwrite chemsrc for biogenic in case of CLM/BVOC option
      ! test if CLM BVOC is activated and overwrite chemsrc.
      ! Below included in order to include CLM-MEGAN biogenic emission
      ! into the gas phase chemistry scheme
      ! bvoc_trmask is used to identify which tracer is also a biogenic
      ! emission species contained in MEGAN.  If it is included then
      ! then include MEGAN emissions into the source value
      ! NOTE:  ibvoc=1 means used MEGAN emissions.  ibvoc is forced to
      ! zero when using BATS
      if ( chemsimtype(1:4) ==  'CB6C') then
        if ( bvoc_trmask(iispr) /= 0 ) then
          do i = ici1, ici2
            chemsrc(j,i,iispr) = cvoc_em(j,i)
          end do
        end if
      else
        if ( bvoc_trmask(iisop) /= 0 ) then
          do i = ici1, ici2
            chemsrc(j,i,iisop) = cvoc_em_clm(j,i,iisop)
          end do
        end if
      end if
#else
#if (defined VOC && defined CLM)
      !overwrite chemsrc for biogenic in case of CLM/BVOC option
      ! test if CLM BVOC is activated and overwrite chemsrc.
      ! Below included in order to include CLM-MEGAN biogenic
      ! emission
      ! NOTE:  ibvoc=1 means used MEGAN emissions.
      ! ibvoc is forced to
      ! zero when using BATS
      if ( chemsimtype(1:4) ==  'CB6C') then
        if ( bvoc_trmask(iispr) /= 0 ) then
          do i = ici1, ici2
            chemsrc(j,i,iispr) = cvoc_em0(j,i)
            if ( ktau == 0 ) cvoc_em0(j,i) = d_zero
            chemsrc(j,i,iispr) = cvoc_em0(j,i)
          end do
        end if
      else
        if ( bvoc_trmask(iisop) /= 0 ) then
          do i = ici1, ici2
            chemsrc(j,i,iisop) = cvoc_em0(j,i)
            if ( ktau == 0 ) cvoc_em0(j,i) = d_zero
            chemsrc(j,i,iisop) = cvoc_em0(j,i)
          end do
        end if
      end if
#else
      ! Modify chemsrc for species that need a dirunal cycle
      if ( chemsimtype(1:4) ==  'CB6C') then
        do i = ici1 , ici2
          dayhr=-tan(declin)*tan(cxlat(j,i)*degrad)
          if(dayhr .lt. -1 .or. dayhr .gt. 1) then
          tmpsrc(j,i,iispr)  = chemsrc(j,i,iispr)
          chemsrc(j,i,iispr) = 0.0
          else
          daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
          daylen = d_two*acos(dayhr)*raddeg
          daylen = daylen*24.0D0/360.0D0
          ! Maximum sun elevation
          maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
          fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
          amp = 12.0D0*mathpi/daylen
          tmpsrc(j,i,iispr)  = chemsrc(j,i,iispr)
          chemsrc(j,i,iispr) = (amp)*chemsrc(j,i,iispr) * &
                              sin(mathpi*fact)*egrav/(dsigma(kz)*1.0D3)
          end if
        end do
      else
        do i = ici1 , ici2
          dayhr=-tan(declin)*tan(cxlat(j,i)*degrad)
          if(dayhr .lt. -1 .or. dayhr .gt. 1) then
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = 0.0
          else
          daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
          daylen = d_two*acos(dayhr)*raddeg
          daylen = daylen*24.0D0/360.0D0
          ! Maximum sun elevation
          maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
          fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
          amp = 12.0D0*mathpi/daylen
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = (amp)*chemsrc(j,i,iisop) * &
                              sin(mathpi*fact)*egrav/(dsigma(kz)*1.0D3)
          end if
        end do
      end if

#endif
#endif
    end if

    if ( iapin > 0 ) then
#if (defined VOC && defined CLM)
      if ( bvoc_trmask(iapin) /= 0 ) then
        do i = ici1, ici2
          if ( ktau == 0 ) cvoc_em1(j,i) = d_zero
          chemsrc(j,i,iapin) = cvoc_em1(j,i)
        end do
      end if
#endif
    end if

    if ( ilimo > 0 ) then
#if (defined VOC && defined CLM)
      if ( bvoc_trmask(ilimo) /= 0 ) then
        do i = ici1, ici2
          if ( ktau == 0 ) cvoc_em2(j,i) = d_zero
          chemsrc(j,i,ilimo) = cvoc_em2(j,i)
        end do
      end if
#endif
    end if
    !
    ! add the source term to tracer tendency
    !
    if ( ichdrdepo /= 2 ) then
      do itr = 1 , ntr
        do i = ici1 , ici2
          if ( chtrname(itr)(1:2) == 'DU' .or. &
               chtrname(itr)(1:4) == 'SSLT' .or. &
               chtrname(itr)(1:6) == 'POLLEN' ) cycle
          chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
              chemsrc(j,i,itr)*egrav/(dsigma(kz)*1.0D3)
          ! diagnostic for source, cumul
          cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
          if ( ichdiag == 1 ) then
            cemisdiag(j,i,kz,itr) = cemisdiag(j,i,kz,itr) + &
                chemsrc(j,i,itr)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
          end if
        end do
      end do
    else if ( ichdrdepo == 2 ) then
      do itr = 1 , ntr
        do i = ici1 , ici2
          if ( chtrname(itr)(1:2).ne.'DU' .or. &
               chtrname(itr)(1:4).ne.'SSLT' .or. &
               chtrname(itr)(1:6).ne.'POLLEN' ) then
            ! if PBL scheme is not UW then calculate emission tendency
            if ( ibltyp /= 2 ) then
              chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                 chemsrc(j,i,itr)*egrav/(dsigma(kz)*1.0D3)
            end if
            ! otherwise emission is injected in the PBL scheme ( together
            ! with dry deposition) for tend calculation
            chifxuw(j,i,itr) = chifxuw(j,i,itr) +  chemsrc(j,i,itr)
            ! diagnostic for source, cumul
            cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
            if ( ichdiag == 1 ) then
              ! in this case we calculate emission tendency diagnostic, but
              ! this term will also be included in BL tendency diagnostic
              ! if UW scheme is used.
              if ( ibltyp /= 2 ) then
                cemisdiag(j,i,kz,itr) = cemisdiag(j,i,kz,itr) + &
                   chemsrc(j,i,itr)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
              end if
            end if
          end if
        end do
      end do
    end if
    ! put back isop source to its nominal value
    if ( chemsimtype(1:4) ==  'CB6C') then
      if ( iispr > 0 ) then
        chemsrc(j,:,iispr) = tmpsrc(j,:,iispr)
      end if
    else
      if ( iisop > 0 ) then
        chemsrc(j,:,iisop) = tmpsrc(j,:,iisop)
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine emis_tend
!
end module mod_che_emission
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
