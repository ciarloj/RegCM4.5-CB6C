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

module mod_che_chemistry

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams , only : iqv , calday
  use mod_che_common
  use mod_che_indices
  use mod_che_species
  use mod_cbmz_Global
  use mod_cbmz_Parameters
  use mod_cbmz_main1
  use mod_cb6_Global
  use mod_cb6_Parameters
  use mod_cb6_main1
  use mod_che_molwg
  use mod_che_soa_cb6
  use mod_che_aerosols_sorgam

  implicit none

  private

  real(rk8) :: dtchsolv
  integer(ik8) :: kchsolv

  integer(ik4) , parameter :: kmin = 2
  real(rk8) , parameter :: kb = 1.380658D-19
  real(rk8) , parameter :: mwa = 28.97D0

  public :: chemistry , chemistry_cb6 , chemistry_srg , dtchsolv , kchsolv 

  contains

    subroutine chemistry(j)
      implicit none
      integer(ik4) , intent(in) :: j
      real(rk8) :: cfactor , pfact
      integer(ik4) :: i , k , ic , n

      time = dtchsolv

      ! Begining of i , k loop
      ! do not solve chemistry anyway for topmost layer
      do k = kmin , kz
        do i = ici1 , ici2
          altmid   = cpb3d(j,i,k)
          ! Skip stratosphere ( ? Should we ? )
          if ( altmid < cptrop(j,i) - d_100 ) cycle
          temp     = ctb3d(j,i,k)
          zenith   = dacos(czen(j,i))*raddeg
          cfactor  = crhob3d(j,i,k) * 1.D-03 * navgdr
          dens     = cfactor / amd
          C_M      = altmid*10.0D0/(kb*temp)
          deptha   = d_zero
          depthb   = d_zero
          altabove = d_zero
          altbelow = d_zero

          if ( ichjphcld == 1 ) then
            deptha   = sum(ctaucld(j,i,1:k-1,8))
            altabove = sum(cdzq(j,i,1:k-1)*ctaucld(j,i,1:k-1,8))
            depthb   = sum(ctaucld(j,i,k+1:kz,8))
            altbelow = sum(cdzq(j,i,k+1:kz)*ctaucld(j,i,k+1:kz,8))
            if ( deptha > d_zero ) altabove = altabove / deptha
            if ( depthb > d_zero ) altbelow = altbelow / depthb
          end if

          ! call the chemistry solver
          xr(:) = d_zero
          xrin(:) = d_zero
          xrout(:) = d_zero
          ! 1 : initialise xrin with the concentrations from
          !     previous chemsolv step
          do ic = 1 , totsp
!           xrin(ic) = chemall(j,i,k,ic)
          end do
          ! 2 : update input concentrations for transported species only
          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
              xrin(trac%indcbmz(n)) = chib3d(j,i,k,n) * cfactor / trac%mw(n)
            end if
          end do
          ! update for water vapor
          xrin(ind_H2O)  = cqxb3d(j,i,k,iqv) * cfactor / amw

          call chemmain(calday,dtchsolv)

          ! save the concentrations of all species for next chemistry step
          do ic = 1 , totsp
            chemall(j,i,k,ic) = xrout(ic)
          end do
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            jphoto(j,i,k,ic) = c_jval(1,ic)
          end do
          !
          ! Now calculate chemical tendencies
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          pfact = cpsb(j,i) / cfactor / dtchsolv

          do n = 1 , ntr
            if ( trac%indcbmz(n) > 0 ) then
              !ashalaby
              if(xrout(trac%indcbmz(n)) .gt. 1.0e20 .or. &
                  xrin(trac%indcbmz(n))  .gt. 1.0e20) then
                write(stderr,*)xrout(trac%indcbmz(n)), & 
                                    xrin(trac%indcbmz(n)),trac%indcbmz(n)
                call fatal(__FILE__,__LINE__,'chemistry')
              else
              !ashalaby_
                chemten(j,i,k,n) = (xrout(trac%indcbmz(n)) - &
                                    xrin(trac%indcbmz(n))) * pfact * trac%mw(n)
              end if
            end if
          end do
        end do ! end i , k loop
      end do
    end subroutine chemistry

    subroutine chemistry_srg(j,cerr) ! subroutine for cb6 with sorgam
      implicit none
      real(rk8), dimension(3*76) :: drog ! size = totch
      integer(ik4) , intent(in) :: j
      integer(ik4) , intent(inout) :: cerr
      real(rk8) :: cfactor , pfact , cfacsrg
      integer(ik4) :: i , k , ic , n , iwt , numchem

      timeb = dtchsolv
      cfacsrg = 1.D+12 / ( navgdr )! * dtchsolv )
      numchem = 3*76

      ! Begining of i , k loop
      ! do not solve chemistry anyway for topmost layer
      do k = kmin , kz
        do i = ici1 , ici2
          altmidb   = cpb3d(j,i,k)
          relhum    = crhb3d(j,i,k)
          ! Skip stratosphere ( ? Should we ? )
          if ( altmidb < cptrop(j,i) - d_100 ) cycle
          tempb     = ctb3d(j,i,k)
          zenithb   = dacos(czen(j,i)) * raddeg
          cfactor   = crhob3d(j,i,k) * 1.D-03 * navgdr
          densb     = cfactor / amd
          C_Mb      = altmidb * 10.0D0 /(kb*tempb)
          depthab   = d_zero
          depthbb   = d_zero
          altaboveb = d_zero
          altbelowb = d_zero

          if ( ichjphcld == 1 ) then
            depthab   = sum(ctaucld(j,i,1:k-1,8))
            altaboveb = sum(cdzq(j,i,1:k-1) * ctaucld(j,i,1:k-1,8))
            depthbb   = sum(ctaucld(j,i,k+1:kz,8))
            altbelowb = sum(cdzq(j,i,k+1:kz) * ctaucld(j,i,k+1:kz,8))
            if ( depthab > d_zero ) altaboveb = altaboveb / depthab
            if ( depthbb > d_zero ) altbelowb = altbelowb / depthbb
          end if

          ! call the chemistry solver
          xrb(:)    = d_zero
          xrinb(:)  = d_zero
          xroutb(:) = d_zero
          drog(:)     = d_zero
          xrin_nu0(1) = ch_nu0(j,i,k)
          xrin_ac0(1) = ch_ac0(j,i,k)
          xrin_nu3(1) = ch_nu3(j,i,k)
          xrin_ac3(1) = ch_ac3(j,i,k)
          ! 1 : initialise xrin with the concentrations from
          !     previous chemsolv step
          do ic = 1 , numchem !totch
            xrinb(ic) = chemall(j,i,k,ic)
          end do
          ! 2 : update input concentrations for transported species only
          do n = 1 , numchem !ntr
            if ( trac%indcb6(n) > 0 ) then
              xrinb(trac%indcb6(n)) = chib3d(j,i,k,n) * cfactor / trac%mw(n)
            end if
          end do
          ! update for water vapor
          wtrin = cqxb3d(j,i,k,iqv) * cfactor / amw

          call chemmain_cb6(calday,dtchsolv)

          ! start process for aerosol conversion molec.cm-3 to ug.m-3.s-1
          do ic = 1 , numchem !totch ! unit conversion necessary for SORGAM
            if ( ic < ngstr+1 ) iwt = ic
            if ( ic > ngstr .and. ic < 2*ngstr+1 ) iwt = ic-ngstr
            if ( ic > 2*ngstr ) iwt = ic-(2*ngstr)
            xrinb(ic)  = xrinb(ic)  * cfacsrg * mw_cb6(iwt)
            xroutb(ic) = xroutb(ic) * cfacsrg * mw_cb6(iwt)
          end do
          wtrout = wtrout * cfacsrg * w_wtr
          do ic = 1 , ngstr
            if ( soa_cat_cb6(ic) == 99 ) then
              drog(ic) = xrinb(ic) - xroutb(ic)
            end if
          end do
          call sorgam_driver(dtchsolv,tempb,relhum,altmidb,xroutb(1:ngstr), &
                 xroutb(ngstr+1:2*ngstr),xroutb(2*ngstr+1:3*ngstr),         &
                 xrin_nu0(1),xrin_ac0(1),xrin_nu3(1),xrin_ac3(1),           &
                 xrinb(ind_SULF),drog(1:ngstr),i,j,k,ngstr,                 &
                 soa_cat_cb6(1:ngstr),pdensn,pdensa,dgnuc,dgacc)
          do ic = 1 , numchem !totch ! convert units back
            if ( ic < ngstr+1 ) iwt = ic
            if ( ic > ngstr .and. ic < 2*ngstr+1 ) iwt = ic-ngstr
            if ( ic > 2*ngstr ) iwt = ic-(2*ngstr)
            xrinb(ic)  = xrinb(ic)  / ( cfacsrg * mw_cb6(iwt) )
            xroutb(ic) = xroutb(ic) / ( cfacsrg * mw_cb6(iwt) )
            dnsn(j,i,k,ic) = pdensn
            dnsa(j,i,k,ic) = pdensa
            dgn(j,i,k,ic)  = dgnuc
            dga(j,i,k,ic)  = dgacc
          end do
          wtrout = wtrout / ( cfacsrg * w_wtr )

          ! save the concentrations of all species for next chemistry step
          do ic = 1 , numchem !totch
            chemall(j,i,k,ic) = xroutb(ic)
          end do
          ch_nu0(j,i,k) = xrin_nu0(1)
          ch_ac0(j,i,k) = xrin_ac0(1)
          ch_nu3(j,i,k) = xrin_nu3(1)
          ch_ac3(j,i,k) = xrin_ac3(1)
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            jphoto(j,i,k,ic) = c_jvalb(1,ic)
          end do

          ! Now calculate chemical tendencies
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          pfact = cpsb(j,i) / cfactor / dtchsolv

          do n = 1 , numchem !ntr
            if ( trac%indcb6(n) > 0 ) then
              !ashalaby
              if(xroutb(trac%indcb6(n)) .gt. 1.0e20 .or. &
                  xrinb(trac%indcb6(n)) .gt. 1.0e20) then
                write(stderr,*)xroutb(trac%indcb6(n)), &
                                    xrinb(trac%indcb6(n)),trac%indcb6(n)
                cerr = 1
                call fatal(__FILE__,__LINE__,'chemistry')
                !ashalaby_
                if ( cerr == 1 ) EXIT
              else
                chemten(j,i,k,n) = (xroutb(trac%indcb6(n)) - &
                                    xrinb(trac%indcb6(n))) * pfact * trac%mw(n)
              end if
            end if
          end do

          chemtmpn(:) = chemten(j,i,k,ngstr+1:2*ngstr)
          chemtmpa(:) = chemten(j,i,k,2*ngstr+1:3*ngstr)
          chemtmpt(:) = chemten(j,i,k,3*ngstr+1:ntr)
          call soa_aggregate(trac%soa_category(1:ngstr),chemtmpn,   &
                             chemtmpa,ngstr,nsoatr,chemtmpt)
          chemten(j,i,k,3*ngstr+1:ntr) = chemtmpt(:)

        if ( cerr == 1 ) EXIT
        end do ! end i loop
      if ( cerr == 1 ) EXIT
      end do ! end k loop

    end subroutine chemistry_srg

    subroutine chemistry_cb6(j,cerr) ! subroutine for cb6 only
      implicit none
      integer(ik4) , intent(in) :: j 
      integer(ik4) , intent(inout) :: cerr
      real(rk8) :: cfactor , pfact , cfacsrg
      integer(ik4) :: i , k , ic , n , iwt , numchem , a

      timeb = dtchsolv
      cfacsrg = 1.D+12 / ( navgdr )! * dtchsolv )
      numchem = totch

      ! Begining of i , k loop
      ! do not solve chemistry anyway for topmost layer
      do k = kmin , kz
        do i = ici1 , ici2
          altmidb   = cpb3d(j,i,k)
          relhum    = crhb3d(j,i,k)
          ! Skip stratosphere ( ? Should we ? )
          if ( altmidb < cptrop(j,i) - d_100 ) cycle
          tempb     = ctb3d(j,i,k)
          zenithb   = dacos(czen(j,i)) * raddeg
          cfactor   = crhob3d(j,i,k) * 1.D-03 * navgdr
          densb     = cfactor / amd
          C_Mb      = altmidb * 10.0D0 /(kb*tempb)
          depthab   = d_zero
          depthbb   = d_zero
          altaboveb = d_zero
          altbelowb = d_zero

          if ( ichjphcld == 1 ) then
            depthab   = sum(ctaucld(j,i,1:k-1,8))
            altaboveb = sum(cdzq(j,i,1:k-1) * ctaucld(j,i,1:k-1,8))
            depthbb   = sum(ctaucld(j,i,k+1:kz,8))
            altbelowb = sum(cdzq(j,i,k+1:kz) * ctaucld(j,i,k+1:kz,8))
            if ( depthab > d_zero ) altaboveb = altaboveb / depthab
            if ( depthbb > d_zero ) altbelowb = altbelowb / depthbb
          end if

          ! call the chemistry solver
          xrb(:)    = d_zero
          xrinb(:)  = d_zero
          xroutb(:) = d_zero
          ! 1 : initialise xrin with the concentrations from
          !     previous chemsolv step
          do ic = 1 , numchem !totch
            xrinb(ic) = chemall(j,i,k,ic)
          end do
          ! 2 : update input concentrations for transported species only
          do n = 1 , ntr !numchem !ntr
            if ( trac%indcb6(n) > 0 ) then
              xrinb(trac%indcb6(n)) = chib3d(j,i,k,n) * cfactor / trac%mw(n)
            end if
          end do
          ! update for water vapor
          wtrin = cqxb3d(j,i,k,iqv) * cfactor / amw

          call chemmain_cb6(calday,dtchsolv)

          ! save the concentrations of all species for next chemistry step
          do ic = 1 , numchem !totch
            chemall(j,i,k,ic) = xroutb(ic)
          end do
          ! Store photolysis rates for diagnostic
          do ic = 1 , nphoto
            jphoto(j,i,k,ic) = c_jvalb(1,ic)
          end do

          ! Now calculate chemical tendencies
          ! mole.cm-3.s-1  to kg.kg-1.s-1.ps (consistency with chiten unit)
          pfact = cpsb(j,i) / cfactor / dtchsolv

          do n = 1 , ntr !numchem !ntr
            if ( trac%indcb6(n) > 0 ) then
              !ashalaby
              if(xroutb(trac%indcb6(n)) .gt. 1.0e20 .or. &
                  xrinb(trac%indcb6(n)) .gt. 1.0e20) then
                write(stderr,*)xroutb(trac%indcb6(n)), &
                                    xrinb(trac%indcb6(n)),trac%indcb6(n)
                cerr = 1
                call fatal(__FILE__,__LINE__,'chemistry')
                !ashalaby_
                if ( cerr == 1 ) EXIT
              else
                chemten(j,i,k,n) = (xroutb(trac%indcb6(n)) - &
                                    xrinb(trac%indcb6(n))) * pfact * trac%mw(n)
              end if
            end if
          end do

        if ( cerr == 1 ) EXIT
        end do ! end i loop
      if ( cerr == 1 ) EXIT
      end do ! end k loop

    end subroutine chemistry_cb6

end module mod_che_chemistry
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
