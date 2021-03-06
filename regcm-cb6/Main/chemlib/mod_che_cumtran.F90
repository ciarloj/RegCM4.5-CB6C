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

module mod_che_cumtran
!
! Tracer convective transport
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_che_common

  implicit none

  private

  public :: cumtran

  contains

  subroutine cumtran
    implicit none
    real(rk8) :: chiabar , chibbar , deltas , cumfrc
    real(rk8) :: chiabar_n , chiabar_a , chibbar_n , chibbar_a
    integer(ik4) :: i , j , k , kctop , n

    if ( ichdiag == 1 ) then
     do j = jci1 , jci2
        do i = ici1 , ici2
         chiten0(j,i,:,:) = chib(j,i,:,:)
        end do
      end do
    end if

    do n = 1 , ntr
      do j = jci1 , jci2
        do i = ici1 , ici2
          if ( kcumtop(j,i) > 0 ) then
            deltas = d_zero
            chiabar = d_zero
            chibbar = d_zero
            kctop = max0(kcumtop(j,i),4)
            do k = kctop , kz
              deltas = deltas + dsigma(k)
              chiabar = chiabar + chia(j,i,k,n)*dsigma(k)
              chibbar = chibbar + chib(j,i,k,n)*dsigma(k)
            end do
            do k = kctop , kz
 !            cumfrc =  ccldfra(j,i,k) - cfcc(j,i,k)
              cumfrc =  convcldfra (j,i,k)
              chia(j,i,k,n) = chia(j,i,k,n)*(d_one-cumfrc)+cumfrc*chiabar/deltas
              chib(j,i,k,n) = chib(j,i,k,n)*(d_one-cumfrc)+cumfrc*chibbar/deltas
            end do
          end if
        end do
      end do
    end do
    ! here calculate a pseudo tendency.
    ! factor 2 is added since we are out of leap frog
    if ( ichdiag == 1 ) then
      do j = jci1 , jci2
        do i = ici1 , ici2
        cconvdiag(j,i,:,:)  = cconvdiag(j,i,:,:) +   &
          ( chib(j,i,:,:) - chiten0(j,i,:,:))/dt  * d_two * cfdout
       end do
     end do
    end if
  end subroutine cumtran
!
end module mod_che_cumtran
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
