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

module mod_diffusion
!
! Diffusion calculations
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_atm_interface
  use mod_runparams
  use mod_mppparam
  use mod_service

  implicit none

  private

  interface diffu_x
    module procedure diffu_x3d
    module procedure diffu_x4d
  end interface

  public :: diffu_d , diffu_x
!
  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     These subroutines computes the diffusion term for decoupled     c
!     variable on constant sigma surface.                             c
!                                                                     c
!     ften    : tendency for variable                                 c
!                                                                     c
!     xkc     : horizontal diffusion coefficient                      c
!                                                                     c
!     f       : variable f on cross/dot points                        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine diffu_d(uten,vten,u,v,press,xkc)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: u , v
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(rk8) , pointer , dimension(:,:) , intent(in) :: press
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: uten , vten
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme for interior:
    !
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jdii1 , jdii2
          uten(j,i,k) = uten(j,i,k) - xkc(j,i,k) *      &
                rdxsq*(u(j+2,i,k)+u(j-2,i,k)+u(j,i+2,k)+u(j,i-2,k)  - &
               d_four*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k)) + &
              d_twelve*u(j,i,k))*press(j,i)
          vten(j,i,k) = vten(j,i,k) - xkc(j,i,k) *      &
                rdxsq*(v(j+2,i,k)+v(j-2,i,k)+v(j,i+2,k)+v(j,i-2,k)  - &
               d_four*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k)) + &
              d_twelve*v(j,i,k))*press(j,i)
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%has_bdyleft ) then
      j = jdi1
      do k = 1 , kz
        do i = idi1 , idi2
          uten(j,i,k) = uten(j,i,k) + xkc(j,i,k) *    &
                rdxsq*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k) - &
                        d_four*u(j,i,k))*press(j,i)
          vten(j,i,k) = vten(j,i,k) + xkc(j,i,k) *    &
                rdxsq*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k) - &
                        d_four*v(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      j = jdi2
      do k = 1 , kz
        do i = idi1 , idi2
          uten(j,i,k) = uten(j,i,k) + xkc(j,i,k) *    &
                rdxsq*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k) - &
                        d_four*u(j,i,k))*press(j,i)
          vten(j,i,k) = vten(j,i,k) + xkc(j,i,k) *    &
                rdxsq*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k) - &
                        d_four*v(j,i,k))*press(j,i)
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%has_bdybottom ) then
      i = idi1
      do k = 1 , kz
        do j = jdi1 , jdi2
          uten(j,i,k) = uten(j,i,k) + xkc(j,i,k) *         &
                rdxsq*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k) - &
                d_four*u(j,i,k))*press(j,i)
          vten(j,i,k) = vten(j,i,k) + xkc(j,i,k) *         &
                rdxsq*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k) - &
                d_four*v(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      i = idi2
      do k = 1 , kz
        do j = jdi1 , jdi2
          uten(j,i,k) = uten(j,i,k) + xkc(j,i,k) *         &
                rdxsq*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k) - &
                d_four*u(j,i,k))*press(j,i)
          vten(j,i,k) = vten(j,i,k) + xkc(j,i,k) *         &
                rdxsq*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k) - &
                d_four*v(j,i,k))*press(j,i)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_d
!
  subroutine diffu_x3d(ften,f,press,xkc,kmax)
!
    implicit none
!
    integer(ik4) , intent(in) :: kmax
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ften
    real(rk8) , pointer , dimension(:,:) , intent(in) :: press
!
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme for interior:
    !
    do k = 1 , kmax
      do i = icii1 , icii2
        do j = jcii1 , jcii2
          ften(j,i,k) = ften(j,i,k) - xkc(j,i,k) *    &
                      rdxsq*(f(j+2,i,k)+f(j-2,i,k) +  &
                             f(j,i+2,k)+f(j,i-2,k) -  &
                     d_four*(f(j+1,i,k)+f(j-1,i,k) +  &
                             f(j,i+1,k)+f(j,i-1,k)) + &
                    d_twelve*f(j,i,k))*press(j,i)
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%has_bdyleft ) then
      j = jci1
      do k = 1 , kmax
        do i = ici1 , ici2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *    &
                        rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                               f(j,i+1,k)+f(j,i-1,k) - &
                        d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      j = jci2
      do k = 1 , kmax
        do i = ici1 , ici2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                        rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                               f(j,i+1,k)+f(j,i-1,k) - &
                        d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%has_bdybottom ) then
      i = ici1
      do k = 1 , kmax
        do j = jci1 , jci2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                      rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                             f(j,i+1,k)+f(j,i-1,k) - &
                      d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      i = ici2
      do k = 1 , kmax
        do j = jci1 , jci2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                      rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                             f(j,i+1,k)+f(j,i-1,k) - &
                      d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x3d
!
  subroutine diffu_x4d(ften,f,press,xkc,kmax,n4)
    implicit none
    integer(ik4) , intent(in) :: kmax
    integer(ik4) , optional , intent(in) :: n4
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: f
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ften
    real(rk8) , pointer , dimension(:,:) , intent(in) :: press

    integer(ik4) :: i , j , k , n , n1 , n2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( present(n4) ) then
      n1 = n4
      n2 = n4
    else
      n1 = lbound(f,4)
      n2 = ubound(f,4)
    end if
    !
    ! fourth-order scheme for interior:
    !
    do n = n1 , n2
      do k = 1 , kmax
        do i = icii1 , icii2
          do j = jcii1 , jcii2
            ften(j,i,k,n) = ften(j,i,k,n) - xkc(j,i,k) *    &
                        rdxsq*(f(j+2,i,k,n)+f(j-2,i,k,n) +  &
                               f(j,i+2,k,n)+f(j,i-2,k,n) -  &
                       d_four*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                               f(j,i+1,k,n)+f(j,i-1,k,n)) + &
                      d_twelve*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%has_bdyleft ) then
      j = jci1
      do n = n1 , n2
        do k = 1 , kmax
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                                 f(j,i+1,k,n)+f(j,i-1,k,n) - &
                          d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      j = jci2
      do n = n1 , n2
        do k = 1 , kmax
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                                 f(j,i+1,k,n)+f(j,i-1,k,n) - &
                          d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%has_bdybottom ) then
      i = ici1
      do n = n1 , n2
        do k = 1 , kmax
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *   &
                        rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                               f(j,i+1,k,n)+f(j,i-1,k,n) - &
                        d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      i = ici2
      do n = n1 , n2
        do k = 1 , kmax
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *   &
                        rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                               f(j,i+1,k,n)+f(j,i-1,k,n) - &
                        d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x4d

end module mod_diffusion
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
