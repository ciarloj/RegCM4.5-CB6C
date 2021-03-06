#include <misc.h>
#include <preproc.h>

module CNDecompMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDecompMod
!
! !DESCRIPTION:
! Module holding routines used in litter and soil decomposition model
! for coupled carbon-nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
   use clm_varpar   , only: nlevsoi
   use clm_varcon   , only: istsoil
   use spmdMod      , only: masterproc
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public:: CNDecompAlloc
!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated to vector data structures
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNDecompAlloc
!
! !INTERFACE:
subroutine CNDecompAlloc (lbc, ubc, num_soilc, filter_soilc, &
   num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use clm_varctl, only: irad
   use CNAllocationMod, only: CNAllocation
   use clm_time_manager, only: get_step_size
   use pft2colMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/15/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! column level
   real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: psisat(:,:)      ! soil water potential at saturation for CN code (MPa)
   real(r8), pointer :: soilpsi(:,:)     ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: cwdc(:)          ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)        ! (kgC/m2) litter labile C
   real(r8), pointer :: litr2c(:)        ! (kgC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)        ! (kgC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)        ! (kgC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)        ! (kgC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)        ! (kgC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)        ! (kgC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: cwdn(:)          ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)        ! (kgN/m2) litter labile N
   real(r8), pointer :: litr2n(:)        ! (kgN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)        ! (kgN/m2) litter lignin N
   integer, pointer :: clandunit(:)      ! index into landunit level quantities
   integer , pointer :: itypelun(:)      ! landunit type
   ! pft level
   real(r8), pointer :: rootfr(:,:)      ! fraction of roots in each soil layer  (nlevsoi)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: fpi(:)           ! fraction of potential immobilization (no units)
   real(r8), pointer :: cwdc_to_litr2c(:)
   real(r8), pointer :: cwdc_to_litr3c(:)
   real(r8), pointer :: litr1_hr(:)
   real(r8), pointer :: litr1c_to_soil1c(:)
   real(r8), pointer :: litr2_hr(:)
   real(r8), pointer :: litr2c_to_soil2c(:)
   real(r8), pointer :: litr3_hr(:)
   real(r8), pointer :: litr3c_to_soil3c(:)
   real(r8), pointer :: soil1_hr(:)
   real(r8), pointer :: soil1c_to_soil2c(:)
   real(r8), pointer :: soil2_hr(:)
   real(r8), pointer :: soil2c_to_soil3c(:)
   real(r8), pointer :: soil3_hr(:)
   real(r8), pointer :: soil3c_to_soil4c(:)
   real(r8), pointer :: soil4_hr(:)
   real(r8), pointer :: cwdn_to_litr2n(:)
   real(r8), pointer :: cwdn_to_litr3n(:)
   real(r8), pointer :: potential_immob(:)
   real(r8), pointer :: litr1n_to_soil1n(:)
   real(r8), pointer :: sminn_to_soil1n_l1(:)
   real(r8), pointer :: litr2n_to_soil2n(:)
   real(r8), pointer :: sminn_to_soil2n_l2(:)
   real(r8), pointer :: litr3n_to_soil3n(:)
   real(r8), pointer :: sminn_to_soil3n_l3(:)
   real(r8), pointer :: soil1n_to_soil2n(:)
   real(r8), pointer :: sminn_to_soil2n_s1(:)
   real(r8), pointer :: soil2n_to_soil3n(:)
   real(r8), pointer :: sminn_to_soil3n_s2(:)
   real(r8), pointer :: soil3n_to_soil4n(:)
   real(r8), pointer :: sminn_to_soil4n_s3(:)
   real(r8), pointer :: soil4n_to_sminn(:)
   real(r8), pointer :: sminn_to_denit_l1s1(:)
   real(r8), pointer :: sminn_to_denit_l2s2(:)
   real(r8), pointer :: sminn_to_denit_l3s3(:)
   real(r8), pointer :: sminn_to_denit_s1s2(:)
   real(r8), pointer :: sminn_to_denit_s2s3(:)
   real(r8), pointer :: sminn_to_denit_s3s4(:)
   real(r8), pointer :: sminn_to_denit_s4(:)
   real(r8), pointer :: sminn_to_denit_excess(:)
   real(r8), pointer :: gross_nmin(:)
   real(r8), pointer :: net_nmin(:)
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,j          !indices
   integer :: fc           !lake filter column index
   integer :: dtime        !time step (s)
   real(r8):: dt           !decomp timestep (seconds)
   real(r8):: dtd          !decomp timestep (days)
   real(r8), pointer:: fr(:,:)   !column-level rooting fraction by soil depth
   real(r8):: frw(lbc:ubc)          !rooting fraction weight
   real(r8):: t_scalar(lbc:ubc)     !soil temperature scalar for decomp
   real(r8):: minpsi, maxpsi        !limits for soil water scalar for decomp
   real(r8):: psi                   !temporary soilpsi for water scalar
   real(r8):: w_scalar(lbc:ubc)     !soil water scalar for decomp
   real(r8):: rate_scalar  !combined rate scalar for decomp
   real(r8):: cn_l1(lbc:ubc)        !C:N for litter 1
   real(r8):: cn_l2(lbc:ubc)        !C:N for litter 2
   real(r8):: cn_l3(lbc:ubc)        !C:N for litter 3
   real(r8):: cn_s1        !C:N for SOM 1
   real(r8):: cn_s2        !C:N for SOM 2
   real(r8):: cn_s3        !C:N for SOM 3
   real(r8):: cn_s4        !C:N for SOM 4
   real(r8):: rf_l1s1      !respiration fraction litter 1 -> SOM 1
   real(r8):: rf_l2s2      !respiration fraction litter 2 -> SOM 2
   real(r8):: rf_l3s3      !respiration fraction litter 3 -> SOM 3
   real(r8):: rf_s1s2      !respiration fraction SOM 1 -> SOM 2
   real(r8):: rf_s2s3      !respiration fraction SOM 2 -> SOM 3
   real(r8):: rf_s3s4      !respiration fraction SOM 3 -> SOM 4
   real(r8):: k_l1         !decomposition rate constant litter 1
   real(r8):: k_l2         !decomposition rate constant litter 2
   real(r8):: k_l3         !decomposition rate constant litter 3
   real(r8):: k_s1         !decomposition rate constant SOM 1
   real(r8):: k_s2         !decomposition rate constant SOM 2
   real(r8):: k_s3         !decomposition rate constant SOM 3
   real(r8):: k_s4         !decomposition rate constant SOM 3
   real(r8):: k_frag       !fragmentation rate constant CWD
   real(r8):: ck_l1        !corrected decomposition rate constant litter 1
   real(r8):: ck_l2        !corrected decomposition rate constant litter 2
   real(r8):: ck_l3        !corrected decomposition rate constant litter 3
   real(r8):: ck_s1        !corrected decomposition rate constant SOM 1
   real(r8):: ck_s2        !corrected decomposition rate constant SOM 2
   real(r8):: ck_s3        !corrected decomposition rate constant SOM 3
   real(r8):: ck_s4        !corrected decomposition rate constant SOM 3
   real(r8):: ck_frag      !corrected fragmentation rate constant CWD
   real(r8):: cwd_fcel     !cellulose fraction of coarse woody debris
   real(r8):: cwd_flig     !lignin fraction of coarse woody debris
   real(r8):: cwdc_loss    !fragmentation rate for CWD carbon (gC/m2/s)
   real(r8):: cwdn_loss    !fragmentation rate for CWD nitrogen (gN/m2/s)
   real(r8):: plitr1c_loss(lbc:ubc) !potential C loss from litter 1
   real(r8):: plitr2c_loss(lbc:ubc) !potential C loss from litter 2
   real(r8):: plitr3c_loss(lbc:ubc) !potential C loss from litter 3
   real(r8):: psoil1c_loss(lbc:ubc) !potential C loss from SOM 1
   real(r8):: psoil2c_loss(lbc:ubc) !potential C loss from SOM 2
   real(r8):: psoil3c_loss(lbc:ubc) !potential C loss from SOM 3
   real(r8):: psoil4c_loss(lbc:ubc) !potential C loss from SOM 4
   real(r8):: pmnf_l1s1(lbc:ubc)    !potential mineral N flux, litter 1 -> SOM 1
   real(r8):: pmnf_l2s2(lbc:ubc)    !potential mineral N flux, litter 2 -> SOM 2
   real(r8):: pmnf_l3s3(lbc:ubc)    !potential mineral N flux, litter 3 -> SOM 3
   real(r8):: pmnf_s1s2(lbc:ubc)    !potential mineral N flux, SOM 1 -> SOM 2
   real(r8):: pmnf_s2s3(lbc:ubc)    !potential mineral N flux, SOM 2 -> SOM 3
   real(r8):: pmnf_s3s4(lbc:ubc)    !potential mineral N flux, SOM 3 -> SOM 4
   real(r8):: pmnf_s4(lbc:ubc)      !potential mineral N flux, SOM 4
   real(r8):: immob(lbc:ubc)        !potential N immobilization
   real(r8):: ratio        !temporary variable
   real(r8):: dnp          !denitrification proportion
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   t_soisno              => clm3%g%l%c%ces%t_soisno
   psisat                => clm3%g%l%c%cps%psisat
   soilpsi               => clm3%g%l%c%cps%soilpsi
   cwdc                  => clm3%g%l%c%ccs%cwdc
   litr1c                => clm3%g%l%c%ccs%litr1c
   litr2c                => clm3%g%l%c%ccs%litr2c
   litr3c                => clm3%g%l%c%ccs%litr3c
   soil1c                => clm3%g%l%c%ccs%soil1c
   soil2c                => clm3%g%l%c%ccs%soil2c
   soil3c                => clm3%g%l%c%ccs%soil3c
   soil4c                => clm3%g%l%c%ccs%soil4c
   cwdn                  => clm3%g%l%c%cns%cwdn
   litr1n                => clm3%g%l%c%cns%litr1n
   litr2n                => clm3%g%l%c%cns%litr2n
   litr3n                => clm3%g%l%c%cns%litr3n
   fpi                   => clm3%g%l%c%cps%fpi
   cwdc_to_litr2c        => clm3%g%l%c%ccf%cwdc_to_litr2c
   cwdc_to_litr3c        => clm3%g%l%c%ccf%cwdc_to_litr3c
   litr1_hr              => clm3%g%l%c%ccf%litr1_hr
   litr1c_to_soil1c      => clm3%g%l%c%ccf%litr1c_to_soil1c
   litr2_hr              => clm3%g%l%c%ccf%litr2_hr
   litr2c_to_soil2c      => clm3%g%l%c%ccf%litr2c_to_soil2c
   litr3_hr              => clm3%g%l%c%ccf%litr3_hr
   litr3c_to_soil3c      => clm3%g%l%c%ccf%litr3c_to_soil3c
   soil1_hr              => clm3%g%l%c%ccf%soil1_hr
   soil1c_to_soil2c      => clm3%g%l%c%ccf%soil1c_to_soil2c
   soil2_hr              => clm3%g%l%c%ccf%soil2_hr
   soil2c_to_soil3c      => clm3%g%l%c%ccf%soil2c_to_soil3c
   soil3_hr              => clm3%g%l%c%ccf%soil3_hr
   soil3c_to_soil4c      => clm3%g%l%c%ccf%soil3c_to_soil4c
   soil4_hr              => clm3%g%l%c%ccf%soil4_hr
   cwdn_to_litr2n        => clm3%g%l%c%cnf%cwdn_to_litr2n
   cwdn_to_litr3n        => clm3%g%l%c%cnf%cwdn_to_litr3n
   potential_immob       => clm3%g%l%c%cnf%potential_immob
   litr1n_to_soil1n      => clm3%g%l%c%cnf%litr1n_to_soil1n
   sminn_to_soil1n_l1    => clm3%g%l%c%cnf%sminn_to_soil1n_l1
   litr2n_to_soil2n      => clm3%g%l%c%cnf%litr2n_to_soil2n
   sminn_to_soil2n_l2    => clm3%g%l%c%cnf%sminn_to_soil2n_l2
   litr3n_to_soil3n      => clm3%g%l%c%cnf%litr3n_to_soil3n
   sminn_to_soil3n_l3    => clm3%g%l%c%cnf%sminn_to_soil3n_l3
   soil1n_to_soil2n      => clm3%g%l%c%cnf%soil1n_to_soil2n
   sminn_to_soil2n_s1    => clm3%g%l%c%cnf%sminn_to_soil2n_s1
   soil2n_to_soil3n      => clm3%g%l%c%cnf%soil2n_to_soil3n
   sminn_to_soil3n_s2    => clm3%g%l%c%cnf%sminn_to_soil3n_s2
   soil3n_to_soil4n      => clm3%g%l%c%cnf%soil3n_to_soil4n
   sminn_to_soil4n_s3    => clm3%g%l%c%cnf%sminn_to_soil4n_s3
   soil4n_to_sminn       => clm3%g%l%c%cnf%soil4n_to_sminn
   sminn_to_denit_l1s1   => clm3%g%l%c%cnf%sminn_to_denit_l1s1
   sminn_to_denit_l2s2   => clm3%g%l%c%cnf%sminn_to_denit_l2s2
   sminn_to_denit_l3s3   => clm3%g%l%c%cnf%sminn_to_denit_l3s3
   sminn_to_denit_s1s2   => clm3%g%l%c%cnf%sminn_to_denit_s1s2
   sminn_to_denit_s2s3   => clm3%g%l%c%cnf%sminn_to_denit_s2s3
   sminn_to_denit_s3s4   => clm3%g%l%c%cnf%sminn_to_denit_s3s4
   sminn_to_denit_s4     => clm3%g%l%c%cnf%sminn_to_denit_s4
   sminn_to_denit_excess => clm3%g%l%c%cnf%sminn_to_denit_excess
   gross_nmin            => clm3%g%l%c%cnf%gross_nmin
   net_nmin              => clm3%g%l%c%cnf%net_nmin
   rootfr                => clm3%g%l%c%p%pps%rootfr
   clandunit             => clm3%g%l%c%landunit
   itypelun              => clm3%g%l%itype

   ! set time steps
   dtime = get_step_size()
   dt = float(irad)*dtime
   dtd = dt/86400.0_r8

   ! set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0)
   cn_s1 = 12.0_r8
   cn_s2 = 12.0_r8
   cn_s3 = 10.0_r8
   cn_s4 = 10.0_r8

   ! set respiration fractions for fluxes between compartments
   ! (from Biome-BGC v4.2.0)
   rf_l1s1 = 0.39_r8
   rf_l2s2 = 0.55_r8
   rf_l3s3 = 0.29_r8
   rf_s1s2 = 0.28_r8
   rf_s2s3 = 0.46_r8
   rf_s3s4 = 0.55

   ! set the cellulose and lignin fractions for coarse woody debris
   cwd_fcel = 0.76_r8
   cwd_flig = 0.24_r8

   ! set initial base rates for decomposition mass loss (1/day)
   ! (from Biome-BGC v4.2.0, using three SOM pools)
   ! Value inside log function is the discrete-time values for a
   ! daily time step model, and the result of the log function is
   ! the corresponding continuous-time decay rate (1/day), following
   ! Olson, 1963.
   k_l1 = -log(1.0_r8-0.7_r8)
   k_l2 = -log(1.0_r8-0.07_r8)
   k_l3 = -log(1.0_r8-0.014_r8)
   k_s1 = -log(1.0_r8-0.07_r8)
   k_s2 = -log(1.0_r8-0.014_r8)
   k_s3 = -log(1.0_r8-0.0014_r8)
   k_s4 = -log(1.0_r8-0.0001_r8)
   k_frag = -log(1.0_r8-0.001_r8)

   ! calculate the new discrete-time decay rate for model timestep
   k_l1 = 1.0_r8-exp(-k_l1*dtd)
   k_l2 = 1.0_r8-exp(-k_l2*dtd)
   k_l3 = 1.0_r8-exp(-k_l3*dtd)
   k_s1 = 1.0_r8-exp(-k_s1*dtd)
   k_s2 = 1.0_r8-exp(-k_s2*dtd)
   k_s3 = 1.0_r8-exp(-k_s3*dtd)
   k_s4 = 1.0_r8-exp(-k_s4*dtd)
   k_frag = 1.0_r8-exp(-k_frag*dtd)
   
   ! The following code implements the acceleration part of the AD spinup
   ! algorithm, by multiplying all of the decomposition base rates by 10.0.

#if (defined AD_SPINUP)
	k_l1 = k_l1 * 10._r8
   k_l2 = k_l2 * 10._r8
   k_l3 = k_l3 * 10._r8
   k_s1 = k_s1 * 10._r8
   k_s2 = k_s2 * 10._r8
   k_s3 = k_s3 * 10._r8
   k_s4 = k_s4 * 10._r8
   k_frag = k_frag * 10._r8
#endif

   ! calculate a weighting function by soil depth that depends on the
   ! fine root distribution per pft and depth and the pft weight on the column.
   ! This will be used to weight the temperature and water potential scalars
   ! for decomposition control.  There is currently no depth
   ! information associated with the litter and soil C and N pools,
   ! but the exponential decrease in fine root distribution with
   ! depth is a useful surrogate for the vertical distribution of
   ! soil organic matter.
   ! This uses the generic pft->column averaging routine, 2d version

   allocate(fr(lbc:ubc,nlevsoi))
   call p2c(nlevsoi,num_soilc,filter_soilc,rootfr,fr)

   ! the following normalizes values in fr so that they
   ! sum to 1.0 across all soil levels on a column

   frw(lbc:ubc) = 0._r8
   do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         frw(c) = frw(c) + fr(c,j)
      end do
   end do

   do j = 1,nlevsoi
!dir$ concurrent
!dir$ prefervector
!cdir nodep
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         if (frw(c) /= 0._r8) then
            fr(c,j) = fr(c,j) / frw(c)
         else
            fr(c,j) = 0._r8
         end if
      end do
   end do

   ! calculate rate constant scalar for soil temperature
   ! assuming that the base rate constants are assigned for non-moisture
   ! limiting conditions at 25 C. The function used here is taken from
   ! Lloyd, J., and J.A. Taylor, 1994. On the temperature dependence of
   ! soil respiration. Functional Ecology, 8:315-323.
   ! This equation is a modification of their eqn. 11, changing the base
   ! temperature from 10 C to 25 C, since most of the microcosm studies
   ! used to get the base decomp rates were controlled at 25 C.

   t_scalar(:) = 0._r8
   do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ! only modify scalar if the soil temperature is >= -10 C
         if (t_soisno(c,j) >= SHR_CONST_TKFRZ-10._r8) then
            t_scalar(c) = t_scalar(c) + exp(308.56_r8*((1.0_r8/71.02_r8)-(1.0_r8/(t_soisno(c,j)-227.13_r8))))*fr(c,j)
         end if
      end do
   end do

   ! calculate the rate constant scalar for soil water content.
   ! Uses the log relationship with water potential given in
   ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
   ! a comparison of models. Ecology, 68(5):1190-1200.
   ! and supported by data in
   ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
   ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

   minpsi = -10.0_r8;
   w_scalar(:) = 0._r8
   do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         maxpsi = psisat(c,j)
         psi = min(soilpsi(c,j),maxpsi)
         ! decomp only if soilpsi is higher than minpsi
         if (psi > minpsi) then
            w_scalar(c) = w_scalar(c) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j);
         end if
      end do
   end do

   ! set initial values for potential C and N fluxes
   plitr1c_loss(:) = 0._r8
   plitr2c_loss(:) = 0._r8
   plitr3c_loss(:) = 0._r8
   psoil1c_loss(:) = 0._r8
   psoil2c_loss(:) = 0._r8
   psoil3c_loss(:) = 0._r8
   psoil4c_loss(:) = 0._r8
   pmnf_l1s1(:) = 0._r8
   pmnf_l2s2(:) = 0._r8
   pmnf_l3s3(:) = 0._r8
   pmnf_s1s2(:) = 0._r8
   pmnf_s2s3(:) = 0._r8
   pmnf_s3s4(:) = 0._r8
   pmnf_s4(:) = 0._r8

   ! column loop to calculate potential decomp rates and total immobilization
   ! demand.
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! calculate litter compartment C:N ratios
      if (litr1n(c) > 0._r8) cn_l1(c) = litr1c(c)/litr1n(c)
      if (litr2n(c) > 0._r8) cn_l2(c) = litr2c(c)/litr2n(c)
      if (litr3n(c) > 0._r8) cn_l3(c) = litr3c(c)/litr3n(c)

      ! calculate the final rate scalar as the product of temperature and water
      ! rate scalars, and correct the base decomp rates

      rate_scalar = t_scalar(c) * w_scalar(c)
      ck_l1 = k_l1 * rate_scalar
      ck_l2 = k_l2 * rate_scalar
      ck_l3 = k_l3 * rate_scalar
      ck_s1 = k_s1 * rate_scalar
      ck_s2 = k_s2 * rate_scalar
      ck_s3 = k_s3 * rate_scalar
      ck_s4 = k_s4 * rate_scalar
      ck_frag = k_frag * rate_scalar

      ! calculate the non-nitrogen-limited fluxes
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.

      ! CWD fragmentation -> litter pools
      cwdc_loss = cwdc(c) * ck_frag / dt
      cwdc_to_litr2c(c) = cwdc_loss * cwd_fcel
      cwdc_to_litr3c(c) = cwdc_loss * cwd_flig
      cwdn_loss = cwdn(c) * ck_frag / dt
      cwdn_to_litr2n(c) = cwdn_loss * cwd_fcel
      cwdn_to_litr3n(c) = cwdn_loss * cwd_flig

      ! litter 1 -> SOM 1
      if (litr1c(c) > 0._r8) then
         plitr1c_loss(c) = litr1c(c) * ck_l1 / dt
         ratio = 0._r8
         if (litr1n(c) > 0._r8) ratio = cn_s1/cn_l1(c)
         pmnf_l1s1(c) = (plitr1c_loss(c) * (1.0_r8 - rf_l1s1 - ratio))/cn_s1
      end if

      ! litter 2 -> SOM 2
      if (litr2c(c) > 0._r8) then
         plitr2c_loss(c) = litr2c(c) * ck_l2 / dt
         ratio = 0._r8
         if (litr2n(c) > 0._r8) ratio = cn_s2/cn_l2(c)
         pmnf_l2s2(c) = (plitr2c_loss(c) * (1.0_r8 - rf_l2s2 - ratio))/cn_s2
      end if

      ! litter 3 -> SOM 3
      if (litr3c(c) > 0._r8) then
         plitr3c_loss(c) = litr3c(c) * ck_l3 / dt
         ratio = 0._r8
         if (litr3n(c) > 0._r8) ratio = cn_s3/cn_l3(c)
         pmnf_l3s3(c) = (plitr3c_loss(c) * (1.0_r8 - rf_l3s3 - ratio))/cn_s3
      end if

      ! SOM 1 -> SOM 2
      if (soil1c(c) > 0._r8) then
         psoil1c_loss(c) = soil1c(c) * ck_s1 / dt
         pmnf_s1s2(c) = (psoil1c_loss(c) * (1.0_r8 - rf_s1s2 - (cn_s2/cn_s1)))/cn_s2
      end if

      ! SOM 2 -> SOM 3
      if (soil2c(c) > 0._r8) then
         psoil2c_loss(c) = soil2c(c) * ck_s2 / dt
         pmnf_s2s3(c) = (psoil2c_loss(c) * (1.0_r8 - rf_s2s3 - (cn_s3/cn_s2)))/cn_s3
      end if

      ! SOM 3 -> SOM 4
      if (soil3c(c) > 0._r8) then
         psoil3c_loss(c) = soil3c(c) * ck_s3 / dt
         pmnf_s3s4(c) = (psoil3c_loss(c) * (1.0_r8 - rf_s3s4 - (cn_s4/cn_s3)))/cn_s4
      end if

      ! Loss from SOM 4 is entirely respiration (no downstream pool)
      if (soil4c(c) > 0._r8) then
         psoil4c_loss(c) = soil4c(c) * ck_s4 / dt
         pmnf_s4(c) = -psoil4c_loss(c)/cn_s4
      end if

      ! Sum up all the potential immobilization fluxes (positive pmnf flux)
      ! and all the mineralization fluxes (negative pmnf flux)

      immob(c) = 0._r8
      ! litter 1 -> SOM 1
      if (pmnf_l1s1(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_l1s1(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 -> SOM 2
      if (pmnf_l2s2(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_l2s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 -> SOM 3
      if (pmnf_l3s3(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_l3s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 -> SOM 2
      if (pmnf_s1s2(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_s1s2(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 -> SOM 3
      if (pmnf_s2s3(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_s2s3(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 -> SOM 4
      if (pmnf_s3s4(c) > 0._r8) then
         immob(c) = immob(c) + pmnf_s3s4(c)
      else
         gross_nmin(c) = gross_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4
      gross_nmin(c) = gross_nmin(c) - pmnf_s4(c)

      potential_immob(c) = immob(c)

   end do ! end column loop

   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.

   call CNAllocation(lbc,ubc,num_soilc,filter_soilc,num_soilp,filter_soilp)

   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N

   dnp = 0.01_r8

!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! upon return from CNAllocation, the fraction of potential immobilization
      ! has been set (cps%fpi). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.

      ! litter 1 fluxes (labile pool)
      if (litr1c(c) > 0._r8) then
         if (pmnf_l1s1(c) > 0._r8) then
            plitr1c_loss(c) = plitr1c_loss(c) * fpi(c)
            pmnf_l1s1(c) = pmnf_l1s1(c) * fpi(c)
            sminn_to_denit_l1s1(c) = 0._r8
         else
            sminn_to_denit_l1s1(c) = -dnp * pmnf_l1s1(c)
         end if
         litr1_hr(c) = rf_l1s1 * plitr1c_loss(c)
         litr1c_to_soil1c(c) = (1._r8 - rf_l1s1) * plitr1c_loss(c)
         if (litr1n(c) > 0._r8) litr1n_to_soil1n(c) = plitr1c_loss(c) / cn_l1(c)
         sminn_to_soil1n_l1(c) = pmnf_l1s1(c)
         net_nmin(c) = net_nmin(c) - pmnf_l1s1(c)
      end if

      ! litter 2 fluxes (cellulose pool)
      if (litr2c(c) > 0._r8) then
         if (pmnf_l2s2(c) > 0._r8) then
            plitr2c_loss(c) = plitr2c_loss(c) * fpi(c)
            pmnf_l2s2(c) = pmnf_l2s2(c) * fpi(c)
            sminn_to_denit_l2s2(c) = 0._r8
         else
            sminn_to_denit_l2s2(c) = -dnp * pmnf_l2s2(c)
         end if
         litr2_hr(c) = rf_l2s2 * plitr2c_loss(c)
         litr2c_to_soil2c(c) = (1._r8 - rf_l2s2) * plitr2c_loss(c)
         if (litr2n(c) > 0._r8) litr2n_to_soil2n(c) = plitr2c_loss(c) / cn_l2(c)
         sminn_to_soil2n_l2(c) = pmnf_l2s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_l2s2(c)
      end if

      ! litter 3 fluxes (lignin pool)
      if (litr3c(c) > 0._r8) then
         if (pmnf_l3s3(c) > 0._r8) then
            plitr3c_loss(c) = plitr3c_loss(c) * fpi(c)
            pmnf_l3s3(c) = pmnf_l3s3(c) * fpi(c)
            sminn_to_denit_l3s3(c) = 0._r8
         else
            sminn_to_denit_l3s3(c) = -dnp * pmnf_l3s3(c)
         end if
         litr3_hr(c) = rf_l3s3 * plitr3c_loss(c)
         litr3c_to_soil3c(c) = (1._r8 - rf_l3s3) * plitr3c_loss(c)
         if (litr3n(c) > 0._r8) litr3n_to_soil3n(c) = plitr3c_loss(c) / cn_l3(c)
         sminn_to_soil3n_l3(c) = pmnf_l3s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_l3s3(c)
      end if

      ! SOM 1 fluxes (fast rate soil organic matter pool)
      if (soil1c(c) > 0._r8) then
         if (pmnf_s1s2(c) > 0._r8) then
            psoil1c_loss(c) = psoil1c_loss(c) * fpi(c)
            pmnf_s1s2(c) = pmnf_s1s2(c) * fpi(c)
            sminn_to_denit_s1s2(c) = 0._r8
         else
            sminn_to_denit_s1s2(c) = -dnp * pmnf_s1s2(c)
         end if
         soil1_hr(c) = rf_s1s2 * psoil1c_loss(c)
         soil1c_to_soil2c(c) = (1._r8 - rf_s1s2) * psoil1c_loss(c)
         soil1n_to_soil2n(c) = psoil1c_loss(c) / cn_s1
         sminn_to_soil2n_s1(c) = pmnf_s1s2(c)
         net_nmin(c) = net_nmin(c) - pmnf_s1s2(c)
      end if

      ! SOM 2 fluxes (medium rate soil organic matter pool)
      if (soil2c(c) > 0._r8) then
         if (pmnf_s2s3(c) > 0._r8) then
            psoil2c_loss(c) = psoil2c_loss(c) * fpi(c)
            pmnf_s2s3(c) = pmnf_s2s3(c) * fpi(c)
            sminn_to_denit_s2s3(c) = 0._r8
         else
            sminn_to_denit_s2s3(c) = -dnp * pmnf_s2s3(c)
         end if
         soil2_hr(c) = rf_s2s3 * psoil2c_loss(c)
         soil2c_to_soil3c(c) = (1._r8 - rf_s2s3) * psoil2c_loss(c)
         soil2n_to_soil3n(c) = psoil2c_loss(c) / cn_s2
         sminn_to_soil3n_s2(c) = pmnf_s2s3(c)
         net_nmin(c) = net_nmin(c) - pmnf_s2s3(c)
      end if

      ! SOM 3 fluxes (slow rate soil organic matter pool)
      if (soil3c(c) > 0._r8) then
         if (pmnf_s3s4(c) > 0._r8) then
            psoil3c_loss(c) = psoil3c_loss(c) * fpi(c)
            pmnf_s3s4(c) = pmnf_s3s4(c) * fpi(c)
            sminn_to_denit_s3s4(c) = 0._r8
         else
            sminn_to_denit_s3s4(c) = -dnp * pmnf_s3s4(c)
         end if
         soil3_hr(c) = rf_s3s4 * psoil3c_loss(c)
         soil3c_to_soil4c(c) = (1._r8 - rf_s3s4) * psoil3c_loss(c)
         soil3n_to_soil4n(c) = psoil3c_loss(c) / cn_s3
         sminn_to_soil4n_s3(c) = pmnf_s3s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s3s4(c)
      end if

      ! SOM 4 fluxes (slowest rate soil organic matter pool)
      if (soil4c(c) > 0._r8) then
         soil4_hr(c) = psoil4c_loss(c)
         soil4n_to_sminn(c) = psoil4c_loss(c) / cn_s4
         sminn_to_denit_s4(c) = -dnp * pmnf_s4(c)
         net_nmin(c) = net_nmin(c) - pmnf_s4(c)
      end if

   end do

   deallocate(fr)

end subroutine CNDecompAlloc
!-----------------------------------------------------------------------
end module CNDecompMod
