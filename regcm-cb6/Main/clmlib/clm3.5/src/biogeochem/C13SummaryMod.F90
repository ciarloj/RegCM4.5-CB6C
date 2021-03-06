#include <misc.h>
#include <preproc.h>

module C13SummaryMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13SummaryMod
!
! !DESCRIPTION:
! Module for isotope carbon summary calculations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon  , only: istsoil
    use spmdMod     , only: masterproc
    use clm_varpar  , only: nlevsoi
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: C13Summary
!
! !REVISION HISTORY:
! 7/13/2005: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Summary
!
! !INTERFACE:
subroutine C13Summary(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform pft and column-level carbon
! summary calculations
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   real(r8), pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
   real(r8), pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: hr(:)            ! (gC/m2/s) total heterotrophic respiration
   real(r8), pointer :: litfire(:)       ! (gC/m2/s) litter fire losses
   real(r8), pointer :: lithr(:)         ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: litr1_hr(:)       
   real(r8), pointer :: litr2_hr(:)        
   real(r8), pointer :: litr3_hr(:)        
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: m_litr1c_to_fire(:)             
   real(r8), pointer :: m_litr2c_to_fire(:)             
   real(r8), pointer :: m_litr3c_to_fire(:)             
   real(r8), pointer :: nee(:)           ! (gC/m2/s) net ecosystem exchange of carbon, includes fire flux, positive for source
   real(r8), pointer :: nep(:)           ! (gC/m2/s) net ecosystem production, excludes fire flux, positive for sink
   real(r8), pointer :: col_ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: col_gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: col_npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: col_pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: col_rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: col_vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: soil1_hr(:)        
   real(r8), pointer :: soil2_hr(:)        
   real(r8), pointer :: soil3_hr(:) 
   real(r8), pointer :: soil4_hr(:) 
   real(r8), pointer :: somfire(:)       ! (gC/m2/s) soil organic matter fire losses
   real(r8), pointer :: somhr(:)         ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: sr(:)            ! (gC/m2/s) total soil respiration (HR + root resp)
   real(r8), pointer :: totfire(:)       ! (gC/m2/s) total ecosystem fire losses
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: col_totpftc(:)        ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: col_totvegc(:)        ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
   real(r8), pointer :: agnpp(:)          ! (gC/m2/s) aboveground NPP
   real(r8), pointer :: ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: bgnpp(:)          ! (gC/m2/s) belowground NPP
   real(r8), pointer :: cpool_deadcroot_gr(:)        
   real(r8), pointer :: cpool_deadcroot_storage_gr(:)
   real(r8), pointer :: cpool_deadstem_gr(:)         
   real(r8), pointer :: cpool_deadstem_storage_gr(:) 
   real(r8), pointer :: cpool_froot_gr(:)            
   real(r8), pointer :: cpool_froot_storage_gr(:)    
   real(r8), pointer :: cpool_leaf_gr(:)             
   real(r8), pointer :: cpool_leaf_storage_gr(:)     
   real(r8), pointer :: cpool_livecroot_gr(:)        
   real(r8), pointer :: cpool_livecroot_storage_gr(:)
   real(r8), pointer :: cpool_livestem_gr(:)         
   real(r8), pointer :: cpool_livestem_storage_gr(:) 
   real(r8), pointer :: cpool_to_deadcrootc(:)        
   real(r8), pointer :: cpool_to_deadstemc(:)         
   real(r8), pointer :: cpool_to_frootc(:)            
   real(r8), pointer :: cpool_to_leafc(:)             
   real(r8), pointer :: cpool_to_livecrootc(:)        
   real(r8), pointer :: cpool_to_livestemc(:)         
   real(r8), pointer :: current_gr(:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:) 
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)       
   real(r8), pointer :: froot_mr(:)     
   real(r8), pointer :: froot_curmr(:)     
   real(r8), pointer :: froot_xsmr(:)     
   real(r8), pointer :: gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: gr(:)             ! (gC/m2/s) total growth respiration
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)         
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: litfall(:)        ! (gC/m2/s) litterfall (leaves and fine roots)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:) 
   real(r8), pointer :: livestem_mr(:)  
   real(r8), pointer :: livestem_curmr(:)  
   real(r8), pointer :: livestem_xsmr(:)  
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:) 
   real(r8), pointer :: m_deadcrootc_to_fire(:)         
   real(r8), pointer :: m_deadcrootc_to_litter(:)           
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)         
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)  
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)  
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)            
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:) 
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:) 
   real(r8), pointer :: m_frootc_storage_to_fire(:)     
   real(r8), pointer :: m_frootc_storage_to_litter(:)     
   real(r8), pointer :: m_frootc_to_fire(:)             
   real(r8), pointer :: m_frootc_to_litter(:)             
   real(r8), pointer :: m_frootc_xfer_to_fire(:)    
   real(r8), pointer :: m_frootc_xfer_to_litter(:)    
   real(r8), pointer :: m_gresp_storage_to_fire(:)      
   real(r8), pointer :: m_gresp_storage_to_litter(:)      
   real(r8), pointer :: m_gresp_xfer_to_fire(:)    
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_fire(:)      
   real(r8), pointer :: m_leafc_storage_to_litter(:)      
   real(r8), pointer :: m_leafc_to_fire(:)             
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_fire(:)     
   real(r8), pointer :: m_leafc_xfer_to_litter(:)    
   real(r8), pointer :: m_livecrootc_storage_to_fire(:) 
   real(r8), pointer :: m_livecrootc_storage_to_litter(:) 
   real(r8), pointer :: m_livecrootc_to_fire(:)         
   real(r8), pointer :: m_livecrootc_to_litter(:)           
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)  
   real(r8), pointer :: m_livestemc_storage_to_litter(:)  
   real(r8), pointer :: m_livestemc_to_fire(:)          
   real(r8), pointer :: m_livestemc_to_litter(:)            
   real(r8), pointer :: m_livestemc_xfer_to_fire(:) 
   real(r8), pointer :: m_livestemc_xfer_to_litter(:) 
   real(r8), pointer :: mr(:)             ! (gC/m2/s) maintenance respiration
   real(r8), pointer :: npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: psnshade_to_cpool(:)
   real(r8), pointer :: psnsun_to_cpool(:) 
   real(r8), pointer :: rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: storage_gr(:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
   real(r8), pointer :: transfer_deadcroot_gr(:)
   real(r8), pointer :: transfer_deadstem_gr(:)      
   real(r8), pointer :: transfer_froot_gr(:)         
   real(r8), pointer :: transfer_gr(:)    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
   real(r8), pointer :: transfer_leaf_gr(:)          
   real(r8), pointer :: transfer_livecroot_gr(:)     
   real(r8), pointer :: transfer_livestem_gr(:)      
   real(r8), pointer :: vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    !(gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    !(gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: tempsum_npp(:)          !temporary annual sum of NPP (gC/m2/yr)
!
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p          ! indices
   integer :: fp,fc        ! lake filter indices

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers
    col_fire_closs                 => clm3%g%l%c%cc13f%col_fire_closs
    er                             => clm3%g%l%c%cc13f%er
    hr                             => clm3%g%l%c%cc13f%hr
    litfire                        => clm3%g%l%c%cc13f%litfire
    lithr                          => clm3%g%l%c%cc13f%lithr
    litr1_hr                       => clm3%g%l%c%cc13f%litr1_hr
    litr2_hr                       => clm3%g%l%c%cc13f%litr2_hr
    litr3_hr                       => clm3%g%l%c%cc13f%litr3_hr
    m_cwdc_to_fire                 => clm3%g%l%c%cc13f%m_cwdc_to_fire
    m_litr1c_to_fire               => clm3%g%l%c%cc13f%m_litr1c_to_fire
    m_litr2c_to_fire               => clm3%g%l%c%cc13f%m_litr2c_to_fire
    m_litr3c_to_fire               => clm3%g%l%c%cc13f%m_litr3c_to_fire
    nee                            => clm3%g%l%c%cc13f%nee
    nep                            => clm3%g%l%c%cc13f%nep
    col_ar                         => clm3%g%l%c%cc13f%pcf_a%ar
    col_gpp                        => clm3%g%l%c%cc13f%pcf_a%gpp
    col_npp                        => clm3%g%l%c%cc13f%pcf_a%npp
    col_pft_fire_closs             => clm3%g%l%c%cc13f%pcf_a%pft_fire_closs
    col_rr                         => clm3%g%l%c%cc13f%pcf_a%rr
    col_vegfire                    => clm3%g%l%c%cc13f%pcf_a%vegfire
    soil1_hr                       => clm3%g%l%c%cc13f%soil1_hr
    soil2_hr                       => clm3%g%l%c%cc13f%soil2_hr
    soil3_hr                       => clm3%g%l%c%cc13f%soil3_hr
    soil4_hr                       => clm3%g%l%c%cc13f%soil4_hr
    somfire                        => clm3%g%l%c%cc13f%somfire
    somhr                          => clm3%g%l%c%cc13f%somhr
    sr                             => clm3%g%l%c%cc13f%sr
    totfire                        => clm3%g%l%c%cc13f%totfire
    cwdc                           => clm3%g%l%c%cc13s%cwdc
    litr1c                         => clm3%g%l%c%cc13s%litr1c
    litr2c                         => clm3%g%l%c%cc13s%litr2c
    litr3c                         => clm3%g%l%c%cc13s%litr3c
    col_totpftc                    => clm3%g%l%c%cc13s%pcs_a%totpftc
    col_totvegc                    => clm3%g%l%c%cc13s%pcs_a%totvegc
    soil1c                         => clm3%g%l%c%cc13s%soil1c
    soil2c                         => clm3%g%l%c%cc13s%soil2c
    soil3c                         => clm3%g%l%c%cc13s%soil3c
    soil4c                         => clm3%g%l%c%cc13s%soil4c
    totcolc                        => clm3%g%l%c%cc13s%totcolc
    totecosysc                     => clm3%g%l%c%cc13s%totecosysc
    totlitc                        => clm3%g%l%c%cc13s%totlitc
    totsomc                        => clm3%g%l%c%cc13s%totsomc
    agnpp                          => clm3%g%l%c%p%pc13f%agnpp
    ar                             => clm3%g%l%c%p%pc13f%ar
    bgnpp                          => clm3%g%l%c%p%pc13f%bgnpp
    cpool_deadcroot_gr             => clm3%g%l%c%p%pc13f%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => clm3%g%l%c%p%pc13f%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => clm3%g%l%c%p%pc13f%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => clm3%g%l%c%p%pc13f%cpool_deadstem_storage_gr
    cpool_froot_gr                 => clm3%g%l%c%p%pc13f%cpool_froot_gr
    cpool_froot_storage_gr         => clm3%g%l%c%p%pc13f%cpool_froot_storage_gr
    cpool_leaf_gr                  => clm3%g%l%c%p%pc13f%cpool_leaf_gr
    cpool_leaf_storage_gr          => clm3%g%l%c%p%pc13f%cpool_leaf_storage_gr
    cpool_livecroot_gr             => clm3%g%l%c%p%pc13f%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => clm3%g%l%c%p%pc13f%cpool_livecroot_storage_gr
    cpool_livestem_gr              => clm3%g%l%c%p%pc13f%cpool_livestem_gr
    cpool_livestem_storage_gr      => clm3%g%l%c%p%pc13f%cpool_livestem_storage_gr
    cpool_to_deadcrootc            => clm3%g%l%c%p%pc13f%cpool_to_deadcrootc
    cpool_to_deadstemc             => clm3%g%l%c%p%pc13f%cpool_to_deadstemc
    cpool_to_frootc                => clm3%g%l%c%p%pc13f%cpool_to_frootc
    cpool_to_leafc                 => clm3%g%l%c%p%pc13f%cpool_to_leafc
    cpool_to_livecrootc            => clm3%g%l%c%p%pc13f%cpool_to_livecrootc
    cpool_to_livestemc             => clm3%g%l%c%p%pc13f%cpool_to_livestemc
    current_gr                     => clm3%g%l%c%p%pc13f%current_gr
    deadcrootc_xfer_to_deadcrootc  => clm3%g%l%c%p%pc13f%deadcrootc_xfer_to_deadcrootc
    deadstemc_xfer_to_deadstemc    => clm3%g%l%c%p%pc13f%deadstemc_xfer_to_deadstemc
    frootc_to_litter               => clm3%g%l%c%p%pc13f%frootc_to_litter
    frootc_xfer_to_frootc          => clm3%g%l%c%p%pc13f%frootc_xfer_to_frootc
    froot_mr                       => clm3%g%l%c%p%pc13f%froot_mr
    froot_curmr                    => clm3%g%l%c%p%pc13f%froot_curmr
    froot_xsmr                     => clm3%g%l%c%p%pc13f%froot_xsmr
    gpp                            => clm3%g%l%c%p%pc13f%gpp
    gr                             => clm3%g%l%c%p%pc13f%gr
    leafc_to_litter                => clm3%g%l%c%p%pc13f%leafc_to_litter
    leafc_xfer_to_leafc            => clm3%g%l%c%p%pc13f%leafc_xfer_to_leafc
    leaf_mr                        => clm3%g%l%c%p%pc13f%leaf_mr
    leaf_curmr                     => clm3%g%l%c%p%pc13f%leaf_curmr
    leaf_xsmr                      => clm3%g%l%c%p%pc13f%leaf_xsmr
    litfall                        => clm3%g%l%c%p%pc13f%litfall
    livecrootc_xfer_to_livecrootc  => clm3%g%l%c%p%pc13f%livecrootc_xfer_to_livecrootc
    livecroot_mr                   => clm3%g%l%c%p%pc13f%livecroot_mr
    livecroot_curmr                => clm3%g%l%c%p%pc13f%livecroot_curmr
    livecroot_xsmr                 => clm3%g%l%c%p%pc13f%livecroot_xsmr
    livestemc_xfer_to_livestemc    => clm3%g%l%c%p%pc13f%livestemc_xfer_to_livestemc
    livestem_mr                    => clm3%g%l%c%p%pc13f%livestem_mr
    livestem_curmr                 => clm3%g%l%c%p%pc13f%livestem_curmr
    livestem_xsmr                  => clm3%g%l%c%p%pc13f%livestem_xsmr
    m_deadcrootc_storage_to_fire   => clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_fire
    m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_fire           => clm3%g%l%c%p%pc13f%m_deadcrootc_to_fire
    m_deadcrootc_to_litter         => clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter
    m_deadcrootc_to_litter_fire    => clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter_fire
    m_deadcrootc_xfer_to_fire      => clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_fire
    m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_fire    => clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_fire
    m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_litter
    m_deadstemc_to_fire            => clm3%g%l%c%p%pc13f%m_deadstemc_to_fire
    m_deadstemc_to_litter          => clm3%g%l%c%p%pc13f%m_deadstemc_to_litter
    m_deadstemc_to_litter_fire     => clm3%g%l%c%p%pc13f%m_deadstemc_to_litter_fire
    m_deadstemc_xfer_to_fire       => clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_fire
    m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_fire       => clm3%g%l%c%p%pc13f%m_frootc_storage_to_fire
    m_frootc_storage_to_litter     => clm3%g%l%c%p%pc13f%m_frootc_storage_to_litter
    m_frootc_to_fire               => clm3%g%l%c%p%pc13f%m_frootc_to_fire
    m_frootc_to_litter             => clm3%g%l%c%p%pc13f%m_frootc_to_litter
    m_frootc_xfer_to_fire          => clm3%g%l%c%p%pc13f%m_frootc_xfer_to_fire
    m_frootc_xfer_to_litter        => clm3%g%l%c%p%pc13f%m_frootc_xfer_to_litter
    m_gresp_storage_to_fire        => clm3%g%l%c%p%pc13f%m_gresp_storage_to_fire
    m_gresp_storage_to_litter      => clm3%g%l%c%p%pc13f%m_gresp_storage_to_litter
    m_gresp_xfer_to_fire           => clm3%g%l%c%p%pc13f%m_gresp_xfer_to_fire
    m_gresp_xfer_to_litter         => clm3%g%l%c%p%pc13f%m_gresp_xfer_to_litter
    m_leafc_storage_to_fire        => clm3%g%l%c%p%pc13f%m_leafc_storage_to_fire
    m_leafc_storage_to_litter      => clm3%g%l%c%p%pc13f%m_leafc_storage_to_litter
    m_leafc_to_fire                => clm3%g%l%c%p%pc13f%m_leafc_to_fire
    m_leafc_to_litter              => clm3%g%l%c%p%pc13f%m_leafc_to_litter
    m_leafc_xfer_to_fire           => clm3%g%l%c%p%pc13f%m_leafc_xfer_to_fire
    m_leafc_xfer_to_litter         => clm3%g%l%c%p%pc13f%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_fire   => clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_fire
    m_livecrootc_storage_to_litter => clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_litter
    m_livecrootc_to_fire           => clm3%g%l%c%p%pc13f%m_livecrootc_to_fire
    m_livecrootc_to_litter         => clm3%g%l%c%p%pc13f%m_livecrootc_to_litter
    m_livecrootc_xfer_to_fire      => clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_fire
    m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_fire    => clm3%g%l%c%p%pc13f%m_livestemc_storage_to_fire
    m_livestemc_storage_to_litter  => clm3%g%l%c%p%pc13f%m_livestemc_storage_to_litter
    m_livestemc_to_fire            => clm3%g%l%c%p%pc13f%m_livestemc_to_fire
    m_livestemc_to_litter          => clm3%g%l%c%p%pc13f%m_livestemc_to_litter
    m_livestemc_xfer_to_fire       => clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_fire
    m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_litter
    mr                             => clm3%g%l%c%p%pc13f%mr
    npp                            => clm3%g%l%c%p%pc13f%npp
    pft_fire_closs                 => clm3%g%l%c%p%pc13f%pft_fire_closs
    psnshade_to_cpool              => clm3%g%l%c%p%pc13f%psnshade_to_cpool
    psnsun_to_cpool                => clm3%g%l%c%p%pc13f%psnsun_to_cpool
    rr                             => clm3%g%l%c%p%pc13f%rr
    storage_gr                     => clm3%g%l%c%p%pc13f%storage_gr
    transfer_deadcroot_gr          => clm3%g%l%c%p%pc13f%transfer_deadcroot_gr
    transfer_deadstem_gr           => clm3%g%l%c%p%pc13f%transfer_deadstem_gr
    transfer_froot_gr              => clm3%g%l%c%p%pc13f%transfer_froot_gr
    transfer_gr                    => clm3%g%l%c%p%pc13f%transfer_gr
    transfer_leaf_gr               => clm3%g%l%c%p%pc13f%transfer_leaf_gr
    transfer_livecroot_gr          => clm3%g%l%c%p%pc13f%transfer_livecroot_gr
    transfer_livestem_gr           => clm3%g%l%c%p%pc13f%transfer_livestem_gr
    vegfire                        => clm3%g%l%c%p%pc13f%vegfire
    cpool                          => clm3%g%l%c%p%pc13s%cpool
    xsmrpool                       => clm3%g%l%c%p%pc13s%xsmrpool
    deadcrootc                     => clm3%g%l%c%p%pc13s%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pc13s%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pc13s%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pc13s%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pc13s%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pc13s%deadstemc_xfer
    dispvegc                       => clm3%g%l%c%p%pc13s%dispvegc
    frootc                         => clm3%g%l%c%p%pc13s%frootc
    frootc_storage                 => clm3%g%l%c%p%pc13s%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pc13s%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pc13s%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pc13s%gresp_xfer
    leafc                          => clm3%g%l%c%p%pc13s%leafc
    leafc_storage                  => clm3%g%l%c%p%pc13s%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pc13s%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pc13s%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pc13s%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pc13s%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pc13s%livestemc
    livestemc_storage              => clm3%g%l%c%p%pc13s%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pc13s%livestemc_xfer
    storvegc                       => clm3%g%l%c%p%pc13s%storvegc
    totpftc                        => clm3%g%l%c%p%pc13s%totpftc
    totvegc                        => clm3%g%l%c%p%pc13s%totvegc
    tempsum_npp                    => clm3%g%l%c%p%pepv%tempsum_npp

   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate pft-level summary carbon fluxes and states

      ! gross primary production (GPP)
      gpp(p) = &
         psnsun_to_cpool(p) + &
         psnshade_to_cpool(p)

      ! maintenance respiration (MR)
      
      leaf_mr(p)      = leaf_curmr(p)      + leaf_xsmr(p)
      froot_mr(p)     = froot_curmr(p)     + froot_xsmr(p)
      livestem_mr(p)  = livestem_curmr(p)  + livestem_xsmr(p)
      livecroot_mr(p) = livecroot_curmr(p) + livecroot_xsmr(p)
      
      mr(p)  = &
         leaf_mr(p)     + &
         froot_mr(p)    + &
         livestem_mr(p) + &
         livecroot_mr(p)

      ! growth respiration (GR)
      ! current GR is respired this time step for new growth displayed in this timestep
      current_gr(p) = &
         cpool_leaf_gr(p)      + &
         cpool_froot_gr(p)     + &
         cpool_livestem_gr(p)  + &
         cpool_deadstem_gr(p)  + &
         cpool_livecroot_gr(p) + &
         cpool_deadcroot_gr(p)

      ! transfer GR is respired this time step for transfer growth displayed in this timestep
      transfer_gr(p) = &
         transfer_leaf_gr(p)      + &
         transfer_froot_gr(p)     + &
         transfer_livestem_gr(p)  + &
         transfer_deadstem_gr(p)  + &
         transfer_livecroot_gr(p) + &
         transfer_deadcroot_gr(p)

      ! storage GR is respired this time step for growth sent to storage for later display
      storage_gr(p) = &
         cpool_leaf_storage_gr(p)      + &
         cpool_froot_storage_gr(p)     + &
         cpool_livestem_storage_gr(p)  + &
         cpool_deadstem_storage_gr(p)  + &
         cpool_livecroot_storage_gr(p) + &
         cpool_deadcroot_storage_gr(p)

      ! GR is the sum of current + transfer + storage GR
      gr(p) = &
         current_gr(p)  + &
         transfer_gr(p) + &
         storage_gr(p)

      ! autotrophic respiration (AR)
      ar(p) = mr(p) + gr(p)

      ! root respiration (RR)
      rr(p) = &
         froot_mr(p) + &
         cpool_froot_gr(p) + &
         cpool_livecroot_gr(p) + &
         cpool_deadcroot_gr(p) + &
         transfer_froot_gr(p) + &
         transfer_livecroot_gr(p) + &
         transfer_deadcroot_gr(p) + &
         cpool_froot_storage_gr(p) + &
         cpool_livecroot_storage_gr(p) + &
         cpool_deadcroot_storage_gr(p)

      ! net primary production (NPP)
      npp(p) = gpp(p) - ar(p)

      ! update the annual NPP accumulator, for use in allocation code
      tempsum_npp(p) = tempsum_npp(p) + npp(p)

      ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
      ! This is supposed to correspond as closely as possible to
      ! field measurements of AGNPP, so it ignores the storage pools
      ! and only treats the fluxes into displayed pools.
      agnpp(p) = &
         cpool_to_leafc(p)                  + &
         leafc_xfer_to_leafc(p)             + &
         cpool_to_livestemc(p)              + &
         livestemc_xfer_to_livestemc(p)     + &
         cpool_to_deadstemc(p)              + &
         deadstemc_xfer_to_deadstemc(p)

     ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
      ! This is supposed to correspond as closely as possible to
      ! field measurements of BGNPP, so it ignores the storage pools
      ! and only treats the fluxes into displayed pools.
      bgnpp(p) = &
         cpool_to_frootc(p)                   + &
         frootc_xfer_to_frootc(p)             + &
         cpool_to_livecrootc(p)               + &
         livecrootc_xfer_to_livecrootc(p)     + &
         cpool_to_deadcrootc(p)               + &
         deadcrootc_xfer_to_deadcrootc(p)

      ! litterfall (LITFALL)
      litfall(p) = &
         leafc_to_litter(p)                 + &
         frootc_to_litter(p)                + &
         m_leafc_to_litter(p)               + &
         m_leafc_storage_to_litter(p)       + &
         m_leafc_xfer_to_litter(p)          + &
         m_frootc_to_litter(p)              + &
         m_frootc_storage_to_litter(p)      + &
         m_frootc_xfer_to_litter(p)         + &
         m_livestemc_to_litter(p)           + &
         m_livestemc_storage_to_litter(p)   + &
         m_livestemc_xfer_to_litter(p)      + &
         m_deadstemc_to_litter(p)           + &
         m_deadstemc_storage_to_litter(p)   + &
         m_deadstemc_xfer_to_litter(p)      + &
         m_livecrootc_to_litter(p)          + &
         m_livecrootc_storage_to_litter(p)  + &
         m_livecrootc_xfer_to_litter(p)     + &
         m_deadcrootc_to_litter(p)          + &
         m_deadcrootc_storage_to_litter(p)  + &
         m_deadcrootc_xfer_to_litter(p)     + &
         m_gresp_storage_to_litter(p)       + &
         m_gresp_xfer_to_litter(p)          + &
         m_deadstemc_to_litter_fire(p)      + &
         m_deadcrootc_to_litter_fire(p)

      ! pft-level fire losses (VEGFIRE)
      vegfire(p) = 0._r8

      ! pft-level carbon losses to fire
      pft_fire_closs(p) = &
         m_leafc_to_fire(p)                + &
         m_leafc_storage_to_fire(p)        + &
         m_leafc_xfer_to_fire(p)           + &
         m_frootc_to_fire(p)               + &
         m_frootc_storage_to_fire(p)       + &
         m_frootc_xfer_to_fire(p)          + &
         m_livestemc_to_fire(p)            + &
         m_livestemc_storage_to_fire(p)    + &
         m_livestemc_xfer_to_fire(p)       + &
         m_deadstemc_to_fire(p)            + &
         m_deadstemc_storage_to_fire(p)    + &
         m_deadstemc_xfer_to_fire(p)       + &
         m_livecrootc_to_fire(p)           + &
         m_livecrootc_storage_to_fire(p)   + &
         m_livecrootc_xfer_to_fire(p)      + &
         m_deadcrootc_to_fire(p)           + &
         m_deadcrootc_storage_to_fire(p)   + &
         m_deadcrootc_xfer_to_fire(p)      + &
         m_gresp_storage_to_fire(p)        + &
         m_gresp_xfer_to_fire(p)

      ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
      dispvegc(p) = &
         leafc(p)      + &
         frootc(p)     + &
         livestemc(p)  + &
         deadstemc(p)  + &
         livecrootc(p) + &
         deadcrootc(p)

      ! stored vegetation carbon, excluding cpool (STORVEGC)
      storvegc(p) = &
      	cpool(p)              + &
         leafc_storage(p)      + &
         frootc_storage(p)     + &
         livestemc_storage(p)  + &
         deadstemc_storage(p)  + &
         livecrootc_storage(p) + &
         deadcrootc_storage(p) + &
         leafc_xfer(p)         + &
         frootc_xfer(p)        + &
         livestemc_xfer(p)     + &
         deadstemc_xfer(p)     + &
         livecrootc_xfer(p)    + &
         deadcrootc_xfer(p)    + &
         gresp_storage(p)      + &
         gresp_xfer(p)

      ! total vegetation carbon, excluding cpool (TOTVEGC)
      totvegc(p) = dispvegc(p) + storvegc(p)

      ! total pft-level carbon, including cpool (TOTPFTC)
      totpftc(p) = totvegc(p) + xsmrpool(p)

   end do  ! end of pfts loop

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c(num_soilc, filter_soilc, gpp, col_gpp)
   call p2c(num_soilc, filter_soilc, ar, col_ar)
   call p2c(num_soilc, filter_soilc, rr, col_rr)
   call p2c(num_soilc, filter_soilc, npp, col_npp)
   call p2c(num_soilc, filter_soilc, vegfire, col_vegfire)
   call p2c(num_soilc, filter_soilc, totvegc, col_totvegc)
   call p2c(num_soilc, filter_soilc, totpftc, col_totpftc)
   call p2c(num_soilc, filter_soilc, pft_fire_closs, col_pft_fire_closs)

   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! litter heterotrophic respiration (LITHR)
      lithr(c) = &
         litr1_hr(c) + &
         litr2_hr(c) + &
         litr3_hr(c)

      ! soil organic matter heterotrophic respiration (SOMHR)
      somhr(c) = &
         soil1_hr(c) + &
         soil2_hr(c) + &
         soil3_hr(c) + &
         soil4_hr(c)

      ! total heterotrophic respiration (HR)
      hr(c) = lithr(c) + somhr(c)

      ! total soil respiration, heterotrophic + root respiration (SR)
      sr(c) = col_rr(c) + hr(c)

      ! total ecosystem respiration, autotrophic + heterotrophic (ER)
      er(c) = col_ar(c) + hr(c)

      ! litter fire losses (LITFIRE)
      litfire(c) = 0._r8

      ! soil organic matter fire losses (SOMFIRE)
      somfire(c) = 0._r8

      ! total ecosystem fire losses (TOTFIRE)
      totfire(c) = &
         litfire(c) + &
         somfire(c) + &
         col_vegfire(c)

      ! column-level carbon losses to fire, including pft losses
      col_fire_closs(c) = &
         m_litr1c_to_fire(c)  + &
         m_litr2c_to_fire(c)  + &
         m_litr3c_to_fire(c)  + &
         m_cwdc_to_fire(c)    + &
         col_pft_fire_closs(c)

      ! net ecosystem production, excludes fire flux, positive for sink (NEP)
      nep(c) = col_gpp(c) - er(c)

      ! net ecosystem exchange of carbon, includes fire flux, positive for source (NEE)
      nee(c) = -nep(c) + col_fire_closs(c)

      ! total litter carbon (TOTLITC)
      totlitc(c) = &
         litr1c(c) + &
         litr2c(c) + &
         litr3c(c)

      ! total soil organic matter carbon (TOTSOMC)
      totsomc(c) = &
         soil1c(c) + &
         soil2c(c) + &
         soil3c(c) + &
         soil4c(c)

      ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
      totecosysc(c) = &
         cwdc(c) + &
         totlitc(c) + &
         totsomc(c) + &
         col_totvegc(c)

      ! total column carbon, including veg and cpool (TOTCOLC)
      totcolc(c) = &
         cwdc(c) + &
         totlitc(c) + &
         totsomc(c) + &
         col_totpftc(c)

   end do ! end of columns loop


end subroutine C13Summary
!-----------------------------------------------------------------------

end module C13SummaryMod
