module mod_che_aerosols_sorgam

  use mod_intkinds
  use mod_realkinds
  use mod_che_molwg
  use mod_che_indices
  use mod_che_data_sorgam

  implicit none

  contains
  subroutine sorgam_driver (dtstep,t_phy,moist,p_phy,chem,chem_n,chem_a,    &
               nu0,ac0,nu3,ac3,vcsulf_old,drog,i,j,k,num_chem,cb6soa_cat,   &
               pdensn,pdensa,dgnuc,dgacc )

    integer(ik4), intent(in) :: num_chem
    real(rk8), dimension(num_chem), intent(inout) :: chem , chem_a , chem_n
    integer(ik4), dimension(num_chem), intent(in) :: cb6soa_cat
    real(rk8), intent(in) :: moist , vcsulf_old 
   
    ! following are aerosol arrays that are not advected
    real(rk8), intent(inout) :: nu0,ac0,nu3,ac3
    ! organic aerosol precursor DeltaROG [ug/m^3]
    real(rk8), dimension(num_chem), intent(in) :: drog 
    real(rk8), intent(in) :: dtstep , t_phy , p_phy

    ! (internal aerosol dynamics)
    real(rk8) :: cblk(num_chem)   ! main array of gas [ug/m^3]
    real(rk8) :: cblk_a(num_chem) ! for accumulation particles
    real(rk8) :: cblk_n(num_chem) ! for Aitken particles
    real(rk8) :: cblk_vnu3 , cblk_vac3 , cblk_vnu0 , cblk_vac0
    real(rk8), intent(out) :: pdensn ! ave particle density in nucl [kg/m^3]
    real(rk8), intent(out) :: pdensa ! ave particle density in accu [kg/m^3]
    real(rk8), intent(out) :: dgnuc  ! nuclei mode  geometric mean diam [m]
    real(rk8), intent(out) :: dgacc  ! accumulation geometric mean diam [m]
    real(rk8) :: so4rat              ! input SO4 formation[ug/m^3/s]
!   real(rk8) :: epmcoarse
!   REAL eorgi_in ! Emission rate of j-mode org. aerosol [ug m**-
!   REAL eorgj_in ! Emission rate of j-mode org. aerosol [ug m**-
    real(rk8) :: p, t, rh, v2
    integer(ik4) :: i,j,k,l

    t   = t_phy
    p   = p_phy
    rh  = max(0.1,moist)
    rh  = min(0.95,rh)
    so4rat = (chem(ind_SULF) - vcsulf_old)/dtstep
            
    cblk(:)    = chem(:)         ! gas phase data
    cblk_n(:)  = max(epsilc,chem_n(:)) ! accumulation data
    cblk_a(:)  = max(epsilc,chem_a(:)) ! Aitken data
    cblk_vnu3  = max(epsilc,nu3) !
    cblk_vac3  = max(epsilc,ac3) !
    cblk_vnu0  = max(1.0E07,nu0) ! 
    cblk_vac0  = max(1.0E07,ac0) ! 

    call rpmmod3(num_chem,dtstep,p,t,rh,cb6soa_cat,so4rat,drog,        &
           cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,cblk_vnu3,cblk_vac3, &
           i,j,k,pdensn,pdensa,dgnuc,dgacc)

    chem(:)   = cblk(:)               ! gas phase data
    chem_a(:) = max(epsilc,cblk_a(:)) ! accumulation data
    chem_n(:) = max(epsilc,cblk_n(:)) ! Aitken data
    nu0 = max(1.0E07,cblk_vnu0)
    ac0 = max(1.0E07,cblk_vac0)
    nu3 = max(epsilc,cblk_vnu3)
    ac3 = max(epsilc,cblk_vac3)
   
    do l = 1, num_chem
      if ( cb6soa_cat(l) /= 0 .and. cb6soa_cat(l) /= 99) then
        chem(l)  =max(epsilc,(chem(l)-chem_a(l)-chem_n(l)))
      end if
    end do

  end subroutine sorgam_driver

  !***********************************************************************

  !       BEGIN OF AEROSOL CALCULATIONS

  !ia*********************************************************************
  !ia                                                                     *
  !ia     MAIN AEROSOL DYNAMICS ROUTINE                                   *
  !ia     based on MODELS3 formulation by FZB                             *
  !ia     Modified by IA in May 97                                        *
  !ia     THIS PROGRAMME IS THE LINK BETWEEN GAS PHASE AND AEROSOL PHASE
  !ia     CALCULATIONS IN THE COLUMN MODEL. IT CONVERTS ALL DATA AND
  !ia     VARIABLES BETWEEN THE TWO PARTS AND DRIVES THE AEROSOL
  !ia     CALCULATIONS.
  !ia     INPUT DATA REQUIRED FOR AEROSOL DYNAMICS ARE SET UP HERE FOR
  !ia     ONE GRID CELL!!!!
  !ia     and passed to dynamics calcs. subroutines.
  !ia                                                                     *
  !ia     Revision history                                                *
  !ia     When    WHO     WHAT                                            *
  !ia     ----    ----    ----                                            *
  !ia     ????    FZB     BEGIN                                           *
  !ia     05/97   IA      Adapted for use in CTM2-S                       *
  !ia                     Modified renaming/bug fixing                    *
  !ia     11/97   IA      Modified for new model version
  !ia                     see comments under iarev02
  !ia     03/98   IA      corrected error on pressure units
  !ia                                                                     *
  !ia     Called BY:      CHEM                                            *
  !ia                                                                     *
  !ia     Calls to:       OUTPUT1,AEROPRC                                 *
  !ia                                                                     *
  !ia*********************************************************************

    subroutine rpmmod3(nspcsda,dtsec,pres,temp,relhum,cb6soa_cat, &
        so4rat,drog,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,       &
        cblk_vnu3,cblk_vac3,igrid,jgrid,kgrid,pdensn,pdensa,dgnuc,dgacc)

      !implicit none
      ! *** inputs
      real(rk8) :: temp      ! Input temperature [K]
      real(rk8) :: relhum    ! Input relative humidity  [fraction]
      real(rk8) :: pres      ! Input pressure [Pa]
      real(rk8) :: numnuc_in ! Input number for Aitken mode [/m^3]
      real(rk8) :: numacc_in ! Input number for accumulation mode [/m^3]
      ! Production rate of soil derived coarse
!     real(rk8) :: eorgi_in   ! Emission rate of i-mode org. aerosol [ug m**-
!     real(rk8) :: eorgj_in   ! Emission rate of j-mode org. aerosol [ug m**-
      real(rk8) :: drog(nspcsda) ! organic aerosol precursor [ug/m^3]
      ! *** Primary emissions rates: [ug/m^3/s]
      ! *** emissions rates for primary organic aerosol
!     real(rk8) :: eorgi  ! Aitken mode
!     real(rk8) :: eorgj  ! Accumululaton mode
      real(rk8) :: dtsec  ! time step [s], PASSED FROM MAIN COLUMN MODE
      real(rk8) :: step   ! synchronization time  [s]
      ! *** arrays for aerosol model codes:
      integer(ik4) :: nspcsda ! number of species in CBLK ciarev02
      integer(ik4), dimension(nspcsda) :: cb6soa_cat
      real(rk8) :: cblk(nspcsda) , cblk_a(nspcsda) , cblk_n(nspcsda)
      real(rk8) :: cblk_vnu0 , cblk_vac0 , cblk_vnu3 , cblk_vac3
      ! *** Meteorological information in blocked arays:
      ! *** Thermodynamic variables:
      ! main array of variables
      real(rk8) :: blkta   ! Air temperature [K]
      real(rk8) :: blkprs  ! Air pressure in [Pa]
      real(rk8) :: blkdens ! Air density  [kg/m^3]
      real(rk8) :: blkrh   ! Fractional relative humidity
      ! *** Chemical production rates [ug/m^-3/s] :
      real(rk8) :: so4rat          ! sulfuric acid vapor-phase production
      real(rk8) :: orgrat(nspcsda) ! organic gas-phase production rate
      ! *** atmospheric properties
      real(rk8) :: xlm ! atmospheric mean free path [m]
      real(rk8) :: amu ! atmospheric dynamic viscosity [kg/m/s]
      ! *** aerosol properties:
      ! *** modal diameters:
      real(rk8) :: dgnuc  ! nuclei mode  geometric mean diam [m]
      real(rk8) :: dgacc  ! accumulation geometric mean diam [m]
      ! *** Modal mass concentrations [ug/m^3] 
      real(rk8) :: pmassn ! mass concentration in Aitken mode
      real(rk8) :: pmassa ! mass concentration in accumulation
      ! *** average modal particle densities  [kg/m^3]
      real(rk8) :: pdensn ! average particle density in nuclei
      real(rk8) :: pdensa ! average particle density in accumu
      ! *** average modal Knudsen numbers
      real(rk8) :: knnuc  ! nuclei mode  Knudsen number
      real(rk8) :: knacc  ! accumulation Knudsen number
     ! *** reciprocal modal condensation rates for sulfuric acid 
      real(rk8) :: fconcn ! reciprocal condensation rate Aitke
      real(rk8) :: fconca ! reciprocal condensation rate acclu
      real(rk8) :: fconcn_org , fconca_org
      ! *** Rates for secondary particle formation:
      ! *** production of new mass concentration [ug/m^3/s]
      real(rk8) :: dmdt   ! by particle formation
      ! *** production of new number concentration [#/m^3/s]
      ! rate of production of new mass concen
      real(rk8) :: dndt   ! by particle formation
      ! *** growth rate for third moment by condensation of precursor
      !      vapor on existing particles [ 3rd mom/m^3/s ]
      ! rate of producton of new particle num
      real(rk8) :: cgrn3  !  Aitken mode [1/s]
      real(rk8) :: cgra3  !  Accumu mode [1/s]
      ! *** Rates for coaglulation: [m^3/s]
      ! *** Unimodal Rates:
      real(rk8) :: urn00  ! Aitken mode 0th moment self-coagulation ra
      real(rk8) :: ura00  ! accumulation mode 0th moment self-coagulat
      ! *** Bimodal Rates:  Aitken mode with accumulation mode ( d( Aitken mod
      real(rk8) :: brna01 ! rate for 0th moment
      real(rk8) :: brna31 ! rate for 3rd moment
      ! *** other processes
      real(rk8) :: deltaso4a ! sulfate aerosol by condensation [ug/m^3]
      integer(ik4) :: l,igrid,jgrid,kgrid

      ! *** set up experimental conditions
      ! *** initialize model variables
      step    = dtsec     ! set time step
      blkta   = temp      ! T in Kelvin
      blkprs  = pres      ! P in  Pa (pres is given in
      blkrh   = relhum    ! fractional RH
      blkdens = blkprs/(rdgas*blkta)

      !...INITIALISE EMISSION RATES
!     eorgi = eorgi_in   ! primary organic
!     eorgj = eorgj_in   ! primary organic
      dgnuc = dginin
      dgacc = dginia

      ! ***  time loop
      call aeroproc(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,            &
        cblk_vnu3,cblk_vac3,step,blkta,blkprs,blkdens,blkrh,so4rat,            &
        orgrat,drog,xlm,amu,dgnuc,dgacc,pmassn,pmassa,pdensn,pdensa,knnuc,     &
        knacc,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3,urn00, &
        ura00,brna01,brna31,deltaso4a,igrid,jgrid,kgrid,cb6soa_cat)
      return

    end subroutine rpmmod3

    ! *********************************************************************** !

    subroutine aeroproc(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,        &
        cblk_vnu3,cblk_vac3,dt,blkta,blkprs,blkdens,blkrh,so4rat,              &
        orgrat,drog,xlm,amu,dgnuc,dgacc,pmassn,pmassa,pdensn,pdensa,knnuc,     &
        knacc,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,cgrn3,cgra3,urn00, &
        ura00,brna01,c30,deltaso4a,igrid,jgrid,kgrid,cb6soa_cat)

      !implicit none
      integer(ik4) :: nspcsda  ! number of species in CBLK
      integer(ik4), dimension(nspcsda) :: cb6soa_cat
      real(rk8) :: cblk(nspcsda) ! main array of variables (INPUT a
      real(rk8) :: cblk_n(nspcsda) , cblk_a(nspcsda)
      real(rk8) :: cblk_vnu0 , cblk_vac0 , cblk_vnu3 , cblk_vac3
      real(rk8) :: dt ! synchronization time  [s]
      ! *** Meteorological information:
      real(rk8) :: blkta   ! Air temperature [K]
      real(rk8) :: blkprs  ! Air pressure in [Pa]
      real(rk8) :: blkdens ! Air density  [kg/m^3]
      real(rk8) :: blkrh   ! Fractional relative humidity
      ! *** Chemical production rates: [ug/m^3/s]
      real(rk8) :: so4rat  ! sulfate gas-phase production rate
      ! organic condensable vapor production rate
      real(rk8) :: drog(nspcsda) ! Delta ROG conc. [ug/m3]
      ! *** organic aerosol mass production rates 
      real(rk8), dimension(nspcsda) :: orgrat
      ! *** Primary emissions rates: [ ug / m**3 s ]
      ! *** emissions rates for primary organic aerosol
!     real(rk8) :: eorgi     ! Aitken mode
!     real(rk8) :: eorgj     ! Accumululaton mode

      ! *** OUTPUT:
      ! *** atmospheric properties
      real(rk8) :: xlm ! atmospheric mean free path [m]
      real(rk8) :: amu ! atmospheric dynamic viscosity [kg/m/s]
      ! *** modal diameters: [m]
      real(rk8) :: dgnuc ! nuclei mode geometric mean diamete
      real(rk8) :: dgacc ! accumulation geometric mean diamet
      ! *** aerosol properties:
      ! *** Modal mass concentrations [ug/m^3]
      real(rk8) :: pmassn ! mass concentration in Aitken mode
      real(rk8) :: pmassa ! mass concentration in accumulation
      ! *** average modal particle densities  [kg/m^3]
      real(rk8) :: pdensn ! average particle density in nuclei
      real(rk8) :: pdensa ! average particle density in accumu
      ! *** average modal Knudsen numbers
      real(rk8) :: knnuc ! nuclei mode  Knudsen number
      real(rk8) :: knacc ! accumulation Knudsen number
      ! ***  modal condensation factors ( see comments in NUCLCOND )
      real(rk8) :: fconcn , fconca , fconcn_org , fconca_org
      ! *** Rates for secondary particle formation:
      ! *** production of new mass concentration [ug/m^3/s]
      real(rk8) :: dmdt ! by particle formation
      ! *** production of new number concentration [#/m^3/s]
      ! rate of production of new mass concen
      real(rk8) :: dndt ! by particle formation
      ! *** growth rate for third moment by condensation of precursor
      !      vapor on existing particles [3rd mom/m^3/s]
      ! rate of producton of new particle num
      real(rk8) :: cgrn3 !  Aitken mode [1/s]
      real(rk8) :: cgra3 !  Accumu mode [1/s]
      ! *** Rates for coaglulation: [m^3/s]
      ! *** Unimodal Rates:
      real(rk8) :: urn00 ! Aitken mode 0th moment self-coagulation ra
      real(rk8) :: ura00 ! accumulation mode 0th moment self-coagulat
      ! *** Bimodal Rates:  Aitken mode with accumulation mode ( d( Aitken mod
      real(rk8) :: brna01
      ! *** 3rd moment intermodal transfer rate replaces coagulation rate ( FS
      ! rate for 0th moment
      real(rk8) :: c30   ! by intermodal c
      ! *** other processes
      ! intermodal 3rd moment transfer r
      real(rk8) :: deltaso4a ! sulfate aerosol by condensation [ug/m^3]
      integer(ik4) :: igrid,jgrid,kgrid

      ! *** get size distribution information:
      call modpar(nspcsda,cblk_n,cblk_a,blkta,blkprs,pmassn,pmassa,           &
        pdensn,pdensa,xlm,amu,dgnuc,dgacc,knnuc,knacc,cblk_vnu0,cblk_vac0,    &
        cblk_vnu3,cblk_vac3,cb6soa_cat)

      ! *** Calculate coagulation rates for fine particles:
      call coagrate(nspcsda,cblk_vnu0,cblk_vac0,blkta,pdensn,                 &
        pdensa,amu,dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)

      ! *** get condensation and particle formation (nucleation) rates:
      call nuclcond(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,           &
        cblk_vnu3,cblk_vac3,dt,blkta,blkprs,blkrh,so4rat,orgrat,              &
        drog,dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt,       &
        deltaso4a,cgrn3,cgra3,cb6soa_cat)

      ! *** advance forward in time  DT seconds:
      call aerostep(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,           &
        cblk_vnu3,cblk_vac3,dt,so4rat,orgrat,dgnuc,dgacc,fconcn,fconca,       &
        fconcn_org,fconca_org,pmassn,pmassa,dmdt,dndt,deltaso4a,              &
        urn00,ura00,brna01,c30,cgrn3,cgra3,igrid,jgrid,kgrid,cb6soa_cat)

      ! *** get new distribution information:
      call modpar(nspcsda,cblk_n,cblk_a,blkta,blkprs,pmassn,pmassa,           &
        pdensn,pdensa,xlm,amu,dgnuc,dgacc,knnuc,knacc,cblk_vnu0,cblk_vac0,    &
        cblk_vnu3,cblk_vac3,cb6soa_cat)

      return

    end subroutine aeroproc

    ! ********************************************************************** !

    subroutine modpar(nspcsda,cblk_n,cblk_a,blkta,blkprs,pmassn,pmassa,       &
        pdensn,pdensa,xlm,amu,dgnuc,dgacc,knnuc,knacc,cblk_vnu0,cblk_vac0,    &
        cblk_vnu3,cblk_vac3,cb6soa_cat)

      !***********************************************************************
      !**    DESCRIPTION:
      !       Calculates modal parameters and derived variables,
      !       log-squared of std deviation, mode mean size, Knudsen number)
      !       based on current values of moments for the modes.
      ! FSB   Now calculates the 3rd moment, mass, and density in all 3 modes.
      !**
      !**    Revision history:
      !       Adapted 3/95 by US and CJC from EAM2's MODPAR and INIT3
      !       Revised  7/23/96 by FSB to use COMMON blocks and small blocks
      !        instead of large 3-d arrays, and to assume a fixed std.
      !       Revised 12/06/96 by FSB to include coarse mode
      !       Revised 1/10/97 by FSB to have arrays passed in call vector
      !**********************************************************************

      ! *** input:
      integer(ik4) :: nspcsda  ! nmber of species in CBLK
      integer(ik4), dimension(nspcsda) :: cb6soa_cat
      real(rk8) :: cblk_n(nspcsda) 
      real(rk8) :: cblk_a(nspcsda)
      real(rk8) :: cblk_vnu3    ! 3rd moment Aitken particles
      real(rk8) :: cblk_vac3    ! 3rd moment accumulation particles
      real(rk8) :: cblk_vnu0    ! 0th moment Aitken particles
      real(rk8) :: cblk_vac0    ! 0th moment accumulation particles
      real(rk8) :: blkta        ! Air temperature [K]
      real(rk8) :: blkprs       ! Air pressure in [Pa]

      ! *** output:
      real(rk8) :: pmassn ! mass concentration in nuclei mode  [ug/m^3]
      real(rk8) :: pmassa ! mass concentration in accumulation [ug/m^3]
      real(rk8) :: pdensn ! average particel density in Aitken [kg/m^3]
      real(rk8) :: pdensa ! average particel density in accumu [kg/m^3]
      real(rk8) :: xlm    ! atmospheric mean free path [m]
      real(rk8) :: amu    ! atmospheric dynamic viscosity [kg/m/s]
      real(rk8) :: dgnuc  ! Aitken mode mean diameter [m]
      real(rk8) :: dgacc  ! accumu mode mean diameter [m]
      real(rk8) :: knnuc  ! Aitken mode Knudsen number
      real(rk8) :: knacc  ! accumulation

      integer (ik4) :: l

      ! *** set up  aerosol 0th moment, 3rd moment, mass, density
      pmassn = 0.
      pmassa = 0.
      cblk_vnu3 = 0.
      cblk_vac3 = 0.
      do l = 1, nspcsda
        if ( cb6soa_cat(l) /= 0 .and. cb6soa_cat(l) /= 99 ) then
          aerfact = orgfac 
          if (cb6soa_cat(l) == 2 ) aerfact = f6dpim9/rhocnyl !carbonyl
          if (cb6soa_cat(l) == 3 ) aerfact = f6dpim9/rhosthc !alkanes + ols
          if (cb6soa_cat(l) == 4 ) aerfact = f6dpim9/rhoolef !olefins
          if (cb6soa_cat(l) == 6 ) aerfact = f6dpim9/rhoarom !aromatics
          if (cb6soa_cat(l) == 9 ) aerfact = f6dpim9/rhontro !organic nitrates
          if (cb6soa_cat(l) == 10) aerfact = f6dpim9/rhomfnc !multifunctional
          if (cb6soa_cat(l) == 22) aerfact = so4fac
          ! *** set up aerosol 3rd moment, mass, density
          cblk_vnu3 = cblk_vnu3 + aerfact*cblk_n(l) ! *** Aitken-mode
          cblk_vac3 = cblk_vac3 + aerfact*cblk_a(l) ! *** Accumu-mode
          ! *** now prepare particle mass and density
          pmassn = pmassn + cblk_n(l) ! *** Aitken-mode
          pmassa = pmassa + cblk_a(l) ! *** Accumulation-mode
        end if
      end do

      ! *** now conclude particle mass and density
      pmassn = max(conmin, pmassn) ! *** Aitken-mode
      pmassa = max(conmin, pmassa) ! *** Accumulation-mode
      cblk_vnu3 = max(conmin,cblk_vnu3)
      cblk_vac3 = max(conmin,cblk_vac3)
      ! *** set up aerosol 0th moment
      cblk_vnu0 = cblk_vnu3/((dginin**3)*esn36)
      cblk_vac0 = cblk_vac3/((dginia**3)*esa36)
      cblk_vnu0 = max(1.0E07,cblk_vnu0)
      cblk_vac0 = max(1.0E07,cblk_vac0)

      ! *** now get particle density, mean free path, and dynamic viscosity
      ! *** density in [ kg m**-3 ]
      pdensn = max(densmin,(f6dpim9*pmassn/cblk_vnu3))
      pdensa = max(densmin,(f6dpim9*pmassa/cblk_vac3))
      ! *** Calculate mean free path [ m ]:
      xlm = 6.6328E-8*pss0*blkta/(tss0*blkprs)
      ! *** 6.6328E-8 is the sea level values given in Table I.2.8
      ! *** on page 10 of U.S. Standard Atmosphere 1962

      ! *** Calculate dynamic viscosity [ kg m**-1 s**-1 ]:
      ! *** U.S. Standard Atmosphere 1962 page 14 expression
      !     for dynamic viscosity is:
      !     dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
      !     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].
      amu = 1.458E-6*blkta*sqrt(blkta)/(blkta+110.4)

      !...............   Standard deviation fixed in both modes, so
      !...............   diagnose diameter from 3rd moment and number concentr

      ! density and mean free path
      ! calculate diameters
      dgnuc = max(dgminn,(cblk_vnu3/(cblk_vnu0*esn36))** one3)
      dgacc = max(dgmina,(cblk_vac3/(cblk_vac0*esa36))** one3)
      ! when running with cloudborne aerosol, apply some very mild bounding
      ! to avoid unrealistic dg values
      !if (cw_phase > 0) then
      !  dgnuc = max( dgnuc, dginin*0.2  )  !  > 0.002 um
      !  dgnuc = min( dgnuc, dginin*10.0 )  !  < 0.10  um
      !  dgacc = max( dgacc, dginia*0.2  )  !  > 0.014 um
      !  dgacc = min( dgacc, dginia*10.0 )  !  < 0.7 um
      !end if
      
      ! Calculate Knudsen numbers
      knnuc = 2.0*xlm/dgnuc
      knacc = 2.0*xlm/dgacc
      return

    end subroutine modpar

    ! *********************************************************************** !

    subroutine coagrate(nspcsda,cblk_vnu0,cblk_vac0,blkta,pdensn,  &
        pdensa,amu,dgnuc,dgacc,knnuc,knacc,urn00,ura00,brna01,c30)
      !***********************************************************************
      !**    DESCRIPTION:  calculates aerosol coagulation rates for unimodal
      !       and bimodal coagulation using E. Whitby 1990's prescription.
      !
      !.......   Rates for coaglulation:
      !.......   Unimodal Rates:
      !.......   URN00:  nuclei       mode 0th moment self-coagulation rate
      !.......   URA00:  accumulation mode 0th moment self-coagulation rate
      !
      !.......   Bimodal Rates:  (only 1st order coeffs appear)
      !.......   NA-- nuclei  with accumulation coagulation rates,
      !.......   AN-- accumulation with nuclei coagulation rates
      !.......   BRNA01:  rate for 0th moment ( d(nuclei mode 0) / dt  term)
      !.......   BRNA31:           3rd        ( d(nuclei mode 3) / dt  term)
      !**
      !**
      !**    Revision history:
      !       prototype 1/95 by Uma and Carlie
      !       Revised   8/95 by US for calculation of density from stmt func
      !                 and collect met variable stmt funcs in one include fil
      !      REVISED 7/25/96 by FSB to use block structure
      !      REVISED 9/13/96 BY FSB for Uma's FIXEDBOTH case only.
      !      REVISED 11/08/96 BY FSB the Whitby Shankar convention on signs
      !                              changed. All coagulation coefficients
      !                              returned with positive signs. Their
      !                              linearization is also abandoned.
      !                              Fixed values are used for the corrections
      !                              to the free-molecular coagulation integra
      !                              The code forces the harmonic means to be
      !                              evaluated in 64 bit arithmetic on 32 bit
      !     REVISED 11/14/96 BY FSB  Internal units are now MKS, moment / unit
      !
      !      REVISED 1/12/98 by FSB   C30 replaces BRNA31 as an array. This wa
      !                              because BRNA31 can become zero on a works
      !                              because of limited precision. With the ch
      !                              aerostep to omit update of the 3rd moment
      !                              C30 is the only variable now needed.
      !                              the logic using ONE88 to force REAL*8 ari
      !                              has been removed and all intermediates ar
      !                              REAL*8.

      !implicit none
      integer(ik4) :: nspcsda  ! nmber of species in CBLK

      real(rk8) :: cblk_vnu0    ! 0th moment Aitken particles
      real(rk8) :: cblk_vac0    ! 0th moment accumu particles
      real(rk8) :: blkta        ! Air temperature [K]
      real(rk8) :: pdensn       ! average particel density in Aitk [kg/m^3]
      real(rk8) :: pdensa       ! average particel density in accu [kg/m^3]
      real(rk8) :: amu          ! atmospheric dynamic viscosity [kg/m/s]
      real(rk8) :: dgnuc        ! Aitken mode mean diameter [m]
      real(rk8) :: dgacc        ! accumu mode mean diameter [m]
      real(rk8) :: knnuc        ! Aitken mode Knudsen number
      real(rk8) :: knacc        ! accumu mode Knudsen number

      ! *** output:
      real(rk8) :: urn00  ! intramodal coagulation rate (Ait
      real(rk8) :: ura00  ! intramodal coagulation rate (acc
      real(rk8) :: brna01 ! intermodal coagulaton rate (numb
      real(rk8) :: c30    ! intermodal 3rd moment transfer r

      ! *** Local variables:
      real(rk8) :: kncnuc, kncacc ! coeffs for unimodal NC coag rate
      real(rk8) :: kfmnuc, kfmacc ! coeffs for unimodal FM coag rate
      real(rk8) :: knc, kfm       ! coeffs for bimodal NC, FM coag rate
      real(rk8) :: bencnn, bencna ! NC 0th moment coag rate (both modes)
      real(rk8) :: bencm3n        ! NC 3rd moment coag rate (nuc mode)
      real(rk8) :: befmnn, befmna ! FM 0th moment coag rate (both modes)
      real(rk8) :: befm3n         ! FM 3rd moment coag rate (nuc mode)
      real(rk8) :: betann, betana ! composite coag rates, mom 0 (both mode
      real(rk8) :: brna31         ! intermodal coagulation rate for 3rd mo
      real(rk8) :: s1             ! scratch subexpression
      real(rk8) :: t1, t2         ! scratch subexpressions
      real(rk8) :: t16, t26       ! T1**6, T2**6
      real(rk8) :: rat, rin       ! ratio of acc to nuc size and its inver
      real(rk8) :: rsqt, rsq4     ! sqrt( rat ), rsqt**4
      real(rk8) :: rsqti, rsqi3   ! sqrt( 1/rat ), sqrt( 1/rat**3 )
      real(rk8) :: dgn3           ! dgnuc**3 [m^3]
      real(rk8) :: dga3           ! dgacc**3 [m^3]

      ! *** Fixed values for correctionss to coagulation
      !      integrals for free-molecular case.
      real(rk8), parameter :: bm0  = 0.8D0
      real(rk8), parameter :: bm0i = 0.9D0
      real(rk8), parameter :: bm3i = 0.9D0
      real(rk8), parameter :: a  = 1.246D0 ! approx Cunningham corr. factor

      !...........   Main computational grid-traversal loops
      !...........   for computing coagulation rates.
      ! *** Both modes have fixed std devs.
      ! *** moment independent factors
      s1 = two3*boltz*blkta/amu
      ! For unimodal coagualtion:
      kncnuc = s1
      kncacc = s1

      kfmnuc = sqrt(3.0*boltz*blkta/pdensn)
      kfmacc = sqrt(3.0*boltz*blkta/pdensa)
      ! For bimodal coagulation:
      knc = s1
      kfm = sqrt(6.0*boltz*blkta/(pdensn+pdensa))
      !...........   Begin unimodal coagulation rate calculations:
      !...........   Near-continuum regime.
      dgn3 = dgnuc**3
      dga3 = dgacc**3

      t1 = sqrt(dgnuc)
      t2 = sqrt(dgacc)
      t16 = dgn3 ! = T1**6
      t26 = dga3 ! = T2**6
      !.......   Note rationalization of fractions and subsequent cancellation
      !.......   from the formulation in  Whitby et al. (1990)
      bencnn = kncnuc*(1.0+esn08+a*knnuc*(esn04+esn20))
      bencna = kncacc*(1.0+esa08+a*knacc*(esa04+esa20))
      !...........   Free molecular regime. Uses fixed value for correction
      !               factor BM0
      befmnn = kfmnuc*t1*(en1+esn25+2.0*esn05)*bm0
      befmna = kfmacc*t2*(ea1+esa25+2.0*esa05)*bm0
      !...........   Calculate half the harmonic mean between unimodal rates
      !...........   free molecular and near-continuum regimes
      ! FSB       64 bit evaluation
      betann = bencnn*befmnn/(bencnn+befmnn)
      betana = bencna*befmna/(bencna+befmna)

      urn00 = betann
      ura00 = betana
      ! *** End of unimodal coagulation calculations.
      !...........   Begin bimodal coagulation rate calculations:
      rat = dgacc/dgnuc
      rin = 1.0D0/rat
      rsqt = sqrt(rat)
      rsq4 = rat**2

      rsqti = 1.0D0/rsqt
      rsqi3 = rin*rsqti
      !...........   Near-continuum coeffs:
      !...........   0th moment nuc mode bimodal coag coefficient
      bencnn = knc*(2.0+a*knnuc*(esn04+rat*esn16*esa04)+a*knacc* &
        (esa04+rin*esa16*esn04)+(rat+rin)*esn04*esa04)
      !...........   3rd moment nuc mode bimodal coag coefficient
      bencm3n = knc*dgn3*(2.0*esn36+a*knnuc*(esn16+rat*esn04*esa04)+a &
        *knacc*(esn36*esa04+rin*esn64*esa16)+rat*esn16*esa04+ &
        rin*esn64*esa04)
      !...........   Free molecular regime coefficients:
      !...........   Uses fixed value for correction
      !               factor BM0I, BM3I
      !...........   0th moment nuc mode coeff
      befmnn = kfm*bm0i*t1*(en1+rsqt*ea1+2.0*rat*en1*esa04+rsq4*esn09*esa16+ &
        rsqi3*esn16*esa09+2.0*rsqti*esn04*ea1)
      !...........   3rd moment nuc mode coeff
      befm3n = kfm*bm3i*t1*t16*(esn49+rsqt*esn36*ea1+2.0*rat*esn25*esa04+ &
        rsq4*esn09*esa16+rsqi3*esn100*esa09+2.0*rsqti*esn64*ea1)
      !...........   Calculate half the harmonic mean between bimodal rates
      !...........   free molecular and near-continuum regimes
      ! FSB       Force 64 bit evaluation
      brna01 = bencnn*befmnn/(bencnn+befmnn)

      brna31 = bencm3n*befm3n/(bencm3n+befm3n) ! BRNA31 now is a scala
      c30 = brna31*cblk_vac0*cblk_vnu0
      ! 3d moment transfer by intermodal coagula
      !         End bimodal coagulation rate.
      return

    end subroutine coagrate
! subroutine  to find the roots of a cubic equation / 3rd order polynomi
! formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   dr. francis binkowski modified the routine on 6/24/91, 8/7/97
! ***

    subroutine nuclcond(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0, &
        cblk_vnu3,cblk_vac3,dt,blkta,blkprs,blkrh,so4rat,orgrat,        &
        drog,dgnuc,dgacc,fconcn,fconca,fconcn_org,fconca_org,dmdt,dndt, &
        deltaso4a,cgrn3,cgra3,cb6soa_cat)
      !***********************************************************************
      !**    DESCRIPTION:  calculates aerosol nucleation and condensational
      !**    growth rates using Binkowski and Shankar (1995) method.
      !
      ! *** In this version, the method od RPM is followed where
      !     the diffusivity, the average molecular ve3locity, and
      !     the accomodation coefficient for sulfuric acid are used for
      !     the organics. This is for consistency.
      !       Future versions will use the correct values.  FSB 12/12/96
      !**
      !**    Revision history:
      !       prototype 1/95 by Uma and Carlie
      !       Corrected 7/95 by Uma for condensation of mass not nucleated
      !       and mass conservation check
      !       Revised   8/95 by US to calculate air density in stmt function
      !                 and collect met variable stmt funcs in one include fil
      !       Revised 7/25/96 by FSB to use block structure.
      !       Revised 9/17/96 by FSB to use Y&K or K&W Nucleation mechanism
      !       Revised 11/15/96 by FSB to use MKS,  and mom m^-3 units.
      !       Revised 1/13/97 by FSB to pass arrays and simplify code.
      !       Added   23/03/99 by BS growth factors for organics
      !**********************************************************************

      ! *** arguments
      ! *** input;
      integer(ik4) :: nspcsda  ! number of species in CBLK
      integer(ik4), dimension(nspcsda) :: cb6soa_cat
      real(rk8) :: cblk(nspcsda) ! main array of variables
      real(rk8) :: cblk_n(nspcsda) , cblk_a(nspcsda)
      real(rk8) :: cblk_vnu0 , cblk_vac0 , cblk_vnu3 , cblk_vac3
      real(rk8) :: dt                    ! model time step in  SECONDS
      real(rk8) :: blkta        ! Air temperature [K]
      real(rk8) :: blkprs       ! Air pressure in [Pa]
      real(rk8) :: blkrh        ! Fractional relative humidity
      real(rk8) :: so4rat       ! sulfate gas-phase production rate [ug/m^3/s]
      !bs * anthropogenic organic condensable vapor production rate
      real(rk8) :: drog(nspcsda) ! Delta ROG conc. [ug/m3]
      real(rk8), dimension(nspcsda) :: orgrat
      real(rk8) :: orgrattot
      ! biogenic organic aerosol production
      real(rk8) :: dgnuc ! nucleation mode   [m]
      real(rk8) :: dgacc ! accumulation mode [m]
      ! *** output:
      ! reciprocal condensation rate [m^3/s --> 1]
      real(rk8) :: fconcn     ! Aitken mode
      real(rk8) :: fconca     ! acclum mode 
      real(rk8) :: fconcn_org ! Aitken mode 
      real(rk8) :: fconca_org ! acclum mode 
      ! rate of production of new mass concent
      real(rk8) :: dmdt       ! by particle formation [ug/m^3/s]
      ! rate of producton of new particle numb
      real(rk8) :: dndt       ! concentration by particle formation [#/m^3/s]
      ! increment of concentration added to
      real(rk8) :: deltaso4a  ! sulfate aerosol by condensation [ug/m^3]
      ! growth rate for 3rd moment for
      real(rk8) :: cgrn3      ! Aitken mode [1/s]
      real(rk8) :: cgra3      ! Accumu mode [1/s]
      !...........    SCRATCH local variables and their descriptions:
      real(rk8) :: chemrat      ! conv rate so2 --> so4  [1/s]
      real(rk8) :: chemrat_org  ! conv rate for organics [1/s]
      real(rk8) :: am1n, am1a   ! 1st mom density (nuc, acc modes) [mom_
      real(rk8) :: am2n, am2a   ! 2nd mom density (nuc, acc modes) [mom_
      real(rk8) :: gnc3n, gnc3a ! near-cont fns (nuc, acc) for mom-3 den
      real(rk8) :: gfm3n, gfm3a ! free-mol  fns (nuc, acc) for mom-3 den
      real(rk8) :: fconc        ! total reciprocal condensation rate 

      real(rk8) :: td      ! d * tinf (cgs)
      real(rk8), parameter :: one88 = 1.0D0 ! Cnstnt to force 64bit eval of

      !  *** variables to set up sulfate and organic condensation rates
      real(rk8) :: vapor1   ! sulfuric acid vapor at current time step
      real(rk8) :: vapor2   ! Sulfuric acid vapor prior to addition from
      real(rk8) :: deltavap ! change to vapor at previous time step
      real(rk8) :: diffcorr
      real(rk8) :: csqt_org
      real(rk8) :: csqt

      integer(ik4) :: l
      !.......................................................................
      !   begin body of subroutine  NUCLCOND

      !...........   Main computational grid-traversal loop nest
      !...........   for computing condensation and nucleation:

      ! *** First moment:
      am1n = cblk_vnu0*dgnuc*esn04
      am1a = cblk_vac0*dgacc*esa04
      !..............   near-continuum factors [ 1 / sec ]
      !bs * adopted from code of FSB
      !bs * correction to DIFFSULF and DIFFORG for temperature and pressure
      diffcorr = (pss0/blkprs)*(blkta/273.16)**1.

      gnc3n = cconc*am1n*diffcorr ! [m^3/s]
      gnc3a = cconc*am1a*diffcorr
      ! *** Second moment:
      am2n = cblk_vnu0*dgnuc*dgnuc*esn16
      am2a = cblk_vac0*dgacc*dgacc*esa16

      csqt = ccofm*sqrt(blkta)
      !...............   free molecular factors [ 1 / sec ]
      ! put in temperature fac
      gfm3n = csqt*am2n ! [m^3/s]
      gfm3a = csqt*am2a
      ! *** Condensation factors in [ s**-1] for h2so4
      ! *** In the future, separate factors for condensing organics will
      !      be included. In this version, the h2so4 values are used.

      !...............   Twice the harmonic mean of fm, nc functions:
      ! *** Force 64 bit evaluation:
      fconcn = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
      fconca = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
      fconc = fconcn + fconca

      ! *** NOTE: FCONCN and FCONCA will be redefined below <<<<<<
      !bs * start modifications for organcis
      gnc3n = cconc_org*am1n*diffcorr
      gnc3a = cconc_org*am1a*diffcorr

      csqt_org = ccofm_org*sqrt(blkta)
      gfm3n = csqt_org*am2n
      gfm3a = csqt_org*am2a

      fconcn_org = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
      fconca_org = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
      !bs * end modifications for organics
      ! *** calculate the total change to sulfuric acid vapor from production
      !                      and condensation
      vapor1 = cblk(ind_SULF) ! curent sulfuric acid vapor
      vapor2 = cblk(ind_SULF) - so4rat*dt 
      vapor2 = max(0.0,vapor2) ! vapor at prev
      deltavap = max(0.0,(so4rat/fconc-vapor2)*(1.0-exp(-fconc*dt)))
      ! *** Calculate increment in total sufate aerosol mass concentration
      ! *** This follows the method of Youngblood & Kreidenweis.
      !bs * allow DELTASO4A to be negative, but the change must not be larger
      !bs * than the amount of vapor available.
      deltaso4a = max(-1.*cblk(ind_SULF),so4rat*dt-deltavap)

      ! *** zero out growth coefficients
      cgrn3 = 0.0
      cgra3 = 0.0

      if (blkta.ge.233.15 .and. blkrh.ge.0.1) then
        call klpnuc(blkta,blkrh,vapor1,dndt,dmdt,so4rat)
!       dndt = dndt/dt
!       dmdt = dmdt/dt
      else
        dndt=0.
        dmdt=0.
      end if
!     if (dndt==0.) dmdt = 0.
!     if (dmdt==0.) dndt = 0.

      !bs * Secondary organic aerosol module (SORGAM)
      ! end of selection of nucleation method
      call sorgam(blkta,blkprs,orgrat,drog,cblk,cblk_n,cblk_a, &
              nspcsda,dt,cb6soa_cat)
      !bs *  Secondary organic aerosol module (SORGAM)

      ! *** redefine FCONCN & FCONCA to be the nondimensional fractionaL
      !     condensation factors
      td = 1.0/(fconcn+fconca)
      fconcn = td*fconcn
      fconca = td*fconca

      td = 1.0/(fconcn_org+fconca_org)
      fconcn_org = td*fconcn_org
      fconca_org = td*fconca_org

      ! *** note CHEMRAT includes  species other than sulfate.
      chemrat = so4fac*so4rat ! [1/s]
      orgrattot = 0.
      do l = 1, nspcsda
        if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
          orgrattot = orgrattot + orgrat(l)
        end if
      end do
      chemrat_org = orgfac*orgrattot
      ! *** Calculate the production rates for new particle
      cgrn3 = so4fac*dmdt
      ! Rate of increase of 3rd
      chemrat = chemrat - cgrn3 !bs 3rd moment production fro
      chemrat = max(chemrat,0.0)
      ! *** Now calculate the rate of condensation on existing particles.

      ! Prevent CHEMRAT from being negativ
      cgrn3 = cgrn3 + chemrat*fconcn + chemrat_org*fconcn_org
      cgra3 = chemrat*fconca + chemrat_org*fconca_org
      return

    end subroutine nuclcond

! *** Time stepping code advances the aerosol moments one timestep;

    subroutine aerostep(nspcsda,cblk,cblk_n,cblk_a,cblk_vnu0,cblk_vac0,   &
       cblk_vnu3,cblk_vac3,dt,so4rat,orgrat,dgnuc,dgacc,fconcn,fconca,    &
       fconcn_org,fconca_org,pmassn,pmassa,dmdt,dndt,deltaso4a,           &
       urn00,ura00,brna01,c30,cgrn3,cgra3,igrid,jgrid,kgrid,cb6soa_cat)

      !***********************************************************************
      ! 
      ! ***  DESCRIPTION: Integrate the Number and Mass equations
      !                   for each mode over the time interval DT.
      !
      !      PRECONDITIONS:
      !       AEROSTEP() must follow calls to all other dynamics routines.
      !
      ! ***   Revision history:
      !       Adapted 3/95 by UAS and CJC from EAM2's code.
      !       Revised 7/29/96 by FSB to use block structure
      !       Revised 11/15/96 by FSB dropped flow-through and cast
      !                        number solver into Riccati equation form.
      !       Revised 8/8/97 by FSB to have mass in Aitken and accumulation mo
      !                        each predicted rather than total mass and
      !                        Aitken mode mass. Also used a local approximati
      !                        the error function. Also added coarse mode.
      !       Revised 9/18/97 by FSB to fix mass transfer from Aitken to
      !                        accumulation mode by coagulation
      !       Revised 10/27/97 by FSB to modify code to use primay emissions
      !                        and to correct 3rd moment updates.
      !                        Also added coarse mode.
      !       Revised 11/4/97 by FSB to fix error in other anthropogenic PM2.5
      !       Revised  11/5/97 by FSB to fix error in MSTRNSFR
      !       Revised  11/6/97 FSB to correct the expression for FACTRANS to
      !                        remove the 6/pi coefficient. UAS found this.
      !       Revised 12/15/97 by FSB to change equations for mass concentrati
      !                        to a chemical production form with analytic
      !                        solutions for the Aitken mode and to remove
      !                        time stepping of the 3rd moments. The mass conc
      !                        in the accumulation mode is updated with a forw
      !                        Euler step.
      !       Revised 1/6/98   by FSB Lowered minimum concentration for
      !                        sulfate aerosol to 0.1 [ ng / m**3 ].
      !       Revised 1/12/98  C30 replaces BRNA31 as a variable. C30 represen
      !                        intermodal transfer rate of 3rd moment in place
      !                        of 3rd moment coagulation rate.
      !       Revised 5/5/98   added new renaming criterion based on diameters
      !       Added   3/23/98  by BS condensational groth factors for organics
      !**********************************************************************

      ! *** ARGUMENTS:
      integer(ik4) :: nspcsda  ! nmber of species in CBLK        
      integer(ik4), dimension(nspcsda) :: cb6soa_cat
      real(rk8) :: cblk(nspcsda) ! main array of variables         
      real(rk8) :: cblk_n(nspcsda) , cblk_a(nspcsda)
      real(rk8) :: cblk_vnu0 , cblk_vac0 , cblk_vnu3 , cblk_vac3
      integer(ik4) :: igrid,jgrid,kgrid
      real(rk8) :: dt ! time step [s]
      ! *** Chemical production rates: [ug/m^3/s]
      real(rk8) :: so4rat          ! sulfate gas-phase production rate
      real(rk8) :: orgrat(nspcsda) ! organic gas-phase production rate
      real(rk8) :: rate            ! gas-phase production rate
      ! *** Primary emissions rates: [ ug / m**3 s ]
      ! *** emissions rates for primary organic aerosol
      ! emission rates set to zero, but kept in case of updates
      real(rk8), parameter :: eorgi = 0.0D0 ! Aitken mode
      real(rk8), parameter :: eorgj = 0.0D0 ! Accumululaton mode
      real(rk8) :: dgnuc     ! nuclei mode mean diameter [m]
      real(rk8) :: dgacc     ! accumulation
      ! reciprocal condensation rate
      real(rk8) :: fconcn    ! Aitken mode [1/s]
      real(rk8) :: fconca    ! acclum mode [1/s]
      real(rk8) :: fca , fcn ! [1/s]
      ! reciprocal condensation rate for organ
      real(rk8) :: fconcn_org! Aitken mode [1/s]
      real(rk8) :: fconca_org! acclum mode [1/s]
      real(rk8) :: dmdt      ! rate of production of new mass concent 
                             ! by particle formation [ug/m^3/s]
      real(rk8) :: dndt      ! rate of producton of new particle numb
                             ! by particle formation [#/m^3/s]
      real(rk8) :: deltaso4a ! sulfate aerosol by condensation [ug/m^3] 
      ! increment of concentration added to   
      real(rk8) :: urn00     ! Aitken intramodal coagulation rate    
      real(rk8) :: ura00     ! Accumulation mode intramodal coagulati
      real(rk8) :: brna01    ! bimodal coagulation rate for number   
      real(rk8) :: c30       ! by intermodal coagulation
      ! intermodal 3rd moment transfer rate by
      real(rk8) :: cgrn3     ! growth rate for 3rd moment for Aitken [1/s] 
      real(rk8) :: cgra3     ! growth rate for 3rd moment for Accumu [1/s]
      ! *** Modal mass concentrations 
      real(rk8) :: pmassn    ! mass concentration in Aitken mode  [ug/m^3] 
      real(rk8) :: pmassa    ! mass concentration in accumulation [ug/m^3]
      real(rk8) :: mini

      ! *** Local Variables
      integer(ik4) :: l 
      ! ** following scratch variables are used for solvers
      ! *** variables needed for modal dynamics solvers:
      ! Loop indices                   
      real(rk8) :: a, b, c
      real(rk8) :: m1, m2, y0, y
      real(rk8) :: dhat, p, pexpdt, expdt
      real(rk8) :: loss, prod, pol, lossinv, aprd
      real(rk8) :: mstrnsfr  ! mass intermodal transfer by coagulation
      real(rk8) :: factrans

      ! *** CODE additions for renaming
      real(rk8) :: getaf2
      real(rk8) :: aaa, xnum, xm3, fnum, fm3, phnum, phm3 ! Defined below
      real(rk8) :: erf, erfc ! Error and complementary error function   
      real(rk8) :: xx
      ! dummy argument for ERF and ERFC
      ! a numerical value for a minimum concentration

      ! *** This value is smaller than any reported tropospheric concentration
      ! *** Statement function given for error function. Source is
      !     Meng, Z., and J.H.Seinfeld (1994) On the source of the submicromet
      !      droplet mode of urban and regional aerosols. Aerosol Sci. and Tec
      !      20:253-265. They cite Reasearch & Education Asociation (REA), (19
      !      Handbook of Mathematical, Scientific, and Engineering Formulas,
      !      Tables, Functions, Graphs, Transforms: REA, Piscataway, NJ. p. 49
      erf(xx)  = sqrt(1.0-exp(-4.0*xx*xx/pi))
      erfc(xx) = 1.0 - erf(xx)

      ! *** set up time-step integration
      ! *** code to move number forward by one time step.
      ! *** solves the Ricatti equation:
      !     Coded 11/21/96 by Dr. Francis S. Binkowski
      ! *** Aitken mode:
      ! *** coefficients
      a = urn00
      b = brna01*cblk_vac0
      c = dndt + factnumn*(orgfac*eorgi) 

      ! includes primary emissions 
      y0 = cblk_vnu0 
      ! ***  trap on C = 0

      ! initial condition                           
      if (c>0.0D0) then
        dhat = sqrt(b*b+4.0D0*a*c)
        m1 = 2.0D0*a*c/(b+dhat)
        m2 = -0.5D0*(b+dhat)
        p = -(m1-a*y0)/(m2-a*y0)
        pexpdt = p*exp(-dhat*dt)
        y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 
        ! solution                       
      else
        ! *** rearrange solution for NUMERICAL stability
        !     note If B << A * Y0, the following form, although
        !     seemingly awkward gives the correct answer.
        expdt = exp(-b*dt)
        if (expdt<1.0D0) then
          y = b*y0*expdt/(b+a*y0*(1.0D0-expdt))
        else
          y = y0
        end if
      end if

      cblk_vnu0 = max(nummin_i,y) 

      ! *** now do accumulation mode number
      ! *** coefficients
      ! update                     
      a = ura00
      b = 0.0D0 ! NOTE B = 0.0                                         
      c = factnuma*(orgfac*eorgj) 
      ! includes primary emissi
      y0 = cblk_vac0 
      ! *** this equation requires special handling, because C can be zero.
      !     if this happens, the form of the equation is different:

      ! initial condition                           
      if (c>0.0D0) then
        dhat = sqrt(4.0D0*a*c)
        m1 = 2.0D0*a*c/dhat
        m2 = -0.5D0*dhat
        p = -(m1-a*y0)/(m2-a*y0)

        pexpdt = p*exp(-dhat*dt)
        y = (m1+m2*pexpdt)/(a*(1.0D0+pexpdt)) 
        ! solution                       
      else
        y = y0/(1.0D0+dt*a*y0) 
        ! correct solution to equatio
      end if

      cblk_vac0 = max(nummin_j,y) 

      ! *** Prepare to advance modal mass concentration one time step.
      ! *** Set up production and and intermodal transfer terms:
      cgrn3 = cgrn3 + orgfac*eorgi 

      ! includes growth from pri
      cgra3 = cgra3 + c30 + orgfac*eorgj ! & transfer of 3rd momen
                                                     ! intermodal coagulation

      ! *** set up transfer coefficients for coagulation between Aitken and ac
      ! *** set up special factors for mass transfer from the Aitken to accumu
      !     intermodal coagulation. The mass transfer rate is proportional to
      !     transfer rate, C30. The proportionality factor is p/6 times the th
      !     density. The average particle density for a species is the species
      !     divided by the particle volume concentration, pi/6 times the 3rd m
      !     The p/6 coefficients cancel.

      ! includes growth from prim
      loss = c30/cblk_vnu3 

      ! Normalized coagulation transfer r
      factrans = loss*dt  ! yields an estimate of the amount of mass t
                          ! the Aitken to the accumulation mode in the

      ! Multiplying this factor by the species con
      expdt = exp(-factrans)  ! decay term is common to all Aitken mode
      ! variable name is re-used here. This expo
      lossinv = 1.0/ loss

      ! *** now advance mass concentrations one time step.
      ! ***  update sulfuric acid vapor concentration by removing mass concent
      !      condensed sulfate and newly produced particles.
      ! *** The method follows Youngblood and Kreidenweis, Further Development
      !     of a Bimodal Aerosol Dynamics Model, Colorado State University Dep
      !     Atmospheric Science Paper Number 550, April,1994, pp 85-89.
      ! set up for multiplication rather than divi
      cblk(ind_SULF) = max(conmin,cblk(ind_SULF)-(deltaso4a+dmdt*dt))

      ! *** Solve Aitken-mode equations of form: dc/dt = P - L*c
      ! *** Solution is:     c(t0 + dt) = p/L + ( c(0) - P/L ) * exp(-L*dt)

      do l = 1, nspcsda
        if (cb6soa_cat(l) /= 0 .and. cb6soa_cat(l) /= 99) then
          if (cb6soa_cat(l) == 22) then
            rate = deltaso4a
            mini = aeroconcmin
            fca  = fconca
            fcn  = fconcn
            prod = rate*fcn + dmdt ! Condensed mass +
            aprd = rate*fca
          else
            rate = orgrat(l)
            mini = conmin
            fca  = fconca_org
            fcn  = fconcn_org
            prod = rate*fcn        ! Condensed mass +
            aprd = rate*fca*dt
          end if
          mstrnsfr = cblk_n(l)*factrans
          pol      = prod*lossinv

          cblk_n(l) = pol + (cblk_n(l)-pol)*expdt
          cblk_n(l) = max(mini,cblk_n(l))
          cblk_a(l) = cblk_a(l) + aprd + mstrnsfr
        end if
      end do
!     ! *** primary anthropogenic organic
!     mstrnsfr = cblk(l,vorgpai)*factrans
!     prod = eorgi(l)
!     pol = prod*lossinv
!
!     cblk(l,vorgpai) = pol + (cblk(l,vorgpai)-pol)*expdt
!     cblk(l,vorgpai) = max(conmin,cblk(l,vorgpai))
!     cblk(l,vorgpaj) = cblk(l,vorgpaj) + eorgj(l)*dt + mstrnsfr

      ! *** Check for mode merging,if Aitken mode is growing faster than j-mod
      !     then merge modes by renaming.
      ! *** use Binkowski-Kreidenweis paradigm, now including emissions
      ! end of time-step loop for total mass                 

!     if (cgrn3>cgra3 .or. dgnuc>.03E-6 .and. cblk_vnu0>cblk_vac0) then
      ! check if mer
      aaa  = getaf(cblk_vnu0,cblk_vac0,dgnuc,dgacc,xxlsgn,xxlsga,sqrt2)
      ! *** AAA is the value of ln( dd / DGNUC ) / ( SQRT2 * XXLSGN ), where
      !     dd is the diameter at which the Aitken-mode and accumulation-mo
      !     distributions intersect (overap).
      xnum = max(aaa,xxm3)     ! this means that no more than one ha
                               ! total Aitken mode number may be tra
                               ! per call.
      ! do not let XNUM become negative bec
      xm3  = xnum - xxm3
      ! set up for 3rd moment and mass tran
      if (xm3>0.0) then
        ! do mode merging if  overlap is corr
        phnum = 0.5*(1.0+erf(xnum))
        phm3  = 0.5*(1.0+erf(xm3))
        fnum  = 0.5*erfc(xnum)
        fm3   = 0.5*erfc(xm3)
        !     In the Aitken mode:
        ! *** FNUM and FM3 are the fractions of the number and 3rd moment
        !     distributions with  diameters greater than dd respectively.
        ! *** PHNUM and PHM3 are the fractions of the number and 3rd moment
        !     distributions with diameters less than dd.
        ! *** rename the  Aitken mode particle number as accumulation mode
        !     particle number
        cblk_vac0 = cblk_vac0 + fnum*cblk_vnu0
        ! *** adjust the Aitken mode number
        cblk_vnu0 = phnum*cblk_vnu0

        do l = 1, nspcsda
          if (cb6soa_cat(l) /= 0 .and. cb6soa_cat(l) /= 99) then
            ! Rename mass from Aitken mode to acumula mode. The mass transfe
            ! to the accumula mode is proportional to the amount of 3rd mome
            ! transferred, therefore FM3 is used for mass transfer.
            cblk_a(l) = cblk_a(l) + cblk_n(l)*fm3
            ! update Aitken mode for mass loss to accumulation mode
            cblk_n(l) = cblk_n(l)*phm3
          end if
        end do
      end if
      ! end check on whether modal overlap is OK             
!     end if
      ! end check on necessity for merging                   
      !     set min value for all concentrations

      ! loop for merging                                       
      do l = 1, nspcsda
        if (cb6soa_cat(l) /= 0 .and. cb6soa_cat(l) /= 99) then
          cblk(l)   = max(cblk(l),conmin)
          cblk_n(l) = max(cblk_n(l),conmin)
          cblk_a(l) = max(cblk_a(l),conmin)
        end if
      end do

      return

    end subroutine aerostep

    ! main box model
    subroutine sorgam(blkta,blkprs,orgrat,drog,cblk,cblk_n,cblk_a, &
        nspcsda,dt,cb6soa_cat)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  SOA_PART calculates the formation and partitioning of secondary   !
      !bs  organic aerosol based on (pseudo-)ideal solution thermodynamics.  !
      !bs                                                                    !
      !bs  This code considers two cases:                                    !
      !bs   i) initil absorbing mass is existend in the aerosol phase        !
      !bs  ii) a threshold has to be exeeded before partitioning (even below !
      !bs      saturation) will take place.                                  !
      !bs                                                                    !
      !bs  The temperature dependence of the saturation concentrations are   !
      !bs  calculated using the Clausius-Clapeyron equation.                 !
      !bs                                                                    !
      !bs  It is assumed that the condensable vapors also evaporate if the   !
      !bs  saturation concentraion lowers e.g. due to temperature effects.   !
      !bs  Therefor negative production rates (= evaporation rates) are      !
      !bs  possible.                                                         !
      !bs                                                                    !
      !bs  If there is no absorbing mass at all the Pandis method is applied !
      !bs  for the first steps.                                              !
      !bs                                                                    !
      !bs  References:                                                       !
      !bs    Pankow (1994):                                                  !
      !bs     An absorption model of the gas/aerosol                         !
      !bs     partitioning involved in the formation of                      !
      !bs     secondary organic aerosol, Atmos. Environ. 28(2),              !
      !bs     189-193.                                                       !
      !bs    Odum et al. (1996):                                             !
      !bs     Gas/particle partitioning and secondary organic                !
      !bs     aerosol yields,  Environ. Sci. Technol. 30(8),                    !
      !bs     2580-2585.                                                     !
      !bs    Bowman et al. (1997):                                           !
      !bs     Mathematical model for gas-particle partitioning               !
      !bs     of secondary organic aerosols, Atmos. Environ.                 !
      !bs     31(23), 3921-3931.                                             !
      !bs    Seinfeld and Pandis (1998):                                     !
      !bs     Atmospheric Chemistry and Physics (0-471-17816-0)              !
      !bs     chapter 13.5.2 Formation of binary ideal solution              !
      !bs     with -- preexisting aerosol                                    !
      !bs          -- other organic vapor                                    !
      !bs * Pandis et al. (1992): Secondary organic aerosol formation and
      !bs *     transport. Atmos Environ. 26A, 2453-2466.
      !bs * STI Report (Sonoma Technology, Inc.) (1998):
      !bs *     Development of gas-phase chemistry, secondary organic aerosol,
      !bs *     and aqueous-phase chemistry modules for PM modeling.
      !bs *     By: R. Strader, C. Gurciullo, S. Pandis, N. Kumar, F. Lurmann
      !bs *     Prepared for: Coordinating Research Council, Atlanta, Aug 24 1
      !bs * Tao and McMurray (1989): Vapor pressures and surface free energies
      !bs *     C14-C18 monocarboxylic acids and C5 and C6 dicarboxylic acids.
      !bs *     Eniron. Sci. Technol. 23, 1519-1523.
      !bs * Pankow (1994): An absorption model of gas/particle partitioning of
      !bs *     organic compounds in the atmosphere. Atmos. Environ. 28, 185-1
      !bs                                                                    !
      !bs--------------------------------------------------------------------!
      !bs                                                                    !
      !bs  History:                                                          !
      !bs   No    Date    Author           Change                            !
      !bs  ____  ______  ________________  _________________________________ !
      !bs   01   170399   B.Schell         Set up                            !
      !bs   02   050499   B.Schell         introduced SR NEWT                !
      !bs   03   040599   B.Schell         include-file sorgam.inc           !
      !jc   04   060416   J.Ciarlo         Adapted for RegCM4.5              !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: layer        ! model layer
      integer(ik4) :: nspcsda      ! number of species in CBLK
      integer(ik4) :: ncv          ! total # of cond. vapors & SOA sp
      integer(ik4) :: nacv         ! # of anthrop. cond. vapors & SOA
      real(rk8) :: cblk(nspcsda)   ! main array of variables
      real(rk8) :: cblk_n(nspcsda)
      real(rk8) :: cblk_a(nspcsda)
      real(rk8) :: dt              ! model time step [s]
      real(rk8) :: blkta           ! Air temperature [K]
      real(rk8) :: blkprs          ! Air pressure in [Pa]
      real(rk8) :: orgrat(nspcsda) ! rates of specific organic vap production
      real(rk8) :: drog(nspcsda)   ! Delta ROG conc. [ug/m3]
      integer(ik4), dimension(nspcsda) :: cb6soa_cat

      !bs * local variable declaration
      real(rk8), parameter :: thrsmin=1.E-19 ! numerical value for a min thresh
      real(rk8), parameter :: tnull=298.     ! reference temperature T0 = 298 K
!     real(rk8), parameter :: rtol=1.E-04    ! relative tolerance for mass check

      integer(ik4) :: l, n    !bs loop index
      real(rk8) :: ttinv      !bs difference of inverse temperatures
      real(rk8) :: minitw     !bs weighted initial organic mass [10^-6
      real(rk8) :: mtotw      !bs weighted total organic mass [10^-6 m
      real(rk8) :: mnonow     !bs weighted inorganic mass [10^-6 mol/m
      real(rk8) :: imtotw     !bs 1. / MTOTW
      real(rk8) :: minit      !bs initial organic mass [ug/m^3]
      real(rk8) :: mnono      !bs inorganic mass [ug/m^3]
      real(rk8) :: mtot       !bs total organic mass [ug/m^3]
      real(rk8) :: thres      !bs threshold for SOA formatio for low M
      real(rk8) :: mcheck     !bs mass check ratio of input/output mas
!     real(rk8) :: msum(nspcsda)  !bs input total mass [ug/m^3]
      real(rk8) :: mwcv(nspcsda)  !bs molecular weight of cond. vapors [g/
      real(rk8) :: imwcv(nspcsda) !bs 1. / MWCV(NCV)
      real(rk8) :: pnull(nspcsda) !bs vapor pres. of pure cond. vapor [Pa]
      real(rk8) :: dhvap(nspcsda) !bs heat of vaporisation of compound i [
      real(rk8) :: pvap(nspcsda)  !bs vapor pressure cond. vapor [Pa]
      real(rk8) :: ctot(nspcsda)  !bs total conc. of cond. vapor aerosol +
      real(rk8) :: cgas(nspcsda)  !bs gasphase concentration of cond. vapo
      real(rk8) :: caer(nspcsda)  !bs aerosolphase concentration of cond.
      real(rk8) :: asav(nspcsda)  !bs saved CAER for iteration
      real(rk8) :: aold(nspcsda)  !bs saved CAER for rate determination
      real(rk8) :: csat(nspcsda)  !bs saturation conc. of cond. vapor [ug/m^3]
      real(rk8) :: alpha(nspcsda) !bs molar yield for condensable vapors
      real(rk8) :: prod(nspcsda)  !bs production of condensable vapor [ug/
      real(rk8) :: prod_all       ! concentration sum all dROG species 
!     real(rk8) :: f(nspcsda)   !bs scaling factor for ind. oxidant
      logical :: check   !bs check convergence of SR NEWT
      integer(ik4) :: its
      !bs * initialisation
      !bs
      !bs * DVAP data: average value calculated from C14-C18 monocarboxylic
      !bs *      acids and C5-C6 dicarboxylic acids. Tao and McMurray (1989):
      !bs *      Eniron. Sci. Technol. 1989, 23, 1519-1523.
      !bs *      average value is 156 kJ/mol
      !bs
      !bs number of iterations in NEWT
!     dhvap(psoaaro1) = 156.0E03
      dhvap(ind_ALKA) = 26.2E03 ! Pentane (paraffin) C5H12
      dhvap(ind_ALDX) = 30.3E03 ! Proprionaldehyde & higher CH3CH2CHO
      dhvap(ind_ISPD) = 29.0E03 ! Methacrolein C4H6O
      dhvap(ind_AALD) = 26.0E03 ! Acetaldehyde CH3CHO
      dhvap(ind_OLE ) = 18.7E03 ! Propene (terminal) CH3CHCH2
      dhvap(ind_FORM) = 23.3E03 ! Formaldehyde HCHO
      dhvap(ind_IOLE) = 23.6E03 ! cis/trans-2-butene (internal) C4H8
      dhvap(ind_ROPN) = 48.4E03 ! Aromatic ring opening product C4H4O2 
                                !! 2-Butynoic acid
      dhvap(ind_MEGY) = 31.3E03 ! Methyl glyoxal CH3COCHO
      dhvap(ind_XOPN) = 44.9E03 ! Unsaturated dicarbonyl C6H8O2 
                                !! (3Z)-hex-3-ene-2,5-dione
      dhvap(ind_GLY ) = 29.4E03 ! Glyoxal CHOCHO
      dhvap(ind_GLYD) = 43.0E03 ! Glycolaldehyde C2H4O3
      dhvap(ind_CRSL) = 62.0E03 ! Cresol C6H4(OH)CH3
      dhvap(ind_AACD) = 42.0E03 ! Acetic acid CH3COOH
      dhvap(ind_CRON) = 52.3E03 ! Nitrocresol C6H3(OH)(NO2)CH3
      dhvap(ind_INTR) = 43.6E03 ! Simple organic nitrate C4H9O3N !! butyl nitrate
      dhvap(ind_ACET) = 32.1E03 ! Acetone CH3COCH3
      dhvap(ind_NTR1) = 83.8E03 ! Organic nitrate  C5H9O4N !! Glutamic acid
      dhvap(ind_FACD) = 29.6E03 ! Formic acid HCOOH
      dhvap(ind_CAT1) = 22.8E03 ! Pyrocatechol C6H3(OH)2CH3
      dhvap(ind_MEOH) = 37.3E03 ! Methanol CH3OH
      dhvap(ind_NTR2) = 83.8E03 ! multifunctional organic nitrates C4H9O4N 
                                !! Glutamic acid
      dhvap(ind_KET ) = 31.3E03 ! Methy ethyl ketone CH3COC2H5
      dhvap(ind_EPOX) = 70.2E03 ! Isoprene epoxide derivative C5H10O4 
                                !! Conduritol B epoxide C6H10O5
!      dhvap(ind_NTR ) = 43.6E03 ! 1/2-butyl nitrate C4H8NO3
     ! Sources of enthalpies given the E03 to match units
     ! webbook.nist.gov/  pubchem.ncbi.nlm.nih.gov/  www.chemicalbook.com/
     ! www.lookchem.com/  www.chemspider.com/ www.chemeo.com/  www.phs.d211.org/

      !bs * aromatic yields from:
      !bs * Odum J.R., T.P.W. Jungkamp, R.J. Griffin, R.C. Flagan, and
      !bs *   J.H. Seinfeld: The atmospheric aerosol-forming potential of whol
      !bs *   gasoline vapor, Science 276, 96-99, 1997.
      !bs * Odum J.R., T.P.W. Jungkamp, R.J. Griffin, H.J.L. Forstner, R.C. Fl
      !bs *   and J.H. Seinfeld: Aromatics, reformulated gasoline, and atmosph
      !bs *   organic aerosol formation, Environ. Sci. Technol. 31, 1890-1897,
      !bs *
      !bs * !! yields provided by Odum are mass-based stoichiometric coefficen
      !bs *    average for high and low yield aromatics
      !bs *    alpha1 = 0.0545  K1 = 0.0475 m^3/ug
      !bs *    alpha2 = 0.1525  K2 = 0.00165 m^3/ug
      !bs *    change to molar yields using the model MW
      !bs *    alpha1 * MW(XYL) / MW(PSOAARO1) = alpha1 * 106 / 150 = 0.0385
      !bs *    alpha2 * MW(XYL) / MW(PSOAARO2) = alpha2 * 106 / 150 = 0.1077
      !bs *   ALPHA(PSOAARO1) = 0.0385; ALPHA(PSOAARO2) = 0.1077
      !bs *
      !bs
      !bs * alkane and alkene yields from:
      !bs * Moucheron M.C. and J. Milford: Development and testing of a proces
      !bs *    model for secondary organic aerosols. Air & Waste Manag. Assoc.
      !bs *    for presentation at the 89th Annual Meeting & Exhibition, Nashv
      !bs *    Tennessee, June 23-28, 96-FA130B.03, 1996.
      !bs *  molar yields used instead of [ ug m^-3 / ppm ], calculation
      !bs *    at T=298K, P=1.0133*10^5 Pa
      !bs *    ALPHA(PSOAALK1) = 0.048; ALPHA(PSOAOLE1) = 0.008
      !bs
      !bs * biogenic yields from:
      !bs * Griffin R.J., D.R. Cocker III, R.C. Flagan, and J.H. Seinfeld:
      !bs *   Organic aerosol formation from the oxidation of biogenic hydro-
      !bs *   carbons, JGR, 1999 in press.
      !bs *   the yields given in Table 3 are mass yields [ ug m^-3 / ug m^-3
      !bs *   change to molar yields via:
      !bs *   molar yield = mass yield * ((R*T/M_soa*p) / (R*T/M_terp*p))
      !bs *               = mass yield * (M_terp / M_soa)
      !bs *               = mass yield * ( M(Terpenes) / M(pinonic acid) )
      !bs *               = mass yield * 136 / 184
      !bs * average for a-Pinene and Limonene, maybe splitted in future versio
      !bs *    0.138 * 0.739 = 0.102; 0.345 * 0.739 = 0.254
      !bs * values for a-Pinene (molar yield) alpha1 = 0.028, alpha2 = 0.241
      !bs * values for limonene (molar yield) alpha1 = 0.163, alpha2 = 0.247
      ! old values for alpha
      !alpha(psoaaro1) = 0.039
      !alpha(psoaaro2) = 0.108
      !alpha(psoaalk1) = 0.048
      !alpha(psoaole1) = 0.008
      !alpha(psoaapi1) = 0.092 !bs API + O3 only Griffin '99
      !alpha(psoaapi2) = 0.075 !bs API + O3 only Griffin '99
      !alpha(psoalim1) = 0.163
      !alpha(psoalim2) = 0.247

      ! alpha(product) = stoich coeff * mw_drog / mw_product
      ! large complexity of chemical system requires calculation out of
      ! another script written in RStudio 
      alpha(ind_ALKA) = 0.214340801705917
      alpha(ind_ALDX) = 0.0419399104108558
      alpha(ind_ISPD) = 0.90879509134125
      alpha(ind_AALD) = 0.110568853979084
      alpha(ind_OLE ) = 0.363475405135388
      alpha(ind_FORM) = 0.558394347807896
      alpha(ind_IOLE) = 0.0290099244974528
      alpha(ind_ROPN) = 0.226111135419519
      alpha(ind_MEGY) = 0.298082305716648
      alpha(ind_XOPN) = 0.220464628656317
      alpha(ind_GLY ) = 0.229109348447759
      alpha(ind_GLYD) = 0.107153527641608
      alpha(ind_CRSL) = 0.310979830722132
      alpha(ind_AACD) = 0.0651273399041627
      alpha(ind_CRON) = 0.124137292513098
      alpha(ind_INTR) = 0.0787471182901095
      alpha(ind_ACET) = 0.0476074658733809
      alpha(ind_NTR1) = 0.0350203213867376
      alpha(ind_FACD) = 0.10848846411915
      alpha(ind_CAT1) = 0.198264689496523
      alpha(ind_MEOH) = 0.0147721483720793
      alpha(ind_NTR2) = 0.0931620192242024
      alpha(ind_KET ) = 0.0182621231518719
      alpha(ind_EPOX) = 0.377920333957553
!      alpha(ind_NTR ) = 0.0993354076482532

      !bs
      !bs * P0 data in Pa for T = 298K:
      !bs *    aromatics: Odum et al. (1997) using R = 8.314 J/(mol*K),
      !bs *         DHvap = 156 kJ/mol, T = 313K, MW = 150 g/mol and averaged
      !bs *         Ki's of high and low aromatics.
      !bs *         T = 313   => PNULL(ARO1) = 1.7E-05, PNULL(ARO2) = 5.1E-04
      !bs *         T = 307.4 => PNULL(ARO1) = 5.7E-05, PNULL(ARO2) = 1.6E-03
      !bs *    biogenics: Hoffmann et al. (1997); Griffin et al. (1999);
      !bs *         using R = 8.314 J/(mol*K),
      !bs *         DHvap = 156 kJ/mol, T = 313, MW = 184 g/mol, and
      !bs *         averaged Ki's of a-pinene and limonene
      !bs *         p1(298K) = 6.1E-06; p2(298K) = 1.5E-04
      !bs *         Ki's for a-pinene p1(298K) = 4.0E-06; p2(298K) = 1.7E-04
      !bs *         Ki's for limonene p1(298K) = 2.5E-05; p2(298K) = 1.2E-04
      !bs *    alkanes and alkenes: no data available, use low value to get cl
      !bs *         to the Pandis yields, 1 ppt = 1*10^-7 Pa.
      !bs ------ original sorgam pnull source ------
      !pnull(psoaaro1) = 5.7E-05
      !pnull(psoaaro2) = 1.6E-03
      !pnull(psoaalk1) = 5.0E-06
      !pnull(psoaole1) = 5.0E-06
      !pnull(psoaapi1) = 2.488E-05 !bs API + O3 only Griffin '99
      !pnull(psoaapi2) = 2.778E-05 !bs API + O3 only Griffin '99
      !pnull(psoalim1) = 2.5E-05
      !pnull(psoalim2) = 1.2E-04
      ! -------- pnull values matched from original sorgam pnull values
      do l = 1, nspcsda
        if (cb6soa_cat(l) == 2 ) pnull(l) = 5.0E-06 !carbonyl         (psoaole1)
        if (cb6soa_cat(l) == 3 ) pnull(l) = 5.0E-06 !alkanes          (psoaalk1)
        if (cb6soa_cat(l) == 4 ) pnull(l) = 5.0E-06 !olefins          (psoaole1)
        if (cb6soa_cat(l) == 6 ) pnull(l) = 5.7E-05 !aromatics 	      (psoaaro1)
        if (cb6soa_cat(l) == 9 ) pnull(l) = 5.0E-06 !organic nitrates (psoaalk1)
        if (cb6soa_cat(l) == 10) pnull(l) = 2.5E-05 !multifunctional  (psoalim1)
      end do
      !bs
      !bs * scaling of contribution of individual oxidants to aerosol formatio
      !bs
!      f(pxyl)  = 1. !bs * XYL + OH
!      f(ptol)  = 1. !bs * TOL + OH
!      f(pcsl1) = 1. !bs * CSL + OH
!      f(pcsl2) = 1. !bs * CSL + NO
!      f(phc8)  = 1. !bs * HC  + OH
!      f(poli1) = 1. !bs * OLI + OH
!      f(poli2) = 1. !bs * OLI + NO
!      f(poli3) = 1. !bs * OLI + O3
!      f(polt1) = 1. !bs * OLT + OH
!      f(polt2) = 1. !bs * OLT + NO
!      f(polt3) = 1. !bs      F(PAPI1) = 0.228          !bs * API + OH
!      f(papi1) = 0. !bs * API + OH
!      f(papi2) = 0. !bs * API + NO
!      f(papi3) = 1. !bs * API + O3
!      f(plim1) = 0.228 !bs * LIM + OH
!      f(plim2) = 0. !bs * LIM + NO
!      f(plim3) = 0.771 !bs !bs * LIM + O3

      !bs * begin code -------------------------------------------------------
!     do l = 1, ldrog
!       drog(l) = f(l)*drog(l)
!     end do
      ttinv = 1./tnull - 1./blkta
      prod_all = 0.
      prod(:)  = 0.
      do l = 1, nspcsda
        if (cb6soa_cat(l) /= 0) then
          cgas(l) = cblk(l)
          caer(l) = cblk_n(l) + cblk_a(l)
          mwcv(l) = mw_cb6(l) 
          if (cb6soa_cat(l) == 99) then
            prod_all = prod_all + drog(l)
          end if
        end if
      end do
      
      prod(ind_ALKA) = drog(ind_ISPR) + drog(ind_TERP) + drog(ind_PRPA)
      prod(ind_ALDX) = prod_all
      prod(ind_ISPD) = drog(ind_ISPR)
      prod(ind_AALD) = prod_all
      prod(ind_OLE ) = drog(ind_ISPR)
      prod(ind_FORM) = prod_all
      prod(ind_IOLE) = drog(ind_ISPR)
      prod(ind_ROPN) = drog(ind_BENZ) + drog(ind_TOLN) + drog(ind_XYLN)
      prod(ind_MEGY) = drog(ind_ISPR) + drog(ind_BENZ) + drog(ind_TOLN) &
                         + drog(ind_XYLN)
      prod(ind_XOPN) = drog(ind_BENZ) + drog(ind_TOLN) + drog(ind_XYLN)
      prod(ind_GLY ) = drog(ind_ETHY) + drog(ind_ETOH) + drog(ind_ETHE) &
                         + drog(ind_ISPR) + drog(ind_BENZ)              &
                         + drog(ind_TOLN) + drog(ind_XYLN) 
      prod(ind_GLYD) = drog(ind_ETOH) + drog(ind_ETHE) + drog(ind_ISPR)
      prod(ind_CRSL) = drog(ind_BENZ) + drog(ind_TOLN) + drog(ind_XYLN)
      prod(ind_AACD) = prod_all
      prod(ind_CRON) = drog(ind_BENZ) + drog(ind_TOLN) + drog(ind_XYLN)
      prod(ind_INTR) = drog(ind_ISPR) 
      prod(ind_ACET) = drog(ind_ISPR) + drog(ind_TERP) + drog(ind_PRPA)
      prod(ind_NTR1) = drog(ind_ETHE) + drog(ind_ISPR) + drog(ind_TERP) & 
                         + drog(ind_PRPA)
      prod(ind_FACD) = prod_all
      prod(ind_CAT1) = drog(ind_BENZ) + drog(ind_TOLN) + drog(ind_XYLN)
      prod(ind_MEOH) = prod_all
      prod(ind_NTR2) = prod_all 
      prod(ind_KET ) = drog(ind_ISPR) + drog(ind_TERP) + drog(ind_PRPA)
      prod(ind_EPOX) = drog(ind_ISPR)
!      prod(ind_NTR ) = drog(ind_TERP)

      !bs * calculate actual production from gasphase reactions [ug/m^3]
      !bs * calculate vapor pressure of pure compound as a liquid
      !bs *   using the Clausius-Clapeyromn equation and the actual
      !bs *   saturation concentration.
      !bs * calculate the threshold for partitioning if no initial mass
      !bs *   is present to partition into.
      thres = 0.
      mtot = 0.
      mtotw = 0.
      do l = 1, nspcsda
        if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
          prod(l)  = alpha(l)*prod(l)
          ctot(l)  = prod(l) + cgas(l) + caer(l) !bs redefined below
          !msum(l) = cgas(l) + caer(l) + prod(l)
          aold(l)  = caer(l)
          imwcv(l) = 1./mwcv(l)
          pvap(l)  = pnull(l)*exp(dhvap(l)/rgas*ttinv)
          csat(l)  = pvap(l)*mwcv(l)*1.0E06/(rgas*blkta)
          thres    = thres + ((cgas(l)+prod(l))/csat(l))
          mtot     = mtot  + caer(l)
          mtotw    = mtotw + caer(l)*imwcv(l)
        end if
      end do

      !bs * small amount of non-volatile absorbing mass is assumed to be
      !bs * present (following Bowman et al. (1997) 0.01% of the inorganic
      !bs * mass in each size section, here mode)
      mnono  = 0.0001*cblk_a(ind_SULF) + 0.0001*cblk_n(ind_SULF)
      mnonow = 0.0001*(cblk_a(ind_SULF)/mwcv(ind_SULF))     &
                + 0.0001*(cblk_n(ind_SULF)/mwcv(ind_SULF))
      mnono  = max(mnono,conmin)
      mnonow = max(mnonow,conmin)

    ! minit  = cblk(lcell,vorgpaj) + cblk(lcell,vorgpai) + mnono
    ! minitw = (cblk(lcell,vorgpaj)+cblk(lcell,vorgpai))/mworg + mnonow

      !bs * If MINIT is set to zero partitioning will occur if the pure
      !bs * saturation concentation is exceeded (Pandis et al. 1992).
      !bs * If some amount of absorbing organic mass is formed gas/particle
      !bs * partitioning will follow the ideal solution approach.
    ! minit = 0.
    ! minitw = 0.

      mtot = mtot + mnono
      mtotw = mtotw + mnonow
      imtotw = 1./mtotw

      !bs * do the gas/particle partitioning
      if ((thres>1 .and. mnonow<thrsmin) .or. (mnonow>thrsmin) .or. &
           (mtot>thrsmin)) then
           do l = 1, nspcsda
             if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
               caer(l) = ctot(l) !bs 'initial' guess
             end if
           end do

        !bs * globally convergent method for nonlinear system of equations
        !bs * adopted from Numerical Recipes 2nd Edition
        call newt(caer,ncv,check,ctot,csat,imwcv,minitw,its)

        do l = 1, nspcsda
          if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
            if (caer(l)<=tolmin) then
              caer(l) = conmin
            end if
            if (caer(l)>ctot(l)) then
              caer(l) = ctot(l)
            end if
!           cgas(l) = ctot(l) - caer(l)
!           cblk(l) = max(cgas(l),conmin)
            orgrat(l) = (caer(l)-aold(l))/dt
          end if
        end do
      else
        !bs do Pandis method
        do l = 1, nspcsda
          if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
            caer(l) = ctot(l) - csat(l)
            caer(l) = max(caer(l),0.)
!           cgas(l) = ctot(l) - caer(l)
!           cblk(l) = cgas(l)
            orgrat(l) = (caer(l)-aold(l))/dt
          end if
        end do
      end if
      !bs * check mass conservation
!      do l = 1, nspcsda
!         if (cb6soa_cat(l) > 0 .and. cb6soa_cat(l) < 20 ) then
!            !rs check is component exits
!            if (cgas(l)==0. .and. caer(l)==0. .and. msum(l)==0) then
!              mcheck = 1.
!            else
!              mcheck = (cgas(l)+caer(l))/msum(l)
!            end if
!            if ((mcheck<1.-rtol) .or. (mcheck>1.+rtol)) then
!90020         FORMAT ('LAYER = ',I2,', L = ',I2,', MCHECK = ',E12.6,', MASS = ',&
!                E12.6)
!            end if
!         end if
!      end do
      return
    end subroutine sorgam  

    subroutine klpnuc(temp,rh,h2so4,ndot1,mdot1,so4rat)
      real(rk8) :: temp   ! ambient temperature [K]
      real(rk8) :: rh     ! fractional relative humidity
      real(rk8) :: h2so4  ! sulfuric acid concentration [ug/m^3]
      real(rk8) :: so4rat ! sulfuric acid production rate [ug/m^3/s]
      real(rk8) :: ndot1  ! particle number production rate [#/m^3/s]
      real(rk8) :: mdot1  ! particle mass production rate [ug/m^3/s]
      real(rk8) :: m2dot  ! particle scnd momnt prod rate [m^2/m^3/s]
      real(rk8) :: ra     ! fractional relative acidity
      real(rk8) :: nav    ! sulfuric acid vaper concentration [1/cm^3]
      real(rk8) :: nwv    ! water vapor concentration [1/cm^3]
      real(rk8) :: nav0   ! equilibrium sulfuric acid vapor conc. [1/cm^3]
                          ! to produce a nucleation rate of 1 [1/cm^3/s]
      real(rk8) :: nac    ! critical sulfuric acid vapor concentration [1/cm^3]
      real(rk8) :: xal    ! mole fractio of the critical nucleus
      real(rk8) :: nsulf, delta
      real(rk8) :: chi    ! factor to calculate Jnuc
      real(rk8) :: jnuc   ! nucleation rate [1/cm^3/s]
      real(rk8) :: tt, rr ! dummy variables for statement functions
      real(rk8) , parameter :: pid6 = pi/6.0
      real(rk8) , parameter :: avo = 6.0221367E23   ! avogadro's constant [1/mol]
      real(rk8) , parameter :: atm = 1013.25E+02    ! 1 atmosphere in pascals
      real(rk8) , parameter :: mwh2so4 = 98.07948   ! mol wt h2so4 [g/mol]
      real(rk8) , parameter :: d35 = 3.5E-07        ! dia of 3.5 nm part [cm]
      real(rk8) , parameter :: d35sq = d35*d35 
      real(rk8) , parameter :: v35 = pid6*d35*d35sq ! vol of 3.5 nm part [cm^3]
      real(rk8) :: mp     ! mass of sulfate in a 3.5 nm particle number/cm3
      ! ***  conversion factors:
      real(rk8) , parameter :: ugm3_ncm3 = (avo/mwh2so4)*1.0E-12 ! ug/m3->mlc/cm3
      real(rk8) , parameter :: nc_ug = (1.0E6)*mwh2so4/avo
      ! *** statement functions **************
      real(rk8) :: pdens, rho_p ! particle density [g/cm^3]
      ! coefficients for density expression
      real(rk8) , parameter :: ad0 =  1.738984
      real(rk8) , parameter :: ad1 = -1.882301
      real(rk8) , parameter :: ad2 =  2.951849
      real(rk8) , parameter :: ad3 = -1.810427

      ! *** Nair and Vohra, Growth of aqueous sulphuric acid droplets
      !     as a function of relative humidity,
      !     J. Aerosol Science, 6, pp 265-271, 1975.

      ! fit to Nair & Vohra data
      real(rk8) :: mp35 ! the mass of sulfate in a 3.5 nm particle
      ! arithmetic statement function to compute
      real(rk8) , parameter :: a0 =  1.961385E2 ! coefficients for cubic in mp35
      real(rk8) , parameter :: a1 = -5.564447E2
      real(rk8) , parameter :: a2 =  8.828801E2
      real(rk8) , parameter :: a3 = -5.231409E2
      real(rk8) :: ph2so4, ph2o     ! for h2so4 and h2o vapor pressures [Pa]

      ! arithmetic statement functions
      pdens(rr)  = ad0 + rr*(ad1+rr*(ad2+rr*ad3))
      ph2o(tt)   = exp(77.34491296-7235.4246512/tt-8.2*log(tt)+tt*5.7113E-03)
      ph2so4(tt) = exp(27.78492066-10156.0/tt)
      ! *** both ph2o and ph2so4 are  as in Kulmala et al.  paper

      ! *** function for the mass of sulfate in   a 3.5 nm sphere
      ! *** obtained from a fit to the number of sulfate monomers in
      !     a 3.5 nm particle. Uses data from Nair & Vohra
      mp35(rr) = nc_ug*(a0+rr*(a1+rr*(a2+rr*a3)))

      !     The 1.0e-6 factor in the following converts from MKS to cgs units
      ! *** get water vapor concentration [ molecles / cm **3 ]
      nwv = rh*ph2o(temp)/(rgas*temp)*avo*1.0E-6

      ! *** calculate the equilibrium h2so4 vapor concentration.
      ! *** use Kulmala corrections:
      nav0 = ph2so4(temp)/(rgas*temp)*avo*1.0E-6

      ! *** convert sulfuric acid vapor concentration from micrograms
      !     per cubic meter to molecules per cubic centimeter.
      nav = ugm3_ncm3*h2so4

      ! *** calculate critical concentration of sulfuric acid vapor
      nac = exp(-14.5125+0.1335*temp-10.5462*rh+1958.4*rh/temp)

      ! *** calculate relative acidity
      ra = nav/nav0

      ! *** calculate temperature correction
      delta = 1.0 + (temp-273.15)/273.14

      ! *** calculate molar fraction
      xal = 1.2233 - 0.0154*ra/(ra+rh) + 0.0102*log(nav) - 0.0415*log(nwv) + &
        0.0016*temp
      ! *** calculate Nsulf
      nsulf = log(nav/nac)

      ! *** calculate particle produtcion rate [ # / cm**3 ]
      chi = 25.1289*nsulf - 4890.8*nsulf/temp - 1743.3/temp - &
        2.2479*delta*nsulf*rh + 7643.4*xal/temp - 1.9712*xal*delta/rh

      jnuc = exp(chi)         ! [ # / cm**3 ]
      ndot1 = (1.0E06)*jnuc

      ! *** calculate particle density
      rho_p = pdens(rh)

      ! *** get the mass of sulfate in a 3.5 nm particle
      mp = mp35(rh)           ! in a 3.5 nm particle at ambient RH

      ! *** calculate mass production rate [ ug / m**3]
      !     assume that the particles are 3.5 nm in diameter.
!     MDOT1 =  (1.0E12) * rho_p * v35 * Jnuc

      ! number of micrograms of sulfate
      mdot1 = mp*ndot1

      if (mdot1>so4rat) then
        mdot1 = so4rat    ! limit nucleated mass by available ma
        ndot1 = mdot1/ mp ! adjust DNDT to this
      end if

      if (mdot1==0.) ndot1 = 0.

     ! *** calculate M2 production rate [ m**2 / (m**3 s)]
      m2dot = 1.0E-04*d35sq*ndot1

      return
    end subroutine klpnuc

    subroutine fdjac(n,x,fjac,ct,cs,imw)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Get the Jacobian of the function                                  !
      !bs                                                                    !
      !bs         ( a1 * X1^2 + b1 * X1 + c1 )                               !
      !bs         ( a2 * X2^2 + b2 * X1 + c2 )                               !
      !bs         ( a3 * X3^2 + b3 * X1 + c3 )                               !
      !bs  F(X) = ( a4 * X4^2 + b4 * X1 + c4 ) = 0.                          !
      !bs         ( a5 * X5^2 + b5 * X1 + c5 )                               !
      !bs         ( a6 * X6^2 + b6 * X1 + c6 )                               !
      !bs                                                                    !
      !bs   a_i = IMW_i                                                      !
      !bs   b_i = SUM(X_j * IMW_j)_j.NE.i + CSAT_i * IMX_i - CTOT_i * IMW_i  !
      !bs   c_i = - CTOT_i * [ SUM(X_j * IMW_j)_j.NE.i + M ]                 !
      !bs                                                                    !
      !bs          delta F_i    ( 2. * a_i * X_i + b_i           if i .EQ. j !
      !bs  J_ij = ----------- = (                                            !
      !bs          delta X_j    ( X_i * IMW_j - CTOT_i * IMW_j   if i .NE. j !
      !bs                                                                    !
      !bs                                                                    !
      !bs  Called by:       NEWT                                             !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      !implicit none
      integer(ik4) :: n   !dimension of problem                   
      real(rk8) :: x(n) 
      real(rk8) :: ct(np) !bs initial guess of CAER               
      real(rk8) :: cs(np)
      real(rk8) :: imw(np)
      real(rk8) :: fjac(n,n)

      integer(ik4) :: i, j !bs loop index  
      real(rk8) :: a(np)
      real(rk8) :: b(np)
      real(rk8) :: b1(np)
      real(rk8) :: b2(np)
      real(rk8) :: sum_jnei

      do i = 1, n
        a(i) = imw(i)
        sum_jnei = 0.
        do j = 1, n
          sum_jnei = sum_jnei + x(j)*imw(j)
        end do
        b1(i) = sum_jnei - (x(i)*imw(i))
        b2(i) = cs(i)*imw(i) - ct(i)*imw(i)
        b(i) = b1(i) + b2(i)
      end do
      do j = 1, n
        do i = 1, n
          if (i==j) then
            fjac(i,j) = 2.*a(i)*x(i) + b(i)
          else
            fjac(i,j) = x(i)*imw(j) - ct(i)*imw(j)
          end if
        end do
      end do

      return
    end subroutine fdjac

    function fmin(x,fvec,n,ct,cs,imw,m)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Adopted from Numerical Recipes in FORTRAN, Chapter 9.7, 2nd ed.   !
      !bs                                                                    !
      !bs  Returns f = 0.5 * F*F at X. SR FUNCV(N,X,F) is a fixed name,      !
      !bs  user-supplied routine that returns the vector of functions at X.  !
      !bs  The common block NEWTV communicates the function values back to   !
      !bs  NEWT.                                                             !
      !bs                                                                    !
      !bs  Called by:       NEWT                                             !
      !bs                                                                    !
      !bs  Calls:           FUNCV                                            !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: n
      real(rk8) :: ct(np)
      real(rk8) :: cs(np)
      real(rk8) :: imw(np)
      real(rk8) :: m,fmin
      real(rk8) :: x(*), fvec(np)

      integer(ik4) :: i
      real(rk8) :: sum

      call funcv(n,x,fvec,ct,cs,imw,m)
      sum = 0.
      do i = 1, n
        sum = sum + fvec(i)**2
      end do
      fmin = 0.5*sum
      return
    end function fmin

    subroutine funcv(n,x,fvec,ct,cs,imw,m)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Called by:       FMIN                                             !
      !bs                                                                    !
      !bs  Calls:           None                                             !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: n
      real(rk8) :: x(*) , fvec(n)
      real(rk8) :: ct(np) , cs(np) , imw(np) , m
      integer(ik4) :: i, j
      real(rk8) :: sum_jnei , a(np) , b(np) , c(np)

      do i = 1, n
        a(i) = imw(i)
        sum_jnei = 0.
        do j = 1, n
          sum_jnei = sum_jnei + x(j)*imw(j)
        end do
        sum_jnei = sum_jnei - (x(i)*imw(i))
        b(i) = sum_jnei + cs(i)*imw(i) - ct(i)*imw(i)
        c(i) = -ct(i)*(sum_jnei+m)
        fvec(i) = a(i)*x(i)**2 + b(i)*x(i) + c(i)
      end do

      return
    end subroutine funcv


    real(rk8) function getaf(ni,nj,dgni,dgnj,xlsgi,xlsgj,sqrt2)
      ! *** set up new processor for renaming of particles from i to j modes
      real(rk8) :: aa, bb, cc, disc, qq, alfa, l, yji
      real(rk8) :: ni, nj, dgni, dgnj, xlsgi, xlsgj, sqrt2

      alfa = xlsgi/xlsgj
      yji = log(dgnj/dgni)/(sqrt2*xlsgi)
      aa = 1.0 - alfa*alfa
      l = log(alfa*nj/ni)
      bb = 2.0*yji*alfa*alfa
      cc = l - yji*yji*alfa*alfa
      disc = bb*bb - 4.0*aa*cc
      if (disc<0.0) then
        getaf = - 5.0 ! error in intersection                     
        return
      end if
      qq = -0.5*(bb+sign(1.0,bb)*sqrt(disc))
      getaf = cc/qq
      return
      ! *** subroutine to implement Kulmala, Laaksonen, Pirjola
    end function getaf
!     Parameterization for sulfuric acid/water
!     nucleation rates, J. Geophys. Research (103), pp 8301-8307,
!     April 20, 1998.

!ia rev01 27.04.99 changes made to calculation of MdoT see RBiV p.2f
!ia rev02 27.04.99 security check on MdoT > SO4RAT

    subroutine lnsrch(ctot,n,xold,fold,g,p,x,f,stpmax,check,func, &
     fvec,ct,cs,imw,m)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Adopted from Numerical Recipes in FORTRAN, Chapter 9.7, 2nd ed.   !
      !bs                                                                    !
      !bs  Given an n-dimensional point XOLD(1:N), the value of the function !
      !bs  and gradient there, FOLD and G(1:N), and a direction P(1:N),      !
      !bs  finds a new point X(1:N) along the direction P from XOLD where    !
      !bs  the function FUNC has decreased 'sufficiently'. The new function  !
      !bs  value is returned in F. STPMAX is an input quantity that limits   !
      !bs  the length of the steps so that you do not try to evaluate the    !
      !bs  function in regions where it is undefined or subject to overflow. !
      !bs  P is usually the Newton direction. The output quantity CHECK is   !
      !bs  false on a normal; exit. It is true when X is too close to XOLD.  !
      !bs  In a minimization algorithm, this usually signals convergence and !
      !bs  can be ignored. However, in a zero-finding algorithm the calling  !
      !bs  program should check whether the convergence is spurious.         !
      !bs                                                                    !
      !bs  Called by:       NEWT                                             !
      !bs                                                                    !
      !bs  Calls:           FUNC                                             !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: n
      logical :: check
      real(rk8) :: f, fold, stpmax
      real(rk8) :: g(n), p(n), x(n), xold(n)
      real(rk8) :: func
      real(rk8) :: ctot(n)
      real(rk8), parameter :: alf=1.E-04
      real(rk8) :: ct(np)
      real(rk8) :: cs(np)
      real(rk8) :: imw(np)
      real(rk8) :: fvec(n)
      real(rk8) :: m

      external func

      integer(ik4) :: i
      real(rk8) :: a, alam, alam2, alamin, b, disc
      real(rk8) :: f2, fold2, rhs1, rhs2, slope
      real(rk8) :: sum, temp, test, tmplam

      check = .FALSE.
      sum = 0.
      do i = 1, n
        sum = sum + p(i)*p(i)
      end do
      sum = sqrt(sum)
      if (sum>stpmax) then
        do i = 1, n
          p(i) = p(i)*stpmax/sum
        end do
      end if
      slope = 0.
      do i = 1, n
        slope = slope + g(i)*p(i)
      end do
      test = 0.
      do i = 1, n
        temp = abs(p(i))/max(abs(xold(i)),1.)
        if (temp>test) test = temp
      end do
      alamin = tolx/test
      alam = 1.

10    continue

      !bs * avoid negative concentrations and set upper limit given by CTOT.
      do i = 1, n
        x(i) = xold(i) + alam*p(i)
        if (x(i)<=0.) x(i) = conmin
        if (x(i)>ctot(i)) x(i) = ctot(i)
      end do
      f = func(x,fvec,n,ct,cs,imw,m)
      if (alam<alamin) then
        do i = 1, n
          x(i) = xold(i)
        end do
        check = .TRUE.
        return
      else if (f<=fold+alf*alam*slope) then
        return
      else
        if (alam==1.) then
          tmplam = -slope/(2.*(f-fold-slope))
        else
          rhs1 = f - fold - alam*slope
          rhs2 = f2 - fold2 - alam2*slope
          a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
          if (a==0.) then
            tmplam = -slope/(2.*b)
          else
            disc = b*b - 3.*a*slope
            tmplam = (-b+sqrt(disc))/(3.*a)
          end if
          if (tmplam>0.5*alam) tmplam = 0.5*alam
        end if
      end if
      alam2 = alam
      f2 = f
      fold2 = fold
      alam = max(tmplam,0.1*alam)
      go to 10

    end subroutine lnsrch

    subroutine lubksb(a,n,np,indx,b)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Adopted from Numerical Recipes in FORTRAN, Chapter 2.3, 2nd ed.   !
      !bs                                                                    !
      !bs  Solves the set of N linear equations A * X = B. Here A is input,  !
      !bs  not as the matrix A but rather as its LU decomposition,           !
      !bs  determined by the routine LUDCMP. B(1:N) is input as the right-   !
      !bs  hand side vector B, and returns with the solution vector X. A, N, !
      !bs  NP, and INDX are not modified by this routine and can be left in  !
      !bs  place for successive calls with different right-hand sides B.     !
      !bs  This routine takes into account the possibilitythat B will begin  !
      !bs  with many zero elements, so it is efficient for use in matrix     !
      !bs  inversion.                                                        !
      !bs                                                                    !
      !bs  Called by:       NEWT                                             !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: n, np, indx(n)
      real(rk8) :: a(np,np), b(n)

      integer(ik4) :: i, ii, j, ll
      real(rk8) :: sum

      ii = 0
      do i = 1, n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii/=0) then
          do j = ii, i - 1
            sum = sum - a(i,j)*b(j)
          end do
        else if (sum/=0) then
          ii = i
        end if
        b(i) = sum
      end do
      do i = n, 1, -1
        sum = b(i)
        do j = i + 1, n
          sum = sum - a(i,j)*b(j)
        end do
        b(i) = sum/a(i,i)
      end do

      return
    end subroutine lubksb

    subroutine ludcmp(a,n,np,indx,d)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Adopted from Numerical Recipes in FORTRAN, Chapter 2.3, 2nd ed.   !
      !bs                                                                    !
      !bs  Equation (2.3.14) Numerical Recipes, p 36:                        !
      !bs   | b_11 b_12 b_13 b_14 |                                          !
      !bs   | a_21 b_22 b_23 b_24 |                                          !
      !bs   | a_31 a_32 b_33 b_34 |                                          !
      !bs   | a_41 a_42 a_43 b_44 |                                          !
      !bs                                                                    !
      !bs  Given a matrix A(1:N,1:N), with physical dimension NP by NP, this !
      !bs  routine replaces it by the LU decomposition of a rowwise          !
      !bs  permutation of itself. A and N are input. A is output arranged as !
      !bs  in equation (2.3.14) above; INDX(1:N) is an output vector that    !
      !bs  records vector that records the row permutation effected by the   !
      !bs  partial pivoting; D is output as +-1 depending on whether the     !
      !bs  number of row interchanges was even or odd, respectively. This    !
      !bs  routine is used in combination with SR LUBKSB to solve linear     !
      !bs  equations or invert a matrix.                                     !
      !bs                                                                    !
      !bs  Called by:       NEWT                                             !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

      integer(ik4) :: n, np, indx(n)
      integer(ik4), parameter :: nmax=10 !largest expected N
      real(rk8) :: d, a(np,np)
      real(rk8), parameter :: tiny=1.0E-20

      integer(ik4) :: i, imax, j, k
      real(rk8) :: aamax, dum, sum, vv(nmax)

      d = 1
      do i = 1, n
        aamax = 0.
        do j = 1, n
          if (abs(a(i,j))>aamax) aamax = abs(a(i,j))
        end do
        if (aamax==0) then
          a(1,1)=epsilc
        end if
        vv(i) = 1./aamax
      end do
      do j = 1, n
        do i = 1, j - 1
          sum = a(i,j)
          do k = 1, i - 1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
        end do
        aamax = 0.
        do i = j, n
          sum = a(i,j)
          do k = 1, j - 1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
          dum = vv(i)*abs(sum)
          if (dum>=aamax) then
            imax = i
            aamax = dum
          end if
        end do
        if (j/=imax) then
          do k = 1, n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vv(imax) = vv(j)
        end if
        indx(j) = imax
        if (a(j,j)==0.) a(j,j) = tiny
        if (j/=n) then
          dum = 1./a(j,j)
          do i = j + 1, n
            a(i,j) = a(i,j)*dum
          end do
        end if
      end do

      return
    end subroutine ludcmp

    subroutine newt(x,n,check,ctot,csat,imwcv,minitw,its)
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      !bs                                                                    !
      !bs  Description:                                                      !
      !bs                                                                    !
      !bs  Adopted from Numerical Recipes in FORTRAN, Chapter 9.7, 2nd ed.   !
      !bs                                                                    !
      !bs  Given an initial guess X(1:N) for a root in N dimensions, find    !
      !bs  the root by globally convergent Newton's method. The vector of    !
      !bs  functions to be zeroed, called FVEC(1:N) in the routine below. is !
      !bs  retuned by a user-supplied function that must be called FUNCV and !
      !bs  have the declaration subroutine FUNCV(NX,FVEC). The output        !
      !bs  quantity CHECK is false on a normal return and true if the        !
      !bs  routine has converged to a local minimum of the function FMIN     !
      !bs  defined below. In this case try restarting from a different       !
      !bs  initial guess.                                                    !
      !bs                                                                    !
      !bs  PARAMETERS                                                        !
      !bs  NP     : maximum expected value of N                              !
      !bs  MAXITS : maximum number of iterations                             !
      !bs  TOLF   : convergence criterion on function values                 !
      !bs  TOLMIN : criterion for decidingwhether spurios convergence to a   !
      !bs           minimum of FMIN has ocurred                              !
      !bs  TOLX   : convergence criterion on delta_X                         !
      !bs  STPMX  : scaled maximum step length allowed in line searches      !
      !bs                                                                    !
      !bs  Called by:       SOA_PART                                         !
      !bs                                                                    !
      !bs  Calls:           FDJAC                                            !
      !bs                   FMIN                                             !
      !bs                   LNSRCH                                           !
      !bs                   LUBKSB                                           !
      !bs                   LUDCMP                                           !
      !bs                                                                    !
      !bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
      
      integer(ik4) :: n     !bs dimension of problem                  
      real(rk8) :: x(n)     !bs initial guess of CAER                 
      logical :: check
      real(rk8) :: ctot(n)  !bs total concentration GAS + AER + PROD  
      real(rk8) :: csat(n)  !bs saturation conc. of cond. vapor [ug/m^
      real(rk8) :: imwcv(n) !bs inverse molecular weights             
      real(rk8) :: minitw   !bs weighted initial mass

      !bs * following Numerical recipes
      integer(ik4) :: nn
      real(rk8) :: fvec(np) !bs vector of functions to be zeroed
      real(rk8) :: ct(np)
      real(rk8) :: cs(np)
      real(rk8) :: imw(np)
      real(rk8) :: m

      integer(ik4) :: i, its, j, indx(np)
      real(rk8) :: d, den, f, fold, stpmax, sum, temp, test
      real(rk8) :: fjac(np,np)
      real(rk8) :: g(np), p(np), xold(np)

      m = minitw
      do i = 1, n
        ct(i) = ctot(i)
        cs(i) = csat(i)
        imw(i) = imwcv(i)
      end do

      nn = n
      f = fmin(x,fvec,nn,ct,cs,imw,m) !The vector FVEC is 
      test = 0.    !Test for initial guess being a root. Us
      do i = 1, n  !stringent test than simply TOLF.       
        if (abs(fvec(i))>test) test = abs(fvec(i))
      end do
      if (test<0.01*tolf) return
      sum = 0.     !Calculate STPMAX for line searches     
      do i = 1, n
        sum = sum + x(i)**2
      end do
      stpmax = stpmx*max(sqrt(sum),float(n))
      do its = 1, maxits !start of iteration loop   
        call fdjac(n,x,fjac,ct,cs,imw) !get Jacobian              
        do i = 1, n !compute Delta f for line search 
          sum = 0.
          do j = 1, n
            sum = sum + fjac(j,i)*fvec(j)
          end do
          g(i) = sum
        end do
        do i = 1, n !store X 
          xold(i) = x(i)
        end do
        fold = f !store F 
        do i = 1, n !right-hand side for linear equations
          p(i) = -fvec(i)
        end do
        call ludcmp(fjac,n,np,indx,d)     !solve linear equations by LU dec
        call lubksb(fjac,n,np,indx,p)
        call lnsrch(ctot,n,xold,fold,g, & !LNSRCH returns new X and F. It a
          p,x,f,stpmax,check,fmin,fvec, & !calculates FVEC at the new X whe
          ct,cs,imw,m)                    !calls FMIN                      
        test = 0.
        do i = 1, n
          if (abs(fvec(i))>test) test = abs(fvec(i))
        end do
        if (test<tolf) then
          check = .FALSE.
          return
        end if
        if (check) then !Check for gradient of F zero,          
          test = 0. !i.e., superious convergence.  
          den = max(f,0.5*n)
          do i = 1, n
            temp = abs(g(i))*max(abs(x(i)),1.)/den
            if (temp>test) test = temp
          end do
          if (test<tolmin) then
            check = .TRUE.
          else
            check = .FALSE.
          end if
          return
        end if
        test = 0. !Test for convergence on delta_x  
        do i = 1, n
          temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.)
          if (temp>test) test = temp
        end do
        if (test<tolx) return
      end do

    end subroutine newt

end module mod_che_aerosols_sorgam
