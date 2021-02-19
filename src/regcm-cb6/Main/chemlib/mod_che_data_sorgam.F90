module mod_che_data_sorgam

  use mod_intkinds
  use mod_realkinds

  implicit none

  real(rk8), parameter :: conmin  = 1.E-16 ! conc lower limit [ug/m^3]
  real(rk8), parameter :: epsilc  = conmin 
! real(rk8), parameter :: dgmin   = 1.E-09 ! lowest particle diameter ( m )
  real(rk8), parameter :: dgminn  = 0.01E-6! init mean diam for nuclei mode [m]
  real(rk8), parameter :: dgmina  = 0.07E-6! init mean diam for accu mode [m]
  real(rk8), parameter :: densmin = 0.626D3! lowest particle density [kg/m^3]
  real(rk8), parameter :: aeroconcmin = 0.0001 ! minimum aerosol sulfate conc

  real(rk8) :: pdensn ! average particle density in nuclei
  real(rk8) :: pdensa ! average particle density in accumu
  real(rk8) :: dgnuc  ! nuclei mode geometric mean diamete
  real(rk8) :: dgacc  ! accumulation geometric mean diamet
  integer(ik4), parameter :: np = 8       !bs maximum expected value of N
  integer(ik4), parameter :: maxits = 100 !bs maximum number of iterations

  real(rk8), parameter :: tolf   = 1.E-09 !convergence criterion on func vals
  real(rk8), parameter :: tolmin = 1.E-12 !criterion whether superio   
                                          !convergence to a min of fmin occurred
  real(rk8), parameter :: tolx   = 1.E-10 !convergence criterion on delta_x
  real(rk8), parameter :: stpmx  = 100.   !scaled maximum step length allowed

  ! //////////////////////////////////////////////////////////////////////
  ! FSB include file
  !
  ! *** declare and set flag for organic aerosol production method
  ! *** Two method are available:
  !
  ! *** The method of Pandis,Harley, Cass, and Seinfeld, 1992,
  !     Secondary aerosol formation and transport, Atmos. Environ., 26A,
  !     pp 2453-2466
  !     Bowman et al. Atmospheric Environment
  !     Vol 29, pp 579-589, 1995.
  ! *** and
  ! *** The method of Odum, Hoffmann, Bowman, Collins, Flagen and
  !     Seinfeld, 1996, Gas/particle partitioning and secondary organic ae
  !     yields, Environ. Sci, Technol, 30, pp 2580-2585.

! integer(ik4) , parameter :: orgaer = 2
  ! 1 = Pandis et al.  1992 method is used
  ! 2 = Pankow 1994/Odum et al. 1996 method is
  ! switch for organic aerosol method

  ! *** set up indices for array  CBLK
  !integer(ik4) , parameter :: vorgpaj=27! Accumulation mode primary anth
  !integer(ik4) , parameter :: vorgpai=28! Aitken mode primary anth

  !BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * !
  !BS *                                                                * !
  !BS *            include file used in SORGAM routines                * !
  !BS *                                                                * !
  !BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * !

  !bs * species pointer for SOA species
  !integer(ik4) , parameter :: psoaaro1=1
  !integer(ik4) , parameter :: psoaaro2=2
  !integer(ik4) , parameter :: psoaalk1=3
  !integer(ik4) , parameter :: psoaole1=4
  !integer(ik4) , parameter :: psoaapi1=5
  !integer(ik4) , parameter :: psoaapi2=6
  !integer(ik4) , parameter :: psoalim1=7
  !integer(ik4) , parameter :: psoalim2=8

  ! *** include file for aerosol routines
  !....................................................................
  !  CONTAINS: Fundamental constants for air quality modeling
  !
  !  DEPendENT UPON:  none
  ! 
  !  REVISION HISTORY:
  !
  !    Adapted 6/92 by CJC from ROM's PI.EXT.
  !
  !    Revised 3/1/93 John McHenry to include constants needed by
  !    LCM aqueous chemistry
  !    Revised 9/93 by John McHenry to include additional constants
  !    needed for FMEM clouds and aqueous chemistry
  !
  !    Revised 3/4/96 by Dr. Francis S. Binkowski to reflect current
  !    Models3 view that MKS units should be used wherever possible,
  !    and that sources be documentated. Some variables have been added
  !    names changed, and values revised.
  !
  !    Revised 3/7/96 to have universal gas constant input and compute
  !    gas constant is chemical form. TWOPI is now calculated rather than
  !
  !    Revised 3/13/96 to group declarations and parameter statements.
  !
  !    Revised 9/13/96 to include more physical constants.
  !    Revised 12/24/96 eliminate silly EPSILON, AMISS
  !
  !    Revised 1/06/97 to eliminate most derived constants
  !
  !    Revised 16/12/15 removed from WRF and introduced to RegCM4
  !
  ! FSB REFERENCES:
  !
  !      CRC76,        CRC Handbook of Chemistry and Physics (76th Ed),
  !                     CRC Press, 1995
  !      Hobbs, P.V.   Basic Physical Chemistry for the Atmospheric Scien
  !                     Cambridge Univ. Press, 206 pp, 1995.
  !      Snyder, J.P., Map Projections-A Working Manual, U.S. Geological
  !                     Paper 1395 U.S.GPO, Washington, DC, 1987.
  !      Stull, R. B., An Introduction to Bounday Layer Meteorology, Klu
  !                     Dordrecht, 1988

  ! Geometric Constants:
  real(rk8) , parameter :: pi   = 4.D+00*datan(1.D+00)! PI
  ! Fundamental Constants: ( Source: CRC76, pp 1-1 to 1-6)
  real(rk8) , parameter :: avo  = 6.0221367E23        ! Avogadro's Const [1/mol]
  real(rk8) , parameter :: rgas = 8.314510            ! univ gas const [J/mol-K]

  ! Atmospheric Constants:
  ! FSB                     78.06%  N2, 21% O2 and 0.943% A on a mole
  ! fraction basis. ( Source : Hobbs, 1995) pp 69-
  real(rk8) , parameter :: mwair   = 28.9628 ! mean molec wght dry air [g/mol]
  real(rk8) , parameter :: rdgas   = 1.0E3*rgas/mwair     ! gas const [J/kg-K]
  real(rk8) , parameter :: f6dpi   = 6.0/pi               ! 6/PI
  real(rk8) , parameter :: f6dpi9  = 1.0E9*f6dpi          ! 1.0e9 * 6/PI
  real(rk8) , parameter :: f6dpim9 = 1.0E-9*f6dpi         ! 1.0e-9 * 6/PI
  real(rk8) , parameter :: sqrt2   = sqrt(2.0)            ! SQRT(2)
  real(rk8) , parameter :: one3    = 1.0/3.0              ! 1/3
  real(rk8) , parameter :: two3    = 2.0/3.0              ! 2/3
  ! *** physical constants:
  real(rk8) , parameter :: boltz  = rgas/avo ! Boltzmann's Const [J/K]
  ! *** bulk density aerosol components [ kg/m^3 ] :
  real(rk8) , parameter :: rhoso4  = 1.8E3   ! sulfate
  real(rk8) , parameter :: rhoorg  = 1.0E3   ! organics
  real(rk8) , parameter :: rhocnyl = 0.791D3 ! carbonyls 
  real(rk8) , parameter :: rhosthc = 0.626D3 ! alkanes and alcohols
  real(rk8) , parameter :: rhoolef = 0.858D3 ! olefins
  real(rk8) , parameter :: rhoarom = 1.371D3 ! aromatics
  real(rk8) , parameter :: rhontro = 1.163D3 ! nitro organic
  real(rk8) , parameter :: rhomfnc = 1.265D3 ! multifunctionals

  ! *** Factors for converting aerosol mass concentration [ug m^-3] to
  !      	  to 3rd moment concentration [m^3 m^-3]
  real(rk8) , parameter :: so4fac = f6dpim9/rhoso4 ![ug^-1 m^3]
  real(rk8) , parameter :: orgfac = f6dpim9/rhoorg
  real(rk8) :: aerfact

  real(rk8) , parameter :: pss0   = 101325.0 ! starting stand surf press [Pa]
  real(rk8) , parameter :: tss0   = 288.15   ! starting stand surf temp [K]
  real(rk8) , parameter :: sginin = 1.70     ! init sigma-G for nucleimode
  real(rk8) , parameter :: sginia = 2.00     ! init sigma-G for accumulation mode
  real(rk8) , parameter :: dginin = dgminn   ! init mean diam for nuclei mode [m]
  real(rk8) , parameter :: dginia = dgmina   ! init mean diam for accu mode [m]

  !     LOGICAL diagnostics
  ! *** Scalar variables for fixed standard deviations.

  ! Flag for writing diagnostics to file       
  real(rk8), parameter :: xxlsgn   = log(sginin)      ! log(sginin)
  real(rk8), parameter :: xxlsga   = log(sginia)      ! log(sginia)
  real(rk8), parameter :: l2sginin = xxlsgn**2        ! log(sginin ) ** 2
  real(rk8), parameter :: l2sginia = xxlsga**2        ! log(sginia ) ** 2
  real(rk8), parameter :: en1   = exp(0.125*l2sginin) ! exp(log^2(sigmag)/8 )
  real(rk8), parameter :: ea1   = exp(0.125*l2sginia) ! exp(log^2(sigmag)
  real(rk8), parameter :: esn04 = en1**4              ! **4 
  real(rk8), parameter :: esa04 = ea1**4  
  real(rk8), parameter :: esn05 = esn04*en1           ! **5
  real(rk8), parameter :: esa05 = esa04*ea1  
  real(rk8), parameter :: esn08 = esn04*esn04         ! **8
  real(rk8), parameter :: esa08 = esa04*esa04  
  real(rk8), parameter :: esn09 = esn04*esn05         ! **9 
  real(rk8), parameter :: esa09 = esa04*esa05 
  real(rk8), parameter :: esn12 = esn04*esn04*esn04   ! **12
  real(rk8), parameter :: esa12 = esa04*esa04*esa04  
  real(rk8), parameter :: esn16 = esn08*esn08         ! **16
  real(rk8), parameter :: esa16 = esa08*esa08  
  real(rk8), parameter :: esn20 = esn16*esn04         ! **20
  real(rk8), parameter :: esa20 = esa16*esa04   
  real(rk8), parameter :: esn24 = esn12*esn12         ! **24
  real(rk8), parameter :: esa24 = esa12*esa12  
  real(rk8), parameter :: esn25 = esn16*esn09         ! **25
  real(rk8), parameter :: esa25 = esa16*esa09  
  real(rk8), parameter :: esn28 = esn20*esn08         ! **28
  real(rk8), parameter :: esa28 = esa20*esa08  
  real(rk8), parameter :: esn32 = esn16*esn16         ! **32
  real(rk8), parameter :: esa32 = esa16*esa16  
  real(rk8), parameter :: esn36 = esn16*esn20         ! **36
  real(rk8), parameter :: esa36 = esa16*esa20  
  real(rk8), parameter :: esn49 = esn25*esn20*esn04   ! **49
  real(rk8), parameter :: esa49 = esa25*esa20*esa04  
  real(rk8), parameter :: esn52 = esn16*esn36         ! **52
  real(rk8), parameter :: esa52 = esa16*esa36  
  real(rk8), parameter :: esn64 = esn32*esn32         ! **64
  real(rk8), parameter :: esa64 = esa32*esa32  
  real(rk8), parameter :: esn100 = esn36*esn64        ! **100
  real(rk8), parameter :: esnm20 = 1.0/esn20          ! **(-20)
  real(rk8), parameter :: esam20 = 1.0/esa20 
  real(rk8), parameter :: esnm32 = 1.0/esn32          ! **(-32)
  real(rk8), parameter :: esam32 = 1.0/esa32 

  ! *** set up COMMON blocks for esg's:
  ! *** SET NUCLEATION FLAG:
  ! Flag for Choice of nucleation Mechanism
  !integer(ik4) , parameter :: inucl=2
  ! INUCL = 0, Kerminen & Wexler Mechanism
  ! INUCL = 1, Youngblood and Kreidenweis mech
  ! INUCL = 2, Kulmala et al. mechanism

  ! special factor to compute mass transfer
  real(rk8), parameter :: dgvem_i = 0.03E-6 ! Aitken [m]
  real(rk8), parameter :: dgvem_j = 0.3E-6  ! accumu [m]

  ! *** factors for getting number emissions rate from mass emissions rate
  real(rk8), parameter :: facatkn_min = 0.04
  real(rk8), parameter :: facacc_min  = 1.0-facatkn_min
  real(rk8), parameter :: factnumn = exp(4.5*log(sginin)**2)/dgvem_i**3 ! Aitken
  real(rk8), parameter :: factnuma = exp(4.5*log(sginia)**2)/dgvem_j**3 ! accumu
  real(rk8), parameter :: xxm3 = 3.0*xxlsgn/ sqrt2
  ! [ ug/m**3 ] ! changed 1/6/98 
  ! factor to set minimum for aerosol mode number
  real(rk8), parameter :: nummin_i = facatkn_min*so4fac*aeroconcmin/ & 
                                       (dginin**3*esn36) ! Aitken  
  real(rk8), parameter :: nummin_j = facacc_min*so4fac*aeroconcmin/ & 
                                       (dginia**3*esa36) ! accumu
      
  !bs      REAL ALPHSULF ! Accommodation coefficient for sulfuric acid
  !bs      PARAMETER ( ALPHSULF = 0.05 ) ! my be set to one in future
  !bs
  !bs      REAL DIFFSULF ! molecular diffusivity for sulfuric acid [ m**2
  !bs      PARAMETER( DIFFSULF = 0.08E-4 ) ! may be changed in future
  !bs
  !bs * 23/03/99 updates of ALPHSULF and DIFFSULF adopted fro new code fro
  !bs * DIFFSULF is calculated from Reid, Prausnitz, and Poling, The prope
  !bs * of gases and liquids, 4th edition, McGraw-Hill, 1987, pp 587-588.
  !bs * Equation (11-4.4) was used.
  !bs * The value is at T = 273.16 K and P = 1.01325E05 Pa
  !bs * Temperature dependence is included for DIFFSULF via DIFFCORR (see

  ! molecular diffusivity for sulfuric acid [m^2/s]
  real(rk8), parameter :: diffsulf   = 9.362223E-06 
  !bs * DIFFORG is calculated from the same formula as DIFFSULF.
  !bs * An average elemental composition of C=8, O=3, N=1, H=17 is asuumed
  !bs * to calculate DIFFORG at T = 273.16K and  P = 1.01325E05 Pa.
  !bs * Temepratur dependence is included below.
  real(rk8) , parameter :: difforg   = 5.151174E-06 !molec diffusivity for organ
  ! *** CCONC is the factor for near-continuum condensation.
  real(rk8) , parameter :: cconc     = 2.0*pi*diffsulf ! ccofm * sqrt( ta )
  !bs * factor for NC condensation for organics
  real(rk8) , parameter :: cconc_org = 2.0*pi*difforg  ! [m^2/s]
  ! Accommodation coefficient
  real(rk8) , parameter :: alphsulf  = 1.0 !for sulfuric
  real(rk8) , parameter :: alphaorg  = 1.0 !for organic
  ! molecular weights
  real(rk8) , parameter :: mwh2so4   = 98.07354E-3 !for sulfuric [kg/mol]
  real(rk8) , parameter :: mworg     = 175.0E-03   !for organic  [kg/mol]
  !bs analogue to CCOFM but for organics ! [m^2/s]  
  real(rk8) , parameter :: ccofm_org = alphaorg*sqrt(pi*rgas/(2.0*mworg))
  ! FSB  CCOFM is  the accommodation coefficient
  !      times the mean molecular velocity for h2so4 without the temperatu
  !      after some algebra
  !bs CCOFM_ORG * sqrt(TA)                
  real(rk8) , parameter :: ccofm = alphsulf*sqrt(pi*rgas/(2.0*mwh2so4))

end module mod_che_data_sorgam
