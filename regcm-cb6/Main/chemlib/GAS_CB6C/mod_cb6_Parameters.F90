! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : mod_cb6_Parameters.f90
! Time                 : Thu Aug 31 17:04:30 2017
! Working directory    : /home/regcm/Code/kpp-2.2.3/Workspace/GAS_CB6r2_mod
! Equation file        : cb6.kpp
! Output root filename : cb6
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mod_cb6_Parameters

  USE mod_cb6_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
! INTEGER, PARAMETER :: NSPEC = 81 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR_CB6 = 76 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 74 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX_CB6 = 5 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 218 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 77 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 938 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 1077 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 77 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 1 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 1 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(indv_spc) = C(ind_spc)

  INTEGER, PARAMETER :: indv_SULF = 1 
  INTEGER, PARAMETER :: indv_CDOX = 2 
  INTEGER, PARAMETER :: indv_SDIO = 3 
  INTEGER, PARAMETER :: indv_OSNG = 4 
  INTEGER, PARAMETER :: indv_ETHY = 5 
  INTEGER, PARAMETER :: indv_DNPO = 6 
  INTEGER, PARAMETER :: indv_BENZ = 7 
  INTEGER, PARAMETER :: indv_EPOX = 8 
  INTEGER, PARAMETER :: indv_ETOH = 9 
  INTEGER, PARAMETER :: indv_KET = 10 
  INTEGER, PARAMETER :: indv_TOLN = 11 
  INTEGER, PARAMETER :: indv_XYLN = 12 
  INTEGER, PARAMETER :: indv_HPLD = 13 
  INTEGER, PARAMETER :: indv_PACN = 14 
  INTEGER, PARAMETER :: indv_PACD = 15 
  INTEGER, PARAMETER :: indv_PNA = 16 
  INTEGER, PARAMETER :: indv_HONO = 17 
  INTEGER, PARAMETER :: indv_MTHA = 18 
  INTEGER, PARAMETER :: indv_NTR2 = 19 
  INTEGER, PARAMETER :: indv_MEPX = 20 
  INTEGER, PARAMETER :: indv_HPOX = 21 
  INTEGER, PARAMETER :: indv_OPAN = 22 
  INTEGER, PARAMETER :: indv_CAT1 = 23 
  INTEGER, PARAMETER :: indv_ISPX = 24 
  INTEGER, PARAMETER :: indv_ETHA = 25 
  INTEGER, PARAMETER :: indv_FACD = 26 
  INTEGER, PARAMETER :: indv_PANX = 27 
  INTEGER, PARAMETER :: indv_HCO3 = 28 
  INTEGER, PARAMETER :: indv_PRPA = 29 
  INTEGER, PARAMETER :: indv_MEOH = 30 
  INTEGER, PARAMETER :: indv_CRER = 31 
  INTEGER, PARAMETER :: indv_RPOX = 32 
  INTEGER, PARAMETER :: indv_NTR1 = 33 
  INTEGER, PARAMETER :: indv_ACET = 34 
  INTEGER, PARAMETER :: indv_INTR = 35 
  INTEGER, PARAMETER :: indv_AACD = 36 
  INTEGER, PARAMETER :: indv_CRON = 37 
  INTEGER, PARAMETER :: indv_ROR = 38 
  INTEGER, PARAMETER :: indv_BZO2 = 39 
  INTEGER, PARAMETER :: indv_ETHE = 40 
  INTEGER, PARAMETER :: indv_CMON = 41 
  INTEGER, PARAMETER :: indv_TOLR = 42 
  INTEGER, PARAMETER :: indv_XYLR = 43 
  INTEGER, PARAMETER :: indv_CRSL = 44 
  INTEGER, PARAMETER :: indv_TERP = 45 
  INTEGER, PARAMETER :: indv_ISPR = 46 
  INTEGER, PARAMETER :: indv_NTRC = 47 
  INTEGER, PARAMETER :: indv_EPX2 = 48 
  INTEGER, PARAMETER :: indv_GLYD = 49 
  INTEGER, PARAMETER :: indv_GLY = 50 
  INTEGER, PARAMETER :: indv_XOPN = 51 
  INTEGER, PARAMETER :: indv_MEGY = 52 
  INTEGER, PARAMETER :: indv_ROPN = 53 
  INTEGER, PARAMETER :: indv_IOLE = 54 
  INTEGER, PARAMETER :: indv_OLE = 55 
  INTEGER, PARAMETER :: indv_FORM = 56 
  INTEGER, PARAMETER :: indv_AALD = 57 
  INTEGER, PARAMETER :: indv_XO2R = 58 
  INTEGER, PARAMETER :: indv_OPO3 = 59 
  INTEGER, PARAMETER :: indv_XO2N = 60 
  INTEGER, PARAMETER :: indv_ISO2 = 61 
  INTEGER, PARAMETER :: indv_XO2H = 62 
  INTEGER, PARAMETER :: indv_ISPD = 63 
  INTEGER, PARAMETER :: indv_MEO2 = 64 
  INTEGER, PARAMETER :: indv_ALKA = 65 
  INTEGER, PARAMETER :: indv_ALDX = 66 
  INTEGER, PARAMETER :: indv_HOX = 67 
  INTEGER, PARAMETER :: indv_ROO = 68 
  INTEGER, PARAMETER :: indv_NTOX = 69 
  INTEGER, PARAMETER :: indv_O = 70 
  INTEGER, PARAMETER :: indv_POX = 71 
  INTEGER, PARAMETER :: indv_CXO3 = 72 
  INTEGER, PARAMETER :: indv_ACOO = 73 
  INTEGER, PARAMETER :: indv_OZN = 74 
  INTEGER, PARAMETER :: indv_NMOX = 75 
  INTEGER, PARAMETER :: indv_NDOX = 76 

! Index declaration for fixed species in C
!   C(indv_spc)

  INTEGER, PARAMETER :: indv_WTR = 77 
  INTEGER, PARAMETER :: indv_DIHY = 78 
  INTEGER, PARAMETER :: indv_O2 = 79 
  INTEGER, PARAMETER :: indv_M = 80 
  INTEGER, PARAMETER :: indv_DUMMY2 = 81 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(indv_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_WTR = 1 
  INTEGER, PARAMETER :: indf_DIHY = 2 
  INTEGER, PARAMETER :: indf_O2 = 3 
  INTEGER, PARAMETER :: indf_M = 4 
  INTEGER, PARAMETER :: indf_DUMMY2 = 5 

END MODULE mod_cb6_Parameters

