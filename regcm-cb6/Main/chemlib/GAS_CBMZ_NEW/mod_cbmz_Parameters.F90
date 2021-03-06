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
! File                 : mod_cbmz_Parameters.f90
! Time                 : Mon Nov 25 13:41:15 2013
! Working directory    : /scratch/ashalaby/kpp-2.2.3/compare/CBMZ
! Equation file        : mod_cbmz.kpp
! Output root filename : mod_cbmz
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mod_cbmz_Parameters

  USE mod_cbmz_Precision
  PUBLIC
  SAVE

! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 58
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 52
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 124
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 59
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 531
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 600
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 59
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 4
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 2
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: indv_CO2 = 1
  INTEGER, PARAMETER :: indv_H2SO4 = 2
  INTEGER, PARAMETER :: indv_HCOOH = 3
  INTEGER, PARAMETER :: indv_RCOOH = 4
  INTEGER, PARAMETER :: indv_MSA = 5
  INTEGER, PARAMETER :: indv_DUMMY = 6
  INTEGER, PARAMETER :: indv_PAN = 7
  INTEGER, PARAMETER :: indv_TOL = 8
  INTEGER, PARAMETER :: indv_O1D = 9
  INTEGER, PARAMETER :: indv_H2O2 = 10
  INTEGER, PARAMETER :: indv_SO2 = 11
  INTEGER, PARAMETER :: indv_XYL = 12
  INTEGER, PARAMETER :: indv_CH4 = 13
  INTEGER, PARAMETER :: indv_C2H6 = 14
  INTEGER, PARAMETER :: indv_CRO = 15
  INTEGER, PARAMETER :: indv_DMS = 16
  INTEGER, PARAMETER :: indv_HNO4 = 17
  INTEGER, PARAMETER :: indv_H2 = 18
  INTEGER, PARAMETER :: indv_TO2 = 19
  INTEGER, PARAMETER :: indv_CH3OH = 20
  INTEGER, PARAMETER :: indv_HNO2 = 21
  INTEGER, PARAMETER :: indv_CH3OOH = 22
  INTEGER, PARAMETER :: indv_ETHOOH = 23
  INTEGER, PARAMETER :: indv_N2O5 = 24
  INTEGER, PARAMETER :: indv_ETH = 25
  INTEGER, PARAMETER :: indv_CRES = 26
  INTEGER, PARAMETER :: indv_O3P = 27
  INTEGER, PARAMETER :: indv_CO = 28
  INTEGER, PARAMETER :: indv_HNO3 = 29
  INTEGER, PARAMETER :: indv_PAR = 30
  INTEGER, PARAMETER :: indv_OPEN = 31
  INTEGER, PARAMETER :: indv_ISOPN = 32
  INTEGER, PARAMETER :: indv_ISOPP = 33
  INTEGER, PARAMETER :: indv_ISOPO2 = 34
  INTEGER, PARAMETER :: indv_H2O = 35
  INTEGER, PARAMETER :: indv_AONE = 36
  INTEGER, PARAMETER :: indv_OLEI = 37
  INTEGER, PARAMETER :: indv_ISOP = 38
  INTEGER, PARAMETER :: indv_HCHO = 39
  INTEGER, PARAMETER :: indv_OLET = 40
  INTEGER, PARAMETER :: indv_XO2 = 41
  INTEGER, PARAMETER :: indv_MGLY = 42
  INTEGER, PARAMETER :: indv_ETHP = 43
  INTEGER, PARAMETER :: indv_NAP = 44
  INTEGER, PARAMETER :: indv_ALD2 = 45
  INTEGER, PARAMETER :: indv_CH3O2 = 46
  INTEGER, PARAMETER :: indv_ISOPRD = 47
  INTEGER, PARAMETER :: indv_ANO2 = 48
  INTEGER, PARAMETER :: indv_ROOH = 49
  INTEGER, PARAMETER :: indv_RO2 = 50
  INTEGER, PARAMETER :: indv_ONIT = 51
  INTEGER, PARAMETER :: indv_HO2 = 52
  INTEGER, PARAMETER :: indv_O3 = 53
  INTEGER, PARAMETER :: indv_OH = 54
  INTEGER, PARAMETER :: indv_NO = 55
  INTEGER, PARAMETER :: indv_NO2 = 56
  INTEGER, PARAMETER :: indv_NO3 = 57
  INTEGER, PARAMETER :: indv_C2O3 = 58

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: indv_O2 = 59
  INTEGER, PARAMETER :: indv_N2 = 60

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_O2 = 1
  INTEGER, PARAMETER :: indf_N2 = 2

END MODULE mod_cbmz_Parameters

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
