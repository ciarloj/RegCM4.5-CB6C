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

module mod_che_molwg

  use mod_intkinds
  use mod_realkinds
  use mod_che_indices

  implicit none

  public

  real(rk8) , parameter :: w_no2 = 46.0D0
  real(rk8) , parameter :: w_no  = 30.0D0

  real(rk8) , parameter :: w_hono = 47.0D0
  real(rk8) , parameter :: w_hno2 = 47.0D0
  real(rk8) , parameter :: w_no3  = 62.0D0
  real(rk8) , parameter :: w_n2o5 = 108.0D0
  real(rk8) , parameter :: w_hno4 = 79.0D0
  real(rk8) , parameter :: w_hno3 = 63.0D0
  real(rk8) , parameter :: w_o3   = 48.0D0
  real(rk8) , parameter :: w_h2o2 = 34.0D0
  real(rk8) , parameter :: w_h2o  = 18.0D0

  real(rk8) , parameter :: w_so2  = 64.0D0
  real(rk8) , parameter :: w_sulf = 98.0D0
  real(rk8) , parameter :: w_h2so4= 98.0D0
  real(rk8) , parameter :: w_co   = 28.0D0
  real(rk8) , parameter :: w_co2  = 44.0D0
  real(rk8) , parameter :: w_h2   = 2.0D0

  real(rk8) , parameter :: w_oh   = 17.0D0
  real(rk8) , parameter :: w_ho2  = 33.0D0
  real(rk8) , parameter :: w_ro2  = 47.0D0
  real(rk8) , parameter :: w_xo2  = 47.0D0

  ! Alkane species in RADM2 and CBMZ
  ! some species are repeate with differen name convention
  real(rk8) , parameter :: w_ch4    = 16.0D0
  real(rk8) , parameter :: w_ethan  = 30.0D0
  real(rk8) , parameter :: w_c2h6   = 30.070D0
  ! assumed molecular wieght for PAR (CBMZ mechanism)
  real(rk8) , parameter :: w_par    = 44.0D0
  real(rk8) , parameter :: w_hc3    = 44.0D0
  real(rk8) , parameter :: w_c3h8   = 44.10D0
  real(rk8) , parameter :: w_hc5    = 72.0D0
  real(rk8) , parameter :: w_hc8    = 114.0D0
  real(rk8) , parameter :: w_alk4   = 58.120D0
  real(rk8) , parameter :: w_alk7   = 100.200D0

  ! Alkene species in RADM2 and CBMZ
  real(rk8) , parameter :: w_ethene = 28.0D0
  real(rk8) , parameter :: w_eth    = 28.0D0
  real(rk8) , parameter :: w_ol2    = 28.0D0
  real(rk8) , parameter :: w_olt    = 42.0D0
  real(rk8) , parameter :: w_oli    = 56.0D0
  real(rk8) , parameter :: w_olet   = 42.0D0
  real(rk8) , parameter :: w_olei   = 56.0D0
  real(rk8) , parameter :: w_prpe   = 42.0D0
  real(rk8) , parameter :: w_bute   = 56.0D0
  real(rk8) , parameter :: w_isop   = 68.0D0

  ! Aromatic
  real(rk8) , parameter :: w_tolu  = 92.0D0
  real(rk8) , parameter :: w_tol   = 92.0D0
  real(rk8) , parameter :: w_csl   = 108.0D0
  real(rk8) , parameter :: w_cres  = 108.0D0
  real(rk8) , parameter :: w_xyle  = 106.0D0
  real(rk8) , parameter :: w_xyl   = 106.0D0
  real(rk8) , parameter :: w_benz    = 78.110D0

  ! Carbonyls
  real(rk8) , parameter :: w_hcho    = 30.0D0
  real(rk8) , parameter :: w_ald2    = 44.0D0
  real(rk8) , parameter :: w_ket     = 72.0D0
  real(rk8) , parameter :: w_aone    = 72.0D0
  real(rk8) , parameter :: w_gly     = 58.0D0
  real(rk8) , parameter :: w_mgly    = 72.0D0

  ! Organic Nitrate
  real(rk8) , parameter :: w_pan     = 121.0D0
  real(rk8) , parameter :: w_tpan    = 147.0D0
  real(rk8) , parameter :: w_onit    = 119.0D0

  ! Organic Acids
  real(rk8), parameter  :: w_hcooh     = 46.0D0   !Formic acid
  real(rk8), parameter  :: w_ch3cooh   = 60.0D0   !Acetic acid
  real(rk8),  parameter :: w_rcooh   = 59.1D0

  ! Alcohol
  real(rk8), parameter  :: w_moh         = 32.0D0 !Methanol
  real(rk8), parameter  :: w_ch3oh       = 32.0D0 !Methanol
  real(rk8), parameter  :: w_eoh         = 46.0D0 !Ethanol
  real(rk8), parameter  :: w_c2h5oh      = 46.0D0 !Ethanol

  ! Organic Peroxid
  real(rk8), parameter  :: w_ch3ooh   = 48.0D0
  real(rk8), parameter  :: w_ethooh   = 74.0D0
  real(rk8) , parameter :: w_rooh    = 48.0D0

  ! Other species
  real(rk8) , parameter :: w_dms     = 62.0D0
  real(rk8) , parameter :: w_msa     = 96.0D0
  real(rk8) , parameter :: w_nh3     = 17.0D0
  real(rk8) , parameter :: w_apin    = 136.230D0
  real(rk8) , parameter :: w_limo    = 136.230D0

  ! intermediate species that do not undergo other process than chemistryi
  ! are assigned an arbitrary molecular weight , as anyway chemistry works
  ! with mol. and the conversion to mass is done for compatibility with
  ! other processes than chem.
  ! if these species are outputed in mass, the unit will be however wrong !

  real(rk8) , parameter :: w_O1D = 16.D0
  real(rk8) , parameter :: w_cro =48.D0
  real(rk8) , parameter :: w_to2 = 32.D0
  real(rk8) , parameter :: w_dummy = 1.D0
  real(rk8) , parameter :: w_open = 1.D0
  real(rk8) , parameter :: w_O3P = 48.D0
  real(rk8) , parameter :: w_isopn   = 68.0D0
  real(rk8) , parameter :: w_isopp   = 68.0D0
  real(rk8) , parameter :: w_isopo2   = 68.0D0
  real(rk8) , parameter :: w_isoprd   = 68.0D0
  real(rk8) , parameter :: w_ethp   = 28.0D0
  real(rk8) , parameter :: w_nap   = 1.0D0
  real(rk8) , parameter :: w_ch3o2   = 47.0D0
  real(rk8) , parameter :: w_ano2   =  46.0D0
  real(rk8) , parameter :: w_c2o3   =  72.0D0


  ! define molecular weight in CB6C 

  real(rk8) , parameter :: w_cdox =  44.0D0
! real(rk8) , parameter :: w_sulf =  98.0D0 ! defined in cbmz
  real(rk8) , parameter :: w_sdio =  64.0D0
  real(rk8) , parameter :: w_osng =  16.0D0
  real(rk8) , parameter :: w_mtha =  16.0D0
  real(rk8) , parameter :: w_etha =  30.0D0
  real(rk8) , parameter :: w_ethy =  26.0D0
  real(rk8) , parameter :: w_dnpo = 108.0D0
! real(rk8) , parameter :: w_benz =  78.0D0
  real(rk8) , parameter :: w_epox = 134.0D0
  real(rk8) , parameter :: w_etoh =  46.0D0
  real(rk8) , parameter :: w_prpa =  44.0D0
! real(rk8) , parameter :: w_ket  =  72.0D0
  real(rk8) , parameter :: w_toln =  92.0D0
  real(rk8) , parameter :: w_xyln = 106.0D0
  real(rk8) , parameter :: w_hpld = 116.0D0
  real(rk8) , parameter :: w_pacn = 121.0D0
  real(rk8) , parameter :: w_pacd =  76.0D0
  real(rk8) , parameter :: w_ntr2 = 135.0D0
  real(rk8) , parameter :: w_pna  =  79.0D0
  real(rk8) , parameter :: w_meoh =  32.0D0
! real(rk8) , parameter :: w_hono =  47.0D0
  real(rk8) , parameter :: w_mepx =  48.0D0
  real(rk8) , parameter :: w_opan = 161.0D0
  real(rk8) , parameter :: w_cat1 = 124.0D0
  real(rk8) , parameter :: w_hpox =  34.0D0
  real(rk8) , parameter :: w_ispx = 118.0D0
  real(rk8) , parameter :: w_facd =  46.0D0
  real(rk8) , parameter :: w_panx = 135.0D0
  real(rk8) , parameter :: w_hco3 =  63.0D0
  real(rk8) , parameter :: w_crer = 107.0D0
  real(rk8) , parameter :: w_rpox =  90.0D0
  real(rk8) , parameter :: w_ntr1 = 119.0D0 
  real(rk8) , parameter :: w_acet =  58.0D0
  real(rk8) , parameter :: w_intr = 119.0D0
  real(rk8) , parameter :: w_bzo2 = 159.0D0
  real(rk8) , parameter :: w_cron = 153.0D0
  real(rk8) , parameter :: w_aacd =  60.0D0
  real(rk8) , parameter :: w_ror  =  71.0D0
  real(rk8) , parameter :: w_tolr = 173.0D0
  real(rk8) , parameter :: w_ethe =  28.0D0
  real(rk8) , parameter :: w_cmon =  28.0D0
  real(rk8) , parameter :: w_xo2r =  87.0D0
  real(rk8) , parameter :: w_terp = 136.0D0
  real(rk8) , parameter :: w_crsl = 108.0D0
  real(rk8) , parameter :: w_ispr =  68.0D0
  real(rk8) , parameter :: w_epx2 = 165.0D0
  real(rk8) , parameter :: w_ntrc =  63.0D0
  real(rk8) , parameter :: w_glyd =  60.0D0 
! real(rk8) , parameter :: w_gly  =  58.0D0
  real(rk8) , parameter :: w_xopn =  98.0D0
  real(rk8) , parameter :: w_megy =  72.0D0
  real(rk8) , parameter :: w_ropn =  84.0D0
  real(rk8) , parameter :: w_iole =  56.0D0
  real(rk8) , parameter :: w_form =  30.0D0
  real(rk8) , parameter :: w_ole  =  42.0D0
  real(rk8) , parameter :: w_aald =  44.0D0
  real(rk8) , parameter :: w_xylr = 187.0D0
  real(rk8) , parameter :: w_opo3 = 115.0D0
  real(rk8) , parameter :: w_xo2n =  87.0D0
  real(rk8) , parameter :: w_iso2 = 117.0D0
  real(rk8) , parameter :: w_xo2h =  87.0D0
  real(rk8) , parameter :: w_ispd =  70.0D0
  real(rk8) , parameter :: w_meo2 =  47.0D0
  real(rk8) , parameter :: w_aldx =  58.0D0
  real(rk8) , parameter :: w_alka =  72.0D0
  real(rk8) , parameter :: w_ozn  =  48.0D0
  real(rk8) , parameter :: w_roo  =  87.0D0
  real(rk8) , parameter :: w_hox  =  17.0D0
  real(rk8) , parameter :: w_pox  =  33.0D0
  real(rk8) , parameter :: w_cxo3 =  89.0D0
  real(rk8) , parameter :: w_ntox =  62.0D0
  real(rk8) , parameter :: w_acoo =  75.0D0
  real(rk8) , parameter :: w_o    =  16.0D0
  real(rk8) , parameter :: w_nmox =  30.0D0
  real(rk8) , parameter :: w_ndox =  46.0D0
  real(rk8) , parameter :: w_o2   =  32.0D0
  real(rk8) , parameter :: w_m    =   1.0D0
  real(rk8) , parameter :: w_dihy =   2.0D0
  real(rk8) , parameter :: w_wtr  =  18.0D0
  real(rk8) , parameter :: w_meth =  16.0D0
  real(rk8) , parameter :: w_dummy2 = 1.0D0

  real(rk8) , parameter :: w_aggrg =  0.0D0

  ! define here a table of molecular weight for the CBMZ species.

  real(rk8),dimension(totsp) ::  mw_cbmz

  data mw_cbmz / W_CO2 ,   & ! 1
                 W_H2SO4,  & ! 2
                 W_HCOOH,  & ! 3
                 W_RCOOH,  & ! 4
                 W_MSA,    & ! 5
                 W_DUMMY,  & ! 6
                 W_PAN,    & ! 7
                 W_TOL,    & ! 8
                 W_O1D,    & ! 9
                 W_H2O2,   & ! 10
                 W_SO2,    & ! 11
                 W_XYL,    & ! 12
                 W_CH4,    & ! 13
                 W_C2H6,   & ! 14
                 W_CRO,    & ! 15
                 W_DMS,    & ! 16
                 W_HNO4,   & ! 17
                 W_H2,     & ! 18
                 W_TO2,    & ! 19
                 W_CH3OH,  & ! 20
                 W_HNO2,   & ! 21
                 W_CH3OOH, & ! 22
                 W_ETHOOH, & ! 23
                 W_N2O5,   & ! 24
                 W_ETH,    & ! 25
                 W_CRES,   & ! 26
                 W_O3P,    & ! 27
                 W_CO,     & ! 28
                 W_HNO3,   & ! 29
                 W_PAR,    & ! 30
                 W_OPEN,   & ! 31
                 W_ISOPN,  & ! 32
                 W_ISOPP,  & ! 33
                 W_ISOPO2, & ! 34
                 W_H2O,    & ! 35
                 W_AONE,   & ! 36
                 W_OLEI,   & ! 37
                 W_ISOP,   & ! 38
                 W_HCHO,   & ! 39
                 W_OLET,   & ! 40
                 W_XO2,    & ! 41
                 W_MGLY,   & ! 42
                 W_ETHP,   & ! 43
                 W_NAP,    & ! 44
                 W_ALD2,   & ! 45
                 W_CH3O2,  & ! 46
                 W_ISOPRD, & ! 47
                 W_ANO2,   & ! 48
                 W_ROOH,   & ! 49
                 W_RO2,    & ! 50
                 W_ONIT,   & ! 51
                 W_HO2,    & ! 52
                 W_O3,     & ! 53
                 W_OH,     & ! 54
                 W_NO,     & ! 55
                 W_NO2,    & ! 56
                 W_NO3,    & ! 57
                 W_C2O3 /    ! 58

  ! define here a table of molecular weight for the CB6C species.

  real(rk8),dimension(nvar_cb6) ::  mw_cb6

  data mw_cb6  / W_CDOX,   & ! 1
                 W_SULF,   & ! 2
                 W_SDIO,   & ! 3
                 W_OSNG,   & ! 4
                 W_MTHA,   & ! 5
                 W_ETHA,   & ! 6
                 W_ETHY,   & ! 7
                 W_DNPO,   & ! 8
                 W_BENZ,   & ! 9
                 W_EPOX,   & ! 10
                 W_ETOH,   & ! 11
                 W_PRPA,   & ! 12
                 W_KET,    & ! 13
                 W_TOLN,   & ! 14
                 W_XYLN,   & ! 15
                 W_HPLD,   & ! 16
                 W_PACN,   & ! 17
                 W_PACD,   & ! 18
                 W_NTR2,   & ! 19
                 W_PNA,    & ! 20
                 W_MEOH,   & ! 21
                 W_HONO,   & ! 22
                 W_MEPX,   & ! 23
                 W_OPAN,   & ! 24
                 W_CAT1,   & ! 25
                 W_HPOX,   & ! 26
                 W_ISPX,   & ! 27
                 W_FACD,   & ! 28
                 W_PANX,   & ! 29
                 W_HCO3,   & ! 30
                 W_CRER,   & ! 31
                 W_RPOX,   & ! 32
                 W_NTR1,   & ! 33
                 W_ACET,   & ! 34
                 W_INTR,   & ! 35
                 W_BZO2,   & ! 36
                 W_CRON,   & ! 37
                 W_AACD,   & ! 38
                 W_ROR,    & ! 39
                 W_TOLR,   & ! 40
                 W_ETHE,   & ! 41
                 W_CMON,   & ! 42
                 W_XO2R,   & ! 43
                 W_TERP,   & ! 44
                 W_CRSL,   & ! 45
                 W_ISPR,   & ! 46
                 W_EPX2,   & ! 47
                 W_NTRC,   & ! 48
                 W_GLYD,   & ! 49
                 W_GLY,    & ! 50
                 W_XOPN,   & ! 51
                 W_MEGY,   & ! 52
                 W_ROPN,   & ! 53
                 W_IOLE,   & ! 54
                 W_FORM,   & ! 55
                 W_OLE,    & ! 56
                 W_AALD,   & ! 57
                 W_XYLR,   & ! 58
                 W_OPO3,   & ! 59
                 W_XO2N,   & ! 60
                 W_ISO2,   & ! 61
                 W_XO2H,   & ! 62
                 W_ISPD,   & ! 63
                 W_MEO2,   & ! 64
                 W_ALDX,   & ! 65
                 W_ALKA,   & ! 66
                 W_OZN,    & ! 67
                 W_ROO,    & ! 68
                 W_HOX,    & ! 69
                 W_POX,    & ! 70
                 W_CXO3,   & ! 71
                 W_NTOX,   & ! 72
                 W_ACOO,   & ! 73
                 W_O,      & ! 74
                 W_NMOX,   & ! 75
                 W_NDOX/     ! 76

  ! define here a table of molecular weight for the CB6C+SORGAM species.

  real(rk8),dimension(3*nvar_cb6+12) ::  mw_srg

  data mw_srg  / W_CDOX,   & ! 1
                 W_SULF,   & ! 2
                 W_SDIO,   & ! 3
                 W_OSNG,   & ! 4
                 W_MTHA,   & ! 5
                 W_ETHA,   & ! 6
                 W_ETHY,   & ! 7
                 W_DNPO,   & ! 8
                 W_BENZ,   & ! 9
                 W_EPOX,   & ! 10
                 W_ETOH,   & ! 11
                 W_PRPA,   & ! 12
                 W_KET,    & ! 13
                 W_TOLN,   & ! 14
                 W_XYLN,   & ! 15
                 W_HPLD,   & ! 16
                 W_PACN,   & ! 17
                 W_PACD,   & ! 18
                 W_NTR2,   & ! 19
                 W_PNA,    & ! 20
                 W_MEOH,   & ! 21
                 W_HONO,   & ! 22
                 W_MEPX,   & ! 23
                 W_OPAN,   & ! 24
                 W_CAT1,   & ! 25
                 W_HPOX,   & ! 26
                 W_ISPX,   & ! 27
                 W_FACD,   & ! 28
                 W_PANX,   & ! 29
                 W_HCO3,   & ! 30
                 W_CRER,   & ! 31
                 W_RPOX,   & ! 32
                 W_NTR1,   & ! 33
                 W_ACET,   & ! 34
                 W_INTR,   & ! 35
                 W_BZO2,   & ! 36
                 W_CRON,   & ! 37
                 W_AACD,   & ! 38
                 W_ROR,    & ! 39
                 W_TOLR,   & ! 40
                 W_ETHE,   & ! 41
                 W_CMON,   & ! 42
                 W_XO2R,   & ! 43
                 W_TERP,   & ! 44
                 W_CRSL,   & ! 45
                 W_ISPR,   & ! 46
                 W_EPX2,   & ! 47
                 W_NTRC,   & ! 48
                 W_GLYD,   & ! 49
                 W_GLY,    & ! 50
                 W_XOPN,   & ! 51
                 W_MEGY,   & ! 52
                 W_ROPN,   & ! 53
                 W_IOLE,   & ! 54
                 W_FORM,   & ! 55
                 W_OLE,    & ! 56
                 W_AALD,   & ! 57
                 W_XYLR,   & ! 58
                 W_OPO3,   & ! 59
                 W_XO2N,   & ! 60
                 W_ISO2,   & ! 61
                 W_XO2H,   & ! 62
                 W_ISPD,   & ! 63
                 W_MEO2,   & ! 64
                 W_ALDX,   & ! 65
                 W_ALKA,   & ! 66
                 W_OZN,    & ! 67
                 W_ROO,    & ! 68
                 W_HOX,    & ! 69
                 W_POX,    & ! 70
                 W_CXO3,   & ! 71
                 W_NTOX,   & ! 72
                 W_ACOO,   & ! 73
                 W_O,      & ! 74
                 W_NMOX,   & ! 75
                 W_NDOX,   & ! 76
                 W_CDOX,   & ! 77 ! nucleation aerosols
                 W_SULF,   & ! 78
                 W_SDIO,   & ! 79
                 W_OSNG,   & ! 80
                 W_MTHA,   & ! 81
                 W_ETHA,   & ! 82
                 W_ETHY,   & ! 83
                 W_DNPO,   & ! 84
                 W_BENZ,   & ! 85
                 W_EPOX,   & ! 86
                 W_ETOH,   & ! 87
                 W_PRPA,   & ! 88
                 W_KET,    & ! 89
                 W_TOLN,   & ! 90
                 W_XYLN,   & ! 91
                 W_HPLD,   & ! 92
                 W_PACN,   & ! 93
                 W_PACD,   & ! 94
                 W_NTR2,   & ! 95
                 W_PNA,    & ! 96
                 W_MEOH,   & ! 97
                 W_HONO,   & ! 98
                 W_MEPX,   & ! 99
                 W_OPAN,   & ! 100
                 W_CAT1,   & ! 101
                 W_HPOX,   & ! 102
                 W_ISPX,   & ! 103
                 W_FACD,   & ! 104
                 W_PANX,   & ! 105
                 W_HCO3,   & ! 106
                 W_CRER,   & ! 107
                 W_RPOX,   & ! 108
                 W_NTR1,   & ! 109
                 W_ACET,   & ! 110
                 W_INTR,   & ! 111
                 W_BZO2,   & ! 112
                 W_CRON,   & ! 113
                 W_AACD,   & ! 114
                 W_ROR,    & ! 115
                 W_TOLR,   & ! 116
                 W_ETHE,   & ! 117
                 W_CMON,   & ! 118
                 W_XO2R,   & ! 119
                 W_TERP,   & ! 120
                 W_CRSL,   & ! 121
                 W_ISPR,   & ! 122
                 W_EPX2,   & ! 123
                 W_NTRC,   & ! 124
                 W_GLYD,   & ! 125
                 W_GLY,    & ! 126
                 W_XOPN,   & ! 127
                 W_MEGY,   & ! 128
                 W_ROPN,   & ! 129
                 W_IOLE,   & ! 130
                 W_FORM,   & ! 131
                 W_OLE,    & ! 132
                 W_AALD,   & ! 133
                 W_XYLR,   & ! 134
                 W_OPO3,   & ! 135
                 W_XO2N,   & ! 136
                 W_ISO2,   & ! 137
                 W_XO2H,   & ! 138
                 W_ISPD,   & ! 139
                 W_MEO2,   & ! 140
                 W_ALDX,   & ! 141
                 W_ALKA,   & ! 142
                 W_OZN,    & ! 143
                 W_ROO,    & ! 144
                 W_HOX,    & ! 145
                 W_POX,    & ! 146
                 W_CXO3,   & ! 147
                 W_NTOX,   & ! 148
                 W_ACOO,   & ! 149
                 W_O,      & ! 150
                 W_NMOX,   & ! 151
                 W_NDOX,   & ! 152
                 W_CDOX,   & ! 153 ! accumulation aerosols
                 W_SULF,   & ! 154
                 W_SDIO,   & ! 155
                 W_OSNG,   & ! 156
                 W_MTHA,   & ! 157
                 W_ETHA,   & ! 158
                 W_ETHY,   & ! 159
                 W_DNPO,   & ! 160
                 W_BENZ,   & ! 161
                 W_EPOX,   & ! 162
                 W_ETOH,   & ! 163
                 W_PRPA,   & ! 164
                 W_KET,    & ! 165
                 W_TOLN,   & ! 166
                 W_XYLN,   & ! 167
                 W_HPLD,   & ! 168
                 W_PACN,   & ! 169
                 W_PACD,   & ! 170
                 W_NTR2,   & ! 171
                 W_PNA,    & ! 172
                 W_MEOH,   & ! 173
                 W_HONO,   & ! 174
                 W_MEPX,   & ! 175
                 W_OPAN,   & ! 176
                 W_CAT1,   & ! 177
                 W_HPOX,   & ! 178
                 W_ISPX,   & ! 179
                 W_FACD,   & ! 180
                 W_PANX,   & ! 181
                 W_HCO3,   & ! 182
                 W_CRER,   & ! 183
                 W_RPOX,   & ! 184
                 W_NTR1,   & ! 185
                 W_ACET,   & ! 186
                 W_INTR,   & ! 187
                 W_BZO2,   & ! 188
                 W_CRON,   & ! 189
                 W_AACD,   & ! 190
                 W_ROR,    & ! 191
                 W_TOLR,   & ! 192
                 W_ETHE,   & ! 193
                 W_CMON,   & ! 194
                 W_XO2R,   & ! 195
                 W_TERP,   & ! 196
                 W_CRSL,   & ! 197
                 W_ISPR,   & ! 198
                 W_EPX2,   & ! 199
                 W_NTRC,   & ! 200
                 W_GLYD,   & ! 201
                 W_GLY,    & ! 202
                 W_XOPN,   & ! 203
                 W_MEGY,   & ! 204
                 W_ROPN,   & ! 205
                 W_IOLE,   & ! 206
                 W_FORM,   & ! 207
                 W_OLE,    & ! 208
                 W_AALD,   & ! 209
                 W_XYLR,   & ! 210
                 W_OPO3,   & ! 211
                 W_XO2N,   & ! 212
                 W_ISO2,   & ! 213
                 W_XO2H,   & ! 214
                 W_ISPD,   & ! 215
                 W_MEO2,   & ! 216
                 W_ALDX,   & ! 217
                 W_ALKA,   & ! 218
                 W_OZN,    & ! 219
                 W_ROO,    & ! 220
                 W_HOX,    & ! 221
                 W_POX,    & ! 222
                 W_CXO3,   & ! 223
                 W_NTOX,   & ! 224
                 W_ACOO,   & ! 225
                 W_O,      & ! 226
                 W_NMOX,   & ! 227
                 W_NDOX,   & ! 228
                 w_aggrg,  & ! 229 ! aggregate aerosols weight not needed
                 w_aggrg,  & ! 230
                 w_aggrg,  & ! 231
                 w_aggrg,  & ! 232
                 w_aggrg,  & ! 233
                 w_aggrg,  & ! 234
                 w_aggrg,  & ! 235
                 w_aggrg,  & ! 236
                 w_aggrg,  & ! 237
                 w_aggrg,  & ! 238
                 w_aggrg,  & ! 239
                 w_aggrg/    ! 240

end module mod_che_molwg

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
