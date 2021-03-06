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

module mod_ch_param

  use mod_intkinds
  use mod_realkinds

  public




!Mozart indcies
 INTEGER(ik4), PARAMETER :: mz_NO      = 1
 INTEGER(ik4), PARAMETER :: mz_NO2     = 2
 INTEGER(ik4), PARAMETER :: mz_N2O5    = 3
 INTEGER(ik4), PARAMETER :: mz_HNO3    = 4
 INTEGER(ik4), PARAMETER :: mz_HO2NO2  = 5
 INTEGER(ik4), PARAMETER :: mz_O3      = 6
 INTEGER(ik4), PARAMETER :: mz_H2O2    = 7
 INTEGER(ik4), PARAMETER :: mz_SO2     = 8
 INTEGER(ik4), PARAMETER :: mz_SO4     = 9
 INTEGER(ik4), PARAMETER :: mz_CH4     = 10
 INTEGER(ik4), PARAMETER :: mz_CH2O    = 11
 INTEGER(ik4), PARAMETER :: mz_CH3OH   = 12
 INTEGER(ik4), PARAMETER :: mz_PAN     = 13
 INTEGER(ik4), PARAMETER :: mz_C2H6    = 14
 INTEGER(ik4), PARAMETER :: mz_C3H8    = 15
 INTEGER(ik4), PARAMETER :: mz_BIGALK  = 16
 INTEGER(ik4), PARAMETER :: mz_C2H4    = 17
 INTEGER(ik4), PARAMETER :: mz_C3H6    = 18
 INTEGER(ik4), PARAMETER :: mz_BIGENE  = 19
 INTEGER(ik4), PARAMETER :: mz_TOLUENE = 20
 INTEGER(ik4), PARAMETER :: mz_ISOP    = 21
 INTEGER(ik4), PARAMETER :: mz_CH3CHO  = 22
 INTEGER(ik4), PARAMETER :: mz_CH3COOH = 23
 INTEGER(ik4), PARAMETER :: mz_GLYALD  = 24
 INTEGER(ik4), PARAMETER :: mz_CH3OOH  = 25
 INTEGER(ik4), PARAMETER :: mz_C2H5OOH = 26
 INTEGER(ik4), PARAMETER :: mz_CH3COCH3= 27
 INTEGER(ik4), PARAMETER :: mz_HYAC    = 28
 INTEGER(ik4), PARAMETER :: mz_CH3COCHO =29
 INTEGER(ik4), PARAMETER :: mz_ONIT    = 30
 INTEGER(ik4), PARAMETER :: mz_MEK     = 31
 INTEGER(ik4), PARAMETER :: mz_MVK     = 32
 INTEGER(ik4), PARAMETER :: mz_MACR    = 33
 INTEGER(ik4), PARAMETER :: mz_HYDRALD = 34
 INTEGER(ik4), PARAMETER :: mz_BIGALD  = 35
 INTEGER(ik4), PARAMETER :: mz_ISOPNO3 = 36
 INTEGER(ik4), PARAMETER :: mz_ONITR   = 37
 INTEGER(ik4), PARAMETER :: mz_CRESOL  = 38
 INTEGER(ik4), PARAMETER :: mz_CO      = 39
 INTEGER(ik4), PARAMETER :: mz_DMS     = 40
 INTEGER(ik4), PARAMETER :: mz_O       = 41 ! These 18 species were added for CB6C
 INTEGER(ik4), PARAMETER :: mz_O1D     = 42
 INTEGER(ik4), PARAMETER :: mz_NO3     = 43
 INTEGER(ik4), PARAMETER :: mz_H2      = 44
 INTEGER(ik4), PARAMETER :: mz_OH      = 45
 INTEGER(ik4), PARAMETER :: mz_HO2     = 46
 INTEGER(ik4), PARAMETER :: mz_CH3O2   = 47
 INTEGER(ik4), PARAMETER :: mz_C2H5OH  = 48
 INTEGER(ik4), PARAMETER :: mz_CH3CO3  = 49
 INTEGER(ik4), PARAMETER :: mz_GLYOXAL = 50
 INTEGER(ik4), PARAMETER :: mz_ENEO2   = 51
 INTEGER(ik4), PARAMETER :: mz_ISOPO2  = 52
 INTEGER(ik4), PARAMETER :: mz_ISOPOOH = 53
 INTEGER(ik4), PARAMETER :: mz_XOOH    = 54
 INTEGER(ik4), PARAMETER :: mz_C10H16  = 55
 INTEGER(ik4), PARAMETER :: mz_TOLO2   = 56

!CBMZ indcies
 INTEGER(ik4), PARAMETER :: cb_O3       = 1
 INTEGER(ik4), PARAMETER :: cb_NO       = 2
 INTEGER(ik4), PARAMETER :: cb_NO2      = 3
 INTEGER(ik4), PARAMETER :: cb_HNO3     = 4
 INTEGER(ik4), PARAMETER :: cb_HNO4     = 5
 INTEGER(ik4), PARAMETER :: cb_N2O5     = 6
 INTEGER(ik4), PARAMETER :: cb_H2O2     = 7
 INTEGER(ik4), PARAMETER :: cb_CH4      = 8
 INTEGER(ik4), PARAMETER :: cb_CO       = 9
 INTEGER(ik4), PARAMETER :: cb_SO2      = 10
 INTEGER(ik4), PARAMETER :: cb_H2SO4    = 11
 INTEGER(ik4), PARAMETER :: cb_DMS      = 12
 INTEGER(ik4), PARAMETER :: cb_PAR      = 13
 INTEGER(ik4), PARAMETER :: cb_C2H6     = 14
 INTEGER(ik4), PARAMETER :: cb_ETH      = 15
 INTEGER(ik4), PARAMETER :: cb_OLET     = 16
 INTEGER(ik4), PARAMETER :: cb_OLEI     = 17
 INTEGER(ik4), PARAMETER :: cb_TOL      = 18
 INTEGER(ik4), PARAMETER :: cb_XYL      = 19
 INTEGER(ik4), PARAMETER :: cb_ISOP     = 20
 INTEGER(ik4), PARAMETER :: cb_CRES     = 21
 INTEGER(ik4), PARAMETER :: cb_OPEN     = 22
 INTEGER(ik4), PARAMETER :: cb_ISOPN    = 23
 INTEGER(ik4), PARAMETER :: cb_ISOPRD   = 24
 INTEGER(ik4), PARAMETER :: cb_ONIT     = 25
 INTEGER(ik4), PARAMETER :: cb_MGLY     = 26
 INTEGER(ik4), PARAMETER :: cb_AONE     = 27
 INTEGER(ik4), PARAMETER :: cb_PAN      = 28
 INTEGER(ik4), PARAMETER :: cb_CH3OOH   = 29
 INTEGER(ik4), PARAMETER :: cb_ETHOOH   = 30
 INTEGER(ik4), PARAMETER :: cb_ALD2     = 31
 INTEGER(ik4), PARAMETER :: cb_HCHO     = 32
 INTEGER(ik4), PARAMETER :: cb_CH3OH    = 33
!CB6C indcies
 INTEGER(ik4), PARAMETER :: cb6_NMOX    = 1
 INTEGER(ik4), PARAMETER :: cb6_NDOX    = 2
 INTEGER(ik4), PARAMETER :: cb6_DNPO    = 3
 INTEGER(ik4), PARAMETER :: cb6_NTRC    = 4
 INTEGER(ik4), PARAMETER :: cb6_PNA     = 5
 INTEGER(ik4), PARAMETER :: cb6_OZN     = 6
 INTEGER(ik4), PARAMETER :: cb6_HPOX    = 7
 INTEGER(ik4), PARAMETER :: cb6_SDIO    = 8
 INTEGER(ik4), PARAMETER :: cb6_SULF    = 9
 INTEGER(ik4), PARAMETER :: cb6_MTHA    = 10
 INTEGER(ik4), PARAMETER :: cb6_FORM    = 11
 INTEGER(ik4), PARAMETER :: cb6_MEOH    = 12
 INTEGER(ik4), PARAMETER :: cb6_PACN    = 13
 INTEGER(ik4), PARAMETER :: cb6_ETHA    = 14
 INTEGER(ik4), PARAMETER :: cb6_PRPA    = 15
 INTEGER(ik4), PARAMETER :: cb6_ALKA    = 16
 INTEGER(ik4), PARAMETER :: cb6_IOLE    = 17
 INTEGER(ik4), PARAMETER :: cb6_OLE     = 18
 INTEGER(ik4), PARAMETER :: cb6_ETHE    = 19
 INTEGER(ik4), PARAMETER :: cb6_TOLN    = 20
 INTEGER(ik4), PARAMETER :: cb6_XYLN    = 21
 INTEGER(ik4), PARAMETER :: cb6_ISPR    = 22
 INTEGER(ik4), PARAMETER :: cb6_AALD    = 23
 INTEGER(ik4), PARAMETER :: cb6_AACD    = 24
 INTEGER(ik4), PARAMETER :: cb6_GLYD    = 25
 INTEGER(ik4), PARAMETER :: cb6_MEPX    = 26
 INTEGER(ik4), PARAMETER :: cb6_ACET    = 27
 INTEGER(ik4), PARAMETER :: cb6_MEGY    = 28
 INTEGER(ik4), PARAMETER :: cb6_NTR2    = 29
 INTEGER(ik4), PARAMETER :: cb6_INTR    = 30
 INTEGER(ik4), PARAMETER :: cb6_KET     = 31
 INTEGER(ik4), PARAMETER :: cb6_ISPD    = 32
 INTEGER(ik4), PARAMETER :: cb6_HPLD    = 33
 INTEGER(ik4), PARAMETER :: cb6_XOPN    = 34
 INTEGER(ik4), PARAMETER :: cb6_NTR1    = 35
 INTEGER(ik4), PARAMETER :: cb6_CRSL    = 36
 INTEGER(ik4), PARAMETER :: cb6_CMON    = 37
 INTEGER(ik4), PARAMETER :: cb6_CRON    = 38
 INTEGER(ik4), PARAMETER :: cb6_OSNG    = 39
 INTEGER(ik4), PARAMETER :: cb6_NTOX    = 40
 INTEGER(ik4), PARAMETER :: cb6_HOX     = 41
 INTEGER(ik4), PARAMETER :: cb6_POX     = 42
 INTEGER(ik4), PARAMETER :: cb6_MEO2    = 43
 INTEGER(ik4), PARAMETER :: cb6_ETOH    = 44
 INTEGER(ik4), PARAMETER :: cb6_ACOO    = 45
 INTEGER(ik4), PARAMETER :: cb6_ROO     = 46
 INTEGER(ik4), PARAMETER :: cb6_ISO2    = 47
 INTEGER(ik4), PARAMETER :: cb6_ISPX    = 48
 INTEGER(ik4), PARAMETER :: cb6_TERP    = 49
 INTEGER(ik4), PARAMETER :: cb6_TOLR    = 50
 INTEGER(ik4), PARAMETER :: cb6_XYLR    = 51
 INTEGER(ik4), PARAMETER :: cb6_GLY     = 52
 INTEGER(ik4), PARAMETER :: cb6_ROPN    = 53
 INTEGER(ik4), PARAMETER :: cb6_BENZ    = 54

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
  real(rk8) , parameter :: w_eth = 28.0D0
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

  ! CB6C species
  real(rk8) , parameter :: w_nmox =  30.0D0
  real(rk8) , parameter :: w_ndox =  46.0D0
  real(rk8) , parameter :: w_dnpo = 108.0D0
  real(rk8) , parameter :: w_ntrc =  63.0D0
  real(rk8) , parameter :: w_pna  =  79.0D0
  real(rk8) , parameter :: w_ozn  =  48.0D0
  real(rk8) , parameter :: w_hpox =  34.0D0
  real(rk8) , parameter :: w_sdio =  64.0D0
  real(rk8) , parameter :: w_mtha =  16.0D0
  real(rk8) , parameter :: w_form =  30.0D0
  real(rk8) , parameter :: w_meoh =  32.0D0
  real(rk8) , parameter :: w_pacn = 121.0D0
  real(rk8) , parameter :: w_etha =  30.0D0
  real(rk8) , parameter :: w_prpa =  44.0D0
  real(rk8) , parameter :: w_alka =  72.0D0
  real(rk8) , parameter :: w_iole =  56.0D0
  real(rk8) , parameter :: w_ole  =  42.0D0
  real(rk8) , parameter :: w_ethe =  28.0D0
  real(rk8) , parameter :: w_toln =  92.0D0
  real(rk8) , parameter :: w_xyln = 106.0D0
  real(rk8) , parameter :: w_ispr =  68.0D0
  real(rk8) , parameter :: w_aald =  44.0D0
  real(rk8) , parameter :: w_aacd =  60.0D0
  real(rk8) , parameter :: w_glyd =  76.0D0
  real(rk8) , parameter :: w_mepx =  48.0D0
  real(rk8) , parameter :: w_acet =  58.0D0
  real(rk8) , parameter :: w_megy =  72.0D0
  real(rk8) , parameter :: w_ntr  = 118.0D0
  real(rk8) , parameter :: w_ntr2 = 135.0D0
  real(rk8) , parameter :: w_intr = 119.0D0
  real(rk8) , parameter :: w_ispd =  70.0D0
  real(rk8) , parameter :: w_hpld = 116.0D0
  real(rk8) , parameter :: w_xopn = 112.0D0
  real(rk8) , parameter :: w_ntr1 = 147.0D0
  real(rk8) , parameter :: w_crsl = 108.0D0
  real(rk8) , parameter :: w_cmon =  28.0D0
  real(rk8) , parameter :: w_o    =  16.0D0
  real(rk8) , parameter :: w_osng =  16.0D0
  real(rk8) , parameter :: w_ntox =  62.0D0
  real(rk8) , parameter :: w_dihy =   2.0D0
  real(rk8) , parameter :: w_hox  =  17.0D0
  real(rk8) , parameter :: w_pox  =  33.0D0
  real(rk8) , parameter :: w_meo2 =  47.0D0
  real(rk8) , parameter :: w_etoh =  46.0D0
  real(rk8) , parameter :: w_acoo =  75.0D0
  real(rk8) , parameter :: w_roo  =  87.0D0
  real(rk8) , parameter :: w_iso2 = 117.0D0
  real(rk8) , parameter :: w_ispx = 118.0D0
  real(rk8) , parameter :: w_terp = 136.0D0
  real(rk8) , parameter :: w_tolr = 173.0D0
  real(rk8) , parameter :: w_xylr = 187.0D0
  real(rk8) , parameter :: w_ropn =  84.0D0
  real(rk8) , parameter :: w_cron = 153.0D0 


end module mod_ch_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
