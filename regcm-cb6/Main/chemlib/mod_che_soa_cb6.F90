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

module mod_che_soa_cb6

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_indices

  implicit none

  public

  ! category parameters for SOAs from CB6C species
  ! soa_cat = SOA Categories 
  ! 0 = no aerosol
  ! 2 = carbonyls / dicarbonyls
  ! 3 = alkanes / alcohols
  ! 4 = olefins
  ! 6 = aromatics
  ! 9 = organic nitrates
  ! 10 = oraganic multifunctionals
  ! 21 = Nitric Acid
  ! 22 = Sulfuric Acid
  ! 99 = ROG - Reactive Organic Gas - Precursor

  integer(ik4), parameter :: soa_cat_NTR  = 9
  integer(ik4), parameter :: soa_cat_SULF = 22
  integer(ik4), parameter :: soa_cat_SDIO = 0
  integer(ik4), parameter :: soa_cat_OSNG = 0
  integer(ik4), parameter :: soa_cat_ECH4 = 0
  integer(ik4), parameter :: soa_cat_ETHA = 99
  integer(ik4), parameter :: soa_cat_ETHY = 99
  integer(ik4), parameter :: soa_cat_DNPO = 0
  integer(ik4), parameter :: soa_cat_BENZ = 99
  integer(ik4), parameter :: soa_cat_EPOX = 10
  integer(ik4), parameter :: soa_cat_ETOH = 99
  integer(ik4), parameter :: soa_cat_PRPA = 99
  integer(ik4), parameter :: soa_cat_KET  = 2
  integer(ik4), parameter :: soa_cat_TOLN = 99
  integer(ik4), parameter :: soa_cat_XYLN = 99
  integer(ik4), parameter :: soa_cat_HPLD = 0
  integer(ik4), parameter :: soa_cat_PACN = 0
  integer(ik4), parameter :: soa_cat_PACD = 0
  integer(ik4), parameter :: soa_cat_NTR2 = 9
  integer(ik4), parameter :: soa_cat_PNA  = 0
  integer(ik4), parameter :: soa_cat_MEOH = 3
  integer(ik4), parameter :: soa_cat_HONO = 0
  integer(ik4), parameter :: soa_cat_MEPX = 0
  integer(ik4), parameter :: soa_cat_OPAN = 0
  integer(ik4), parameter :: soa_cat_CAT1 = 6
  integer(ik4), parameter :: soa_cat_HPOX = 0
  integer(ik4), parameter :: soa_cat_ISPX = 0
  integer(ik4), parameter :: soa_cat_FACD = 2
  integer(ik4), parameter :: soa_cat_PANX = 0
  integer(ik4), parameter :: soa_cat_HCO3 = 0
  integer(ik4), parameter :: soa_cat_CRER = 0
  integer(ik4), parameter :: soa_cat_RPOX = 0
  integer(ik4), parameter :: soa_cat_NTR1 = 9
  integer(ik4), parameter :: soa_cat_ACET = 2
  integer(ik4), parameter :: soa_cat_INTR = 9
  integer(ik4), parameter :: soa_cat_BZO2 = 0
  integer(ik4), parameter :: soa_cat_CRON = 6
  integer(ik4), parameter :: soa_cat_AACD = 2
  integer(ik4), parameter :: soa_cat_ROR  = 0
  integer(ik4), parameter :: soa_cat_TOLR = 0
  integer(ik4), parameter :: soa_cat_ETHE = 99
  integer(ik4), parameter :: soa_cat_CMON = 0
  integer(ik4), parameter :: soa_cat_XLO2 = 0
  integer(ik4), parameter :: soa_cat_TERP = 99
  integer(ik4), parameter :: soa_cat_CRSL = 6
  integer(ik4), parameter :: soa_cat_ISPR = 99
  integer(ik4), parameter :: soa_cat_EPX2 = 0
  integer(ik4), parameter :: soa_cat_NTRC = 0 
  integer(ik4), parameter :: soa_cat_GLYD = 10
  integer(ik4), parameter :: soa_cat_GLY  = 10
  integer(ik4), parameter :: soa_cat_XOPN = 10
  integer(ik4), parameter :: soa_cat_MEGY = 10
  integer(ik4), parameter :: soa_cat_ROPN = 10
  integer(ik4), parameter :: soa_cat_IOLE = 4
  integer(ik4), parameter :: soa_cat_FORM = 2
  integer(ik4), parameter :: soa_cat_OLE  = 4
  integer(ik4), parameter :: soa_cat_AALD = 2
  integer(ik4), parameter :: soa_cat_XYLR = 0
  integer(ik4), parameter :: soa_cat_OPO3 = 0
  integer(ik4), parameter :: soa_cat_XO2N = 0
  integer(ik4), parameter :: soa_cat_ISO2 = 0
  integer(ik4), parameter :: soa_cat_XO2H = 0
  integer(ik4), parameter :: soa_cat_ISPD = 10
  integer(ik4), parameter :: soa_cat_MEO2 = 0
  integer(ik4), parameter :: soa_cat_ALDX = 2
  integer(ik4), parameter :: soa_cat_ALKA = 3
  integer(ik4), parameter :: soa_cat_OZN  = 0
  integer(ik4), parameter :: soa_cat_ROO  = 0
  integer(ik4), parameter :: soa_cat_HOX  = 0
  integer(ik4), parameter :: soa_cat_POX  = 0
  integer(ik4), parameter :: soa_cat_CXO3 = 0
  integer(ik4), parameter :: soa_cat_NTOX = 0
  integer(ik4), parameter :: soa_cat_ACOO = 0
  integer(ik4), parameter :: soa_cat_O    = 0
  integer(ik4), parameter :: soa_cat_NMOX = 0
  integer(ik4), parameter :: soa_cat_NDOX = 0

  ! define here a table of SOA selectors for the CB6C species.

  integer(ik4), dimension(nvar_cb6) ::  soa_cat_cb6

  data soa_cat_cb6 / soa_cat_NTR,    & ! 1
                     soa_cat_SULF,   & ! 2
                     soa_cat_SDIO,   & ! 3
                     soa_cat_OSNG,   & ! 4
                     soa_cat_ECH4,   & ! 5
                     soa_cat_ETHA,   & ! 6
                     soa_cat_ETHY,   & ! 7
                     soa_cat_DNPO,   & ! 8
                     soa_cat_BENZ,   & ! 9
                     soa_cat_EPOX,   & ! 10
                     soa_cat_ETOH,   & ! 11
                     soa_cat_PRPA,   & ! 12
                     soa_cat_KET,    & ! 13
                     soa_cat_TOLN,   & ! 14
                     soa_cat_XYLN,   & ! 15
                     soa_cat_HPLD,   & ! 16
                     soa_cat_PACN,   & ! 17
                     soa_cat_PACD,   & ! 18
                     soa_cat_NTR2,   & ! 19
                     soa_cat_PNA,    & ! 20
                     soa_cat_MEOH,   & ! 21
                     soa_cat_HONO,   & ! 22
                     soa_cat_MEPX,   & ! 23
                     soa_cat_OPAN,   & ! 24
                     soa_cat_CAT1,   & ! 25
                     soa_cat_HPOX,   & ! 26
                     soa_cat_ISPX,   & ! 27
                     soa_cat_FACD,   & ! 28
                     soa_cat_PANX,   & ! 29
                     soa_cat_HCO3,   & ! 30
                     soa_cat_CRER,   & ! 31
                     soa_cat_RPOX,   & ! 32
                     soa_cat_NTR1,   & ! 33
                     soa_cat_ACET,   & ! 34
                     soa_cat_INTR,   & ! 35
                     soa_cat_BZO2,   & ! 36
                     soa_cat_CRON,   & ! 37
                     soa_cat_AACD,   & ! 38
                     soa_cat_ROR,    & ! 39
                     soa_cat_TOLR,   & ! 40
                     soa_cat_ETHE,   & ! 41
                     soa_cat_CMON,   & ! 42
                     soa_cat_XLO2,   & ! 43
                     soa_cat_TERP,   & ! 44
                     soa_cat_CRSL,   & ! 45
                     soa_cat_ISPR,   & ! 46
                     soa_cat_EPX2,   & ! 47
                     soa_cat_NTRC,   & ! 48
                     soa_cat_GLYD,   & ! 49
                     soa_cat_GLY,    & ! 50
                     soa_cat_XOPN,   & ! 51
                     soa_cat_MEGY,   & ! 52
                     soa_cat_ROPN,   & ! 53
                     soa_cat_IOLE,   & ! 54
                     soa_cat_FORM,   & ! 55
                     soa_cat_OLE,    & ! 56
                     soa_cat_AALD,   & ! 57
                     soa_cat_XYLR,   & ! 58
                     soa_cat_OPO3,   & ! 59
                     soa_cat_XO2N,   & ! 60
                     soa_cat_ISO2,   & ! 61
                     soa_cat_XO2H,   & ! 62
                     soa_cat_ISPD,   & ! 63
                     soa_cat_MEO2,   & ! 64
                     soa_cat_ALDX,   & ! 65
                     soa_cat_ALKA,   & ! 66
                     soa_cat_OZN,    & ! 67
                     soa_cat_ROO,    & ! 68
                     soa_cat_HOX,    & ! 69
                     soa_cat_POX,    & ! 70
                     soa_cat_CXO3,   & ! 71
                     soa_cat_NTOX,   & ! 72
                     soa_cat_ACOO,   & ! 73
                     soa_cat_O,      & ! 74
                     soa_cat_NMOX,   & ! 75
                     soa_cat_NDOX/     ! 76

  integer(ik4) , parameter :: nsoa  = 12
! integer(ik4) , parameter :: prnuc = 1    ! first Aitken index
! integer(ik4) , parameter :: ltnuc = 6    ! last Aitken index
! integer(ik4) , parameter :: pracc = 7    ! first accumulation index
! integer(ik4) , parameter :: ltacc = nsoa ! last accum index
! integer(ik4) , parameter :: nacc  = nsoa-ltnuc ! number of accum SOAs
! integer(ik4) , parameter :: nnuc  = nsoa-nacc  ! number of Aitke SOAs

  !!! NOTE all these optical properties are dummies, need to be replaced!!!!
  ! Define all k values and organise appropriately in matrix
  real(rk8) , dimension(nspi) , parameter :: k_nu_CNYL = &
                      (/2.35534D-01, 1.19018D-01, 9.14269D-02, 7.77163D-02, &
                        6.65251D-02, 5.73101D-02, 3.98556D-02, 9.61258D-03, &
                        2.00215D-03, 2.71931D-02, 2.71923D-02, 2.71923D-02, &
                        2.71923D-02, 2.71923D-02, 2.71914D-02, 2.71914D-02, &
                        3.17975D-02, 3.27297D-03, 3.27297D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_nu_STHC = &
                      (/1.69419D-01, 9.45272D-02, 7.46227D-02, 6.43585D-02, &
                        5.57891D-02, 4.85933D-02, 3.45639D-02, 8.69899D-03, &
                        1.84669D-03, 9.86224D-03, 9.81814D-03, 9.81814D-03, &
                        9.81814D-03, 9.81814D-03, 9.77402D-03, 9.77402D-03, &
                        7.81669D-04, 1.53541D-03, 1.53541D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_nu_OLEF = &
                      (/2.17141D-01, 1.09724D-01, 8.42875D-02, 7.16476D-02, &
                        6.13303D-02, 5.28348D-02, 3.67433D-02, 8.86194D-03, &
                        1.84580D-03, 2.50697D-02, 2.50689D-02, 2.50689D-02, &
                        2.50689D-02, 2.50689D-02, 2.50681D-02, 2.50681D-02, &
                        2.93145D-02, 3.01739D-03, 3.01739D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_nu_AROM = &
                      (/1.35181D-01, 7.51237D-02, 5.92252D-02, 5.10390D-02, &
                        4.42111D-02, 3.84912D-02, 2.74580D-02, 6.97449D-03, &
                        1.46065D-03, 1.15217D-03, 1.15065D-03, 1.15065D-03, &
                        1.15065D-03, 1.15065D-03, 1.14913D-03, 1.14913D-03, &
                        2.78195D-04, 4.41930D-04, 4.41930D-04/)
  real(rk8) , dimension(nspi) , parameter :: k_nu_NTRO = &
                      (/1.60196D-01, 8.09486D-02, 6.21829D-02, 5.28578D-02, &
                        4.52462D-02, 3.89787D-02, 2.71073D-02, 6.53787D-03, &
                        1.36174D-03, 1.84951D-02, 1.84945D-02, 1.84945D-02, &
                        1.84945D-02, 1.84945D-02, 1.84939D-02, 1.84939D-02, &
                        2.16267D-02, 2.22607D-03, 2.22607D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_nu_MFNC = &
                      (/1.47279D-01, 7.44215D-02, 5.71689D-02, 4.85957D-02, &
                        4.15979D-02, 3.58358D-02, 2.49216D-02, 6.01071D-03, &
                        1.25194D-03, 1.70038D-02, 1.70032D-02, 1.70032D-02, &
                        1.70032D-02, 1.70032D-02, 1.70027D-02, 1.70027D-02, &
                        1.98829D-02, 2.04658D-03, 2.04658D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_CNYL = &
                      (/8.21925D+00, 1.22037D+01, 1.30340D+01, 1.44094D+01, &
                        1.57298D+01, 1.57553D+01, 1.62409D+01, 1.15205D+01, &
                        6.54074D+00, 5.60093D-01, 5.58608D-01, 5.58608D-01, &
                        5.58608D-01, 5.58608D-01, 5.57182D-01, 5.57182D-01, &
                        8.89227D-02, 1.36137D-02, 1.36137D-02/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_STHC = &
                      (/1.74331D+01, 1.93958D+01, 2.02185D+01, 2.01725D+01, &
                        1.96052D+01, 1.93418D+01, 1.86515D+01, 1.12016D+01, &
                        5.88181D+00, 4.72530D-01, 4.71177D-01, 4.71177D-01, &
                        4.71177D-01, 4.71177D-01, 4.69828D-01, 4.69828D-01, &
                        5.39505D-02, 1.03960D-02, 1.03960D-02/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_OLEF = &
                      (/7.57742D+00, 1.12507D+01, 1.20162D+01, 1.32842D+01, &
                        1.45015D+01, 1.45250D+01, 1.49727D+01, 1.06209D+01, &
                        6.02998D+00, 5.16356D-01, 5.14987D-01, 5.14987D-01, &
                        5.14987D-01, 5.14987D-01, 5.13672D-01, 5.13672D-01, &
                        8.19788D-02, 1.25507D-02, 1.25507D-02/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_AROM = &
                      (/4.59196D+00, 6.34892D+00, 6.82500D+00, 7.27739D+00, &
                        8.19759D+00, 8.81775D+00, 9.10538D+00, 7.82871D+00, &
                        4.74797D+00, 4.02416D-01, 4.01394D-01, 4.01394D-01, &
                        4.01394D-01, 4.01394D-01, 4.00374D-01, 4.00374D-01, &
                        4.64314D-02, 7.72708D-03, 7.72708D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_NTRO = &
                      (/5.59022D+00, 8.30017D+00, 8.86488D+00, 9.80036D+00, &
                        1.06984D+01, 1.07158D+01, 1.10460D+01, 7.83552D+00, &
                        4.44860D+00, 3.80940D-01, 3.79930D-01, 3.79930D-01, &
                        3.79930D-01, 3.79930D-01, 3.78960D-01, 3.78960D-01, &
                        6.04796D-02, 9.25921D-03, 9.25921D-03/)
  real(rk8) , dimension(nspi) , parameter :: k_ac_MFNC = &
                      (/5.13947D+00, 7.63091D+00, 8.15008D+00, 9.01013D+00, &
                        9.83579D+00, 9.85175D+00, 1.01554D+01, 7.20372D+00, &
                        4.08990D+00, 3.50224D-01, 3.49296D-01, 3.49296D-01, &
                        3.49296D-01, 3.49296D-01, 3.48404D-01, 3.48404D-01, &
                        5.56030D-02, 8.51262D-03, 8.51262D-03/)
  real(rk8) , dimension(nsoa, nspi) , parameter :: ka_soa = &
         reshape((/ k_nu_CNYL , k_nu_STHC , k_nu_OLEF , k_nu_AROM ,   &
                    k_nu_NTRO , k_nu_MFNC , k_ac_CNYL , k_ac_STHC ,   &
                    k_ac_OLEF , k_ac_AROM , k_ac_NTRO , k_ac_MFNC /)  &
                 ,(/nsoa, nspi/), order =(/2,1/) )

  ! Define all w values and organise appropriately in matrix
  real(rk8) , dimension(nspi) , parameter :: w_nu_CNYL = &
                      (/9.97599D-01, 9.99593D-01, 9.99567D-01, 9.99546D-01, &
                        9.99565D-01, 9.99553D-01, 9.99861D-01, 9.99804D-01, &
                        9.99428D-01, 1.54770D-01, 1.54620D-01, 1.54620D-01, &
                        1.54620D-01, 1.54620D-01, 1.54469D-01, 1.54469D-01, &
                        8.57653D-04, 4.92071D-04, 4.92071D-04/)
  real(rk8) , dimension(nspi) , parameter :: w_nu_STHC = &
                      (/9.99998D-01, 9.99997D-01, 9.99997D-01, 9.99996D-01, &
                        9.99996D-01, 9.99996D-01, 9.99994D-01, 9.99966D-01, &
                        9.99815D-01, 1.61179D-01, 1.60977D-01, 1.60977D-01, &
                        1.60977D-01, 1.60977D-01, 1.60784D-01, 1.60784D-01, &
                        1.01407D-02, 6.72980D-04, 6.72980D-04/)
  real(rk8) , dimension(nspi) , parameter :: w_nu_OLEF = &
                      (/9.97599D-01, 9.99593D-01, 9.99567D-01, 9.99546D-01, &
                        9.99565D-01, 9.99553D-01, 9.99861D-01, 9.99804D-01, &
                        9.99428D-01, 1.54770D-01, 1.54620D-01, 1.54620D-01, &
                        1.54620D-01, 1.54620D-01, 1.54469D-01, 1.54469D-01, &
                        8.57653D-04, 4.92071D-04, 4.92071D-04/)
  real(rk8) , dimension(nspi) , parameter :: w_nu_AROM = &
                      (/9.99997D-01, 9.99995D-01, 9.99994D-01, 9.99993D-01, &
                        9.99992D-01, 9.99992D-01, 9.99989D-01, 9.99966D-01, &
                        9.99920D-01, 2.55058D-01, 2.54922D-01, 2.54922D-01, &
                        2.54922D-01, 2.54922D-01, 2.54781D-01, 2.54781D-01, &
                        1.93985D-02, 1.94241D-03, 1.94241D-03/)
  real(rk8) , dimension(nspi) , parameter :: w_nu_NTRO = &
                      (/9.97599D-01, 9.99593D-01, 9.99567D-01, 9.99546D-01, &
                        9.99565D-01, 9.99553D-01, 9.99861D-01, 9.99804D-01, &
                        9.99428D-01, 1.54770D-01, 1.54620D-01, 1.54620D-01, &
                        1.54620D-01, 1.54620D-01, 1.54469D-01, 1.54469D-01, &
                        8.57653D-04, 4.92071D-04, 4.92071D-04/)
  real(rk8) , dimension(nspi) , parameter :: w_nu_MFNC = &
                      (/9.97599D-01, 9.99593D-01, 9.99567D-01, 9.99546D-01, &
                        9.99565D-01, 9.99553D-01, 9.99861D-01, 9.99804D-01, &
                        9.99428D-01, 1.54770D-01, 1.54620D-01, 1.54620D-01, &
                        1.54620D-01, 1.54620D-01, 1.54469D-01, 1.54469D-01, &
                        8.57653D-04, 4.92071D-04, 4.92071D-04/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_CNYL = &
                      (/9.99720D-01, 9.99989D-01, 9.99993D-01, 9.99994D-01, &
                        9.99995D-01, 9.99996D-01, 9.99999D-01, 1.00000D+00, &
                        1.00000D+00, 7.39413D-01, 7.39381D-01, 7.39381D-01, &
                        7.39381D-01, 7.39381D-01, 7.39349D-01, 7.39349D-01, &
                        7.03119D-01, 7.50131D-01, 7.50131D-01/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_STHC = &
                      (/1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, &
                        1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, &
                        1.00000D+00, 8.82593D-01, 8.82558D-01, 8.82558D-01, &
                        8.82558D-01, 8.82558D-01, 8.82523D-01, 8.82523D-01, &
                        9.84732D-01, 8.48367D-01, 8.48367D-01/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_OLEF = &
                      (/9.99720D-01, 9.99989D-01, 9.99993D-01, 9.99994D-01, &
                        9.99995D-01, 9.99996D-01, 9.99999D-01, 1.00000D+00, &
                        1.00000D+00, 7.39413D-01, 7.39381D-01, 7.39381D-01, &
                        7.39381D-01, 7.39381D-01, 7.39349D-01, 7.39349D-01, &
                        7.03119D-01, 7.50131D-01, 7.50131D-01/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_AROM = &
                      (/1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, &
                        1.00000D+00, 1.00000D+00, 1.00000D+00, 1.00000D+00, &
                        1.00000D+00, 9.53201D-01, 9.53228D-01, 9.53228D-01, &
                        9.53228D-01, 9.53228D-01, 9.53255D-01, 9.53255D-01, &
                        9.93281D-01, 9.40317D-01, 9.40317D-01/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_NTRO = &
                      (/9.99720D-01, 9.99989D-01, 9.99993D-01, 9.99994D-01, &
                        9.99995D-01, 9.99996D-01, 9.99999D-01, 1.00000D+00, &
                        1.00000D+00, 7.39413D-01, 7.39381D-01, 7.39381D-01, &
                        7.39381D-01, 7.39381D-01, 7.39349D-01, 7.39349D-01, &
                        7.03119D-01, 7.50131D-01, 7.50131D-01/)
  real(rk8) , dimension(nspi) , parameter :: w_ac_MFNC = &
                      (/9.99720D-01, 9.99989D-01, 9.99993D-01, 9.99994D-01, &
                        9.99995D-01, 9.99996D-01, 9.99999D-01, 1.00000D+00, &
                        1.00000D+00, 7.39413D-01, 7.39381D-01, 7.39381D-01, &
                        7.39381D-01, 7.39381D-01, 7.39349D-01, 7.39349D-01, &
                        7.03119D-01, 7.50131D-01, 7.50131D-01/)
  real(rk8) , dimension(nsoa, nspi) , parameter :: wa_soa = &
         reshape((/ w_nu_CNYL , w_nu_STHC , w_nu_OLEF , w_nu_AROM ,   &
                    w_nu_NTRO , w_nu_MFNC , w_ac_CNYL , w_ac_STHC ,   &
                    w_ac_OLEF , w_ac_AROM , w_ac_NTRO , w_ac_MFNC /)  &
                 ,(/nsoa, nspi/), order =(/2,1/) )

  ! Define all g values and organise appropriately in matrix
  real(rk8) , dimension(nspi) , parameter :: g_nu_CNYL = &
                      (/2.11953D-02, 1.58107D-02, 1.40388D-02, 1.30301D-02, &
                        1.21276D-02, 1.13166D-02, 9.51234D-03, 4.51790D-03, &
                        2.24077D-03, 3.02016D-04, 3.01551D-04, 3.01551D-04, &
                        3.01551D-04, 3.01551D-04, 3.01089D-04, 3.01089D-04, &
                        1.30998D-04, 5.30108D-05, 5.30108D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_nu_STHC = &
                      (/1.99149D-02, 1.50304D-02, 1.33905D-02, 1.24507D-02, &
                        1.16062D-02, 1.08447D-02, 9.14107D-03, 4.36654D-03, &
                        2.17034D-03, 2.92418D-04, 2.91969D-04, 2.91969D-04, &
                        2.91969D-04, 2.91969D-04, 2.91522D-04, 2.91522D-04, &
                        1.28266D-04, 5.10933D-05, 5.10933D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_nu_OLEF = &
                      (/2.11953D-02, 1.58107D-02, 1.40388D-02, 1.30301D-02, &
                        1.21276D-02, 1.13166D-02, 9.51234D-03, 4.51790D-03, &
                        2.24077D-03, 3.02016D-04, 3.01551D-04, 3.01551D-04, &
                        3.01551D-04, 3.01551D-04, 3.01089D-04, 3.01089D-04, &
                        1.30998D-04, 5.30108D-05, 5.30108D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_nu_AROM = &
                      (/2.12054D-02, 1.60054D-02, 1.42586D-02, 1.32573D-02, &
                        1.23574D-02, 1.15462D-02, 9.73719D-03, 4.65567D-03, &
                        2.30752D-03, 3.10716D-04, 3.10238D-04, 3.10238D-04, &
                        3.10238D-04, 3.10238D-04, 3.09763D-04, 3.09763D-04, &
                        1.36925D-04, 5.43943D-05, 5.43943D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_nu_NTRO = &
                      (/2.11953D-02, 1.58107D-02, 1.40388D-02, 1.30301D-02, &
                        1.21276D-02, 1.13166D-02, 9.51234D-03, 4.51790D-03, &
                        2.24077D-03, 3.02016D-04, 3.01551D-04, 3.01551D-04, &
                        3.01551D-04, 3.01551D-04, 3.01089D-04, 3.01089D-04, &
                        1.30998D-04, 5.30108D-05, 5.30108D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_nu_MFNC = &
                      (/2.11953D-02, 1.58107D-02, 1.40388D-02, 1.30301D-02, &
                        1.21276D-02, 1.13166D-02, 9.51234D-03, 4.51790D-03, &
                        2.24077D-03, 3.02016D-04, 3.01551D-04, 3.01551D-04, &
                        3.01551D-04, 3.01551D-04, 3.01089D-04, 3.01089D-04, &
                        1.30998D-04, 5.30108D-05, 5.30108D-05/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_CNYL = &
                      (/5.02088D-01, 6.57026D-01, 6.90758D-01, 7.28837D-01, &
                        7.43852D-01, 7.46850D-01, 7.68431D-01, 7.46481D-01, &
                        6.54579D-01, 1.29052D-01, 1.28932D-01, 1.28932D-01, &
                        1.28932D-01, 1.28932D-01, 1.28812D-01, 1.28812D-01, &
                        5.67142D-02, 2.30474D-02, 2.30474D-02/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_STHC = &
                      (/7.83208D-01, 8.10112D-01, 8.25102D-01, 8.27398D-01, &
                        8.24953D-01, 8.23172D-01, 8.24616D-01, 7.64935D-01, &
                        6.75659D-01, 1.26210D-01, 1.26084D-01, 1.26084D-01, &
                        1.26084D-01, 1.26084D-01, 1.25958D-01, 1.25958D-01, &
                        5.56438D-02, 2.22554D-02, 2.22554D-02/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_OLEF = &
                      (/5.02088D-01, 6.57026D-01, 6.90758D-01, 7.28837D-01, &
                        7.43852D-01, 7.46850D-01, 7.68431D-01, 7.46481D-01, &
                        6.54579D-01, 1.29052D-01, 1.28932D-01, 1.28932D-01, &
                        1.28932D-01, 1.28932D-01, 1.28812D-01, 1.28812D-01, &
                        5.67142D-02, 2.30474D-02, 2.30474D-02/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_AROM = &
                      (/4.81473D-01, 6.08802D-01, 6.17712D-01, 6.67687D-01, &
                        6.95697D-01, 6.98189D-01, 7.24141D-01, 7.25906D-01, &
                        6.38508D-01, 1.31356D-01, 1.31239D-01, 1.31239D-01, &
                        1.31239D-01, 1.31239D-01, 1.31122D-01, 1.31122D-01, &
                        5.90147D-02, 2.36228D-02, 2.36228D-02/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_NTRO = &
                      (/5.02088D-01, 6.57026D-01, 6.90758D-01, 7.28837D-01, &
                        7.43852D-01, 7.46850D-01, 7.68431D-01, 7.46481D-01, &
                        6.54579D-01, 1.29052D-01, 1.28932D-01, 1.28932D-01, &
                        1.28932D-01, 1.28932D-01, 1.28812D-01, 1.28812D-01, &
                        5.67142D-02, 2.30474D-02, 2.30474D-02/)
  real(rk8) , dimension(nspi) , parameter :: g_ac_MFNC = &
                      (/5.02088D-01, 6.57026D-01, 6.90758D-01, 7.28837D-01, &
                        7.43852D-01, 7.46850D-01, 7.68431D-01, 7.46481D-01, &
                        6.54579D-01, 1.29052D-01, 1.28932D-01, 1.28932D-01, &
                        1.28932D-01, 1.28932D-01, 1.28812D-01, 1.28812D-01, &
                        5.67142D-02, 2.30474D-02, 2.30474D-02/)
  real(rk8) , dimension(nsoa, nspi) , parameter :: ga_soa = &
         reshape((/ g_nu_CNYL , g_nu_STHC , g_nu_OLEF , g_nu_AROM ,   &
                    g_nu_NTRO , g_nu_MFNC , g_ac_CNYL , g_ac_STHC ,   &
                    g_ac_OLEF , g_ac_AROM , g_ac_NTRO , g_ac_MFNC /)  &
                 ,(/nsoa, nspi/), order =(/2,1/) )

  contains

  subroutine soa_aggregate(cat,chiin_n,chiin_a,numtr,numsoa,chiout)

    implicit none

    integer(ik4) , intent(in) :: numsoa , numtr
    integer(ik4) :: n , i , j , k
    integer(ik4) , dimension(numtr) , intent(in)  :: cat
    real(rk8) , dimension(numtr ) , intent(in)    :: chiin_n , chiin_a
    real(rk8) , dimension(numsoa) , intent(inout) :: chiout
    real(rk8) , dimension(numsoa) :: chivar

    chivar(:) = d_zero
    do n=1, numtr
      if ( cat(n) == 2 ) then       ! 2 = carbonyls / dicarbonyls
        chivar(icnylnu-3*numtr) = chivar(icnylnu-3*numtr) + chiin_n(n)
        chivar(icnylac-3*numtr) = chivar(icnylac-3*numtr) + chiin_a(n)
      else if ( cat(n) == 3 ) then  ! 3 = alkanes / alcohols
        chivar(isthcnu-3*numtr) = chivar(isthcnu-3*numtr) + chiin_n(n)
        chivar(isthcac-3*numtr) = chivar(isthcac-3*numtr) + chiin_a(n)
      else if ( cat(n) == 4 ) then  ! 4 = olefins
        chivar(iolefnu-3*numtr) = chivar(iolefnu-3*numtr) + chiin_n(n)
        chivar(iolefac-3*numtr) = chivar(iolefac-3*numtr) + chiin_a(n)
      else if ( cat(n) == 6 ) then  ! 6 = aromatics
        chivar(iaromnu-3*numtr) = chivar(iaromnu-3*numtr) + chiin_n(n)
        chivar(iaromac-3*numtr) = chivar(iaromac-3*numtr) + chiin_a(n)
      else if ( cat(n) == 9 ) then  ! 9 = organic nitrates
        chivar(intronu-3*numtr) = chivar(intronu-3*numtr) + chiin_n(n)
        chivar(introac-3*numtr) = chivar(introac-3*numtr) + chiin_a(n)
      else if ( cat(n) == 10 ) then ! 10 = oraganic multifunctionals
        chivar(imfncnu-3*numtr) = chivar(imfncnu-3*numtr) + chiin_n(n)
        chivar(imfncac-3*numtr) = chivar(imfncac-3*numtr) + chiin_a(n)
      end if
    end do

    do n=1, numsoa
      chiout(n) = chivar(n)
    end do

  end subroutine soa_aggregate 


end module mod_che_soa_cb6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
