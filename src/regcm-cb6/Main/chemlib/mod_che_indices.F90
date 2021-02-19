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

module mod_che_indices

  use mod_intkinds
  use mod_cbmz_parameters , only : nvar
  use mod_cb6_parameters , only : nvar_cb6

  implicit none

  public
! declarartoin of usefull chemical indices for species
! IMPORTANT : "INTERFACE SPECIES" indices
!    ibchl , ibchb , iochl , iochb , iisop, ianh4, iano3
!    now declared in mod_runparam for surface/chem  interface compatibility !!
  integer(ik4) :: iso2 , iso4 , idms
  integer(ik4) :: imsa
  integer(ik4) :: io3 , ino , ino2 , ino3 , ioh , iho2 , ih2o2
  integer(ik4) :: ihno2 , ihno3 , ihno4
  integer(ik4) :: isulf , ih2so4 , ihono , in2o5 , ihc , ihcr , ic2h4
  integer(ik4) :: ico , ihcho , iald2 , ieth , ic2h6 , ic3h8,ic3h6
  integer(ik4) :: itol , ixyl , inh3 , ipan , in2o
  integer(ik4) :: irooh , iaone , ibenz , ich4 , ico2
  integer(ik4) :: inox , ihox , isox , ieoh , ich3oh , iaco2 , ircooh,ihcooh
  integer(ik4) :: ipar , iolet , iolei , imgly , icres , iopen , iisoprd,iisopn
  integer(ik4) :: iethooh , ixo2 , iro2
  integer(ik4) :: iapin , ilimo
  integer(ik4) :: ialk4, ialk7
! integer(ik4) :: ianh4, iano3 : now declared in mod_runparam

  !CB6C tracers
  integer(ik4) :: intr,  isdio, ietha, iethy, idnpo, iepox, ietoh
  integer(ik4) :: iprpa, iket,  itoln, ixyln, ihpld, ipacn, ipacd, intr2
  integer(ik4) :: ipna,  imeoh, imepx, iopan, icat1, ihpox, iispx
  integer(ik4) :: ifacd, irpox, intr1, iacet, iintr, icron, iaacd, iethe
  integer(ik4) :: icmon, iterp, icrsl, iispr, intrc, iglyd, igly,  ixopn
  integer(ik4) :: imegy, iropn, iiole, iform, iole,  iaald, iispd, ialdx
  integer(ik4) :: ialka, iozn,  inmox, indox, intox, ipanx, iech4, ixo2h
  integer(ik4) :: icdox, ixo2r, imtha
  ! transporting all species !
  integer(ik4) :: iosng, ihco3, icrer, ibzo2, iror,  itolr, ixlo2
  integer(ik4) :: iepx2, ixylr, iopo3, ixo2n, iiso2, imeo2, iroo
  integer(ik4) :: ipox,  icxo3, iacoo, ioxy
  ! SORGAM - nucleation species
  integer(ik4) :: inntr,  insdio, inetha, inethy, indnpo, inepox, inetoh
  integer(ik4) :: inprpa, inket,  intoln, inxyln, inhpld, inpacn, inpacd, inntr2
  integer(ik4) :: inpna,  inmeoh, inmepx, inopan, incat1, inhpox, inispx
  integer(ik4) :: infacd, inrpox, inntr1, inacet, inintr, incron, inaacd, inethe
  integer(ik4) :: incmon, interp, incrsl, inispr, inntrc, inglyd, ingly,  inxopn
  integer(ik4) :: inmegy, inropn, iniole, inform, inole,  inaald, inispd, inaldx
  integer(ik4) :: inalka, inozn,  innmox, inndox, inntox, inpanx, inech4, inxo2h
  integer(ik4) :: inosng, inhco3, increr, inbzo2, inror,  intolr, inxlo2
  integer(ik4) :: inepx2, inxylr, inopo3, inxo2n, iniso2, inmeo2, inroo
  integer(ik4) :: inpox,  incxo3, inacoo, inoxy,  insulf, inbenz, inhono, inhox
  integer(ik4) :: incdox, inmtha, inxo2r
  ! SORGAM - accumulation species
  integer(ik4) :: iantr,  iasdio, iaetha, iaethy, iadnpo, iaepox, iaetoh
  integer(ik4) :: iaprpa, iaket,  iatoln, iaxyln, iahpld, iapacn, iapacd, iantr2
  integer(ik4) :: iapna,  iameoh, iamepx, iaopan, iacat1, iahpox, iaispx
  integer(ik4) :: iafacd, iarpox, iantr1, iaacet, iaintr, iacron, iaaacd, iaethe
  integer(ik4) :: iacmon, iaterp, iacrsl, iaispr, iantrc, iaglyd, iagly,  iaxopn
  integer(ik4) :: iamegy, iaropn, iaiole, iaform, iaole,  iaaald, iaispd, iaaldx
  integer(ik4) :: iaalka, iaozn,  ianmox, iandox, iantox, iapanx, iaech4, iaxo2h
  integer(ik4) :: iaosng, iahco3, iacrer, iabzo2, iaror,  iatolr, iaxlo2
  integer(ik4) :: iaepx2, iaxylr, iaopo3, iaxo2n, iaiso2, iameo2, iaroo
  integer(ik4) :: iapox,  iacxo3, iaacoo, iaoxy,  iasulf, iabenz, iahono, iahox
  integer(ik4) :: iacdox, iamtha, iaxo2r
  ! SORGAM - aggragates
  integer(ik4) :: icnylnu, isthcnu, iolefnu, iaromnu, intronu, imfncnu
  integer(ik4) :: icnylac, isthcac, iolefac, iaromac, introac, imfncac

  !*** abt added from wetdep scheme
  integer(ik4) :: iisopno3 , ich3ooh , ihydrald , ihyac , ipooh
  integer(ik4) :: ic3h7ooh , ic2h5ooh
  integer(ik4) :: iisopooh , imacrooh , ipb , ionit , ich3coooh
  integer(ik4) :: ich3cocho , ixooh
  integer(ik4) :: ionitr , iglyald , imvk , imacr , isoa , inh4
  integer(ik4) :: inh4no3 , ich3cooh
  integer(ik4) :: iterpooh , itolooh , imekooh , ialkooh

  !
  integer(ik4) :: ipollen

  ! list and name of cbmz species : must be absolutely consistant with
  ! mod_cbmz_Parameters

  integer(ik4) , parameter :: totsp = nvar
  character(len=8),target, dimension(totsp) :: cbmzspec
  data  cbmzspec /'CO2',     & ! 1
                  'H2SO4',   & ! 2
                  'HCOOH',   & ! 3
                  'RCOOH',   & ! 4
                  'MSA',     & ! 5
                  'DUMMY',   & ! 6
                  'PAN',     & ! 7
                  'TOL',     & ! 8
                  'O1D',     & ! 9
                  'H2O2',    & ! 10
                  'SO2',     & ! 11
                  'XYL',     & ! 12
                  'CH4',     & ! 13
                  'C2H6',    & ! 14
                  'CRO',     & ! 15
                  'DMS',     & ! 16
                  'HNO4',    & ! 17
                  'H2',      & ! 18
                  'TO2',     & ! 19
                  'CH3OH',   & ! 20
                  'HNO2',    & ! 21
                  'CH3OOH',  & ! 22
                  'ETHOOH',  & ! 23
                  'N2O5',    & ! 24
                  'ETH',     & ! 25
                  'CRES',    & ! 26
                  'O3P',     & ! 27
                  'CO',      & ! 28
                  'HNO3',    & ! 29
                  'PAR',     & ! 30
                  'OPEN',    & ! 31
                  'ISOPN',   & ! 32
                  'ISOPP',   & ! 33
                  'ISOPO2',  & ! 34
                  'H2O',     & ! 35
                  'AONE',    & ! 36
                  'OLEI',    & ! 37
                  'ISOP',    & ! 38
                  'HCHO',    & ! 39
                  'OLET',    & ! 40
                  'XO2',     & ! 41
                  'MGLY',    & ! 42
                  'ETHP',    & ! 43
                  'NAP',     & ! 44
                  'ALD2',    & ! 45
                  'CH3O2',   & ! 46
                  'ISOPRD',  & ! 47
                  'ANO2',    & ! 48
                  'ROOH',    & ! 49
                  'RO2',     & ! 50
                  'ONIT',    & ! 51
                  'HO2',     & ! 52
                  'O3',      & ! 53
                  'OH',      & ! 54
                  'NO',      & ! 55
                  'NO2',     & ! 56
                  'NO3',     & ! 57
                  'C2O3'/      ! 58

  integer(ik4) , parameter :: ind_CO2 = 1
  integer(ik4) , parameter :: ind_H2SO4 = 2
  integer(ik4) , parameter :: ind_HCOOH = 3
  integer(ik4) , parameter :: ind_RCOOH = 4
  integer(ik4) , parameter :: ind_MSA = 5
  integer(ik4) , parameter :: ind_DUMMY = 6
  integer(ik4) , parameter :: ind_PAN = 7
  integer(ik4) , parameter :: ind_TOL = 8
  integer(ik4) , parameter :: ind_O1D = 9
  integer(ik4) , parameter :: ind_H2O2 = 10
  integer(ik4) , parameter :: ind_SO2 = 11
  integer(ik4) , parameter :: ind_XYL = 12
  integer(ik4) , parameter :: ind_CH4 = 13
  integer(ik4) , parameter :: ind_C2H6 = 14
  integer(ik4) , parameter :: ind_CRO = 15
  integer(ik4) , parameter :: ind_DMS = 16
  integer(ik4) , parameter :: ind_HNO4 = 17
  integer(ik4) , parameter :: ind_H2 = 18
  integer(ik4) , parameter :: ind_TO2 = 19
  integer(ik4) , parameter :: ind_CH3OH = 20
  integer(ik4) , parameter :: ind_HNO2 = 21
  integer(ik4) , parameter :: ind_CH3OOH = 22
  integer(ik4) , parameter :: ind_ETHOOH = 23
  integer(ik4) , parameter :: ind_N2O5 = 24
  integer(ik4) , parameter :: ind_ETH = 25
  integer(ik4) , parameter :: ind_CRES = 26
  integer(ik4) , parameter :: ind_O3P = 27
  integer(ik4) , parameter :: ind_CO = 28
  integer(ik4) , parameter :: ind_HNO3 = 29
  integer(ik4) , parameter :: ind_PAR = 30
  integer(ik4) , parameter :: ind_OPEN = 31
  integer(ik4) , parameter :: ind_ISOPN = 32
  integer(ik4) , parameter :: ind_ISOPP = 33
  integer(ik4) , parameter :: ind_ISOPO2 = 34
  integer(ik4) , parameter :: ind_H2O = 35
  integer(ik4) , parameter :: ind_AONE = 36
  integer(ik4) , parameter :: ind_OLEI = 37
  integer(ik4) , parameter :: ind_ISOP = 38
  integer(ik4) , parameter :: ind_HCHO = 39
  integer(ik4) , parameter :: ind_OLET = 40
  integer(ik4) , parameter :: ind_XO2 = 41
  integer(ik4) , parameter :: ind_MGLY = 42
  integer(ik4) , parameter :: ind_ETHP = 43
  integer(ik4) , parameter :: ind_NAP = 44
  integer(ik4) , parameter :: ind_ALD2 = 45
  integer(ik4) , parameter :: ind_CH3O2 = 46
  integer(ik4) , parameter :: ind_ISOPRD = 47
  integer(ik4) , parameter :: ind_ANO2 = 48
  integer(ik4) , parameter :: ind_ROOH = 49
  integer(ik4) , parameter :: ind_RO2 = 50
  integer(ik4) , parameter :: ind_ONIT = 51
  integer(ik4) , parameter :: ind_HO2 = 52
  integer(ik4) , parameter :: ind_O3 = 53
  integer(ik4) , parameter :: ind_OH = 54
  integer(ik4) , parameter :: ind_NO = 55
  integer(ik4) , parameter :: ind_NO2 = 56
  integer(ik4) , parameter :: ind_NO3 = 57
  integer(ik4) , parameter :: ind_C2O3 = 58

  ! indices for jvalues
  integer(ik4) , parameter :: jvO2 = 1
  integer(ik4) , parameter :: jvO3a = 2
  integer(ik4) , parameter :: jvO3b = 3
  integer(ik4) , parameter :: jvNO2 = 4
  integer(ik4) , parameter :: jvNO3a = 5
  integer(ik4) , parameter :: jvNO3b = 6
  integer(ik4) , parameter :: jvN2O5a = 7
  integer(ik4) , parameter :: jvN2O5b = 8
  integer(ik4) , parameter :: jvN2O = 9
  integer(ik4) , parameter :: jvHO2 = 10
  integer(ik4) , parameter :: jvH2O2 = 11
  integer(ik4) , parameter :: jvHNO2 = 12
  integer(ik4) , parameter :: jvHNO3 = 13
  integer(ik4) , parameter :: jvHNO4 = 14
  integer(ik4) , parameter :: jvCH2Oa = 15
  integer(ik4) , parameter :: jvCH2Ob = 16
  integer(ik4) , parameter :: jvCH3CHOa = 17
  integer(ik4) , parameter :: jvCH3CHOb = 18
  integer(ik4) , parameter :: jvCH3CHOc = 19
  integer(ik4) , parameter :: jvC2H5CHO = 20
  integer(ik4) , parameter :: jvCHOCHO = 21
  integer(ik4) , parameter :: jvCH3COCHO = 22
  integer(ik4) , parameter :: jvCH3COCH3 = 23
  integer(ik4) , parameter :: jvCH3OOH = 24
  integer(ik4) , parameter :: jvCH3ONO2 = 25
  integer(ik4) , parameter :: jvPAN = 26

  integer(ik4) , parameter :: totch = nvar_cb6 !76
  character(len=6),target, dimension(nvar_cb6) :: cb6spec
  data  cb6spec  /'SULF',   & ! 1
                  'CDOX',   & ! 2
                  'SDIO',   & ! 3
                  'OSNG',   & ! 4
                  'ETHY',   & ! 5
                  'DNPO',   & ! 6
                  'BENZ',   & ! 7
                  'EPOX',   & ! 8
                  'ETOH',   & ! 9
                  'KET',    & ! 10
                  'TOLN',   & ! 11
                  'XYLN',   & ! 12
                  'HPLD',   & ! 13
                  'PACN',   & ! 14
                  'PACD',   & ! 15
                  'PNA',    & ! 16
                  'HONO',   & ! 17
                  'MTHA',   & ! 18
                  'NTR2',   & ! 19
                  'MEPX',   & ! 20
                  'HPOX',   & ! 21
                  'OPAN',   & ! 22
                  'CAT1',   & ! 23
                  'ISPX',   & ! 24
                  'ETHA',   & ! 25
                  'FACD',   & ! 26
                  'PANX',   & ! 27
                  'HCO3',   & ! 28
                  'PRPA',   & ! 29
                  'MEOH',   & ! 30
                  'CRER',   & ! 31
                  'RPOX',   & ! 32
                  'NTR1',   & ! 33
                  'ACET',   & ! 34
                  'INTR',   & ! 35
                  'AACD',   & ! 36
                  'CRON',   & ! 37
                  'ROR',    & ! 38
                  'BZO2',   & ! 39
                  'ETHE',   & ! 40
                  'CMON',   & ! 41
                  'TOLR',   & ! 42
                  'XYLR',   & ! 43
                  'CRSL',   & ! 44
                  'TERP',   & ! 45
                  'ISPR',   & ! 46
                  'NTRC',   & ! 47
                  'EPX2',   & ! 48
                  'GLYD',   & ! 49
                  'GLY',    & ! 50
                  'XOPN',   & ! 51
                  'MEGY',   & ! 52
                  'ROPN',   & ! 53
                  'IOLE',   & ! 54
                  'OLE',    & ! 55
                  'FORM',   & ! 56
                  'AALD',   & ! 57
                  'XO2R',   & ! 58
                  'OPO3',   & ! 59
                  'XO2N',   & ! 60
                  'ISO2',   & ! 61
                  'XO2H',   & ! 62
                  'ISPD',   & ! 63
                  'MEO2',   & ! 64
                  'ALKA',   & ! 65
                  'ALDX',   & ! 66
                  'HOX',    & ! 67
                  'ROO',    & ! 68
                  'NTOX',   & ! 69
                  'O',      & ! 70
                  'POX',    & ! 71
                  'CXO3',   & ! 72
                  'ACOO',   & ! 73
                  'OZN',    & ! 74
                  'NMOX',   & ! 75
                  'NDOX'/     ! 76

  integer(ik4) , parameter :: ind_SULF = 1
  integer(ik4) , parameter :: ind_CDOX = 2
  integer(ik4) , parameter :: ind_SDIO = 3
  integer(ik4) , parameter :: ind_OSNG = 4
  integer(ik4) , parameter :: ind_ETHY = 5
  integer(ik4) , parameter :: ind_DNPO = 6
  integer(ik4) , parameter :: ind_BENZ = 7
  integer(ik4) , parameter :: ind_EPOX = 8
  integer(ik4) , parameter :: ind_ETOH = 9
  integer(ik4) , parameter :: ind_KET  = 10
  integer(ik4) , parameter :: ind_TOLN = 11
  integer(ik4) , parameter :: ind_XYLN = 12
  integer(ik4) , parameter :: ind_HPLD = 13
  integer(ik4) , parameter :: ind_PACN = 14
  integer(ik4) , parameter :: ind_PACD = 15
  integer(ik4) , parameter :: ind_PNA  = 16
  integer(ik4) , parameter :: ind_HONO = 17
  integer(ik4) , parameter :: ind_MTHA = 18
  integer(ik4) , parameter :: ind_NTR2 = 19
  integer(ik4) , parameter :: ind_MEPX = 20
  integer(ik4) , parameter :: ind_HPOX = 21
  integer(ik4) , parameter :: ind_OPAN = 22
  integer(ik4) , parameter :: ind_CAT1 = 23
  integer(ik4) , parameter :: ind_ISPX = 24
  integer(ik4) , parameter :: ind_ETHA = 25
  integer(ik4) , parameter :: ind_FACD = 26
  integer(ik4) , parameter :: ind_PANX = 27
  integer(ik4) , parameter :: ind_HCO3 = 28
  integer(ik4) , parameter :: ind_PRPA = 29
  integer(ik4) , parameter :: ind_MEOH = 30
  integer(ik4) , parameter :: ind_CRER = 31
  integer(ik4) , parameter :: ind_RPOX = 32
  integer(ik4) , parameter :: ind_NTR1 = 33
  integer(ik4) , parameter :: ind_ACET = 34
  integer(ik4) , parameter :: ind_INTR = 35
  integer(ik4) , parameter :: ind_AACD = 36
  integer(ik4) , parameter :: ind_CRON = 37
  integer(ik4) , parameter :: ind_ROR  = 38
  integer(ik4) , parameter :: ind_BZO2 = 39
  integer(ik4) , parameter :: ind_ETHE = 40
  integer(ik4) , parameter :: ind_CMON = 41
  integer(ik4) , parameter :: ind_TOLR = 42
  integer(ik4) , parameter :: ind_XYLR = 43
  integer(ik4) , parameter :: ind_CRSL = 44
  integer(ik4) , parameter :: ind_TERP = 45
  integer(ik4) , parameter :: ind_ISPR = 46
  integer(ik4) , parameter :: ind_NTRC = 47
  integer(ik4) , parameter :: ind_EPX2 = 48
  integer(ik4) , parameter :: ind_GLYD = 49
  integer(ik4) , parameter :: ind_GLY  = 50
  integer(ik4) , parameter :: ind_XOPN = 51
  integer(ik4) , parameter :: ind_MEGY = 52
  integer(ik4) , parameter :: ind_ROPN = 53
  integer(ik4) , parameter :: ind_IOLE = 54
  integer(ik4) , parameter :: ind_OLE  = 55
  integer(ik4) , parameter :: ind_FORM = 56
  integer(ik4) , parameter :: ind_AALD = 57
  integer(ik4) , parameter :: ind_XO2R = 58
  integer(ik4) , parameter :: ind_OPO3 = 59
  integer(ik4) , parameter :: ind_XO2N = 60
  integer(ik4) , parameter :: ind_ISO2 = 61
  integer(ik4) , parameter :: ind_XO2H = 62
  integer(ik4) , parameter :: ind_ISPD = 63
  integer(ik4) , parameter :: ind_MEO2 = 64
  integer(ik4) , parameter :: ind_ALKA = 65
  integer(ik4) , parameter :: ind_ALDX = 66
  integer(ik4) , parameter :: ind_HOX  = 67
  integer(ik4) , parameter :: ind_ROO  = 68
  integer(ik4) , parameter :: ind_NTOX = 69
  integer(ik4) , parameter :: ind_O    = 70
  integer(ik4) , parameter :: ind_POX  = 71
  integer(ik4) , parameter :: ind_CXO3 = 72
  integer(ik4) , parameter :: ind_ACOO = 73
  integer(ik4) , parameter :: ind_OZN  = 74
  integer(ik4) , parameter :: ind_NMOX = 75
  integer(ik4) , parameter :: ind_NDOX = 76

  ! indices for jvalues
  integer(ik4) , parameter :: jvO31D  = 2
  integer(ik4) , parameter :: jvO33P  = 3
  integer(ik4) , parameter :: jvNDOX  = 4
  integer(ik4) , parameter :: jvNTOXa = 5
  integer(ik4) , parameter :: jvNTOXb = 6
  integer(ik4) , parameter :: jvDNPOb = 8
  integer(ik4) , parameter :: jvHPOX  = 11
  integer(ik4) , parameter :: jvHONO  = 12
  integer(ik4) , parameter :: jvNTRC  = 13
  integer(ik4) , parameter :: jvPNA   = 14
  integer(ik4) , parameter :: jvFORM  = 15
  integer(ik4) , parameter :: jvAALD  = 17
  integer(ik4) , parameter :: jvISPD  = 18
  integer(ik4) , parameter :: jvALDX  = 20
  integer(ik4) , parameter :: jvGLY   = 21
  integer(ik4) , parameter :: jvMEGY  = 22
  integer(ik4) , parameter :: jvACET  = 23
  integer(ik4) , parameter :: jvMEPX  = 24
  integer(ik4) , parameter :: jvNTR   = 25
  integer(ik4) , parameter :: jvPACN  = 26
  integer(ik4) , parameter :: jvPANX  = 29
  integer(ik4) , parameter :: jvRPOX  = 30
  integer(ik4) , parameter :: jvGLYD  = 31
  integer(ik4) , parameter :: jvKET   = 32
  integer(ik4) , parameter :: jvHPLD  = 33
  integer(ik4) , parameter :: jvCRON  = 34
  integer(ik4) , parameter :: jvXOPN  = 35
  integer(ik4) , parameter :: jvROPN  = 36

  character(len=6),target, dimension(3*nvar_cb6+12) :: srgspec
  data  srgspec  /'CDOX',   & ! 1  ! start of gas species
                  'SULF',   & ! 2
                  'SDIO',   & ! 3
                  'OSNG',   & ! 4
                  'MTHA',   & ! 5
                  'ETHA',   & ! 6
                  'ETHY',   & ! 7
                  'DNPO',   & ! 8
                  'BENZ',   & ! 9
                  'EPOX',   & ! 10
                  'ETOH',   & ! 11
                  'PRPA',   & ! 12
                  'KET',    & ! 13
                  'TOLN',   & ! 14
                  'XYLN',   & ! 15
                  'HPLD',   & ! 16
                  'PACN',   & ! 17
                  'PACD',   & ! 18
                  'NTR2',   & ! 19
                  'PNA',    & ! 20
                  'MEOH',   & ! 21
                  'HONO',   & ! 22
                  'MEPX',   & ! 23
                  'OPAN',   & ! 24
                  'CAT1',   & ! 25
                  'HPOX',   & ! 26
                  'ISPX',   & ! 27
                  'FACD',   & ! 28
                  'PANX',   & ! 29
                  'HCO3',   & ! 30
                  'CRER',   & ! 31
                  'RPOX',   & ! 32
                  'NTR1',   & ! 33
                  'ACET',   & ! 34
                  'INTR',   & ! 35
                  'BZO2',   & ! 36
                  'CRON',   & ! 37
                  'AACD',   & ! 38
                  'ROR',    & ! 39
                  'TOLR',   & ! 40
                  'ETHE',   & ! 41
                  'CMON',   & ! 42
                  'XO2R',   & ! 43
                  'TERP',   & ! 44
                  'CRSL',   & ! 45
                  'ISPR',   & ! 46
                  'EPX2',   & ! 47
                  'NTRC',   & ! 48
                  'GLYD',   & ! 49
                  'GLY',    & ! 50
                  'XOPN',   & ! 51
                  'MEGY',   & ! 52
                  'ROPN',   & ! 53
                  'IOLE',   & ! 54
                  'FORM',   & ! 55
                  'OLE',    & ! 56
                  'AALD',   & ! 57
                  'XYLR',   & ! 58
                  'OPO3',   & ! 59
                  'XO2N',   & ! 60
                  'ISO2',   & ! 61
                  'XO2H',   & ! 62
                  'ISPD',   & ! 63
                  'MEO2',   & ! 64
                  'ALDX',   & ! 65
                  'ALKA',   & ! 66
                  'OZN',    & ! 67
                  'ROO',    & ! 68
                  'HOX',    & ! 69
                  'POX',    & ! 70
                  'CXO3',   & ! 71
                  'NTOX',   & ! 72
                  'ACOO',   & ! 73
                  'O',      & ! 74
                  'NMOX',   & ! 75
                  'NDOX',   & ! 76 ! end of gas species
                  'CDOXnu', & ! 77 ! start of nucleation aerosol species
                  'SULFnu', & ! 78
                  'SDIOnu', & ! 79
                  'OSNGnu', & ! 80
                  'MTHAnu', & ! 81
                  'ETHAnu', & ! 82
                  'ETHYnu', & ! 83
                  'DNPOnu', & ! 84
                  'BENZnu', & ! 85
                  'EPOXnu', & ! 86
                  'ETOHnu', & ! 87
                  'PRPAnu', & ! 88
                  'KETnu',  & ! 89
                  'TOLNnu', & ! 90
                  'XYLNnu', & ! 91
                  'HPLDnu', & ! 92
                  'PACNnu', & ! 93
                  'PACDnu', & ! 94
                  'NTR2nu', & ! 95
                  'PNAnu',  & ! 96
                  'MEOHnu', & ! 97
                  'HONOnu', & ! 98
                  'MEPXnu', & ! 99 
                  'OPANnu', & ! 100
                  'CAT1nu', & ! 101
                  'HPOXnu', & ! 102
                  'ISPXnu', & ! 103
                  'FACDnu', & ! 104
                  'PANXnu', & ! 105
                  'HCO3nu', & ! 106
                  'CRERnu', & ! 107
                  'RPOXnu', & ! 108
                  'NTR1nu', & ! 109
                  'ACETnu', & ! 110
                  'INTRnu', & ! 111
                  'BZO2nu', & ! 112
                  'CRONnu', & ! 113
                  'AACDnu', & ! 114
                  'RORnu',  & ! 115
                  'TOLRnu', & ! 116
                  'ETHEnu', & ! 117
                  'CMONnu', & ! 118
                  'XO2Rnu', & ! 119
                  'TERPnu', & ! 120
                  'CRSLnu', & ! 121
                  'ISPRnu', & ! 122
                  'EPX2nu', & ! 123
                  'NTRCnu', & ! 124
                  'GLYDnu', & ! 125
                  'GLYnu',  & ! 126
                  'XOPNnu', & ! 127
                  'MEGYnu', & ! 128
                  'ROPNnu', & ! 129
                  'IOLEnu', & ! 130
                  'FORMnu', & ! 131
                  'OLEnu',  & ! 132
                  'AALDnu', & ! 133
                  'XYLRnu', & ! 134
                  'OPO3nu', & ! 135
                  'XO2Nnu', & ! 136
                  'ISO2nu', & ! 137
                  'XO2Hnu', & ! 138
                  'ISPDnu', & ! 139
                  'MEO2nu', & ! 140
                  'ALDXnu', & ! 141
                  'ALKAnu', & ! 142
                  'OZNnu',  & ! 143
                  'ROOnu',  & ! 144
                  'HOXnu',  & ! 145
                  'POXnu',  & ! 146
                  'CXO3nu', & ! 147
                  'NTOXnu', & ! 148
                  'ACOOnu', & ! 149
                  'Onu',    & ! 150
                  'NMOXnu', & ! 151
                  'NDOXnu', & ! 152! end of nucleation aerosol species
                  'CDOXac', & ! 153! start of accumulation aerosol species
                  'SULFac', & ! 154
                  'SDIOac', & ! 155
                  'OSNGac', & ! 156
                  'MTHAac', & ! 157
                  'ETHAac', & ! 158
                  'ETHYac', & ! 159
                  'DNPOac', & ! 160
                  'BENZac', & ! 161
                  'EPOXac', & ! 162
                  'ETOHac', & ! 163
                  'PRPAac', & ! 164
                  'KETac',  & ! 165
                  'TOLNac', & ! 166
                  'XYLNac', & ! 167
                  'HPLDac', & ! 168
                  'PACNac', & ! 169
                  'PACDac', & ! 170
                  'NTR2ac', & ! 171
                  'PNAac',  & ! 172
                  'MEOHac', & ! 173
                  'HONOac', & ! 174
                  'MEPXac', & ! 175
                  'OPANac', & ! 176
                  'CAT1ac', & ! 177
                  'HPOXac', & ! 178
                  'ISPXac', & ! 179
                  'FACDac', & ! 180
                  'PANXac', & ! 181
                  'HCO3ac', & ! 182
                  'CRERac', & ! 183
                  'RPOXac', & ! 184
                  'NTR1ac', & ! 185
                  'ACETac', & ! 186
                  'INTRac', & ! 187
                  'BZO2ac', & ! 188
                  'CRONac', & ! 189
                  'AACDac', & ! 190
                  'RORac',  & ! 191
                  'TOLRac', & ! 192
                  'ETHEac', & ! 193
                  'CMONac', & ! 194
                  'XO2Rac', & ! 195
                  'TERPac', & ! 196
                  'CRSLac', & ! 197
                  'ISPRac', & ! 198
                  'EPX2ac', & ! 199
                  'NTRCac', & ! 200
                  'GLYDac', & ! 201
                  'GLYac',  & ! 202
                  'XOPNac', & ! 203
                  'MEGYac', & ! 204
                  'ROPNac', & ! 205
                  'IOLEac', & ! 206
                  'FORMac', & ! 207
                  'OLEac',  & ! 208
                  'AALDac', & ! 209
                  'XYLRac', & ! 210
                  'OPO3ac', & ! 211
                  'XO2Nac', & ! 212
                  'ISO2ac', & ! 213
                  'XO2Hac', & ! 214
                  'ISPDac', & ! 215
                  'MEO2ac', & ! 216
                  'ALDXac', & ! 217
                  'ALKAac', & ! 218
                  'OZNac',  & ! 219
                  'ROOac',  & ! 220
                  'HOXac',  & ! 221
                  'POXac',  & ! 222
                  'CXO3ac', & ! 223
                  'NTOXac', & ! 224
                  'ACOOac', & ! 225
                  'Oac',    & ! 226
                  'NMOXac', & ! 227
                  'NDOXac', & ! 228! end of accumulation aerosol species
                  'CNYLnu', & ! 229! 1  ! start of nuc aggregates
                  'STHCnu', & ! 230! 2
                  'OLEFnu', & ! 231! 3
                  'AROMnu', & ! 232! 4
                  'NTROnu', & ! 233! 5
                  'MFNCnu', & ! 234! 6  
                  'CNYLac', & ! 235! 7  ! start of acc aggregates
                  'STHCac', & ! 236! 8
                  'OLEFac', & ! 237! 9
                  'AROMac', & ! 238! 10
                  'NTROac', & ! 239! 11
                  'MFNCac'/   ! 240! 12

! nucleation aerosol indices
  integer(ik4) , parameter :: idn_CDOX = nvar_cb6+ind_CDOX
  integer(ik4) , parameter :: idn_SULF = nvar_cb6+ind_SULF
  integer(ik4) , parameter :: idn_SDIO = nvar_cb6+ind_SDIO
  integer(ik4) , parameter :: idn_OSNG = nvar_cb6+ind_OSNG
  integer(ik4) , parameter :: idn_MTHA = nvar_cb6+ind_MTHA
  integer(ik4) , parameter :: idn_ETHA = nvar_cb6+ind_ETHA
  integer(ik4) , parameter :: idn_ETHY = nvar_cb6+ind_ETHY
  integer(ik4) , parameter :: idn_DNPO = nvar_cb6+ind_DNPO
  integer(ik4) , parameter :: idn_BENZ = nvar_cb6+ind_BENZ
  integer(ik4) , parameter :: idn_EPOX = nvar_cb6+ind_EPOX
  integer(ik4) , parameter :: idn_ETOH = nvar_cb6+ind_ETOH
  integer(ik4) , parameter :: idn_PRPA = nvar_cb6+ind_PRPA
  integer(ik4) , parameter :: idn_KET  = nvar_cb6+ind_KET
  integer(ik4) , parameter :: idn_TOLN = nvar_cb6+ind_TOLN
  integer(ik4) , parameter :: idn_XYLN = nvar_cb6+ind_XYLN
  integer(ik4) , parameter :: idn_HPLD = nvar_cb6+ind_HPLD
  integer(ik4) , parameter :: idn_PACN = nvar_cb6+ind_PACN
  integer(ik4) , parameter :: idn_PACD = nvar_cb6+ind_PACD
  integer(ik4) , parameter :: idn_NTR2 = nvar_cb6+ind_NTR2
  integer(ik4) , parameter :: idn_PNA  = nvar_cb6+ind_PNA
  integer(ik4) , parameter :: idn_MEOH = nvar_cb6+ind_MEOH
  integer(ik4) , parameter :: idn_HONO = nvar_cb6+ind_HONO
  integer(ik4) , parameter :: idn_MEPX = nvar_cb6+ind_MEPX
  integer(ik4) , parameter :: idn_OPAN = nvar_cb6+ind_OPAN
  integer(ik4) , parameter :: idn_CAT1 = nvar_cb6+ind_CAT1
  integer(ik4) , parameter :: idn_HPOX = nvar_cb6+ind_HPOX
  integer(ik4) , parameter :: idn_ISPX = nvar_cb6+ind_ISPX
  integer(ik4) , parameter :: idn_FACD = nvar_cb6+ind_FACD
  integer(ik4) , parameter :: idn_PANX = nvar_cb6+ind_PANX
  integer(ik4) , parameter :: idn_HCO3 = nvar_cb6+ind_HCO3
  integer(ik4) , parameter :: idn_CRER = nvar_cb6+ind_CRER
  integer(ik4) , parameter :: idn_RPOX = nvar_cb6+ind_RPOX
  integer(ik4) , parameter :: idn_NTR1 = nvar_cb6+ind_NTR1
  integer(ik4) , parameter :: idn_ACET = nvar_cb6+ind_ACET
  integer(ik4) , parameter :: idn_INTR = nvar_cb6+ind_INTR
  integer(ik4) , parameter :: idn_BZO2 = nvar_cb6+ind_BZO2
  integer(ik4) , parameter :: idn_CRON = nvar_cb6+ind_CRON
  integer(ik4) , parameter :: idn_AACD = nvar_cb6+ind_AACD
  integer(ik4) , parameter :: idn_ROR  = nvar_cb6+ind_ROR
  integer(ik4) , parameter :: idn_TOLR = nvar_cb6+ind_TOLR
  integer(ik4) , parameter :: idn_ETHE = nvar_cb6+ind_ETHE
  integer(ik4) , parameter :: idn_CMON = nvar_cb6+ind_CMON
  integer(ik4) , parameter :: idn_XO2R = nvar_cb6+ind_XO2R
  integer(ik4) , parameter :: idn_TERP = nvar_cb6+ind_TERP
  integer(ik4) , parameter :: idn_CRSL = nvar_cb6+ind_CRSL
  integer(ik4) , parameter :: idn_ISPR = nvar_cb6+ind_ISPR
  integer(ik4) , parameter :: idn_EPX2 = nvar_cb6+ind_EPX2
  integer(ik4) , parameter :: idn_NTRC = nvar_cb6+ind_NTRC
  integer(ik4) , parameter :: idn_GLYD = nvar_cb6+ind_GLYD
  integer(ik4) , parameter :: idn_GLY  = nvar_cb6+ind_GLY
  integer(ik4) , parameter :: idn_XOPN = nvar_cb6+ind_XOPN
  integer(ik4) , parameter :: idn_MEGY = nvar_cb6+ind_MEGY
  integer(ik4) , parameter :: idn_ROPN = nvar_cb6+ind_ROPN
  integer(ik4) , parameter :: idn_IOLE = nvar_cb6+ind_IOLE
  integer(ik4) , parameter :: idn_FORM = nvar_cb6+ind_FORM
  integer(ik4) , parameter :: idn_OLE  = nvar_cb6+ind_OLE
  integer(ik4) , parameter :: idn_AALD = nvar_cb6+ind_AALD
  integer(ik4) , parameter :: idn_XYLR = nvar_cb6+ind_XYLR
  integer(ik4) , parameter :: idn_OPO3 = nvar_cb6+ind_OPO3
  integer(ik4) , parameter :: idn_XO2N = nvar_cb6+ind_XO2N
  integer(ik4) , parameter :: idn_ISO2 = nvar_cb6+ind_ISO2
  integer(ik4) , parameter :: idn_XO2H = nvar_cb6+ind_XO2H
  integer(ik4) , parameter :: idn_ISPD = nvar_cb6+ind_ISPD
  integer(ik4) , parameter :: idn_MEO2 = nvar_cb6+ind_MEO2
  integer(ik4) , parameter :: idn_ALDX = nvar_cb6+ind_ALDX
  integer(ik4) , parameter :: idn_ALKA = nvar_cb6+ind_ALKA
  integer(ik4) , parameter :: idn_OZN  = nvar_cb6+ind_OZN
  integer(ik4) , parameter :: idn_ROO  = nvar_cb6+ind_ROO
  integer(ik4) , parameter :: idn_HOX  = nvar_cb6+ind_HOX
  integer(ik4) , parameter :: idn_POX  = nvar_cb6+ind_POX
  integer(ik4) , parameter :: idn_CXO3 = nvar_cb6+ind_CXO3
  integer(ik4) , parameter :: idn_NTOX = nvar_cb6+ind_NTOX
  integer(ik4) , parameter :: idn_ACOO = nvar_cb6+ind_ACOO
  integer(ik4) , parameter :: idn_O    = nvar_cb6+ind_O
  integer(ik4) , parameter :: idn_NMOX = nvar_cb6+ind_NMOX
  integer(ik4) , parameter :: idn_NDOX = nvar_cb6+ind_NDOX
! accumulation aerosol indices
  integer(ik4) , parameter :: ida_CDOX = 2*nvar_cb6+ind_CDOX
  integer(ik4) , parameter :: ida_SULF = 2*nvar_cb6+ind_SULF
  integer(ik4) , parameter :: ida_SDIO = 2*nvar_cb6+ind_SDIO
  integer(ik4) , parameter :: ida_OSNG = 2*nvar_cb6+ind_OSNG
  integer(ik4) , parameter :: ida_MTHA = 2*nvar_cb6+ind_MTHA
  integer(ik4) , parameter :: ida_ETHA = 2*nvar_cb6+ind_ETHA
  integer(ik4) , parameter :: ida_ETHY = 2*nvar_cb6+ind_ETHY
  integer(ik4) , parameter :: ida_DNPO = 2*nvar_cb6+ind_DNPO
  integer(ik4) , parameter :: ida_BENZ = 2*nvar_cb6+ind_BENZ
  integer(ik4) , parameter :: ida_EPOX = 2*nvar_cb6+ind_EPOX
  integer(ik4) , parameter :: ida_ETOH = 2*nvar_cb6+ind_ETOH
  integer(ik4) , parameter :: ida_PRPA = 2*nvar_cb6+ind_PRPA
  integer(ik4) , parameter :: ida_KET  = 2*nvar_cb6+ind_KET
  integer(ik4) , parameter :: ida_TOLN = 2*nvar_cb6+ind_TOLN
  integer(ik4) , parameter :: ida_XYLN = 2*nvar_cb6+ind_XYLN
  integer(ik4) , parameter :: ida_HPLD = 2*nvar_cb6+ind_HPLD
  integer(ik4) , parameter :: ida_PACN = 2*nvar_cb6+ind_PACN
  integer(ik4) , parameter :: ida_PACD = 2*nvar_cb6+ind_PACD
  integer(ik4) , parameter :: ida_NTR2 = 2*nvar_cb6+ind_NTR2
  integer(ik4) , parameter :: ida_PNA  = 2*nvar_cb6+ind_PNA
  integer(ik4) , parameter :: ida_MEOH = 2*nvar_cb6+ind_MEOH
  integer(ik4) , parameter :: ida_HONO = 2*nvar_cb6+ind_HONO
  integer(ik4) , parameter :: ida_MEPX = 2*nvar_cb6+ind_MEPX
  integer(ik4) , parameter :: ida_OPAN = 2*nvar_cb6+ind_OPAN
  integer(ik4) , parameter :: ida_CAT1 = 2*nvar_cb6+ind_CAT1
  integer(ik4) , parameter :: ida_HPOX = 2*nvar_cb6+ind_HPOX
  integer(ik4) , parameter :: ida_ISPX = 2*nvar_cb6+ind_ISPX
  integer(ik4) , parameter :: ida_FACD = 2*nvar_cb6+ind_FACD
  integer(ik4) , parameter :: ida_PANX = 2*nvar_cb6+ind_PANX
  integer(ik4) , parameter :: ida_HCO3 = 2*nvar_cb6+ind_HCO3
  integer(ik4) , parameter :: ida_CRER = 2*nvar_cb6+ind_CRER
  integer(ik4) , parameter :: ida_RPOX = 2*nvar_cb6+ind_RPOX
  integer(ik4) , parameter :: ida_NTR1 = 2*nvar_cb6+ind_NTR1
  integer(ik4) , parameter :: ida_ACET = 2*nvar_cb6+ind_ACET
  integer(ik4) , parameter :: ida_INTR = 2*nvar_cb6+ind_INTR
  integer(ik4) , parameter :: ida_BZO2 = 2*nvar_cb6+ind_BZO2
  integer(ik4) , parameter :: ida_CRON = 2*nvar_cb6+ind_CRON
  integer(ik4) , parameter :: ida_AACD = 2*nvar_cb6+ind_AACD
  integer(ik4) , parameter :: ida_ROR  = 2*nvar_cb6+ind_ROR
  integer(ik4) , parameter :: ida_TOLR = 2*nvar_cb6+ind_TOLR
  integer(ik4) , parameter :: ida_ETHE = 2*nvar_cb6+ind_ETHE
  integer(ik4) , parameter :: ida_CMON = 2*nvar_cb6+ind_CMON
  integer(ik4) , parameter :: ida_XO2R = 2*nvar_cb6+ind_XO2R
  integer(ik4) , parameter :: ida_TERP = 2*nvar_cb6+ind_TERP
  integer(ik4) , parameter :: ida_CRSL = 2*nvar_cb6+ind_CRSL
  integer(ik4) , parameter :: ida_ISPR = 2*nvar_cb6+ind_ISPR
  integer(ik4) , parameter :: ida_EPX2 = 2*nvar_cb6+ind_EPX2
  integer(ik4) , parameter :: ida_NTRC = 2*nvar_cb6+ind_NTRC
  integer(ik4) , parameter :: ida_GLYD = 2*nvar_cb6+ind_GLYD
  integer(ik4) , parameter :: ida_GLY  = 2*nvar_cb6+ind_GLY
  integer(ik4) , parameter :: ida_XOPN = 2*nvar_cb6+ind_XOPN
  integer(ik4) , parameter :: ida_MEGY = 2*nvar_cb6+ind_MEGY
  integer(ik4) , parameter :: ida_ROPN = 2*nvar_cb6+ind_ROPN
  integer(ik4) , parameter :: ida_IOLE = 2*nvar_cb6+ind_IOLE
  integer(ik4) , parameter :: ida_FORM = 2*nvar_cb6+ind_FORM
  integer(ik4) , parameter :: ida_OLE  = 2*nvar_cb6+ind_OLE
  integer(ik4) , parameter :: ida_AALD = 2*nvar_cb6+ind_AALD
  integer(ik4) , parameter :: ida_XYLR = 2*nvar_cb6+ind_XYLR
  integer(ik4) , parameter :: ida_OPO3 = 2*nvar_cb6+ind_OPO3
  integer(ik4) , parameter :: ida_XO2N = 2*nvar_cb6+ind_XO2N
  integer(ik4) , parameter :: ida_ISO2 = 2*nvar_cb6+ind_ISO2
  integer(ik4) , parameter :: ida_XO2H = 2*nvar_cb6+ind_XO2H
  integer(ik4) , parameter :: ida_ISPD = 2*nvar_cb6+ind_ISPD
  integer(ik4) , parameter :: ida_MEO2 = 2*nvar_cb6+ind_MEO2
  integer(ik4) , parameter :: ida_ALDX = 2*nvar_cb6+ind_ALDX
  integer(ik4) , parameter :: ida_ALKA = 2*nvar_cb6+ind_ALKA
  integer(ik4) , parameter :: ida_OZN  = 2*nvar_cb6+ind_OZN
  integer(ik4) , parameter :: ida_ROO  = 2*nvar_cb6+ind_ROO
  integer(ik4) , parameter :: ida_HOX  = 2*nvar_cb6+ind_HOX
  integer(ik4) , parameter :: ida_POX  = 2*nvar_cb6+ind_POX
  integer(ik4) , parameter :: ida_CXO3 = 2*nvar_cb6+ind_CXO3
  integer(ik4) , parameter :: ida_NTOX = 2*nvar_cb6+ind_NTOX
  integer(ik4) , parameter :: ida_ACOO = 2*nvar_cb6+ind_ACOO
  integer(ik4) , parameter :: ida_O    = 2*nvar_cb6+ind_O
  integer(ik4) , parameter :: ida_NMOX = 2*nvar_cb6+ind_NMOX
  integer(ik4) , parameter :: ida_NDOX = 2*nvar_cb6+ind_NDOX
  integer(ik4) , parameter :: indnu_CNYL = 3*nvar_cb6+1
  integer(ik4) , parameter :: indnu_STHC = 3*nvar_cb6+2
  integer(ik4) , parameter :: indnu_OLEF = 3*nvar_cb6+3
  integer(ik4) , parameter :: indnu_AROM = 3*nvar_cb6+4
  integer(ik4) , parameter :: indnu_NTRO = 3*nvar_cb6+5
  integer(ik4) , parameter :: indnu_MFNC = 3*nvar_cb6+6
  integer(ik4) , parameter :: indac_CNYL = 3*nvar_cb6+7
  integer(ik4) , parameter :: indac_STHC = 3*nvar_cb6+8
  integer(ik4) , parameter :: indac_OLEF = 3*nvar_cb6+9
  integer(ik4) , parameter :: indac_AROM = 3*nvar_cb6+10
  integer(ik4) , parameter :: indac_NTRO = 3*nvar_cb6+11
  integer(ik4) , parameter :: indac_MFNC = 3*nvar_cb6+12

end module mod_che_indices
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
