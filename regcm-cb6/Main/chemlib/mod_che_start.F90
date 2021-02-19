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

module mod_che_start

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_common
  use mod_che_indices
  use mod_che_bdyco
  use mod_che_wetdep
  use mod_che_chemistry
  use mod_che_carbonaer
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_bdyco
  use mod_che_dust
  use mod_che_sox
  use mod_che_seasalt
  use mod_mppparam
  use mod_che_hvread
  use mod_che_molwg
  use mod_che_bionit
  use mod_che_soa_cb6

  implicit none

  private

  real(rk8) , parameter :: solso4 = 0.9D0

  public  :: start_chem

  contains

  !--------------------------------------------------------------------------

  subroutine start_chem
    implicit none
    integer(ik4) :: i , j , k , itr , ibin , jbin , kbin, n 

    ! A : Intialise chemistry tracer indices

    dtchsolv = dtche
    kchsolv = ntche

    iso2  = 0
    iso4  = 0
    idms  = 0
    imsa  = 0
    ibchl = 0
    ibchb = 0
    iochl = 0
    iochb = 0
    idust = 0
    isslt = 0
    icarb = 0
    ianh4 = 0
    iano3 = 0


    ino     = 0
    ino2    = 0
    in2o5   = 0
    ihno2   = 0
    ihno3   = 0
    ihno4   = 0
    io3     = 0
    ih2o2   = 0
    ico     = 0
    iso2    = 0
    idms    = 0
    ih2so4  = 0
    ich4    = 0
    ic2h6   = 0
    ipar    = 0
    ich3oh  = 0
    ihcho   = 0
    iald2   = 0
    iaone   = 0
    ieth    = 0
    iolet   = 0
    iolei   = 0
    itol    = 0
    ixyl    = 0
    iisop   = 0
    ionit   = 0
    ipan    = 0
    ihcooh  = 0
    ircooh  = 0
    ich3ooh = 0
    iethooh = 0
    irooh   = 0
    imgly   = 0
    iisoprd = 0
    iisopn  = 0
    iopen   = 0
    icres   = 0
    ipollen = 0

! CBSIX tracer initialization
    intr    = 0
    isulf   = 0
    isdio   = 0
    ietha   = 0
    iethy   = 0
    idnpo   = 0
    ibenz   = 0
    iepox   = 0
    ietoh   = 0
    iprpa   = 0
    iket    = 0
    itoln   = 0
    ixyln   = 0
    ihpld   = 0
    ipacn   = 0
    ipacd   = 0
    intr2   = 0
    ipna    = 0
    imeoh   = 0
    ihono   = 0
    imepx   = 0
    iopan   = 0
    icat1   = 0
    ihpox   = 0
    iispx   = 0
    ifacd   = 0
    irpox   = 0
    intr1   = 0
    iacet   = 0
    iintr   = 0
    icron   = 0
    iaacd   = 0
    iethe   = 0
    icmon   = 0
    iterp   = 0
    icrsl   = 0
    iispr   = 0
    intrc   = 0
    iglyd   = 0
    igly    = 0
    ixopn   = 0
    imegy   = 0
    iropn   = 0
    iiole   = 0
    iform   = 0
    iole    = 0
    iaald   = 0
    iispd   = 0
    ialdx   = 0
    ialka   = 0
    iozn    = 0
    inmox   = 0
    indox   = 0
    intox   = 0
    ipanx   = 0
    icdox   = 0
    imtha   = 0
    ixo2r   = 0
    ! transporting all species !
    iosng   = 0
    iech4   = 0 ! to be removed later with srg
    ihco3   = 0
    icrer   = 0
    ibzo2   = 0
    iror    = 0
    itolr   = 0
    ixlo2   = 0
    iepx2   = 0
    ixylr   = 0
    iopo3   = 0
    ixo2n   = 0
    iiso2   = 0
    ixo2h   = 0
    imeo2   = 0
    iroo    = 0
    ihox    = 0
    ipox    = 0
    icxo3   = 0
    iacoo   = 0
    ioxy    = 0
    ! nucleation species
    inntr   = 0
    insulf  = 0
    insdio  = 0
    inetha  = 0
    inethy  = 0
    indnpo  = 0
    inbenz  = 0
    inepox  = 0
    inetoh  = 0
    inprpa  = 0
    inket   = 0
    intoln  = 0
    inxyln  = 0
    inhpld  = 0
    inpacn  = 0
    inpacd  = 0
    inntr2  = 0
    inpna   = 0
    inmeoh  = 0
    inhono  = 0
    inmepx  = 0
    inopan  = 0
    incat1  = 0
    inhpox  = 0
    inispx  = 0
    infacd  = 0
    inrpox  = 0
    inntr1  = 0
    inacet  = 0
    inintr  = 0
    incron  = 0
    inaacd  = 0
    inethe  = 0
    incmon  = 0
    interp  = 0
    incrsl  = 0
    inispr  = 0
    inntrc  = 0
    inglyd  = 0
    ingly   = 0
    inxopn  = 0
    inmegy  = 0
    inropn  = 0
    iniole  = 0
    inform  = 0
    inole   = 0
    inaald  = 0
    inispd  = 0
    inaldx  = 0
    inalka  = 0
    inozn   = 0
    innmox  = 0
    inndox  = 0
    inntox  = 0
    inpanx  = 0
    inosng  = 0
    inech4  = 0
    inhco3  = 0
    increr  = 0
    inbzo2  = 0
    inror   = 0
    intolr  = 0
    inxlo2  = 0
    inepx2  = 0
    inxylr  = 0
    inopo3  = 0
    inxo2n  = 0
    iniso2  = 0
    inxo2h  = 0
    inmeo2  = 0
    inroo   = 0
    inhox   = 0
    inpox   = 0
    incxo3  = 0
    inacoo  = 0
    inoxy   = 0
    incdox  = 0
    inmtha  = 0
    inxo2r  = 0
    ! accumulation species
    iantr   = 0
    iasulf  = 0
    iasdio  = 0
    iaetha  = 0
    iaethy  = 0
    iadnpo  = 0
    iabenz  = 0
    iaepox  = 0
    iaetoh  = 0
    iaprpa  = 0
    iaket   = 0
    iatoln  = 0
    iaxyln  = 0
    iahpld  = 0
    iapacn  = 0
    iapacd  = 0
    iantr2  = 0
    iapna   = 0
    iameoh  = 0
    iahono  = 0
    iamepx  = 0
    iaopan  = 0
    iacat1  = 0
    iahpox  = 0
    iaispx  = 0
    iafacd  = 0
    iarpox  = 0
    iantr1  = 0
    iaacet  = 0
    iaintr  = 0
    iacron  = 0
    iaaacd  = 0
    iaethe  = 0
    iacmon  = 0
    iaterp  = 0
    iacrsl  = 0
    iaispr  = 0
    iantrc  = 0
    iaglyd  = 0
    iagly   = 0
    iaxopn  = 0
    iamegy  = 0
    iaropn  = 0
    iaiole  = 0
    iaform  = 0
    iaole   = 0
    iaaald  = 0
    iaispd  = 0
    iaaldx  = 0
    iaalka  = 0
    iaozn   = 0
    ianmox  = 0
    iandox  = 0
    iantox  = 0
    iapanx  = 0
    iaosng  = 0
    iaech4  = 0
    iahco3  = 0
    iacrer  = 0
    iabzo2  = 0
    iaror   = 0
    iatolr  = 0
    iaxlo2  = 0
    iaepx2  = 0
    iaxylr  = 0
    iaopo3  = 0
    iaxo2n  = 0
    iaiso2  = 0
    iaxo2h  = 0
    iameo2  = 0
    iaroo   = 0
    iahox   = 0
    iapox   = 0
    iacxo3  = 0
    iaacoo  = 0
    iaoxy   = 0
    iacdox  = 0
    iamtha  = 0
    iaxo2r  = 0

! SORGAM tracer initialization
    icnylnu = 0
    isthcnu = 0
    iolefnu = 0
    iaromnu = 0
    intronu = 0
    imfncnu = 0
    icnylac = 0
    isthcac = 0
    iolefac = 0
    iaromac = 0
    introac = 0
    imfncac = 0

    !abt *** For initializing megan tracer biogenic voc mask
    !    *** Mask not equal to zero when any MEGAN species is
    !    *** defined as a tracer within regcm.in
    !    *** If not equal to zero then that compound will be
    !    *** used as a surface emission from MEGAN and not
    !    *** from inventory (see surftnd.F for use)

#if defined CLM45 || (defined VOC && defined CLM)
    if ( igaschem == 1 ) then
      bvoc_trmask(:) = 0
    end if
#endif

    ibin = 0
    jbin = 0
    kbin = 0

    trac%indcbmz(:) = -1
    trac%indcb6(:)  = -1
    trac%mw(:)  = d_zero
    do itr = 1 , ntr
      if ( chtrname(itr) == 'SO2' ) iso2 = itr
      if ( chtrname(itr) == 'DMS' ) idms = itr
!
! for sulfuric acid and sulfate aer, we define here the same tracer index for compatibility with gas-phase
! chemistry options, and simple sulfate aer option. this might change with adding explicit sulfate aq chemistry.
! we consider here that all h2so4 partition in aerosol phase
      if ( chtrname(itr) == 'SO4' .or. chtrname(itr) == 'H2SO4'    ) then
        ! sulfate index is added to carb vector for treatment in drydep
        ! and wetdep sulfate effective diameter and bin is taken equal to ochl
        kbin = kbin + 1
        ih2so4 = itr
        iso4 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl
        chtrsol(iso4) = solso4
      end if
      if ( chtrname(itr) == 'ANO3' ) then
        kbin = kbin + 1
        iano3 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl
        chtrsol(iano3) = solso4
      end if
      if ( chtrname(itr) == 'ANH4' ) then
        kbin = kbin + 1
        ianh4 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl
        chtrsol(ianh4) = solso4
      end if
      if ( chtrname(itr) == 'BC_HL' ) then
        kbin = kbin + 1
        ibchl = itr
        icarb(kbin) = itr
        carbed(kbin) = reffbchl
        chtrsol(itr) = solbchl
      end if
      if ( chtrname(itr) == 'BC_HB' ) then
        kbin = kbin + 1
        ibchb = itr
        icarb(kbin) = itr
        carbed(kbin) = reffbc
        chtrsol(itr) = solbc
      end if
      if ( chtrname(itr) == 'OC_HL' ) then
        kbin = kbin + 1
        iochl = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl
        chtrsol(itr) = soloc
      end if
      if ( chtrname(itr) == 'OC_HB' ) then
        kbin = kbin + 1
        iochb = itr
        icarb(kbin) = itr
        carbed(kbin) = reffoc
        chtrsol(itr) = solochl
      end if
      if ( chtrname(itr)(1:4) ==  'DUST') then
        ibin = ibin + 1
        idust(ibin) = itr
        chtrsol(itr) = soldust(ibin)
      end if
      if ( chtrname(itr)(1:4) ==  'SSLT') then
        jbin = jbin + 1
        isslt(jbin) = itr
        chtrsol(itr) = solsslt(jbin)
      end if

      if ( chemsimtype(1:4) ==  'CBMZ') then
      ! gas phas species (CBMZ),
      ! max configuration : number of tracer = number of species

      do n = 1,totsp
        if ( chtrname(itr) == cbmzspec(n) )then
          ! index of the tracer in the CBMZ list of species
          trac%indcbmz(itr) = n
          ! correponding molecular weight
          trac%mw(itr) = mw_cbmz(n)
        end if
      end do

      if ( myid == italk ) then
        if ( itr == 1 )  write(stdout,*) 'tracer', ' cbmz index',' molw'
        write(stdout,*) chtrname(itr),  trac%indcbmz(itr),  trac%mw(itr)
      end if

      !!$ Define also some specific indices for practical purpose
      !  Not all the possible cbmz species have a tracer index:
      !  however this information is also contained in trac%indcbmz table
      !
      !CBMZ mechanims
      if ( chtrname(itr) == 'NO'     ) ino         = itr
      if ( chtrname(itr) == 'NO2'    ) ino2        = itr
      if ( chtrname(itr) == 'N2O5'   ) in2o5       = itr
      if ( chtrname(itr) == 'HNO2'   ) ihno2       = itr
      if ( chtrname(itr) == 'HNO3'   ) ihno3       = itr
      if ( chtrname(itr) == 'HNO4'   ) ihno4       = itr
      if ( chtrname(itr) == 'O3'     ) io3         = itr
      if ( chtrname(itr) == 'H2O2'   ) ih2o2       = itr
      if ( chtrname(itr) == 'CO'     ) ico         = itr
      if ( chtrname(itr) == 'SO2'    ) iso2        = itr
      if ( chtrname(itr) == 'DMS'    ) idms        = itr
      if ( chtrname(itr) == 'H2SO4'  ) ih2so4      = itr
      if ( chtrname(itr) == 'CH4'    ) ich4        = itr
      if ( chtrname(itr) == 'C2H6'   ) ic2h6       = itr
      if ( chtrname(itr) == 'PAR'    ) ipar        = itr
      if ( chtrname(itr) == 'CH3OH'  ) ich3oh      = itr
      if ( chtrname(itr) == 'HCHO'   ) ihcho       = itr
      if ( chtrname(itr) == 'ALD2'   ) iald2       = itr
      if ( chtrname(itr) == 'AONE'   ) iaone       = itr
      if ( chtrname(itr) == 'ETH'    ) ieth        = itr
      if ( chtrname(itr) == 'OLET'   ) iolet       = itr
      if ( chtrname(itr) == 'OLEI'   ) iolei       = itr
      if ( chtrname(itr) == 'TOL'    ) itol        = itr
      if ( chtrname(itr) == 'XYL'    ) ixyl        = itr
      if ( chtrname(itr) == 'ISOP'   ) iisop       = itr
      if ( chtrname(itr) == 'ONIT'   ) ionit       = itr
      if ( chtrname(itr) == 'PAN'    ) ipan        = itr
      if ( chtrname(itr) == 'HCOOH'  ) ihcooh      = itr
      if ( chtrname(itr) == 'RCOOH'  ) ircooh      = itr
      if ( chtrname(itr) == 'CH3OOH' ) ich3ooh     = itr
      if ( chtrname(itr) == 'ETHOOH' ) iethooh     = itr
      if ( chtrname(itr) == 'ROOH'   ) irooh       = itr
      if ( chtrname(itr) == 'MGLY'   ) imgly       = itr
      if ( chtrname(itr) == 'ISOPRD' ) iisoprd     = itr
      if ( chtrname(itr) == 'ISOPN'  ) iisopn      = itr
      if ( chtrname(itr) == 'OPEN'   ) iopen       = itr
      if ( chtrname(itr) == 'CRES'   ) icres       = itr
      !!$
      end if

      if ( chemsimtype(1:5) ==  'CB6C') then
      ! gas phas species (CB6r2 = CBSIX)

      if ( isorgam == 1 ) then
        do n = 1,3*ngstr+nsoatr
          if ( chtrname(itr) == srgspec(n))then
            trac%indcb6(itr) = n
            trac%mw(itr) = mw_srg(n)
          end if
        end do

        if ( myid == italk ) then
          if ( itr == 1 )  write(stdout,*)'associating cb6 index and molw values'
        end if

      else
        do n = 1,totch
          if ( chtrname(itr) == cb6spec(n))then
            trac%indcb6(itr) = n  
            ! index of the tracer in the CB6r2 list of species
            trac%mw(itr) = mw_cb6(n) 
            ! correponding molecular weight
          end if
        end do

        if ( myid == italk ) then
          if ( itr == 1 )  write(stdout,*) 'tracer', ' cb6 index',' molw'
          write(stdout,*) chtrname(itr),  trac%indcb6(itr),  trac%mw(itr)
        end if

      end if

      !!$ Define also some specific indices for practical purpose
      !  Not all the possible cb6 species have a tracer index:
      !  however this information is also contained in trac%indcb6 table
      !
      !CBSIX mechanims
      if ( chtrname(itr) == 'CDOX'    ) icdox       = itr
      if ( chtrname(itr) == 'SULF'    ) isulf       = itr
      if ( chtrname(itr) == 'SDIO'    ) isdio       = itr
      if ( chtrname(itr) == 'ETHA'    ) ietha       = itr
      if ( chtrname(itr) == 'ETHY'    ) iethy       = itr
      if ( chtrname(itr) == 'DNPO'    ) idnpo       = itr
      if ( chtrname(itr) == 'BENZ'    ) ibenz       = itr
      if ( chtrname(itr) == 'EPOX'    ) iepox       = itr
      if ( chtrname(itr) == 'ETOH'    ) ietoh       = itr
      if ( chtrname(itr) == 'PRPA'    ) iprpa       = itr
      if ( chtrname(itr) == 'KET'     ) iket        = itr
      if ( chtrname(itr) == 'TOLN'    ) itoln       = itr
      if ( chtrname(itr) == 'XYLN'    ) ixyln       = itr
      if ( chtrname(itr) == 'HPLD'    ) ihpld       = itr
      if ( chtrname(itr) == 'PACN'    ) ipacn       = itr
      if ( chtrname(itr) == 'PACD'    ) ipacd       = itr
      if ( chtrname(itr) == 'NTR2'    ) intr2       = itr
      if ( chtrname(itr) == 'PNA'     ) ipna        = itr
      if ( chtrname(itr) == 'MEOH'    ) imeoh       = itr
      if ( chtrname(itr) == 'HONO'    ) ihono       = itr
      if ( chtrname(itr) == 'MEPX'    ) imepx       = itr
      if ( chtrname(itr) == 'OPAN'    ) iopan       = itr
      if ( chtrname(itr) == 'CAT1'    ) icat1       = itr
      if ( chtrname(itr) == 'HPOX'    ) ihpox       = itr
      if ( chtrname(itr) == 'ISPX'    ) iispx       = itr
      if ( chtrname(itr) == 'FACD'    ) ifacd       = itr
      if ( chtrname(itr) == 'RPOX'    ) irpox       = itr
      if ( chtrname(itr) == 'NTR1'    ) intr1       = itr
      if ( chtrname(itr) == 'ACET'    ) iacet       = itr
      if ( chtrname(itr) == 'INTR'    ) iintr       = itr
      if ( chtrname(itr) == 'CRON'    ) icron       = itr
      if ( chtrname(itr) == 'AACD'    ) iaacd       = itr
      if ( chtrname(itr) == 'ETHE'    ) iethe       = itr
      if ( chtrname(itr) == 'CMON'    ) icmon       = itr
      if ( chtrname(itr) == 'TERP'    ) iterp       = itr
      if ( chtrname(itr) == 'CRSL'    ) icrsl       = itr
      if ( chtrname(itr) == 'ISPR'    ) iispr       = itr
      if ( chtrname(itr) == 'NTRC'    ) intrc       = itr
      if ( chtrname(itr) == 'GLYD'    ) iglyd       = itr
      if ( chtrname(itr) == 'GLY'     ) igly        = itr
      if ( chtrname(itr) == 'XOPN'    ) ixopn       = itr
      if ( chtrname(itr) == 'MEGY'    ) imegy       = itr
      if ( chtrname(itr) == 'ROPN'    ) iropn       = itr
      if ( chtrname(itr) == 'IOLE'    ) iiole       = itr
      if ( chtrname(itr) == 'FORM'    ) iform       = itr
      if ( chtrname(itr) == 'OLE'     ) iole        = itr
      if ( chtrname(itr) == 'AALD'    ) iaald       = itr
      if ( chtrname(itr) == 'ISPD'    ) iispd       = itr
      if ( chtrname(itr) == 'ALDX'    ) ialdx       = itr
      if ( chtrname(itr) == 'ALKA'    ) ialka       = itr
      if ( chtrname(itr) == 'OZN'     ) iozn        = itr
      if ( chtrname(itr) == 'NMOX'    ) inmox       = itr
      if ( chtrname(itr) == 'NDOX'    ) indox       = itr
      if ( chtrname(itr) == 'NTOX'    ) intox       = itr
      if ( chtrname(itr) == 'PANX'    ) ipanx       = itr
      ! transporting all species !
!      if ( chtrname(itr) == 'OSNG'    ) iosng       = itr
      if ( chtrname(itr) == 'MTHA'    ) imtha       = itr
!      if ( chtrname(itr) == 'HCO3'    ) ihco3       = itr
!      if ( chtrname(itr) == 'CRER'    ) icrer       = itr
!      if ( chtrname(itr) == 'BZO2'    ) ibzo2       = itr
!      if ( chtrname(itr) == 'ROR'     ) iror        = itr
!      if ( chtrname(itr) == 'TOLR'    ) itolr       = itr
!      if ( chtrname(itr) == 'XO2R'    ) ixo2r       = itr
!      if ( chtrname(itr) == 'EPX2'    ) iepx2       = itr
!      if ( chtrname(itr) == 'XYLR'    ) ixylr       = itr
!      if ( chtrname(itr) == 'OPO3'    ) iopo3       = itr
!      if ( chtrname(itr) == 'XO2N'    ) ixo2n       = itr
!      if ( chtrname(itr) == 'ISO2'    ) iiso2       = itr
!      if ( chtrname(itr) == 'XO2H'    ) ixo2h       = itr
!      if ( chtrname(itr) == 'MEO2'    ) imeo2       = itr
!      if ( chtrname(itr) == 'ROO'     ) iroo        = itr
!      if ( chtrname(itr) == 'HOX'     ) ihox        = itr
!      if ( chtrname(itr) == 'POX'     ) ipox        = itr
!      if ( chtrname(itr) == 'CXO3'    ) icxo3       = itr
!      if ( chtrname(itr) == 'ACOO'    ) iacoo       = itr
!      if ( chtrname(itr) == 'O'       ) ioxy        = itr
      ! nucleation aerosols
      if ( chtrname(itr) == 'CDOXnu'  ) incdox      = itr
      if ( chtrname(itr) == 'SULFnu'  ) insulf      = itr
      if ( chtrname(itr) == 'SDIOnu'  ) insdio      = itr
      if ( chtrname(itr) == 'ETHAnu'  ) inetha      = itr
      if ( chtrname(itr) == 'ETHYnu'  ) inethy      = itr
      if ( chtrname(itr) == 'DNPOnu'  ) indnpo      = itr
      if ( chtrname(itr) == 'BENZnu'  ) inbenz      = itr
      if ( chtrname(itr) == 'EPOXnu'  ) inepox      = itr
      if ( chtrname(itr) == 'ETOHnu'  ) inetoh      = itr
      if ( chtrname(itr) == 'PRPAnu'  ) inprpa      = itr
      if ( chtrname(itr) == 'KETnu'   ) inket       = itr
      if ( chtrname(itr) == 'TOLNnu'  ) intoln      = itr
      if ( chtrname(itr) == 'XYLNnu'  ) inxyln      = itr
      if ( chtrname(itr) == 'HPLDnu'  ) inhpld      = itr
      if ( chtrname(itr) == 'PACNnu'  ) inpacn      = itr
      if ( chtrname(itr) == 'PACDnu'  ) inpacd      = itr
      if ( chtrname(itr) == 'NTR2nu'  ) inntr2      = itr
      if ( chtrname(itr) == 'PNAnu'   ) inpna       = itr
      if ( chtrname(itr) == 'MEOHnu'  ) inmeoh      = itr
      if ( chtrname(itr) == 'HONOnu'  ) inhono      = itr
      if ( chtrname(itr) == 'MEPXnu'  ) inmepx      = itr
      if ( chtrname(itr) == 'OPANnu'  ) inopan      = itr
      if ( chtrname(itr) == 'CAT1nu'  ) incat1      = itr
      if ( chtrname(itr) == 'HPOXnu'  ) inhpox      = itr
      if ( chtrname(itr) == 'ISPXnu'  ) inispx      = itr
      if ( chtrname(itr) == 'FACDnu'  ) infacd      = itr
      if ( chtrname(itr) == 'RPOXnu'  ) inrpox      = itr
      if ( chtrname(itr) == 'NTR1nu'  ) inntr1      = itr
      if ( chtrname(itr) == 'ACETnu'  ) inacet      = itr
      if ( chtrname(itr) == 'INTRnu'  ) inintr      = itr
      if ( chtrname(itr) == 'CRONnu'  ) incron      = itr
      if ( chtrname(itr) == 'AACDnu'  ) inaacd      = itr
      if ( chtrname(itr) == 'ETHEnu'  ) inethe      = itr
      if ( chtrname(itr) == 'CMONnu'  ) incmon      = itr
      if ( chtrname(itr) == 'TERPnu'  ) interp      = itr
      if ( chtrname(itr) == 'CRSLnu'  ) incrsl      = itr
      if ( chtrname(itr) == 'ISPRnu'  ) inispr      = itr
      if ( chtrname(itr) == 'NTRCnu'  ) inntrc      = itr
      if ( chtrname(itr) == 'GLYDnu'  ) inglyd      = itr
      if ( chtrname(itr) == 'GLYnu'   ) ingly       = itr
      if ( chtrname(itr) == 'XOPNnu'  ) inxopn      = itr
      if ( chtrname(itr) == 'MEGYnu'  ) inmegy      = itr
      if ( chtrname(itr) == 'ROPNnu'  ) inropn      = itr
      if ( chtrname(itr) == 'IOLEnu'  ) iniole      = itr
      if ( chtrname(itr) == 'FORMnu'  ) inform      = itr
      if ( chtrname(itr) == 'OLEnu'   ) inole       = itr
      if ( chtrname(itr) == 'AALDnu'  ) inaald      = itr
      if ( chtrname(itr) == 'ISPDnu'  ) inispd      = itr
      if ( chtrname(itr) == 'ALDXnu'  ) inaldx      = itr
      if ( chtrname(itr) == 'ALKAnu'  ) inalka      = itr
      if ( chtrname(itr) == 'OZNnu'   ) inozn       = itr
      if ( chtrname(itr) == 'NMOXnu'  ) innmox      = itr
      if ( chtrname(itr) == 'NDOXnu'  ) inndox      = itr
      if ( chtrname(itr) == 'NTOXnu'  ) inntox      = itr
      if ( chtrname(itr) == 'PANXnu'  ) inpanx      = itr
      if ( chtrname(itr) == 'OSNGnu'  ) inosng      = itr
      if ( chtrname(itr) == 'MTHAnu'  ) inmtha      = itr
      if ( chtrname(itr) == 'HCO3nu'  ) inhco3      = itr
      if ( chtrname(itr) == 'CRERnu'  ) increr      = itr
      if ( chtrname(itr) == 'BZO2nu'  ) inbzo2      = itr
      if ( chtrname(itr) == 'RORnu'   ) inror       = itr
      if ( chtrname(itr) == 'TOLRnu'  ) intolr      = itr
      if ( chtrname(itr) == 'XO2Rnu'  ) inxo2r      = itr
      if ( chtrname(itr) == 'EPX2nu'  ) inepx2      = itr
      if ( chtrname(itr) == 'XYLRnu'  ) inxylr      = itr
      if ( chtrname(itr) == 'OPO3nu'  ) inopo3      = itr
      if ( chtrname(itr) == 'XO2Nnu'  ) inxo2n      = itr
      if ( chtrname(itr) == 'ISO2nu'  ) iniso2      = itr
      if ( chtrname(itr) == 'XO2Hnu'  ) inxo2h      = itr
      if ( chtrname(itr) == 'MEO2nu'  ) inmeo2      = itr
      if ( chtrname(itr) == 'ROOnu'   ) inroo       = itr
      if ( chtrname(itr) == 'HOXnu'   ) inhox       = itr
      if ( chtrname(itr) == 'POXnu'   ) inpox       = itr
      if ( chtrname(itr) == 'CXO3nu'  ) incxo3      = itr
      if ( chtrname(itr) == 'ACOOnu'  ) inacoo      = itr
      if ( chtrname(itr) == 'Onu'     ) inoxy       = itr
      ! accumulation aerosols
      if ( chtrname(itr) == 'CDOXac'  ) iacdox      = itr
      if ( chtrname(itr) == 'SULFac'  ) iasulf      = itr
      if ( chtrname(itr) == 'SDIOac'  ) iasdio      = itr
      if ( chtrname(itr) == 'ETHAac'  ) iaetha      = itr
      if ( chtrname(itr) == 'ETHYac'  ) iaethy      = itr
      if ( chtrname(itr) == 'DNPOac'  ) iadnpo      = itr
      if ( chtrname(itr) == 'BENZac'  ) iabenz      = itr
      if ( chtrname(itr) == 'EPOXac'  ) iaepox      = itr
      if ( chtrname(itr) == 'ETOHac'  ) iaetoh      = itr
      if ( chtrname(itr) == 'PRPAac'  ) iaprpa      = itr
      if ( chtrname(itr) == 'KETac'   ) iaket       = itr
      if ( chtrname(itr) == 'TOLNac'  ) iatoln      = itr
      if ( chtrname(itr) == 'XYLNac'  ) iaxyln      = itr
      if ( chtrname(itr) == 'HPLDac'  ) iahpld      = itr
      if ( chtrname(itr) == 'PACNac'  ) iapacn      = itr
      if ( chtrname(itr) == 'PACDac'  ) iapacd      = itr
      if ( chtrname(itr) == 'NTR2ac'  ) iantr2      = itr
      if ( chtrname(itr) == 'PNAac'   ) iapna       = itr
      if ( chtrname(itr) == 'MEOHac'  ) iameoh      = itr
      if ( chtrname(itr) == 'HONOac'  ) iahono      = itr
      if ( chtrname(itr) == 'MEPXac'  ) iamepx      = itr
      if ( chtrname(itr) == 'OPANac'  ) iaopan      = itr
      if ( chtrname(itr) == 'CAT1ac'  ) iacat1      = itr
      if ( chtrname(itr) == 'HPOXac'  ) iahpox      = itr
      if ( chtrname(itr) == 'ISPXac'  ) iaispx      = itr
      if ( chtrname(itr) == 'FACDac'  ) iafacd      = itr
      if ( chtrname(itr) == 'RPOXac'  ) iarpox      = itr
      if ( chtrname(itr) == 'NTR1ac'  ) iantr1      = itr
      if ( chtrname(itr) == 'ACETac'  ) iaacet      = itr
      if ( chtrname(itr) == 'INTRac'  ) iaintr      = itr
      if ( chtrname(itr) == 'CRONac'  ) iacron      = itr
      if ( chtrname(itr) == 'AACDac'  ) iaaacd      = itr
      if ( chtrname(itr) == 'ETHEac'  ) iaethe      = itr
      if ( chtrname(itr) == 'CMONac'  ) iacmon      = itr
      if ( chtrname(itr) == 'TERPac'  ) iaterp      = itr
      if ( chtrname(itr) == 'CRSLac'  ) iacrsl      = itr
      if ( chtrname(itr) == 'ISPRac'  ) iaispr      = itr
      if ( chtrname(itr) == 'NTRCac'  ) iantrc      = itr
      if ( chtrname(itr) == 'GLYDac'  ) iaglyd      = itr
      if ( chtrname(itr) == 'GLYac'   ) iagly       = itr
      if ( chtrname(itr) == 'XOPNac'  ) iaxopn      = itr
      if ( chtrname(itr) == 'MEGYac'  ) iamegy      = itr
      if ( chtrname(itr) == 'ROPNac'  ) iaropn      = itr
      if ( chtrname(itr) == 'IOLEac'  ) iaiole      = itr
      if ( chtrname(itr) == 'FORMac'  ) iaform      = itr
      if ( chtrname(itr) == 'OLEac'   ) iaole       = itr
      if ( chtrname(itr) == 'AALDac'  ) iaaald      = itr
      if ( chtrname(itr) == 'ISPDac'  ) iaispd      = itr
      if ( chtrname(itr) == 'ALDXac'  ) iaaldx      = itr
      if ( chtrname(itr) == 'ALKAac'  ) iaalka      = itr
      if ( chtrname(itr) == 'OZNac'   ) iaozn       = itr
      if ( chtrname(itr) == 'NMOXac'  ) ianmox      = itr
      if ( chtrname(itr) == 'NDOXac'  ) iandox      = itr
      if ( chtrname(itr) == 'NTOXac'  ) iantox      = itr
      if ( chtrname(itr) == 'PANXac'  ) iapanx      = itr
      if ( chtrname(itr) == 'OSNGac'  ) iaosng      = itr
      if ( chtrname(itr) == 'MTHAac'  ) iamtha      = itr
      if ( chtrname(itr) == 'HCO3ac'  ) iahco3      = itr
      if ( chtrname(itr) == 'CRERac'  ) iacrer      = itr
      if ( chtrname(itr) == 'BZO2ac'  ) iabzo2      = itr
      if ( chtrname(itr) == 'RORac'   ) iaror       = itr
      if ( chtrname(itr) == 'TOLRac'  ) iatolr      = itr
      if ( chtrname(itr) == 'XO2Rac'  ) iaxo2r      = itr
      if ( chtrname(itr) == 'EPX2ac'  ) iaepx2      = itr
      if ( chtrname(itr) == 'XYLRac'  ) iaxylr      = itr
      if ( chtrname(itr) == 'OPO3ac'  ) iaopo3      = itr
      if ( chtrname(itr) == 'XO2Nac'  ) iaxo2n      = itr
      if ( chtrname(itr) == 'ISO2ac'  ) iaiso2      = itr
      if ( chtrname(itr) == 'XO2Hac'  ) iaxo2h      = itr
      if ( chtrname(itr) == 'MEO2ac'  ) iameo2      = itr
      if ( chtrname(itr) == 'ROOac'   ) iaroo       = itr
      if ( chtrname(itr) == 'HOXac'   ) iahox       = itr
      if ( chtrname(itr) == 'POXac'   ) iapox       = itr
      if ( chtrname(itr) == 'CXO3ac'  ) iacxo3      = itr
      if ( chtrname(itr) == 'ACOOac'  ) iaacoo      = itr
      if ( chtrname(itr) == 'Oac'     ) iaoxy       = itr
      ! sorgam aggregates
      if ( chtrname(itr) == 'CNYLnu'  ) icnylnu     = itr
      if ( chtrname(itr) == 'STHCnu'  ) isthcnu     = itr
      if ( chtrname(itr) == 'OLEFnu'  ) iolefnu     = itr
      if ( chtrname(itr) == 'AROMnu'  ) iaromnu     = itr
      if ( chtrname(itr) == 'NTROnu'  ) intronu     = itr
      if ( chtrname(itr) == 'MFNCnu'  ) imfncnu     = itr
      if ( chtrname(itr) == 'CNYLac'  ) icnylac     = itr
      if ( chtrname(itr) == 'STHCac'  ) isthcac     = itr
      if ( chtrname(itr) == 'OLEFac'  ) iolefac     = itr
      if ( chtrname(itr) == 'AROMac'  ) iaromac     = itr
      if ( chtrname(itr) == 'NTROac'  ) introac     = itr
      if ( chtrname(itr) == 'MFNCac'  ) imfncac     = itr

      end if


      if ( chtrname(itr) == 'NH3'   ) inh3       = itr

      if ( chtrname(itr) == 'POLLEN') ipollen   = itr

! special case of biogenic options
#if defined CLM45 || (defined VOC && defined CLM)
      !abt *** Added below to determine which MEGAN biogenic emission species
      !    *** will be passed to the gas phase mechanism
      !    *** commented out lines correspond to species not advected but
      !    *** potentially used in chemistry mechanism.
      !    *** Uncomment to give potential to advect
      if ( igaschem == 1 ) then
        if ( chtrname(itr) == 'ISOP'  ) bvoc_trmask(itr) = 1
        if ( chtrname(itr) == 'APIN'  ) bvoc_trmask(itr) = 1
        if ( chtrname(itr) == 'LIMO'  ) bvoc_trmask(itr) = 1
      end if
#endif

    end do

    if ( isorgam == 1 ) then
      trac%soa_category(:) = -1
      chtrsol(:) = soloc ! all tracers given soloc for wetdep
      do n = 1,ngstr
        trac%soa_category(n) = soa_cat_cb6(trac%indcb6(n))
        ! argument for soa from cb6 gas
      end do
      trac%soa_category(ngstr+1:2*ngstr)   = trac%soa_category(1:ngstr)
      trac%soa_category(2*ngstr+1:3*ngstr) = trac%soa_category(1:ngstr)
      do itr = 1 , ntr
        if ( myid == italk ) then
          if ( itr == 1         ) write(stdout,*) 'tracer', ' cb6 index',' molw'
          if ( itr < ngstr+1    ) write(stdout,*) & 
             chtrname(itr),  trac%indcb6(itr),  trac%mw(itr)
          if ( itr == ngstr+1   ) write(stdout,*) 'nu-aer', ' cb6 index',' molw'
          if ( itr > ngstr   .and. itr < 2*ngstr+1 ) write(stdout,*) & 
             chtrname(itr),  trac%indcb6(itr),  trac%mw(itr)
          if ( itr == 2*ngstr+1 ) write(stdout,*) 'ac-aer', ' cb6 index',' molw'
          if ( itr > 2*ngstr .and. itr < 3*ngstr+1 ) write(stdout,*) & 
             chtrname(itr),  trac%indcb6(itr),  trac%mw(itr)
          if ( itr == 3*ngstr+1 ) write(stdout,*) 'tr-soa', ' cb6 index: sorgam'
          if ( itr > 3*ngstr    ) write(stdout,*) & 
             chtrname(itr),  trac%indcb6(itr)
        end if
      end do
    end if
    !
    ! define now correspndance between boundary species indices and
    ! determine tracer indices corresponding to ch boundary conditions
    !
    ichbdy2trac(:) = 0
    itr = 1
    if( igaschem == 1 ) then
      if ( chemsimtype(1:5) ==  'CB6C') then
        do n = 1 , n_chbcvar6
          do i = 1 , ntr
            if (chbcname6(n) == chtrname(i)) then
              ichbdy2trac(itr) = i
              itr = itr + 1
            end if
          end do
        end do
      else
        do n = 1 , n_chbcvar
          do i = 1 , ntr
            if (chbcname(n) == chtrname(i)) then
              ichbdy2trac(itr) = i
              itr = itr + 1
            end if
          end do
        end do
      end if
    end if
    !
    ! look also in aerosol bc and pile them after.
    !
    if ( iaerosol == 1 .and. isorgam == 0 ) then
      do n = 1 , size(aeaero)
        do i = 1 , ntr
          if ( aeaero(n) == chtrname(i) ) then
            ichbdy2trac(itr) = i
            itr = itr + 1
          end if
        end do
      end do
    end if
    if ( myid == italk ) then
      write(stdout,*) 'tracer index coreesponding to bdy species '
      call iprntv(ichbdy2trac,size(ichbdy2trac),'ichbdy2trac')
    end if
!!$  FAB : work on that later
!!$    do itr = 1,ntr
!!$       do n = 1,n_chbcvar
!!$       if ( chtrname(itr) == chbcname(n))then
!!$        trac%indchbdy(itr) = n  ! index of the tracer in the chbdy list
!!$       end if
!!$      end do
!!$        print*,'test', itr, chtrname(itr), trac%indchbdy(itr)
!!$     end do
    if ( idust(1) > 0 .or. ichbion==1) then
      ! activate dust initialization
      if ( myid == italk ) write(stdout,*) 'Calling inidust'
      call inidust
    end if

    if ( ichbion == 1 ) call ini_bionit

    if ( igaschem == 1 ) then
      !if ( chemsimtype(1:5) ==  'CB6C') then
      !  open(26,file='TUVGRID3', status='old', err=903)
      !else
        open(26,file='TUVGRID2', status='old', err=900)
      !end if
      ! not used in KPP
!       open(25,file='REACTION.DAT_CBMZ', status='old', err=901)
!       open(27,file='cbmz_chemmech.out', status='replace', err=902)
! 902   continue
!       call chemread
      call hvread
!       call cheminit
    end if

    call init_mod_che_ncio(chemsimtype)
    call che_init_bdy

    ! Finally initialise chia and chib to chib0 over the whole domain

    if ( .not. ifrest ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            chia(j,i,k,:) = chib0(j,i,k,:)
            chib(j,i,k,:) = chib0(j,i,k,:)
          end do
        end do
      end do
    end if

    return

900 write(stderr,*) 'Cannot open required file TUVGRID2.'
    call fatal(__FILE__,__LINE__,'TUVGRID2 NOT FOUND')
!901 write(stderr,*) 'Cannot open required file REACTION.DAT_CBMZ.'
!    call fatal(__FILE__,__LINE__,'REACTION.DAT_CBMZ NOT FOUND')
903 write(stderr,*) 'Cannot open required file TUVGRID3.'
    call fatal(__FILE__,__LINE__,'TUVGRID3 NOT FOUND')

  end subroutine start_chem

end module mod_che_start
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
