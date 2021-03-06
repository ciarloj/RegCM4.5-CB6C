! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Utility Data Module File
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
! File                 : mod_cb6_Monitor.f90
! Time                 : Thu Aug 31 17:04:30 2017
! Working directory    : /home/regcm/Code/kpp-2.2.3/Workspace/GAS_CB6r2_mod
! Equation file        : cb6.kpp
! Output root filename : cb6
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mod_cb6_Monitor


  CHARACTER(LEN=15), PARAMETER, DIMENSION(81) :: SPC_NAMES = (/ &
     'SULF           ','CDOX           ','SDIO           ', &
     'OSNG           ','ETHY           ','DNPO           ', &
     'BENZ           ','EPOX           ','ETOH           ', &
     'KET            ','TOLN           ','XYLN           ', &
     'HPLD           ','PACN           ','PACD           ', &
     'PNA            ','HONO           ','MTHA           ', &
     'NTR2           ','MEPX           ','HPOX           ', &
     'OPAN           ','CAT1           ','ISPX           ', &
     'ETHA           ','FACD           ','PANX           ', &
     'HCO3           ','PRPA           ','MEOH           ', &
     'CRER           ','RPOX           ','NTR1           ', &
     'ACET           ','INTR           ','AACD           ', &
     'CRON           ','ROR            ','BZO2           ', &
     'ETHE           ','CMON           ','TOLR           ', &
     'XYLR           ','CRSL           ','TERP           ', &
     'ISPR           ','NTRC           ','EPX2           ', &
     'GLYD           ','GLY            ','XOPN           ', &
     'MEGY           ','ROPN           ','IOLE           ', &
     'OLE            ','FORM           ','AALD           ', &
     'XO2R           ','OPO3           ','XO2N           ', &
     'ISO2           ','XO2H           ','ISPD           ', &
     'MEO2           ','ALKA           ','ALDX           ', &
     'HOX            ','ROO            ','NTOX           ', &
     'O              ','POX            ','CXO3           ', &
     'ACOO           ','OZN            ','NMOX           ', &
     'NDOX           ','WTR            ','DIHY           ', &
     'O2             ','M              ','DUMMY2         ' /)

  INTEGER, PARAMETER, DIMENSION(1) :: LOOKAT = (/ &
      74 /)

  INTEGER, PARAMETER, DIMENSION(1) :: MONITOR = (/ &
      74 /)

  CHARACTER(LEN=15), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_0 = (/ &
     '             NDOX --> O + NMOX                                                                      ', &
     '       O + O2 + M --> OZN + M                                                                       ', &
     '       OZN + NMOX --> NDOX                                                                          ', &
     '     O + NMOX + M --> NDOX + M                                                                      ', &
     '         O + NDOX --> NMOX                                                                          ', &
     '         O + NDOX --> NTOX                                                                          ', &
     '          O + OZN --> DUMMY2                                                                        ', &
     '              OZN --> O                                                                             ', &
     '              OZN --> OSNG                                                                          ', &
     '         OSNG + M --> O + M                                                                         ', &
     '       OSNG + WTR --> 2 HOX                                                                         ', &
     '        HOX + OZN --> POX                                                                           ', &
     '        POX + OZN --> HOX                                                                           ', &
     '          HOX + O --> POX                                                                           ', &
     '          O + POX --> HOX                                                                           ', &
     '            2 HOX --> O                                                                             ', &
     '            2 HOX --> HPOX                                                                          ', &
     '        HOX + POX --> DUMMY2                                                                        ', &
     '            2 POX --> HPOX                                                                          ', &
     '      2 POX + WTR --> HPOX                                                                          ', &
     '             HPOX --> 2 HOX                                                                         ', &
     '       HPOX + HOX --> POX                                                                           ', &
     '         HPOX + O --> HOX + POX                                                                     ', &
     '      2 NMOX + O2 --> 2 NDOX                                                                        ', &
     '       POX + NMOX --> HOX + NDOX                                                                    ', &
     '       OZN + NDOX --> NTOX                                                                          ', &
     '             NTOX --> O + NDOX                                                                      ', &
     '             NTOX --> NMOX                                                                          ', &
     '      NTOX + NMOX --> 2 NDOX                                                                        ', &
     '      NTOX + NDOX --> NMOX + NDOX                                                                   ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_1 = (/ &
     '         NTOX + O --> NDOX                                                                          ', &
     '       HOX + NTOX --> POX + NDOX                                                                    ', &
     '       NTOX + POX --> HOX + NDOX                                                                    ', &
     '       NTOX + OZN --> NDOX                                                                          ', &
     '           2 NTOX --> 2 NDOX                                                                        ', &
     '      NTOX + NDOX --> DNPO                                                                          ', &
     '             DNPO --> NTOX + NDOX                                                                   ', &
     '             DNPO --> NTOX + NDOX                                                                   ', &
     '       DNPO + WTR --> 2 NTRC                                                                        ', &
     '       HOX + NMOX --> HONO                                                                          ', &
     'NMOX + NDOX + WTR --> 2 HONO                                                                        ', &
     '           2 HONO --> NMOX + NDOX                                                                   ', &
     '             HONO --> HOX + NMOX                                                                    ', &
     '       HONO + HOX --> NDOX                                                                          ', &
     '       HOX + NDOX --> NTRC                                                                          ', &
     '       NTRC + HOX --> NTOX                                                                          ', &
     '             NTRC --> HOX + NDOX                                                                    ', &
     '       POX + NDOX --> PNA                                                                           ', &
     '              PNA --> POX + NDOX                                                                    ', &
     '              PNA --> 0.41 HOX + 0.41 NTOX + 0.59 POX + 0.59 NDOX                                   ', &
     '        PNA + HOX --> NDOX                                                                          ', &
     '       SDIO + HOX --> SULF + POX                                                                    ', &
     '      ACOO + NMOX --> MEO2 + ROO + NDOX                                                             ', &
     '      ACOO + NDOX --> PACN                                                                          ', &
     '             PACN --> ACOO + NDOX                                                                   ', &
     '             PACN --> 0.4 MEO2 + 0.4 ROO + 0.4 NTOX + 0.6 ACOO + 0.6 NDOX ... etc.                  ', &
     '       POX + ACOO --> 0.41 PACD + 0.15 AACD + 0.44 MEO2 + 0.44 HOX + 0.44 ROO ... etc.              ', &
     '       ROO + ACOO --> ACOO                                                                          ', &
     '           2 ACOO --> 2 MEO2 + 2 ROO                                                                ', &
     '      CXO3 + ACOO --> AALD + XO2H + MEO2 + 2 ROO                                                    ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_2 = (/ &
     '      CXO3 + NMOX --> AALD + XO2H + ROO + NDOX                                                      ', &
     '      CXO3 + NDOX --> PANX                                                                          ', &
     '             PANX --> CXO3 + NDOX                                                                   ', &
     '             PANX --> 0.4 AALD + 0.4 XO2H + 0.4 ROO + 0.4 NTOX + 0.6 CXO3 ... etc.                  ', &
     '       POX + CXO3 --> 0.41 PACD + 0.15 AACD + 0.44 AALD + 0.44 XO2H + 0.44 HOX ... etc.             ', &
     '       ROO + CXO3 --> 0.8 AALD + 0.8 XO2H + 0.8 ROO                                                 ', &
     '           2 CXO3 --> 2 AALD + 2 XO2H + 2 ROO                                                       ', &
     '       ROO + NMOX --> NMOX                                                                          ', &
     '        ROO + POX --> POX                                                                           ', &
     '            2 ROO --> DUMMY2                                                                        ', &
     '      MEO2 + NMOX --> FORM + POX + NDOX                                                             ', &
     '       MEO2 + POX --> 0.9 MEPX + 0.1 FORM                                                           ', &
     '      MEO2 + ACOO --> 0.1 AACD + FORM + 0.9 MEO2 + 0.9 ROO + 0.9 POX                                ', &
     '       MEO2 + ROO --> 0.315 MEOH + 0.685 FORM + ROO + 0.37 POX                                      ', &
     '      XO2H + NMOX --> POX + NDOX                                                                    ', &
     '       XO2H + POX --> RPOX                                                                          ', &
     '      XO2H + ACOO --> 0.2 AACD + 0.8 MEO2 + 0.8 ROO + 0.8 POX                                       ', &
     '       XO2H + ROO --> ROO + 0.6 POX                                                                 ', &
     '      XO2R + NMOX --> NDOX                                                                          ', &
     '       XO2R + POX --> RPOX                                                                          ', &
     '      XO2R + ACOO --> 0.2 AACD + 0.8 MEO2 + 0.8 ROO                                                 ', &
     '       XO2R + ROO --> ROO                                                                           ', &
     '      XO2N + NMOX --> 0.5 NTR2 + 0.5 NTR1                                                           ', &
     '       XO2N + POX --> RPOX                                                                          ', &
     '      XO2N + ACOO --> 0.2 AACD + 0.8 MEO2 + 0.8 ROO + 0.8 POX                                       ', &
     '       XO2N + ROO --> ROO                                                                           ', &
     '       MEPX + HOX --> 0.4 FORM + 0.6 MEO2 + 0.4 HOX + 0.6 ROO                                       ', &
     '             MEPX --> MEO2 + HOX + ROO                                                              ', &
     '       RPOX + HOX --> 0.06 XO2N + 0.54 XO2H + 0.4 HOX + 0.6 ROO                                     ', &
     '             RPOX --> HOX + POX                                                                     ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_3 = (/ &
     '       NTR1 + HOX --> NTR2                                                                          ', &
     '             NTR1 --> NDOX                                                                          ', &
     '       FACD + HOX --> POX                                                                           ', &
     '       AACD + HOX --> MEO2 + ROO                                                                    ', &
     '       PACD + HOX --> ACOO                                                                          ', &
     '       FORM + HOX --> CMON + POX                                                                    ', &
     '             FORM --> CMON + 2 POX                                                                  ', &
     '             FORM --> CMON + DIHY                                                                   ', &
     '         FORM + O --> CMON + HOX + POX                                                              ', &
     '      FORM + NTOX --> CMON + NTRC + POX                                                             ', &
     '       FORM + POX --> HCO3                                                                          ', &
     '             HCO3 --> FORM + POX                                                                    ', &
     '      HCO3 + NMOX --> FACD + POX + NDOX                                                             ', &
     '       HCO3 + POX --> 0.5 MEPX + 0.5 FACD + 0.2 HOX + 0.2 POX                                       ', &
     '         AALD + O --> HOX + ACOO                                                                    ', &
     '       AALD + HOX --> ACOO                                                                          ', &
     '      AALD + NTOX --> NTRC + ACOO                                                                   ', &
     '             AALD --> CMON + MEO2 + ROO + POX                                                       ', &
     '         ALDX + O --> HOX + CXO3                                                                    ', &
     '       ALDX + HOX --> CXO3                                                                          ', &
     '      ALDX + NTOX --> NTRC + CXO3                                                                   ', &
     '             ALDX --> CMON + AALD + XO2H + ROO + POX                                                ', &
     '       GLYD + HOX --> 0.2 GLY + 0.2 POX + 0.8 ACOO                                                  ', &
     '             GLYD --> 0.15 MEOH + 0.89 CMON + 0.11 GLY + 0.74 FORM + 0.11 XO2H ... etc.             ', &
     '      GLYD + NTOX --> NTRC + ACOO                                                                   ', &
     '        GLY + HOX --> 1.8 CMON + 0.2 XO2R + 0.2 ROO + POX                                           ', &
     '              GLY --> 2 CMON + 2 POX                                                                ', &
     '       GLY + NTOX --> 1.5 CMON + NTRC + 0.5 XO2R + 0.5 ROO + POX                                    ', &
     '             MEGY --> CMON + POX + ACOO                                                             ', &
     '      MEGY + NTOX --> NTRC + XO2R + ROO + ACOO                                                      ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_4 = (/ &
     '       MEGY + HOX --> CMON + ACOO                                                                   ', &
     '       HOX + DIHY --> POX                                                                           ', &
     '       CMON + HOX --> POX                                                                           ', &
     '       MTHA + HOX --> MEO2 + ROO                                                                    ', &
     '       ETHA + HOX --> 0.991 AALD + 0.009 XO2N + 0.991 XO2H + ROO                                    ', &
     '       MEOH + HOX --> FORM + POX                                                                    ', &
     '       ETOH + HOX --> 0.011 GLYD + 0.078 FORM + 0.95 AALD + 0.1 XO2H + 0.1 ROO ... etc.             ', &
     '              KET --> 0.5 AALD + 0.5 XO2H + 0.5 MEO2 - -2.5 ALKA + ROO + 0.5 CXO3 ... etc.          ', &
     '             ACET --> 0.38 CMON + 1.38 MEO2 + 1.38 ROO + 0.62 ACOO                                  ', &
     '       ACET + HOX --> FORM + XO2R + ROO + ACOO                                                      ', &
     '       PRPA + HOX --> 0.71 ACET + 0.03 XO2N + 0.97 XO2H + 0.26 ALKA + 0.26 ALDX ... etc.            ', &
     '       ALKA + HOX --> 0.76 ROR + 0.76 XO2R + 0.13 XO2N + 0.11 XO2H - -0.11 ALKA ... etc.            ', &
     '              ROR --> 0.2 KET + 0.42 ACET + 0.02 ROR + 0.74 AALD + 0.04 XO2N ... etc.               ', &
     '         ROR + O2 --> KET + POX                                                                     ', &
     '       ROR + NDOX --> NTR1                                                                          ', &
     '       ETHY + HOX --> 0.3 FACD + 0.3 CMON + 0.7 GLY + 0.7 HOX + 0.3 POX                             ', &
     '         ETHE + O --> CMON + FORM + 0.7 XO2H + 0.3 HOX + 0.7 ROO + POX                              ', &
     '       ETHE + HOX --> 0.22 GLYD + 1.56 FORM + XO2H + ROO                                            ', &
     '       ETHE + OZN --> 0.37 FACD + 0.51 CMON + FORM + 0.16 HOX + 0.16 POX                            ', &
     '      ETHE + NTOX --> 0.5 NTR1 + 1.125 FORM + 0.5 XO2R + 0.5 XO2H + ROO + 0.5 NDOX ... etc.         ', &
     '          OLE + O --> 0.2 CMON + 0.2 FORM + 0.2 AALD + 0.01 XO2N + 0.2 XO2H ... etc.                ', &
     '        OLE + HOX --> 0.781 FORM + 0.488 AALD + 0.195 XO2R + 0.024 XO2N + 0.976 XO2H ... etc.       ', &
     '        OLE + OZN --> 0.04 HPOX + 0.09 FACD + 0.13 AACD + 0.378 CMON + 0.075 GLY ... etc.           ', &
     '       OLE + NTOX --> 0.5 NTR1 + 0.5 FORM + 0.25 AALD + 0.48 XO2R + 0.04 XO2N ... etc.              ', &
     '         IOLE + O --> 0.1 CMON + 1.24 AALD + 0.1 XO2H + 0.1 ALKA + 0.66 ALDX ... etc.               ', &
     '       IOLE + HOX --> 1.3 AALD + XO2H + 0.7 ALDX + ROO                                              ', &
     '       IOLE + OZN --> 0.08 HPOX + 0.08 AACD + 0.245 CMON + 0.24 GLY + 0.06 MEGY ... etc.            ', &
     '      IOLE + NTOX --> 0.5 NTR1 + 0.5 AALD + 0.48 XO2R + 0.04 XO2N + 0.48 XO2H ... etc.              ', &
     '       ISPR + HOX --> ISO2 + ROO                                                                    ', &
     '         ISPR + O --> 0.5 FORM + 0.25 XO2R + 0.75 ISPD + 0.25 ALKA + 0.25 ROO ... etc.              ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_5 = (/ &
     '      ISO2 + NMOX --> 0.1 INTR + 0.673 FORM + 0.082 XO2H + 0.9 ISPD + 0.082 ROO ... etc.            ', &
     '       ISO2 + POX --> 0.88 ISPX + 0.12 FORM + 0.12 ISPD + 0.12 HOX + 0.12 POX ... etc.              ', &
     '      ISO2 + ACOO --> 0.2 AACD + 0.598 FORM + 0.072 XO2H + ISPD + 0.8 MEO2 ... etc.                 ', &
     '       ISO2 + ROO --> 0.598 FORM + 0.072 XO2H + ISPD + 0.072 ROO + 0.728 POX ... etc.               ', &
     '             ISO2 --> HPLD + POX                                                                    ', &
     '       ISPR + OZN --> 0.066 CMON + 0.6 FORM + 0.2 XO2R + 0.65 ISPD + 0.35 ALKA ... etc.             ', &
     '      ISPR + NTOX --> 0.65 NTR2 + 0.35 FORM + 0.33 XO2R + 0.03 XO2N + 0.64 XO2H ... etc.            ', &
     '       ISPD + HOX --> 0.137 ACET + 0.137 CMON + 0.269 GLYD + 0.115 MEGY + 0.521 XO2R ... etc.       ', &
     '       ISPD + OZN --> 0.15 FACD + 0.17 ACET + 0.543 CMON + 0.17 GLY + 0.531 MEGY ... etc.           ', &
     '      ISPD + NTOX --> 0.142 NTR2 + 0.717 NTRC + 0.113 GLYD + 0.113 MEGY + 0.142 XO2R ... etc.       ', &
     '             ISPD --> 0.17 ACET + 0.128 GLYD + 0.24 OLE + 0.26 FORM + 0.16 XO2R ... etc.            ', &
     '       ISPX + HOX --> 0.904 EPOX + 0.029 IOLE + 0.067 ISO2 + 0.029 ALDX + 0.933 HOX ... etc.        ', &
     '             HPLD --> ISPD + HOX                                                                    ', &
     '      HPLD + NTOX --> NTRC + ISPD                                                                   ', &
     '       EPOX + HOX --> EPX2 + ROO                                                                    ', &
     '       EPX2 + POX --> 0.074 FACD + 0.251 CMON + 0.275 GLYD + 0.275 GLY + 0.275 MEGY ... etc.        ', &
     '      EPX2 + NMOX --> 0.251 CMON + 0.275 GLYD + 0.275 GLY + 0.275 MEGY + 0.375 FORM ... etc.        ', &
     '      EPX2 + ACOO --> 0.2 AACD + 0.2 CMON + 0.22 GLYD + 0.22 GLY + 0.22 MEGY ... etc.               ', &
     '       EPX2 + ROO --> 0.251 CMON + 0.275 GLYD + 0.275 GLY + 0.275 MEGY + 0.375 FORM ... etc.        ', &
     '       INTR + HOX --> 0.266 NTR2 + 0.185 FACD + 0.104 INTR + 0.331 GLYD + 0.098 OLE ... etc.        ', &
     '         TERP + O --> 5.12 ALKA + 0.15 ALDX                                                         ', &
     '       TERP + HOX --> 0.28 FORM + 0.5 XO2R + 0.25 XO2N + 0.75 XO2H + 1.66 ALKA ... etc.             ', &
     '       TERP + OZN --> 0.001 CMON + 0.24 FORM + 0.69 XO2R + 0.18 XO2N + 0.07 XO2H ... etc.           ', &
     '      TERP + NTOX --> 0.53 NTR2 + 0.75 XO2R + 0.25 XO2N + 0.28 XO2H + 0.47 ALDX ... etc.            ', &
     '       BENZ + HOX --> 0.352 BZO2 + 0.53 CRSL + 0.118 ROPN + 0.118 HOX + 0.352 ROO ... etc.          ', &
     '      BZO2 + NMOX --> 0.082 NTR1 + 0.918 GLY + 0.918 ROPN + 0.918 POX + 0.918 NDOX ... etc.         ', &
     '      BZO2 + ACOO --> GLY + ROPN + MEO2 + ROO + POX                                                 ', &
     '       BZO2 + POX --> DUMMY2                                                                        ', &
     '       BZO2 + ROO --> GLY + ROPN + ROO + POX                                                        ', &
     '       TOLN + HOX --> 0.65 TOLR + 0.18 CRSL + 0.1 ROPN + 0.07 XO2H + 0.1 HOX ... etc.               ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_6 = (/ &
     '      TOLR + NMOX --> 0.14 NTR2 + 0.417 GLY + 0.2 XOPN + 0.443 MEGY + 0.66 ROPN ... etc.            ', &
     '      TOLR + ACOO --> 0.48 GLY + 0.23 XOPN + 0.52 MEGY + 0.77 ROPN + MEO2 ... etc.                  ', &
     '       TOLR + POX --> DUMMY2                                                                        ', &
     '       TOLR + ROO --> 0.48 GLY + 0.23 XOPN + 0.52 MEGY + 0.77 ROPN + ROO + POX ... etc.             ', &
     '       XYLN + HOX --> 0.544 XYLR + 0.155 CRSL + 0.244 XOPN + 0.058 XO2H + 0.244 HOX ... etc.        ', &
     '      XYLR + NMOX --> 0.14 NTR2 + 0.221 GLY + 0.56 XOPN + 0.675 MEGY + 0.3 ROPN ... etc.            ', &
     '       XYLR + POX --> DUMMY2                                                                        ', &
     '      XYLR + ACOO --> 0.26 GLY + 0.65 XOPN + 0.77 MEGY + 0.35 ROPN + MEO2 ... etc.                  ', &
     '       XYLR + ROO --> 0.26 GLY + 0.65 XOPN + 0.77 MEGY + 0.35 ROPN + ROO + POX ... etc.             ', &
     '       CRSL + HOX --> 0.732 CAT1 + 0.2 CRER + 0.025 GLY + 0.025 ROPN + 0.02 XO2N ... etc.           ', &
     '      CRSL + NTOX --> 0.3 CRER + NTRC + 0.24 GLY + 0.24 MEGY + 0.48 XO2R + 0.48 OPO3 ... etc.       ', &
     '      CRER + NDOX --> CRON                                                                          ', &
     '       CRER + POX --> CRSL                                                                          ', &
     '       CRON + HOX --> NTR2 + 0.5 CRER                                                               ', &
     '      CRON + NTOX --> NTR2 + 0.5 CRER + NTRC                                                        ', &
     '             CRON --> HONO + ROPN + FORM + POX                                                      ', &
     '             XOPN --> 0.7 CMON + 0.4 GLY + XO2H + 0.7 POX + 0.3 ACOO                                ', &
     '       XOPN + HOX --> 0.4 GLY + MEGY + 2 XO2H + 2 ROO                                               ', &
     '       XOPN + OZN --> 0.5 CMON + 1.2 MEGY + 0.1 AALD + 0.3 XO2H + 0.5 HOX ... etc.                  ', &
     '      XOPN + NTOX --> 0.5 NTR2 + 0.25 MEGY + 0.25 ROPN + 0.45 XO2R + 0.1 XO2N ... etc.              ', &
     '             ROPN --> CMON + OPO3 + POX                                                             ', &
     '       ROPN + HOX --> 0.4 GLY + 0.6 OPO3 + 0.4 XO2H + 0.4 ROO                                       ', &
     '       ROPN + OZN --> 1.98 CMON + 1.4 GLY + 0.24 MEGY + 0.08 FORM + 0.02 AALD ... etc.              ', &
     '      ROPN + NTOX --> NTRC + OPO3                                                                   ', &
     '       CAT1 + HOX --> 0.5 CRER + 0.14 FORM + 0.2 POX                                                ', &
     '      CAT1 + NTOX --> CRER + NTRC                                                                   ', &
     '      OPO3 + NMOX --> 0.5 CMON + 0.5 GLY + 0.8 POX + 0.2 CXO3 + NDOX                                ', &
     '      OPO3 + NDOX --> OPAN                                                                          ', &
     '             OPAN --> OPO3 + NDOX                                                                   ', &
     '       OPO3 + POX --> 0.41 PACD + 0.15 AACD + 0.44 XO2H + 0.44 ALDX + 0.44 HOX ... etc.             ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(8) :: EQN_NAMES_7 = (/ &
     '      OPO3 + ACOO --> XO2R + MEO2 + ALDX + 2 ROO                                                    ', &
     '       OPO3 + ROO --> 0.2 AACD + 0.8 XO2H + 0.8 ALDX + 1.8 ROO                                      ', &
     '       OPAN + HOX --> 0.5 NTR2 + CMON + 0.5 GLY + 0.5 NDOX                                          ', &
     '       PANX + HOX --> AALD + NDOX                                                                   ', &
     '             NTR2 --> NTRC                                                                          ', &
     '       ETHE + OZN --> 0.24 CDOX + 0.52 FACD + 0.24 CMON + FORM + 0.12 HOX ... etc.                  ', &
     '        OLE + OZN --> 0.22 CDOX + 0.06 MTHA + 0.01 ETHA + 0.03 MEOH + 0.31 CMON ... etc.            ', &
     '       IOLE + OZN --> 0.18 CDOX + 0.08 MTHA + 0.01 ETHA - -2.26 PRPA + 0.04 MEOH ... etc.           ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(218) :: EQN_NAMES = (/&
    EQN_NAMES_0, EQN_NAMES_1, EQN_NAMES_2, EQN_NAMES_3, EQN_NAMES_4, &
    EQN_NAMES_5, EQN_NAMES_6, EQN_NAMES_7 /)

! INLINED global variables

! End INLINED global variables


END MODULE mod_cb6_Monitor
