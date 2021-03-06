! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The ODE Function of Chemical Model File
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
! File                 : mod_cbmz_Function.f90
! Time                 : Mon Nov 25 13:41:15 2013
! Working directory    : /scratch/ashalaby/kpp-2.2.3/compare/CBMZ
! Equation file        : mod_cbmz.kpp
! Output root filename : mod_cbmz
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mod_cbmz_Function

  USE mod_cbmz_Parameters
  IMPLICIT NONE

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)


! Computation of equation rates
  A(1) = RCT(1)*V(56)
  A(2) = RCT(2)*V(57)
  A(3) = RCT(3)*V(21)
  A(4) = RCT(4)*V(29)
  A(5) = RCT(5)*V(17)
  A(6) = RCT(6)*V(53)
  A(7) = RCT(7)*V(53)
  A(8) = RCT(8)*V(10)
  A(9) = RCT(9)*V(9)*F(1)
  A(10) = RCT(10)*V(9)*F(2)
  A(11) = 2.2e-10*V(9)*V(35)
  A(12) = RCT(12)*V(27)*F(1)
  A(13) = RCT(13)*V(27)*V(53)
  A(14) = RCT(14)*V(27)*V(56)
  A(15) = RCT(15)*V(27)*V(56)
  A(16) = RCT(16)*V(27)*V(55)
  A(17) = RCT(17)*V(53)*V(55)
  A(18) = RCT(18)*V(53)*V(56)
  A(19) = RCT(19)*V(53)*V(54)
  A(20) = RCT(20)*V(52)*V(53)
  A(21) = RCT(21)*V(18)*V(54)
  A(22) = RCT(22)*V(54)*V(55)
  A(23) = RCT(23)*V(54)*V(56)
  A(24) = 2.2e-11*V(54)*V(57)
  A(25) = RCT(25)*V(21)*V(54)
  A(26) = RCT(26)*V(29)*V(54)
  A(27) = RCT(27)*V(17)*V(54)
  A(28) = RCT(28)*V(52)*V(54)
  A(29) = RCT(29)*V(10)*V(54)
  A(30) = RCT(30)*V(52)*V(52)
  A(31) = RCT(31)*V(35)*V(52)*V(52)
  A(32) = RCT(32)*V(52)*V(55)
  A(33) = RCT(33)*V(52)*V(56)
  A(34) = 5e-16*V(52)*V(56)
  A(35) = RCT(35)*V(17)
  A(36) = RCT(36)*V(55)*V(57)
  A(37) = RCT(37)*V(56)*V(57)
  A(38) = RCT(38)*V(56)*V(57)
  A(39) = RCT(39)*V(57)*V(57)
  A(40) = 3.5e-12*V(52)*V(57)
  A(41) = 2e-21*V(24)*V(35)
  A(42) = RCT(42)*V(24)
  A(43) = RCT(43)*V(55)*V(55)*F(1)
  A(44) = RCT(44)*V(28)*V(54)
  A(45) = RCT(45)*V(11)*V(54)
  A(46) = RCT(46)*V(13)*V(54)
  A(47) = RCT(47)*V(14)*V(54)
  A(48) = 8.1e-13*V(30)*V(54)
  A(49) = RCT(49)*V(20)*V(54)
  A(50) = RCT(50)*V(39)
  A(51) = RCT(51)*V(39)
  A(52) = 1e-11*V(39)*V(54)
  A(53) = RCT(53)*V(39)*V(57)
  A(54) = RCT(54)*V(45)
  A(55) = RCT(55)*V(45)*V(54)
  A(56) = RCT(56)*V(45)*V(57)
  A(57) = RCT(57)*V(36)
  A(58) = RCT(58)*V(36)*V(54)
  A(59) = RCT(59)*V(42)
  A(60) = 1.7e-11*V(42)*V(54)
  A(61) = RCT(61)*V(42)*V(57)
  A(62) = RCT(62)*V(25)*V(53)
  A(63) = RCT(63)*V(25)*V(54)
  A(64) = RCT(64)*V(40)*V(53)
  A(65) = RCT(65)*V(37)*V(53)
  A(66) = RCT(66)*V(40)*V(54)
  A(67) = RCT(67)*V(37)*V(54)
  A(68) = RCT(68)*V(40)*V(57)
  A(69) = 2.5e-12*V(37)*V(57)
  A(70) = RCT(70)*V(8)*V(54)
  A(71) = RCT(71)*V(12)*V(54)
  A(72) = 8.1e-12*V(19)*V(55)
  A(73) = 4.1e-11*V(26)*V(54)
  A(74) = 2.2e-11*V(26)*V(57)
  A(75) = 1.4e-11*V(15)*V(56)
  A(76) = 3e-11*V(31)*V(54)
  A(77) = RCT(77)*V(31)
  A(78) = RCT(78)*V(31)*V(53)
  A(79) = RCT(79)*V(38)*V(54)
  A(80) = RCT(80)*V(38)*V(53)
  A(81) = RCT(81)*V(38)*V(57)
  A(82) = 3.3e-11*V(47)*V(54)
  A(83) = 7e-18*V(47)*V(53)
  A(84) = RCT(84)*V(47)
  A(85) = 1e-15*V(47)*V(57)
  A(86) = RCT(86)*V(22)
  A(87) = RCT(87)*V(23)
  A(88) = RCT(88)*V(49)
  A(89) = RCT(89)*V(22)*V(54)
  A(90) = RCT(90)*V(23)*V(54)
  A(91) = RCT(91)*V(49)*V(54)
  A(92) = RCT(92)*V(51)*V(54)
  A(93) = RCT(93)*V(51)
  A(94) = RCT(94)*V(56)*V(58)
  A(95) = RCT(95)*V(7)
  A(96) = RCT(96)*V(46)*V(55)
  A(97) = RCT(97)*V(43)*V(55)
  A(98) = 4e-12*V(50)*V(55)
  A(99) = RCT(99)*V(55)*V(58)
  A(100) = 4e-12*V(48)*V(55)
  A(101) = 4e-12*V(44)*V(55)
  A(102) = 4e-12*V(33)*V(55)
  A(103) = 4e-12*V(32)*V(55)
  A(104) = 4e-12*V(34)*V(55)
  A(105) = 4e-12*V(41)*V(55)
  A(106) = 1.1e-12*V(46)*V(57)
  A(107) = 2.5e-12*V(43)*V(57)
  A(108) = 2.5e-12*V(50)*V(57)
  A(109) = 4e-12*V(57)*V(58)
  A(110) = 1.2e-12*V(48)*V(57)
  A(111) = 4e-12*V(44)*V(57)
  A(112) = 2.5e-12*V(41)*V(57)
  A(113) = RCT(113)*V(46)*V(52)
  A(114) = RCT(114)*V(43)*V(52)
  A(115) = RCT(115)*V(50)*V(52)
  A(116) = RCT(116)*V(52)*V(58)
  A(117) = RCT(117)*V(48)*V(52)
  A(118) = RCT(118)*V(44)*V(52)
  A(119) = RCT(119)*V(33)*V(52)
  A(120) = RCT(120)*V(32)*V(52)
  A(121) = RCT(121)*V(34)*V(52)
  A(122) = RCT(122)*V(41)*V(52)
  A(123) = RCT(123)*V(16)*V(54)
  A(124) = RCT(124)*V(16)*V(57)

! Aggregate function
  Vdot(1) = 0.24*A(62)+0.22*A(64)+0.18*A(65)+A(99)
  Vdot(2) = A(45)
  Vdot(3) = 0.52*A(62)+0.22*A(64)
  Vdot(4) = 0.09*A(64)+0.16*A(65)+0.39*A(80)+0.46*A(83)+0.4*A(116)
  Vdot(5) = 0.6*A(123)
  Vdot(6) = A(122)
  Vdot(7) = A(94)-A(95)
  Vdot(8) = -A(70)
  Vdot(9) = A(7)-A(9)-A(10)-A(11)
  Vdot(10) = -A(8)-A(29)+A(30)+A(31)
  Vdot(11) = -A(45)+0.4*A(123)+A(124)
  Vdot(12) = -A(71)
  Vdot(13) = -A(46)+0.06*A(64)+0.08*A(65)
  Vdot(14) = -A(47)+0.01*A(64)+0.01*A(65)
  Vdot(15) = 0.4*A(73)+A(74)-A(75)
  Vdot(16) = -A(123)-A(124)
  Vdot(17) = -A(5)-A(27)+A(33)-A(35)
  Vdot(18) = -A(21)+0.08*A(64)
  Vdot(19) = 0.8*A(70)+0.45*A(71)-A(72)
  Vdot(20) = -A(49)+0.03*A(64)+0.04*A(65)
  Vdot(21) = -A(3)+A(22)-A(25)+A(34)
  Vdot(22) = -A(86)-A(89)+A(113)
  Vdot(23) = -A(87)-A(90)+A(114)
  Vdot(24) = A(38)-A(41)-A(42)
  Vdot(25) = -A(62)-A(63)
  Vdot(26) = 0.12*A(70)+0.05*A(71)-A(73)-A(74)
  Vdot(27) = A(1)+0.89*A(2)+A(6)+A(9)+A(10)-A(12)-A(13)-A(14)-A(15)-A(16)
  Vdot(28) = -A(44)+A(50)+A(51)+A(52)+A(53)+A(54)+A(59)+A(61)+0.24*A(62)+0.31*A(64)+0.3*A(65)+2*A(76)+A(77)+0.69*A(78)&
               &+0.07*A(80)+0.1*A(83)+0.33*A(84)+0.64*A(85)+0.59*A(104)
  Vdot(29) = -A(4)+A(23)-A(26)+0.3*A(40)+2*A(41)+A(53)+A(56)+A(61)+A(74)+0.07*A(85)+A(124)
  Vdot(30) = -A(48)-1.06*A(64)-2.26*A(65)-A(66)-2.23*A(67)+1.1*A(71)+1.86*A(85)-1.98*A(88)+0.42*A(91)-1.98*A(93)-1.68&
               &*A(98)-A(101)+0.18*A(102)+1.6*A(103)-1.98*A(108)-A(111)+2*A(120)
  Vdot(31) = 0.95*A(72)+0.3*A(73)-A(76)-A(77)-A(78)
  Vdot(32) = A(81)-A(103)-A(120)
  Vdot(33) = A(79)-A(102)-A(119)
  Vdot(34) = 0.5*A(82)-A(104)-A(121)
  Vdot(35) = -A(11)+A(21)+A(28)-A(31)-A(41)
  Vdot(36) = -A(57)-A(58)+0.07*A(65)+0.23*A(67)+0.09*A(83)+0.03*A(84)+0.74*A(88)+0.74*A(93)+0.62*A(98)+0.63*A(104)+0.74&
               &*A(108)
  Vdot(37) = -A(65)-A(67)-A(69)
  Vdot(38) = -A(79)-A(80)-A(81)
  Vdot(39) = A(49)-A(50)-A(51)-A(52)-A(53)+A(62)+1.56*A(63)+0.57*A(64)+A(66)+A(76)+0.7*A(78)+0.6*A(80)+0.15*A(83)+0.2&
               &*A(84)+0.28*A(85)+A(86)+0.3*A(89)+A(96)+A(100)+0.5*A(101)+0.63*A(102)+0.25*A(104)+A(106)+A(110)+0.5*A(111)
  Vdot(40) = -A(64)-A(66)-A(68)
  Vdot(41) = A(60)+A(63)+A(66)+A(67)+0.08*A(70)+0.5*A(71)+0.6*A(73)+A(76)+0.03*A(78)+0.08*A(79)+0.2*A(80)+0.2*A(82)+0.07&
               &*A(83)+0.93*A(85)+0.4*A(88)+0.41*A(93)+0.34*A(98)-A(105)+0.4*A(108)-A(112)-A(122)
  Vdot(42) = -A(59)-A(60)-A(61)+0.04*A(64)+0.07*A(65)+0.8*A(71)+0.2*A(78)+0.85*A(83)+0.19*A(91)+0.34*A(104)
  Vdot(43) = A(47)+0.06*A(64)+0.05*A(65)+0.1*A(88)+0.7*A(90)+0.1*A(93)-A(97)+0.08*A(98)-A(107)+0.1*A(108)-A(114)
  Vdot(44) = A(68)+A(69)+A(92)-A(101)-A(111)-A(118)
  Vdot(45) = -A(54)-A(55)-A(56)+0.22*A(63)+0.47*A(64)+1.03*A(65)+A(66)+1.77*A(67)+0.03*A(78)+0.15*A(80)+0.02*A(83)+0.07&
               &*A(84)+0.28*A(85)+A(87)+0.3*A(88)+0.3*A(90)+0.04*A(91)+0.3*A(93)+A(97)+0.25*A(98)+0.5*A(101)+0.8*A(103)+0.55&
               &*A(104)+A(107)+0.3*A(108)+0.5*A(111)
  Vdot(46) = A(46)+A(54)+A(57)+0.07*A(64)+0.1*A(65)+0.05*A(83)+0.7*A(84)+0.7*A(89)-A(96)+A(99)-A(106)+A(109)-A(113)
  Vdot(47) = 0.65*A(80)-A(82)-A(83)-A(84)-A(85)+0.91*A(102)+0.2*A(103)
  Vdot(48) = A(58)+0.11*A(65)-A(100)-A(110)-A(117)
  Vdot(49) = -A(88)-A(91)+A(115)+A(117)+A(119)+A(121)
  Vdot(50) = A(48)+0.03*A(64)+0.09*A(65)+0.77*A(91)-A(98)-A(108)-A(115)
  Vdot(51) = 0.05*A(72)+A(75)+0.93*A(85)-A(92)-A(93)+0.16*A(98)+0.5*A(101)+0.09*A(102)+0.8*A(103)+0.5*A(111)+A(118)&
               &+A(120)
  Vdot(52) = A(5)+A(19)-A(20)+A(21)+A(24)-A(28)+A(29)-2*A(30)-2*A(31)-A(32)-A(33)-A(34)+A(35)-A(40)+A(44)+A(45)+A(49)+2&
               &*A(50)+A(52)+A(53)+A(54)+A(59)+0.2*A(62)+A(63)+0.26*A(64)+0.22*A(65)+A(66)+A(67)+0.2*A(70)+0.55*A(71)+0.95&
               &*A(72)+0.6*A(73)+2*A(76)+A(77)+0.76*A(78)+0.07*A(80)+0.1*A(83)+0.33*A(84)+0.93*A(85)+A(86)+A(87)+0.9*A(88)&
               &+0.9*A(93)+A(96)+A(97)+0.76*A(98)+0.5*A(101)+0.91*A(102)+0.8*A(103)+A(104)+A(106)+A(107)+0.9*A(108)+0.5&
               &*A(111)-A(113)-A(114)-A(115)-A(116)-A(117)-A(118)-A(119)-A(120)-A(121)-A(122)
  Vdot(53) = -A(6)-A(7)+A(12)-A(13)-A(17)-A(18)-A(19)-A(20)-A(62)-A(64)-A(65)-A(78)-A(80)-A(83)+0.4*A(116)
  Vdot(54) = A(3)+A(4)+2*A(8)+2*A(11)-A(19)+A(20)-A(21)-A(22)-A(23)-A(24)-A(25)-A(26)-A(27)-A(28)-A(29)+A(32)+0.7*A(40)&
               &-A(44)-A(45)-A(46)-A(47)-A(48)-A(49)-A(52)-A(55)-A(58)-A(60)+0.12*A(62)-A(63)+0.33*A(64)+0.6*A(65)-A(66)&
               &-A(67)-A(70)-A(71)-A(73)-A(76)+0.08*A(78)-A(79)+0.27*A(80)-A(82)+0.27*A(83)+A(86)+A(87)+A(88)-0.7*A(89)-0.7&
               &*A(90)-0.77*A(91)-A(92)-A(123)
  Vdot(55) = A(1)+0.11*A(2)+A(3)+A(14)-A(16)-A(17)-A(22)-A(32)-A(36)+A(37)-2*A(43)-A(72)-A(96)-A(97)-A(98)-A(99)-A(100)&
               &-A(101)-A(102)-A(103)-A(104)-A(105)
  Vdot(56) = -A(1)+0.89*A(2)+A(4)+A(5)-A(14)-A(15)+A(16)+A(17)-A(18)-A(23)+A(24)+A(25)+A(27)+A(32)-A(33)-A(34)+A(35)+2&
               &*A(36)-A(38)+2*A(39)+0.7*A(40)+A(42)+2*A(43)+0.95*A(72)-A(75)+A(93)-A(94)+A(95)+A(96)+A(97)+0.84*A(98)+A(99)&
               &+A(100)+1.5*A(101)+0.91*A(102)+1.2*A(103)+A(104)+A(105)+A(106)+A(107)+A(108)+A(109)+A(110)+1.5*A(111)+A(112)
  Vdot(57) = -A(2)+A(15)+A(18)-A(24)+A(26)-A(36)-A(37)-A(38)-2*A(39)-A(40)+A(42)-A(53)-A(56)-A(61)-A(68)-A(69)-A(74)&
               &-A(81)-A(85)-A(106)-A(107)-A(108)-A(109)-A(110)-A(111)-A(112)-A(124)
  Vdot(58) = A(55)+A(56)+A(57)+A(59)+A(60)+A(61)+0.13*A(64)+0.19*A(65)+A(76)+A(77)+0.62*A(78)+0.2*A(80)+0.5*A(82)+0.11&
               &*A(83)+0.97*A(84)+0.07*A(85)-A(94)+A(95)-A(99)+A(100)-A(109)+A(110)-A(116)

END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE mod_cbmz_Function

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
