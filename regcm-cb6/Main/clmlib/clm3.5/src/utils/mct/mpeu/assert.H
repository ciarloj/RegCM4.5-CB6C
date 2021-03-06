!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id: assert.H 18 2005-12-12 17:49:42Z mvr $
! CVS $Name$  
!BOP -------------------------------------------------------------------
!
! !MODULE: assert.H - an #include section of ASSERT() macro for Fortran
!
! !DESCRIPTION:
!
! !INTERFACE:
!
!	#define NDEBUG
!	#include "assert.H"
!	...
!		use m_die,only : assert_
!	...
!		ASSERT( <Fortran expression> )
!		ALWAYS_ASSERT( <Fortran expression> )
!
! !BUGS
!	This macro requires Fortran friendly cpp() for macro processing.
!
! !REVISION HISTORY:
! 	28Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- modified
!		- added the prolog for a brief documentation
!	before  - Tom Clune
!		- Created for MP PSAS
!EOP ___________________________________________________________________

	! This implementation allows multi-"#include" in a single file

#ifndef ALWAYS_ASSERT

#define ALWAYS_ASSERT(EX) If (.not. (EX) ) call assert_("EX",__FILE__,__LINE__)
#endif


#ifndef ASSERT

#ifdef NDEBUG

#define ASSERT(EX) ! Skip assertion: EX

#else

#define ASSERT(EX) ALWAYS_ASSERT(EX)

#endif

#endif
