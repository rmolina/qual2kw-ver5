!element data structure
!
!
!INCLUDE "NRTYPE.F90"
! classreach.f90

!reach data structure
MODULE Class_Reach
	USE nrtype

	IMPLICIT NONE

	TYPE Reach_type
		CHARACTER(LEN=30) :: rlab='', rname =''
		REAL(DP) :: xrup, xrdn
		REAL(DP) :: xpm =0			!xpm - center of reach
		LOGICAL(LGT) :: IsJct=.FALSE.		!junction reach (with tributary flow comes in)?
		LOGICAL(LGT) :: IsHw=.FALSE.		!headwater reach (its immediate upstream is headwater)?
		INTEGER(I2B) :: upID=0, dwID=0		!upstream and downstream element ID
																			!immediate upstream reach id in the same tributary
!		REAL(DP) xup,xdown, xl, xreach		!xl - reach lengh (m), xreach - reach length (km)
!		REAL(DP) elemLen			!element length
!		INTEGER(I2B) elems			!number of elements in this reach
!		REAL(DP) elev1, elev2					
!		REAL(DP) latd, latm, lats
!		REAL(DP) lond, lonm, lons
!		REAL(DP) latr, lonr			!latitude, longitude radius

	END TYPE Reach_type

END MODULE Class_Reach

