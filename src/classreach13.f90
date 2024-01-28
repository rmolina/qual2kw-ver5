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
        REAL(DP) :: xrdn=0.0_DP
!		REAL(DP) xup,xdown, xl, xreach			!xl - reach lengh (m), xreach - reach length (km)
!		REAL(DP) elemLen							!element length
!		INTEGER(I4B) elems						!number of elements in this reach
!		REAL(DP) elev1, elev2
!		REAL(DP) latd, latm, lats
!		REAL(DP) lond, lonm, lons
!		REAL(DP) latr, lonr						!latitude, longitude radius
        !xpm - center of reach,				!
        REAL(DP) :: xpm =0
        INTEGER(I4B) :: upID=0, dwnID=0			!upstream and downstream element ID
    END TYPE Reach_type

END MODULE Class_Reach

