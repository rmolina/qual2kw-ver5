!classwaterquality.f90
!water quality constituent type
MODULE Class_WaterQuality
    USE nrtype
    IMPLICIT NONE

    TYPE WaterQuality_type
        REAL(DP) :: Te =0.0
        REAL(DP) :: c(1:nv) =0.0
        REAL(DP) :: pH = 7.0
    END TYPE

END MODULE
