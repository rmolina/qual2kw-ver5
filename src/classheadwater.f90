!classhadwater.f90

MODULE Class_Headwater
    USE m_water_quality
    USE Class_SystemParams	!, ONLY:steadystate
    USE Class_RiverTopo	!, ONLY: nHw
    USE Class_Phsolve, ONLY: cT
    IMPLICIT NONE

!	PRIVATE
!	PUBLIC :: Headwater_,

    TYPE Headwater_type
        INTEGER(I4B) :: NumdatPnt = 24					!Number of data points per day
        TYPE(t_water_quality), POINTER :: dat(:)			!headwater data(time:headwaterID)
    END TYPE Headwater_type

CONTAINS

    !HEADWATER DATA STRUCUTRE CONSTRUCTOR
    FUNCTION Headwater_(HwFileIn) RESULT(Hw)
        TYPE(Headwater_type) Hw
        TYPE(t_water_quality), INTENT(IN) :: HwFileIn(:)
        INTEGER(I4B) i, status

!		IF (steadystate()) THEN
        ALLOCATE (Hw%dat(0:Hw%NumdatPnt-1),STAT=status)

        IF (status==1) STOP 'ERROR: Class_Headwater(Headwater_) Insufficient memory for dyanmic allocation'

        Hw%dat = HwFileIn

!		ELSE		!dyanmic simulation
        !ToDo: implement
!			STOP 'dynamic simulation has not been implemented yet!'
!		END IF
    END FUNCTION Headwater_


! Interpolate Instanteneous Headwater data
    SUBROUTINE InstanteousHeadwater(Hw, t, Te, c, pH)
        USE Class_Phsolve, ONLY: cT

        TYPE(Headwater_type), INTENT(IN):: Hw
        REAL(DP), INTENT(IN) :: t
        REAL(DP):: Te, c(:), pH
        REAL(DP) t_hr
        INTEGER(I4B) t0_hr, t1_hr

        t_hr = (t - Int(t)) * Hw%NumdatPnt
        IF (t_hr < Hw%NumdatPnt-1) Then
            t0_hr = Int(t_hr)
            t1_hr = t0_hr + 1
        ELSE
            t0_hr = Hw%NumdatPnt-1
            t1_hr = 0
        END IF

        !gp 12-Jan-06
        !WRITE(10,'(4F13.4)')	Hw%dat(t0_hr)%Te + (t_hr - t0_hr) * (Hw%dat(t1_hr)%Te - Hw%dat(t0_hr)%Te), &
        !						Hw%dat(t0_hr)%c + (t_hr - t0_hr) 	* (Hw%dat(t1_hr)%c - Hw%dat(t0_hr)%c), &
        !						Hw%dat(t0_hr)%pH + (t_hr - t0_hr) * (Hw%dat(t1_hr)%pH - Hw%dat(t0_hr)%pH), &
        !						cT(pH, c(nv-2), Te, c(1))

        Te = Hw%dat(t0_hr)%Te + (t_hr - t0_hr) * (Hw%dat(t1_hr)%Te - Hw%dat(t0_hr)%Te)
        c = Hw%dat(t0_hr)%c + (t_hr - t0_hr) 	* (Hw%dat(t1_hr)%c - Hw%dat(t0_hr)%c)
        pH = Hw%dat(t0_hr)%pH + (t_hr - t0_hr) * (Hw%dat(t1_hr)%pH - Hw%dat(t0_hr)%pH)
        !solve Total Carbon concentration
        c(nv-1)=cT(pH, c(nv-2), Te, c(1))

        !gp 12-Jan-06
        !WRITE(10,'(4F13.4)') Te, c, pH, c(nv-1)

    END SUBROUTINE InstanteousHeadwater


    Function Interpolate_hourly(t, c)
        !for interpolation of instantaneous values from hourly point estimates

        REAL(DP) interpolate_hourly
        REAL(DP),INTENT(IN) :: c(:), t
        REAL(DP) t_hr
        INTEGER(I4B) t0_hr, t1_hr

        t_hr = (t - Int(t)) * 24

        IF (t_hr < 23) Then
            t0_hr = Int(t_hr)
            t1_hr = t0_hr + 1
        ELSE
            t0_hr = 23
            t1_hr = 0
        END IF

        Interpolate_hourly = c(t0_hr) + (t_hr - t0_hr) * (c(t1_hr) - c(t0_hr))

    END FUNCTION


END MODULE Class_Headwater
