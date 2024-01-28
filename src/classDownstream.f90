
MODULE Class_Downstream
    USE nrtype
    USE Class_WaterQuality
!	USE Class_SystemParams	!, ONLY:steadystate
!	USE Class_RiverTopo	!, ONLY: nHw
    IMPLICIT NONE


    TYPE Downstream_type
        LOGICAL(LGT) :: downstreamBound =.FALSE.			!specify downsteam boundary
        INTEGER(I4B) :: NumdatPnt = 24						!Number of data points per day
        TYPE(WaterQuality_type), POINTER :: Dat(:)			!downstream boundary time series data
    END TYPE

CONTAINS

    FUNCTION Downstream_(DBFilein, downstreamBound) RESULT(DB)
        TYPE(WaterQuality_type), INTENT(IN) :: DBFilein(0:)
        LOGICAL(LGT), INTENT(IN) :: downstreamBound
        TYPE(Downstream_type) DB
        INTEGER(I4B) :: status

        IF (downstreamBound) THEN
            ALLOCATE (DB%dat(0:DB%NumdatPnt-1), STAT=status)
            IF (status==1) THEN
                STOP 'QUAL2K ERROR: Class_Headwater(Downstream_) Insufficient memory for dyanmic allocation'
            END IF
            DB%Dat=DBFilein
        ELSE
        END IF
        DB%downstreamBound= downstreamBound
    END FUNCTION Downstream_


    ! interpolate instanteous downstream boundary condition
    SUBROUTINE InstanteousDownstreamboundary(DB, t, lastTe, lastc, Te, c)
        USE Class_Phsolve, ONLY: cT
        TYPE(Downstream_type), INTENT(IN) :: DB

!		TYPE(HeadwaterDownstream_type), INTENT(IN):: HDboundary
        REAL(DP), INTENT(IN) :: lastTe, lastc(:), t		!last reach temperature
        REAL(DP) Te, c(:)
        REAL(DP) t_hr, pH
        INTEGER(I4B) t0_hr, t1_hr

        !gp 12-Jan-06
        !WRITE(10,'(4F13.4)')	DB%dat(t0_hr)%Te + (t_hr - t0_hr) * (DB%dat(t1_hr)%Te - DB%dat(t0_hr)%Te), &
        !						DB%dat(t0_hr)%c + (t_hr - t0_hr)   * (DB%dat(t1_hr)%c - DB%dat(t0_hr)%c), &
        !						DB%dat(t0_hr)%pH + (t_hr - t0_hr) * (DB%dat(t1_hr)%pH - DB%dat(t0_hr)%pH), &
        !						cT(pH, c(nv - 2), Te, c(1))

        IF (DB%downstreamBound) Then						!specified downstream boundary
            t_hr = (t - Int(t)) * DB%NumdatPnt
            IF (t_hr < DB%NumdatPnt -1) THEN
                t0_hr = INT(t_hr)
                t1_hr = t0_hr + 1
            ELSE

                !gp 17-Jan-06 debug
                !t0_hr = DB%NumdatPnt
                t0_hr = DB%NumdatPnt-1

                t1_hr = 0
            END IF

            !Te = Hw%dat(t0_hr)%Te + (t_hr - t0_hr) * (Hw%dat(t1_hr)%Te - Hw%dat(t0_hr)%Te)
            Te = DB%dat(t0_hr)%Te + (t_hr - t0_hr) * (DB%dat(t1_hr)%Te - DB%dat(t0_hr)%Te)

            !c = Hw%dat(t0_hr)%c + (t_hr - t0_hr) 	* (Hw%dat(t1_hr)%c - Hw%dat(t0_hr)%c)
            c  = DB%dat(t0_hr)%c + (t_hr - t0_hr)   * (DB%dat(t1_hr)%c - DB%dat(t0_hr)%c)

            !pH = Hw%dat(t0_hr)%pH + (t_hr - t0_hr) * (Hw%dat(t1_hr)%pH - Hw%dat(t0_hr)%pH)
            pH = DB%dat(t0_hr)%pH + (t_hr - t0_hr) * (DB%dat(t1_hr)%pH - DB%dat(t0_hr)%pH)

            !c(nv-1)=cT(pH, c(nv-2), Te, c(1))
            c(nv - 1) = cT(pH, c(nv - 2), Te, c(1))

        ELSE																!same as last reach's
            Te=lastTe
            c=lastc
        END IF
    END SUBROUTINE InstanteousDownstreamboundary

END MODULE Class_Downstream
