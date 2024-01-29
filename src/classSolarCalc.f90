! solarRadiation.f90
! calculate solar radiation

MODULE Class_SolarCalc
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    USE nrtype
    use m_date, only: date_t
    USE Class_SolarPosition, only: calcJD, degToRad, radToDeg, solarposition, sunrise, sunset
    USE Class_SystemParams
    USE m_RiverTopo
    IMPLICIT NONE
    PRIVATE SolarCalcHelper

!	PUBLIC :: Solar, SunriseSunset, SolarCalc, sitesolar_

    TYPE solar_type
        !REAL(DP) xyear, xmon, xday, Julday				!Julian day
        REAL(DP) :: nfacBras =2.0_DP, atcRyanStolz=0.8_DP		!Bras parameter, Ryan-Stolz parameter

        !GP 23-Nov-09
        !INTEGER(I4B) timezoneHour, dlstime				!timezone hour, daylight saving time in hours
        REAL(DP) timezoneHour, dlstime				!timezone hour, daylight saving time in hours

        CHARACTER(LEN=30) :: solarMethod = "Bras"
        REAL(DP), DIMENSION(:), POINTER :: ff, sunrs, sunss, Jsnt
    END TYPE

!	REAL(DP), ALLOCATABLE ::ff(:), sunrs(:), sunss(:)
!	REAL(DP), ALLOCATABLE ::Jsnt(:)

CONTAINS

    FUNCTION sitesolar_(nr, timezone, solarMethod, fBras, fRyan, dlstime) RESULT(Solar)
        !Data constructor for Solar type

        TYPE(solar_type) Solar

        INTEGER(I4B), INTENT(IN) :: nr
        REAL(DP), INTENT(IN) :: fBras, fRyan

        !GP 23-Nov-09
        !INTEGER(I4B), INTENT(IN) :: dlstime
        REAL(DP), INTENT(IN) :: dlstime

        !GP 23-Nov-09
        !CHARACTER (*), INTENT(IN) :: timezone, solarMethod
        CHARACTER (*), INTENT(IN) :: solarMethod
        REAL(DP), INTENT(IN) :: timezone

        INTEGER(I4B) status(4), i

        ALLOCATE (Solar%ff(0:nr), STAT=status(1))
        ALLOCATE (Solar%sunrs(0:nr), STAT=status(2))
        ALLOCATE (Solar%sunss(0:nr), STAT=status(3))
        ALLOCATE (Solar%Jsnt(0:nr), STAT=status(4))

        DO i=1, 4
            IF (status(i)==1) THEN
                STOP 'Class_SolarCalc:sitesolar_ failed. Insufficient Memory!'
            END IF
        END DO

        !GP 23-Nov-09
        !Select Case (timezone)
        !	Case ("Atlantic")
        !		Solar%timezoneHour = -4
        !	Case ("Eastern")
        !		Solar%timezoneHour = -5
        !	Case ("Central")
        !		Solar%timezoneHour = -6
        !	Case ("Mountain")
        !		Solar%timezoneHour = -7
        !	Case ("Pacific")
        !		Solar%timezoneHour = -8
        !	Case ("Alaska")
        !		Solar%timezoneHour = -9
        !	Case ("Hawaii-Aleutian")
        !		Solar%timezoneHour = -10
        !	Case ("Samoa")
        !		Solar%timezoneHour = -11
        !	Case Default		!gp 17-Nov-04 allow user to enter any integer hour time zone (e.g. PST=-8, GMT/UTC=0, etc)
        !		IF (LEN_TRIM(timezone)==0) THEN
        !			Solar%timezoneHour = 0		!time zone is GMT/UTC if it is left blank in the xls file
        !		ELSE
        !			OPEN (unit=9, File='c:\qual2kw5\scratch.q2k', status='SCRATCH', ACTION='READWRITE')
        !			WRITE(9,*) timezone				!write character value to scratch file
        !			REWIND(9)
        !			READ(9,*) Solar%timezoneHour	!read integer value from scratch file
        !			CLOSE (9)
        !			IF (Solar%timezoneHour < -12 .or. Solar%timezoneHour > 14) THEN
        !				STOP		!Invalid input for time zone. Select from the list or enter an integer hour (e.g. PST=-8, GMT=0, etc)
        !			END IF
        !		END IF
        !End Select
        Solar%timezoneHour = timezone

        Solar%solarMethod=solarMethod
        Solar%nfacBras=fBras
        Solar%atcRyanStolz = fRyan
        Solar%dlstime=dlstime

    END FUNCTION sitesolar_

    SUBROUTINE SunriseSunset(nr, Solar, hydrau, today)
        !gp calculate sunrise, sunset, and photoperid
        USE m_hydraulics, ONLY: RiverHydraulics_type
        USE m_date
        IMPLICIT NONE

        INTEGER(I4B), INTENT(IN) :: nr
        TYPE(date_t), INTENT(IN) :: today						!11/16/04
        TYPE(solar_type), INTENT(INOUT) :: Solar
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau

        REAL(DP) photo, tset, tsun
        INTEGER(I4B) i

        DO i = 1, nr ! ne() for dynamic simulation

            tsun = sunrise(hydrau%reach(i)%latr, -hydrau%reach(i)%lonr, today%year, today%month, &
                today%day, Solar%timezoneHour, Solar%dlstime)				!sunrise in days
            tset = sunset(hydrau%reach(i)%latr, -hydrau%reach(i)%lonr, today%year, today%month, &
                today%day, Solar%timezoneHour, Solar%dlstime)				!sunset in days

            photo = tset - tsun                !photoperiod in days
            !tnoon = solarnoon(latr, -lonr, Solar%xyear, Solar%xmon, &
            !	Solar%xday, Solar%timezoneHour, Solar%dlstime)				!time of solar noon in days

            !gp note: the next three variables are not used for solar elevation/radiation calculations
            !because solar elevation for each segment at each time step is calculated with the new NOAA functions
            !instead of the half-sine approximation
            Solar%ff(i) = photo
            Solar%sunrs(i) = tsun * 24.0_DP
            Solar%sunss(i) = tset * 24.0_DP
        END DO
    End SUBROUTINE SunriseSunset

    SUBROUTINE SolarCalc(nr, Solar, siteMeteo, hydrau, system)
        USE Class_SystemParams
        USE m_hydraulics
        USE m_meteorology

        INTEGER(I4B), INTENT(IN) :: nr
        TYPE(solar_type), INTENT(INOUT) :: Solar
        TYPE(meteorology_t) siteMeteo								!meteology information
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau
        TYPE(SystemParams), INTENT(IN) :: system

        INTEGER(I4B) i, j

        DO i=1, nr
            IF ((system%tday>= Solar%sunrs(i)/24.0_DP - 0.01_DP).AND. &
                (system%tday <= Solar%sunss(i)/24.0_DP + 0.01_DP)) THEN

                !gp 16-Jul-08
                !Solar%Jsnt(i)=SolarCalcHelper(system%today,hydrau%reach(i)%latr, &
                !			hydrau%reach(i)%lonr, hydrau%reach(i)%elev, siteMeteo%cc(i), &
                !			siteMeteo%shadeT(i), system%tday, Solar)
                Solar%Jsnt(i)=SolarCalcHelper(system%today,hydrau%reach(i)%latr, &
                    hydrau%reach(i)%lonr, hydrau%reach(i)%elev, siteMeteo%cc(i), &
                    siteMeteo%shadeT(i), siteMeteo%solarT(i), system%tday, Solar)

            ELSE
                Solar%Jsnt(i) = 0
            END IF
        END DO
    END SUBROUTINE SolarCalc


    !gp 16-Jul-08
    !FUNCTION SolarCalcHelper(today, latr, lonr, elev, cc, shadeT, tday, &
    !						 Solar) RESULT(Jsnt)
    FUNCTION SolarCalcHelper(today, latr, lonr, elev, cc, shadeT, solarT, tday, &
        Solar) RESULT(Jsnt)

        !solar radiation at the current reach at this time step
!		USE  Class_Hydraulics, ONLY: channel
        USE m_date

        !gp 23-Jun-09
        USE Class_LightHeat

        IMPLICIT NONE
        TYPE(solar_type), INTENT(IN) :: Solar
        TYPE(date_t), INTENT(IN) :: today							!11/16/04

        !gp 23-Jun-09
        !TYPE(LightHeat_type), INTENT(IN) :: lightheat

        !gp 16-Jul-08
        !REAL(DP), INTENT(IN):: latr, lonr, elev, cc, shadeT, tday	!cc cloudy cover
        REAL(DP), INTENT(IN):: latr, lonr, elev, cc, shadeT, solarT, tday	!cc cloud cover

!		TYPE(Solar), INTENT(IN) :: Solar
        REAL(DP) Jsnt
        REAL(DP) curdayfrac      !current time as fraction of day from 0-1
        REAL(DP) curHourFrac     !intermediate calc
        REAL(DP) curMinFrac      !intermediate calc
        INTEGER(I4B) curhh        !current hour during integration
        INTEGER(I4B) curmm        !current minute during integration
        REAL(DP) curss           !current second during integration
        REAL(DP) el              !solar elevation (deg from horizon)
        REAL(DP) ALDO            !reflection of solar rad from water surface
        REAL(DP) Iclear
        REAL(DP) az							!solar azimuth in deg from N
        REAL(DP) erv						!earth radius vector distance to sun in AU

        REAL(DP) Julday

        curdayfrac = tday
        curHourFrac = curdayfrac * 24.0_DP - Int(curdayfrac * 24.0_DP)     !intermediate calc
        curMinFrac = curHourFrac * 60.0_DP - Int(curHourFrac * 60.0_DP)    !intermediate calc
        curhh = Int(curdayfrac * 24.0_DP)                             !current hour
        curmm = Int(curHourFrac * 60.0_DP)                            !current minute
        curss = curMinFrac * 60.0_DP                                  !current second

        !gp For i = 1 To nr

        !gp calculate solar elevation

        !el = solarelevation(latr, -lonr, today%year, today%month, today%day, curhh, curmm, curss, &
        !										Solar%timezoneHour, Solar%dlstime)
        Call solarposition(latr, -lonr, today%year, today%month, today%day, curhh, curmm, curss, &
            Solar%timezoneHour, Solar%dlstime, az, el, erv)
!gp 16-Jul-08
        If (Solar%solarMethod == "Observed") Then
            Jsnt = solarT / (4.183076 * 100 * 100 / 86400)    !'convert from W/m^2 to cal/cm^2/d
        Else

            !gp optional solar radiation codes for clear sky
            If (el <= 0) Then
                Jsnt = 0
            Else
                Julday = calcJD(today%year, today%month, today%day)
                SELECT CASE (Solar%solarMethod)
                  CASE ("Bras")							!clear sky Iclear units of W/m^2
                    Iclear = BrasSolar(Julday, today%year, tday, el, Solar%nfacBras, erv)
                  CASE DEFAULT
                    ! ("Ryan-Stolzenbach")
                    Iclear = RyanStolzSolar(Julday, today%year, tday, el, Solar%atcRyanStolz, elev, erv)
                End SELECT
                Jsnt = Iclear / (4.183076_DP * 100.0_DP * 100.0_DP / 86400.0_DP)   !convert from W/m^2 to cal/cm^2/d
            End If

            !cloud cover correction for solar radiation

            !gp 24-Jun-09
            !Jsnt = Jsnt * (1.0_DP - 0.65_DP * cc ** 2)
            Jsnt = Jsnt * (1.0_DP - lightheat%KCL1 * cc ** 2)

!gp 16-Jul-08
        End If

        !adjustment of Jsnt for effective shading from vegetation and topography
        !note: shadeT(i) at time=t is interpolated from hourly input values in Sub Derivs
        Jsnt = Jsnt * (1.0_DP - shadeT)

        !gp calculate fraction of solar radiation reflected from the water surface
        !using Anderson (1954) as reported by Brown and Barnwell (1987) for QUAL2e
        If (el > 0) Then
            If (cc < 0.1_DP) Then
                ALDO = 1.18_DP * el ** (-0.77_DP)
                !	ElseIf ((cc(i) >= 0.1) .And. (cc(i) < 0.5)) Then
            ElseIf (cc < 0.5_DP) Then
                ALDO = 2.2_DP * el ** (-0.97_DP)
                !	ElseIf ((cc(i) >= 0.5) .And. (cc(i) < 0.9)) Then
            ElseIf (cc < 0.9_DP) Then
                ALDO = 0.95_DP * el ** (-0.75_DP)
                !	ElseIf (cc(i) >= 0.9) Then
            Else
                ALDO = 0.35_DP * el ** (-0.45_DP)
            End If
        Else
            ALDO = 0
        End If
        If (ALDO < 0) ALDO = 0
        If (ALDO > 1.0_DP)	ALDO = 1.0_DP
        !adjust for reflection
        Jsnt = Jsnt * (1 - ALDO)          !units of cal/cm^2/day

    End FUNCTION SolarCalcHelper


    FUNCTION BrasSolar(jd, year, dayfrac, el, nfac, R) RESULT(Iclear)

        !**********************************************************************
        !*inputs:
        !*jd = julian day (Jan 1=1, etc.)
        !*year = current year
        !*dayfrac = current time of day as a fraction of the day (0-1)
        !*el = solar elevation (deg from horizon)
        !*nfac = atmospheric turbidity paramter (2=clear, 5=smoggy"
        !*
        !*output:
        !*Iclear = clear-sky solar radiation at input solar elevation (W/m^2)
        !**********************************************************************
        USE nrtype
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: jd, year, dayfrac, el, nfac
        REAL(DP) Iclear
        REAL(DP) R, I0, m, a1

        !NREL solar constant (W/m^2)
        !	W0 = 1367

        !ratio of actual earth-sun distance to mean earth-sun distance (Bras eqn 2.10
        !R=Ri(year,jd, dayfrac)

        !solar radiation on horizontal surface at top of atmosphere (Bras eqn 2.9)
        I0 = (W0 / R ** 2) * SIN(degToRad(el))

        !optical air mass (Bras eqn 2.22)
        m = (SIN(degToRad(el)) + 0.15_DP * (el + 3.885_DP) ** (-1.253_DP)) ** (-1)    !Bras eqn 2.22

        !molecular scattering coeff (Bras eqn 2.26)
        a1 = 0.128_DP - 0.054_DP * Log(m) / e

        !clear-sky solar radiation at earth surface on horizontal surface (W/m^2) (Bras eqn 2.25)
        Iclear = I0 * Exp(-nfac * a1 * m)

    End FUNCTION BrasSolar


    FUNCTION RyanStolzSolar(jd, year, dayfrac, el, atc, z, R) RESULT(rs)
        !**************************************************************************
        ! input variables
        !   jd      julian day (Jan 1=1 etc)
        !   year = current year
        !   dayfrac = current time of day as a fraction of the day (0-1)
        !   el      solar elevation deg from horizon
        !   atc     atmospheric transmission coefficient (0.70-0.91, default 0.8)
        !   z       elevation, metres -- required if imthd=2

        ! output variable
        !   rs      clear-sky solar radiation, W m-2
        !***************************************************************************
        USE nrtype
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: jd, year, dayfrac, el, atc, z
        REAL(DP) rs
        REAL(DP) at
        REAL(DP) sinal , al , a0
        REAL(DP) rs_toa , R , rm

        !NREL solar constant, W m-2
        !W0 = 1367#
        ! atmospheric transmission coefficient (0.70-0.91)
        ! from ryan et al. mit publication
        at = atc
        sinal = SIN(degToRad(el))						!Sine of the solar elevation angle

        If (sinal < 0) Then
            rs = 0.0
        Else
            al = Asin(sinal)
            a0 = radToDeg(al)                 !convert the radians to degree
            !ratio of actual earth-sun distance to mean earth-sun distance (Bras eqn 2.10
            !R=Ri(year,jd, dayfrac)

            rm = (((288.0_DP - 0.0065_DP * z) / 288.0_DP) ** 5.256_DP) / (sinal + 0.15_DP * (a0 + 3.885_DP) ** (-1.253_DP))
            rs_toa = W0 * sinal / R ** 2       !RS on the top of atmosphere
            rs = rs_toa * (at ** rm)           !RS on the ground
        End If
    End FUNCTION RyanStolzSolar

    PURE FUNCTION Ri(year,jd, dayfrac)
        USE	nrtype
        IMPLICIT NONE

        REAL(DP), INTENT(IN) :: year, jd, dayfrac
        REAL(DP) Ri
        INTEGER(I4B) days

        If (year / 4 - Int(year / 4) == 0) Then
            days=366
        Else
            days=365
        End If
        Ri = 1 + 0.017_DP * Cos((2 * PII / days) * (186.0_DP - jd - 1 + dayfrac))

    END FUNCTION Ri

END MODULE Class_SolarCalc
