! solarRadiation.f90
! calculate solar radiation

MODULE Class_SolarCalc
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    USE m_constants, only: pii, w0, e
    use m_date, only: date_t
    USE m_solar_position, only: calcJD, degToRad, radToDeg, solarposition, sunrise, sunset
    USE m_RiverTopo, only: t_rivertopo
    use m_hydraulics, only: riverhydraulics_type
    use class_lightheat, only: lightheat
    IMPLICIT NONE
    PRIVATE
	PUBLIC :: Solar_type, SunriseSunset, SolarCalc, sitesolar_

    TYPE solar_type
        !REAL(R64) xyear, xmon, xday, Julday				!Julian day
        REAL(R64) :: nfacBras =2.0_R64, atcRyanStolz=0.8_R64		!Bras parameter, Ryan-Stolz parameter

        !GP 23-Nov-09
        !INTEGER(I32) timezoneHour, dlstime				!timezone hour, daylight saving time in hours
        REAL(R64) timezoneHour, dlstime				!timezone hour, daylight saving time in hours

        CHARACTER(LEN=30) :: solarMethod = "Bras"
        REAL(R64), DIMENSION(:), POINTER :: ff, sunrs, sunss, Jsnt
    END TYPE

!	REAL(R64), ALLOCATABLE ::ff(:), sunrs(:), sunss(:)
!	REAL(R64), ALLOCATABLE ::Jsnt(:)

CONTAINS

    FUNCTION sitesolar_(nr, timezone, solarMethod, fBras, fRyan, dlstime) RESULT(Solar)
        !Data constructor for Solar type

        TYPE(solar_type) Solar

        INTEGER(I32), INTENT(IN) :: nr
        REAL(R64), INTENT(IN) :: fBras, fRyan

        !GP 23-Nov-09
        !INTEGER(I32), INTENT(IN) :: dlstime
        REAL(R64), INTENT(IN) :: dlstime

        !GP 23-Nov-09
        !CHARACTER (*), INTENT(IN) :: timezone, solarMethod
        CHARACTER (*), INTENT(IN) :: solarMethod
        REAL(R64), INTENT(IN) :: timezone

        INTEGER(I32) status(4), i

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
        INTEGER(I32), INTENT(IN) :: nr
        TYPE(date_t), INTENT(IN) :: today						!11/16/04
        TYPE(solar_type), INTENT(INOUT) :: Solar
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau

        REAL(R64) photo, tset, tsun
        INTEGER(I32) i

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
            Solar%sunrs(i) = tsun * 24.0_R64
            Solar%sunss(i) = tset * 24.0_R64
        END DO
    End SUBROUTINE SunriseSunset

    SUBROUTINE SolarCalc(nr, Solar, siteMeteo, hydrau, system)
        USE m_system_params
        USE m_hydraulics
        USE m_meteorology

        INTEGER(I32), INTENT(IN) :: nr
        TYPE(solar_type), INTENT(INOUT) :: Solar
        TYPE(meteorology_t) siteMeteo								!meteology information
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau
        TYPE(system_params_t), INTENT(IN) :: system

        INTEGER(I32) i, j

        DO i=1, nr
            IF ((system%tday>= Solar%sunrs(i)/24.0_R64 - 0.01_R64).AND. &
                (system%tday <= Solar%sunss(i)/24.0_R64 + 0.01_R64)) THEN

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
        TYPE(solar_type), INTENT(IN) :: Solar
        TYPE(date_t), INTENT(IN) :: today							!11/16/04

        !gp 23-Jun-09
        !TYPE(LightHeat_type), INTENT(IN) :: lightheat

        !gp 16-Jul-08
        !REAL(R64), INTENT(IN):: latr, lonr, elev, cc, shadeT, tday	!cc cloudy cover
        REAL(R64), INTENT(IN):: latr, lonr, elev, cc, shadeT, solarT, tday	!cc cloud cover

!		TYPE(Solar), INTENT(IN) :: Solar
        REAL(R64) Jsnt
        REAL(R64) curdayfrac      !current time as fraction of day from 0-1
        REAL(R64) curHourFrac     !intermediate calc
        REAL(R64) curMinFrac      !intermediate calc
        INTEGER(I32) curhh        !current hour during integration
        INTEGER(I32) curmm        !current minute during integration
        REAL(R64) curss           !current second during integration
        REAL(R64) el              !solar elevation (deg from horizon)
        REAL(R64) ALDO            !reflection of solar rad from water surface
        REAL(R64) Iclear
        REAL(R64) az							!solar azimuth in deg from N
        REAL(R64) erv						!earth radius vector distance to sun in AU

        REAL(R64) Julday

        curdayfrac = tday
        curHourFrac = curdayfrac * 24.0_R64 - Int(curdayfrac * 24.0_R64)     !intermediate calc
        curMinFrac = curHourFrac * 60.0_R64 - Int(curHourFrac * 60.0_R64)    !intermediate calc
        curhh = Int(curdayfrac * 24.0_R64)                             !current hour
        curmm = Int(curHourFrac * 60.0_R64)                            !current minute
        curss = curMinFrac * 60.0_R64                                  !current second

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
                Jsnt = Iclear / (4.183076_R64 * 100.0_R64 * 100.0_R64 / 86400.0_R64)   !convert from W/m^2 to cal/cm^2/d
            End If

            !cloud cover correction for solar radiation

            !gp 24-Jun-09
            !Jsnt = Jsnt * (1.0_R64 - 0.65_R64 * cc ** 2)
            Jsnt = Jsnt * (1.0_R64 - lightheat%KCL1 * cc ** 2)

!gp 16-Jul-08
        End If

        !adjustment of Jsnt for effective shading from vegetation and topography
        !note: shadeT(i) at time=t is interpolated from hourly input values in Sub Derivs
        Jsnt = Jsnt * (1.0_R64 - shadeT)

        !gp calculate fraction of solar radiation reflected from the water surface
        !using Anderson (1954) as reported by Brown and Barnwell (1987) for QUAL2e
        If (el > 0) Then
            If (cc < 0.1_R64) Then
                ALDO = 1.18_R64 * el ** (-0.77_R64)
                !	ElseIf ((cc(i) >= 0.1) .And. (cc(i) < 0.5)) Then
            ElseIf (cc < 0.5_R64) Then
                ALDO = 2.2_R64 * el ** (-0.97_R64)
                !	ElseIf ((cc(i) >= 0.5) .And. (cc(i) < 0.9)) Then
            ElseIf (cc < 0.9_R64) Then
                ALDO = 0.95_R64 * el ** (-0.75_R64)
                !	ElseIf (cc(i) >= 0.9) Then
            Else
                ALDO = 0.35_R64 * el ** (-0.45_R64)
            End If
        Else
            ALDO = 0
        End If
        If (ALDO < 0) ALDO = 0
        If (ALDO > 1.0_R64)	ALDO = 1.0_R64
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
        REAL(R64), INTENT(IN) :: jd, year, dayfrac, el, nfac
        REAL(R64) Iclear
        REAL(R64) R, I0, m, a1

        !NREL solar constant (W/m^2)
        !	W0 = 1367

        !ratio of actual earth-sun distance to mean earth-sun distance (Bras eqn 2.10
        !R=Ri(year,jd, dayfrac)

        !solar radiation on horizontal surface at top of atmosphere (Bras eqn 2.9)
        I0 = (W0 / R ** 2) * SIN(degToRad(el))

        !optical air mass (Bras eqn 2.22)
        m = (SIN(degToRad(el)) + 0.15_R64 * (el + 3.885_R64) ** (-1.253_R64)) ** (-1)    !Bras eqn 2.22

        !molecular scattering coeff (Bras eqn 2.26)
        a1 = 0.128_R64 - 0.054_R64 * Log(m) / e

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
        REAL(R64), INTENT(IN) :: jd, year, dayfrac, el, atc, z
        REAL(R64) rs
        REAL(R64) at
        REAL(R64) sinal , al , a0
        REAL(R64) rs_toa , R , rm

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

            rm = (((288.0_R64 - 0.0065_R64 * z) / 288.0_R64) ** 5.256_R64) / (sinal + 0.15_R64 * (a0 + 3.885_R64) ** (-1.253_R64))
            rs_toa = W0 * sinal / R ** 2       !RS on the top of atmosphere
            rs = rs_toa * (at ** rm)           !RS on the ground
        End If
    End FUNCTION RyanStolzSolar

    PURE FUNCTION Ri(year,jd, dayfrac)
        REAL(R64), INTENT(IN) :: year, jd, dayfrac
        REAL(R64) Ri
        INTEGER(I32) days

        If (year / 4 - Int(year / 4) == 0) Then
            days=366
        Else
            days=365
        End If
        Ri = 1 + 0.017_R64 * Cos((2 * PII / days) * (186.0_R64 - jd - 1 + dayfrac))

    END FUNCTION Ri

END MODULE Class_SolarCalc
