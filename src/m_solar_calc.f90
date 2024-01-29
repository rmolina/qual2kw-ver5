! solarradiation.f90
! calculate solar radiation

module m_solar_calc
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_constants, only: pii, w0, e
    use m_date, only: date_t
    use m_solar_position, only: calcjd, degtorad, radtodeg, solarposition, sunrise, sunset
    use m_rivertopo, only: t_rivertopo
    use m_hydraulics, only: riverhydraulics_type
    use m_light_heat, only: lightheat
    implicit none
    private
    public :: solar_type, sunrisesunset, solarcalc, sitesolar_

    type solar_type
        !real(r64) xyear, xmon, xday, julday !julian day
        real(r64) :: nfacbras =2.0_r64, atcryanstolz=0.8_r64 !bras parameter, ryan-stolz parameter

        !gp 23-nov-09
        !integer(i32) timezonehour, dlstime !timezone hour, daylight saving time in hours
        real(r64) timezonehour, dlstime !timezone hour, daylight saving time in hours

        character(len=30) :: solarmethod = "Bras"
        real(r64), dimension(:), pointer :: ff, sunrs, sunss, jsnt
    end type

! real(r64), allocatable ::ff(:), sunrs(:), sunss(:)
! real(r64), allocatable ::jsnt(:)

contains

    function sitesolar_(nr, timezone, solarmethod, fbras, fryan, dlstime) result(solar)
        !data constructor for solar type

        type(solar_type) solar

        integer(i32), intent(in) :: nr
        real(r64), intent(in) :: fbras, fryan

        !gp 23-nov-09
        !integer(i32), intent(in) :: dlstime
        real(r64), intent(in) :: dlstime

        !gp 23-nov-09
        !character (*), intent(in) :: timezone, solarmethod
        character (*), intent(in) :: solarmethod
        real(r64), intent(in) :: timezone

        integer(i32) status(4), i

        allocate (solar%ff(0:nr), stat=status(1))
        allocate (solar%sunrs(0:nr), stat=status(2))
        allocate (solar%sunss(0:nr), stat=status(3))
        allocate (solar%jsnt(0:nr), stat=status(4))

        do i=1, 4
            if (status(i)==1) then
                stop 'Class_SolarCalc:sitesolar_ failed. Insufficient Memory!'
            end if
        end do

        !gp 23-nov-09
        !select case (timezone)
        ! case ("atlantic")
        ! solar%timezonehour = -4
        ! case ("eastern")
        ! solar%timezonehour = -5
        ! case ("central")
        ! solar%timezonehour = -6
        ! case ("mountain")
        ! solar%timezonehour = -7
        ! case ("pacific")
        ! solar%timezonehour = -8
        ! case ("alaska")
        ! solar%timezonehour = -9
        ! case ("hawaii-aleutian")
        ! solar%timezonehour = -10
        ! case ("samoa")
        ! solar%timezonehour = -11
        ! case default !gp 17-nov-04 allow user to enter any integer hour time zone (e.g. pst=-8, gmt/utc=0, etc)
        ! if (len_trim(timezone)==0) then
        ! solar%timezonehour = 0 !time zone is gmt/utc if it is left blank in the xls file
        ! else
        ! open (unit=9, file='c:\qual2kw5\scratch.q2k', status='scratch', action='readwrite')
        ! write(9,*) timezone !write character value to scratch file
        ! rewind(9)
        ! read(9,*) solar%timezonehour !read integer value from scratch file
        ! close (9)
        ! if (solar%timezonehour < -12 .or. solar%timezonehour > 14) then
        ! stop !invalid input for time zone. select from the list or enter an integer hour (e.g. pst=-8, gmt=0, etc)
        ! end if
        ! end if
        !end select
        solar%timezonehour = timezone

        solar%solarmethod=solarmethod
        solar%nfacbras=fbras
        solar%atcryanstolz = fryan
        solar%dlstime=dlstime

    end function sitesolar_

    subroutine sunrisesunset(nr, solar, hydrau, today)
        !gp calculate sunrise, sunset, and photoperid
        integer(i32), intent(in) :: nr
        type(date_t), intent(in) :: today !11/16/04
        type(solar_type), intent(inout) :: solar
        type(riverhydraulics_type), intent(in) :: hydrau

        real(r64) photo, tset, tsun
        integer(i32) i

        do i = 1, nr ! ne() for dynamic simulation

            tsun = sunrise(hydrau%reach(i)%latr, -hydrau%reach(i)%lonr, today%year, today%month, &
                today%day, solar%timezonehour, solar%dlstime) !sunrise in days
            tset = sunset(hydrau%reach(i)%latr, -hydrau%reach(i)%lonr, today%year, today%month, &
                today%day, solar%timezonehour, solar%dlstime) !sunset in days

            photo = tset - tsun !photoperiod in days
            !tnoon = solarnoon(latr, -lonr, solar%xyear, solar%xmon, &
            ! solar%xday, solar%timezonehour, solar%dlstime) !time of solar noon in days

            !gp note: the next three variables are not used for solar elevation/radiation calculations
            !because solar elevation for each segment at each time step is calculated with the new noaa functions
            !instead of the half-sine approximation
            solar%ff(i) = photo
            solar%sunrs(i) = tsun * 24.0_r64
            solar%sunss(i) = tset * 24.0_r64
        end do
    end subroutine sunrisesunset

    subroutine solarcalc(nr, solar, sitemeteo, hydrau, system)
        use m_system_params
        use m_hydraulics
        use m_meteorology

        integer(i32), intent(in) :: nr
        type(solar_type), intent(inout) :: solar
        type(meteorology_t) sitemeteo !meteology information
        type(riverhydraulics_type), intent(in) :: hydrau
        type(system_params_t), intent(in) :: system

        integer(i32) i, j

        do i=1, nr
            if ((system%tday>= solar%sunrs(i)/24.0_r64 - 0.01_r64).and. &
                (system%tday <= solar%sunss(i)/24.0_r64 + 0.01_r64)) then

                !gp 16-jul-08
                !solar%jsnt(i)=solarcalchelper(system%today,hydrau%reach(i)%latr, &
                ! hydrau%reach(i)%lonr, hydrau%reach(i)%elev, sitemeteo%cc(i), &
                ! sitemeteo%shadet(i), system%tday, solar)
                solar%jsnt(i)=solarcalchelper(system%today,hydrau%reach(i)%latr, &
                    hydrau%reach(i)%lonr, hydrau%reach(i)%elev, sitemeteo%cc(i), &
                    sitemeteo%shadet(i), sitemeteo%solart(i), system%tday, solar)

            else
                solar%jsnt(i) = 0
            end if
        end do
    end subroutine solarcalc


    !gp 16-jul-08
    !function solarcalchelper(today, latr, lonr, elev, cc, shadet, tday, &
    ! solar) result(jsnt)
    function solarcalchelper(today, latr, lonr, elev, cc, shadet, solart, tday, &
        solar) result(jsnt)

        !solar radiation at the current reach at this time step
! use class_hydraulics, only: channel
        type(solar_type), intent(in) :: solar
        type(date_t), intent(in) :: today !11/16/04

        !gp 23-jun-09
        !type(lightheat_type), intent(in) :: lightheat

        !gp 16-jul-08
        !real(r64), intent(in):: latr, lonr, elev, cc, shadet, tday !cc cloudy cover
        real(r64), intent(in):: latr, lonr, elev, cc, shadet, solart, tday !cc cloud cover

! type(solar), intent(in) :: solar
        real(r64) jsnt
        real(r64) curdayfrac !current time as fraction of day from 0-1
        real(r64) curhourfrac !intermediate calc
        real(r64) curminfrac !intermediate calc
        integer(i32) curhh !current hour during integration
        integer(i32) curmm !current minute during integration
        real(r64) curss !current second during integration
        real(r64) el !solar elevation (deg from horizon)
        real(r64) aldo !reflection of solar rad from water surface
        real(r64) iclear
        real(r64) az !solar azimuth in deg from n
        real(r64) erv !earth radius vector distance to sun in au

        real(r64) julday

        curdayfrac = tday
        curhourfrac = curdayfrac * 24.0_r64 - int(curdayfrac * 24.0_r64) !intermediate calc
        curminfrac = curhourfrac * 60.0_r64 - int(curhourfrac * 60.0_r64) !intermediate calc
        curhh = int(curdayfrac * 24.0_r64) !current hour
        curmm = int(curhourfrac * 60.0_r64) !current minute
        curss = curminfrac * 60.0_r64  !current second

        !gp for i = 1 to nr

        !gp calculate solar elevation

        !el = solarelevation(latr, -lonr, today%year, today%month, today%day, curhh, curmm, curss, &
        ! solar%timezonehour, solar%dlstime)
        call solarposition(latr, -lonr, today%year, today%month, today%day, curhh, curmm, curss, &
            solar%timezonehour, solar%dlstime, az, el, erv)
!gp 16-jul-08
        if (solar%solarmethod == "Observed") then
            jsnt = solart / (4.183076 * 100 * 100 / 86400) !'convert from w/m^2 to cal/cm^2/d
        else

            !gp optional solar radiation codes for clear sky
            if (el <= 0) then
                jsnt = 0
            else
                julday = calcjd(today%year, today%month, today%day)
                select case (solar%solarmethod)
                  case ("Bras") !clear sky iclear units of w/m^2
                    iclear = brassolar(julday, today%year, tday, el, solar%nfacbras, erv)
                  case default
                    ! ("ryan-stolzenbach")
                    iclear = ryanstolzsolar(julday, today%year, tday, el, solar%atcryanstolz, elev, erv)
                end select
                jsnt = iclear / (4.183076_r64 * 100.0_r64 * 100.0_r64 / 86400.0_r64) !convert from w/m^2 to cal/cm^2/d
            end if

            !cloud cover correction for solar radiation

            !gp 24-jun-09
            !jsnt = jsnt * (1.0_r64 - 0.65_r64 * cc ** 2)
            jsnt = jsnt * (1.0_r64 - lightheat%kcl1 * cc ** 2)

!gp 16-jul-08
        end if

        !adjustment of jsnt for effective shading from vegetation and topography
        !note: shadet(i) at time=t is interpolated from hourly input values in sub derivs
        jsnt = jsnt * (1.0_r64 - shadet)

        !gp calculate fraction of solar radiation reflected from the water surface
        !using anderson (1954) as reported by brown and barnwell (1987) for qual2e
        if (el > 0) then
            if (cc < 0.1_r64) then
                aldo = 1.18_r64 * el ** (-0.77_r64)
                ! elseif ((cc(i) >= 0.1) .and. (cc(i) < 0.5)) then
            elseif (cc < 0.5_r64) then
                aldo = 2.2_r64 * el ** (-0.97_r64)
                ! elseif ((cc(i) >= 0.5) .and. (cc(i) < 0.9)) then
            elseif (cc < 0.9_r64) then
                aldo = 0.95_r64 * el ** (-0.75_r64)
                ! elseif (cc(i) >= 0.9) then
            else
                aldo = 0.35_r64 * el ** (-0.45_r64)
            end if
        else
            aldo = 0
        end if
        if (aldo < 0) aldo = 0
        if (aldo > 1.0_r64) aldo = 1.0_r64
        !adjust for reflection
        jsnt = jsnt * (1 - aldo) !units of cal/cm^2/day

    end function solarcalchelper


    function brassolar(jd, year, dayfrac, el, nfac, r) result(iclear)

        !**********************************************************************
        !*inputs:
        !*jd = julian day (jan 1=1, etc.)
        !*year = current year
        !*dayfrac = current time of day as a fraction of the day (0-1)
        !*el = solar elevation (deg from horizon)
        !*nfac = atmospheric turbidity paramter (2=clear, 5=smoggy"
        !*
        !*output:
        !*iclear = clear-sky solar radiation at input solar elevation (w/m^2)
        !**********************************************************************
        real(r64), intent(in) :: jd, year, dayfrac, el, nfac
        real(r64) iclear
        real(r64) r, i0, m, a1

        !nrel solar constant (w/m^2)
        ! w0 = 1367

        !ratio of actual earth-sun distance to mean earth-sun distance (bras eqn 2.10
        !r=ri(year,jd, dayfrac)

        !solar radiation on horizontal surface at top of atmosphere (bras eqn 2.9)
        i0 = (w0 / r ** 2) * sin(degtorad(el))

        !optical air mass (bras eqn 2.22)
        m = (sin(degtorad(el)) + 0.15_r64 * (el + 3.885_r64) ** (-1.253_r64)) ** (-1) !bras eqn 2.22

        !molecular scattering coeff (bras eqn 2.26)
        a1 = 0.128_r64 - 0.054_r64 * log(m) / e

        !clear-sky solar radiation at earth surface on horizontal surface (w/m^2) (bras eqn 2.25)
        iclear = i0 * exp(-nfac * a1 * m)

    end function brassolar


    function ryanstolzsolar(jd, year, dayfrac, el, atc, z, r) result(rs)
        !**************************************************************************
        ! input variables
        ! jd julian day (jan 1=1 etc)
        ! year = current year
        ! dayfrac = current time of day as a fraction of the day (0-1)
        ! el solar elevation deg from horizon
        ! atc atmospheric transmission coefficient (0.70-0.91, default 0.8)
        ! z elevation, metres -- required if imthd=2

        ! output variable
        ! rs clear-sky solar radiation, w m-2
        !***************************************************************************
        real(r64), intent(in) :: jd, year, dayfrac, el, atc, z
        real(r64) rs
        real(r64) at
        real(r64) sinal , al , a0
        real(r64) rs_toa , r , rm

        !nrel solar constant, w m-2
        !w0 = 1367#
        ! atmospheric transmission coefficient (0.70-0.91)
        ! from ryan et al. mit publication
        at = atc
        sinal = sin(degtorad(el)) !sine of the solar elevation angle

        if (sinal < 0) then
            rs = 0.0
        else
            al = asin(sinal)
            a0 = radtodeg(al) !convert the radians to degree
            !ratio of actual earth-sun distance to mean earth-sun distance (bras eqn 2.10
            !r=ri(year,jd, dayfrac)

            rm = (((288.0_r64 - 0.0065_r64 * z) / 288.0_r64) ** 5.256_r64) / (sinal + 0.15_r64 * (a0 + 3.885_r64) ** (-1.253_r64))
            rs_toa = w0 * sinal / r ** 2 !rs on the top of atmosphere
            rs = rs_toa * (at ** rm) !rs on the ground
        end if
    end function ryanstolzsolar

    pure function ri(year,jd, dayfrac)
        real(r64), intent(in) :: year, jd, dayfrac
        real(r64) ri
        integer(i32) days

        if (year / 4 - int(year / 4) == 0) then
            days=366
        else
            days=365
        end if
        ri = 1 + 0.017_r64 * cos((2 * pii / days) * (186.0_r64 - jd - 1 + dayfrac))

    end function ri

end module m_solar_calc
