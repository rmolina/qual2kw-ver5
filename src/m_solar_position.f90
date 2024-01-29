module m_solar_position
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_constants, only: pii
    implicit none
    private
    public solarelevation, radtodeg, degtorad, sunrise, sunset, solarnoon, solarposition, calcjd

    !the sunrise/sunset and solar position functions are a vba translation of
    !noaa!s sunrise/sunset calculator and noaa!s solar position calculator
    !at the following web pages:
    !
    !http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html
    !http://www.srrb.noaa.gov/highlights/sunrise/azel.html
    !
    !the calculations in the noaa sunrise/sunset and solar position calculators
    !are based on equations from astronomical algorithms, by jean meeus.
    !the sunrise and sunset results have been verified by noaa to be accurate
    !to within a minute for locations between +/- 72� latitude,
    !and within 10 minutes outside of those latitudes.
    !
    !five main functions are included for use from excel worksheets or vba programs:
    !
    !sunrise(lat, lon, year, month, day, timezone, dlstime)
    ! calculates the local time of sunrise for a location and date
    !
    !solarnoon(lat, lon, year, month, day, timezone, dlstime)
    ! calculates the local time of solar noon for a location and date
    ! (the time when the sun crosses the meridian)
    !
    !sunset(lat, lon, year, month, day, timezone, dlstime)
    ! calculates the local time of sunset for a location and date
    !
    !solarazimuth(lat, lon, year, month, day, hour, minute, second, timezone, dlstime)
    ! calculates the solar azimuth for a location, date, and time
    ! (degrees clockwise from north to the point on the horizon directly below the sun)
    !
    !solarelevation(lat, lon, year, month, day, hour, minute, second, timezone, dlstime)
    ! calculates the solar elevation for a location, date, and time
    ! (degrees vertically from horizon to the sun)
    !
    !a subroutine is also provided that calculates solar azimuth (az), solar elevation (el):
    !
    !solarposition(lat, lon, year, month, day, hour, minute, second, &
    ! timezone, dlstime, az, el, earthradiusvector)
    !
    !the sign convention for the main functions and subroutine is:
    !
    ! positive latitude decimal degrees for northern hemisphere
    ! negative longitude degrees for western hemisphere
    ! negative time zone hours for western hemisphere
    !
contains
    pure function radtodeg(anglerad)
    !// convert radian angle to degrees
        real(r64) radtodeg
        real(r64), intent(in) :: anglerad
        radtodeg = (180.0_r64 * anglerad /pii)

    end function radtodeg


    pure function degtorad(angledeg)
        real(r64) degtorad
        real(r64), intent(in) :: angledeg
        !// convert degree angle to radians
        degtorad = (pii * angledeg / 180.0)

    end function


    !***********************************************************************/
    !* name: calcjd
    !* type: function
    !* purpose: julian day from calendar day
    !* arguments:
    !* year : 4 digit year
    !* month: january = 1
    !* day : 1 - 31
    !* return value:
    !* the julian day corresponding to the date
    !* note:
    !* number is returned for start of day. fractional days should be
    !* added later.
    !***********************************************************************/
    pure function calcjd(year, month, day)
        real(r64) calcjd
        real(r64), intent(in) :: year, month, day
        real(r64) year1, month1
        real(r64) a, b, jd

        year1=year ; month1=month

        if (month1 <= 2) then
            year1 = year1 - 1.0_r64
            month1 = month1 + 12.0_r64
        end if

        a = floor(year1 / 100.0_r64)
        b = 2.0_r64 - a + floor(a / 4.0_r64)

        jd = floor(365.25_r64 * (year1 + 4716.0_r64)) + &
            floor(30.6001_r64 * (month1 + 1.0_r64)) + day + b - 1524.5_r64
        calcjd = jd

    end function calcjd


    pure function calctimejuliancent(jd)

        !***********************************************************************/
        !* name: calctimejuliancent
        !* type: function
        !* purpose: convert julian day to centuries since j2000.0.
        !* arguments:
        !* jd : the julian day to convert
        !* return value:
        !* the t value corresponding to the julian day
        !***********************************************************************/

    real(r64) calctimejuliancent
        real(r64), intent(in) :: jd
        real(r64) t

        t = (jd - 2451545.0_r64) / 36525.0_r64
        calctimejuliancent = t

    end function

    pure function calcjdfromjuliancent(t)

        !***********************************************************************/
        !* name: calcjdfromjuliancent
        !* type: function
        !* purpose: convert centuries since j2000.0 to julian day.
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* the julian day corresponding to the t value
        !***********************************************************************/
        real(r64) calcjdfromjuliancent
        real(r64), intent(in) :: t
        real(r64) jd

        jd = t * 36525.0_r64 + 2451545.0_r64
        calcjdfromjuliancent = jd

    end function


    pure function calcgeommeanlongsun(t)

        !***********************************************************************/
        !* name: calgeommeanlongsun
        !* type: function
        !* purpose: calculate the geometric mean longitude of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* the geometric mean longitude of the sun in degrees
        !***********************************************************************/
        real(r64) calcgeommeanlongsun
        real(r64), intent(in) :: t
        real(r64) l0

        l0 = 280.46646_r64 + t * (36000.76983_r64 + 0.0003032_r64 * t)
        do
            if (l0 > 360.0_r64) then
                l0 = l0 - 360.0_r64
            else if (l0 < 0) then
                l0 = l0 + 360.0_r64
            else
                exit
            end if
        end do

        calcgeommeanlongsun = l0

    end function


    pure function calcgeommeananomalysun(t)

        !***********************************************************************/
        !* name: calgeomanomalysun
        !* type: function
        !* purpose: calculate the geometric mean anomaly of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* the geometric mean anomaly of the sun in degrees
        !***********************************************************************/
        real(r64) calcgeommeananomalysun
        real(r64), intent(in) :: t
        real(r64) m

        m = 357.52911_r64 + t * (35999.05029_r64 - 0.0001537_r64 * t)
        calcgeommeananomalysun = m

    end function


    pure function calceccentricityearthorbit(t)

        !***********************************************************************/
        !* name: calceccentricityearthorbit
        !* type: function
        !* purpose: calculate the eccentricity of earth!s orbit
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* the unitless eccentricity
        !***********************************************************************/
        real(r64) calceccentricityearthorbit
        real(r64), intent(in) :: t
        real(r64) e0

        e0 = 0.016708634_r64 - t * (0.000042037_r64 + 0.0000001267_r64 * t)
        calceccentricityearthorbit = e0

    end function


    pure function calcsuneqofcenter(t)


        !***********************************************************************/
        !* name: calcsuneqofcenter
        !* type: function
        !* purpose: calculate the equation of center for the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* in degrees
        !***********************************************************************/
        real(r64) calcsuneqofcenter
        real(r64), intent(in) :: t
        real(r64) m, mrad, sinm, sin2m, sin3m
        real(r64) c

        m = calcgeommeananomalysun(t)

        mrad = degtorad(m)
        sinm = sin(mrad)
        sin2m = sin(mrad + mrad)
        sin3m = sin(mrad + mrad + mrad)

        c = sinm * (1.914602_r64 - t * (0.004817_r64 + 0.000014_r64 * t)) &
            + sin2m * (0.019993_r64 - 0.000101_r64 * t) + sin3m * 0.000289_r64

        calcsuneqofcenter = c

    end function


    pure function calcsuntruelong(t)

        !***********************************************************************/
        !* name: calcsuntruelong
        !* type: function
        !* purpose: calculate the true longitude of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun!s true longitude in degrees
        !***********************************************************************/
        real(r64) calcsuntruelong
        real(r64), intent(in) :: t

        real(r64) l0, c, o

        l0 = calcgeommeanlongsun(t)
        c = calcsuneqofcenter(t)

        o = l0 + c
        calcsuntruelong = o

    end function


    pure function calcsuntrueanomaly(t)

        !***********************************************************************/
        !* name: calcsuntrueanomaly (not used by sunrise, solarnoon, sunset)
        !* type: function
        !* purpose: calculate the true anamoly of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun!s true anamoly in degrees
        !***********************************************************************/
        real(r64) calcsuntrueanomaly
        real(r64), intent(in) :: t

        real(r64) m, c, v

        m = calcgeommeananomalysun(t)
        c = calcsuneqofcenter(t)

        v = m + c
        calcsuntrueanomaly = v

    end function


    pure function calcsunradvector(t)

        !***********************************************************************/
        !* name: calcsunradvector (not used by sunrise, solarnoon, sunset)
        !* type: function
        !* purpose: calculate the distance to the sun in au
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun radius vector in aus
        !***********************************************************************/
        real(r64) calcsunradvector
        real(r64), intent(in) :: t
        real(r64) v, e0, r

        v = calcsuntrueanomaly(t)
        e0 = calceccentricityearthorbit(t)

        r = (1.000001018_r64 * (1.0_r64 - e0 * e0)) / (1 + e0 * cos(degtorad(v)))
        calcsunradvector = r

    end function


    pure function calcsunapparentlong(t)

        !***********************************************************************/
        !* name: calcsunapparentlong (not used by sunrise, solarnoon, sunset)
        !* type: function
        !* purpose: calculate the apparent longitude of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun!s apparent longitude in degrees
        !***********************************************************************/
        real(r64) calcsunapparentlong
        real(r64), intent(in) :: t
        real(r64) o, omega, lambda

        o = calcsuntruelong(t)

        omega = 125.04_r64 - 1934.136_r64 * t
        lambda = o - 0.00569_r64 - 0.00478_r64 * sin(degtorad(omega))
        calcsunapparentlong = lambda

    end function


    pure function calcmeanobliquityofecliptic(t)

        !***********************************************************************/
        !* name: calcmeanobliquityofecliptic
        !* type: function
        !* purpose: calculate the mean obliquity of the ecliptic
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* mean obliquity in degrees
        !***********************************************************************/
        real(r64) calcmeanobliquityofecliptic
        real(r64), intent(in) :: t
        real(r64) seconds, e0

        seconds = 21.448_r64 - t * (46.815_r64 + t * (0.00059_r64 - t * (0.001813_r64)))
        e0 = 23.0_r64 + (26.0_r64 + (seconds / 60.0_r64)) / 60.0_r64
        calcmeanobliquityofecliptic = e0

    end function


    pure function calcobliquitycorrection(t)

        !***********************************************************************/
        !* name: calcobliquitycorrection
        !* type: function
        !* purpose: calculate the corrected obliquity of the ecliptic
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* corrected obliquity in degrees
        !***********************************************************************/
        real(r64) calcobliquitycorrection
        real(r64), intent(in) :: t
        real(r64) e0, omega, e1

        e0 = calcmeanobliquityofecliptic(t)

        omega = 125.04_r64 - 1934.136_r64 * t
        e1 = e0 + 0.00256_r64 * cos(degtorad(omega))
        calcobliquitycorrection = e1

    end function calcobliquitycorrection


    pure function calcsunrtascension(t)

        !***********************************************************************/
        !* name: calcsunrtascension (not used by sunrise, solarnoon, sunset)
        !* type: function
        !* purpose: calculate the right ascension of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun!s right ascension in degrees
        !***********************************************************************/
        real(r64) calcsunrtascension
        real(r64), intent(in) :: t
        real(r64) e0, lambda, tananum, tanadenom, alpha

        e0 = calcobliquitycorrection(t)
        lambda = calcsunapparentlong(t)

        tananum = (cos(degtorad(e0)) * sin(degtorad(lambda)))
        tanadenom = (cos(degtorad(lambda)))

        !original noaa code using javascript math.atan2(y,x) convention:
        ! var alpha = radtodeg(math.atan2(tananum, tanadenom));
        ! alpha = radtodeg(atan2(tananum, tanadenom))

        !translated using excel vba atan2(x,y) convention:
        alpha = radtodeg(atan2(tanadenom, tananum))

        calcsunrtascension = alpha

    end function calcsunrtascension


    pure function calcsundeclination(t)

        !***********************************************************************/
        !* name: calcsundeclination
        !* type: function
        !* purpose: calculate the declination of the sun
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* sun!s declination in degrees
        !***********************************************************************/
        real(r64) calcsundeclination
        real(r64), intent(in) :: t
        real(r64) e0, lambda, sint, theta

        e0 = calcobliquitycorrection(t)
        lambda = calcsunapparentlong(t)

        sint = sin(degtorad(e0)) * sin(degtorad(lambda))
        theta = radtodeg(asin(sint))
        calcsundeclination = theta

    end function calcsundeclination


    pure function calcequationoftime(t)

        !***********************************************************************/
        !* name: calcequationoftime
        !* type: function
        !* purpose: calculate the difference between true solar time and mean
        !* solar time
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* return value:
        !* equation of time in minutes of time
        !***********************************************************************/
        real(r64) calcequationoftime
        real(r64), intent(in) :: t

        real(r64) epsilon, l0, e0, m
        real(r64) y, sin2l0, sinm
        real(r64) cos2l0, sin4l0, sin2m, etime

        epsilon = calcobliquitycorrection(t)
        l0 = calcgeommeanlongsun(t)
        e0 = calceccentricityearthorbit(t)
        m = calcgeommeananomalysun(t)

        y = tan(degtorad(epsilon) / 2.0_r64)
        y = y * y

        sin2l0 = sin(2.0_r64 * degtorad(l0))
        sinm = sin(degtorad(m))
        cos2l0 = cos(2.0_r64 * degtorad(l0))
        sin4l0 = sin(4.0_r64 * degtorad(l0))
        sin2m = sin(2.0_r64 * degtorad(m))

        etime = y * sin2l0 - 2.0_r64 * e0 * sinm + 4.0_r64 * e0 * y * sinm * cos2l0 &
            - 0.5_r64 * y * y * sin4l0 - 1.25_r64 * e0 * e0 * sin2m

        calcequationoftime = radtodeg(etime) * 4.0_r64

    end function calcequationoftime


    pure function calchouranglesunrise(lat, solardec)

        !***********************************************************************/
        !* name: calchouranglesunrise
        !* type: function
        !* purpose: calculate the hour angle of the sun at sunrise for the
        !* latitude
        !* arguments:
        !* lat : latitude of observer in degrees
        !* solardec : declination angle of sun in degrees
        !* return value:
        !* hour angle of sunrise in radians
        !***********************************************************************/
        real(r64) calchouranglesunrise
        real(r64), intent(in) :: lat, solardec

        real(r64) latrad, sdrad, haarg, ha

        latrad = degtorad(lat)
        sdrad = degtorad(solardec)

        haarg = (cos(degtorad(90.833_r64)) / (cos(latrad) * cos(sdrad)) - tan(latrad) * tan(sdrad))

        ha = (acos(cos(degtorad(90.833_r64)) &
            / (cos(latrad) * cos(sdrad)) - tan(latrad) * tan(sdrad)))

        calchouranglesunrise = ha

    end function


    function calchouranglesunset(lat, solardec)

        !***********************************************************************/
        !* name: calchouranglesunset
        !* type: function
        !* purpose: calculate the hour angle of the sun at sunset for the
        !* latitude
        !* arguments:
        !* lat : latitude of observer in degrees
        !* solardec : declination angle of sun in degrees
        !* return value:
        !* hour angle of sunset in radians
        !***********************************************************************/
        real(r64) calchouranglesunset
        real(r64), intent(in) :: lat, solardec

        real(r64) latrad, sdrad, haarg, ha

        latrad = degtorad(lat)
        sdrad = degtorad(solardec)

        haarg = (cos(degtorad(90.833_r64)) / (cos(latrad) * cos(sdrad)) - tan(latrad) * tan(sdrad))

        ha = (acos(cos(degtorad(90.833_r64)) &
            / (cos(latrad) * cos(sdrad)) - tan(latrad) * tan(sdrad)))

        calchouranglesunset = -ha

    end function


    function calcsunriseutc(jd, latitude, longitude)

        !***********************************************************************/
        !* name: calcsunriseutc
        !* type: function
        !* purpose: calculate the universal coordinated time (utc) of sunrise
        !* for the given day at the given location on earth
        !* arguments:
        !* jd : julian day
        !* latitude : latitude of observer in degrees
        !* longitude : longitude of observer in degrees
        !* return value:
        !* time in minutes from zero z
        !***********************************************************************/
        real(r64) calcsunriseutc
        real(r64), intent(in) :: jd, latitude, longitude

        real(r64) t, eqtime, solardec, hourangle
        real(r64) delta, timediff, timeutc
        real(r64) newt

        t = calctimejuliancent(jd)

        ! // *** first pass to approximate sunrise

        eqtime = calcequationoftime(t)
        solardec = calcsundeclination(t)
        hourangle = calchouranglesunrise(latitude, solardec)

        delta = longitude - radtodeg(hourangle)
        timediff = 4.0_r64 * delta
        ! in minutes of time
        timeutc = 720.0_r64 + timediff - eqtime
        ! in minutes

        ! *** second pass includes fractional jday in gamma calc

        newt = calctimejuliancent(calcjdfromjuliancent(t) + timeutc / 1440.0_r64)
        eqtime = calcequationoftime(newt)
        solardec = calcsundeclination(newt)
        hourangle = calchouranglesunrise(latitude, solardec)
        delta = longitude - radtodeg(hourangle)
        timediff = 4.0_r64 * delta
        timeutc = 720.0_r64 + timediff - eqtime
        ! in minutes

        calcsunriseutc = timeutc

    end function


    pure function calcsolnoonutc(t, longitude)

        !***********************************************************************/
        !* name: calcsolnoonutc
        !* type: function
        !* purpose: calculate the universal coordinated time (utc) of solar
        !* noon for the given day at the given location on earth
        !* arguments:
        !* t : number of julian centuries since j2000.0
        !* longitude : longitude of observer in degrees
        !* return value:
        !* time in minutes from zero z
        !***********************************************************************/
        real(r64) calcsolnoonutc
        real(r64), intent(in) :: t, longitude
        real(r64) newt, eqtime, solarnoondec, solnoonutc

        newt = calctimejuliancent(calcjdfromjuliancent(t) + 0.5_r64 + longitude / 360.0_r64)
        eqtime = calcequationoftime(newt)
        solarnoondec = calcsundeclination(newt)
        solnoonutc = 720 + (longitude * 4) - eqtime

        calcsolnoonutc = solnoonutc

    end function


    function calcsunsetutc(jd, latitude, longitude)

        !***********************************************************************/
        !* name: calcsunsetutc
        !* type: function
        !* purpose: calculate the universal coordinated time (utc) of sunset
        !* for the given day at the given location on earth
        !* arguments:
        !* jd : julian day
        !* latitude : latitude of observer in degrees
        !* longitude : longitude of observer in degrees
        !* return value:
        !* time in minutes from zero z
        !***********************************************************************/
        real(r64) calcsunsetutc
        real(r64), intent(in) :: jd, latitude, longitude
        real(r64) t, eqtime, solardec, hourangle
        real(r64) delta, timediff, timeutc
        real(r64) newt

        t = calctimejuliancent(jd)

        ! // first calculates sunrise and approx length of day

        eqtime = calcequationoftime(t)
        solardec = calcsundeclination(t)
        hourangle = calchouranglesunset(latitude, solardec)

        delta = longitude - radtodeg(hourangle)
        timediff = 4 * delta
        timeutc = 720 + timediff - eqtime

        ! // first pass used to include fractional day in gamma calc

        newt = calctimejuliancent(calcjdfromjuliancent(t) + timeutc / 1440.0_r64)
        eqtime = calcequationoftime(newt)
        solardec = calcsundeclination(newt)
        hourangle = calchouranglesunset(latitude, solardec)

        delta = longitude - radtodeg(hourangle)
        timediff = 4 * delta
        timeutc = 720 + timediff - eqtime
        ! // in minutes

        calcsunsetutc = timeutc

    end function calcsunsetutc


    function sunrise(lat, lon, year, month, day, timezone, dlstime)

        !***********************************************************************/
        !* name: sunrise
        !* type: main function called by spreadsheet
        !* purpose: calculate time of sunrise for the entered date
        !* and location.
        !* for latitudes greater than 72 degrees n and s, calculations are
        !* accurate to within 10 minutes. for latitudes less than +/- 72�
        !* accuracy is approximately one minute.
        !* arguments:
        ! latitude = latitude (decimal degrees)
        ! longitude = longitude (decimal degrees)
        ! note: longitude is negative for western hemisphere for input cells
        ! in the spreadsheet for calls to the functions named
        ! sunrise, solarnoon, and sunset. those functions convert the
        ! longitude to positive for the western hemisphere for calls to
        ! other functions using the original sign convention
        ! from the noaa javascript code.
        ! year = year
        ! month = month
        ! day = day
        ! timezone = time zone hours relative to gmt/utc (hours)
        ! dlstime = daylight savings time (0 = no, 1 = yes) (hours)
        !* return value:
        !* sunrise time in local time (days)
        !***********************************************************************/
        real(r64) sunrise
        real(r64), intent(in) :: lat, lon, year, month, day

        !gp 23-nov-09
        !integer(i32), intent(in) :: timezone, dlstime
        real(r64), intent(in) :: timezone, dlstime

        real(r64) longitude, latitude, jd
        real(r64) risetimegmt, risetimelst

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8_r64) then
            latitude = 89.8_r64
        else if (latitude < -89.8_r64) then
            latitude = -89.8_r64
        end if
        jd = calcjd(year, month, day)

        ! // calculate sunrise for this date
        risetimegmt = calcsunriseutc(jd, latitude, longitude)

        ! // adjust for time zone and daylight savings time in minutes
        risetimelst = risetimegmt + (60 * timezone) + (dlstime * 60)

        ! // convert to days
        sunrise = risetimelst / 1440.0_r64

    end function


    pure function solarnoon(lat, lon, year, month, day, timezone, dlstime)

        !***********************************************************************/
        !* name: solarnoon
        !* type: main function called by spreadsheet
        !* purpose: calculate the universal coordinated time (utc) of solar
        !* noon for the given day at the given location on earth
        !* arguments:
        ! year
        ! month
        ! day
        !* longitude : longitude of observer in degrees
        ! note: longitude is negative for western hemisphere for input cells
        ! in the spreadsheet for calls to the functions named
        ! sunrise, solarnoon, and sunset. those functions convert the
        ! longitude to positive for the western hemisphere for calls to
        ! other functions using the original sign convention
        ! from the noaa javascript code.
        !* return value:
        !* time of solar noon in local time days
        !***********************************************************************/
        real(r64) solarnoon
        real(r64), intent(in) :: lat, lon, year, month, day

        !gp 23-nov-09
        !integer(i32), intent(in) :: timezone, dlstime
        real(r64), intent(in) :: timezone, dlstime

        real(r64) longitude, latitude, jd
        real(r64) t, newt, eqtime
        real(r64) solarnoondec, solnoonutc

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8) then
            latitude = 89.8
        else if (latitude < -89.8) then
            latitude = -89.8
        end if
        jd = calcjd(year, month, day)
        t = calctimejuliancent(jd)

        newt = calctimejuliancent(calcjdfromjuliancent(t) + 0.5_r64 + longitude / 360.0_r64)

        eqtime = calcequationoftime(newt)
        solarnoondec = calcsundeclination(newt)
        solnoonutc = 720 + (longitude * 4) - eqtime

! // adjust for time zone and daylight savings time in minutes
        solarnoon = solnoonutc + (60 * timezone) + (dlstime * 60)

! // convert to days
        solarnoon = solarnoon / 1440.0_r64

    end function


    function sunset(lat, lon, year, month, day, timezone, dlstime)

        !***********************************************************************/
        !* name: sunset
        !* type: main function called by spreadsheet
        !* purpose: calculate time of sunrise and sunset for the entered date
        !* and location.
        !* for latitudes greater than 72 degrees n and s, calculations are
        !* accurate to within 10 minutes. for latitudes less than +/- 72�
        !* accuracy is approximately one minute.
        !* arguments:
        ! latitude = latitude (decimal degrees)
        ! longitude = longitude (decimal degrees)
        ! note: longitude is negative for western hemisphere for input cells
        ! in the spreadsheet for calls to the functions named
        ! sunrise, solarnoon, and sunset. those functions convert the
        ! longitude to positive for the western hemisphere for calls to
        ! other functions using the original sign convention
        ! from the noaa javascript code.
        ! year = year
        ! month = month
        ! day = day
        ! timezone = time zone hours relative to gmt/utc (hours)
        ! dlstime = daylight savings time (0 = no, 1 = yes) (hours)
        !* return value:
        !* sunset time in local time (days)
        !***********************************************************************/

        real(r64) sunset
        real(r64), intent(in) :: lat, lon, year, month, day

        !gp 23-nov-09
        !integer(i32), intent(in) :: timezone, dlstime
        real(r64), intent(in) :: timezone, dlstime

        real(r64) longitude, latitude, jd
        real(r64) settimegmt, settimelst

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8) then
            latitude = 89.8
        else if (latitude < -89.8) then
            latitude = -89.8
        end if
        jd = calcjd(year, month, day)

        ! // calculate sunset for this date
        settimegmt = calcsunsetutc(jd, latitude, longitude)

        ! // adjust for time zone and daylight savings time in minutes
        settimelst = settimegmt + (60 * timezone) + (dlstime * 60)

        ! // convert to days
        sunset = settimelst / 1440.0_r64

    end function


    function solarazimuth(lat, lon, year, month, day, &
        hours, minutes, seconds, timezone, dlstime)

        !***********************************************************************/
        !* name: solarazimuth
        !* type: main function
        !* purpose: calculate solar azimuth (deg from north) for the entered
        !* date, time and location. returns -999999 if darker than twilight
        !*
        !* arguments:
        !* latitude, longitude, year, month, day, hour, minute, second,
        !* timezone, daylightsavingstime
        !* return value:
        !* solar azimuth in degrees from north
        !*
        !* note: solarelevation and solarazimuth functions are identical
        !* and could be converted to a vba subroutine that would return
        !* both values.
        !*
        !***********************************************************************/
        real(r64) solarazimuth
        real(r64), intent(in) ::lat, lon, year, month, day, &
            hours, minutes, seconds, timezone, dlstime
        real(r64) longitude, latitude
        real(r64) zone, daysavings
        real(r64) hh, mm, ss, timenow
        real(r64) jd, t, r
        real(r64) alpha, theta, etime, eqtime
        real(r64) solardec, earthradvec, solartimefix
        real(r64) truesolartime, hourangle, harad
        real(r64) csz, zenith, azdenom, azrad
        real(r64) azimuth, exoatmelevation
        real(r64) step1, step2, step3
        real(r64) refractioncorrection, te, solarzen

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8) then
            latitude = 89.8
        else if (latitude < -89.8) then
            latitude = -89.8_r64
        end if

        !change time zone to ppositive hours in western hemisphere
        zone = timezone * (-1)
        daysavings = dlstime * 60
        hh = hours - (daysavings / 60.0_r64)
        mm = minutes
        ss = seconds

        !// timenow is gmt time for calculation in hours since 0z
        timenow = hh + mm / 60.0_r64 + ss / 3600.0_r64 + zone

        jd = calcjd(year, month, day)
        t = calctimejuliancent(jd + timenow / 24.0)
        r = calcsunradvector(t)
        alpha = calcsunrtascension(t)
        theta = calcsundeclination(t)
        etime = calcequationoftime(t)

        eqtime = etime
        solardec = theta !// in degrees
        earthradvec = r

        solartimefix = eqtime - 4.0 * longitude + 60.0 * zone
        truesolartime = hh * 60.0 + mm + ss / 60.0 + solartimefix
        !// in minutes

        do while (truesolartime > 1440)
            truesolartime = truesolartime - 1440
        end do

        hourangle = truesolartime / 4.0 - 180.0
        !// thanks to louis schwarzmayr for the next line:
        if (hourangle < -180) hourangle = hourangle + 360.0

        harad = degtorad(hourangle)

        csz = sin(degtorad(latitude)) * &
            sin(degtorad(solardec)) + &
            cos(degtorad(latitude)) * &
            cos(degtorad(solardec)) * cos(harad)

        if (csz > 1.0) then
            csz = 1.0
        elseif (csz < -1.0) then
            csz = -1.0
        end if

        zenith = radtodeg(acos(csz))

        azdenom = (cos(degtorad(latitude)) * sin(degtorad(zenith)))

        if (abs(azdenom) > 0.001) then
            azrad = ((sin(degtorad(latitude)) * &
                cos(degtorad(zenith))) - &
                sin(degtorad(solardec))) / azdenom
            if (abs(azrad) > 1.0) then
                if (azrad < 0) then
                    azrad = -1.0
                else
                    azrad = 1.0
                end if
            end if

            azimuth = 180.0 - radtodeg(acos(azrad))

            if (hourangle > 0.0) then
                azimuth = -azimuth
            end if
        else
            if (latitude > 0.0) then
                azimuth = 180.0
            else
                azimuth = 0.0
            end if
        end if
        if (azimuth < 0.0) then
            azimuth = azimuth + 360.0
        end if

        exoatmelevation = 90.0 - zenith

        !beginning of complex expression commented out
        ! if (exoatmelevation > 85.0) then
        ! refractioncorrection = 0.0
        ! else
        ! te = tan(degtorad(exoatmelevation))
        ! if (exoatmelevation > 5.0) then
        ! refractioncorrection = 58.1 / te - 0.07 / (te * te * te) + &
        ! 0.000086 / (te * te * te * te * te)
        ! elseif (exoatmelevation > -0.575) then
        ! refractioncorrection = 1735.0 + exoatmelevation * &
        ! (-518.2 + exoatmelevation * (103.4 + &
        ! exoatmelevation * (-12.79 + &
        ! exoatmelevation * 0.711)))
        ! else
        ! refractioncorrection = -20.774 / te
        ! end if
        ! refractioncorrection = refractioncorrection / 3600.0
        ! end if
        !end of complex expression

        !beginning of simplified expression
        if (exoatmelevation > 85.0) then
            refractioncorrection = 0.0
        else
            te = tan(degtorad(exoatmelevation))
            if (exoatmelevation > 5.0_r64) then
                refractioncorrection = 58.1 / te - 0.07_r64 / (te * te * te) + &
                    0.000086_r64 / (te * te * te * te * te)
            elseif (exoatmelevation > -0.575) then
                step1 = (-12.79_r64 + exoatmelevation * 0.711_r64)
                step2 = (103.4_r64 + exoatmelevation * (step1))
                step3 = (-518.2_r64 + exoatmelevation * (step2))
                refractioncorrection = 1735.0_r64 + exoatmelevation * (step3)
            else
                refractioncorrection = -20.774 / te
            end if
            refractioncorrection = refractioncorrection / 3600.0
        end if
        !end of simplified expression

        solarzen = zenith - refractioncorrection

        ! if (solarzen < 108.0) then
        solarazimuth = azimuth
        ! solarelevation = 90.0 - solarzen
        ! if (solarzen < 90.0) then
        ! coszen = cos(degtorad(solarzen))
        ! else
        ! coszen = 0.0
        ! end if
        ! else !// do not report az & el after astro twilight
        ! solarazimuth = -999999
        ! solarelevation = -999999
        ! coszen = -999999
        ! end if

    end function


    pure function solarelevation(lat, lon, year, month, day, &
        hours, minutes, seconds, timezone, dlstime)

        !***********************************************************************/
        !* name: solarazimuth
        !* type: main function
        !* purpose: calculate solar azimuth (deg from north) for the entered
        !* date, time and location. returns -999999 if darker than twilight
        !*
        !* arguments:
        !* latitude, longitude, year, month, day, hour, minute, second,
        !* timezone, daylightsavingstime
        !* return value:
        !* solar azimuth in degrees from north
        !*
        !* note: solarelevation and solarazimuth functions are identical
        !* and could converted to a vba subroutine that would return
        !* both values.
        !*
        !***********************************************************************/
        real(r64) solarelevation
        real(r64), intent(in) ::lat, lon, year, month, day, seconds

        !gp 23-nov-09
        !integer(i32), intent(in) :: hours, minutes, timezone, dlstime
        integer(i32), intent(in) :: hours, minutes
        real(r64), intent(in) :: timezone, dlstime

        real(r64) longitude, latitude
        real(r64) zone, daysavings
        real(r64) hh, mm, ss, timenow
        real(r64) jd, t, r
        real(r64) alpha, theta, etime, eqtime
        real(r64) solardec, earthradvec, solartimefix
        real(r64) truesolartime, hourangle, harad
        real(r64) csz, zenith, azdenom, azrad
        real(r64) azimuth, exoatmelevation
        real(r64) step1, step2, step3
        real(r64) refractioncorrection, te, solarzen

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8) latitude = 89.8
        if (latitude < -89.8) latitude = -89.8

        !change time zone to ppositive hours in western hemisphere
        zone = timezone * (-1)
        daysavings = dlstime * 60.0
        hh = hours - (daysavings / 60.0)
        mm = minutes
        ss = seconds

        !// timenow is gmt time for calculation in hours since 0z
        timenow = hh + mm / 60.0 + ss / 3600.0 + zone

        jd = calcjd(year, month, day)
        t = calctimejuliancent(jd + timenow / 24.0)
        r = calcsunradvector(t)
        alpha = calcsunrtascension(t)
        theta = calcsundeclination(t)
        etime = calcequationoftime(t)

        eqtime = etime
        solardec = theta !// in degrees
        earthradvec = r

        solartimefix = eqtime - 4.0 * longitude + 60.0 * zone
        truesolartime = hh * 60.0 + mm + ss / 60.0 + solartimefix
        !// in minutes

        do while (truesolartime > 1440)
            truesolartime = truesolartime - 1440
        end do

        hourangle = truesolartime / 4.0 - 180.0
        !// thanks to louis schwarzmayr for the next line:
        if (hourangle < -180) hourangle = hourangle + 360.0

        harad = degtorad(hourangle)

        csz = sin(degtorad(latitude)) * &
            sin(degtorad(solardec)) + &
            cos(degtorad(latitude)) * &
            cos(degtorad(solardec)) * cos(harad)

        if (csz > 1.0) then
            csz = 1.0
        elseif (csz < -1.0) then
            csz = -1.0
        end if

        zenith = radtodeg(acos(csz))

        azdenom = (cos(degtorad(latitude)) * sin(degtorad(zenith)))

        if (abs(azdenom) > 0.001) then
            azrad = ((sin(degtorad(latitude)) * &
                cos(degtorad(zenith))) - &
                sin(degtorad(solardec))) / azdenom
            if (abs(azrad) > 1.0) then
                if (azrad < 0) then
                    azrad = -1.0
                else
                    azrad = 1.0
                end if
            end if

            azimuth = 180.0 - radtodeg(acos(azrad))

            if (hourangle > 0.0) then
                azimuth = -azimuth
            end if
        else
            if (latitude > 0.0) then
                azimuth = 180.0
            else
                azimuth = 0.0
            end if
        end if
        if (azimuth < 0.0) then
            azimuth = azimuth + 360.0
        end if

        exoatmelevation = 90.0 - zenith

        !beginning of complex expression commented out
        ! if (exoatmelevation > 85.0) then
        ! refractioncorrection = 0.0
        ! else
        ! te = tan(degtorad(exoatmelevation))
        ! if (exoatmelevation > 5.0) then
        ! refractioncorrection = 58.1 / te - 0.07 / (te * te * te) + &
        ! 0.000086 / (te * te * te * te * te)
        ! elseif (exoatmelevation > -0.575) then
        ! refractioncorrection = 1735.0 + exoatmelevation * &
        ! (-518.2 + exoatmelevation * (103.4 + &
        ! exoatmelevation * (-12.79 + &
        ! exoatmelevation * 0.711)))
        ! else
        ! refractioncorrection = -20.774 / te
        ! end if
        ! refractioncorrection = refractioncorrection / 3600.0
        ! end if
        !end of complex expression

        !beginning of simplified expression
        if (exoatmelevation > 85.0) then
            refractioncorrection = 0.0
        else
            te = tan(degtorad(exoatmelevation))
            if (exoatmelevation > 5.0) then
                refractioncorrection = 58.1 / te - 0.07 / (te * te * te) + &
                    0.000086 / (te * te * te * te * te)
            elseif (exoatmelevation > -0.575) then
                step1 = (-12.79 + exoatmelevation * 0.711)
                step2 = (103.4 + exoatmelevation * (step1))
                step3 = (-518.2 + exoatmelevation * (step2))
                refractioncorrection = 1735.0 + exoatmelevation * (step3)
            else
                refractioncorrection = -20.774 / te
            end if
            refractioncorrection = refractioncorrection / 3600.0
        end if
        !end of simplified expression

        solarzen = zenith - refractioncorrection

        ! if (solarzen < 108.0) then
        ! solarazimuth = azimuth
        solarelevation = 90.0 - solarzen
        ! if (solarzen < 90.0) then
        ! coszen = cos(degtorad(solarzen))
        ! else
        ! coszen = 0.0
        ! end if
        ! else !// do not report az & el after astro twilight
        ! solarazimuth = -999999
        ! solarelevation = -999999
        ! coszen = -999999
        ! end if

    end function solarelevation


    subroutine solarposition(lat, lon, year, month, day, &
        hours, minutes, seconds, timezone, dlstime, &
        solarazimuth0, solarelevation0, earthradvec)

        !***********************************************************************/
        !* name: solarposition
        !* type: subroutine
        !* purpose: calculate solar azimuth (deg from north)
        !* and elevation (deg from horizeon) for the entered
        !* date, time and location.
        !*
        !* arguments:
        !* latitude, longitude, year, month, day, hour, minute, second,
        !* timezone, daylightsavingstime
        !* return value:
        !* solar azimuth in degrees from north
        !* solar elevation in degrees from horizon
        !* earth radius vector (distance to the sun in au)
        !*
        !* note: solarelevation and solarazimuth functions are identical
        !* and could converted to a vba subroutine that would return
        !* both values.
        !*
        !***********************************************************************/
        real(r64) lat, lon, year, month, day

        !gp 23-nov-09
        !integer(i32) hours, minutes, timezone, dlstime
        integer(i32) hours, minutes
        real(r64) timezone, dlstime

        real(r64) seconds, solarazimuth0, solarelevation0, earthradvec
        real(r64) longitude, latitude
        real(r64) zone, daysavings
        real(r64) hh, mm, ss, timenow
        real(r64) jd, t, r
        real(r64) alpha, theta, etime, eqtime
        !real(r64) solardec, earthradvec, solartimefix
        real(r64) solardec, solartimefix
        real(r64) truesolartime, hourangle, harad
        real(r64) csz, zenith, azdenom, azrad
        real(r64) azimuth, exoatmelevation
        real(r64) step1, step2, step3
        real(r64) refractioncorrection, te, solarzen

        ! change sign convention for longitude from negative to positive in western hemisphere
        longitude = lon * (-1)
        latitude = lat
        if (latitude > 89.8) then
            latitude = 89.8
        else if (latitude < -89.8) then
            latitude = -89.8
        end if


        !change time zone to ppositive hours in western hemisphere
        zone = timezone * (-1)
        daysavings = dlstime * 60
        hh = hours - (daysavings / 60)
        mm = minutes
        ss = seconds

        !// timenow is gmt time for calculation in hours since 0z
        timenow = hh + mm / 60 + ss / 3600 + zone

        jd = calcjd(year, month, day)
        t = calctimejuliancent(jd + timenow / 24.0)
        r = calcsunradvector(t)
        alpha = calcsunrtascension(t)
        theta = calcsundeclination(t)
        etime = calcequationoftime(t)

        eqtime = etime
        solardec = theta !// in degrees
        earthradvec = r

        solartimefix = eqtime - 4.0 * longitude + 60.0 * zone
        truesolartime = hh * 60.0 + mm + ss / 60.0 + solartimefix
        !// in minutes

        do while (truesolartime > 1440)
            truesolartime = truesolartime - 1440
        end do

        hourangle = truesolartime / 4.0 - 180.0
        !// thanks to louis schwarzmayr for the next line:
        if (hourangle < -180) hourangle = hourangle + 360.0

        harad = degtorad(hourangle)

        csz = sin(degtorad(latitude)) * &
            sin(degtorad(solardec)) + &
            cos(degtorad(latitude)) * &
            cos(degtorad(solardec)) * cos(harad)

        if (csz > 1.0) then
            csz = 1.0
        elseif (csz < -1.0) then
            csz = -1.0
        end if

        zenith = radtodeg(acos(csz))

        azdenom = (cos(degtorad(latitude)) * sin(degtorad(zenith)))

        if (abs(azdenom) > 0.001) then
            azrad = ((sin(degtorad(latitude)) * &
                cos(degtorad(zenith))) - &
                sin(degtorad(solardec))) / azdenom
            if (abs(azrad) > 1.0) then
                if (azrad < 0) then
                    azrad = -1.0
                else
                    azrad = 1.0
                end if
            end if

            azimuth = 180.0 - radtodeg(acos(azrad))

            if (hourangle > 0.0) then
                azimuth = -azimuth
            end if
        else
            if (latitude > 0.0) then
                azimuth = 180.0
            else
                azimuth = 0.0
            end if
        end if
        if (azimuth < 0.0) then
            azimuth = azimuth + 360.0
        end if

        exoatmelevation = 90.0 - zenith

        !beginning of complex expression commented out
        ! if (exoatmelevation > 85.0) then
        ! refractioncorrection = 0.0
        ! else
        ! te = tan(degtorad(exoatmelevation))
        ! if (exoatmelevation > 5.0) then
        ! refractioncorrection = 58.1 / te - 0.07 / (te * te * te) + &
        ! 0.000086 / (te * te * te * te * te)
        ! elseif (exoatmelevation > -0.575) then
        ! refractioncorrection = 1735.0 + exoatmelevation * &
        ! (-518.2 + exoatmelevation * (103.4 + &
        ! exoatmelevation * (-12.79 + &
        ! exoatmelevation * 0.711)))
        ! else
        ! refractioncorrection = -20.774 / te
        ! end if
        ! refractioncorrection = refractioncorrection / 3600.0
        ! end if
        !end of complex expression


        !beginning of simplified expression
        if (exoatmelevation > 85.0) then
            refractioncorrection = 0.0
        else
            te = tan(degtorad(exoatmelevation))
            if (exoatmelevation > 5.0) then
                refractioncorrection = 58.1 / te - 0.07 / (te * te * te) + &
                    0.000086 / (te * te * te * te * te)
            elseif (exoatmelevation > -0.575) then
                step1 = (-12.79 + exoatmelevation * 0.711)
                step2 = (103.4 + exoatmelevation * (step1))
                step3 = (-518.2 + exoatmelevation * (step2))
                refractioncorrection = 1735.0 + exoatmelevation * (step3)
            else
                refractioncorrection = -20.774 / te
            end if
            refractioncorrection = refractioncorrection / 3600.0
        end if
        !end of simplified expression


        solarzen = zenith - refractioncorrection

        ! if (solarzen < 108.0) then
        solarazimuth0 = azimuth
        solarelevation0 = 90.0 - solarzen
        ! if (solarzen < 90.0) then
        ! coszen = cos(degtorad(solarzen))
        ! else
        ! coszen = 0.0
        ! end if
        ! else !// do not report az & el after astro twilight
        ! solarazimuth = -999999
        ! solarelevation = -999999
        ! coszen = -999999
        ! end if

    end subroutine

end module m_solar_position







