!  solarposition.f90 
!
MODULE Class_SolarPosition
	USE nrtype
	IMPLICIT NONE
	PRIVATE															!All private, unless declare public
	PUBLIC solarelevation, radToDeg, degToRad, sunrise, sunset, solarnoon, solarposition

	!The sunrise/sunset and solar position functions are a VBA translation of
	!NOAA!s sunrise/sunset calculator and NOAA!s solar position calculator
	!at the following web pages:
	!
	!http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html
	!http://www.srrb.noaa.gov/highlights/sunrise/azel.html
	!
	!The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators
	!are based on equations from Astronomical Algorithms, by Jean Meeus.
	!The sunrise and sunset results have been verified by NOAA to be accurate
	!to within a minute for locations between +/- 72° latitude,
	!and within 10 minutes outside of those latitudes.
	!
	!Five main functions are included for use from Excel worksheets or VBA programs:
	!
	!sunrise(lat, lon, year, month, day, timezone, dlstime)
	!     calculates the local time of sunrise for a location and date
	!
	!solarnoon(lat, lon, year, month, day, timezone, dlstime)
	!     calculates the local time of solar noon for a location and date
	!     (the time when the sun crosses the meridian)
	!
	!sunset(lat, lon, year, month, day, timezone, dlstime)
	!     calculates the local time of sunset for a location and date
	!
	!solarazimuth(lat, lon, year, month, day, hour, minute, second, timezone, dlstime)
	!     calculates the solar azimuth for a location, date, and time
	!     (degrees clockwise from north to the point on the horizon directly below the sun)
	!
	!solarelevation(lat, lon, year, month, day, hour, minute, second, timezone, dlstime)
	!     calculates the solar elevation for a location, date, and time
	!     (degrees vertically from horizon to the sun)
	!
	!A subroutine is also provided that calculates solar azimuth (az), solar elevation (el):
	!
	!solarposition(lat, lon, year, month, day, hour, minute, second, &
	!                      timezone, dlstime, az, el, earthRadiusVector)
	!
	!The sign convention for the main functions and subroutine is:
	!
	!     positive latitude decimal degrees for northern hemisphere
	!     negative longitude degrees for western hemisphere
	!     negative time zone hours for western hemisphere
	!
	CONTAINS
		PURE Function radToDeg(angleRad)
			USE nrtype
			IMPLICIT NONE

			!// Convert radian angle to degrees
			REAL(DP) radToDeg
			REAL(DP), INTENT(IN) :: angleRad
			radToDeg = (180.0_dp * angleRad /PII)

		End Function radToDeg


		PURE Function degToRad(angleDeg)
			USE nrtype
			IMPLICIT NONE
			REAL(DP) degToRad
			REAL(DP), INTENT(IN) :: angleDeg
			!// Convert degree angle to radians
			degToRad = (PII * angleDeg / 180.0)
      
		End Function


		!***********************************************************************/
		!* Name:    calcJD
		!* Type:    Function
		!* Purpose: Julian day from calendar day
		!* Arguments:
		!*   year : 4 digit year
		!*   month: January = 1
		!*   day  : 1 - 31
		!* Return value:
		!*   The Julian day corresponding to the date
		!* Note:
		!*   Number is returned for start of day.  Fractional days should be
		!*   added later.
		!***********************************************************************/
		PURE Function calcJD(year, month, day)
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcJD
			REAL(DP), INTENT(IN) ::  year, month, day
			REAL(DP)  year1, month1
			REAL(DP) a, b, jd

			year1=year ; month1=month

			If (month1 <= 2) Then
			 year1 = year1 - 1.0_dp
			 month1 = month1 + 12.0_dp
			End If

			a = Floor(year1 / 100.0_dp)
			b = 2.0_dp - a + Floor(a / 4.0_dp)

			jd = Floor(365.25_dp * (year1 + 4716.0_dp)) + &
					 Floor(30.6001_dp * (month1 + 1.0_dp)) + day + b - 1524.5_dp
			calcJD = jd
        
		End Function calcJD


		PURE Function calcTimeJulianCent(jd)

		!***********************************************************************/
		!* Name:    calcTimeJulianCent
		!* Type:    Function
		!* Purpose: convert Julian Day to centuries since J2000.0.
		!* Arguments:
		!*   jd : the Julian Day to convert
		!* Return value:
		!*   the T value corresponding to the Julian Day
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE

			REAL(DP) calcTimeJulianCent
			REAL(DP), INTENT(IN) ::  jd
			REAL(DP) t

			t = (jd - 2451545.0_dp) / 36525.0_dp
			calcTimeJulianCent = t

		End Function

		PURE Function calcJDFromJulianCent(t)

		!***********************************************************************/
		!* Name:    calcJDFromJulianCent
		!* Type:    Function
		!* Purpose: convert centuries since J2000.0 to Julian Day.
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   the Julian Day corresponding to the t value
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcJDFromJulianCent
		REAL(DP), INTENT(IN) ::  t
		REAL(DP) jd

		jd = t * 36525.0_dp + 2451545.0_dp
		calcJDFromJulianCent = jd

		End Function


		PURE Function calcGeomMeanLongSun(t)

		!***********************************************************************/
		!* Name:    calGeomMeanLongSun
		!* Type:    Function
		!* Purpose: calculate the Geometric Mean Longitude of the Sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   the Geometric Mean Longitude of the Sun in degrees
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcGeomMeanLongSun
		REAL(DP), INTENT(IN) ::  t
		REAL(DP) l0

		l0 = 280.46646_dp + t * (36000.76983_dp + 0.0003032_dp * t)
		Do
			 IF (l0 > 360.0_dp) Then
				l0 = l0 - 360.0_dp
			 ELSE IF (l0 < 0) Then
				l0 = l0 + 360.0_dp
			 ELSE
				 Exit
			 END IF
		END DO

		calcGeomMeanLongSun = l0

		End Function


		PURE Function calcGeomMeanAnomalySun(t)

		!***********************************************************************/
		!* Name:    calGeomAnomalySun
		!* Type:    Function
		!* Purpose: calculate the Geometric Mean Anomaly of the Sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   the Geometric Mean Anomaly of the Sun in degrees
		!***********************************************************************/
		!USE nrtype
		!IMPLICIT NONE
		REAL(DP) calcGeomMeanAnomalySun
		REAL(DP), INTENT(IN) :: t
		REAL(DP) m
  
		m = 357.52911_dp + t * (35999.05029_dp - 0.0001537_dp * t)
		calcGeomMeanAnomalySun = m
      
		End Function
      

		PURE Function calcEccentricityEarthOrbit(t)

		!***********************************************************************/
		!* Name:    calcEccentricityEarthOrbit
		!* Type:    Function
		!* Purpose: calculate the eccentricity of earth!s orbit
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   the unitless eccentricity
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcEccentricityEarthOrbit
			REAL(DP), INTENT(IN) :: t
			REAL(DP) e0

			e0 = 0.016708634_dp - t * (0.000042037_dp + 0.0000001267_dp * t)
			calcEccentricityEarthOrbit = e0
      
		End Function


		PURE Function calcSunEqOfCenter(t)
		

		!***********************************************************************/
		!* Name:    calcSunEqOfCenter
		!* Type:    Function
		!* Purpose: calculate the equation of center for the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunEqOfCenter
			REAL(DP), INTENT(IN) :: t
			REAL(DP) e0
			REAL(DP) m, mrad, sinm, sin2m, sin3m
			REAL(DP) c

					m = calcGeomMeanAnomalySun(t)

					mrad = degToRad(m)
					sinm = Sin(mrad)
					sin2m = Sin(mrad + mrad)
					sin3m = Sin(mrad + mrad + mrad)

					c = sinm * (1.914602_dp - t * (0.004817_dp + 0.000014_dp * t)) &
							+ sin2m * (0.019993_dp - 0.000101_dp * t) + sin3m * 0.000289_dp
      
					calcSunEqOfCenter = c
      
		End Function


		PURE Function calcSunTrueLong(t)

		!***********************************************************************/
		!* Name:    calcSunTrueLong
		!* Type:    Function
		!* Purpose: calculate the true longitude of the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun!s true longitude in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunTrueLong
			REAL(DP), INTENT(IN) :: t
			REAL(DP) e0

			REAL(DP) l0, c, O

			l0 = calcGeomMeanLongSun(t)
			c = calcSunEqOfCenter(t)

			O = l0 + c
			calcSunTrueLong = O
      
		End Function


		PURE Function calcSunTrueAnomaly(t)
		
		!***********************************************************************/
		!* Name:    calcSunTrueAnomaly (not used by sunrise, solarnoon, sunset)
		!* Type:    Function
		!* Purpose: calculate the true anamoly of the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun!s true anamoly in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunTrueAnomaly
			REAL(DP), INTENT(IN) :: t
			REAL(DP) e0

			REAL(DP) m, c, v

			m = calcGeomMeanAnomalySun(t)
			c = calcSunEqOfCenter(t)

			v = m + c
			calcSunTrueAnomaly = v
      
		End Function


		PURE Function calcSunRadVector(t)

		!***********************************************************************/
		!* Name:    calcSunRadVector (not used by sunrise, solarnoon, sunset)
		!* Type:    Function
		!* Purpose: calculate the distance to the sun in AU
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun radius vector in AUs
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunRadVector
			REAL(DP), INTENT(IN) :: t
			REAL(DP) v, e0, R

			v = calcSunTrueAnomaly(t)
			e0 = calcEccentricityEarthOrbit(t)

			R = (1.000001018_dp * (1.0_dp - e0 * e0)) / (1 + e0 * Cos(degToRad(v)))
			calcSunRadVector = R
      
		End Function


		PURE Function calcSunApparentLong(t)

		!***********************************************************************/
		!* Name:    calcSunApparentLong (not used by sunrise, solarnoon, sunset)
		!* Type:    Function
		!* Purpose: calculate the apparent longitude of the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun!s apparent longitude in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunApparentLong
			REAL(DP), INTENT(IN) :: t
			REAL(DP) O, omega, lambda

			O = calcSunTrueLong(t)

			omega = 125.04_dp - 1934.136_dp * t
			lambda = O - 0.00569_dp - 0.00478_dp * Sin(degToRad(omega))
			calcSunApparentLong = lambda

		End Function


		PURE Function calcMeanObliquityOfEcliptic(t)

		!***********************************************************************/
		!* Name:    calcMeanObliquityOfEcliptic
		!* Type:    Function
		!* Purpose: calculate the mean obliquity of the ecliptic
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   mean obliquity in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcMeanObliquityOfEcliptic
			REAL(DP), INTENT(IN) :: t
			REAL(DP) seconds, e0

			seconds = 21.448_dp - t * (46.815_dp + t * (0.00059_dp - t * (0.001813_dp)))
			e0 = 23.0_dp + (26.0_dp + (seconds / 60.0_dp)) / 60.0_dp
			calcMeanObliquityOfEcliptic = e0
      
		End Function
  

		PURE Function calcObliquityCorrection(t)

		!***********************************************************************/
		!* Name:    calcObliquityCorrection
		!* Type:    Function
		!* Purpose: calculate the corrected obliquity of the ecliptic
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   corrected obliquity in degrees
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcObliquityCorrection
			REAL(DP), INTENT(IN) :: t
			REAL(DP) e0, omega, e1

			e0 = calcMeanObliquityOfEcliptic(t)

			omega = 125.04_dp - 1934.136_dp * t
			e1 = e0 + 0.00256_dp * Cos(degToRad(omega))
			calcObliquityCorrection = e1
      
		End Function calcObliquityCorrection
      

		PURE Function calcSunRtAscension(t)

		!***********************************************************************/
		!* Name:    calcSunRtAscension (not used by sunrise, solarnoon, sunset)
		!* Type:    Function
		!* Purpose: calculate the right ascension of the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun!s right ascension in degrees
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcSunRtAscension
		REAL(DP), INTENT(IN) :: t
		REAL(DP) e0, lambda, tananum, tanadenom, alpha

		e0 = calcObliquityCorrection(t)
		lambda = calcSunApparentLong(t)

		tananum = (Cos(degToRad(e0)) * Sin(degToRad(lambda)))
		tanadenom = (Cos(degToRad(lambda)))

		!original NOAA code using javascript Math.Atan2(y,x) convention:
		!        var alpha = radToDeg(Math.atan2(tananum, tanadenom));
		!        alpha = radToDeg(Atan2(tananum, tanadenom))

		!translated using Excel VBA Atan2(x,y) convention:
		alpha = radToDeg(Atan2(tanadenom, tananum))

		calcSunRtAscension = alpha

		End Function calcSunRtAscension


		PURE Function calcSunDeclination(t)

		!***********************************************************************/
		!* Name:    calcSunDeclination
		!* Type:    Function
		!* Purpose: calculate the declination of the sun
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   sun!s declination in degrees
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcSunDeclination
		REAL(DP), INTENT(IN) :: t
		REAL(DP) e0, lambda, sint, theta

		e0 = calcObliquityCorrection(t)
		lambda = calcSunApparentLong(t)

		sint = Sin(degToRad(e0)) * Sin(degToRad(lambda))
		theta = radToDeg(Asin(sint))
		calcSunDeclination = theta

		End Function calcSunDeclination


		PURE Function calcEquationOfTime(t)

		!***********************************************************************/
		!* Name:    calcEquationOfTime
		!* Type:    Function
		!* Purpose: calculate the difference between true solar time and mean
		!*     solar time
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!* Return value:
		!*   equation of time in minutes of time
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcEquationOfTime
		REAL(DP), INTENT(IN) :: t

		REAL(DP) epsilon, l0, e0, m
		REAL(DP) y, sin2l0, sinm
		REAL(DP) cos2l0, sin4l0, sin2m, Etime

		epsilon = calcObliquityCorrection(t)
		l0 = calcGeomMeanLongSun(t)
		e0 = calcEccentricityEarthOrbit(t)
		m = calcGeomMeanAnomalySun(t)

		y = Tan(degToRad(epsilon) / 2.0_dp)
		y = y * y

		sin2l0 = Sin(2.0_dp * degToRad(l0))
		sinm = Sin(degToRad(m))
		cos2l0 = Cos(2.0_dp * degToRad(l0))
		sin4l0 = Sin(4.0_dp * degToRad(l0))
		sin2m = Sin(2.0_dp * degToRad(m))

		Etime = y * sin2l0 - 2.0_dp * e0 * sinm + 4.0_dp * e0 * y * sinm * cos2l0 &
						- 0.5_dp * y * y * sin4l0 - 1.25_dp * e0 * e0 * sin2m

		calcEquationOfTime = radToDeg(Etime) * 4.0_dp
      
		End Function calcEquationOfTime
  
  
		PURE Function calcHourAngleSunrise(lat, SolarDec)

		!***********************************************************************/
		!* Name:    calcHourAngleSunrise
		!* Type:    Function
		!* Purpose: calculate the hour angle of the sun at sunrise for the
		!*         latitude
		!* Arguments:
		!*   lat : latitude of observer in degrees
		!* solarDec : declination angle of sun in degrees
		!* Return value:
		!*   hour angle of sunrise in radians
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcHourAngleSunrise
			REAL(DP), INTENT(IN) :: lat, SolarDec

			REAL(DP) latrad, sdRad, HAarg, HA

			latrad = degToRad(lat)
			sdRad = degToRad(SolarDec)

			HAarg = (Cos(degToRad(90.833_DP)) / (Cos(latrad) * Cos(sdRad)) - Tan(latrad) * Tan(sdRad))

			HA = (Acos(Cos(degToRad(90.833_DP)) &
						/ (Cos(latrad) * Cos(sdRad)) - Tan(latrad) * Tan(sdRad)))

			calcHourAngleSunrise = HA
      
		End Function


		Function calcHourAngleSunset(lat, SolarDec)

		!***********************************************************************/
		!* Name:    calcHourAngleSunset
		!* Type:    Function
		!* Purpose: calculate the hour angle of the sun at sunset for the
		!*         latitude
		!* Arguments:
		!*   lat : latitude of observer in degrees
		!* solarDec : declination angle of sun in degrees
		!* Return value:
		!*   hour angle of sunset in radians
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcHourAngleSunset
		REAL(DP), INTENT(IN) :: lat, SolarDec
		REAL(DP) e0

		REAL(DP) latrad, sdRad, HAarg, HA

					latrad = degToRad(lat)
					sdRad = degToRad(SolarDec)

					HAarg = (Cos(degToRad(90.833_DP)) / (Cos(latrad) * Cos(sdRad)) - Tan(latrad) * Tan(sdRad))

					HA = (Acos(Cos(degToRad(90.833_DP)) &
								 / (Cos(latrad) * Cos(sdRad)) - Tan(latrad) * Tan(sdRad)))

					calcHourAngleSunset = -HA
      
		End Function


		Function calcSunriseUTC(jd, Latitude, longitude)

		!***********************************************************************/
		!* Name:    calcSunriseUTC
		!* Type:    Function
		!* Purpose: calculate the Universal Coordinated Time (UTC) of sunrise
		!*         for the given day at the given location on earth
		!* Arguments:
		!*   JD  : julian day
		!*   latitude : latitude of observer in degrees
		!*   longitude : longitude of observer in degrees
		!* Return value:
		!*   time in minutes from zero Z
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) calcSunriseUTC
		REAL(DP), INTENT(IN) :: jd, Latitude, longitude

		REAL(DP) t, eqtime, SolarDec, hourangle
		REAL(DP) delta, timeDiff, timeUTC
		REAL(DP) newt

					t = calcTimeJulianCent(jd)

		!        // *** First pass to approximate sunrise

					eqtime = calcEquationOfTime(t)
					SolarDec = calcSunDeclination(t)
					hourangle = calcHourAngleSunrise(Latitude, SolarDec)

					delta = longitude - radToDeg(hourangle)
					timeDiff = 4.0_dp * delta
		! in minutes of time
					timeUTC = 720.0_dp + timeDiff - eqtime
		! in minutes

		! *** Second pass includes fractional jday in gamma calc

					newt = calcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC / 1440.0_dp)
					eqtime = calcEquationOfTime(newt)
					SolarDec = calcSunDeclination(newt)
					hourangle = calcHourAngleSunrise(Latitude, SolarDec)
					delta = longitude - radToDeg(hourangle)
					timeDiff = 4.0_dp * delta
					timeUTC = 720.0_dp + timeDiff - eqtime
		! in minutes

					calcSunriseUTC = timeUTC

		End Function


		PURE Function calcSolNoonUTC(t, longitude)

		!***********************************************************************/
		!* Name:    calcSolNoonUTC
		!* Type:    Function
		!* Purpose: calculate the Universal Coordinated Time (UTC) of solar
		!*     noon for the given day at the given location on earth
		!* Arguments:
		!*   t : number of Julian centuries since J2000.0
		!*   longitude : longitude of observer in degrees
		!* Return value:
		!*   time in minutes from zero Z
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSolNoonUTC
			REAL(DP), INTENT(IN) :: t, longitude
			REAL(DP) newt, eqtime, solarNoonDec, solNoonUTC
      
			newt = calcTimeJulianCent(calcJDFromJulianCent(t) + 0.5_dp + longitude / 360.0_dp)
			eqtime = calcEquationOfTime(newt)
			solarNoonDec = calcSunDeclination(newt)
			solNoonUTC = 720 + (longitude * 4) - eqtime
  
			calcSolNoonUTC = solNoonUTC

		End Function


		Function calcSunsetUTC(jd, Latitude, longitude)

		!***********************************************************************/
		!* Name:    calcSunsetUTC
		!* Type:    Function
		!* Purpose: calculate the Universal Coordinated Time (UTC) of sunset
		!*         for the given day at the given location on earth
		!* Arguments:
		!*   JD  : julian day
		!*   latitude : latitude of observer in degrees
		!*   longitude : longitude of observer in degrees
		!* Return value:
		!*   time in minutes from zero Z
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) calcSunsetUTC
			REAL(DP), INTENT(IN) :: jd, Latitude, longitude
			REAL(DP) t, eqtime, SolarDec, hourangle
			REAL(DP) delta, timeDiff, timeUTC
			REAL(DP) newt
              
					t = calcTimeJulianCent(jd)

		!        // First calculates sunrise and approx length of day

					eqtime = calcEquationOfTime(t)
					SolarDec = calcSunDeclination(t)
					hourangle = calcHourAngleSunset(Latitude, SolarDec)

					delta = longitude - radToDeg(hourangle)
					timeDiff = 4 * delta
					timeUTC = 720 + timeDiff - eqtime

		!        // first pass used to include fractional day in gamma calc

					newt = calcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC / 1440.0_dp)
					eqtime = calcEquationOfTime(newt)
					SolarDec = calcSunDeclination(newt)
					hourangle = calcHourAngleSunset(Latitude, SolarDec)

					delta = longitude - radToDeg(hourangle)
					timeDiff = 4 * delta
					timeUTC = 720 + timeDiff - eqtime
		!        // in minutes

					calcSunsetUTC = timeUTC

		End Function calcSunsetUTC


		Function sunrise(lat, lon, year, month, day, timezone, dlstime)
  
		!***********************************************************************/
		!* Name:    sunrise
		!* Type:    Main Function called by spreadsheet
		!* Purpose: calculate time of sunrise  for the entered date
		!*     and location.
		!* For latitudes greater than 72 degrees N and S, calculations are
		!* accurate to within 10 minutes. For latitudes less than +/- 72°
		!* accuracy is approximately one minute.
		!* Arguments:
		!   latitude = latitude (decimal degrees)
		!   longitude = longitude (decimal degrees)
		!    NOTE: longitude is negative for western hemisphere for input cells
		!          in the spreadsheet for calls to the functions named
		!          sunrise, solarnoon, and sunset. Those functions convert the
		!          longitude to positive for the western hemisphere for calls to
		!          other functions using the original sign convention
		!          from the NOAA javascript code.
		!   year = year
		!   month = month
		!   day = day
		!   timezone = time zone hours relative to GMT/UTC (hours)
		!   dlstime = daylight savings time (0 = no, 1 = yes) (hours)
		!* Return value:
		!*   sunrise time in local time (days)
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) sunrise
		REAL(DP), INTENT(IN) :: lat, lon, year, month, day

		!GP 23-Nov-09
		!INTEGER(I4B), INTENT(IN) :: timezone, dlstime
		REAL(DP), INTENT(IN) :: timezone, dlstime

		REAL(DP) longitude, Latitude, jd
		REAL(DP) riseTimeGMT, riseTimeLST

		! change sign convention for longitude from negative to positive in western hemisphere
							longitude = lon * -1
							Latitude = lat
							If (Latitude > 89.8_dp) Then
								Latitude = 89.8_dp
							ELSE If (Latitude < -89.8_dp) THEN
								Latitude = -89.8_dp
							END IF
							jd = calcJD(year, month, day)

		!            // Calculate sunrise for this date
							riseTimeGMT = calcSunriseUTC(jd, Latitude, longitude)
           
		!            //  adjust for time zone and daylight savings time in minutes
							riseTimeLST = riseTimeGMT + (60 * timezone) + (dlstime * 60)

		!            //  convert to days
							sunrise = riseTimeLST / 1440.0_dp

		End Function


		PURE Function solarnoon(lat, lon, year, month, day, timezone, dlstime)

		!***********************************************************************/
		!* Name:    solarnoon
		!* Type:    Main Function called by spreadsheet
		!* Purpose: calculate the Universal Coordinated Time (UTC) of solar
		!*     noon for the given day at the given location on earth
		!* Arguments:
		!    year
		!    month
		!    day
		!*   longitude : longitude of observer in degrees
		!    NOTE: longitude is negative for western hemisphere for input cells
		!          in the spreadsheet for calls to the functions named
		!          sunrise, solarnoon, and sunset. Those functions convert the
		!          longitude to positive for the western hemisphere for calls to
		!          other functions using the original sign convention
		!          from the NOAA javascript code.
		!* Return value:
		!*   time of solar noon in local time days
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE
			REAL(DP) solarnoon
			REAL(DP), INTENT(IN) :: lat, lon, year, month, day

			!GP 23-Nov-09
			!INTEGER(I4B), INTENT(IN) :: timezone, dlstime
			REAL(DP), INTENT(IN) :: timezone, dlstime
      
			REAL(DP) longitude, Latitude, jd
			REAL(DP) t, newt, eqtime
			REAL(DP) solarNoonDec, solNoonUTC
              
		! change sign convention for longitude from negative to positive in western hemisphere
			longitude = lon * -1
			Latitude = lat
			If (Latitude > 89.8) Then
				Latitude = 89.8
			ELSE If (Latitude < -89.8) Then
			  Latitude = -89.8
			END IF
			jd = calcJD(year, month, day)
			t = calcTimeJulianCent(jd)
  
			newt = calcTimeJulianCent(calcJDFromJulianCent(t) + 0.5_dp + longitude / 360.0_dp)

			eqtime = calcEquationOfTime(newt)
			solarNoonDec = calcSunDeclination(newt)
			solNoonUTC = 720 + (longitude * 4) - eqtime
  
!            //  adjust for time zone and daylight savings time in minutes
			solarnoon = solNoonUTC + (60 * timezone) + (dlstime * 60)

!            //  convert to days
			solarnoon = solarnoon / 1440.0_dp

		End Function


		Function sunset(lat, lon, year, month, day, timezone, dlstime)
  
		!***********************************************************************/
		!* Name:    sunset
		!* Type:    Main Function called by spreadsheet
		!* Purpose: calculate time of sunrise and sunset for the entered date
		!*     and location.
		!* For latitudes greater than 72 degrees N and S, calculations are
		!* accurate to within 10 minutes. For latitudes less than +/- 72°
		!* accuracy is approximately one minute.
		!* Arguments:
		!   latitude = latitude (decimal degrees)
		!   longitude = longitude (decimal degrees)
		!    NOTE: longitude is negative for western hemisphere for input cells
		!          in the spreadsheet for calls to the functions named
		!          sunrise, solarnoon, and sunset. Those functions convert the
		!          longitude to positive for the western hemisphere for calls to
		!          other functions using the original sign convention
		!          from the NOAA javascript code.
		!   year = year
		!   month = month
		!   day = day
		!   timezone = time zone hours relative to GMT/UTC (hours)
		!   dlstime = daylight savings time (0 = no, 1 = yes) (hours)
		!* Return value:
		!*   sunset time in local time (days)
		!***********************************************************************/
          
			REAL(DP) sunset
			REAL(DP), INTENT(IN) :: lat, lon, year, month, day

			!GP 23-Nov-09
			!INTEGER(I4B), INTENT(IN) :: timezone, dlstime
			REAL(DP), INTENT(IN) :: timezone, dlstime

			REAL(DP) longitude, Latitude, jd
			REAL(DP) setTimeGMT, setTimeLST

		! change sign convention for longitude from negative to positive in western hemisphere
							longitude = lon * -1
							Latitude = lat
							If (Latitude > 89.8) Then
								 Latitude = 89.8
							ELSE If (Latitude < -89.8) Then
							 Latitude = -89.8
							END IF
							jd = calcJD(year, month, day)

		!           // Calculate sunset for this date
							setTimeGMT = calcSunsetUTC(jd, Latitude, longitude)
          
		!            //  adjust for time zone and daylight savings time in minutes
							setTimeLST = setTimeGMT + (60 * timezone) + (dlstime * 60)

		!            //  convert to days
							sunset = setTimeLST / 1440.0_dp

		End Function


		Function solarazimuth(lat, lon, year, month, day, &
												hours, minutes, seconds, timezone, dlstime)

		!***********************************************************************/
		!* Name:    solarazimuth
		!* Type:    Main Function
		!* Purpose: calculate solar azimuth (deg from north) for the entered
		!*          date, time and location. Returns -999999 if darker than twilight
		!*
		!* Arguments:
		!*   latitude, longitude, year, month, day, hour, minute, second,
		!*   timezone, daylightsavingstime
		!* Return value:
		!*   solar azimuth in degrees from north
		!*
		!* Note: solarelevation and solarazimuth functions are identical
		!*       and could be converted to a VBA subroutine that would return
		!*       both values.
		!*
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) solarazimuth
		REAL(DP), INTENT(IN) ::lat, lon, year, month, day, &
													 hours, minutes, seconds, timezone, dlstime
		REAL(DP) longitude, Latitude
		REAL(DP) Zone, daySavings
		REAL(DP) hh, mm, ss, timenow
		REAL(DP) jd, t, R
		REAL(DP) alpha, theta, Etime, eqtime
		REAL(DP) SolarDec, earthRadVec, solarTimeFix
		REAL(DP) trueSolarTime, hourangle, harad
		REAL(DP) csz, zenith, azDenom, azRad
		REAL(DP) azimuth, exoatmElevation
		REAL(DP) step1, step2, step3
		REAL(DP) refractionCorrection, Te, solarzen

		! change sign convention for longitude from negative to positive in western hemisphere
							longitude = lon * -1
							Latitude = lat
							If (Latitude > 89.8) Then
								 Latitude = 89.8
							ELSE If (Latitude < -89.8) Then
							 Latitude = -89.8_dp
							END IF
          
		!change time zone to ppositive hours in western hemisphere
							Zone = timezone * -1
							daySavings = dlstime * 60
							hh = hours - (daySavings / 60.0_dp)
							mm = minutes
							ss = seconds

		!//    timenow is GMT time for calculation in hours since 0Z
							timenow = hh + mm / 60.0_dp + ss / 3600.0_dp + Zone

							jd = calcJD(year, month, day)
							t = calcTimeJulianCent(jd + timenow / 24.0)
							R = calcSunRadVector(t)
							alpha = calcSunRtAscension(t)
							theta = calcSunDeclination(t)
							Etime = calcEquationOfTime(t)

							eqtime = Etime
							SolarDec = theta !//    in degrees
							earthRadVec = R

							solarTimeFix = eqtime - 4.0 * longitude + 60.0 * Zone
							trueSolarTime = hh * 60.0 + mm + ss / 60.0 + solarTimeFix
							!//    in minutes

							Do While (trueSolarTime > 1440)
									trueSolarTime = trueSolarTime - 1440
							END DO
          
							hourangle = trueSolarTime / 4.0 - 180.0
							!//    Thanks to Louis Schwarzmayr for the next line:
							If (hourangle < -180) hourangle = hourangle + 360.0

							harad = degToRad(hourangle)

							csz = Sin(degToRad(Latitude)) * &
										Sin(degToRad(SolarDec)) + &
										Cos(degToRad(Latitude)) * &
										Cos(degToRad(SolarDec)) * Cos(harad)

							If (csz > 1.0) Then
									csz = 1.0
							ElseIf (csz < -1.0) Then
									csz = -1.0
							End If
          
							zenith = radToDeg(Acos(csz))

							azDenom = (Cos(degToRad(Latitude)) * Sin(degToRad(zenith)))
          
							If (Abs(azDenom) > 0.001) Then
									azRad = ((Sin(degToRad(Latitude)) * &
											Cos(degToRad(zenith))) - &
											Sin(degToRad(SolarDec))) / azDenom
									If (Abs(azRad) > 1.0) Then
											If (azRad < 0) Then
													azRad = -1.0
											Else
													azRad = 1.0
											End If
									End If

									azimuth = 180.0 - radToDeg(Acos(azRad))

									If (hourangle > 0.0) Then
											azimuth = -azimuth
									End If
							Else
									If (Latitude > 0.0) Then
											azimuth = 180.0
									Else
											azimuth = 0.0
									End If
							End If
							If (azimuth < 0.0) Then
									azimuth = azimuth + 360.0
							End If
                      
							exoatmElevation = 90.0 - zenith

		!beginning of complex expression commented out
		!            If (exoatmElevation > 85.0) Then
		!                refractionCorrection = 0.0
		!            Else
		!                te = Tan(degToRad(exoatmElevation))
		!                If (exoatmElevation > 5.0) Then
		!                    refractionCorrection = 58.1 / te - 0.07 / (te * te * te) + &
		!                        0.000086 / (te * te * te * te * te)
		!                ElseIf (exoatmElevation > -0.575) Then
		!                    refractionCorrection = 1735.0 + exoatmElevation * &
		!                        (-518.2 + exoatmElevation * (103.4 + &
		!                        exoatmElevation * (-12.79 + &
		!                        exoatmElevation * 0.711)))
		!                Else
		!                    refractionCorrection = -20.774 / te
		!                End If
		!                refractionCorrection = refractionCorrection / 3600.0
		!            End If
		!end of complex expression

		!beginning of simplified expression
							If (exoatmElevation > 85.0) Then
									refractionCorrection = 0.0
							Else
									Te = Tan(degToRad(exoatmElevation))
									If (exoatmElevation > 5.0_dp) Then
											refractionCorrection = 58.1 / Te - 0.07_dp / (Te * Te * Te) + &
													0.000086_dp / (Te * Te * Te * Te * Te)
									ElseIf (exoatmElevation > -0.575) Then
											step1 = (-12.79_dp + exoatmElevation * 0.711_dp)
											step2 = (103.4_dp + exoatmElevation * (step1))
											step3 = (-518.2_dp + exoatmElevation * (step2))
											refractionCorrection = 1735.0_dp + exoatmElevation * (step3)
									Else
											refractionCorrection = -20.774 / Te
									End If
									refractionCorrection = refractionCorrection / 3600.0
							End If
		!end of simplified expression
          
							solarzen = zenith - refractionCorrection
                   
		!            If (solarZen < 108.0) Then
								solarazimuth = azimuth
		!              solarelevation = 90.0 - solarZen
		!              If (solarZen < 90.0) Then
		!                coszen = Cos(degToRad(solarZen))
		!              Else
		!                coszen = 0.0
		!              End If
		!            Else    !// do not report az & el after astro twilight
		!              solarazimuth = -999999
		!              solarelevation = -999999
		!              coszen = -999999
		!            End If

		End Function


		PURE Function solarelevation(lat, lon, year, month, day, &
												hours, minutes, seconds, timezone, dlstime)

		!***********************************************************************/
		!* Name:    solarazimuth
		!* Type:    Main Function
		!* Purpose: calculate solar azimuth (deg from north) for the entered
		!*          date, time and location. Returns -999999 if darker than twilight
		!*
		!* Arguments:
		!*   latitude, longitude, year, month, day, hour, minute, second,
		!*   timezone, daylightsavingstime
		!* Return value:
		!*   solar azimuth in degrees from north
		!*
		!* Note: solarelevation and solarazimuth functions are identical
		!*       and could converted to a VBA subroutine that would return
		!*       both values.
		!*
		!***********************************************************************/
		USE nrtype
		IMPLICIT NONE
		REAL(DP) solarelevation
		REAL(DP), INTENT(IN) ::lat, lon, year, month, day, seconds

		!GP 23-Nov-09
		!INTEGER(I4B), INTENT(IN) :: hours, minutes, timezone, dlstime
		INTEGER(I4B), INTENT(IN) :: hours, minutes
		REAL(DP), INTENT(IN) :: timezone, dlstime

		REAL(DP) longitude, Latitude
		REAL(DP) Zone, daySavings
		REAL(DP) hh, mm, ss, timenow
		REAL(DP) jd, t, R
		REAL(DP) alpha, theta, Etime, eqtime
		REAL(DP) SolarDec, earthRadVec, solarTimeFix
		REAL(DP) trueSolarTime, hourangle, harad
		REAL(DP) csz, zenith, azDenom, azRad
		REAL(DP) azimuth, exoatmElevation
		REAL(DP) step1, step2, step3
		REAL(DP) refractionCorrection, Te, solarzen

		! change sign convention for longitude from negative to positive in western hemisphere
		longitude = lon * -1
		Latitude = lat
		If (Latitude > 89.8) Latitude = 89.8
		If (Latitude < -89.8) Latitude = -89.8
          
		!change time zone to ppositive hours in western hemisphere
		Zone = timezone * -1
		daySavings = dlstime * 60.0
		hh = hours - (daySavings / 60.0)
		mm = minutes
		ss = seconds

		!//    timenow is GMT time for calculation in hours since 0Z
		timenow = hh + mm / 60.0 + ss / 3600.0 + Zone

		jd = calcJD(year, month, day)
		t = calcTimeJulianCent(jd + timenow / 24.0)
		R = calcSunRadVector(t)
		alpha = calcSunRtAscension(t)
		theta = calcSunDeclination(t)
		Etime = calcEquationOfTime(t)

		eqtime = Etime
		SolarDec = theta !//    in degrees
		earthRadVec = R

		solarTimeFix = eqtime - 4.0 * longitude + 60.0 * Zone
		trueSolarTime = hh * 60.0 + mm + ss / 60.0 + solarTimeFix
		!//    in minutes

		Do While (trueSolarTime > 1440)
			trueSolarTime = trueSolarTime - 1440
		END DO
          
		hourangle = trueSolarTime / 4.0 - 180.0
			 !//    Thanks to Louis Schwarzmayr for the next line:
		If (hourangle < -180) hourangle = hourangle + 360.0

		harad = degToRad(hourangle)

		csz = Sin(degToRad(Latitude)) * &
					Sin(degToRad(SolarDec)) + &
					Cos(degToRad(Latitude)) * &
					Cos(degToRad(SolarDec)) * Cos(harad)

		If (csz > 1.0) Then
				csz = 1.0
		ElseIf (csz < -1.0) Then
				csz = -1.0
		End If

		zenith = radToDeg(Acos(csz))

		azDenom = (Cos(degToRad(Latitude)) * Sin(degToRad(zenith)))

		If (Abs(azDenom) > 0.001) Then
			azRad = ((Sin(degToRad(Latitude)) * &
									Cos(degToRad(zenith))) - &
									Sin(degToRad(SolarDec))) / azDenom
				If (Abs(azRad) > 1.0) Then
						If (azRad < 0) Then
								azRad = -1.0
						Else
								azRad = 1.0
						End If
				End If

				azimuth = 180.0 - radToDeg(Acos(azRad))

				If (hourangle > 0.0) Then
						azimuth = -azimuth
				End If
		Else
				If (Latitude > 0.0) Then
						azimuth = 180.0
				Else
						azimuth = 0.0
				End If
		End If
		If (azimuth < 0.0) Then
				azimuth = azimuth + 360.0
		End If
            
		exoatmElevation = 90.0 - zenith

		!beginning of complex expression commented out
		!            If (exoatmElevation > 85.0) Then
		!                refractionCorrection = 0.0
		!            Else
		!                te = Tan(degToRad(exoatmElevation))
		!                If (exoatmElevation > 5.0) Then
		!                    refractionCorrection = 58.1 / te - 0.07 / (te * te * te) + &
		!                        0.000086 / (te * te * te * te * te)
		!                ElseIf (exoatmElevation > -0.575) Then
		!                    refractionCorrection = 1735.0 + exoatmElevation * &
		!                        (-518.2 + exoatmElevation * (103.4 + &
		!                        exoatmElevation * (-12.79 + &
		!                        exoatmElevation * 0.711)))
		!                Else
		!                    refractionCorrection = -20.774 / te
		!                End If
		!                refractionCorrection = refractionCorrection / 3600.0
		!            End If
		!end of complex expression

		!beginning of simplified expression
							If (exoatmElevation > 85.0) Then
									refractionCorrection = 0.0
							Else
									Te = Tan(degToRad(exoatmElevation))
									If (exoatmElevation > 5.0) Then
											refractionCorrection = 58.1 / Te - 0.07 / (Te * Te * Te) + &
													0.000086 / (Te * Te * Te * Te * Te)
									ElseIf (exoatmElevation > -0.575) Then
											step1 = (-12.79 + exoatmElevation * 0.711)
											step2 = (103.4 + exoatmElevation * (step1))
											step3 = (-518.2 + exoatmElevation * (step2))
											refractionCorrection = 1735.0 + exoatmElevation * (step3)
									Else
											refractionCorrection = -20.774 / Te
									End If
									refractionCorrection = refractionCorrection / 3600.0
							End If
		!end of simplified expression
          
							solarzen = zenith - refractionCorrection
                   
		!            If (solarZen < 108.0) Then
		!              solarazimuth = azimuth
								solarelevation = 90.0 - solarzen
		!              If (solarZen < 90.0) Then
		!                coszen = Cos(degToRad(solarZen))
		!              Else
		!                coszen = 0.0
		!              End If
		!            Else    !// do not report az & el after astro twilight
		!              solarazimuth = -999999
		!              solarelevation = -999999
		!              coszen = -999999
		!            End If

		End Function solarelevation


		SUBROUTINE solarposition(lat, lon, year, month, day, &
					hours, minutes, seconds, timezone, dlstime, &
					solarazimuth, solarelevation, earthRadVec)

		!***********************************************************************/
		!* Name:    solarposition
		!* Type:    Subroutine
		!* Purpose: calculate solar azimuth (deg from north)
		!*          and elevation (deg from horizeon) for the entered
		!*          date, time and location.
		!*
		!* Arguments:
		!*   latitude, longitude, year, month, day, hour, minute, second,
		!*   timezone, daylightsavingstime
		!* Return value:
		!*   solar azimuth in degrees from north
		!*   solar elevation in degrees from horizon
		!*   earth radius vector (distance to the sun in AU)
		!*
		!* Note: solarelevation and solarazimuth functions are identical
		!*       and could converted to a VBA subroutine that would return
		!*       both values.
		!*
		!***********************************************************************/
			USE nrtype
			IMPLICIT NONE

			REAL(DP) lat, lon, year, month, day

			!GP 23-Nov-09
			!INTEGER(I4B)	hours, minutes, timezone, dlstime
			INTEGER(I4B)	hours, minutes
			REAL(DP)	timezone, dlstime

			REAL(DP) seconds, solarazimuth, solarelevation, earthRadVec
			REAL(DP) longitude, Latitude
			REAL(DP) Zone, daySavings
			REAL(DP) hh, mm, ss, timenow
			REAL(DP) jd, t, R
			REAL(DP) alpha, theta, Etime, eqtime
			!REAL(DP) SolarDec, earthRadVec, solarTimeFix
			REAL(DP) SolarDec, solarTimeFix
			REAL(DP) trueSolarTime, hourangle, harad
			REAL(DP) csz, zenith, azDenom, azRad
			REAL(DP) azimuth, exoatmElevation
			REAL(DP) step1, step2, step3
			REAL(DP) refractionCorrection, Te, solarzen

			! change sign convention for longitude from negative to positive in western hemisphere
							longitude = lon * -1
							Latitude = lat
							If (Latitude > 89.8) Then
								 Latitude = 89.8
							ELSE If (Latitude < -89.8) Then
							 Latitude = -89.8
							END IF

          
			!change time zone to ppositive hours in western hemisphere
							Zone = timezone * -1
							daySavings = dlstime * 60
							hh = hours - (daySavings / 60)
							mm = minutes
							ss = seconds

		!//    timenow is GMT time for calculation in hours since 0Z
							timenow = hh + mm / 60 + ss / 3600 + Zone

							jd = calcJD(year, month, day)
							t = calcTimeJulianCent(jd + timenow / 24.0)
							R = calcSunRadVector(t)
							alpha = calcSunRtAscension(t)
							theta = calcSunDeclination(t)
							Etime = calcEquationOfTime(t)

							eqtime = Etime
							SolarDec = theta !//    in degrees
							earthRadVec = R

							solarTimeFix = eqtime - 4.0 * longitude + 60.0 * Zone
							trueSolarTime = hh * 60.0 + mm + ss / 60.0 + solarTimeFix
							!//    in minutes

							Do While (trueSolarTime > 1440)
									trueSolarTime = trueSolarTime - 1440
							END DO
          
							hourangle = trueSolarTime / 4.0 - 180.0
							!//    Thanks to Louis Schwarzmayr for the next line:
							If (hourangle < -180) hourangle = hourangle + 360.0

							harad = degToRad(hourangle)

							csz = Sin(degToRad(Latitude)) * &
										Sin(degToRad(SolarDec)) + &
										Cos(degToRad(Latitude)) * &
										Cos(degToRad(SolarDec)) * Cos(harad)

							If (csz > 1.0) Then
									csz = 1.0
							ElseIf (csz < -1.0) Then
									csz = -1.0
							End If
          
							zenith = radToDeg(Acos(csz))

							azDenom = (Cos(degToRad(Latitude)) * Sin(degToRad(zenith)))
          
							If (Abs(azDenom) > 0.001) Then
									azRad = ((Sin(degToRad(Latitude)) * &
											Cos(degToRad(zenith))) - &
											Sin(degToRad(SolarDec))) / azDenom
									If (Abs(azRad) > 1.0) Then
											If (azRad < 0) Then
													azRad = -1.0
											Else
													azRad = 1.0
											End If
									End If

									azimuth = 180.0 - radToDeg(Acos(azRad))

									If (hourangle > 0.0) Then
											azimuth = -azimuth
									End If
							Else
									If (Latitude > 0.0) Then
											azimuth = 180.0
									Else
											azimuth = 0.0
									End If
							End If
							If (azimuth < 0.0) Then
									azimuth = azimuth + 360.0
							End If
                      
							exoatmElevation = 90.0 - zenith

		!beginning of complex expression commented out
		!            If (exoatmElevation > 85.0) Then
		!                refractionCorrection = 0.0
		!            Else
		!                te = Tan(degToRad(exoatmElevation))
		!                If (exoatmElevation > 5.0) Then
		!                    refractionCorrection = 58.1 / te - 0.07 / (te * te * te) + &
		!                        0.000086 / (te * te * te * te * te)
		!                ElseIf (exoatmElevation > -0.575) Then
		!                    refractionCorrection = 1735.0 + exoatmElevation * &
		!                        (-518.2 + exoatmElevation * (103.4 + &
		!                        exoatmElevation * (-12.79 + &
		!                        exoatmElevation * 0.711)))
		!                Else
		!                    refractionCorrection = -20.774 / te
		!                End If
		!                refractionCorrection = refractionCorrection / 3600.0
		!            End If
		!end of complex expression


		!beginning of simplified expression
							If (exoatmElevation > 85.0) Then
									refractionCorrection = 0.0
							Else
									Te = Tan(degToRad(exoatmElevation))
									If (exoatmElevation > 5.0) Then
											refractionCorrection = 58.1 / Te - 0.07 / (Te * Te * Te) + &
													0.000086 / (Te * Te * Te * Te * Te)
									ElseIf (exoatmElevation > -0.575) Then
											step1 = (-12.79 + exoatmElevation * 0.711)
											step2 = (103.4 + exoatmElevation * (step1))
											step3 = (-518.2 + exoatmElevation * (step2))
											refractionCorrection = 1735.0 + exoatmElevation * (step3)
									Else
											refractionCorrection = -20.774 / Te
									End If
									refractionCorrection = refractionCorrection / 3600.0
							End If
		!end of simplified expression
          
          
							solarzen = zenith - refractionCorrection
                   
		!            If (solarZen < 108.0) Then
								solarazimuth = azimuth
								solarelevation = 90.0 - solarzen
		!              If (solarZen < 90.0) Then
		!                coszen = Cos(degToRad(solarZen))
		!              Else
		!                coszen = 0.0
		!              End If
		!            Else    !// do not report az & el after astro twilight
		!              solarazimuth = -999999
		!              solarelevation = -999999
		!              coszen = -999999
		!            End If

		End SUBROUTINE

END MODULE Class_SolarPosition







