MODULE Class_Date
	USE nrtype

	!Derived data types
	TYPE SimDate																!11/16/04
		REAL(DP) year
		REAL(DP) month
		REAL(DP) day
		REAL(DP) Julday															!Julian day
	END TYPE

	CONTAINS

	!Data structure constructor
	PURE FUNCTION Date_(year, month, day)

		REAL(DP), INTENT(IN) :: year, month, day
		TYPE(SimDate) Date_														!11/16/04

		Date_%year = year
		Date_%month= month
		Date_%day = day
		Date_%Julday = Julcvt(INT(month, I2B), INT(day, I2B), INT(year, I2B))
	END FUNCTION Date_

!Calculate Julian day, Jan 1st =1, 
	PURE FUNCTION Julcvt(month, day, year)
	
		REAL(DP) Julcvt
		INTEGER(I4B), INTENT(IN):: month, day, year
		INTEGER(I4B) leap
		Integer(I4B) modtest

!gp edit for Absoft syntax
!		IF (MOD(year, 4) == 0) leap = 1
		IF (MOD(year, 4) .eq. 0) leap = 1
!!		IF MOD(year, 4) = 0 then
!!                modtest = MOD(year,4)
!		modtest = year - int(year/4)*4
!!		IF (modtest == 0) then
!		IF (modtest .eq. 0) then
!		 leap = 1
!		end if

		SELECT CASE (month)
			CASE (1)
				Julcvt = 0
			CASE (2)
				Julcvt = 31
			CASE (3)
				Julcvt = 59 + leap
			CASE (4)
				Julcvt = 90 + leap
			CASE (5)
				Julcvt = 120 + leap
			CASE (6)
				Julcvt = 151 + leap
			CASE (7)
				Julcvt = 181 + leap
			CASE (8)
				Julcvt = 212 + leap
			CASE (9)
				Julcvt = 243 + leap
			CASE (10)
				Julcvt = 273 + leap
			CASE (11)
				Julcvt = 304 + leap
			CASE (12)
				Julcvt = 334 + leap
		END SELECT
		Julcvt = Julcvt + day
	END FUNCTION

END MODULE


