!Handle meteology data interpolation
MODULE Class_Meteo
	USE nrtype
	USE Class_RiverTopo
	USE Class_SystemParams, ONLY: steadystate
	IMPLICIT NONE

	PRIVATE !all variables and functions are private, unless declare public 
					!only the following subroutines are public 
	PUBLIC MeteoData_, InstanteousMeteo, Meteo_type
	
	TYPE Meteo_type

		!gp 16-Jul-08
		!REAL(DP), DIMENSION(:,:), POINTER :: shadeHH, TaHH, TdHH, UwHH, ccHH
		REAL(DP), DIMENSION(:,:), POINTER :: shadeHH, TaHH, TdHH, UwHH, ccHH, solarHH

		!instantenous Meteorology data

		!gp 16-Jul-08
		!REAL(DP), DIMENSION(:),   POINTER :: shadeT, Ta, Td, Uw, cc
		REAL(DP), DIMENSION(:),   POINTER :: shadeT, Ta, Td, Uw, cc, solarT

		INTEGER(I4B) :: NumdatPnt= 24		!use to identify the size of time series
																		!hardwire to (0-23 hour) for now
	END TYPE Meteo_type

	CONTAINS
	
	!gp 16-Jul-08
	!FUNCTION MeteoData_(nr, shadeHHin, TaHHin, TdHHin, UwHHin, ccHHin) RESULT(Met)
	FUNCTION MeteoData_(nr, shadeHHin, TaHHin, TdHHin, UwHHin, ccHHin, solarHHin) RESULT(Met)
	
	TYPE(Meteo_type) Met
	INTEGER(I4B), INTENT(IN) :: nr

	!gp 16-Jul-08
	!REAL(DP), DIMENSION(0:,:), INTENT(IN) :: shadeHHin, TaHHin, TdHHin, UwHHin, ccHHin
	REAL(DP), DIMENSION(0:,:), INTENT(IN) :: shadeHHin, TaHHin, TdHHin, UwHHin, ccHHin, solarHHin

	!gp 16-Jul-08
	!INTEGER(I4B) status(10), i	 
	INTEGER(I4B) status(12), i	 

	!IF (steadystate()) THEN
			
	!END IF

	IF (nr > 0) THEN
		ALLOCATE (Met%shadeHH(0:Met%NumdatPnt-1 ,nr), STAT=status(1))
		ALLOCATE (Met%TaHH(0:Met%NumdatPnt-1,nr), STAT=status(2))
		ALLOCATE (Met%TdHH(0:Met%NumdatPnt-1,nr), STAT=status(3))
		ALLOCATE (Met%UwHH(0:Met%NumdatPnt-1,nr), STAT=status(4))
		ALLOCATE (Met%ccHH(0:Met%NumdatPnt-1,nr), STAT=status(5))

		ALLOCATE (Met%shadeT(nr), STAT=status(6))
		ALLOCATE (Met%Ta(nr), STAT=status(7))
		ALLOCATE (Met%Td(nr), STAT=status(8))
		ALLOCATE (Met%Uw(nr), STAT=status(9))
		ALLOCATE (Met%cc(nr), STAT=status(10))

		!gp 16-Jul-08
		ALLOCATE (Met%solarHH(0:Met%NumdatPnt-1,nr), STAT=status(11))
		ALLOCATE (Met%solarT(nr), STAT=status(12))

		DO i=1, 10
			IF (status(i)==1) STOP 'ERROR: Class_Meteo:AllocateMeteoDataArray.Allocation Failed'
		END DO

		Met%shadeHH = shadeHHin
		Met%TaHH = TaHHin
		Met%TdHH = TdHHin
		Met%UwHH = UwHHin
		Met%ccHH = ccHHin

		!gp 16-Jul-08
		Met%solarHH = solarHHin

	ELSE
		PRINT *, 'ERROR:element number must be great than 0'
		STOP 'Class_Meteo:AllocateMeteoDataArray failed'
	END IF

	END FUNCTION MeteoData_


!interpolate Hourly Meteology data
	SUBROUTINE InstanteousMeteo(nr, t, Met)

		INTEGER(I4B), INTENT(IN) :: nr
		TYPE(Meteo_type), INTENT(INOUT) :: Met
		REAL(DP), INTENT(IN) :: t

		CALL InterpolateHelper(nr, t, Met%Ta, Met%TaHH, Met%NumdatPnt)
		CALL InterpolateHelper(nr, t, Met%Td, Met%TdHH, Met%NumdatPnt)		
		CALL InterpolateHelper(nr, t, Met%Uw, Met%UwHH, Met%NumdatPnt)
		CALL InterpolateHelper(nr, t, Met%cc, Met%ccHH, Met%NumdatPnt)

		!gp 16-Jul-08
		CALL InterpolateHelper(nr, t, Met%solarT, Met%solarHH, Met%NumdatPnt)

		!shade use different method
		CALL InterpolateShade(nr, t, Met%shadeT, Met%shadeHH, Met%NumdatPnt)		
	END SUBROUTINE InstanteousMeteo


	!private subroutine
	SUBROUTINE InterpolateHelper(nr, t, Interp, c, upbound) 
	!for interpolation of instantaneous values from hourly point estimates of 2D array
	!input arrays have dimensions of (ihour, k) where ihour is hour (0 to 23)
	!and k is variable number (1 to 15) or reach number (1-1000)
	INTEGER(I4B), INTENT(IN) :: nr
	REAL(DP), INTENT(OUT) :: Interp(:)
	REAL(DP), INTENT(IN) :: t, c (0:,:)
	INTEGER(I4B), INTENT(IN):: upBound
	REAL(DP) t_hr
	INTEGER(I4B) t0_hr, t1_hr, k

	t_hr = (t - Int(t)) * (upBound)

	If (t_hr < upBound-1) Then
	 t0_hr = Int(t_hr)
	 t1_hr = t0_hr + 1
	Else
	 t0_hr = upBound-1
	 t1_hr = 0
	End If

	DO k=1, nr
		Interp(k) = c(t0_hr, k) + (t_hr - t0_hr) * (c(t1_hr, k) - c(t0_hr, k))
	END DO

	End SUBROUTINE


	SUBROUTINE InterpolateShade(nr, t, hourlyShade, shadeHH, upbound)
	!input values of hourly shade are integrated values that apply to the entire hour
	!therefore this function assigns the constant hourly values for each hourly period of the simulation

	!input arrays have dimensions of (ihour, k) where ihour is hour (0 to 23)
	!and k is variable number (1 to 15) or reach number (1-1000)
	INTEGER(I4B), INTENT(IN) :: nr
	REAL(DP), INTENT(OUT) :: HourlyShade(:)
	REAL(DP), INTENT(IN):: t, shadeHH(0:,:)
	INTEGER(I4B), INTENT(IN):: upbound
	REAL(DP) t_hr
	INTEGER(I4B) t0_hr, t1_hr, k

	t_hr = (t - Int(t)) * (upBound)

	If (t_hr < upBound-1) Then
	 t0_hr = Int(t_hr)
	 t1_hr = t0_hr + 1
	Else
	 t0_hr = upBound-1
	 t1_hr = 0
	End If
	
	DO k=1, nr
		hourlyShade(k) = shadeHH(t0_hr, k)
	END DO
	End SUBROUTINE InterpolateShade

END MODULE Class_Meteo