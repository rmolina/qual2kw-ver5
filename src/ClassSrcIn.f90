!sourceinbackup.f90

MODULE Class_SourceIn
	USE nrtype
	USE Class_RiverTopo					!only num of reach and num of element
	USE Class_WaterQuality
	USE Class_Hydraulics
	IMPLICIT NONE
	
	PRIVATE PointIn_, NonPointIn_
!	PUBLIC :: Qpt, Qpta 

	!derived type for diffusion load/source
	TYPE Diffusion_type					!contain the raw data for diffusin load and source
		CHARACTER(LEN=30) name
		REAL(DP) xdup, xddn
		REAL(DP) :: Q =0.0				!if load, then >0; source, then <0
		REAL(DP) Te
		REAL(DP) c(nv-1)
		REAL(DP) pH
		INTEGER(I4B) beginRch			!beginning reach number
	END TYPE Diffusion_type	
	!derived type for point load/source
	TYPE Point_type							!
		CHARACTER(LEN=30) name
		REAL(DP) x
		REAL(DP) :: Q	=0.0				!if load, then >0; source, then <0
		REAL(DP) :: TeMean=0, TeAmp=0, TeMaxTime =0
		REAL(DP) cMean(nv-2), cAmp(nv-2), cMaxTime(nv-2)			
		REAL(DP) pHMean, pHAmp, phMaxTime
		INTEGER(I4B) beginRch			!beginning reach number
	END TYPE Point_type
	
	INTEGER(I4B) npt, ndiff

	!contains the original point and diffusion data
	TYPE(Point_type), POINTER :: point(:)
	TYPE(Diffusion_type), POINTER :: diffu(:)
	TYPE(WaterQuality_type), ALLOCATABLE :: load(:)				!combined load from both point and diffusion
	REAL(DP), ALLOCATABLE :: HeatDiff(:), loadDiff(:,:) 	!diffusion load not vary by time
	REAL(DP), ALLOCATABLE:: Qpta(:), Qpt (:) 

	CONTAINS

	! Read in Point and diffusion load and abastraction data
	SUBROUTINE sourceIn_(nr, nptIn, ndiffIn, flag, Topo, hydrau, PtName, xptt, Qptta, Qptt, TepttMean, &
										 TepttAmp, TepttMaxTime, cpttMean, cpttAmp, cpttMaxTime, phpttMean, &
										 phpttAmp, phpttMaxTime, DiffName, xdup, xddn, Qdifa, Qdif, Tedif, cdif, pHind)

		INTEGER(I4B), INTENT(IN) :: nr, nptIn, ndiffIn, flag
		TYPE(RiverTopo_type), INTENT(IN) :: Topo									!river topology
		TYPE(RiverHydraulics_type) hydrau					!channel dimensions, hydraulics, physical characters
		REAL(DP) xptt(:), Qptta(:), Qptt(:), TepttMean(:), TepttAmp(:), TepttMaxTime(:)
		REAL(DP) cpttMean(:,:), cpttAmp(:,:), cpttMaxTime(:,:)
		REAL(DP) phpttMean(:), phpttAmp(:), phpttMaxTime(:)
		CHARACTER(LEN=30), INTENT(IN) :: PtName(:)
		REAL(DP), INTENT(IN) :: xdup(:), xddn(:), Qdifa(:), Qdif(:), pHind(:), Tedif(:), cdif(:,:)
		CHARACTER(LEN=30), INTENT(IN) :: DiffName(:)
		INTEGER(I4B) status(5), i
		
		ALLOCATE(Qpta(nr), STAT=status(1))
		ALLOCATE(Qpt (nr), STAT=status(2))
		ALLOCATE(load(nr), STAT=status(3))
		ALLOCATE(HeatDiff(nr), STAT=status(4))
		ALLOCATE(loadDiff(nr,nv), STAT=status(5))
						
		DO i=1, 5
			IF (status(i)==1) THEN
				STOP 'ERROR:Class_SourceIn(sourceIN) Insufficient memory for dyanmic allocation'
			END IF
		END DO
		Qpta = 0
		Qpt  = 0
	
		CALL PointIn_(nr, nptIn, topo, flag, PtName, xptt, Qptta, Qptt, TepttMean, TepttAmp, TepttMaxTime, &
										   cpttMean, cpttAmp, cpttMaxTime, phpttMean, phpttAmp, phpttMaxTime)
		CALL NonPointIn_(nr, topo, flag, ndiffIn, DiffName, xdup, xddn, Qdifa, Qdif, Tedif, cdif, pHind)
		!generate average reach flows for hydraulics (m^3/s)
		DO i = 1, nr
			If (Qpt(i) < 0) Qpt(i) = 0
			hydrau%reach(i)%Q = hydrau%reach(i-1)%Q + Qpt(i) - Qpta(i)
			hydrau%reach(i)%Qpt=Qpt(i)
			hydrau%reach(i)%Qpta = Qpta(i)
		END DO
	
	END SUBROUTINE sourceIn_
	

	SUBROUTINE PointIn_(nr, nptIn, topo, flag, PtName, xptt, Qptta, Qptt, TepttMean, TepttAmp, TepttMaxTime, &
										   cpttMean, cpttAmp, cpttMaxTime, phpttMean, phpttAmp, phpttMaxTime)
	
		INTEGER(I4B), INTENT(IN) :: nr, nptIn, flag
		TYPE(RiverTopo_type), INTENT(IN) :: Topo									!river topology
		REAL(DP) xptt(:), Qptta(:), Qptt(:), TepttMean(:), TepttAmp(:), TepttMaxTime(:)
		REAL(DP) cpttMean(:,:), cpttAmp(:,:), cpttMaxTime(:,:)
		REAL(DP) phpttMean(:), phpttAmp(:), phpttMaxTime(:)
		CHARACTER(LEN=30), INTENT(IN) :: PtName(:)
		
		INTEGER(I4B) i, j, status
		LOGICAL(2) cond1
		npt=nptIn

		IF (npt>0) THEN

			ALLOCATE(point(npt), STAT=status)
			IF (status==1) THEN
				STOP 'Class_SourceIn:PointIn_, dynamic allocation failed!'
			END IF
			DO i=1, npt
				point(i)%name = PtName(i)
				point(i)%x = xptt(i) 
				IF (Qptta(i)>0) THEN
					point(i)%Q = -Qptta(i)					!if load, then >0; source, then <0
				ELSE
					point(i)%Q = Qptt(i)						!if load, then >0; source, then <0
				END IF
				point(i)%TeMean   = TepttMean(i)
				point(i)%TeAmp = TepttAmp(i)
				point(i)%TeMaxTime = TepttMaxTime(i)

				DO j=1, nv-2
					point(i)%cMean(j) =cpttMean(i,j)
					point(i)%cAmp(j) =cpttAmp(i,j)
					point(i)%cMaxTime(j) = cpttMaxTime(i,j)
				END DO
				point(i)%phMean = phpttMean(i)
				point(i)%phAmp = phpttAmp(i)
				point(i)%phMaxTime = phpttMaxTime(i)
				!distribute point flows to elements for hydraulics
				DO j=1, nr
					IF (flag == 1) THEN
						cond1 = (xptt(i) >= topo%reach(j-1)%xrdn) .AND. &
										(xptt(i) < topo%reach(j)%xrdn)
					ELSE
						cond1 = (xptt(i) <= topo%reach(j-1)%xrdn) .AND. &
										(xptt(i) > topo%reach(j)%xrdn)
					END IF
					IF (cond1) THEN
						IF (point(i)%Q < 0) THEN
							Qpta(j) = Qpta(j) -  point(i)%Q
						ELSE
							Qpt(j) = Qpt(j) + point(i)%Q
						END IF
						point(i)%beginRch =j
						Exit
					END IF
				END DO 
			END DO
		END IF

	END SUBROUTINE PointIn_


	!Diffusion -- nonpointer source

	SUBROUTINE NonPointIn_(nr, topo, flag, ndiffIn, DiffName, xdup, xddn, Qdifa, Qdif, Tedif, cdif, pHind)
	USE Class_Phsolve, ONLY: cT

	TYPE(RiverTopo_type) Topo									!river topology
	INTEGER(I4B), INTENT(IN) ::	ndiffIn, nr, flag
	REAL(DP), INTENT(IN) :: xdup(:), xddn(:), Qdifa(:), Qdif(:), pHind(:), Tedif(:), cdif(:,:)
	CHARACTER(LEN=30), INTENT(IN) :: DiffName(:)
	INTEGER(I4B) i, j, k, status
	LOGICAL(2) cond1, cond2, cond3, cond4, cond5
	REAL(DP) Qd, Lend 
	ndiff=ndiffIn							!number of diffusion source and abstraction

	HeatDiff=0
	LoadDiff=0

	IF (ndiff>0) THEN

		ALLOCATE(diffu(ndiff), STAT=status)
		IF (status==1) THEN
			STOP 'Class_SourceIn:PointIn_, dynamic allocation failed!'
		END IF

		DO i=1, ndiff
			diffu(i)%name = DiffName(i)
			diffu(i)%xdup = xdup(i)
			diffu(i)%xddn = xddn(i)
			IF (Qdifa(i)>0) THEN
				diffu(i)%Q	= -Qdifa(i)					!if load, then >0; source, then <0
			ELSE
				diffu(i)%Q	= Qdif(i)
			END IF
		  diffu(i)%Te   = Tedif(i)
			DO j=1, nv-2
				diffu(i)%c(j) =cdif(i,j)
			END DO

			IF (diffu(i)%c(nv - 2) == 0) diffu(i)%c(nv - 2) = 100.0
			IF (pHind(i)==0) THEN
				diffu(i)%pH = 7.0_DP
			ELSE
				diffu(i)%pH = pHind(i)
			END IF
			!Total inorganic carbon
			diffu(i)%c(nv - 1) = cT(diffu(i)%pH, diffu(i)%c(nv - 2), Tedif(i), diffu(i)%c(1))
			!distribute nonpoint flows to elements for hydraulics

			Qd = diffu(i)%Q / (xddn(i) - xdup(i)) * flag
			
			DO j=1, nr
				If (flag == 1) THEN
					cond1 = topo%reach(j)%xrdn < xdup(i) .OR. topo%reach(j-1)%xrdn > xddn(i)
					cond2 = topo%reach(j-1)%xrdn <= xdup(i) .AND. xddn(i) <= topo%reach(j)%xrdn
					cond3 = xdup(i) <= topo%reach(j-1)%xrdn .AND. xddn(i) >= topo%reach(j)%xrdn
					cond4 = topo%reach(j-1)%xrdn >= xdup(i) .AND. topo%reach(j)%xrdn >= xddn(i)
					cond5 = topo%reach(j)%xrdn >= xdup(i) .AND. topo%reach(j)%xrdn <= xddn(i)
				Else
					cond1 = topo%reach(j)%xrdn > xdup(i) .OR. topo%reach(j-1)%xrdn < xddn(i)
					cond2 = topo%reach(j-1)%xrdn >= xdup(i) .AND. xddn(i) >= topo%reach(j)%xrdn
					cond3 = xdup(i) >= topo%reach(j-1)%xrdn .AND. xddn(i) <= topo%reach(j)%xrdn
					cond4 = topo%reach(j-1)%xrdn <= xdup(i) .AND. topo%reach(j)%xrdn <= xddn(i)
					cond5 = topo%reach(j)%xrdn <= xdup(i) .AND. topo%reach(j)%xrdn >= xddn(i)
				End If
				IF (cond1) THEN
					Lend = 0
				ElseIf (cond2) THEN
					Lend = flag * (xddn(i) - xdup(i))
				ElseIf (cond3) THEN
					Lend = flag * (topo%reach(j)%xrdn - topo%reach(j-1)%xrdn)
				ElseIf (cond4) THEN
					Lend = flag * (xddn(i) - topo%reach(j-1)%xrdn)
				ElseIf (cond5) THEN
					Lend = flag * (topo%reach(j)%xrdn - xdup(i))
				End If
				If (Qd <= 0) THEN
					Qpta(j) = Qpta(j) - Lend * Qd
				Else
					Qpt(j) = Qpt(j) + Lend * Qd
				End If

				IF (Lend>0 .AND. Qd>0) THEN
!					Qpt(i) = Qpt(i) + Lend * Qd
					HeatDiff(j) = HeatDiff(j) + Lend * Qd * diffu(i)%Te
					DO k = 1, nv - 1
						LoadDiff(j, k) = LoadDiff(j, k) + Lend * Qd * diffu(i)%c(k)
					END DO
				END IF

			END DO 
		END DO
	END IF
  !distribute point flows to elements for hydraulics
!  DO i = 1, nr
!    Qpt(i) = 0
!    Qpta(i) = 0
!  END DO

	END SUBROUTINE NonPointIn_


!Calculate instanteneous sources for time t
	SUBROUTINE SourcesCalc(t, nr, flag)
	USE Class_Phsolve, ONLY: cT
	!gp new sub to evaluate point source sine functions and distribute loads to reaches at time t

	REAL(DP), INTENT(IN) :: t
	INTEGER(I4B), INTENT(IN) :: nr, flag
!	TYPE(RiverTopo_type) Topo									!river topology
	INTEGER(I4B) i, j, k, kk
	REAL(DP) Teptt(nr)
	REAL(DP) :: Heat(nr)
	REAL(DP) :: Loadi(nr, nv)
	LOGICAL(2) cond1, cond2, cond3, cond4, cond5 
	TYPE(WaterQuality_type) :: ptt	
	REAL(DP) Qd, Lend

	Heat = 0;	Loadi =0

	IF (npt > 0) THEN
		!gp evaluate the point source diel sine functions for the current time step
		DO i = 1, npt

			IF (point(i)%Q<=0) CONTINUE ! no need to process for Abastraction
			 
			ptt%Te = sinday(t, point(i)%TeMean, point(i)%TeAmp, point(i)%TeMaxTime)
			IF (ptt%Te < 0) ptt%Te = 0
			DO j = 1, nv - 2
				ptt%c(j) = sinday(t, point(i)%cMean(j), point(i)%cAmp(j), point(i)%cMaxTime(j))
				IF (ptt%c(j) < 0) ptt%c(j) = 0
			END DO
			IF (ptt%c(nv - 2) == 0) ptt%c(nv - 2) = 100
			ptt%pH = sinday(t, point(i)%pHMean, point(i)%pHAmp, point(i)%pHMaxTime)
			IF (ptt%pH < 0.01) THEN
				ptt%pH = 0.01_DP
			ELSEIF (ptt%pH > 13.99) THEN
			 ptt%pH = 13.99_DP
			END IF

			IF (point(i)%phMean == 0) ptt%pH = 7.0_DP
			ptt%c(nv - 1) = cT(ptt%pH, ptt%c(nv - 2), ptt%Te, ptt%c(1))

			j=point(i)%beginRch
			IF (point(i)%Q <= 0) THEN	!abstraction	
!						Qpta(j) = Qpta(j) - point(i)%Q
			ELSE												!load
!				Qpt(j) = Qpt(j) + point(i)%Q
				Heat(j) = Heat(j) + point(i)%Q * rhow * cpw * ptt%Te
				DO k = 1, nv - 1
					Loadi(j, k) = Loadi(j, k) + point(i)%Q * ptt%c(k)
				END DO
			END IF
		END DO

	ELSE		!no point sources
	END IF

	!generate average reach input temperatures and concentrations
	loadi=loadi+loadDiff
	Heat=Heat+HeatDiff
	DO i = 1, nr
		IF (Qpt(i) > 0) THEN
			load(i)%Te = Heat(i) / Qpt(i)
			DO j = 1, nv - 1
				load(i)%c(j) = Loadi(i, j) / Qpt(i)
			END DO
		ELSE
			Qpt(i) = 0
			load(i)%Te = 0
			DO j = 1, nv - 1
				load(i)%c(j) = 0
			END DO
		END IF
!		Q(i) = Q(i - 1) + Qpt(i) - Qpta(i)
	END DO
  
	END SUBROUTINE


	PURE FUNCTION sinday(t, xMean, xAmp, xMaxTime)
	!gp new function sinday to calculate a constituent
	!   at a particular time of day given daily mean, amplitude=(max-min)/2, and time of max (days)
		REAL(DP) sinday
		REAL(DP), INTENT(IN) :: t, xMean, xAmp, xMaxTime
		sinday = xMean + xAmp * COS(2.0_DP * PII * (t - xMaxTime))

	END FUNCTION


	PURE FUNCTION sinday2(t, xMin, xMax, xMaxTime)
	!gp new function sinday2 to calculate a constituent
	!   at a particular time of day given daily min, max, and time of max
		REAL(DP) sinday2
		REAL(DP), INTENT(IN) :: t, xMin, xMax, xMaxTime
		REAL(DP) xMean, xAmp, xTheta, xOmega

		xMean = (xMax + xMin) / 2.0_DP
		xAmp = (xMax - xMin) / 2.0_DP
		xTheta = (xMaxTime - 0.25_DP) * 2.0_DP * PII
		xOmega = 2.0_DP * PII
		sinday2 = xMean + xAmp * SIN(xOmega * t - xTheta)

	END FUNCTION

END MODULE Class_SourceIn
