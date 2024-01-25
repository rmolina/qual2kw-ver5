! classrivertopo.f90
! data structure for reach and elements
MODULE Class_RiverTopo
	USE nrtype
!	USE Class_Element
	USE Class_Reach
	IMPLICIT NONE

!	PRIVATE															!unless declared public
!	PUBLIC :: ne, nr, nHw, River_

	!/*module variables */
	TYPE RiverTopo_type
		INTEGER(I4B) nr									!number of reaches
!		INTEGER(I4B) ne									!number of elements
		INTEGER(I4B) :: nHw =1					!Hardwired to 1 for mainstem only
!		TYPE(Element_type), POINTER :: elems(:) !element array, containing topology information
		TYPE(Reach_type), DIMENSION(:),  POINTER :: reach !reaches array, 

		!gp 17-Nov-04
		CHARACTER(LEN=30) geoMethod		!gp 17-Nov-04 Depth or Width for col T and U of 'Reach' sheet

	END TYPE RiverTopo_type

!	TYPE(RiverTopo_type) Topo
	!steady state data types
	
	CONTAINS
	!/* public functions */

	!PURE FUNCTION ne()
	!	INTEGER(I4B) ne
	!	ne = Topo%ne
	!END FUNCTION ne

	!PURE FUNCTION nr()
	!	INTEGER(I4B) nr
	!	nr = Topo%nr
	!END FUNCTION nr

	!PURE FUNCTION nHw()
	!	INTEGER(I4B) nHw
	!	nHw= Topo%nHw
	!END FUNCTION nHw
 
	!data structure constructor	
	!gp 17-Nov-04 FUNCTION RiverTopo_(nRch, rlab2, rname, xrdn) RESULT(Topo)
	FUNCTION RiverTopo_(nRch, rlab2, rname, xrdn, geoMethod) RESULT(Topo)		!gp 17-Nov-04

		TYPE(RiverTopo_type) Topo	
		INTEGER(I4B), INTENT(IN) :: nRch
		REAL(DP), INTENT(IN) :: xrdn(0:)
		CHARACTER(LEN=30), INTENT(IN) :: rlab2(0:), rname(0:)
		INTEGER(I4B) i, j, nElem, status

		!gp 17-Nov-04
		CHARACTER(LEN=30), INTENT(IN) :: geoMethod		!gp 17-Nov-04 Depth or Width for col T and U of 'Reach' sheet

		!
		!CALL AllocateReachArray(Topo, nRch)
		Topo%nr = nRch
		
		!gp 17-Nov-04
		Topo%geoMethod = geoMethod		!gp 17-Nov-04 Depth or Width for col T and U of 'Reach' sheet
		
		ALLOCATE (Topo%reach(0:Topo%nr), STAT=status)
		IF (status==1) THEN 
			STOP 'Class_River:AllocateReachArray failed. Insufficient Memory!'
		END IF

		DO i=0, nRch
			Topo%reach(i)%rlab  = rlab2(i)
			Topo%reach(i)%rname = rname(i)
			Topo%reach(i)%xrdn = xrdn(i)
			IF (i==0) THEN
				Topo%reach(i)%xpm=xrdn(i)
			ELSE	
				Topo%reach(i)%xpm = (xrdn(i-1)+xrdn(i))/2
			END IF
	!		Topo%reach(i)%xup	 = xrdn (i-1)
	!		Topo%reach(i)%xdown = xrdn (i)
	!		Topo%reach(i)%elems = 1											!// hardwire to 1 for test
	!		Topo%reach(i)%elemLen = (xrdn(i)-xrdn(i-1))/Topo%reach(i)%elems
	!		nElem = nElem + Topo%reach(i)%elems					!calculate total elements
		END DO
	!	CALL AllocateElementArray(nElem)
	END FUNCTION RiverTopo_

	SUBROUTINE AllocateReachArray(Topo, nRch)

		TYPE(RiverTopo_type) Topo
		INTEGER(I4B), INTENT(IN) :: nRch
		INTEGER(I4B) status
		IF (.NOT. ASSOCIATED(Topo%reach)) THEN
			IF (nRch>0)THEN
				Topo%nr = nRch
				ALLOCATE (Topo%reach(0:Topo%nr), STAT=status)
				IF (status==1) THEN 
					STOP 'Class_River:AllocateReachArray failed. Insufficient Memory!'
				END IF

			ELSE
				PRINT *, 'ERROR:element and reach number must be great than 0'
				STOP 'Class_River:AllocateReachArray failed'
			END IF
		ELSE
			PRINT *, 'Warning: River array can only allocate once. Allocation failed!'
		END IF

	END SUBROUTINE AllocateReachArray

!	SUBROUTINE AllocateElementArray(nElem)
!		INTEGER(I4B), INTENT(IN) :: nElem
!		INTEGER(I4B) status

!		IF (.NOT. ASSOCIATED(Topo%elems)) THEN
!			IF (nElem>0)THEN
!				Topo%ne = nElem
!				ALLOCATE (Topo%elems (0:Topo%ne), STAT=status)
!				IF (status==1) THEN 
!					STOP 'Class_River:AllocateReachArray failed. Insufficient Memory!'
!				END IF
!			ELSE
!				PRINT *, 'ERROR:element and reach number must be great than 0'
!				STOP 'Class_River:AllocateElementArray failed'
!			END IF
!		ELSE
!			PRINT *, 'Warning: Topo array can only allocate once. Allocation failed!'
!		END IF

!	END SUBROUTINE AllocateElementArray

	!dynamic data types

	SUBROUTINE BuildElemTopology()
		
		
		
	END SUBROUTINE

END MODULE Class_RiverTopo