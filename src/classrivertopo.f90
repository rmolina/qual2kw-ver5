! classrivertopo.f90
! data structure for river topology

MODULE Class_RiverTopo
	USE nrtype
!	USE Class_Element
	USE Class_Reach
	IMPLICIT NONE

	PRIVATE :: buildTopoUpstream, buildTopoDownstream
!	PUBLIC :: ne, nr, nHw, RiverTopo_

	TYPE tributary_type
		INTEGER(I2B) begID					!headwater reach ID
		INTEGER(I2B) endID					!last reach ID in this tributary
	END TYPE
	!/*module variables */
	TYPE RiverTopo_type
		INTEGER(I2B) nr						!number of reaches
!		INTEGER(I2B) ne						!number of elements
		INTEGER(I2B) :: nHw =1					!number of headwater
!		TYPE(Element_type), POINTER :: elems(:) 		!element array, containing topology information
		TYPE(Reach_type), DIMENSION(:),  POINTER :: reach	!reaches array, 
		TYPE(tributary_type), DIMENSION(:), POINTER :: trib	!	
	END TYPE RiverTopo_type

!	TYPE(RiverTopo_type) Topo
	!steady state data types
	
	CONTAINS
	!/* public functions */
 
	!data structure constructor	
	FUNCTION RiverTopo_(nRch, nHw, rlab2, rname, xrup, xrdn) RESULT(Topo)

		TYPE(RiverTopo_type) Topo	
		INTEGER(I2B), INTENT(IN) :: nRch, nHw
		REAL(DP), INTENT(IN) :: xrup(:), xrdn(:)
		CHARACTER(LEN=30), INTENT(IN) :: rlab2(:), rname(:)
		INTEGER(I2B) i, j, nElem, status(2)

		!
		!CALL AllocateReachArray(Topo, nRch)
		Topo%nr = nRch; 		Topo%nHw=nHw
		ALLOCATE (Topo%reach(0:Topo%nr), STAT=status(1))
		ALLOCATE (Topo%trib(nHw), STAT=status(2))
		!
		DO i=1,2 
			IF (status(i)==1) THEN 
				STOP 'Class_River:AllocateReachArray failed. Insufficient Memory!'
			END IF
		END DO

		DO i=1, nRch
			Topo%reach(i)%rlab  = rlab2(i)
			Topo%reach(i)%rname = rname(i)
			Topo%reach(i)%xrup = xrup(i)
			Topo%reach(i)%xrdn = xrdn(i)
			Topo%reach(i)%xpm = (xrup(i)+ xrdn(i))/2
	!		Topo%reach(i)%xup	 = xrdn (i-1)
	!		Topo%reach(i)%xdown = xrdn (i)
	!		Topo%reach(i)%elems = 1											!// hardwire to 1 for test
	!		Topo%reach(i)%elemLen = (xrdn(i)-xrdn(i-1))/Topo%reach(i)%elems
	!		nElem = nElem + Topo%reach(i)%elems					!calculate total elements
		END DO
	!	CALL AllocateElementArray(nElem)
	END FUNCTION RiverTopo_


	SUBROUTINE BuildRiverTopology(Topo)
	TYPE(RiverTopo_type) Topo
	INTEGER(I2B) i
	INTEGER(I2B) curRchID, curTribID					                      !current reach ID
																																	!current tributary ID
	IF ((Topo%nr <=0) .OR. (Topo%nHw<=0)) THEN
		WRITE(*,*) 'QUAL2K error: Class_River(StreamTopology) nr/nHw not initialized!'
		STOP 
	END IF

	DO i = 1,Topo%nHw                                               !set headwater signs
		Topo%reach(Topo%trib(i)%begID)%IsHw = .TRUE.
	END DO

	curRchID = 1                                                    !first reach
	curTribID = 1																										!first tributary
	!build river topology
	Call buildTopoUpstream(Topo, curRchID, curTribID)
	Call buildTopoDownstream(Topo)

	DO i = 1, Topo%nr                                               !set junction reach signs
	!not a headwater and not consecutive ID with upstream
		IF ((Topo%reach(i)%upID > 0) .AND. (i - Topo%reach(i)%upID > 1)) THEN
			Topo%reach(i)%IsJct = .TRUE.
		ELSE
			Topo%reach(i)%IsJct = .FALSE.
		END IF
	END DO

	END SUBROUTINE

!----------------------------------------------------------------------------------------------
!buildTopoUpstream()
!so each reach knows who's my immediate upstream reach
!purpose: test connectivity of main stem and tributaries, beginning from the headwater reach of mainstem to
!         downstream of mainstem. IF a tributary is encountered, then go to the headwater of that tributary.
!         When finishing testing that tributary, then program resumes the test of mainstem from where it is
!         interrupted. This subroutine has the ability to recursivly test as many level of tributaries as
!         possible
!notice: headwater reach use "0" as its upstream reach ID
!----------------------------------------------------------------------------------------------
RECURSIVE SUBROUTINE buildTopoUpstream(Topo, curRchID, curTribID)
	TYPE(RiverTopo_type) Topo
	INTEGER(I2B), INTENT(INOUT) :: curRchID, curTribID
  INTEGER(I2B) :: preRchID =0																					!previous Reach ID
																																			!0: beginning of tributary
                                                                      !is always a Headwater Reach
  DO WHILE (.TRUE.)
    
    IF (Topo%reach(curRchID)%xrup <= Topo%reach(curRchID)%xrdn) THEN    !upstream > downstream ?
                                                                      !error message, exit program
      WRITE(*,*) "QUAL2K Error:Upstream location should be always greater than the downstream" &
			        ," Reach ", curRchID
			STOP
      
    END IF
    
    !check connectivity with previous reach
    IF (preRchID > 0) THEN                                             !not a headwater reach
    
      IF (Topo%reach(preRchID)%xrdn /= Topo%reach(curRchID)%xrup) THEN
				!upstream end = current start?
        WRITE(*,*) "QUAL2K Error:upstream location should equal downstream of previous reach!", &
            " Reach ", curRchID
      END IF
    END IF
    Topo%reach(curRchID)%upID = preRchID                                     !set immediate upstream reach
    
    IF ((curRchID >= Topo%nr) .OR. (Topo%reach(curRchID)%xrdn == 0)) THEN
																																		!last reach of the mainstem
                                                                     !or current tributary
      Topo%trib(curTribID)%endID = curRchID                         !assign last reach ID
      curRchID = curRchID + 1                                        !next reach ID
      curTribID = curTribID - 1                                      !parent tributary's id
      RETURN	                                                       !exit the recursion
                                                                     !and go back to resume parent
                                                                     !tributary where was interrupted
    END IF
    !else next reach in current tributary
    !housing keeping next reach
    preRchID = curRchID
    curRchID = curRchID + 1
    
    IF (Topo%reach(curRchID)%IsHw) THEN                              !headwater reach of a new tributary ?
      curTribID = curTribID + 1
      Call buildTopoUpstream(Topo, curRchID, curTribID)              !then start testing my child
                                                                     !tributary, recursion starts!!!
      !check anything wrong?
      IF (curRchID > Topo%nr) THEN                                   !last reach is not in mainstem
				WRITE(*,*) "QUAL2K Error:last reach should be a mainstem reach!", " Reach ", curRchID
        STOP
      ELSEIF (Topo%reach(curRchID)%IsHw) THEN                        !next is immediately a new tributary ?
        !two tributary flow into the same location of their main stream
				WRITE(*,*) "QUAL2K Error:two tributaries are not allowed to enter the same location of their mainstem", &
					" Reach ", curRchID
        STOP 
      END IF
    END IF
  END DO
END SUBROUTINE

!build river topology for downstream
!so each reach knows who's my immediate downstream reach
SUBROUTINE buildTopoDownstream(Topo)
	TYPE(RiverTopo_type) Topo
	INTEGER(I2B) i

	!reverse array index and upstream reach ID array to get downstream reach ID array
	DO i = 1, Topo%nr
		IF (Topo%reach(i)%upID > 0) THEN
			Topo%reach(Topo%reach(i)%upID)%dwID = i
		END IF
	END DO

	!if no reach ID set, then it must be last reach in its tributary
	DO i = 1,	Topo%nr
		IF (Topo%reach(i)%dwID <= 0) THEN
			Topo%reach(i)%dwID = i + 1                            !downstream ID is alway 1 greater
		END IF
	END DO
END SUBROUTINE


END MODULE Class_RiverTopo
