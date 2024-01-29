MODULE Class_SystemParams
    USE nrtype
    USE m_date, only: date_t
    IMPLICIT NONE

! PRIVATE !/* All data and subroutines are private unless declare as public */
! PUBLIC :: SystemParams_, NextTimeStep, steadystate

    TYPE SystemParams

        !GP 23-Nov-09
        !CHARACTER(LEN=30) BASINNAME, FILENAME, PATH, TITLE, timezone
        CHARACTER(LEN=30) BASINNAME, FILENAME, PATH, TITLE
        REAL(DP) timezone

        REAL(DP) dtuser
        TYPE(date_t) today !11/16/04 current date
        REAL(DP) :: dt, tday =0 !timestep, time of the day
        LOGICAL(LGT) :: steadystate =.TRUE. !identify the simulation type
        TYPE(date_t) LastDay !11/16/04 last day of simulation, for dynamic
        INTEGER(I4B) days !final time for steady state only

        !gp 29-Oct-09
        !INTEGER(I4B) np, nc !
        INTEGER(I4B) np, nc !

        INTEGER(I4B) stepCount !count time steps
        INTEGER(I4B) dayCount !count day simulated
        CHARACTER(LEN=30) :: IMeth = 'Euler' ! integration method
        CHARACTER(LEN=30) :: IMethpH = 'Newton-Raphson' !pH method

        !gp 29-Oct-04
        CHARACTER(LEN=30) :: simHyporheicWQ !simulate hyporheic pore water quality

        !03-Feb-05
        CHARACTER(LEN=30) :: showDielResults !Yes or No (only used in Excel VBA)
        CHARACTER(LEN=30) :: stateVariables !All or Temperature (used to bypass WQ derivs unless 'All' is selected)

        !gp 11-Jan-06
        CHARACTER(LEN=30) :: calcSedFlux !Yes or No

        !gp 26-Oct-07
        CHARACTER(LEN=30) :: simAlk !Yes or No

        !gp 24-Jun-96
        CHARACTER(LEN=30) :: writeDynamic !Yes or No

    END TYPE SystemParams

    !declare the system parameter variables
! TYPE(SystemParams) system

CONTAINS
    !/* external functions */

!sytemParams data type constructor

    !gp 24-Jun-09
    !gp 28-Oct-04 FUNCTION SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
    !gp timezone, dtuser, tf, IMeth, IMethpH) RESULT(system)
    !FUNCTION SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
    ! timezone, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, &
    ! showDielResults, stateVariables, calcSedFlux, simAlk) RESULT(system) !gp 26-Oct-07
    FUNCTION SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
        timezone, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, &
        showDielResults, stateVariables, calcSedFlux, simAlk, writeDynamic) RESULT(system)
        IMPLICIT NONE

        TYPE(SystemParams) system

        !GP 23-Nov-09
        !CHARACTER(LEN=30), INTENT(IN) :: BASINNAME, FILENAME, PATH, TITLE, timezone
        CHARACTER(LEN=30), INTENT(IN) :: BASINNAME, FILENAME, PATH, TITLE
        REAL(DP), INTENT(IN) :: timezone

        REAL(DP), INTENT(IN) :: year, month, day
        REAL(DP), INTENT(IN) :: dtuser, tf
        CHARACTER(LEN=30), INTENT(IN) :: IMethpH !pH method
        CHARACTER(LEN=30), INTENT(IN) :: IMeth !integration method
        CHARACTER(LEN=30), INTENT(IN) :: simHyporheicWQ !gp 28-Oct-04 Yes or No to simulate hyoprheic WQ
        CHARACTER(LEN=30), INTENT(IN) :: showDielResults !gp 03-Feb-05 Yes or No (only used in Excel VBA)
        CHARACTER(LEN=30), INTENT(IN) :: stateVariables !gp 03-Feb-05 All or Temperature
        CHARACTER(LEN=30), INTENT(IN) :: calcSedFlux !gp 11-Jan-06 Yes or No
        CHARACTER(LEN=30), INTENT(IN) :: simAlk !gp 26-Oct-07 Yes or No

        !gp 24-Jun-09
        CHARACTER(LEN=30), INTENT(IN) :: writeDynamic !Yes or No

        REAL(DP) dtmax

        system%BASINNAME=BASINNAME ; system%FILENAME=FILENAME; system%PATH=PATH
        system%TITLE = TITLE;
        system%today = date_t(year,month,day); system%timezone=timezone
        system%dtuser = dtuser ; system%days = tf

        !time-step control

        system%dt = 2.0 ** (Int(Log(system%dtuser) / Log(2.0)))
        !//added by th, since int() implement differently in VB
        IF (system%dt > system%dtuser) system%dt=system%dt/2.0
        !//end th
        dtmax = 2.0 ** (Int(Log(4.0 / 24.0) / Log(2.0)))
        If (system%dt > dtmax) system%dt = dtmax

        system%np = tf
        system%nc = 1.0_dp /system%dt

        !integration method
        SELECT CASE (IMeth)
          CASE ('Adaptive step') !Adaptive step method (
            system%IMeth = IMeth
          CASE ('Runge-Kutta')
            system%IMeth = IMeth !Runge-Kutta 4th order
          CASE DEFAULT
            system%IMeth = 'Euler' !default: Euler method
        END SELECT

        !pH solver method
        SELECT CASE (IMethpH)
          CASE ('Bisection')
            system%IMethpH = IMethpH

            !gp 10-Dec-09
            !CASE DEFAULT
            ! system%IMethpH = 'Newton-Raphson'
          CASE ('Newton-Raphson')
            system%IMethpH = IMethpH
          CASE DEFAULT
            system%IMethpH = 'Brent'

        END SELECT

        !gp 29-Oct-04
        system%simHyporheicWQ = simHyporheicWQ

        !gp 03-Feb-05
        system%showDielResults = showDielResults
        system%stateVariables = stateVariables

        !gp 11-Jan-06
        system%calcSedFlux = calcSedFlux

        !gp 26-Oct-07
        system%simAlk = simAlk

        !gp 24-Jun-09
        system%writeDynamic = writeDynamic


! system%steadystate= steadystate
! system%dt = dt

! system%tday =0 ; system%stepCount=0; system%dayCount = 0

! IF (steadystate .AND. PRESENT(days)) THEN
! system%days= days
! ELSE IF (.NOT. steadystate .AND. PRESENT(LastDay)) THEN
! system%LastDay= LastDay
! ELSE
        !PRINT *, !Warning: Wrong inputs
! END IF

    END FUNCTION SystemParams_

    SUBROUTINE NextTimeStep(system)
        IMPLICIT NONE

        TYPE(SystemParams) system
! TYPE(SystemParams), INTENT(INOUT) :: system
        system%tday= system%tday + system%dt
        system%stepCount= system%stepCount+1


        IF (system%tday>=1.0) THEN !END of day
            system%tday= 0.0
            system%dayCount=system%dayCount + 1
            IF (.NOT.system%steadystate) THEN !steadystate simulation
                !dynamic simulation, the next day

            END IF
        END IF
    END SUBROUTINE NextTimeStep

    FUNCTION steadystate(system)

        TYPE(SystemParams) system
        LOGICAL(2) steadystate
        steadystate = system%steadystate
    END FUNCTION

END MODULE Class_SystemParams
