module m_system_params
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use nrtype, only: lgt
    use m_date, only: date_t
    implicit none
    private
    public :: system_params_t

    type system_params_t

        !gp 23-nov-09
        !character(len=30) basinname, filename, path, title, timezone
        character(len=30) basinname, filename, path, title
        real(r64) timezone

        real(r64) dtuser
        type(date_t) today !11/16/04 current date
        real(r64) :: dt, tday =0 !timestep, time of the day
        logical(lgt) :: steadystate =.true. !identify the simulation type
        type(date_t) lastday !11/16/04 last day of simulation, for dynamic
        integer(i32) days !final time for steady state only

        !gp 29-oct-09
        !integer(i32) np, nc !
        integer(i32) np, nc !

        integer(i32) stepcount !count time steps
        integer(i32) daycount !count day simulated
        character(len=30) :: imeth = 'Euler' ! integration method
        character(len=30) :: imethph = 'Newton-Raphson' !pH method

        !gp 29-oct-04
        character(len=30) :: simhyporheicwq !simulate hyporheic pore water quality

        !03-feb-05
        character(len=30) :: showdielresults !yes or no (only used in excel vba)
        character(len=30) :: statevariables !all or temperature (used to bypass wq derivs unless 'all' is selected)

        !gp 11-jan-06
        character(len=30) :: calcsedflux !yes or no

        !gp 26-oct-07
        character(len=30) :: simalk !yes or no

        !gp 24-jun-96
        character(len=30) :: writedynamic !yes or no

    end type system_params_t

    interface system_params_t
        procedure :: system_params_ctor
    end interface system_params_t

contains
    !/* external functions */

!sytemparams data type constructor

    !gp 24-jun-09
    !gp 28-oct-04 function systemparams_(basinname, filename, path, title, year, month, day, &
    !gp timezone, dtuser, tf, imeth, imethph) result(system)
    !function systemparams_(basinname, filename, path, title, year, month, day, &
    ! timezone, dtuser, tf, imeth, imethph, simhyporheicwq, &
    ! showdielresults, statevariables, calcsedflux, simalk) result(system) !gp 26-oct-07
    function system_params_ctor(basinname, filename, path, title, year, month, day, &
        timezone, dtuser, tf, imeth, imethph, simhyporheicwq, &
        showdielresults, statevariables, calcsedflux, simalk, writedynamic) result(system)
        implicit none

        type(system_params_t) system

        !gp 23-nov-09
        !character(len=30), intent(in) :: basinname, filename, path, title, timezone
        character(len=30), intent(in) :: basinname, filename, path, title
        real(r64), intent(in) :: timezone

        real(r64), intent(in) :: year, month, day
        real(r64), intent(in) :: dtuser, tf
        character(len=30), intent(in) :: imethph !ph method
        character(len=30), intent(in) :: imeth !integration method
        character(len=30), intent(in) :: simhyporheicwq !gp 28-oct-04 yes or no to simulate hyoprheic wq
        character(len=30), intent(in) :: showdielresults !gp 03-feb-05 yes or no (only used in excel vba)
        character(len=30), intent(in) :: statevariables !gp 03-feb-05 all or temperature
        character(len=30), intent(in) :: calcsedflux !gp 11-jan-06 yes or no
        character(len=30), intent(in) :: simalk !gp 26-oct-07 yes or no

        !gp 24-jun-09
        character(len=30), intent(in) :: writedynamic !yes or no

        real(r64) dtmax

        system%basinname=basinname ; system%filename=filename; system%path=path
        system%title = title;
        system%today = date_t(year,month,day); system%timezone=timezone
        system%dtuser = dtuser ; system%days = int(tf, i32)

        !time-step control

        system%dt = 2.0 ** (int(log(system%dtuser) / log(2.0)))
        !//added by th, since int() implement differently in vb
        if (system%dt > system%dtuser) system%dt=system%dt/2.0
        !//end th
        dtmax = 2.0 ** (int(log(4.0 / 24.0) / log(2.0)))
        if (system%dt > dtmax) system%dt = dtmax

        system%np = int(tf, i32)
        system%nc = int(1.0_r64 / system%dt, i32)

        !integration method
        select case (imeth)
          case ('Adaptive step') !adaptive step method (
            system%imeth = imeth
          case ('Runge-Kutta')
            system%imeth = imeth !runge-kutta 4th order
          case default
            system%imeth = 'Euler' !default: euler method
        end select

        !ph solver method
        select case (imethph)
          case ('Bisection')
            system%imethph = imethph

            !gp 10-dec-09
            !case default
            ! system%imethph = 'newton-raphson'
          case ('Newton-Raphson')
            system%imethph = imethph
          case default
            system%imethph = 'Brent'

        end select

        !gp 29-oct-04
        system%simhyporheicwq = simhyporheicwq

        !gp 03-feb-05
        system%showdielresults = showdielresults
        system%statevariables = statevariables

        !gp 11-jan-06
        system%calcsedflux = calcsedflux

        !gp 26-oct-07
        system%simalk = simalk

        !gp 24-jun-09
        system%writedynamic = writedynamic


! system%steadystate= steadystate
! system%dt = dt

! system%tday =0 ; system%stepcount=0; system%daycount = 0

! if (steadystate .and. present(days)) then
! system%days= days
! else if (.not. steadystate .and. present(lastday)) then
! system%lastday= lastday
! else
        !print *, !warning: wrong inputs
! end if

    end function system_params_ctor

    subroutine nexttimestep(system)
        implicit none

        type(system_params_t) system
! type(systemparams), intent(inout) :: system
        system%tday= system%tday + system%dt
        system%stepcount= system%stepcount+1


        if (system%tday>=1.0) then !end of day
            system%tday= 0.0
            system%daycount=system%daycount + 1
            if (.not.system%steadystate) then !steadystate simulation
                !dynamic simulation, the next day

            end if
        end if
    end subroutine nexttimestep

    function steadystate(system)

        type(system_params_t) system
        logical(2) steadystate
        steadystate = system%steadystate
    end function

end module m_system_params
