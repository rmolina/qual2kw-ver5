module m_upstream_boundary
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_constants, only: nv
    use m_water_quality, only: water_quality_t
    use m_phsolve, only: ct_function
    implicit none
    private
    public :: upstream_boundary_t, instanteousheadwater

    type upstream_boundary_t
        integer(i32) :: numdatpnt = 24 !number of data points per day
        type(water_quality_t), pointer :: dat(:) !headwater data(time:headwaterid)
    end type upstream_boundary_t

    interface upstream_boundary_t
        procedure :: upstream_boundary_ctor
    end interface upstream_boundary_t


contains

    !headwater data strucutre constructor
    function upstream_boundary_ctor(hwfilein) result(hw)
        type(upstream_boundary_t) hw
        type(water_quality_t), intent(in) :: hwfilein(:)
        integer(i32) status

! if (steadystate()) then
        allocate (hw%dat(0:hw%numdatpnt-1),stat=status)

        if (status==1) stop 'ERROR: Class_Headwater(Headwater_) Insufficient memory for dyanmic allocation'

        hw%dat = hwfilein

! else !dyanmic simulation
        !todo: implement
! stop 'dynamic simulation has not been implemented yet!'
! end if
    end function upstream_boundary_ctor


! interpolate instanteneous headwater data
    subroutine instanteousheadwater(hw, t, te, c, ph)

        type(upstream_boundary_t), intent(in):: hw
        real(r64), intent(in) :: t
        real(r64):: te, c(:), ph
        real(r64) t_hr
        integer(i32) t0_hr, t1_hr

        t_hr = (t - int(t)) * hw%numdatpnt
        if (t_hr < hw%numdatpnt-1) then
            t0_hr = int(t_hr)
            t1_hr = t0_hr + 1
        else
            t0_hr = hw%numdatpnt-1
            t1_hr = 0
        end if

        !gp 12-jan-06
        !write(10,'(4f13.4)') hw%dat(t0_hr)%te + (t_hr - t0_hr) * (hw%dat(t1_hr)%te - hw%dat(t0_hr)%te), &
        ! hw%dat(t0_hr)%c + (t_hr - t0_hr) * (hw%dat(t1_hr)%c - hw%dat(t0_hr)%c), &
        ! hw%dat(t0_hr)%ph + (t_hr - t0_hr) * (hw%dat(t1_hr)%ph - hw%dat(t0_hr)%ph), &
        ! ct(ph, c(nv-2), te, c(1))

        te = hw%dat(t0_hr)%te + (t_hr - t0_hr) * (hw%dat(t1_hr)%te - hw%dat(t0_hr)%te)
        c = hw%dat(t0_hr)%c + (t_hr - t0_hr) * (hw%dat(t1_hr)%c - hw%dat(t0_hr)%c)
        ph = hw%dat(t0_hr)%ph + (t_hr - t0_hr) * (hw%dat(t1_hr)%ph - hw%dat(t0_hr)%ph)
        !solve total carbon concentration
        c(nv-1)=ct_function(ph, c(nv-2), te, c(1))

        !gp 12-jan-06
        !write(10,'(4f13.4)') te, c, ph, c(nv-1)

    end subroutine instanteousheadwater


    function interpolate_hourly(t, c)
        !for interpolation of instantaneous values from hourly point estimates

        real(r64) interpolate_hourly
        real(r64),intent(in) :: c(:), t
        real(r64) t_hr
        integer(i32) t0_hr, t1_hr

        t_hr = (t - int(t)) * 24

        if (t_hr < 23) then
            t0_hr = int(t_hr)
            t1_hr = t0_hr + 1
        else
            t0_hr = 23
            t1_hr = 0
        end if

        interpolate_hourly = c(t0_hr) + (t_hr - t0_hr) * (c(t1_hr) - c(t0_hr))

    end function


end module m_upstream_boundary
