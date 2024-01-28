!handle meteology data interpolation
module m_meteorology
    use, intrinsic :: iso_fortran_env, only: dp=>real64, i4b=>int32
    !use class_rivertopo
    !use class_systemparams, only: steadystate
    implicit none

    private !all variables and functions are private, unless declare public
    !only the following subroutines are public
    public meteodata_, instanteousmeteo, t_meteorolody

    type t_meteorolody

        !gp 16-jul-08
        !real(dp), dimension(:,:), pointer :: shadehh, tahh, tdhh, uwhh, cchh
        real(dp), dimension(:,:), pointer :: shadehh, tahh, tdhh, uwhh, cchh, solarhh

        !instantenous meteorology data

        !gp 16-jul-08
        !real(dp), dimension(:),   pointer :: shadet, ta, td, uw, cc
        real(dp), dimension(:),   pointer :: shadet, ta, td, uw, cc, solart

        integer(i4b) :: numdatpnt= 24 !use to identify the size of time series
        !hardwire to (0-23 hour) for now
    end type t_meteorolody

contains

    !gp 16-jul-08
    !function meteodata_(nr, shadehhin, tahhin, tdhhin, uwhhin, cchhin) result(met)
    function meteodata_(nr, shadehhin, tahhin, tdhhin, uwhhin, cchhin, solarhhin) result(met)

        type(t_meteorolody) met
        integer(i4b), intent(in) :: nr

        !gp 16-jul-08
        !real(dp), dimension(0:,:), intent(in) :: shadehhin, tahhin, tdhhin, uwhhin, cchhin
        real(dp), dimension(0:,:), intent(in) :: shadehhin, tahhin, tdhhin, uwhhin, cchhin, solarhhin

        !gp 16-jul-08
        !integer(i4b) status(10), i
        integer(i4b) status(12), i

        !if (steadystate()) then

        !end if

        if (nr > 0) then
            allocate (met%shadehh(0:met%numdatpnt-1 ,nr), stat=status(1))
            allocate (met%tahh(0:met%numdatpnt-1,nr), stat=status(2))
            allocate (met%tdhh(0:met%numdatpnt-1,nr), stat=status(3))
            allocate (met%uwhh(0:met%numdatpnt-1,nr), stat=status(4))
            allocate (met%cchh(0:met%numdatpnt-1,nr), stat=status(5))

            allocate (met%shadet(nr), stat=status(6))
            allocate (met%ta(nr), stat=status(7))
            allocate (met%td(nr), stat=status(8))
            allocate (met%uw(nr), stat=status(9))
            allocate (met%cc(nr), stat=status(10))

            !gp 16-jul-08
            allocate (met%solarhh(0:met%numdatpnt-1,nr), stat=status(11))
            allocate (met%solart(nr), stat=status(12))

            do i=1, 10
                if (status(i)==1) stop 'ERROR: Class_Meteo:AllocateMeteoDataArray.Allocation Failed'
            end do

            met%shadehh = shadehhin
            met%tahh = tahhin
            met%tdhh = tdhhin
            met%uwhh = uwhhin
            met%cchh = cchhin

            !gp 16-jul-08
            met%solarhh = solarhhin

        else
            print *, 'ERROR:element number must be great than 0'
            stop 'Class_Meteo:AllocateMeteoDataArray failed'
        end if

    end function meteodata_


!interpolate hourly meteology data
    subroutine instanteousmeteo(nr, t, met)

        integer(i4b), intent(in) :: nr
        type(t_meteorolody), intent(inout) :: met
        real(dp), intent(in) :: t

        call interpolatehelper(nr, t, met%ta, met%tahh, met%numdatpnt)
        call interpolatehelper(nr, t, met%td, met%tdhh, met%numdatpnt)
        call interpolatehelper(nr, t, met%uw, met%uwhh, met%numdatpnt)
        call interpolatehelper(nr, t, met%cc, met%cchh, met%numdatpnt)

        !gp 16-jul-08
        call interpolatehelper(nr, t, met%solart, met%solarhh, met%numdatpnt)

        !shade use different method
        call interpolateshade(nr, t, met%shadet, met%shadehh, met%numdatpnt)
    end subroutine instanteousmeteo


    !private subroutine
    subroutine interpolatehelper(nr, t, interp, c, upbound)
        !for interpolation of instantaneous values from hourly point estimates of 2d array
        !input arrays have dimensions of (ihour, k) where ihour is hour (0 to 23)
        !and k is variable number (1 to 15) or reach number (1-1000)
        integer(i4b), intent(in) :: nr
        real(dp), intent(out) :: interp(:)
        real(dp), intent(in) :: t, c (0:,:)
        integer(i4b), intent(in):: upbound
        real(dp) t_hr
        integer(i4b) t0_hr, t1_hr, k

        t_hr = (t - int(t)) * (upbound)

        if (t_hr < upbound-1) then
            t0_hr = int(t_hr)
            t1_hr = t0_hr + 1
        else
            t0_hr = upbound-1
            t1_hr = 0
        end if

        do k=1, nr
            interp(k) = c(t0_hr, k) + (t_hr - t0_hr) * (c(t1_hr, k) - c(t0_hr, k))
        end do

    end subroutine


    subroutine interpolateshade(nr, t, hourlyshade, shadehh, upbound)
        !input values of hourly shade are integrated values that apply to the entire hour
        !therefore this function assigns the constant hourly values for each hourly period of the simulation

        !input arrays have dimensions of (ihour, k) where ihour is hour (0 to 23)
        !and k is variable number (1 to 15) or reach number (1-1000)
        integer(i4b), intent(in) :: nr
        real(dp), intent(out) :: hourlyshade(:)
        real(dp), intent(in):: t, shadehh(0:,:)
        integer(i4b), intent(in):: upbound
        real(dp) t_hr
        integer(i4b) t0_hr, t1_hr, k

        t_hr = (t - int(t)) * (upbound)

        if (t_hr < upbound-1) then
            t0_hr = int(t_hr)
            t1_hr = t0_hr + 1
        else
            t0_hr = upbound-1
            t1_hr = 0
        end if

        do k=1, nr
            hourlyshade(k) = shadehh(t0_hr, k)
        end do
    end subroutine interpolateshade

end module m_meteorology
