
module m_downstream_boundary
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_constants, only: LGT, nv
    use m_water_quality, only: water_quality_t
    use m_phsolve, only: ct
    implicit none
    private
    public :: downstream_boundary_t, instanteousdownstreamboundary

    type downstream_boundary_t
        logical(lgt) :: downstreambound =.false. !specify downsteam boundary
        integer(i32) :: numdatpnt = 24 !number of data points per day
        type(water_quality_t), pointer :: dat(:) !downstream boundary time series data
    end type

    interface downstream_boundary_t
        procedure :: downstream_boundary_ctor
    end interface downstream_boundary_t

contains

    function downstream_boundary_ctor(dbfilein, downstreambound) result(db)
        type(water_quality_t), intent(in) :: dbfilein(0:)
        logical(lgt), intent(in) :: downstreambound
        type(downstream_boundary_t) db
        integer(i32) :: status

        if (downstreambound) then
            allocate (db%dat(0:db%numdatpnt-1), stat=status)
            if (status==1) then
                stop 'QUAL2K ERROR: Class_Headwater(Downstream_) Insufficient memory for dyanmic allocation'
            end if
            db%dat=dbfilein
        else
        end if
        db%downstreambound= downstreambound
    end function downstream_boundary_ctor


    ! interpolate instanteous downstream boundary condition
    subroutine instanteousdownstreamboundary(db, t, lastte, lastc, te, c)
        type(downstream_boundary_t), intent(in) :: db

! type(headwaterdownstream_type), intent(in):: hdboundary
        real(r64), intent(in) :: lastte, lastc(:), t !last reach temperature
        real(r64) te, c(:)
        real(r64) t_hr, ph
        integer(i32) t0_hr, t1_hr

        !gp 12-jan-06
        !write(10,'(4f13.4)') db%dat(t0_hr)%te + (t_hr - t0_hr) * (db%dat(t1_hr)%te - db%dat(t0_hr)%te), &
        ! db%dat(t0_hr)%c + (t_hr - t0_hr) * (db%dat(t1_hr)%c - db%dat(t0_hr)%c), &
        ! db%dat(t0_hr)%ph + (t_hr - t0_hr) * (db%dat(t1_hr)%ph - db%dat(t0_hr)%ph), &
        ! ct(ph, c(nv - 2), te, c(1))

        if (db%downstreambound) then !specified downstream boundary
            t_hr = (t - int(t)) * db%numdatpnt
            if (t_hr < db%numdatpnt -1) then
                t0_hr = int(t_hr)
                t1_hr = t0_hr + 1
            else

                !gp 17-jan-06 debug
                !t0_hr = db%numdatpnt
                t0_hr = db%numdatpnt-1

                t1_hr = 0
            end if

            !te = hw%dat(t0_hr)%te + (t_hr - t0_hr) * (hw%dat(t1_hr)%te - hw%dat(t0_hr)%te)
            te = db%dat(t0_hr)%te + (t_hr - t0_hr) * (db%dat(t1_hr)%te - db%dat(t0_hr)%te)

            !c = hw%dat(t0_hr)%c + (t_hr - t0_hr) * (hw%dat(t1_hr)%c - hw%dat(t0_hr)%c)
            c = db%dat(t0_hr)%c + (t_hr - t0_hr) * (db%dat(t1_hr)%c - db%dat(t0_hr)%c)

            !ph = hw%dat(t0_hr)%ph + (t_hr - t0_hr) * (hw%dat(t1_hr)%ph - hw%dat(t0_hr)%ph)
            ph = db%dat(t0_hr)%ph + (t_hr - t0_hr) * (db%dat(t1_hr)%ph - db%dat(t0_hr)%ph)

            !c(nv-1)=ct(ph, c(nv-2), te, c(1))
            c(nv - 1) = ct(ph, c(nv - 2), te, c(1))

        else !same as last reach's
            te=lastte
            c=lastc
        end if
    end subroutine instanteousdownstreamboundary

end module m_downstream_boundary
