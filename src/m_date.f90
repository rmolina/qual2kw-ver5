module m_date
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    private
    public :: date_t

    type date_t !11/16/04
        real (r64) year
        real (r64) month
        real (r64) day
        ! real (r64) julday !julian day
    end type date_t

!     interface date_t
!         procedure :: date_ctor
!     end interface date_t


! contains

!     !data structure constructor
!     pure function date_ctor(year, month, day)

!         real (r64), intent(in) :: year, month, day
!         type(date_t) date_ctor !11/16/04

!         date_ctor%year = year
!         date_ctor%month= month
!         date_ctor%day = day
!         date_ctor%julday = julcvt(int(month, i32), int(day, i32), int(year, i32))
!     end function date_ctor

! !calculate julian day, jan 1st =1,
!     pure function julcvt(month, day, year)

!         real (r64) julcvt
!         integer (i32), intent(in):: month, day, year
!         integer (i32) leap
!         integer (i32) modtest

! !gp edit for absoft syntax
! ! if (mod(year, 4) == 0) leap = 1
!         if (mod(year, 4) .eq. 0) leap = 1
! !! if mod(year, 4) = 0 then
! !! modtest = mod(year,4)
! ! modtest = year - int(year/4)*4
! !! if (modtest == 0) then
! ! if (modtest .eq. 0) then
! ! leap = 1
! ! end if

!         select case (month)
!           case (1)
!             julcvt = 0
!           case (2)
!             julcvt = 31
!           case (3)
!             julcvt = 59 + leap
!           case (4)
!             julcvt = 90 + leap
!           case (5)
!             julcvt = 120 + leap
!           case (6)
!             julcvt = 151 + leap
!           case (7)
!             julcvt = 181 + leap
!           case (8)
!             julcvt = 212 + leap
!           case (9)
!             julcvt = 243 + leap
!           case (10)
!             julcvt = 273 + leap
!           case (11)
!             julcvt = 304 + leap
!           case (12)
!             julcvt = 334 + leap
!         end select
!         julcvt = julcvt + day
!     end function

end module


