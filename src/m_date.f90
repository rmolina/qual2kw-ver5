module m_date
    use, intrinsic :: iso_fortran_env, only: r64 => real64
    implicit none
    private
    public :: date_t

    type date_t
        real (r64) year
        real (r64) month
        real (r64) day
    end type date_t

end module
