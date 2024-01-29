module m_water_quality
    use, intrinsic :: iso_fortran_env, only: r64 => real64
    use nrtype, only: nv
    implicit none
    private
    public :: water_quality_t

    type water_quality_t
        real(r64) :: te = 0.0
        real(r64) :: c(1:nv) = 0.0
        real(r64) :: ph = 7.0
    end type

end module
