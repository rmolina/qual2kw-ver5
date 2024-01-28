module m_water_quality
    use nrtype, only: dp, nv
    implicit none
    private
    public :: water_quality_t

    type water_quality_t
        real(dp) :: te = 0.0
        real(dp) :: c(1:nv) = 0.0
        real(dp) :: ph = 7.0
    end type

end module
