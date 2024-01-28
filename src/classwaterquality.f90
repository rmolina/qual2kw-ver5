module m_water_quality
    use nrtype, only: dp, nv
    implicit none
    private
    public :: t_water_quality

    type t_water_quality
        real(dp) :: te = 0.0
        real(dp) :: c(1:nv) = 0.0
        real(dp) :: ph = 7.0
    end type

end module
