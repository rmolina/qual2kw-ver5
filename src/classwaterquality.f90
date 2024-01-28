module class_waterquality
    use nrtype
    implicit none

    type waterquality_type
        real(dp) :: te = 0.0
        real(dp) :: c(1:nv) = 0.0
        real(dp) :: ph = 7.0
    end type

end module
