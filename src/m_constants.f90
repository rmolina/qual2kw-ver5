module m_constants
    !symbolic names for kind types of 4-, 2-, and 1-byte integers:
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: i2b = selected_int_kind(4)
    integer, parameter :: i1b = selected_int_kind(2)
    !symbolic names for kind types of single- and double-precision reals:
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
! nrtype.f90
! definition of all data types and constants

    !symbolic names for kind types of single- and double-precision complex:
    integer, parameter :: spc = kind((1.0,1.0))
    integer, parameter :: dpc = kind((1.0d0,1.0d0))
    !symbolic name for kind type of default logical:
    integer, parameter :: lgt = kind(.true.)
    !frequently used mathematical constants (with precision to spare):
    !real(sp), parameter :: pi=3.141592653589793238462643383279502884197_sp
    !real(sp), parameter :: pio2=1.57079632679489661923132169163975144209858_sp
    !real(sp), parameter :: twopi=6.283185307179586476925286766559005768394_sp
    !real(sp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_sp
    !real(sp), parameter :: euler=0.5772156649015328606065120900824024310422_sp

    ! constants used in q2k model
    !real(dp), parameter :: pii=3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: pii=3.14159265358979_dp
    real(dp), parameter :: pio2_d=1.57079632679489661923132169163975144209858_dp
    real(dp), parameter :: twopi_d=6.283185307179586476925286766559005768394_dp

    real(dp), parameter :: con = 1.74532925199433e-02_dp
    real(dp), parameter :: rhow = 1.0_dp, cpw = 1.0_dp
    real(dp), parameter :: acoeff = 0.6_dp, rl = 0.03_dp, bowen = 0.47_dp
    real(dp), parameter :: eps = 0.97_dp, sigma = 0.000000117_dp
    real(dp), parameter :: alphas = 0.0035_dp * 86400.0_dp / 10000.0_dp, hsed = 0.1_dp
    real(dp), parameter :: adam = 1.25_dp, bdam= 0.90_dp
    real(dp), parameter :: grav = 9.81_dp
    real(dp), parameter :: es = 0.001_dp, e=2.302585093_dp
    real(dp), parameter :: w0 = 1367.0_dp

    integer(i4b), parameter :: imax = 13
    integer(i4b), parameter ::hrsday = 24 !total hours per day
    !gp 30-nov-04 integer(i4b), parameter :: nv = 16, nl = 2
    integer(i4b), parameter :: nv = 17, nl = 2 !gp 30-nov-04 add generic constituent

    !gp 30-jan-06 define null value for identifying missing reach-specific rates
    real(dp), parameter :: null_val = -999_dp

    type outputsummary_type
        real(dp) min, max, avg
    end type outputsummary_type
end module m_constants
