module m_reach
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    implicit none
    private
    public :: t_reach

    type t_reach
        character(len=30) :: rlab='', rname =''
        real(r64) :: xrdn = 0
        real(r64) :: xpm = 0
        integer(i32) :: upid = 0, dwnid = 0 !upstream and downstream element id
    end type t_reach

end module m_reach
