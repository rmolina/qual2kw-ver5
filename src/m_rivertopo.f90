module m_rivertopo
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_reach, only: t_reach
    implicit none
    private
    public :: t_rivertopo

!	private															!unless declared public
!	public :: ne, nr, nhw, river_

    !/*module variables */
    type t_rivertopo
        integer(i32) nr									!number of reaches
!		integer(i32) ne									!number of elements
        integer(i32) :: nhw =1					!hardwired to 1 for mainstem only
!		type(element_type), pointer :: elems(:) !element array, containing topology information
        type(t_reach), dimension(:),  pointer :: reach !reaches array,

        !gp 17-nov-04
        character(len=30) geomethod		!gp 17-nov-04 depth or width for col t and u of 'reach' sheet

    end type t_rivertopo

    interface t_rivertopo
        procedure :: t_rivertopo_ctor
    end interface t_rivertopo

!	type(rivertopo_type) topo
    !steady state data types

contains
    !/* public functions */

    !pure function ne()
    !	integer(i32) ne
    !	ne = topo%ne
    !end function ne

    !pure function nr()
    !	integer(i32) nr
    !	nr = topo%nr
    !end function nr

    !pure function nhw()
    !	integer(i32) nhw
    !	nhw= topo%nhw
    !end function nhw

    !data structure constructor
    !gp 17-nov-04 function rivertopo_(nrch, rlab2, rname, xrdn) result(topo)
    function t_rivertopo_ctor(nrch, rlab2, rname, xrdn, geomethod) result(topo)		!gp 17-nov-04

        type(t_rivertopo) topo
        integer(i32), intent(in) :: nrch
        real(r64), intent(in) :: xrdn(0:)
        character(len=30), intent(in) :: rlab2(0:), rname(0:)
        integer(i32) i, j, nelem, status

        !gp 17-nov-04
        character(len=30), intent(in) :: geomethod		!gp 17-nov-04 depth or width for col t and u of 'reach' sheet

        !
        !call allocatereacharray(topo, nrch)
        topo%nr = nrch

        !gp 17-nov-04
        topo%geomethod = geomethod		!gp 17-nov-04 depth or width for col t and u of 'reach' sheet

        allocate (topo%reach(0:topo%nr), stat=status)
        if (status==1) then
            stop 'Class_River:AllocateReachArray failed. Insufficient Memory!'
        end if

        do i=0, nrch
            topo%reach(i)%rlab  = rlab2(i)
            topo%reach(i)%rname = rname(i)
            topo%reach(i)%xrdn = xrdn(i)
            if (i==0) then
                topo%reach(i)%xpm=xrdn(i)
            else
                topo%reach(i)%xpm = (xrdn(i-1)+xrdn(i))/2
            end if
            !		topo%reach(i)%xup	 = xrdn (i-1)
            !		topo%reach(i)%xdown = xrdn (i)
            !		topo%reach(i)%elems = 1											!// hardwire to 1 for test
            !		topo%reach(i)%elemlen = (xrdn(i)-xrdn(i-1))/topo%reach(i)%elems
            !		nelem = nelem + topo%reach(i)%elems					!calculate total elements
        end do
        !	call allocateelementarray(nelem)
    end function t_rivertopo_ctor

    subroutine allocatereacharray(topo, nrch)

        type(t_rivertopo) topo
        integer(i32), intent(in) :: nrch
        integer(i32) status
        if (.not. associated(topo%reach)) then
            if (nrch>0)then
                topo%nr = nrch
                allocate (topo%reach(0:topo%nr), stat=status)
                if (status==1) then
                    stop 'Class_River:AllocateReachArray failed. Insufficient Memory!'
                end if

            else
                print *, 'ERROR:element and reach number must be great than 0'
                stop 'Class_River:AllocateReachArray failed'
            end if
        else
            print *, 'Warning: River array can only allocate once. Allocation failed!'
        end if

    end subroutine allocatereacharray

!	subroutine allocateelementarray(nelem)
!		integer(i32), intent(in) :: nelem
!		integer(i32) status

!		if (.not. associated(topo%elems)) then
!			if (nelem>0)then
!				topo%ne = nelem
!				allocate (topo%elems (0:topo%ne), stat=status)
!				if (status==1) then
!					stop 'class_river:allocatereacharray failed. insufficient memory!'
!				end if
!			else
!				print *, 'error:element and reach number must be great than 0'
!				stop 'class_river:allocateelementarray failed'
!			end if
!		else
!			print *, 'warning: topo array can only allocate once. allocation failed!'
!		end if

!	end subroutine allocateelementarray

    !dynamic data types

    subroutine buildelemtopology()



    end subroutine

end module m_rivertopo