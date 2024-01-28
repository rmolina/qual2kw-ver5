! classrivertopo.f90
! data structure for reach and elements
module m_rivertopo
    use nrtype
!	use class_element
    use m_reach
    implicit none
    private
    public :: rivertopo_type

!	private															!unless declared public
!	public :: ne, nr, nhw, river_

    !/*module variables */
    type rivertopo_type
        integer(i4b) nr									!number of reaches
!		integer(i4b) ne									!number of elements
        integer(i4b) :: nhw =1					!hardwired to 1 for mainstem only
!		type(element_type), pointer :: elems(:) !element array, containing topology information
        type(t_reach), dimension(:),  pointer :: reach !reaches array,

        !gp 17-nov-04
        character(len=30) geomethod		!gp 17-nov-04 depth or width for col t and u of 'reach' sheet

    end type rivertopo_type

    interface rivertopo_type
        procedure :: rivertopo_
    end interface rivertopo_type

!	type(rivertopo_type) topo
    !steady state data types

contains
    !/* public functions */

    !pure function ne()
    !	integer(i4b) ne
    !	ne = topo%ne
    !end function ne

    !pure function nr()
    !	integer(i4b) nr
    !	nr = topo%nr
    !end function nr

    !pure function nhw()
    !	integer(i4b) nhw
    !	nhw= topo%nhw
    !end function nhw

    !data structure constructor
    !gp 17-nov-04 function rivertopo_(nrch, rlab2, rname, xrdn) result(topo)
    function rivertopo_(nrch, rlab2, rname, xrdn, geomethod) result(topo)		!gp 17-nov-04

        type(rivertopo_type) topo
        integer(i4b), intent(in) :: nrch
        real(dp), intent(in) :: xrdn(0:)
        character(len=30), intent(in) :: rlab2(0:), rname(0:)
        integer(i4b) i, j, nelem, status

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
    end function rivertopo_

    subroutine allocatereacharray(topo, nrch)

        type(rivertopo_type) topo
        integer(i4b), intent(in) :: nrch
        integer(i4b) status
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
!		integer(i4b), intent(in) :: nelem
!		integer(i4b) status

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
