!sourceinbackup.f90

module class_sourcein
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use nrtype, only: nv, pii, cpw, rhow
    use m_rivertopo, only: t_rivertopo
    use class_phsolve, only: ct
    use m_water_quality, only: t_water_quality
    use class_hydraulics, only: riverhydraulics_type
    implicit none
    private
    public sourcein_, load, ndiff, npt, sourcescalc

    !derived type for diffusion load/source
    type diffusion_type !contain the raw data for diffusin load and source
        character(len=30) name
        real(r64) xdup, xddn
        real(r64) :: q =0.0 !if load, then >0; source, then <0
        real(r64) te
        real(r64) c(nv-1)
        real(r64) ph
        integer(i32) beginrch !beginning reach number
    end type diffusion_type
    !derived type for point load/source
    type point_type !
        character(len=30) name
        real(r64) x
        real(r64) :: q =0.0 !if load, then >0; source, then <0
        real(r64) :: temean=0, teamp=0, temaxtime =0
        real(r64) cmean(nv-2), camp(nv-2), cmaxtime(nv-2)
        real(r64) phmean, phamp, phmaxtime
        integer(i32) beginrch !beginning reach number
    end type point_type

    integer(i32) npt, ndiff

    !contains the original point and diffusion data
    type(point_type), pointer :: point(:)
    type(diffusion_type), pointer :: diffu(:)
    type(t_water_quality), allocatable :: load(:) !combined load from both point and diffusion
    real(r64), allocatable :: heatdiff(:), loaddiff(:,:) !diffusion load not vary by time
    real(r64), allocatable:: qpta(:), qpt (:)

contains

    ! read in point and diffusion load and abastraction data
    subroutine sourcein_(nr, nptin, ndiffin, flag, topo, hydrau, ptname, xptt, qptta, qptt, tepttmean, &
        tepttamp, tepttmaxtime, cpttmean, cpttamp, cpttmaxtime, phpttmean, &
        phpttamp, phpttmaxtime, diffname, xdup, xddn, qdifa, qdif, tedif, cdif, phind)

        integer(i32), intent(in) :: nr, nptin, ndiffin, flag
        type(t_rivertopo), intent(in) :: topo !river topology
        type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
        real(r64) xptt(:), qptta(:), qptt(:), tepttmean(:), tepttamp(:), tepttmaxtime(:)
        real(r64) cpttmean(:,:), cpttamp(:,:), cpttmaxtime(:,:)
        real(r64) phpttmean(:), phpttamp(:), phpttmaxtime(:)
        character(len=30), intent(in) :: ptname(:)
        real(r64), intent(in) :: xdup(:), xddn(:), qdifa(:), qdif(:), phind(:), tedif(:), cdif(:,:)
        character(len=30), intent(in) :: diffname(:)
        integer(i32) status(5), i

        allocate(qpta(nr), stat=status(1))
        allocate(qpt (nr), stat=status(2))
        allocate(load(nr), stat=status(3))
        allocate(heatdiff(nr), stat=status(4))
        allocate(loaddiff(nr,nv), stat=status(5))

        do i=1, 5
            if (status(i)==1) then
                stop 'ERROR:Class_SourceIn(sourceIN) Insufficient memory for dyanmic allocation'
            end if
        end do
        qpta = 0
        qpt = 0

        call pointin_(nr, nptin, topo, flag, ptname, xptt, qptta, qptt, tepttmean, tepttamp, tepttmaxtime, &
            cpttmean, cpttamp, cpttmaxtime, phpttmean, phpttamp, phpttmaxtime)
        call nonpointin_(nr, topo, flag, ndiffin, diffname, xdup, xddn, qdifa, qdif, tedif, cdif, phind)
        !generate average reach flows for hydraulics (m^3/s)
        do i = 1, nr
            if (qpt(i) < 0) qpt(i) = 0
            hydrau%reach(i)%q = hydrau%reach(i-1)%q + qpt(i) - qpta(i)
            hydrau%reach(i)%qpt=qpt(i)
            hydrau%reach(i)%qpta = qpta(i)
        end do

    end subroutine sourcein_


    subroutine pointin_(nr, nptin, topo, flag, ptname, xptt, qptta, qptt, tepttmean, tepttamp, tepttmaxtime, &
        cpttmean, cpttamp, cpttmaxtime, phpttmean, phpttamp, phpttmaxtime)

        integer(i32), intent(in) :: nr, nptin, flag
        type(t_rivertopo), intent(in) :: topo !river topology
        real(r64) xptt(:), qptta(:), qptt(:), tepttmean(:), tepttamp(:), tepttmaxtime(:)
        real(r64) cpttmean(:,:), cpttamp(:,:), cpttmaxtime(:,:)
        real(r64) phpttmean(:), phpttamp(:), phpttmaxtime(:)
        character(len=30), intent(in) :: ptname(:)

        integer(i32) i, j, status
        logical(2) cond1
        npt=nptin

        if (npt>0) then

            allocate(point(npt), stat=status)
            if (status==1) then
                stop 'Class_SourceIn:PointIn_, dynamic allocation failed!'
            end if
            do i=1, npt
                point(i)%name = ptname(i)
                point(i)%x = xptt(i)
                if (qptta(i)>0) then
                    point(i)%q = -qptta(i) !if load, then >0; source, then <0
                else
                    point(i)%q = qptt(i) !if load, then >0; source, then <0
                end if
                point(i)%temean = tepttmean(i)
                point(i)%teamp = tepttamp(i)
                point(i)%temaxtime = tepttmaxtime(i)

                do j=1, nv-2
                    point(i)%cmean(j) =cpttmean(i,j)
                    point(i)%camp(j) =cpttamp(i,j)
                    point(i)%cmaxtime(j) = cpttmaxtime(i,j)
                end do
                point(i)%phmean = phpttmean(i)
                point(i)%phamp = phpttamp(i)
                point(i)%phmaxtime = phpttmaxtime(i)
                !distribute point flows to elements for hydraulics
                do j=1, nr
                    if (flag == 1) then
                        cond1 = (xptt(i) >= topo%reach(j-1)%xrdn) .and. &
                            (xptt(i) < topo%reach(j)%xrdn)
                    else
                        cond1 = (xptt(i) <= topo%reach(j-1)%xrdn) .and. &
                            (xptt(i) > topo%reach(j)%xrdn)
                    end if
                    if (cond1) then
                        if (point(i)%q < 0) then
                            qpta(j) = qpta(j) - point(i)%q
                        else
                            qpt(j) = qpt(j) + point(i)%q
                        end if
                        point(i)%beginrch =j
                        exit
                    end if
                end do
            end do
        end if

    end subroutine pointin_


    !diffusion -- nonpointer source

    subroutine nonpointin_(nr, topo, flag, ndiffin, diffname, xdup, xddn, qdifa, qdif, tedif, cdif, phind)

        type(t_rivertopo) topo !river topology
        integer(i32), intent(in) :: ndiffin, nr, flag
        real(r64), intent(in) :: xdup(:), xddn(:), qdifa(:), qdif(:), phind(:), tedif(:), cdif(:,:)
        character(len=30), intent(in) :: diffname(:)
        integer(i32) i, j, k, status
        logical(2) cond1, cond2, cond3, cond4, cond5
        real(r64) qd, lend
        ndiff=ndiffin !number of diffusion source and abstraction

        heatdiff=0
        loaddiff=0

        if (ndiff>0) then

            allocate(diffu(ndiff), stat=status)
            if (status==1) then
                stop 'Class_SourceIn:PointIn_, dynamic allocation failed!'
            end if

            do i=1, ndiff
                diffu(i)%name = diffname(i)
                diffu(i)%xdup = xdup(i)
                diffu(i)%xddn = xddn(i)
                if (qdifa(i)>0) then
                    diffu(i)%q = -qdifa(i) !if load, then >0; source, then <0
                else
                    diffu(i)%q = qdif(i)
                end if
                diffu(i)%te = tedif(i)
                do j=1, nv-2
                    diffu(i)%c(j) =cdif(i,j)
                end do

                if (diffu(i)%c(nv - 2) == 0) diffu(i)%c(nv - 2) = 100.0
                if (phind(i)==0) then
                    diffu(i)%ph = 7.0_r64
                else
                    diffu(i)%ph = phind(i)
                end if
                !total inorganic carbon
                diffu(i)%c(nv - 1) = ct(diffu(i)%ph, diffu(i)%c(nv - 2), tedif(i), diffu(i)%c(1))
                !distribute nonpoint flows to elements for hydraulics

                qd = diffu(i)%q / (xddn(i) - xdup(i)) * flag

                do j=1, nr
                    if (flag == 1) then
                        cond1 = topo%reach(j)%xrdn < xdup(i) .or. topo%reach(j-1)%xrdn > xddn(i)
                        cond2 = topo%reach(j-1)%xrdn <= xdup(i) .and. xddn(i) <= topo%reach(j)%xrdn
                        cond3 = xdup(i) <= topo%reach(j-1)%xrdn .and. xddn(i) >= topo%reach(j)%xrdn
                        cond4 = topo%reach(j-1)%xrdn >= xdup(i) .and. topo%reach(j)%xrdn >= xddn(i)
                        cond5 = topo%reach(j)%xrdn >= xdup(i) .and. topo%reach(j)%xrdn <= xddn(i)
                    else
                        cond1 = topo%reach(j)%xrdn > xdup(i) .or. topo%reach(j-1)%xrdn < xddn(i)
                        cond2 = topo%reach(j-1)%xrdn >= xdup(i) .and. xddn(i) >= topo%reach(j)%xrdn
                        cond3 = xdup(i) >= topo%reach(j-1)%xrdn .and. xddn(i) <= topo%reach(j)%xrdn
                        cond4 = topo%reach(j-1)%xrdn <= xdup(i) .and. topo%reach(j)%xrdn <= xddn(i)
                        cond5 = topo%reach(j)%xrdn <= xdup(i) .and. topo%reach(j)%xrdn >= xddn(i)
                    end if
                    if (cond1) then
                        lend = 0
                    elseif (cond2) then
                        lend = flag * (xddn(i) - xdup(i))
                    elseif (cond3) then
                        lend = flag * (topo%reach(j)%xrdn - topo%reach(j-1)%xrdn)
                    elseif (cond4) then
                        lend = flag * (xddn(i) - topo%reach(j-1)%xrdn)
                    elseif (cond5) then
                        lend = flag * (topo%reach(j)%xrdn - xdup(i))
                    end if
                    if (qd <= 0) then
                        qpta(j) = qpta(j) - lend * qd
                    else
                        qpt(j) = qpt(j) + lend * qd
                    end if

                    if (lend>0 .and. qd>0) then
! qpt(i) = qpt(i) + lend * qd
                        heatdiff(j) = heatdiff(j) + lend * qd * diffu(i)%te
                        do k = 1, nv - 1
                            loaddiff(j, k) = loaddiff(j, k) + lend * qd * diffu(i)%c(k)
                        end do
                    end if

                end do
            end do
        end if
        !distribute point flows to elements for hydraulics
! do i = 1, nr
! qpt(i) = 0
! qpta(i) = 0
! end do

    end subroutine nonpointin_


!calculate instanteneous sources for time t
    subroutine sourcescalc(t, nr, flag)
        !gp new sub to evaluate point source sine functions and distribute loads to reaches at time t

        real(r64), intent(in) :: t
        integer(i32), intent(in) :: nr, flag
! type(rivertopo_type) topo !river topology
        integer(i32) i, j, k, kk
        real(r64) teptt(nr)
        real(r64) :: heat(nr)
        real(r64) :: loadi(nr, nv)
        logical(2) cond1, cond2, cond3, cond4, cond5
        type(t_water_quality) :: ptt
        real(r64) qd, lend

        heat = 0; loadi =0

        if (npt > 0) then
            !gp evaluate the point source diel sine functions for the current time step
            do i = 1, npt

                if (point(i)%q<=0) continue ! no need to process for abastraction

                ptt%te = sinday(t, point(i)%temean, point(i)%teamp, point(i)%temaxtime)
                if (ptt%te < 0) ptt%te = 0
                do j = 1, nv - 2
                    ptt%c(j) = sinday(t, point(i)%cmean(j), point(i)%camp(j), point(i)%cmaxtime(j))
                    if (ptt%c(j) < 0) ptt%c(j) = 0
                end do
                if (ptt%c(nv - 2) == 0) ptt%c(nv - 2) = 100
                ptt%ph = sinday(t, point(i)%phmean, point(i)%phamp, point(i)%phmaxtime)
                if (ptt%ph < 0.01) then
                    ptt%ph = 0.01_r64
                elseif (ptt%ph > 13.99) then
                    ptt%ph = 13.99_r64
                end if

                if (point(i)%phmean == 0) ptt%ph = 7.0_r64
                ptt%c(nv - 1) = ct(ptt%ph, ptt%c(nv - 2), ptt%te, ptt%c(1))

                j=point(i)%beginrch
                if (point(i)%q <= 0) then !abstraction
! qpta(j) = qpta(j) - point(i)%q
                else !load
! qpt(j) = qpt(j) + point(i)%q
                    heat(j) = heat(j) + point(i)%q * rhow * cpw * ptt%te
                    do k = 1, nv - 1
                        loadi(j, k) = loadi(j, k) + point(i)%q * ptt%c(k)
                    end do
                end if
            end do

        else !no point sources
        end if

        !generate average reach input temperatures and concentrations
        loadi=loadi+loaddiff
        heat=heat+heatdiff
        do i = 1, nr
            if (qpt(i) > 0) then
                load(i)%te = heat(i) / qpt(i)
                do j = 1, nv - 1
                    load(i)%c(j) = loadi(i, j) / qpt(i)
                end do
            else
                qpt(i) = 0
                load(i)%te = 0
                do j = 1, nv - 1
                    load(i)%c(j) = 0
                end do
            end if
! q(i) = q(i - 1) + qpt(i) - qpta(i)
        end do

    end subroutine


    pure function sinday(t, xmean, xamp, xmaxtime)
        !gp new function sinday to calculate a constituent
        ! at a particular time of day given daily mean, amplitude=(max-min)/2, and time of max (days)
        real(r64) sinday
        real(r64), intent(in) :: t, xmean, xamp, xmaxtime
        sinday = xmean + xamp * cos(2.0_r64 * pii * (t - xmaxtime))

    end function


    pure function sinday2(t, xmin, xmax, xmaxtime)
        !gp new function sinday2 to calculate a constituent
        ! at a particular time of day given daily min, max, and time of max
        real(r64) sinday2
        real(r64), intent(in) :: t, xmin, xmax, xmaxtime
        real(r64) xmean, xamp, xtheta, xomega

        xmean = (xmax + xmin) / 2.0_r64
        xamp = (xmax - xmin) / 2.0_r64
        xtheta = (xmaxtime - 0.25_r64) * 2.0_r64 * pii
        xomega = 2.0_r64 * pii
        sinday2 = xmean + xamp * sin(xomega * t - xtheta)

    end function

end module class_sourcein
