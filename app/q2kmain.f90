! Q2KMAIN.f90
!
!****************************************************************************
!
! PROGRAM: Q2KFORTRANTEST
!
! PURPOSE: Entry point for the console application.
!
! The following water quality state variables are simulated in QUAL2Kw ver5.x
! (i=reach, j=layer [1=water column, 2=sediment temperature/hyporheic pore water])
!
! Variable Variable name Units
! ---------- ------------- -----
! c(i, 1, j) Conductivity umho/cm
! c(i, 2, j) Inorganic Suspended Solids mgD/L (D=dry weight)
! c(i, 3, j) Dissolved Oxygen mgO/L (O=oxygen)
! c(i, 4, j) Slow CBOD mgO/L
! c(i, 5, j) Fast CBOD mgO/L
! c(i, 6, j) Organic N (dissolved and particulate) ugN/L
! c(i, 7, j) Ammonia N ugN/L
! c(i, 8, j) Nitrate + Nitrite N ugN/L
! c(i, 9, j) Organic P (dissolved and particulate) ugP/L
! c(i, 10, j) Soluble Reactive P ugP/L
! c(i, 11, j) Phytoplankton ugA/L (A=chlorophyll a)
! c(i, 12, j) Particulate Organic Mater (POM) mgD/L
! c(i, 13, j) Pathogen indicator organism number/volume
! c(i, 14, j) Generic constituent user defined
! c(i, nv-2, j) Alkalinity mgCaCO3/L
! c(i, nv-1, j) Total Inorganic C moles/L
! c(i, nv, 1) bottom algae gD/m^2
! c(i, nv, 2) hyporheic biofilm gD/m^2
! Te(i, j) temperature deg C
! INb(i) Internal N in bottom algae cells mgN/m^2
! IPb(i) Internal P in bottom algae cells mgP/m^2
!
!****************************************************************************

program q2kmain
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_hydraulics, only: riverhydraulics_type
    use m_meteorology, only: meteorology_t
    use m_upstream_boundary, only: upstream_boundary_t
    use m_downstream_boundary, only: downstream_boundary_t
    use m_rates, only: rates_t
    use m_solar_calc, only: solar_type
    use m_system_params, only: system_params_t
    use m_rivertopo, only: t_rivertopo
    use m_output, only: outdata_t, output
    use m_readfile, only: readinputfile
    use m_integration, only: integration


!gp 04-feb-05
! use dfport !11/16/04


    implicit none




    ! variables
    character(len=260) :: infile, outfile !input & output file name
    !gp long file names are limited to 255 characters (260 for full paths)
    type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
    type(meteorology_t) meteo !meteology information
    type(upstream_boundary_t) hw !headwater
    type(downstream_boundary_t) db !downstream boundary
    type(rates_t) rates !stoch, reaction, temperature and all other rate
    type(solar_type) :: solar !solar radiation
    type(system_params_t) sys !declare the system parameter variables
    type(t_rivertopo) topo !river topology
    type(outdata_t) prout
    integer(i32) begintime, endtime


!gp 04-feb-05
! character(len=260) ::msgfile, dirname !11/16/04
! integer(i32) istat !11/16/04
    character(len=260) ::msgfile !gp long file names are limited to 255 characters (260 for full paths)


    call system_clock(begintime)

    write(*,*)
    WRITE(*,'(50A)') ' QUAL2Kw version 5.1'
    WRITE(*,'(50A)') ' Department of Ecology and Tufts University'
    WRITE(*,*)
    WRITE(*,'(50A)') ' G.J. Pelletier, S.C. Chapra, and Hua Tao'
    !WRITE(*,'(50A)') ' modified from QUAL2K ver 1.4 by Hua Tao and S.C. Chapra'
    WRITE(*,*)
    WRITE(*,'(50A)') ' Program is running, please wait...'
    WRITE(*,*)


    msgFile = 'message.dat'


    open (unit=8, file=msgfile, status='old', action='read')
    read(8,*) infile, outfile

    close (8)


!gp debug 04-feb-05
!open (unit=11, file='debug.out', status='replace', action='readwrite')
!write(11,*) sdir
!write(11,*) msgfile
!write(11,*) infile, outfile
!close (11)


    !infile='c:\research\qual2k\q2k\input\bc092187.q2k'
    !infile='c:\research\qual2k\input\bc092187thackp1.q2k'
! infile='c:\research\qual2k\input\bc092187v1_2fortran.q2k'
! infile='c:\research\qual2k\input\srdummyfortran.q2k'
    open (unit=8, file=infile, status='old', action='read')
    !read in file

!gp 12-jan-06, 21-nov-06
! open (unit=10, file='debug.txt', status='replace', action='write')
! open (unit=11, file='debug2.txt', status='replace', action='write')

    call readinputfile(sys, hydrau, meteo, hw, db, rates, topo, solar)


    !finishing reading
    close (8)

    !do simulation
    !gp 17-nov-04 prout= outdata_(topo%nr)
    prout= outdata_t(topo%nr, sys) !gp 17-nov-04 pass sys to allocate dynamic diel arrays

! outfile= 'c:\research\qual2k\input\bc092187v1_2fortran.out'
! outfile= 'c:\research\qual2k\input\srdummyfortran.out'
    open (unit=8, file=outfile, status='replace', action='write')

    !gp 23-jun-09
    if (sys%writedynamic == "Yes") then
        open (unit=12, file='dynamic.txt', status='replace', action='write')
        write(12,'(44a13)') 'reach', 't', 'tempw', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'norg', &
            'nh4n', 'no23n', 'porg', 'srp', 'phyto', 'detritus', 'pathogen', 'generic', 'alk', 'ph', &
            'bot algae', 'bot algae', 'tss', 'tp', 'tn', 'dosat', 'nh3ui', 'int n', 'int p', &
            'phitemp', 'philight', 'phinitr', 'phiphos', 'phicarb', 'phitotal', &
            'botalgphoto', 'botalgresp', 'botalgdeath', 'botalggrow', &
            'botalgphoto', 'botalgresp', 'botalgdeath', 'botalggrow'
    end if


    call integration(sys, rates, meteo, solar, hw, db, &
        hydrau, prout, topo%nr)

    !gp 23-jun-09
    if (sys%writedynamic == "Yes") then
        close (12)
    end if

    call output(prout, topo%nr, topo, hydrau, rates, sys)

    call system_clock(endtime)

    write(*,*) 'elapsed time: ', (endtime-begintime)/1000.0, ' seconds'
    close (8)

! call output(temn, temx, teav, osav, phsav, cmn, cmx, cav, &
! phmn, phmx, phav, tnmn, tnmx, tnav, &
! tpmn, tpmx, tpav, nh3mn, nh3mx, nh3av, time1, nj)

    ! body of q2kfortrantest

!gp 12-jan-06, 21-nov-06
! close (10)
! close (11)

end program q2kmain




