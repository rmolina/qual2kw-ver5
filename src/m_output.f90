module m_output
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use nrtype, only: nv, nl
    use class_hydraulics, only: riverhydraulics_type
    use class_phsolve, only: ph_solver
    use class_sourcein, only: load, npt, ndiff, sourcescalc
    use class_systemparams, only: SystemParams
    use m_rates, only: rates_t
    use m_rivertopo, only: t_rivertopo
    use m_oxygen, only: oxygen_saturation
    implicit none
    private
    public :: outdata_t, output

    type outdata_t
        !gp 24-oct-04 integer(i32) nj
        integer(i32) nj !gp need long integer for nj
        real(r64), dimension(:), pointer :: tdy
        !gp 27-oct-04 add dimension for nl
        !gp real(r64), dimension(:,:), pointer :: cmn, cmx, cav !constituents concentrations
        !gp real(r64), dimension(:), pointer :: temn, temx, teav, osav !temperature, saturated do
        !gp real(r64), dimension(:), pointer :: phmn, phmx, phav, phsav !ph
        !gp real(r64), dimension(:), pointer :: tnmn, tnmx, tnav !total nitrogen
        !gp real(r64), dimension(:), pointer :: tpmn, tpmx, tpav !total phosphorus
        !gp real(r64), dimension(:), pointer :: nh3mn, nh3mx, nh3av !nh3
        !gp real(r64), dimension(:,:), pointer :: phpr, ninbpr, nipbpr !bottom algae luxury update
        !gp real(r64), dimension(:,:), pointer :: phitotalsavepr, phitsavepr, philsavepr, &
        !gp phinsavepr, phipsavepr, phicsavepr !gp 20-oct-04 bottom algae growth limitation factors
        !gp real(r64), dimension(:,:,:), pointer :: cpr, tepr !print out
        real(r64), dimension(:,:,:), pointer :: cmn, cmx, cav !constituents concentrations (nr, nv, nl)
        real(r64), dimension(:,:), pointer :: temn, temx, teav, osav !temperature, saturated do (nr, nl)
        real(r64), dimension(:,:), pointer :: phmn, phmx, phav, phsav !ph (nr, nl)
        real(r64), dimension(:,:), pointer :: tnmn, tnmx, tnav !total nitrogen (nr, nl)
        real(r64), dimension(:,:), pointer :: tpmn, tpmx, tpav !total phosphorus (nr, nl)
        real(r64), dimension(:,:), pointer :: nh3mn, nh3mx, nh3av !nh3 (nr, nl)
        real(r64), dimension(:,:,:), pointer :: phpr !bottom algae luxury update (nr, nj, nl)
        real(r64), dimension(:,:), pointer :: ninbpr, nipbpr !bottom algae luxury update (nr, nj)

        !gp 20-oct-04 growth limitation factors for bottom algae (nr, nj)
        real(r64), dimension(:,:), pointer :: phitotalsavepr, phitsavepr, philsavepr, phinsavepr, phipsavepr, phicsavepr

        !gp 28-oct-04 diagenesis flux between sediment/water (nr, nj)
        real(r64), dimension(:,:), pointer :: diagfluxdopr, diagfluxcbodpr, diagfluxnh4pr, diagfluxno3pr, &
            diagfluxsrppr, diagfluxicpr

        !gp 28-oct-04 hyporheic exchange flux between sediment/water (nr, nj)
        real(r64), dimension(:,:), pointer :: hypofluxdopr, hypofluxcbodpr, hypofluxnh4pr, hypofluxno3pr, hypofluxsrppr, &
            hypofluxicpr

        !gp 15-nov-04 reach-average daily-average flux between sediment/water (nr)
        real(r64), dimension(:), pointer :: diagfluxdoav, diagfluxcbodav, diagfluxnh4av, diagfluxno3av, diagfluxsrpav
        real(r64), dimension(:), pointer :: hypofluxdoav, hypofluxcbodav, hypofluxnh4av, hypofluxno3av, hypofluxsrpav

        !gp 11-jan-05 reach-min/max/mean cell quota mgn/gd and mgp/gd (nr)
        real(r64), dimension(:), pointer :: ninbmn, ninbmx, ninbav, nipbmn, nipbmx, nipbav

        real(r64), dimension(:,:,:), pointer :: tepr !print out (nr, nj, nl)
        real(r64), dimension(:,:,:,:), pointer :: cpr !print out (nr, nv, nj, nl)

        !gp 05-jul-05 heat/do/co2 fluxes
        real(r64), dimension(:,:), pointer :: pr_saveheatfluxjsnt, pr_saveheatfluxlongat, pr_saveheatfluxback, pr_saveheatfluxconv
        real(r64), dimension(:,:), pointer :: pr_saveheatfluxevap, pr_saveheatfluxjsed, pr_saveheatfluxjhyporheic, &
            pr_saveheatfluxtribs
        real(r64), dimension(:,:), pointer :: pr_saveheatfluxadvecdisp
        real(r64), dimension(:,:), pointer :: pr_savedofluxreaer, pr_savedofluxcbodfast, pr_savedofluxcbodslow, pr_savedofluxnitrif
        real(r64), dimension(:,:), pointer :: pr_savedofluxphytoresp, pr_savedofluxphytophoto, pr_savedofluxbotalgresp, &
            pr_savedofluxbotalgphoto
        real(r64), dimension(:,:), pointer :: pr_savedofluxsod, pr_savedofluxcod, pr_savedofluxhyporheic, pr_savedofluxtribs
        real(r64), dimension(:,:), pointer :: pr_savedofluxadvecdisp
        real(r64), dimension(:,:), pointer :: pr_saveco2fluxreaer, pr_saveco2fluxcbodfast, pr_saveco2fluxcbodslow
        real(r64), dimension(:,:), pointer :: pr_saveco2fluxphytoresp, pr_saveco2fluxphytophoto, pr_saveco2fluxbotalgresp, &
            pr_saveco2fluxbotalgphoto
        real(r64), dimension(:,:), pointer :: pr_saveco2fluxsod, pr_saveco2fluxhyporheic, pr_saveco2fluxtribs
        real(r64), dimension(:,:), pointer :: pr_saveco2fluxadvecdisp

        !gp 25-jun-09
        real(r64), dimension(:), pointer :: av_botalgphoto, av_botalgresp, av_botalgdeath, av_botalgnetgrowth

    end type outdata_t

    interface outdata_t
        procedure :: outdata_ctor
    end interface outdata_t

contains

!gp 17-nov-04 function outdata_(nr) result(pr)
    function outdata_ctor(nr, sys) result(pr) !gp 17-nov-04 pass sys to minimize size of dynamic diel arrays

        !gp 17-nov-04

        integer(i32), intent(in) :: nr
        type(outdata_t) pr

        !gp 05-jul-05 integer(i32) i,status(0:60) !gp 11-jan-05
        !gp 25-jun-09 integer(i32) i,status(0:93) !gp 05-jul-05
        integer(i32) i,status(0:97) !gp 25-jun-09

        !gp 17-nov-04
        type(systemparams) sys
        integer(i32) nsteps
        if (sys%imeth == "Adaptive step") then
            nsteps = 2400
        else
            nsteps = sys%nc !minimizes array sizes for euler and rk4 integration
        end if

        status=0

        !gp 17-nov-04 allocate(pr%tdy(0:2400), stat=status(0))
        allocate(pr%tdy(0:nsteps), stat=status(0)) !gp 17-nov-04 replace 2400 with nsteps for all dynamic diel output arrays below

        !gp allocate(pr%cpr(0:nr, nv, 0:2400), stat=status(1))
        !gp allocate(pr%tepr(0:nr, 0:2400, 2), stat=status(2))
        !gp allocate(pr%phpr(0:nr, 0:2400), stat=status(3))
        allocate(pr%cpr(0:nr, nv, 0:nsteps, nl), stat=status(1))
        allocate(pr%tepr(0:nr, 0:nsteps, nl), stat=status(2))
        allocate(pr%phpr(0:nr, 0:nsteps, nl), stat=status(3)) !gp end new block

        allocate(pr%ninbpr(0:nr, 0:nsteps), stat=status(4))
        allocate(pr%nipbpr(0:nr, 0:nsteps), stat=status(5))

        !gp allocate(pr%temn(0:nr), stat=status(7))
        !gp allocate(pr%temx(0:nr), stat=status(8))
        !gp allocate(pr%teav(0:nr), stat=status(9))
        !gp allocate(pr%osav(0:nr), stat=status(10))
        !gp allocate(pr%phsav(0:nr), stat=status(11))
        !gp allocate(pr%cmn(0:nr,nv), stat=status(12))
        !gp allocate(pr%cmx(0:nr,nv), stat=status(13))
        !gp allocate(pr%cav(0:nr,nv), stat=status(14))
        !gp allocate(pr%phmn(0:nr), stat=status(15))
        !gp allocate(pr%phmx(0:nr), stat=status(16))
        !gp allocate(pr%phav(0:nr), stat=status(17))
        !gp allocate(pr%tnmn(0:nr), stat=status(18))
        !gp allocate(pr%tnmx(0:nr), stat=status(19))
        !gp allocate(pr%tnav(0:nr), stat=status(20))
        !gp allocate(pr%tpmn(0:nr), stat=status(21))
        !gp allocate(pr%tpmx(0:nr), stat=status(22))
        !gp allocate(pr%tpav(0:nr), stat=status(23))
        !gp allocate(pr%nh3mn(0:nr), stat=status(24))
        !gp allocate(pr%nh3mx(0:nr), stat=status(25))
        !gp allocate(pr%nh3av(0:nr), stat=status(26))
        allocate(pr%temn(0:nr, nl), stat=status(7))
        allocate(pr%temx(0:nr, nl), stat=status(8))
        allocate(pr%teav(0:nr, nl), stat=status(9))
        allocate(pr%osav(0:nr, nl), stat=status(10))
        allocate(pr%phsav(0:nr, nl), stat=status(11))
        allocate(pr%cmn(0:nr,nv, nl), stat=status(12))
        allocate(pr%cmx(0:nr,nv, nl), stat=status(13))
        allocate(pr%cav(0:nr,nv, nl), stat=status(14))
        allocate(pr%phmn(0:nr, nl), stat=status(15))
        allocate(pr%phmx(0:nr, nl), stat=status(16))
        allocate(pr%phav(0:nr, nl), stat=status(17))
        allocate(pr%tnmn(0:nr, nl), stat=status(18))
        allocate(pr%tnmx(0:nr, nl), stat=status(19))
        allocate(pr%tnav(0:nr, nl), stat=status(20))
        allocate(pr%tpmn(0:nr, nl), stat=status(21))
        allocate(pr%tpmx(0:nr, nl), stat=status(22))
        allocate(pr%tpav(0:nr, nl), stat=status(23))
        allocate(pr%nh3mn(0:nr, nl), stat=status(24))
        allocate(pr%nh3mx(0:nr, nl), stat=status(25))
        allocate(pr%nh3av(0:nr, nl), stat=status(26)) !gp end new block

        !gp 20-oct-04 growth limitation factors for bottom algae
        allocate(pr%phitotalsavepr(0:nr, 0:nsteps), stat=status(27))
        allocate(pr%phitsavepr(0:nr, 0:nsteps), stat=status(28))
        allocate(pr%philsavepr(0:nr, 0:nsteps), stat=status(29))
        allocate(pr%phinsavepr(0:nr, 0:nsteps), stat=status(30))
        allocate(pr%phipsavepr(0:nr, 0:nsteps), stat=status(31))
        allocate(pr%phicsavepr(0:nr, 0:nsteps), stat=status(32)) !gp 20-oct-04 end new block

        !gp 28-oct-04 diagenesis flux between sediment/water
        allocate(pr%diagfluxdopr(0:nr, 0:nsteps), stat=status(33))
        allocate(pr%diagfluxcbodpr(0:nr, 0:nsteps), stat=status(34))
        allocate(pr%diagfluxnh4pr(0:nr, 0:nsteps), stat=status(35))
        allocate(pr%diagfluxno3pr(0:nr, 0:nsteps), stat=status(36))
        allocate(pr%diagfluxsrppr(0:nr, 0:nsteps), stat=status(37))
        allocate(pr%diagfluxicpr(0:nr, 0:nsteps), stat=status(38)) !gp end new block

        !gp 28-oct-04 diagenesis flux between sediment/water
        allocate(pr%hypofluxdopr(0:nr, 0:nsteps), stat=status(39))
        allocate(pr%hypofluxcbodpr(0:nr, 0:nsteps), stat=status(40))
        allocate(pr%hypofluxnh4pr(0:nr, 0:nsteps), stat=status(41))
        allocate(pr%hypofluxno3pr(0:nr, 0:nsteps), stat=status(42))
        allocate(pr%hypofluxsrppr(0:nr, 0:nsteps), stat=status(43))
        allocate(pr%hypofluxicpr(0:nr, 0:nsteps), stat=status(44)) !gp end new block

        !gp 15-nov-04 reach-average daily-average flux between sediment/water
        allocate(pr%diagfluxdoav(0:nr), stat=status(45))
        allocate(pr%diagfluxcbodav(0:nr), stat=status(46))
        allocate(pr%diagfluxnh4av(0:nr), stat=status(47))
        allocate(pr%diagfluxno3av(0:nr), stat=status(48))
        allocate(pr%diagfluxsrpav(0:nr), stat=status(49))
        allocate(pr%hypofluxdoav(0:nr), stat=status(50))
        allocate(pr%hypofluxcbodav(0:nr), stat=status(51))
        allocate(pr%hypofluxnh4av(0:nr), stat=status(52))
        allocate(pr%hypofluxno3av(0:nr), stat=status(53))
        allocate(pr%hypofluxsrpav(0:nr), stat=status(54))

        !gp 11-jan-05 cell quota mgn/gd and mgp/gd
        allocate(pr%ninbmn(0:nr), stat=status(55))
        allocate(pr%ninbmx(0:nr), stat=status(56))
        allocate(pr%ninbav(0:nr), stat=status(57))
        allocate(pr%nipbmn(0:nr), stat=status(58))
        allocate(pr%nipbmx(0:nr), stat=status(59))
        allocate(pr%nipbav(0:nr), stat=status(60))

        !gp 05-jul-05 heat/do/co2 fluxes
        allocate(pr%pr_saveheatfluxjsnt(0:nr, 0:nsteps), stat=status(61))
        allocate(pr%pr_saveheatfluxlongat(0:nr, 0:nsteps), stat=status(62))
        allocate(pr%pr_saveheatfluxback(0:nr, 0:nsteps), stat=status(63))
        allocate(pr%pr_saveheatfluxconv(0:nr, 0:nsteps), stat=status(64))
        allocate(pr%pr_saveheatfluxevap(0:nr, 0:nsteps), stat=status(65))
        allocate(pr%pr_saveheatfluxjsed(0:nr, 0:nsteps), stat=status(66))
        allocate(pr%pr_saveheatfluxjhyporheic(0:nr, 0:nsteps), stat=status(67))
        allocate(pr%pr_saveheatfluxtribs(0:nr, 0:nsteps), stat=status(68))
        allocate(pr%pr_saveheatfluxadvecdisp(0:nr, 0:nsteps), stat=status(69))
        allocate(pr%pr_savedofluxreaer(0:nr, 0:nsteps), stat=status(70))
        allocate(pr%pr_savedofluxcbodfast(0:nr, 0:nsteps), stat=status(71))
        allocate(pr%pr_savedofluxcbodslow(0:nr, 0:nsteps), stat=status(72))
        allocate(pr%pr_savedofluxcod(0:nr, 0:nsteps), stat=status(73))
        allocate(pr%pr_savedofluxnitrif(0:nr, 0:nsteps), stat=status(74))
        allocate(pr%pr_savedofluxphytoresp(0:nr, 0:nsteps), stat=status(75))
        allocate(pr%pr_savedofluxphytophoto(0:nr, 0:nsteps), stat=status(76))
        allocate(pr%pr_savedofluxbotalgresp(0:nr, 0:nsteps), stat=status(77))
        allocate(pr%pr_savedofluxbotalgphoto(0:nr, 0:nsteps), stat=status(78))
        allocate(pr%pr_savedofluxsod(0:nr, 0:nsteps), stat=status(79))
        allocate(pr%pr_savedofluxhyporheic(0:nr, 0:nsteps), stat=status(80))
        allocate(pr%pr_savedofluxtribs(0:nr, 0:nsteps), stat=status(81))
        allocate(pr%pr_savedofluxadvecdisp(0:nr, 0:nsteps), stat=status(82))
        allocate(pr%pr_saveco2fluxreaer(0:nr, 0:nsteps), stat=status(83))
        allocate(pr%pr_saveco2fluxcbodfast(0:nr, 0:nsteps), stat=status(84))
        allocate(pr%pr_saveco2fluxcbodslow(0:nr, 0:nsteps), stat=status(85))
        allocate(pr%pr_saveco2fluxphytoresp(0:nr, 0:nsteps), stat=status(86))
        allocate(pr%pr_saveco2fluxphytophoto(0:nr, 0:nsteps), stat=status(87))
        allocate(pr%pr_saveco2fluxbotalgresp(0:nr, 0:nsteps), stat=status(88))
        allocate(pr%pr_saveco2fluxbotalgphoto(0:nr, 0:nsteps), stat=status(89))
        allocate(pr%pr_saveco2fluxsod(0:nr, 0:nsteps), stat=status(90))
        allocate(pr%pr_saveco2fluxhyporheic(0:nr, 0:nsteps), stat=status(91))
        allocate(pr%pr_saveco2fluxtribs(0:nr, 0:nsteps), stat=status(92))
        allocate(pr%pr_saveco2fluxadvecdisp(0:nr, 0:nsteps), stat=status(93))

        !gp 25-jun-09
        allocate(pr%av_botalgphoto(0:nr), stat=status(94))
        allocate(pr%av_botalgresp(0:nr), stat=status(95))
        allocate(pr%av_botalgdeath(0:nr), stat=status(96))
        allocate(pr%av_botalgnetgrowth(0:nr), stat=status(97))

        pr%temn=0; pr%temx=0; pr%teav=0
        pr%osav=0; pr%phsav=0; pr%cmn=0
        pr%cmx=0; pr%cav=0; pr%phmn=0
        pr%phmx=0; pr%phav=0; pr%tnmn=0
        pr%tnmx=0; pr%tnav=0; pr%tpmn=0
        pr%tpmx=0; pr%tpav=0; pr%nh3mn=0
        pr%nh3mx=0; pr%nh3av=0;
        pr%cpr=0; pr%tepr=0; pr%phpr=0
        pr%ninbpr=0; pr%nipbpr=0

        !gp 20-oct-04
        pr%phitotalsavepr=0; pr%phitsavepr=0; pr%philsavepr=0; pr%phinsavepr=0; pr%phipsavepr=0; pr%phicsavepr=0

        !gp 28-oct-04 sed fluxes at each calc step
        pr%diagfluxdopr=0; pr%diagfluxcbodpr=0; pr%diagfluxnh4pr=0; pr%diagfluxno3pr=0; pr%diagfluxsrppr=0; pr%diagfluxicpr=0
        pr%hypofluxdopr=0; pr%hypofluxcbodpr=0; pr%hypofluxnh4pr=0; pr%hypofluxno3pr=0; pr%hypofluxsrppr=0; pr%hypofluxicpr=0

        !gp 15-nov-04 reach-average daily average sed fluxes
        pr%diagfluxdoav=0; pr%diagfluxcbodav=0; pr%diagfluxnh4av=0; pr%diagfluxno3av=0; pr%diagfluxsrpav=0
        pr%hypofluxdoav=0; pr%hypofluxcbodav=0; pr%hypofluxnh4av=0; pr%hypofluxno3av=0; pr%hypofluxsrpav=0

        !gp 11-jan-05 cell quota mgn/gd and mgp/gd
        pr%ninbmn=0; pr%ninbmx=0; pr%ninbav=0; pr%nipbmn=0; pr%nipbmx=0; pr%nipbav=0


        !gp 05-jul-05
        pr%pr_saveheatfluxjsnt=0; pr%pr_saveheatfluxlongat=0; pr%pr_saveheatfluxback=0; pr%pr_saveheatfluxconv=0
        pr%pr_saveheatfluxevap=0; pr%pr_saveheatfluxjsed=0; pr%pr_saveheatfluxjhyporheic=0; pr%pr_saveheatfluxtribs=0
        pr%pr_saveheatfluxadvecdisp=0
        pr%pr_savedofluxreaer=0; pr%pr_savedofluxcbodfast=0; pr%pr_savedofluxcbodslow=0; pr%pr_savedofluxnitrif=0
        pr%pr_savedofluxphytoresp=0; pr%pr_savedofluxphytophoto=0; pr%pr_savedofluxbotalgresp=0; pr%pr_savedofluxbotalgphoto=0
        pr%pr_savedofluxsod=0; pr%pr_savedofluxcod=0; pr%pr_savedofluxhyporheic=0; pr%pr_savedofluxtribs=0
        pr%pr_savedofluxadvecdisp=0
        pr%pr_saveco2fluxreaer=0; pr%pr_saveco2fluxcbodfast=0; pr%pr_saveco2fluxcbodslow=0
        pr%pr_saveco2fluxphytoresp=0; pr%pr_saveco2fluxphytophoto=0; pr%pr_saveco2fluxbotalgresp=0; pr%pr_saveco2fluxbotalgphoto=0
        pr%pr_saveco2fluxsod=0; pr%pr_saveco2fluxhyporheic=0; pr%pr_saveco2fluxtribs=0
        pr%pr_saveco2fluxadvecdisp=0

        !gp 25-jun-09
        pr%av_botalgphoto=0; pr%av_botalgresp=0; pr%av_botalgdeath=0; pr%av_botalgnetgrowth=0

        !gp 15-nov-04 do i=0, 26
        !gp 05-jul-05 do i=0, 60
        do i=0, 93
            if (status(i)==1) then
                write(8,*) '** Class_Output:outData_ failed. Insufficient memory for dynamic diel output arrays. **'
                close (8) !gp 17-nov-04
                stop !class_integration:integration_ failed. insufficient memory!'
            end if

            !gp debug
            !open (unit=9, file='debug2.out', status='replace', action='write')
            !write(9,*) 'status(', i, ') =', status(i)
            !close (9)

        end do

    end function outdata_ctor

    subroutine output(pr, nr, topo, hydrau, rates, system)

        type(outdata_t), intent(in) :: pr
        type(t_rivertopo) topo
        type(riverhydraulics_type), intent(in) :: hydrau
        type(rates_t) rates
        type(systemparams) system
        integer(i32), intent(in) :: nr

        !gp integer(i32) i, j, nrp
        integer(i32) i, j, nrp, k !gp
        real(r64) toc, tkn, tss, tp, tn, bottomalgae, dosat, nh3
        integer(i32) ihour
        real(r64) t, ph, cbodu
        !gp output of hydraulics
        !sheets("hydraulics summary")

        call hydraulics_summary(topo, hydrau, rates, nr)

        !output of hourly summary of loads
        !sheets("source summary")
        write(8,*)
        write(8,*) '** Source summary **'
        !gp 30-Nov-04
        !gp write(8,'(A6, A2, 2A15, 20A12)') 'Time', '', 'Reach', 'Downstream', 'UpDist', 'Down Dist', 'Abstraction', 'Inflow', &
        !gp 'Temp', 'Cond', 'Iss', 'Oxygen', 'CBODs', 'CBODf', 'No', 'NH4', &
        !gp 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogens', 'Alk', 'pH'
        !gp write(8,'(A8, 2A15, 20A12)') '', 'Label', 'Label', 'x(km)', 'x(km)', 'cms', 'cms', 'C', 'umhos', 'mgD/L', &
        !gp 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', 'ugN/L', 'ugN/L', &
        !gp 'ugN/L', 'ugN/L', 'ugN/L', 'mgD/L', 'cfu/100mL', 'mgCaCO3/L'
        write(8,'(A6, A2, 2A15, 21A12)') 'Time', '', 'Reach', 'Downstream', 'UpDist', 'Down Dist', 'Abstraction', 'Inflow', &
            'Temp', 'Cond', 'Iss', 'Oxygen', 'CBODs', 'CBODf', 'No', 'NH4', &
            'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogens', 'Generic', 'Alk', 'pH'
        write(8,'(A8, 2A15, 21A12)') '', 'Label', 'Label', 'x(km)', 'x(km)', 'cms', 'cms', 'C', 'umhos', 'mgD/L', &
            'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', 'ugN/L', 'ugN/L', &
            'ugN/L', 'ugN/L', 'ugN/L', 'mgD/L', 'cfu/100mL', 'user define', 'mgCaCO3/L' !gp 30-Nov-04 end new block
        if (npt /= 0 .or. ndiff /= 0) then
            ! rlab2(0) = rlab1(1)
            do ihour = 0, 23 !loop through output of hourly sources
                t = ihour / 24.0_r64 !current output time in days
                !evaluate sine functions and distribute loads to reaches at time t
                call sourcescalc(t, nr, hydrau%flag)

                !output of simulation date+time
                !activecell.value = dateserial(xyear, xmon, xday) + t

                do i = 1, nr

                    if (load(i)%c(nv - 2) > 0 .and. load(i)%c(nv - 1) > 0) then !gp 03-dec-04
                        call ph_solver(system%imethph, &
                            ph, load(i)%c(nv - 1), load(i)%te, load(i)%c(nv - 2), load(i)%c(1))
                    else
                        ph=0.0
                    end if
                    !gp 30-nov-04
                    !gp write(8,'(f6.5, a2, 2a15, 20f12.4)') t,'', topo%reach(i)%rname, topo%reach(i)%rlab, topo%reach(i-1)%xrdn, &
                    !gp topo%reach(i)%xrdn, hydrau%reach(i)%qpta, hydrau%reach(i)%qpt, &
                    !gp load(i)%te, (load(i)%c(j), j=1, nv-2), ph
                    write(8,'(F6.5, A2, 2A15, 21F12.4)') t,'', topo%reach(i)%rname, topo%reach(i)%rlab, topo%reach(i-1)%xrdn, &
                        topo%reach(i)%xrdn, hydrau%reach(i)%qpta, hydrau%reach(i)%qpt, &
                        load(i)%te, (load(i)%c(j), j=1, nv-2), ph !gp 30-nov-04 add generic constituent
                end do
            end do
        end if

! output temperature for the water column
! sheets("temperature output")
        write(8,*)
        write(8,*) '** Temperature summary (water column temperature) **'
        write(8,'(A5, A10, 4A10)') 'Reach', '', 'Distance','Temp(C)', 'Temp(C)', 'Temp(C)'
        write(8,'(A5, A10, 4A10)') 'Label', '', 'x(km)','Average', 'Minimum', 'Maximum'
        j = 1 !gp 27-oct-04 water column is layer 1
        do i=0, nr
            !gp 27-oct-04 write (8,'(a15, 4f10.4)') topo%reach(i)%rname, topo%reach(i)%xpm, pr%teav(i), &
            !gp pr%temn(i), pr%temx(i)
            write (8,'(A15, 4F10.4)') topo%reach(i)%rname, topo%reach(i)%xpm, pr%teav(i, j), &
                pr%temn(i, j), pr%temx(i, j) !gp add nl dimension
        end do

! output concentrations for the water column
        write(8,*)
        write(8,*) '** Daily average water quality summary (water column constituents) **'
        !gp 30-nov-04
        !gp write(8,'(a5, a10, 27a13)') 'reach', '', 'x', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'no', 'nh4', 'no3', &
        !gp 'po', 'inorgp', 'phyto', 'detritus', 'pathogen', 'alk', 'ph', &
        !gp 'bot alg', 'toc', 'tn', 'tp', 'tkn', 'tss', 'cbodu', 'bot alg', &
        !gpwrite(8,'(a5, a10, 27a13)') 'label', '', 'km', 'umhos', 'mgd/l', 'mgo2/l', 'mgo2/l', 'mgo2/l', 'ugn/l', &
        !gp 'ugn/l', 'ugn/l', &
        !gp 'ugp/l', 'ugp/l', 'uga/l', 'mgd/l', '', 'alk', 'ph', &
        !gp 'gd/m2', '', '', '', '', 'mgd/l', '', 'mga/m2', &
        !gp '', '', ''

        !gp 25-jun-09
        !write(8,'(a5, a10, 33a13)') 'reach', '', 'x', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'norg', 'nh4', 'no3', &
        ! 'porg', 'inorgp', 'phyto', 'detritus', 'pathogen', 'generic', 'alk', 'ph', &
        ! 'bot alg', 'toc', 'tn', 'tp', 'tkn', 'tss', 'cbodu', 'bot alg', &
        ! 'nh3', 'do sat', 'ph sat', 'hypo biofilm', &
        ! 'ninbav', 'nipbav', 'ninbav', 'nipbav'
        !write(8,'(a5, a10, 33a13)') 'label', '', 'km', 'umhos', 'mgd/l', 'mgo2/l', 'mgo2/l', 'mgo2/l', 'ugn/l', &
        ! 'ugn/l', 'ugn/l', &
        ! 'ugp/l', 'ugp/l', 'uga/l', 'mgd/l', '', 'user defined', 'mgcaco3/l', 's.u.', &
        ! 'gd/m2', '', '', '', '', 'mgd/l', '', 'mga/m2', &
        ! '', '', '', 'gd/m2', 'mgn/mga', 'mgp/mga', 'mgn/gd', 'mgp/gd' !gp 11-jan-05 end new block
        write(8,'(A5, A10, 41A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
            'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
            'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
            'NH3', 'DO sat', 'pH sat', 'Hypo biofilm', &
            'NINbav', 'NIPbav', 'NINbav', 'NIPbav', &
            'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow', &
            'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow'
        write(8,'(A5, A10, 41A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
            'ugN/L', 'ugN/L', &
            'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
            'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', &
            '', '', '', 'gD/m2', 'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD', &
            'gD/m2/d', 'gD/m2/d', 'gD/m2/d', 'gD/m2/d', &
            'mgA/m2/d', 'mgA/m2/d', 'mgA/m2/d', 'mgA/m2/d'

        do i = 0, nr
            !gp 27-oct-04 add nl dimension
            !gptoc = (pr%cav(i, 4) + pr%cav(i, 5)) / rates%roc + &
            !gp rates%aca * pr%cav(i, 11) + rates%aca / rates%ada * pr%cav(i, 12)
            !gptkn = pr%cav(i, 6) + pr%cav(i, 7) + rates%ana * pr%cav(i, 11)
            !gptss = rates%ada * pr%cav(i, 11) + pr%cav(i, 2) + pr%cav(i, 12)
            !gpcbodu = pr%cav(i, 4) + pr%cav(i, 5) + rates%roa * pr%cav(i, 11) + &
            !gp rates%roc * rates%aca / rates%ada * pr%cav(i, 12)
            !gp!bottom algae as chl a
            !gpbottomalgae= pr%cav(i, 16) / (rates%adc * rates%aca)
            !gpwrite(8,'(a15, 27f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, j), j=1, nv-2), pr%phav(i) , &
            !gp pr%cav(i, nv), toc, pr%tnav(i), pr%tpav(i), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3av(i), pr%osav(i), &
            !gp pr%phsav(i)
            j = 1 !gp water column is layer 1
            toc = (pr%cav(i, 4, j) + pr%cav(i, 5, j)) / rates%roc + &
                rates%aca * pr%cav(i, 11, j) + rates%aca / rates%ada * pr%cav(i, 12, j)
            tkn = pr%cav(i, 6, j) + pr%cav(i, 7, j) + rates%ana * pr%cav(i, 11, j)
            tss = rates%ada * pr%cav(i, 11, j) + pr%cav(i, 2, j) + pr%cav(i, 12, j)
            cbodu = pr%cav(i, 4, j) + pr%cav(i, 5, j) + rates%roa * pr%cav(i, 11, j) + &
                rates%roc * rates%aca / rates%ada * pr%cav(i, 12, j)
            !bottom algae as chl a
            bottomalgae= pr%cav(i, nv, j) / (rates%adc * rates%aca)
            !gp 30-nov-04
            !gp write(8,'(a15, 27f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
            !gp pr%cav(i, nv, j), toc, pr%tnav(i, j), pr%tpav(i, j), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3av(i, j), pr%osav(i, j), &
            !gp pr%phsav(i, j)

            !gp 25-jun-09
            !if (i == 0) then !make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
            ! write(8,'(a15, 33f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
            ! pr%cav(1, nv, j), toc, pr%tnav(i, j), pr%tpav(i, j), &
            ! tkn, tss, cbodu, pr%cav(1, nv, j) / (rates%adc * rates%aca), pr%nh3av(i, j), pr%osav(i, j), &
            ! pr%phsav(i, j), pr%cav(1, nv, 2), &
            ! pr%ninbav(1) * rates%mgd / rates%mga / 1000, pr%nipbav(1) * rates%mgd / rates%mga / 1000, &
            ! pr%ninbav(1), pr%nipbav(1) !gp 11-jan-05 end new block
            !else
            ! write(8,'(a15, 33f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
            ! pr%cav(i, nv, j), toc, pr%tnav(i, j), pr%tpav(i, j), &
            ! tkn, tss, cbodu, bottomalgae, pr%nh3av(i, j), pr%osav(i, j), &
            ! pr%phsav(i, j), pr%cav(i, nv, 2), &
            ! pr%ninbav(i) * rates%mgd / rates%mga / 1000, pr%nipbav(i) * rates%mgd / rates%mga / 1000, &
            ! pr%ninbav(i), pr%nipbav(i) !gp 11-jan-05 end new block
            !end if
            if (i == 0) then !make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
                write(8,'(A15, 41F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
                    pr%cav(1, nv, j), toc, pr%tnav(i, j), pr%tpav(i, j), &
                    tkn, tss, cbodu, pr%cav(1, nv, j) / (rates%adc * rates%aca), pr%nh3av(i, j), pr%osav(i, j), &
                    pr%phsav(i, j), pr%cav(1, nv, 2), &
                    pr%ninbav(1) * rates%mgd / rates%mga / 1000, pr%nipbav(1) * rates%mgd / rates%mga / 1000, &
                    pr%ninbav(1), pr%nipbav(1), &
                    pr%av_botalgphoto(1), pr%av_botalgresp(1), pr%av_botalgdeath(1), pr%av_botalgnetgrowth(1), &
                    pr%av_botalgphoto(1)/rates%ada, pr%av_botalgresp(1)/rates%ada, &
                    pr%av_botalgdeath(1)/rates%ada, pr%av_botalgnetgrowth(1)/rates%ada
            else
                write(8,'(A15, 41F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
                    pr%cav(i, nv, j), toc, pr%tnav(i, j), pr%tpav(i, j), &
                    tkn, tss, cbodu, bottomalgae, pr%nh3av(i, j), pr%osav(i, j), &
                    pr%phsav(i, j), pr%cav(i, nv, 2), &
                    pr%ninbav(i) * rates%mgd / rates%mga / 1000, pr%nipbav(i) * rates%mgd / rates%mga / 1000, &
                    pr%ninbav(i), pr%nipbav(i), &
                    pr%av_botalgphoto(i), pr%av_botalgresp(i), pr%av_botalgdeath(i), pr%av_botalgnetgrowth(i), &
                    pr%av_botalgphoto(i)/rates%ada, pr%av_botalgresp(i)/rates%ada, &
                    pr%av_botalgdeath(i)/rates%ada, pr%av_botalgnetgrowth(i)/rates%ada
            end if

        end do

! output min concentrations
        write(8,*)
        write(8,*) '** Daily minimum water quality summary (water column constituents) **'
        !gp 15-nov-04
        !gp write(8,'(a5, a10, 27a13)') 'reach', '', 'x', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'no', 'nh4', 'no3', &
        !gp 'po', 'inorgp', 'phyto', 'detritus', 'pathogen', 'alk', 'ph', &
        !gp 'bot alg', 'toc', 'tn', 'tp', 'tkn', 'tss', 'cbodu', 'bot alg', &
        !gp 'nh3'
        !gp write(8,'(a5, a10, 27a13)') 'label', '', 'km', 'umhos', 'mgd/l', 'mgo2/l', 'mgo2/l', 'mgo2/l', 'ugn/l', &
        !gp 'ugn/l', 'ugn/l', &
        !gp 'ugp/l', 'ugp/l', 'uga/l', 'mgd/l', '', 'alk', 'ph', &
        !gp 'gd/m2', '', '', '', '', 'mgd/l', '', 'mga/m2', ''
        !write(8,'(a5, a10, 31a13)') 'reach', '', 'x', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'no', 'nh4', 'no3', &
        write(8,'(A5, A10, 31A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
            'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
            'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
            'NH3', 'Hypo biofilm', 'NINbmn', 'NIPbmn', 'NINbmn', 'NIPbmn'
        !write(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
        write(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
            'ugN/L', 'ugN/L', &
            'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'Alk', 'pH', &
            'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', '', 'gD/m2', &
            'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD' !30-Nov-04 add hypo biofilm and generic const
        do i = 0, nr
            !gp 27-oct-04 add dimension for nl
            !gp toc = (pr%cmn(i, 4) + pr%cmn(i, 5)) / rates%roc + &
            !gp rates%aca * pr%cmn(i, 11) + rates%aca / rates%ada * pr%cmn(i, 12)
            !gp tkn = pr%cmn(i, 6) + pr%cmn(i, 7) + rates%ana * pr%cmn(i, 11)
            !gp tss = rates%ada * pr%cmn(i, 11) + pr%cmn(i, 2) + pr%cmn(i, 12)
            !gp cbodu = pr%cmn(i, 4) + pr%cmn(i, 5) + rates%roa * pr%cmn(i, 11) + &
            !gp rates%roc * rates%aca / rates%ada * pr%cmn(i, 12)
            !gp !bottom algae as chl a
            !gp bottomalgae= pr%cmn(i, 16) / (rates%adc * rates%aca)
            !gp write(8,'(a15, 25f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmn(i, j), j=1, nv-2), pr%phmn(i) , &
            !gp pr%cmn(i, nv), toc, pr%tnmn(i), pr%tpmn(i), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3mn(i)
            j = 1 !gp water column is layer 1
            toc = (pr%cmn(i, 4, j) + pr%cmn(i, 5, j)) / rates%roc + &
                rates%aca * pr%cmn(i, 11, j) + rates%aca / rates%ada * pr%cmn(i, 12, j)
            tkn = pr%cmn(i, 6, j) + pr%cmn(i, 7, j) + rates%ana * pr%cmn(i, 11, j)
            tss = rates%ada * pr%cmn(i, 11, j) + pr%cmn(i, 2, j) + pr%cmn(i, 12, j)
            cbodu = pr%cmn(i, 4, j) + pr%cmn(i, 5, j) + rates%roa * pr%cmn(i, 11, j) + &
                rates%roc * rates%aca / rates%ada * pr%cmn(i, 12, j)
            !bottom algae as chl a
            bottomalgae= pr%cmn(i, nv, j) / (rates%adc * rates%aca)
            !gp 30-nov-04
            !gp write(8,'(a15, 25f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%phmn(i, j) , &
            !gp pr%cmn(i, nv, j), toc, pr%tnmn(i, j), pr%tpmn(i, j), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3mn(i, j)
            if (i == 0) then !make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
                write(8,'(A15, 31F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%phmn(i, j) , &
                    pr%cmn(1, nv, j), toc, pr%tnmn(i, j), pr%tpmn(i, j), &
                    tkn, tss, cbodu, pr%cmn(1, nv, j) / (rates%adc * rates%aca), pr%nh3mn(i, j), pr%cmn(1, nv, 2), & !gp 15-nov-04 end new block
                    pr%ninbmn(1) * rates%mgd / rates%mga / 1000, pr%nipbmn(1) * rates%mgd / rates%mga / 1000, &
                    pr%ninbmn(1), pr%nipbmn(1) !gp 11-jan-05 end new block
            else
                write(8,'(A15, 31F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%phmn(i, j) , &
                    pr%cmn(i, nv, j), toc, pr%tnmn(i, j), pr%tpmn(i, j), &
                    tkn, tss, cbodu, bottomalgae, pr%nh3mn(i, j), pr%cmn(i, nv, 2), &
                    pr%ninbmn(i) * rates%mgd / rates%mga / 1000, pr%nipbmn(i) * rates%mgd / rates%mga / 1000, &
                    pr%ninbmn(i), pr%nipbmn(i) !gp 11-jan-05 end new block
            end if
        end do

! Output max concentrations
        write(8,*)
        write(8,*) '** Daily maximum water quality summary (water column constituents) **'
        !gp 30-Nov-04
        !gp write(8,'(A5, A10, 27A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
        !gp 'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Alk', 'pH', &
        !gp 'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
        !gp 'NH3'
        !gp write(8,'(A5, A10, 27A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
        !gp 'ugN/L', 'ugN/L', &
        !gp 'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'Alk', 'pH', &
        !gp 'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', ''
        write(8,'(A5, A10, 31A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
            'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
            'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
            'NH3', 'Hypo biofilm', 'NINbmx', 'NIPbmx', 'NINbmx', 'NIPbmx'
        write(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
            'ugN/L', 'ugN/L', &
            'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
            'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', '', 'gD/m2', &
            'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD' !gp 11-Jan-05 end new block
        do i = 0, nr
            !gp 27-oct-04 add dimension for nl
            !gp toc = (pr%cmx(i, 4) + pr%cmx(i, 5)) / rates%roc + &
            !gp rates%aca * pr%cmx(i, 11) + rates%aca / rates%ada * pr%cmx(i, 12)
            !gp tkn = pr%cmx(i, 6) + pr%cmx(i, 7) + rates%ana * pr%cmx(i, 11)
            !gp tss = rates%ada * pr%cmx(i, 11) + pr%cmx(i, 2) + pr%cmx(i, 12)
            !gp cbodu = pr%cmx(i, 4) + pr%cmx(i, 5) + rates%roa * pr%cmx(i, 11) + &
            !gp rates%roc * rates%aca / rates%ada * pr%cmx(i, 12)
            !gp !bottom algae as chl a
            !gp bottomalgae= pr%cmx(i, 16) / (rates%adc * rates%aca)
            !gp write(8,'(a15, 25f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmx(i, j), j=1, nv-2), pr%phmx(i) , &
            !gp pr%cmx(i, nv), toc, pr%tnmx(i), pr%tpmx(i), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3mx(i)
            j = 1 !gp water column is layer 1
            toc = (pr%cmx(i, 4, j) + pr%cmx(i, 5, j)) / rates%roc + &
                rates%aca * pr%cmx(i, 11, j) + rates%aca / rates%ada * pr%cmx(i, 12, j)
            tkn = pr%cmx(i, 6, j) + pr%cmx(i, 7, j) + rates%ana * pr%cmx(i, 11, j)
            tss = rates%ada * pr%cmx(i, 11, j) + pr%cmx(i, 2, j) + pr%cmx(i, 12, j)
            cbodu = pr%cmx(i, 4, j) + pr%cmx(i, 5, j) + rates%roa * pr%cmx(i, 11, j) + &
                rates%roc * rates%aca / rates%ada * pr%cmx(i, 12, j)
            !bottom algae as chl a
            bottomalgae= pr%cmx(i, nv, j) / (rates%adc * rates%aca)
            !gp 30-nov-04
            !gp write(8,'(a15, 25f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%phmx(i, j) , &
            !gp pr%cmx(i, nv, j), toc, pr%tnmx(i, j), pr%tpmx(i, j), &
            !gp tkn, tss, cbodu, bottomalgae, pr%nh3mx(i, j)
            if (i == 0) then !make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
                write(8,'(A15, 31F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%phmx(i, j) , &
                    pr%cmx(1, nv, j), toc, pr%tnmx(i, j), pr%tpmx(i, j), &
                    tkn, tss, cbodu, pr%cmx(1, nv, j) / (rates%adc * rates%aca), pr%nh3mx(i, j), pr%cmx(1, nv, 2), &
                    pr%ninbmx(1) * rates%mgd / rates%mga / 1000, pr%nipbmx(1) * rates%mgd / rates%mga / 1000, &
                    pr%ninbmx(1), pr%nipbmx(1) !gp 11-jan-05 end new block
            else
                write(8,'(A15, 31F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%phmx(i, j), &
                    pr%cmx(i, nv, j), toc, pr%tnmx(i, j), pr%tpmx(i, j), &
                    tkn, tss, cbodu, bottomalgae, pr%nh3mx(i, j), pr%cmx(i, nv, 2), &
                    pr%ninbmx(i) * rates%mgd / rates%mga / 1000, pr%nipbmx(i) * rates%mgd / rates%mga / 1000, &
                    pr%ninbmx(i), pr%nipbmx(i) !gp 11-jan-05 end new block
            end if
        end do

! output sediment fluxes !gp 15-nov-04 reach-averaged and daily-averaged and include hyporheic and total flux
        write(8,*)
        write(8,*) '** Sediment fluxes (reach-average daily-average) **'
        !gp 15-nov-04
        !gp write(8,'(a5, a10, 6a13)') 'reach', '', 'distance', 'sod', 'flux ch4', 'flux nh4', &
        !gp 'flux inorgp', 'flux no3'
        !gp write(8,'(a5, a10, 6a13)') 'label', '', 'x(km)', 'go2/m2/d', 'go2/m2/d', 'mgn/m2/d', &
        !gp 'mgp/m2/d', 'mgn/m2/d'
        write(8,'(A5, A10, 16A13)') 'Reach', '', 'Distance', &
            'DiagFluxDO', 'DiagFluxCBOD', 'DiagFluxNH4', 'DiagFluxSRP', 'DiagFluxNO3', &
            'HypoFluxDO', 'HypoFluxCBOD', 'HypoFluxNH4', 'HypoFluxSRP', 'HypoFluxNO3', &
            'TotFluxDO', 'TotFluxCBOD', 'TotFluxNH4', 'TotFluxSRP', 'TotFluxNO3'
        write(8,'(A5, A10, 16A13)') 'Label', '', 'x(km)', &
            'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d', &
            'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d', &
            'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d' !gp 15-Nov-04 end new block
        do i=1, nr
            !gp 15-nov-04 reach-average daily-average sediment fluxes
            !gp write(8,'(a15, 6f13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, sodpr(i), &
            !gp jch4pr(i), jnh4pr(i), jsrppr(i), jno3pr(i)
            write(8,'(A15, 16F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, &
                pr%diagfluxdoav(i), pr%diagfluxcbodav(i), pr%diagfluxnh4av(i), pr%diagfluxsrpav(i), pr%diagfluxno3av(i), &
                pr%hypofluxdoav(i), pr%hypofluxcbodav(i), pr%hypofluxnh4av(i), pr%hypofluxsrpav(i), pr%hypofluxno3av(i), &
                pr%diagfluxdoav(i) + pr%hypofluxdoav(i), &
                pr%diagfluxcbodav(i) + pr%hypofluxcbodav(i), &
                pr%diagfluxnh4av(i) + pr%hypofluxnh4av(i), &
                pr%diagfluxsrpav(i) + pr%hypofluxsrpav(i), &
                pr%diagfluxno3av(i) + pr%hypofluxno3av(i) !gp 15-nov-04 end new block
        end do


!gp 17-feb-05
!only output diel if showdielresults = "yes"
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        if (system%showdielresults == "Yes") then

!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

!diel output result for water column
            write(8,*)
            write(8,*) '** Diel water quality in the water column **'

            !gp 01-nov-04 write(8,'(32a13)') 't', 'tempw', ' temps', 'cond', 'iss', 'do', 'cbods', 'cbodf', 'no', &
            write(8,'(32A13)') 't', 'Tempw', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', & !gp 30-Nov-04
                'NH4', 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', & !gp 30-Nov-04
                'Bot Algae', 'TSS', 'TP', 'TN', 'DOsat', 'NH3', 'IntN', 'Int P', &
                'phiTemp', 'phiLight', 'phiNitr', 'phiPhos', 'phiCarb', 'phiTotal' !gp 20-Oct-04

            !gp 01-Nov-04 write(8,'(32A13)') 'hr', 'c', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', &
            write(8,'(32A13)') 'hr', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', & !gp 32-Nov-04
                'ugN/L', 'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', '', '', '', '', 'gD/m2', 'mgD/L', & !gp 30-Nov-04
                'ugP/L', 'ugN/L', 'mg/L', 'ugN/L', 'mgN/mgA', 'mgP/mgA', &
                'frac', 'frac', 'frac', 'frac', 'frac', 'frac' !gp 20-Oct-04

            !gp write(8,*) pr%nj
            write(8,'(I13)') pr%nj
            do nrp=0, nr
                !tdy=0
                do i=0, pr%nj
                    !gp 27-oct-04 add dimension for nl
                    !gp tss = pr%cpr(nrp, 11, i) * rates%ada + pr%cpr(nrp, 2, i) + &
                    !gp pr%cpr(nrp, 12, i)
                    !gp tp = pr%cpr(nrp, 11, i) * rates%apa + pr%cpr(nrp, 9, i) + &
                    !gp pr%cpr(nrp, 10, i)
                    !gp tn = pr%cpr(nrp, 11, i) * rates%ana + pr%cpr(nrp, 6, i) + &
                    !gp pr%cpr(nrp, 7, i) + pr%cpr(nrp, 8, i)
                    !gp dosat = oxsat(pr%tepr(nrp, i, 1), hydrau%reach(nrp)%elev)
                    !gp nh3 = 1.0_r64/(1 + 10.0_r64 ** (-pr%phpr(nrp, i))/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                    !gp (pr%tepr(nrp, i, 1) + 273.15_r64))) * pr%cpr(nrp, 7, i)
                    !gp write(8, '(32f13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, 1), pr%tepr(nrp, i, 2), &
                    !gp (pr%cpr(nrp, j, i), j=1, nv-2), pr%phpr(nrp, i), &
                    !gp pr%cpr(nrp, nv, i)* rates%mga / rates%mgd * 1000, & !scc 08/09/2004
                    !gp tss, tp, tn , dosat, nh3, &
                    !gp pr%ninbpr(nrp, i)* rates%mgd / rates%mga / 1000, & !scc 08/09/2004
                    !gp pr%nipbpr(nrp, i)* rates%mgd / rates%mga / 1000, & !scc 08/09/2004
                    !gp pr%phitsavepr(nrp, i), pr%philsavepr(nrp, i), & !gp 20-oct-04
                    !gp pr%phinsavepr(nrp, i), pr%phipsavepr(nrp, i), & !gp 20-oct-04
                    !gp pr%phicsavepr(nrp, i), pr%phitotalsavepr(nrp, i) !gp 20-oct-04
                    !gp ! tdy=tdy+ system%dt * 24
                    j = 1 !gp water column is layer 1
                    tss = pr%cpr(nrp, 11, i, j) * rates%ada + pr%cpr(nrp, 2, i, j) + &
                        pr%cpr(nrp, 12, i, j)
                    tp = pr%cpr(nrp, 11, i, j) * rates%apa + pr%cpr(nrp, 9, i, j) + &
                        pr%cpr(nrp, 10, i, j)
                    tn = pr%cpr(nrp, 11, i, j) * rates%ana + pr%cpr(nrp, 6, i, j) + &
                        pr%cpr(nrp, 7, i, j) + pr%cpr(nrp, 8, i, j)
                    dosat = oxygen_saturation(pr%tepr(nrp, i, j), hydrau%reach(nrp)%elev)
                    nh3 = 1.0_r64/(1 + 10.0_r64 ** (-pr%phpr(nrp, i, j))/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                        (pr%tepr(nrp, i, j) + 273.15_r64))) * pr%cpr(nrp, 7, i, j)

                    !gp 01-nov-04 write(8, '(32f13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, 1), pr%tepr(nrp, i, 2), &
                    !gp 05-jul-05 write(8, '(33f13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, j), &
                    write(8, '(32F13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, j), & !gp 05-jul-05
                        (pr%cpr(nrp, k, i, j), k=1, nv-2), pr%phpr(nrp, i, j), &
                        pr%cpr(nrp, nv, i, j)* rates%mga / rates%mgd * 1000, & !scc 08/09/2004
                        tss, tp, tn , dosat, nh3, &
                        pr%ninbpr(nrp, i)* rates%mgd / rates%mga / 1000, & !scc 08/09/2004
                        pr%nipbpr(nrp, i)* rates%mgd / rates%mga / 1000, & !scc 08/09/2004
                        pr%phitsavepr(nrp, i), pr%philsavepr(nrp, i), & !gp 20-oct-04
                        pr%phinsavepr(nrp, i), pr%phipsavepr(nrp, i), & !gp 20-oct-04
                        pr%phicsavepr(nrp, i), pr%phitotalsavepr(nrp, i) !gp 20-oct-04

                end do
            end do

            !gp 27-Oct-04 (all code below here is new)
            !
            ! ---------------------------------------------------
            ! --- output of hyporheic pore water constituents ---
            ! ---------------------------------------------------
            !

            !gp diel output result for hyporheic pore water
            write(8,*)
            write(8,*) '** Diel hyporheic pore water and sediment flux **'

            write(8,'(41A13)') 't', ' Temps', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', &
                'NH4', 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
                'TSS', 'TP', 'TN', 'DOsat', 'NH3', &
                'DiagFluxDO', 'DiagFluxCBOD', 'DiagFluxNH4', 'DiagFluxNO3', 'DiagFluxSRP', 'DiagFluxIC', &
                'HypoFluxDO', 'HypoFluxCBOD', 'HypoFluxNH4', 'HypoFluxNO3', 'HypoFluxSRP', 'DiagFluxIC', &
                'TotFluxDO', 'TotFluxCBOD', 'TotFluxNH4', 'TotFluxNO3', 'TotFluxSRP', 'TotFluxIC'
            write(8,'(41A13)') 'hr', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', &
                'ugN/L', 'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', '', '', '', '', 'mgD/L', &
                'ugP/L', 'ugN/L', 'mg/L', 'ugN/L', &
                'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d', &
                'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d', &
                'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d'
            write(8,*) pr%nj
            do nrp=0, nr
                do i=0, pr%nj
                    select case (system%simhyporheicwq)
                      case ('Level 1', 'Level 2')
                        j = 2 !gp hyporheic pore water is layer 2
                        tss = pr%cpr(nrp, 11, i, j) * rates%ada + pr%cpr(nrp, 2, i, j) + &
                            pr%cpr(nrp, 12, i, j)
                        tp = pr%cpr(nrp, 11, i, j) * rates%apa + pr%cpr(nrp, 9, i, j) + &
                            pr%cpr(nrp, 10, i, j)
                        tn = pr%cpr(nrp, 11, i, j) * rates%ana + pr%cpr(nrp, 6, i, j) + &
                            pr%cpr(nrp, 7, i, j) + pr%cpr(nrp, 8, i, j)
                        nh3 = 1.0_r64/(1 + 10.0_r64 ** (-pr%phpr(nrp, i, j))/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                            (pr%tepr(nrp, i, j) + 273.15_r64))) * pr%cpr(nrp, 7, i, j)
                        dosat = oxygen_saturation(pr%tepr(nrp, i, j), hydrau%reach(nrp)%elev)
                        write(8, '(41F13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, j), &
                            (pr%cpr(nrp, k, i, j), k=1, nv-2), pr%phpr(nrp, i, j), &
                            tss, tp, tn , dosat, nh3, &
                            pr%diagfluxdopr(nrp, i), pr%diagfluxcbodpr(nrp, i), pr%diagfluxnh4pr(nrp, i), &
                            pr%diagfluxno3pr(nrp, i), pr%diagfluxsrppr(nrp, i), pr%diagfluxicpr(nrp, i), &
                            pr%hypofluxdopr(nrp, i), pr%hypofluxcbodpr(nrp, i), pr%hypofluxnh4pr(nrp, i), &
                            pr%hypofluxno3pr(nrp, i), pr%hypofluxsrppr(nrp, i), pr%hypofluxicpr(nrp, i), &
                            pr%diagfluxdopr(nrp, i) + pr%hypofluxdopr(nrp, i), &
                            pr%diagfluxcbodpr(nrp, i) + pr%hypofluxcbodpr(nrp, i), &
                            pr%diagfluxnh4pr(nrp, i) + pr%hypofluxnh4pr(nrp, i), &
                            pr%diagfluxno3pr(nrp, i) + pr%hypofluxno3pr(nrp, i), &
                            pr%diagfluxsrppr(nrp, i) + pr%hypofluxsrppr(nrp, i), &
                            pr%diagfluxicpr(nrp, i) + pr%hypofluxicpr(nrp, i)
                      case default !gp only write sediment temperatures and diagenesis fluxes if hyporheic wq is not being simulated
                        write(8, '(41F13.4)') pr%tdy(i)*24, pr%tepr(nrp, i, 2), (0.0*k,k=3,23), &
                            pr%diagfluxdopr(nrp, i), pr%diagfluxcbodpr(nrp, i), pr%diagfluxnh4pr(nrp, i), &
                            pr%diagfluxno3pr(nrp, i), pr%diagfluxsrppr(nrp, i), pr%diagfluxicpr(nrp, i), &
                            0.0,0.0,0.0,0.0,0.0,0.0, &
                            pr%diagfluxdopr(nrp, i), pr%diagfluxcbodpr(nrp, i), pr%diagfluxnh4pr(nrp, i), &
                            pr%diagfluxno3pr(nrp, i), pr%diagfluxsrppr(nrp, i), pr%diagfluxicpr(nrp, i)
                    end select
                end do
            end do


!gp 05-jul-05 diel fluxes for heat/do/co2
            write(8,*)
            write(8,*) '** Diel fluxes of heat (W/m^2), DO (gO2/m^2/d), and CO2 (gC/m^2/d) **'

            write(8,'(34A13)') 't', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', &
                'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', &
                'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2'

            write(8,'(34A13)') 'hr', 'solar', 'longat', 'back', 'air conv', 'evap', 'sed cond', 'hyporheic', 'tribs/GW', &
                'advec/disp', 'reaeration', 'fast CBOD', 'slow CBOD', 'COD', 'nitrif', &
                'phyto resp', 'phyto photo', 'botalg resp', 'botalg photo', &
                'SOD', 'hyporheic', 'tribs/GW', 'advec/disp', &
                'reaeration', 'fast CBOD', 'slow CBOD', &
                'phyto resp', 'phyto photo', 'botalg resp', 'botalg photo', &
                'SOD', 'hyporheic', 'tribs/GW', 'advec/disp'

            write(8,'(I13)') pr%nj
            do nrp=0, nr
                do i=0, pr%nj
                    write(8, '(34F13.4)') pr%tdy(i)*24, &
                        pr%pr_saveheatfluxjsnt(nrp, i), & !'heat fluxes
                        pr%pr_saveheatfluxlongat(nrp, i), &
                        pr%pr_saveheatfluxback(nrp, i), &
                        pr%pr_saveheatfluxconv(nrp, i), &
                        pr%pr_saveheatfluxevap(nrp, i), &
                        pr%pr_saveheatfluxjsed(nrp, i), &
                        pr%pr_saveheatfluxjhyporheic(nrp, i), &
                        pr%pr_saveheatfluxtribs(nrp, i), &
                        pr%pr_saveheatfluxadvecdisp(nrp, i), &
                        pr%pr_savedofluxreaer(nrp, i), & !'do fluxes
                        pr%pr_savedofluxcbodfast(nrp, i), &
                        pr%pr_savedofluxcbodslow(nrp, i), &
                        pr%pr_savedofluxcod(nrp, i), &
                        pr%pr_savedofluxnitrif(nrp, i), &
                        pr%pr_savedofluxphytoresp(nrp, i), &
                        pr%pr_savedofluxphytophoto(nrp, i), &
                        pr%pr_savedofluxbotalgresp(nrp, i), &
                        pr%pr_savedofluxbotalgphoto(nrp, i), &
                        pr%pr_savedofluxsod(nrp, i), &
                        pr%pr_savedofluxhyporheic(nrp, i), &
                        pr%pr_savedofluxtribs(nrp, i), &
                        pr%pr_savedofluxadvecdisp(nrp, i), &
                        pr%pr_saveco2fluxreaer(nrp, i), & !'co2 fluxes
                        pr%pr_saveco2fluxcbodfast(nrp, i), &
                        pr%pr_saveco2fluxcbodslow(nrp, i), &
                        pr%pr_saveco2fluxphytoresp(nrp, i), &
                        pr%pr_saveco2fluxphytophoto(nrp, i), &
                        pr%pr_saveco2fluxbotalgresp(nrp, i), &
                        pr%pr_saveco2fluxbotalgphoto(nrp, i), &
                        pr%pr_saveco2fluxsod(nrp, i), &
                        pr%pr_saveco2fluxhyporheic(nrp, i), &
                        pr%pr_saveco2fluxtribs(nrp, i), &
                        pr%pr_saveco2fluxadvecdisp(nrp, i)

                end do
            end do


!gp 17-feb-05
!only output diel if showdielresults = "yes"
!
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

        end if

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        !hyporheic pore water
        call hyporheic_pore_water(pr, topo, system, rates, nr)


    end subroutine



    Subroutine hydraulics_summary(topo, hydrau, rates, nr)

        type(riverhydraulics_type), intent(in) :: hydrau
        type(t_rivertopo), intent(in) :: topo
        type(rates_t), intent(in) :: rates
        integer(i32), intent(in) :: nr

        character(len=30) reaformular
        integer(i32) i

        !calculate ka(0) for output, not used in calculations
        !kawind = 0.728_r64 * uw(0) ** 0.5_r64 - 0.317_r64 * uw(0) + 0.0372_r64 * uw(0) ** 2
        !ka(0) = kau(0) + kawind / depth(0)
        write(8,*) '** Hydraulics Summary **'
        write (8,'(A9, 9A12, A25)')'Downstream', 'Hydraulics', "E'", 'H', 'B', 'Ac', 'U', 'trav time', &
            'slope', 'Reaeration', 'Reaeration formulas'
        write (8,'(A9, 9A12, A25)') 'distance', 'Q,m3/s', 'm3/s', 'm', 'm', 'm2', 'mps', 'd', '', &
            'ka,20,/d', 'water/wind'
        do i = 0, nr
            if (i==0) then
                reaformular =''
            else
                if (rates%kawindmethod=='None') then
                    reaformular = trim(hydrau%reach(i)%kaf) // '/No wind'
                else
                    reaformular = trim(hydrau%reach(i)%kaf) // rates%kawindmethod
                end if
            end if
            write(8,'(F9.4, 9F12.5, 1A35)') topo%reach(i)%xrdn, hydrau%reach(i)%q, hydrau%reach(i)%epout, &
                hydrau%reach(i)%depth, hydrau%reach(i)%b, hydrau%reach(i)%ac, &
                hydrau%reach(i)%u, hydrau%reach(i)%trav, hydrau%reach(i)%s, &
                hydrau%reach(i)%ka, reaformular
        end do

    end subroutine hydraulics_summary







    subroutine hyporheic_pore_water(pr, topo, system, rates, nr)

        type(outdata_t), intent(in) :: pr
        type(t_rivertopo), intent(in) :: topo
        type(systemparams), intent(in) :: system
        type(rates_t), intent(in) :: rates
        integer(i32), intent(in) :: nr

        real(r64) cbodu, tkn, toc, tss
        integer(i32) i, j, k

        ! Output temperature for the hyporheic pore water
        write(8,*)
        write(8,*) '** Temperature summary (hyporheic pore water temperature) **'
        write(8,'(A5, A10, 4A10)') 'Reach', '', 'Distance','Temp(C)', 'Temp(C)', 'Temp(C)'
        write(8,'(A5, A10, 4A10)') 'Label', '', 'x(km)','Average', 'Minimum', 'Maximum'
        j = 2 !gp 27-oct-04 hypoprheic pore water is layer 2
        do i=0, nr
            write (8,'(A15, 4F10.4)') topo%reach(i)%rname, topo%reach(i)%xpm, pr%teav(i, j), &
                pr%temn(i, j), pr%temx(i, j) !gp add nl dimension
        end do

        select case (system%simhyporheicwq)
          case ('Level 1', 'Level 2')

            ! Output concentrations for the hyporheic pore water constituents
            write(8,*)
            write(8,*) '** Daily average water quality summary (hyporheic pore water constituents) **'
            write(8,'(A5, A10, 26A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
                'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
                'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
                'NH3', 'DO sat', 'pH sat'
            write(8,'(A5, A10, 26A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
                'ugN/L', 'ugN/L', &
                'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
                '', '', '', '', 'mgD/L', '', &
                '', '', ''
            do i = 0, nr
                j = 2 !gp hyporheic pore water is layer 2
                toc = (pr%cav(i, 4, j) + pr%cav(i, 5, j)) / rates%roc + &
                    rates%aca * pr%cav(i, 11, j) + rates%aca / rates%ada * pr%cav(i, 12, j)
                tkn = pr%cav(i, 6, j) + pr%cav(i, 7, j) + rates%ana * pr%cav(i, 11, j)
                tss = rates%ada * pr%cav(i, 11, j) + pr%cav(i, 2, j) + pr%cav(i, 12, j)
                cbodu = pr%cav(i, 4, j) + pr%cav(i, 5, j) + rates%roa * pr%cav(i, 11, j) + &
                    rates%roc * rates%aca / rates%ada * pr%cav(i, 12, j)
                write(8,'(A15, 26F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%phav(i, j) , &
                    toc, pr%tnav(i, j), pr%tpav(i, j), &
                    tkn, tss, cbodu, pr%nh3av(i, j), pr%osav(i, j), &
                    pr%phsav(i, j)
            end do

            ! output min concentrations
            write(8,*)
            write(8,*) '** Daily minimum water quality summary (hyporheic pore water constituents) **'
            write(8,'(A5, A10, 24A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
                'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
                'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
                'NH3'
            write(8,'(A5, A10, 24A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
                'ugN/L', 'ugN/L', &
                'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', 'cfu/100mL', 'user defined', 'mgCaCO3/L', 's.u.', &
                '', '', '', '', 'mgD/L', '', ''
            do i = 0, nr
                j = 2 !gp hyporheic pore water is layer 2
                toc = (pr%cmn(i, 4, j) + pr%cmn(i, 5, j)) / rates%roc + &
                    rates%aca * pr%cmn(i, 11, j) + rates%aca / rates%ada * pr%cmn(i, 12, j)
                tkn = pr%cmn(i, 6, j) + pr%cmn(i, 7, j) + rates%ana * pr%cmn(i, 11, j)
                tss = rates%ada * pr%cmn(i, 11, j) + pr%cmn(i, 2, j) + pr%cmn(i, 12, j)
                cbodu = pr%cmn(i, 4, j) + pr%cmn(i, 5, j) + rates%roa * pr%cmn(i, 11, j) + &
                    rates%roc * rates%aca / rates%ada * pr%cmn(i, 12, j)
                write(8,'(A15, 24F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%phmn(i, j) , &
                    toc, pr%tnmn(i, j), pr%tpmn(i, j), &
                    tkn, tss, cbodu, pr%nh3mn(i, j)
            end do

            ! Output max concentrations
            write(8,*)
            write(8,*) '** Daily maximum water quality summary (hyporheic pore water constituents) **'
            write(8,'(A5, A10, 24A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
                'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
                'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
                'NH3'
            write(8,'(A5, A10, 24A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
                'ugN/L', 'ugN/L', &
                'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', 'cfu/100mL', 'user defined', 'mgCaCO3/L', 's.u.', &
                '', '', '', '', 'mgD/L', '', ''
            do i = 0, nr
                j = 2 !gp hyporheic pore water is layer 2
                toc = (pr%cmx(i, 4, j) + pr%cmx(i, 5, j)) / rates%roc + &
                    rates%aca * pr%cmx(i, 11, j) + rates%aca / rates%ada * pr%cmx(i, 12, j)
                tkn = pr%cmx(i, 6, j) + pr%cmx(i, 7, j) + rates%ana * pr%cmx(i, 11, j)
                tss = rates%ada * pr%cmx(i, 11, j) + pr%cmx(i, 2, j) + pr%cmx(i, 12, j)
                cbodu = pr%cmx(i, 4, j) + pr%cmx(i, 5, j) + rates%roa * pr%cmx(i, 11, j) + &
                    rates%roc * rates%aca / rates%ada * pr%cmx(i, 12, j)
                write(8,'(A15, 24F13.6)')topo%reach(i)%rname, topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%phmx(i, j) , &
                    toc, pr%tnmx(i, j), pr%tpmx(i, j), &
                    tkn, tss, cbodu, pr%nh3mx(i, j)
            end do

        end select

    end subroutine hyporheic_pore_water


end module
