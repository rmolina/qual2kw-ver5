module m_integration_data

    use m_constants

    !gp 01-nov-04
    implicit none
    private

    public :: integration_data_t, &
        saveheatfluxtribs, saveheatfluxadvecdisp, saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, &
        saveheatfluxconv, saveheatfluxevap, saveheatfluxjsed, saveheatfluxjhyporheic, &
        os, phitsave, savedofluxheadwater, saveco2fluxheadwater, savedofluxtribs, savedofluxadvecdisp, &
        saveco2fluxtribs, saveco2fluxadvecdisp, phinsave, phipsave, phicsave, philsave, phitotalsave , &
        savebotalgphoto, savebotalgresp, savebotalgdeath, savebotalgnetgrowth, &
        sodpr, jnh4pr, jno3pr, jch4pr, jsrppr, csodpr, &
        diagfluxdo, diagfluxcbod, diagfluxnh4, diagfluxno3, diagfluxsrp, diagfluxic, savedofluxreaer, &
        savedofluxcbodfast, savedofluxcbodslow, savedofluxnitrif, savedofluxphytoresp, &
        savedofluxphytophoto, savedofluxbotalgresp, savedofluxbotalgphoto, savedofluxsod, savedofluxcod, &
        saveco2fluxreaer, saveco2fluxcbodfast, saveco2fluxcbodslow, saveco2fluxphytoresp, saveco2fluxphytophoto, &
        saveco2fluxbotalgresp, saveco2fluxbotalgphoto, saveco2fluxsod, hypofluxdo, hypofluxcbod, hypofluxnh4, &
        hypofluxno3, hypofluxsrp, hypofluxic

    type integration_data_t
! integer(i4b) np, nc !days, step in each day
        !dependent variables
        real(dp), dimension(:), pointer :: inb, ipb
        !gp 27-oct-04 real(dp), dimension(:,:), pointer :: te, c
        real(dp), dimension(:,:), pointer :: te !gp
        real(dp), dimension(:,:,:), pointer :: c !gp add dimension for nl
    end type integration_data_t

    interface integration_data_t
        procedure :: integration_data_ctor
    end interface integration_data_t


    !sediment nutrient flux
    real(dp), dimension(:), pointer :: sodpr, jnh4pr, jno3pr, jch4pr, jsrppr

    !gp real(dp), dimension(:), pointer :: os
    real(dp), dimension(:,:), pointer :: os !gp

    !gp 20-oct-04 bottom algae limitation factors (fraction of maximum potential growth)
    real(dp), dimension(:), pointer :: phitotalsave, phitsave, philsave, phinsave, phipsave, phicsave !gp 20-oct-04

    !gp 28-oct-04 diagenesis and hyporheic sediment/water fluxes (positive is source to water from sediment)
    real(dp), dimension(:), pointer :: diagfluxdo, diagfluxcbod, diagfluxnh4, diagfluxno3, diagfluxsrp, diagfluxic !gp 28-oct-04
    real(dp), dimension(:), pointer :: hypofluxdo, hypofluxcbod, hypofluxnh4, hypofluxno3, hypofluxsrp, hypofluxic !gp 28-oct-04

    !gp 29-oct-04
    real(dp), dimension(:), pointer :: csodpr

    !gp 05-jul-05 heat/do/co2 fluxes
    real(dp), dimension(:), pointer :: saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, saveheatfluxconv
    real(dp), dimension(:), pointer :: saveheatfluxevap, saveheatfluxjsed, saveheatfluxjhyporheic, saveheatfluxtribs
    real(dp), dimension(:), pointer :: saveheatfluxadvecdisp
    real(dp), dimension(:), pointer :: savedofluxreaer, savedofluxcbodfast, savedofluxcbodslow, savedofluxnitrif
    real(dp), dimension(:), pointer :: savedofluxphytoresp, savedofluxphytophoto, savedofluxbotalgresp, savedofluxbotalgphoto
    real(dp), dimension(:), pointer :: savedofluxsod, savedofluxcod, savedofluxhyporheic, savedofluxtribs
    real(dp), dimension(:), pointer :: savedofluxadvecdisp, savedofluxheadwater
    real(dp), dimension(:), pointer :: saveco2fluxreaer, saveco2fluxcbodfast, saveco2fluxcbodslow
    real(dp), dimension(:), pointer :: saveco2fluxphytoresp, saveco2fluxphytophoto, saveco2fluxbotalgresp, saveco2fluxbotalgphoto
    real(dp), dimension(:), pointer :: saveco2fluxsod, saveco2fluxhyporheic, saveco2fluxtribs
    real(dp), dimension(:), pointer :: saveco2fluxadvecdisp, saveco2fluxheadwater

    !gp 25-jun-09
    real(dp), dimension(:), pointer :: savebotalgphoto, savebotalgresp, savebotalgdeath, savebotalgnetgrowth
    !real(dp), dimension(:), pointer :: av_botalgphoto, av_botalgresp, av_botalgdeath, av_botalgnetgrowth

contains

    !user-defined type 'integral_type' constructor
    function integration_data_ctor(nr) result(integral)

! use class_rivertopo
        type(integration_data_t) integral
        integer(i4b), intent(in) :: nr

        !gp 05-jul-05 integer(i4b) i, status(29) !gp 28-oct-04
        !gp 25-jun-09 integer(i4b) i, status(64) !gp 05-jul-05
        integer(i4b) i, status(68) !gp 25-jun-09

        !gp allocate(os(0:nr), stat=status(1))
        allocate(os(0:nr, nl), stat=status(1)) !gp

        allocate(integral%inb(0:nr), stat=status(2))
        allocate(integral%ipb(0:nr), stat=status(3))

        !gp allocate(integral%c(0:nr+1, nv), stat=status(4))
        allocate(integral%c(0:nr+1, nv, nl), stat=status(4))

        allocate(integral%te(0:nr+1, nl), stat=status(5))
        allocate(jnh4pr(0:nr), stat=status(6))
        allocate(jno3pr(0:nr), stat=status(7))
        allocate(jch4pr(0:nr), stat=status(8))
        allocate(jsrppr(0:nr), stat=status(9))
        allocate(sodpr(0:nr), stat=status(10))

        !gp 20-oct-04 growth limitation factors for bottom algae
        allocate(phitotalsave(0:nr), stat=status(11))
        allocate(phitsave(0:nr), stat=status(12))
        allocate(philsave(0:nr), stat=status(13))
        allocate(phinsave(0:nr), stat=status(14))
        allocate(phipsave(0:nr), stat=status(15))
        allocate(phicsave(0:nr), stat=status(16)) !gp end new block

        !gp 28-oct-04 diagenesis fluxes between sediment/water
        allocate(diagfluxdo(0:nr), stat=status(17))
        allocate(diagfluxcbod(0:nr), stat=status(18))
        allocate(diagfluxnh4(0:nr), stat=status(19))
        allocate(diagfluxno3(0:nr), stat=status(20))
        allocate(diagfluxsrp(0:nr), stat=status(21))
        allocate(diagfluxic(0:nr), stat=status(22)) !gp end new block

        !gp 28-oct-04 hyporheic exchange fluxes between sediment/water
        allocate(hypofluxdo(0:nr), stat=status(23))
        allocate(hypofluxcbod(0:nr), stat=status(24))
        allocate(hypofluxnh4(0:nr), stat=status(25))
        allocate(hypofluxno3(0:nr), stat=status(26))
        allocate(hypofluxsrp(0:nr), stat=status(27))
        allocate(hypofluxic(0:nr), stat=status(28)) !gp end new block

        !gp 29-oct-04
        allocate(csodpr(0:nr), stat=status(29))

        !gp 05-jul-05 heat/do/co2 fluxes
        allocate(saveheatfluxjsnt(0:nr), stat=status(30))
        allocate(saveheatfluxlongat(0:nr), stat=status(31))
        allocate(saveheatfluxback(0:nr), stat=status(32))
        allocate(saveheatfluxconv(0:nr), stat=status(33))
        allocate(saveheatfluxevap(0:nr), stat=status(34))
        allocate(saveheatfluxjsed(0:nr), stat=status(35))
        allocate(saveheatfluxjhyporheic(0:nr), stat=status(36))
        allocate(saveheatfluxtribs(0:nr), stat=status(37))
        allocate(saveheatfluxadvecdisp(0:nr), stat=status(38))
        allocate(savedofluxreaer(0:nr), stat=status(39))
        allocate(savedofluxcbodfast(0:nr), stat=status(40))
        allocate(savedofluxcbodslow(0:nr), stat=status(41))
        allocate(savedofluxcod(0:nr), stat=status(42))
        allocate(savedofluxnitrif(0:nr), stat=status(43))
        allocate(savedofluxphytoresp(0:nr), stat=status(44))
        allocate(savedofluxphytophoto(0:nr), stat=status(45))
        allocate(savedofluxbotalgresp(0:nr), stat=status(46))
        allocate(savedofluxbotalgphoto(0:nr), stat=status(47))
        allocate(savedofluxsod(0:nr), stat=status(48))
        allocate(savedofluxhyporheic(0:nr), stat=status(49))
        allocate(savedofluxtribs(0:nr), stat=status(50))
        allocate(savedofluxadvecdisp(0:nr), stat=status(51))
        allocate(savedofluxheadwater(0:nr), stat=status(52))
        allocate(saveco2fluxreaer(0:nr), stat=status(53))
        allocate(saveco2fluxcbodfast(0:nr), stat=status(54))
        allocate(saveco2fluxcbodslow(0:nr), stat=status(55))
        allocate(saveco2fluxphytoresp(0:nr), stat=status(56))
        allocate(saveco2fluxphytophoto(0:nr), stat=status(57))
        allocate(saveco2fluxbotalgresp(0:nr), stat=status(58))
        allocate(saveco2fluxbotalgphoto(0:nr), stat=status(59))
        allocate(saveco2fluxsod(0:nr), stat=status(60))
        allocate(saveco2fluxhyporheic(0:nr), stat=status(61))
        allocate(saveco2fluxtribs(0:nr), stat=status(62))
        allocate(saveco2fluxadvecdisp(0:nr), stat=status(63))
        allocate(saveco2fluxheadwater(0:nr), stat=status(64))

        !gp 25-jun-09
        allocate(savebotalgphoto(0:nr), stat=status(65))
        allocate(savebotalgresp(0:nr), stat=status(66))
        allocate(savebotalgdeath(0:nr), stat=status(67))
        allocate(savebotalgnetgrowth(0:nr), stat=status(68))
        !allocate(av_botalgphoto(0:nr), stat=status(69))
        !allocate(av_botalgresp(0:nr), stat=status(70))
        !allocate(av_botalgdeath(0:nr), stat=status(71))
        !allocate(av_botalgnetgrowth(0:nr), stat=status(72))

        !initialize
        os=0;
        integral%inb=0; integral%ipb=0; integral%c=0; integral%te=0
        sodpr=0; jnh4pr=0; jno3pr=0; jch4pr=0; jsrppr=0

        !gp 20-oct-04
        phitotalsave=0; phitsave=0; philsave=0; phinsave=0; phipsave=0; phicsave=0 !gp end new block

        !gp 25-jun-09
        savebotalgphoto=0; savebotalgresp=0; savebotalgdeath=0; savebotalgnetgrowth=0
        !av_botalgphoto=0; av_botalgresp=0; av_botalgdeath=0; av_botalgnetgrowth=0

        !28-oct-04
        diagfluxdo=0; diagfluxcbod=0; diagfluxnh4=0; diagfluxno3=0; diagfluxsrp=0; diagfluxic=0
        hypofluxdo=0; hypofluxcbod=0; hypofluxnh4=0; hypofluxno3=0; hypofluxsrp=0; hypofluxic=0 !gp end new block

        !gp 29-oct-04
        csodpr=0

        !gp 05-jul-05
        saveheatfluxjsnt=0; saveheatfluxlongat=0; saveheatfluxback=0; saveheatfluxconv=0
        saveheatfluxevap=0; saveheatfluxjsed=0; saveheatfluxjhyporheic=0; saveheatfluxtribs=0
        saveheatfluxadvecdisp=0
        savedofluxreaer=0; savedofluxcbodfast=0; savedofluxcbodslow=0; savedofluxnitrif=0
        savedofluxphytoresp=0; savedofluxphytophoto=0; savedofluxbotalgresp=0; savedofluxbotalgphoto=0
        savedofluxsod=0; savedofluxcod=0; savedofluxhyporheic=0; savedofluxtribs=0
        savedofluxadvecdisp=0; savedofluxheadwater=0
        saveco2fluxreaer=0; saveco2fluxcbodfast=0; saveco2fluxcbodslow=0
        saveco2fluxphytoresp=0; saveco2fluxphytophoto=0; saveco2fluxbotalgresp=0; saveco2fluxbotalgphoto=0
        saveco2fluxsod=0; saveco2fluxhyporheic=0; saveco2fluxtribs=0
        saveco2fluxadvecdisp=0; saveco2fluxheadwater=0

        !gp 20-oct-04 do i=2, 10
        !gp 05-jul-05 do i=2, 29
        do i=2, 64
            if (status(i)==1) then
                write(8,*) '** Class_Integration:Integration_ failed. Insufficient Memory **'
                stop !class_integration:integration_ failed. insufficient memory!'
            end if
        end do

    end function integration_data_ctor

end module m_integration_data
