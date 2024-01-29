!classintegration

module class_integration
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_derivs, only: derivs
    use class_hydraulics, only: riverhydraulics_type
    use class_integrationdata, only: integral_type, integration_, &
        saveheatfluxtribs, saveheatfluxadvecdisp, saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, &
        saveheatfluxconv, saveheatfluxevap, saveheatfluxjsed, saveheatfluxjhyporheic, &
        os, phitsave, savedofluxtribs, savedofluxadvecdisp, &
        saveco2fluxtribs, saveco2fluxadvecdisp, phinsave, phipsave, phicsave, philsave, phitotalsave , &
        savebotalgphoto, savebotalgresp, savebotalgdeath, savebotalgnetgrowth, &
        sodpr, &
        diagfluxdo, diagfluxcbod, diagfluxnh4, diagfluxno3, diagfluxsrp, diagfluxic, savedofluxreaer, &
        savedofluxcbodfast, savedofluxcbodslow, savedofluxnitrif,  savedofluxphytoresp, &
        savedofluxphytophoto, savedofluxbotalgresp, savedofluxbotalgphoto, savedofluxsod, savedofluxcod, &
        saveco2fluxreaer, saveco2fluxcbodfast, saveco2fluxcbodslow, saveco2fluxphytoresp, saveco2fluxphytophoto, &
        saveco2fluxbotalgresp, saveco2fluxbotalgphoto, saveco2fluxsod, hypofluxdo, hypofluxcbod, hypofluxnh4, &
        hypofluxno3, hypofluxsrp, hypofluxic
    use class_lightheat, only: lightheat, heatbalance, lightextinction
    use class_phsolve, only: chemrates, modfp2, phsolnewton, phsolbisect, phsolbrent, ct
    use class_solarcalc, only: solar_type, solarcalc
    use class_sourcein, only: load, sourcescalc
    use class_systemparams, only: systemparams
    use m_downstream_boundary, only: downstream_boundary_t, instanteousdownstreamboundary
    use m_meteorology, only: meteorology_t, instanteousmeteo
    use m_output, only: outdata_t, output
    use m_oxygen, only: oxygen_inhibition_and_enhancement, oxygen_saturation
    use m_rates, only: rates_t
    use m_tempadjust, only: temp_adjust
    use m_upstream_boundary, only: upstream_boundary_t, instanteousheadwater
    use nrtype, only: nv, nl, adam, bdam, cpw, rhow

    implicit none

contains
!----------------------------------------------------------------------------------
    subroutine integration(sys, rates, meteo, solar, hw, db, hydrau, pr, nr)


        type(integral_type) intg !integral data structure
        type(systemparams) sys
        type(rates_t), intent(in) :: rates
        type(meteorology_t) meteo
        type(solar_type) solar !solar radiation
        type(upstream_boundary_t) hw !headwater
        type(downstream_boundary_t) db !downstream boundary
        type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
        type(outdata_t), intent(out) :: pr !print out data structure
        integer(i32), intent(in) :: nr !number of reach
        real(r64) t
        !heat related values

        integer(i32) ip, ic, i, j, k
        real(r64) ph
        !gp 27-oct-04 real(r64) dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64) dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp add dim to dc for nl

        !gp 24-jun-09
        real(r64) tss, tp, tn, dosat, nh3
        integer(i32) iskip, nskip
        nskip = sys%nc / 32 !for output of dynamic calculations every 45 minutes (32 times per day)
        iskip = 0 !initialize counter for skipping dynamic output

        !gp 24-jun-09
        !character(len=30), intent(in) :: imethph ! ph solve method

        intg=integration_(nr)

        !set initial conditions
        dc=0; dinb=0; dipb=0; dte=0

        !initial time
        t=0

        !initial temperature and concentrations
        !call instanteousheadwater(hdboundary, sys%tday, 1, intg%te(0,1), intg%c(0,:), ph)
        !gp call instanteousheadwater(hw, t, intg%te(0,1), intg%c(0,:), ph)
        call instanteousheadwater(hw, t, intg%te(0,1), intg%c(0,:,1), ph) !gp

        do i = 1, nr

            !gp 21-nov-06
            !do j = 1, nl
            ! intg%te(i, j) = intg%te(0, 1)
            !!gp 27-oct-04 end do
            ! do k = 1, nv - 1
            ! intg%c(i, k, j) = intg%c(0, k, 1)
            ! end do
            !end do !gp
            !!gp 11-jan-05 if (rates%typef == "first-order") intg%c(i, nv) = 1.0
            !intg%c(i, nv, 1) = hydrau%reach(i)%botalg0 * rates%adc * rates%aca !'convert mga/m^2 to gd/m^2
            !if (rates%typef == "first-order" .and. hydrau%reach(i)%botalg0==0) intg%c(i, nv, 1) = 1.0
            do j = 1, nl
                ! if (hydrau%reach(i)%te_ini < 0) then
                !     intg%te(i, j) = intg%te(0, 1)
                ! else
                !     intg%te(i, j) = hydrau%reach(i)%te_ini
                ! end if
                intg%te(i, j) = ini_or_prev(hydrau%reach(i)%te_ini, intg%te(0, 1))

                ! if (hydrau%reach(i)%c01_ini < 0) then
                !     intg%c(i, 1, j) = intg%c(0, 1, 1)
                ! else
                !     intg%c(i, 1, j) = hydrau%reach(i)%c01_ini
                ! end if
                intg%c(i, 1, j) = ini_or_prev(hydrau%reach(i)%c01_ini, intg%c(0, 1, 1))

                ! if (hydrau%reach(i)%c02_ini < 0) then
                !     intg%c(i, 2, j) = intg%c(0, 2, 1)
                ! else
                !     intg%c(i, 2, j) = hydrau%reach(i)%c02_ini
                ! end if
                intg%c(i, 2, j) = ini_or_prev(hydrau%reach(i)%c02_ini, intg%c(0, 2, 1))

                ! if (hydrau%reach(i)%c03_ini < 0) then
                !     intg%c(i, 3, j) = intg%c(0, 3, 1)
                ! else
                !     intg%c(i, 3, j) = hydrau%reach(i)%c03_ini
                ! end if
                intg%c(i, 3, j) = ini_or_prev(hydrau%reach(i)%c03_ini, intg%c(0, 3, 1))

                ! if (hydrau%reach(i)%c04_ini < 0) then
                !     intg%c(i, 4, j) = intg%c(0, 4, 1)
                ! else
                !     intg%c(i, 4, j) = hydrau%reach(i)%c04_ini
                ! end if
                intg%c(i, 4, j) = ini_or_prev(hydrau%reach(i)%c04_ini, intg%c(0, 4, 1))

                ! if (hydrau%reach(i)%c05_ini < 0) then
                !     intg%c(i, 5, j) = intg%c(0, 5, 1)
                ! else
                !     intg%c(i, 5, j) = hydrau%reach(i)%c05_ini
                ! end if
                intg%c(i, 5, j) = ini_or_prev(hydrau%reach(i)%c05_ini, intg%c(0, 5, 1))

                ! if (hydrau%reach(i)%c06_ini < 0) then
                !     intg%c(i, 6, j) = intg%c(0, 6, 1)
                ! else
                !     intg%c(i, 6, j) = hydrau%reach(i)%c06_ini
                ! end if
                intg%c(i, 6, j) = ini_or_prev(hydrau%reach(i)%c06_ini, intg%c(0, 6, 1))

                ! if (hydrau%reach(i)%c07_ini < 0) then
                !     intg%c(i, 7, j) = intg%c(0, 7, 1)
                ! else
                !     intg%c(i, 7, j) = hydrau%reach(i)%c07_ini
                ! end if
                intg%c(i, 7, j) = ini_or_prev(hydrau%reach(i)%c07_ini, intg%c(0, 7, 1))

                ! if (hydrau%reach(i)%c08_ini < 0) then
                !     intg%c(i, 8, j) = intg%c(0, 8, 1)
                ! else
                !     intg%c(i, 8, j) = hydrau%reach(i)%c08_ini
                ! end if
                intg%c(i, 8, j) = ini_or_prev(hydrau%reach(i)%c08_ini, intg%c(0, 8, 1))

                ! if (hydrau%reach(i)%c09_ini < 0) then
                !     intg%c(i, 9, j) = intg%c(0, 9, 1)
                ! else
                !     intg%c(i, 9, j) = hydrau%reach(i)%c09_ini
                ! end if
                intg%c(i, 9, j) = ini_or_prev(hydrau%reach(i)%c09_ini, intg%c(0, 9, 1))

                ! if (hydrau%reach(i)%c10_ini < 0) then
                !     intg%c(i, 10, j) = intg%c(0, 10, 1)
                ! else
                !     intg%c(i, 10, j) = hydrau%reach(i)%c10_ini
                ! end if
                intg%c(i, 10, j) = ini_or_prev(hydrau%reach(i)%c10_ini, intg%c(0, 10, 1))

                ! if (hydrau%reach(i)%c11_ini < 0) then
                !     intg%c(i, 11, j) = intg%c(0, 11, 1)
                ! else
                !     intg%c(i, 11, j) = hydrau%reach(i)%c11_ini
                ! end if
                intg%c(i, 11, j) = ini_or_prev(hydrau%reach(i)%c11_ini, intg%c(0, 11, 1))

                ! if (hydrau%reach(i)%c12_ini < 0) then
                !     intg%c(i, 12, j) = intg%c(0, 12, 1)
                ! else
                !     intg%c(i, 12, j) = hydrau%reach(i)%c12_ini
                ! end if
                intg%c(i, 12, j) = ini_or_prev(hydrau%reach(i)%c12_ini, intg%c(0, 12, 1))

                ! if (hydrau%reach(i)%c13_ini < 0) then
                !     intg%c(i, 13, j) = intg%c(0, 13, 1)
                ! else
                !     intg%c(i, 13, j) = hydrau%reach(i)%c13_ini
                ! end if
                intg%c(i, 13, j) = ini_or_prev(hydrau%reach(i)%c13_ini, intg%c(0, 13, 1))

                ! if (hydrau%reach(i)%c14_ini < 0) then
                !     intg%c(i, 14, j) = intg%c(0, 14, 1)
                ! else
                !     intg%c(i, 14, j) = hydrau%reach(i)%c14_ini
                ! end if
                intg%c(i, 14, j) = ini_or_prev(hydrau%reach(i)%c14_ini, intg%c(0, 14, 1))

                ! if (hydrau%reach(i)%c15_ini < 0) then
                !     intg%c(i, nv-2, j) = intg%c(0, nv-2, 1)
                ! else
                !     intg%c(i, nv-2, j) = hydrau%reach(i)%c15_ini
                ! end if
                intg%c(i, nv-2, j) = ini_or_prev(hydrau%reach(i)%c15_ini, intg%c(0, nv-2, 1))

                if (hydrau%reach(i)%ph_ini < 0) then
                    intg%c(i, nv-1, j) = intg%c(0, nv-1, 1)
                else
                    intg%c(i, nv-1, j) = ct(hydrau%reach(i)%ph_ini, intg%c(i, nv-2, j), intg%te(i, j), intg%c(i, 1, j))
                end if

            end do
            if (hydrau%reach(i)%c17_ini < 0) then !bottom algae gd/m^2
                intg%c(i, nv, 1) = 0
            else
                intg%c(i, nv, 1) = hydrau%reach(i)%c17_ini
            end if
            if (rates%typef == "First-order" .and. intg%c(i, nv, 1) < 1) intg%c(i, nv, 1) = 1.0

            !gp 02-nov-07
            !if (hydrau%reach(i)%ninb_ini >= 0) then !intracellular n mgn/gd
            !intg%inb(i) = hydrau%reach(i)%ninb_ini * intg%c(i, nv, j)
            if (hydrau%reach(i)%ninb_ini >= 0) then !intracellular n mgn/gd
                intg%inb(i) = hydrau%reach(i)%ninb_ini * intg%c(i, nv, 1)

            end if

            !gp 02-nov-07
            !if (hydrau%reach(i)%nipb_ini >= 0) then !intracellular p mgp/gd
            !intg%ipb(i) = hydrau%reach(i)%nipb_ini * intg%c(i, nv, j)
            if (hydrau%reach(i)%nipb_ini >= 0) then !intracellular p mgp/gd
                intg%ipb(i) = hydrau%reach(i)%nipb_ini * intg%c(i, nv, 1)

            end if

            if (sys%simhyporheicwq == "Level 2" .and. rates%typeh == "First-order") intg%c(i, nv, 2) = 1.0 !gp 15-Nov-04
        END DO

        call print_report_integration_method(sys%imeth)

        !integration
        do ip = 1, sys%np !np

            if (ip >= sys%np) then

                !gp 17-feb-05
                !call save_init_step(nr, intg, pr, sys%dt, sys%imethph, rates)
                call save_init_step(nr, intg, pr, sys%dt, sys%imethph, sys%showdielresults, rates)

            end if

            !euler method
            if (sys%imeth == "Euler") then
                write(*,*) 'day', ip
                do ic = 1, sys%nc !time step in each day
                    call derivs(nr, meteo, solar,hw, db, hydrau, sys, intg%te, intg%c, intg%inb, &
                        intg%ipb, rates, dte, dc, dinb, dipb, t)
                    !finite method: new= old + old * derive
                    do i = 1, nr
                        do j = 1, nl

                            !gp 08-jan-10
                            !intg%te(i, j) = intg%te(i, j) + dte(i, j) * sys%dt
                            if (sys%statevariables /= "All except temperature") then
                                intg%te(i, j) = intg%te(i, j) + dte(i, j) * sys%dt

                                !gp 12-jan-10
                            else
                                if (hydrau%reach(i)%te_ini < 0) then
                                    !use diel headwater as diel temp for each reach during intg
                                    !call instanteousheadwater(hw, t, intg%te(0,1), intg%c(0,:,1), ph)
                                    intg%te(i, j) = intg%te(0, 1)
                                else
                                    !use initial temp for each reach as const temp during intg
                                    intg%te(i, j) = hydrau%reach(i)%te_ini
                                end if

                            end if

                            !gp 27-oct-04 end do
                            ! write (8, *) i, intg%c(i,3), dc(i,3)
                            ! if (i==1) then
                            ! write(8, '(i3, 16f26.18)') i, (intg%c(i, k), k=1, nv)
                            ! write(8, '(i3, 16f26.18)') i, (dc(i, k), k=1, nv)
                            ! write(8, '(i3, 4f26.18)') i, intg%inb(i), dinb(i), intg%ipb(i), dipb(i)
                            ! end if
                            do k = 1, nv
                                intg%c(i, k, j) = intg%c(i, k, j) + dc(i, k, j) * sys%dt
                                if (intg%c(i, k, j) < 0) intg%c(i, k, j) = 1.0e-6_r64
                            end do
                        end do !gp 27-oct-04
                        intg%inb(i) = intg%inb(i) + dinb(i) * sys%dt
                        intg%ipb(i) = intg%ipb(i) + dipb(i) * sys%dt
                        if (intg%inb(i) < 0) intg%inb(i) = 1.0e-6_r64
                        if (intg%ipb(i) < 0) intg%ipb(i) = 1.0e-6_r64
                    end do


                    !gp 24-jun-09 write dynamic output every 45 minutes
                    if (sys%writedynamic == "Yes") then
                        if (iskip == 0) then
                            do i=0, nr
                                j = 1
                                tss = intg%c(i, 11, j) * rates%ada + intg%c(i, 2, j) + &
                                    intg%c(i, 12, j)
                                tp = intg%c(i, 11, j) * rates%apa + intg%c(i, 9, j) + &
                                    intg%c(i, 10, j)
                                tn = intg%c(i, 11, j) * rates%ana + intg%c(i, 6, j) + &
                                    intg%c(i, 7, j) + intg%c(i, 8, j)
                                dosat = oxygen_saturation(intg%te(i, j), hydrau%reach(i)%elev)

                                !gp 23-nov-09
                                !if (sys%imethph == "newton-raphson") then
                                !call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !else
                                !call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !end if
                                if (sys%imethph == "Newton-Raphson") then
                                    call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                elseif (sys%IMethpH == "Bisection") then
                                    call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                else
                                    call phsolbrent(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                end if

                                nh3 = 1.0_r64/(1 + 10.0_r64 ** (-ph)/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                                    (intg%te(i, j) + 273.15_r64))) * intg%c(i, 7, j)
                                write(12, '(I13, 41F13.4)') i, t, intg%te(i, j), &
                                    (intg%c(i, k, j), k=1, nv-2), ph, &
                                    intg%c(i, nv, j), intg%c(i, nv, j)* rates%mga / rates%mgd * 1000, &
                                    tss, tp, tn , dosat, nh3, &
                                    intg%inb(i), &
                                    intg%inb(i), &
                                    phitsave(i), philsave(i), &
                                    phinsave(i), phipsave(i), &
                                    phicsave(i), phitotalsave(i), &
                                    savebotalgphoto(i), savebotalgresp(i), savebotalgdeath(i), savebotalgnetgrowth(i), &
                                    savebotalgphoto(i)/rates%ada, savebotalgresp(i)/rates%ada, &
                                    savebotalgdeath(i)/rates%ada, savebotalgnetgrowth(i)/rates%ada
                            end do
                        end if
                        iskip = iskip + 1
                        if (iskip == nskip) then
                            iskip = 0
                        end if
                    end if

                    t = t + sys%dt

                    !save intemediate steps ?
                    if (ip >= sys%np) then

                        !gp 17-feb-05
                        !call save_a_step(nr, intg,pr, sys%dt, sys%imethph, rates)
                        call save_a_step(nr, intg,pr, sys%dt, sys%imethph, sys%showdielresults, rates)

                    end if
                end do
            else if (sys%IMeth == 'Runge-Kutta') then
                write(*,*) 'day', ip
                do ic = 1, sys%nc !time step in each day

                    !gp 08-jan-10
                    !call rk4(nr, meteo, solar,hw, db, hydrau, sys, intg, &
                    ! rates, dte, dc, dinb, dipb, t, sys%dt)
                    call rk4(nr, meteo, solar,hw, db, hydrau, sys, intg, &
                        rates, dte, dc, dinb, dipb, t, sys%dt, sys%statevariables)

                    !gp 24-jun-09 write dynamic output every 45 minutes
                    if (sys%writedynamic == "Yes") then
                        if (iskip == 0) then
                            do i=0, nr
                                j = 1
                                tss = intg%c(i, 11, j) * rates%ada + intg%c(i, 2, j) + &
                                    intg%c(i, 12, j)
                                tp = intg%c(i, 11, j) * rates%apa + intg%c(i, 9, j) + &
                                    intg%c(i, 10, j)
                                tn = intg%c(i, 11, j) * rates%ana + intg%c(i, 6, j) + &
                                    intg%c(i, 7, j) + intg%c(i, 8, j)
                                dosat = oxygen_saturation(intg%te(i, j), hydrau%reach(i)%elev)

                                !gp 23-nov-09
                                !if (sys%imethph == "newton-raphson") then
                                !call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !else
                                !call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !end if
                                if (sys%imethph == "Newton-Raphson") then
                                    call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                elseif (sys%imethph == "Bisection") then
                                    call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                else
                                    call phsolbrent(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                end if

                                nh3 = 1.0_r64/(1 + 10.0_r64 ** (-ph)/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                                    (intg%te(i, j) + 273.15_r64))) * intg%c(i, 7, j)
                                write(12, '(I13, 41F13.4)') i, t, intg%te(i, j), &
                                    (intg%c(i, k, j), k=1, nv-2), ph, &
                                    intg%c(i, nv, j), intg%c(i, nv, j)* rates%mga / rates%mgd * 1000, &
                                    tss, tp, tn , dosat, nh3, &
                                    intg%inb(i), &
                                    intg%inb(i), &
                                    phitsave(i), philsave(i), &
                                    phinsave(i), phipsave(i), &
                                    phicsave(i), phitotalsave(i), &
                                    savebotalgphoto(i), savebotalgresp(i), savebotalgdeath(i), savebotalgnetgrowth(i), &
                                    savebotalgphoto(i)/rates%ada, savebotalgresp(i)/rates%ada, &
                                    savebotalgdeath(i)/rates%ada, savebotalgnetgrowth(i)/rates%ada
                            end do
                        end if
                        iskip = iskip + 1
                        if (iskip == nskip) then
                            iskip = 0
                        end if
                    end if

                    t = t + sys%dt

                    !save intemediate steps ?
                    if (ip >= sys%np) then

                        !gp 17-feb-05
                        !call save_a_step(nr, intg,pr, sys%dt, sys%imethph, rates)
                        call save_a_step(nr, intg,pr, sys%dt, sys%imethph, sys%showdielresults, rates)

                    end if
                end do
            else if (sys%imeth == "Adaptive step") then
                write(*,*) 'day', ip
                if (ip >= sys%np) then

                    !gp 08-jan-10
                    !call odeint(intg, sys, rates, meteo, solar, &
                    ! hw, db, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .true.)
                    call odeint(intg, sys, rates, meteo, solar, &
                        hw, db, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .true., sys%statevariables)

                else

                    !gp 08-jan-10
                    !call odeint(intg, sys, rates, meteo, solar, &
                    ! hw, db, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .false.)
                    call odeint(intg, sys, rates, meteo, solar, &
                        hw, db, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .false., sys%statevariables)

                end if
            else
                stop 'Please specify integration method'
            end if
        end do

    end subroutine

    !-----------------------------------------------------------------------------------------
    !save the initial step of the last simulation day
    !called only once

    !gp 17-feb-05
    !subroutine save_init_step(nr, intg, pr, dt, imethph, rates)
    subroutine save_init_step(nr, intg, pr, dt, imethph, showdielresults, rates)

        ! use class_output, only: outdata_type

        integer(i32), intent(in) :: nr

        real(r64), intent(in) :: dt
        character(len=30), intent(in) :: imethph ! ph solve method

        !gp 17-feb-05
        character(len=30), intent(in) :: showdielresults !yes or no (only used in excel vba)

        type(rates_t), intent(in) :: rates
        type(integral_type), intent(in) :: intg
        type(outdata_t), intent(out) :: pr

        integer(i32) i, j, k
        !gp 27-oct-04 real(r64) :: phss(0:nr), kamm, ph
        real(r64) :: phss(0:nr, 0:nl), kamm, ph !gp add dim for nl


        ! !gp debug
        ! open (unit=9, file='debug.out', status='replace', action='write')
        ! write(9,*) 'program is now just before type(riverhydraulics_type) hydrau'
        ! close (9)

        !gp 29-oct-04
        type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters

        ! !gp debug
        ! open (unit=9, file='debug.out', status='replace', action='write')
        ! write(9,'(a64)') 'program is now just after type(riverhydraulics_type) hydrau'
        ! close (9)

        pr%tdy(0)=0
        do i = 0, nr !nr public
            !gp 27-oct-04 add nl dimension
            !gp pr%temn(i) = intg%te(i, 1)
            !gp pr%temx(i) = intg%te(i, 1)
            !gp pr%teav(i) = 0
            !gp pr%osav(i) = 0
            !gp pr%phsav(i) = 0
            !gp do k = 1, nv
            !gp pr%cmn(i, k) = intg%c(i, k)
            !gp pr%cmx(i, k) = intg%c(i, k)
            !gp end do
            !gp if (imethph == "newton-raphson") then !public imethph as string
            !gp call phsolnewton(ph, intg%c(i, 15), intg%te(i, 1), intg%c(i, 14), &
            !gp intg%c(i, 1))
            !gp else
            !gp call phsolbisect(ph, intg%c(i, 15), intg%te(i, 1), intg%c(i, 14), &
            !gp intg%c(i, 1))
            !gp end if
            !gp pr%phav(i) = 0
            !gp pr%phmn(i) = ph
            !gp pr%phmx(i) = ph
            !gp phss(i) = ph
            !gp kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%te(i, 1) + 273.15_r64)))
            !gp pr%nh3av(i) = 0
            !gp pr%nh3mn(i) = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-ph) / kamm) * intg%c(i, 7)
            !gp pr%nh3mx(i) = pr%nh3mn(i)
            !gp pr%tpav(i) = 0
            !gp pr%tpmn(i) = intg%c(i, 9) + intg%c(i, 10) + intg%c(i, 11) * rates%apa
            !gp pr%tpmx(i) = pr%tpmn(i)
            !gp pr%tnav(i) = 0
            !gp pr%tnmn(i) = intg%c(i, 6) + intg%c(i, 7) + intg%c(i, 8) + &
            !gp intg%c(i, 11) * rates%ana
            !gp pr%tnmx(i) = pr%tnmn(i)
            do j = 1, nl
                pr%temn(i, j) = intg%te(i, j)
                pr%temx(i, j) = intg%te(i, j)
                pr%teav(i, j) = 0
                pr%osav(i, j) = 0
                pr%phsav(i, j) = 0
                do k = 1, nv
                    pr%cmn(i, k, j) = intg%c(i, k, j)
                    pr%cmx(i, k, j) = intg%c(i, k, j)
                end do

                !gp 23-nov-09
                !if (imethph == "newton-raphson") then !public imethph as string
                ! call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !else
                !call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !end if
                if (imethph == "Newton-Raphson") then
                    call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                elseif (imethph == "Bisection") then
                    call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                else
                    call phsolbrent(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                end if

                pr%phav(i, j) = 0
                pr%phmn(i, j) = ph
                pr%phmx(i, j) = ph
                phss(i, j) = ph
                kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%te(i, j) + 273.15_r64)))
                pr%nh3av(i, j) = 0
                pr%nh3mn(i, j) = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-ph) / kamm) * intg%c(i, 7, j)
                pr%nh3mx(i, j) = pr%nh3mn(i, j)
                pr%tpav(i, j) = 0
                pr%tpmn(i, j) = intg%c(i, 9, j) + intg%c(i, 10, j) + intg%c(i, 11, j) * rates%apa
                pr%tpmx(i, j) = pr%tpmn(i, j)
                pr%tnav(i, j) = 0
                pr%tnmn(i, j) = intg%c(i, 6, j) + intg%c(i, 7, j) + intg%c(i, 8, j) + &
                    intg%c(i, 11, j) * rates%ana
                pr%tnmx(i, j) = pr%tnmn(i, j) !gp end of new block of code 27-oct-04
            end do
            !'gp 15-nov-04 reach average daily average sediment fluxes
            pr%hypofluxdoav(i) = 0; pr%hypofluxcbodav(i) = 0; pr%hypofluxnh4av(i) = 0
            pr%hypofluxno3av(i) = 0; pr%hypofluxsrpav(i) = 0
            pr%diagfluxdoav(i) = 0; pr%diagfluxcbodav(i) = 0; pr%diagfluxnh4av(i) = 0
            pr%diagfluxno3av(i) = 0; pr%diagfluxsrpav(i) = 0
            !'gp 11-jan-05 min/max/mean cell quotas mgn/gd and mgp/gd
            if (intg%c(i, nv, 1) > 0) then
                pr%ninbav(i) = 0
                pr%ninbmn(i) = intg%inb(i) / intg%c(i, nv, 1)
                pr%ninbmx(i) = intg%inb(i) / intg%c(i, nv, 1)
                pr%nipbav(i) = 0
                pr%nipbmn(i) = intg%ipb(i) / intg%c(i, nv, 1)
                pr%nipbmx(i) = intg%ipb(i) / intg%c(i, nv, 1)
            else
                pr%ninbav(i) = 0; pr%ninbmn(i) = 0; pr%ninbmx(i) = 0
                pr%nipbav(i) = 0; pr%nipbmn(i) = 0; pr%nipbmx(i) = 0
            end if

            !gp 25-jun-09
            pr%av_botalgphoto(i)=0; pr%av_botalgresp(i)=0; pr%av_botalgdeath(i)=0; pr%av_botalgnetgrowth(i)=0

        end do

        !gp 17-feb-05
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if (showdielresults == "Yes") then
            !x x x x x x x x x x x x x x x x x x x

            pr%nj = 0
            do i = 0, nr
                do j = 1, nl
                    pr%tepr(i, pr%nj, j) = intg%te(i, j)
                    !gp 27-oct-04 end do
                    !gp do k = 1, nv
                    !gp pr%cpr(i, k, pr%nj) = intg%c(i, k)
                    !gp end do
                    !gp pr%phpr(i, pr%nj) = phss(i)
                    do k = 1, nv
                        pr%cpr(i, k, pr%nj, j) = intg%c(i, k, j)
                    end do
                    pr%phpr(i, pr%nj, j) = phss(i, j)
                end do !gp end of new block of code 27-oct-04
            end do
            do i = 1, nr
                !gp 02-nov-04 pr%ninbpr(i, pr%nj) = intg%inb(i) / intg%c(i, nv) !?
                !gp pr%nipbpr(i, pr%nj) = intg%ipb(i) / intg%c(i, nv) !?
                pr%ninbpr(i, pr%nj) = intg%inb(i) / intg%c(i, nv, 1)
                pr%nipbpr(i, pr%nj) = intg%ipb(i) / intg%c(i, nv, 1) !gp 02-nov-04 end new block
                !gp 20-oct-04 output growth limitation factors for bottom algae
                pr%phitotalsavepr(i, pr%nj) = phitotalsave(i)
                pr%phitsavepr(i, pr%nj) = phitsave(i)
                pr%philsavepr(i, pr%nj) = philsave(i)
                pr%phinsavepr(i, pr%nj) = phinsave(i)
                pr%phipsavepr(i, pr%nj) = phipsave(i)
                pr%phicsavepr(i, pr%nj) = phicsave(i)
                !gp 28-oct-04 output hyporheic fluxes of constituents (total surface area)
                pr%hypofluxdopr(i, pr%nj) = hypofluxdo(i) !'dissolved oxygen go2/m^2/d
                pr%hypofluxcbodpr(i, pr%nj) = hypofluxcbod(i) !'fast cbod go2/m^2/d
                pr%hypofluxnh4pr(i, pr%nj) = hypofluxnh4(i) !'ammonia mgn/m^2/d
                pr%hypofluxno3pr(i, pr%nj) = hypofluxno3(i) !'nitrate mgn/m^2/d
                pr%hypofluxsrppr(i, pr%nj) = hypofluxsrp(i) !'srp mgp/m^2/d
                pr%hypofluxicpr(i, pr%nj) = hypofluxic(i) !'inorganic c gc/m^2/d
                !gp 01-nov-04 output diagenesis fluxes of constituents (total surface area)
                pr%diagfluxdopr(i, pr%nj) = diagfluxdo(i) !'dissolved oxygen go2/m^2/d
                pr%diagfluxcbodpr(i, pr%nj) = diagfluxcbod(i) !'fast cbod go2/m^2/d
                pr%diagfluxnh4pr(i, pr%nj) = diagfluxnh4(i) !'ammonia mgn/m^2/d
                pr%diagfluxno3pr(i, pr%nj) = diagfluxno3(i) !'nitrate mgn/m^2/d
                pr%diagfluxsrppr(i, pr%nj) = diagfluxsrp(i) !'srp mgp/m^2/d
                pr%diagfluxicpr(i, pr%nj) = diagfluxic(i) !'inorganic c gc/m^2/d


                !'gp 05-jul-05
                !'heat fluxes (converted from cal/cm^2/d to w/m^2)
                pr%pr_saveheatfluxjsnt(i, pr%nj) = saveheatfluxjsnt(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxlongat(i, pr%nj) = saveheatfluxlongat(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxback(i, pr%nj) = saveheatfluxback(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxconv(i, pr%nj) = saveheatfluxconv(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxevap(i, pr%nj) = saveheatfluxevap(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxjsed(i, pr%nj) = saveheatfluxjsed(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxjhyporheic(i, pr%nj) = saveheatfluxjhyporheic(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxtribs(i, pr%nj) = saveheatfluxtribs(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxadvecdisp(i, pr%nj) = saveheatfluxadvecdisp(i) * (4.183076 * 100 * 100 / 86400)
                !'do fluxes (go2/m^2/d)
                pr%pr_savedofluxreaer(i, pr%nj) = savedofluxreaer(i)
                pr%pr_savedofluxcbodfast(i, pr%nj) = savedofluxcbodfast(i)
                pr%pr_savedofluxcbodslow(i, pr%nj) = savedofluxcbodslow(i)
                pr%pr_savedofluxnitrif(i, pr%nj) = savedofluxnitrif(i)
                pr%pr_savedofluxphytoresp(i, pr%nj) = savedofluxphytoresp(i)
                pr%pr_savedofluxphytophoto(i, pr%nj) = savedofluxphytophoto(i)
                pr%pr_savedofluxbotalgresp(i, pr%nj) = savedofluxbotalgresp(i)
                pr%pr_savedofluxbotalgphoto(i, pr%nj) = savedofluxbotalgphoto(i)
                pr%pr_savedofluxsod(i, pr%nj) = savedofluxsod(i)
                pr%pr_savedofluxcod(i, pr%nj) = savedofluxcod(i)
                pr%pr_savedofluxhyporheic(i, pr%nj) = hypofluxdo(i)
                pr%pr_savedofluxtribs(i, pr%nj) = savedofluxtribs(i)
                pr%pr_savedofluxadvecdisp(i, pr%nj) = savedofluxadvecdisp(i)
                !'co2 fluxes (gc/m^2/d)
                pr%pr_saveco2fluxreaer(i, pr%nj) = saveco2fluxreaer(i)
                pr%pr_saveco2fluxcbodfast(i, pr%nj) = saveco2fluxcbodfast(i)
                pr%pr_saveco2fluxcbodslow(i, pr%nj) = saveco2fluxcbodslow(i)
                pr%pr_saveco2fluxphytoresp(i, pr%nj) = saveco2fluxphytoresp(i)
                pr%pr_saveco2fluxphytophoto(i, pr%nj) = saveco2fluxphytophoto(i)
                pr%pr_saveco2fluxbotalgresp(i, pr%nj) = saveco2fluxbotalgresp(i)
                pr%pr_saveco2fluxbotalgphoto(i, pr%nj) = saveco2fluxbotalgphoto(i)
                pr%pr_saveco2fluxsod(i, pr%nj) = saveco2fluxsod(i)
                pr%pr_saveco2fluxhyporheic(i, pr%nj) = hypofluxic(i)
                pr%pr_saveco2fluxtribs(i, pr%nj) = saveco2fluxtribs(i)
                pr%pr_saveco2fluxadvecdisp(i, pr%nj) = saveco2fluxadvecdisp(i)


            end do

            !gp 17-feb-05
            !x x x x x x x x x x x x x x x x x x x
        end if
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    end subroutine save_init_step

!--------------------------------------------------------------------------------------------
    !save an intermediate step during the last day of the simulation

    !gp 17-feb-05
    !subroutine save_a_step(nr, intg, pr, dt, imethph, rates)
    subroutine save_a_step(nr, intg, pr, dt, imethph, showdielresults, rates)


        integer(i32), intent(in) :: nr
        real(r64), intent(in) :: dt
        type(rates_t), intent(in) :: rates
        !character(len=30) :: imeth = 'euler' ! integration method
        character(len=30), intent(in) :: imethph ! ph solve method

        !gp 17-feb-05
        character(len=30), intent(in) :: showdielresults !yes or no (only used in excel vba)

        type(integral_type), intent(in) :: intg
        type(outdata_t), intent(inout) :: pr
        integer(i32) i, j, k
        real(r64) kamm, nh3, ph, totp, totn
        real(r64) k1, k2, kw, kh, co2sat
        !gp 28-oct-04 real(r64) :: phss(0:nr)
        real(r64) :: phss(0:nr, nl) !gp

        !gp 29-oct-04
        type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters

        do i = 0, nr

            !gp 27-oct-04
            !gp !update maximum and minimum water temperature
            do j = 1, nl
                if (intg%te(i, j) < pr%temn(i, j)) then
                    pr%temn(i, j) = intg%te(i, j)
                else if (intg%te(i, j) > pr%temx(i, j)) then
                    pr%temx(i, j) = intg%te(i, j)
                end if
                !update average water temperature
                pr%teav(i, j) = pr%teav(i, j) + intg%te(i, j) * dt
                !update average saturation do
                pr%osav(i, j) = pr%osav(i, j) + os(i, j) * dt
                call chemrates(intg%te(i, j), k1, k2, kw, kh, intg%c(i, 1, j))
                co2sat = kh * rates%pco2 !co2 saturation concentration
                call modfp2(4.0_r64, 12.0_r64, ph, co2sat, intg%te(i, j), intg%c(i, nv - 2, j), &
                    intg%c(i, 1, j)) !solving ph under co2 saturation condition
                !update average ph under co2 saturation condition
                pr%phsav(i, j) = pr%phsav(i, j) + ph * dt
                !
                do k = 1, nv
                    !update constituents max and min concentrations or values
                    if (intg%c(i, k, j) < pr%cmn(i, k, j)) then
                        pr%cmn(i, k, j) = intg%c(i, k, j)
                    else if (intg%c(i, k, j) > pr%cmx(i, k, j)) then
                        pr%cmx(i, k, j) = intg%c(i, k, j)
                    end if
                    !update average constituents concentrations or values
                    pr%cav(i, k, j) = pr%cav(i, k, j) + intg%c(i, k, j) * dt
                end do

                !solve ph value for now

                !gp 23-nov-09
                !if (imethph == "newton-raphson") then
                ! call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !else
                ! call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !end if
                if (imethph == "Newton-Raphson") then
                    call phsolnewton(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                elseif (imethph == "Bisection") then
                    call phsolbisect(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                else
                    call phsolbrent(ph, intg%c(i, nv - 1, j), intg%te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                end if

                phss(i, j) = ph
                if (ph < pr%phmn(i, j)) then
                    pr%phmn(i, j) = ph
                else if (ph > pr%phmx(i, j)) then
                    pr%phmx(i, j) = ph
                end if
                !update average ph value
                pr%phav(i, j) = pr%phav(i, j) + ph * dt

                !// nh3 //
                kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%te(i, j) + 273.15_r64)))
                nh3 = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-ph) / kamm) * intg%c(i, 7, j)
                if (nh3 < pr%nh3mn(i, j)) then
                    pr%nh3mn(i, j) = nh3
                else if (nh3 > pr%nh3mx(i, j)) then
                    pr%nh3mx(i, j) = nh3
                end if
                !average value of nh3
                pr%nh3av(i, j) = pr%nh3av(i, j) + nh3 * dt

                !// total phosphorus //
                totp = intg%c(i, 9, j) + intg%c(i, 10, j) + intg%c(i, 11, j) * rates%apa
                pr%tpav(i, j) = pr%tpav(i, j) + totp * dt
                if (totp < pr%tpmn(i, j)) then
                    pr%tpmn(i, j) = totp
                else if (totp > pr%tpmx(i, j)) then
                    pr%tpmx(i, j) = totp
                end if

                !// total nitrogen //
                totn = intg%c(i, 6, j) + intg%c(i, 7, j) + intg%c(i, 8, j) + &
                    intg%c(i, 11, j) * rates%ana
                pr%tnav(i, j) = pr%tnav(i, j) + totn * dt
                if (totn < pr%tnmn(i, j)) then
                    pr%tnmn(i, j) = totn
                elseif (totn > pr%tnmx(i, j)) then
                    pr%tnmx(i, j) = totn
                end if
            end do !gp 27-oct-04 end block of new code
            !'gp 15-nov-04 reach average daily average sediment fluxes
            if (i > 0) then
                pr%hypofluxdoav(i) = pr%hypofluxdoav(i) + hypofluxdo(i) * dt !'dissolved oxygen go2/m^2/d
                pr%hypofluxcbodav(i) = pr%hypofluxcbodav(i) + hypofluxcbod(i) * dt !'fast cbod go2/m^2/d
                pr%hypofluxnh4av(i) = pr%hypofluxnh4av(i) + hypofluxnh4(i) * dt !'ammonia mgn/m^2/d
                pr%hypofluxno3av(i) = pr%hypofluxno3av(i) + hypofluxno3(i) * dt !'nitrate mgn/m^2/d
                pr%hypofluxsrpav(i) = pr%hypofluxsrpav(i) + hypofluxsrp(i) * dt !'srp mgp/m^2/d
                pr%diagfluxdoav(i) = pr%diagfluxdoav(i) + diagfluxdo(i) * dt !'dissolved oxygen go2/m^2/d
                pr%diagfluxcbodav(i) = pr%diagfluxcbodav(i) + diagfluxcbod(i) * dt !'fast cbod go2/m^2/d
                pr%diagfluxnh4av(i) = pr%diagfluxnh4av(i) + diagfluxnh4(i) * dt !'ammonia mgn/m^2/d
                pr%diagfluxno3av(i) = pr%diagfluxno3av(i) + diagfluxno3(i) * dt !'nitrate mgn/m^2/d
                pr%diagfluxsrpav(i) = pr%diagfluxsrpav(i) + diagfluxsrp(i) * dt !'srp mgp/m^2/d
                !'gp 11-jan-05 mgn/gd and mgp/gd
                if (intg%c(i, nv, 1) > 0) then
                    pr%ninbav(i) = pr%ninbav(i) + dt * intg%inb(i) / intg%c(i, nv, 1)
                    if (intg%inb(i) / intg%c(i, nv, 1) < pr%ninbmn(i)) pr%ninbmn(i) = intg%inb(i) / intg%c(i, nv, 1)
                    if (intg%inb(i) / intg%c(i, nv, 1) > pr%ninbmx(i)) pr%ninbmx(i) = intg%inb(i) / intg%c(i, nv, 1)
                    pr%nipbav(i) = pr%nipbav(i) + dt * intg%ipb(i) / intg%c(i, nv, 1)
                    if (intg%ipb(i) / intg%c(i, nv, 1) < pr%nipbmn(i)) pr%nipbmn(i) = intg%ipb(i) / intg%c(i, nv, 1)
                    if (intg%ipb(i) / intg%c(i, nv, 1) > pr%nipbmx(i)) pr%nipbmx(i) = intg%ipb(i) / intg%c(i, nv, 1)
                end if

                !gp 25-jun-09
                pr%av_botalgphoto(i) = pr%av_botalgphoto(i) + savebotalgphoto(i) * dt !average botalgphoto/asb gd/m2/d
                pr%av_botalgresp(i) = pr%av_botalgresp(i) + savebotalgresp(i) * dt !average botalgreso/asb gd/m2/d
                pr%av_botalgdeath(i) = pr%av_botalgdeath(i) + savebotalgdeath(i) * dt !average botalgdeath/asb gd/m2/d
                pr%av_botalgnetgrowth(i) = pr%av_botalgnetgrowth(i) + savebotalgnetgrowth(i) * dt !average botalg net growth (photo-resp-death) gd/m2/d

            end if !'gp 15-nov-04 end new block
        end do

        !gp 17-feb-05
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if (showdielresults == "Yes") then
            !x x x x x x x x x x x x x x x x x x x

            pr%nj = pr%nj + 1
            pr%tdy(pr%nj)= pr%tdy(pr%nj-1) + dt
            do i = 0, nr
                do j = 1, nl
                    pr%tepr(i, pr%nj, j) = intg%te(i, j)
                    do k = 1, nv
                        pr%cpr(i, k, pr%nj, j) = intg%c(i, k, j)
                    end do
                    pr%phpr(i, pr%nj, j) = phss(i, j)
                end do !gp 27-oct-04
            end do
            do i = 1, nr
                !gp 02-nov-04 pr%ninbpr(i, pr%nj) = intg%inb(i) / intg%c(i, nv)
                !gp pr%nipbpr(i, pr%nj) = intg%ipb(i) / intg%c(i, nv)
                pr%ninbpr(i, pr%nj) = intg%inb(i) / intg%c(i, nv, 1)
                pr%nipbpr(i, pr%nj) = intg%ipb(i) / intg%c(i, nv, 1) !gp 02-nov-04 end new block
                pr%phitotalsavepr(i, pr%nj) = phitotalsave(i) !gp 20-oct-04
                pr%phitsavepr(i, pr%nj) = phitsave(i) !gp 20-oct-04
                pr%philsavepr(i, pr%nj) = philsave(i) !gp 20-oct-04
                pr%phinsavepr(i, pr%nj) = phinsave(i) !gp 20-oct-04
                pr%phipsavepr(i, pr%nj) = phipsave(i) !gp 20-oct-04
                pr%phicsavepr(i, pr%nj) = phicsave(i) !gp 20-oct-04
                !gp 28-oct-04 output hyporheic fluxes of constituents (total surface area)
                pr%hypofluxdopr(i, pr%nj) = hypofluxdo(i) !'dissolved oxygen go2/m^2/d
                pr%hypofluxcbodpr(i, pr%nj) = hypofluxcbod(i) !'fast cbod go2/m^2/d
                pr%hypofluxnh4pr(i, pr%nj) = hypofluxnh4(i) !'ammonia mgn/m^2/d
                pr%hypofluxno3pr(i, pr%nj) = hypofluxno3(i) !'nitrate mgn/m^2/d
                pr%hypofluxsrppr(i, pr%nj) = hypofluxsrp(i) !'srp mgp/m^2/d
                pr%hypofluxicpr(i, pr%nj) = hypofluxic(i) !'inorganic c gc/m^2/d
                !gp 01-nov-04 output diagenesis fluxes of constituents (total surface area)
                pr%diagfluxdopr(i, pr%nj) = diagfluxdo(i) !'dissolved oxygen go2/m^2/d
                pr%diagfluxcbodpr(i, pr%nj) = diagfluxcbod(i) !'fast cbod go2/m^2/d
                pr%diagfluxnh4pr(i, pr%nj) = diagfluxnh4(i) !'ammonia mgn/m^2/d
                pr%diagfluxno3pr(i, pr%nj) = diagfluxno3(i) !'nitrate mgn/m^2/d
                pr%diagfluxsrppr(i, pr%nj) = diagfluxsrp(i) !'srp mgp/m^2/d
                pr%diagfluxicpr(i, pr%nj) = diagfluxic(i) !'inorganic c gc/m^2/d

                !'gp 05-jul-05
                !'heat fluxes (converted from cal/cm^2/d to w/m^2)
                pr%pr_saveheatfluxjsnt(i, pr%nj) = saveheatfluxjsnt(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxlongat(i, pr%nj) = saveheatfluxlongat(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxback(i, pr%nj) = saveheatfluxback(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxconv(i, pr%nj) = saveheatfluxconv(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxevap(i, pr%nj) = saveheatfluxevap(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxjsed(i, pr%nj) = saveheatfluxjsed(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxjhyporheic(i, pr%nj) = saveheatfluxjhyporheic(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxtribs(i, pr%nj) = saveheatfluxtribs(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveheatfluxadvecdisp(i, pr%nj) = saveheatfluxadvecdisp(i) * (4.183076 * 100 * 100 / 86400)
                !'do fluxes (go2/m^2/d)
                pr%pr_savedofluxreaer(i, pr%nj) = savedofluxreaer(i)
                pr%pr_savedofluxcbodfast(i, pr%nj) = savedofluxcbodfast(i)
                pr%pr_savedofluxcbodslow(i, pr%nj) = savedofluxcbodslow(i)
                pr%pr_savedofluxnitrif(i, pr%nj) = savedofluxnitrif(i)
                pr%pr_savedofluxphytoresp(i, pr%nj) = savedofluxphytoresp(i)
                pr%pr_savedofluxphytophoto(i, pr%nj) = savedofluxphytophoto(i)
                pr%pr_savedofluxbotalgresp(i, pr%nj) = savedofluxbotalgresp(i)
                pr%pr_savedofluxbotalgphoto(i, pr%nj) = savedofluxbotalgphoto(i)
                pr%pr_savedofluxsod(i, pr%nj) = savedofluxsod(i)
                pr%pr_savedofluxcod(i, pr%nj) = savedofluxcod(i)
                pr%pr_savedofluxhyporheic(i, pr%nj) = hypofluxdo(i)
                pr%pr_savedofluxtribs(i, pr%nj) = savedofluxtribs(i)
                pr%pr_savedofluxadvecdisp(i, pr%nj) = savedofluxadvecdisp(i)
                !'co2 fluxes (gc/m^2/d)
                pr%pr_saveco2fluxreaer(i, pr%nj) = saveco2fluxreaer(i)
                pr%pr_saveco2fluxcbodfast(i, pr%nj) = saveco2fluxcbodfast(i)
                pr%pr_saveco2fluxcbodslow(i, pr%nj) = saveco2fluxcbodslow(i)
                pr%pr_saveco2fluxphytoresp(i, pr%nj) = saveco2fluxphytoresp(i)
                pr%pr_saveco2fluxphytophoto(i, pr%nj) = saveco2fluxphytophoto(i)
                pr%pr_saveco2fluxbotalgresp(i, pr%nj) = saveco2fluxbotalgresp(i)
                pr%pr_saveco2fluxbotalgphoto(i, pr%nj) = saveco2fluxbotalgphoto(i)
                pr%pr_saveco2fluxsod(i, pr%nj) = saveco2fluxsod(i)
                pr%pr_saveco2fluxhyporheic(i, pr%nj) = hypofluxic(i)
                pr%pr_saveco2fluxtribs(i, pr%nj) = saveco2fluxtribs(i)
                pr%pr_saveco2fluxadvecdisp(i, pr%nj) = saveco2fluxadvecdisp(i)

            end do

            !gp 17-feb-05
            !x x x x x x x x x x x x x x x x x x x
        end if
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    end subroutine save_a_step

    !----------------------------------------------------------------------------------------------
    !make appropriate step
    !Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust
    !stepsize.

    !gp 08-jan-10
    !subroutine odeint (begin, sys, rates, meteo, solar, &
    ! hw, db, hydrau, pr, nr, t1, t2, dt1, savesteps)
    subroutine odeint (begin, sys, rates, meteo, solar, &
        hw, db, hydrau, pr, nr, t1, t2, dt1, savesteps, statevariables)

        !constants
        real(r64), parameter :: tiny=1.0e-30
        integer(i32), parameter :: maxstp=1000
        !dummy variables
        integer(i32), intent(in) :: nr

        !gp 08-jan-10
        character(len=30), intent(in) :: statevariables !to test for 'all except temperature'

        type(meteorology_t) :: meteo
        type(solar_type) solar !solar radiation
        ! type(headwaterdownstream_type), intent(in) :: hdboundary
        type(upstream_boundary_t) hw !headwater
        type(downstream_boundary_t) db !downstream boundary
        type(riverhydraulics_type), intent(in) :: hydrau !channel dimensions, hydraulics, physical characters
        type(systemparams) sys
        type(integral_type), intent(inout) :: begin
        type(rates_t), intent(in) :: rates
        type(outdata_t), intent(out) :: pr !output data
        real(r64), intent(in) :: t1, t2, dt1
        logical(4), intent(in) :: savesteps
        !local variables
        real(r64) t, dt !initial time
        integer(i32) nstp, i, j
        real(r64) dtnext, dtdid
        type(integral_type) scal, now
        real(r64), dimension(nr,nl) :: tenow
        real(r64), dimension(nr) :: inbnow, ipbnow
        !gp 28-oct-04 real(r64) dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64) dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp
        integer(i32) kmax, kout

        t=t1 !initial time
        dt=dt1

        !constructor
        now = integration_(nr)
        scal= integration_(nr)
        !initialize
        now = begin

        do nstp = 1, maxstp
            ! write(*,*) t
            call derivs(nr, meteo, solar, hw, db, hydrau, sys, now%te, now%c, now%inb, now%ipb, &
                rates, dte, dc, dinb, dipb, t)

            !scaling to monitor accuracy.
            !gp 28-oct-04 scal%c(1:nr, 1:nv) = abs(now%c(1:nr, 1:nv)) + abs(dc(1:nr, 1:nv) * dt) + tiny
            scal%c(1:nr, 1:nv, 1:nl) = abs(now%c(1:nr, 1:nv, 1:nl)) + abs(dc(1:nr, 1:nv, 1:nl) * dt) + tiny !gp
            scal%te(1:nr, 1:nl) = abs(now%te(1:nr, 1:nl)) +abs(dte(1:nr, 1:nl) *dt) + tiny
            scal%inb(1:nr) = abs(now%inb(1:nr)) +abs(dinb(1:nr) *dt) + tiny
            scal%ipb(1:nr) = abs(now%ipb(1:nr)) +abs(dipb(1:nr) *dt) + tiny
            ! if stepsize can overshoot, decrease
            if (t + dt > t2) dt = t2 - t
            !adaptive step

            !gp 08-jan-10
            !call rkqs(nr, meteo, solar, hw, db, hydrau, sys, now, &
            ! scal, rates, dte, dc, dinb, dipb, t, dt, dtdid,dtnext)
            call rkqs(nr, meteo, solar, hw, db, hydrau, sys, now, &
                scal, rates, dte, dc, dinb, dipb, t, dt, dtdid, dtnext, statevariables)

            if (savesteps) then

                !gp 17-feb-05
                !call save_a_step(nr, now, pr, dtdid, sys%imethph, rates)
                call save_a_step(nr, now, pr, dtdid, sys%imethph, sys%showdielresults, rates)

            end if
            if (t - t2 >= 0) then !are we done?
                begin= now
                return
            end if
            dt = dtnext
        end do

    end subroutine odeint

!-----------------------------------------------------------------------------------------

    !gp 08-jan-10
    !subroutine rkqs(nr, meteo, solar, hw, db, hydrau, sys, intg,&
    ! scal, rates, dte, dc, dinb, dipb, t, dttry, dtdid, dtnext)
    subroutine rkqs(nr, meteo, solar, hw, db, hydrau, sys, intg,&
        scal, rates, dte, dc, dinb, dipb, t, dttry, dtdid, dtnext, statevariables)


        real(r64), parameter :: safety=0.9
        real(r64), parameter :: pgrow=-0.2
        real(r64), parameter :: pshrnk=-0.25
        real(r64), parameter :: errcon=0.000189
        real(r64), parameter :: epstol=0.0001

        integer(i32), intent(in) :: nr

        !gp 08-jan-10
        character(len=30), intent(in) :: statevariables !to test for 'all except temperature'

        type(meteorology_t) :: meteo
        type(solar_type) solar !solar radiation
        type(upstream_boundary_t), intent(in) :: hw !headwater
        type(downstream_boundary_t), intent(in) :: db !downstream boundary
        type(riverhydraulics_type), intent(in) :: hydrau !channel dimensions, hydraulics, physical characters
        type(systemparams) sys
        type(integral_type), intent(inout) :: intg
        type(integral_type), intent(in) :: scal
        type(rates_t), intent(in) :: rates
        !gp 28-oct-04 real(r64), intent(in):: dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64), intent(in):: dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp
        real(r64), intent(inout) :: t !initial time, double
        real(r64), intent(in) :: dttry
        real(r64), intent(out) :: dtdid, dtnext
        integer(i32) i, j, k
        real(r64) errmax, dt, dttemp, tnew

        type(integral_type) tmp, err !current value and incremented value

        tmp = integration_(nr)
        err = integration_(nr)

        !set the stepsize to the initial trial value
        dt = dttry
        do while (.true.)
            !take a step

            !gp 08-jan-10
            !call rkck(nr, meteo, solar, hw, db, hydrau, sys, intg, &
            ! rates, dte, dc, dinb, dipb, t, dt, tmp, err)
            call rkck(nr, meteo, solar, hw, db, hydrau, sys, intg, &
                rates, dte, dc, dinb, dipb, t, dt, tmp, err, statevariables)

            !gp 28-oct-04 errmax = 0.0
            !gp !evaluate accuracy
            !gp do i = 1 , nr
            !gp do j = 1 , nv -1
            !gp errmax = max(errmax, abs(err%c(i, j) / scal%c(i, j)))
            !gp ! write(*,*) i, j, errmax
            !gp end do
            !gp end do
            !gp
            !gp ! write (*,*) 'temperature'
            !gp do i = 1 , nr
            !gp do k=1, nl
            !gp errmax= max(errmax, abs(err%te(i,k)/scal%te(i,k)))
            !gp ! write(*,*) i,k, errmax
            !gp end do
            !gp ! errmax= max(errmax, abs(err%inb(i)/scal%inb(i)))
            !gp ! errmax= max(errmax, abs(err%ipb(i)/scal%ipb(i)))
            !gp end do
            !evaluate accuracy
            errmax = 0.0
            do i = 1 , nr
                do j = 1, nl
                    errmax= max(errmax, abs(err%te(i, j)/scal%te(i, j)))
                    do k = 1 , nv -1
                        errmax = max(errmax, abs(err%c(i, k, j) / scal%c(i, k, j)))
                    end do
                end do
            end do !gp end new block

            !scale relative to required tolerance
            errmax = errmax / epstol
            if (errmax <= 1.0) then
                exit !step succeeded, compute size of next step
            end if
            dttemp = safety * dt * (errmax ** pshrnk)

            !truncation error too large, reduce stepsize
            dt = max(dttemp, 0.1 * dt) !no more than a factor of 10.

            tnew = t + dt
            if (tnew == t) then
                ! write(*,*) 'dt = ', dt, 'l = ', l
                stop "Stepsize underflow in rkqs"
            end if
        end do

        ! l=l+1
        if (errmax > errcon) then
            dtnext = safety * dt * (errmax ** pgrow)
        else
            dtnext = 5.0 * dt !no more than a factor of 5 increase
        end if

        dtdid = dt
        t = t + dt
        if (dt < 0) then
            stop 'Step size underflow in rkqs'
        end if

        ! off next 4 lines by hua 09/09/04
        ! c(1:nr,1:nv) = tmp%c(1:nr,1:nv)
        ! te(1:nr,1:nl) = tmp%te(1:nr,1:nl)
        ! inb(1:nr) = tmp%inb(1:nr)
        ! ipb(1:nr) = tmp%ipb(1:nr)
        ! os=tmp%os
        ! hua then add
        intg=tmp
        ! write(*,*) t
    end subroutine

!---------------------------------------------------------------------------------------------

    !gp 08-jan-10
    !subroutine rkck(nr, meteo, solar, hw, db, hydrau, sys, intg,&
    ! rates, dte, dc, dinb, dipb, t, dt, outintg, err)
    subroutine rkck(nr, meteo, solar, hw, db, hydrau, sys, intg,&
        rates, dte, dc, dinb, dipb, t, dt, outintg, err, statevariables)

        integer(i32), intent(in) :: nr

        !gp 08-jan-10
        character(len=30), intent(in) :: statevariables !to test for 'all except temperature'

        type(meteorology_t) :: meteo
        type(solar_type) solar !solar radiation
        ! type(headwaterdownstream_type), intent(in) :: hdboundary
        type(upstream_boundary_t), intent(in) :: hw !headwater
        type(downstream_boundary_t), intent(in) :: db !downstream boundary
        type(riverhydraulics_type), intent(in) :: hydrau !channel dimensions, hydraulics, physical characters
        type(systemparams) sys
        type(integral_type), intent(in) :: intg
        type(integral_type), intent(out) :: outintg, err !current value and incremented value
        type(rates_t), intent(in) :: rates

        !gp 28-oct-04 real(r64), intent(in):: dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64), intent(in):: dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp
        real(r64), intent(in) :: t, dt !initial time, and step size
        integer(i32) i,j,k
        real(r64), parameter :: a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, &
            b21 = 0.2, b31 = 3./40., b32 = 9./ 40., b41 = 0.3, &
            b42 = -0.9, b43 = 1.2, b51 = -11./54., b52 = 2.5, &
            b53 = -70./27., b54 = 35./27., b61 = 1631./55296., b62 = 175./512., &
            b63 = 575./13824., b64 = 44275./110592., b65 = 253./4096., c1 = 37./378., &
            c3 = 250./621., c4 = 125./594., c6 = 512./1771., dc5 = -277./14336., &
            dc1 =c1 - 2825./ 27648., dc3 = c3 - 18575./48384., dc4 = c4 - 13525./55296., &
            dc6 = c6 - 0.25
        type(integral_type) temp
        !runge-kutta derivatives
        !gp 28-oct-04 real(r64) ak2dc(nr, nv), ak2dte(nr, nl), ak2dinb(nr), ak2dipb(nr)
        !gp real(r64) ak3dc(nr, nv), ak3dte(nr, nl), ak3dinb(nr), ak3dipb(nr)
        !gp real(r64) ak4dc(nr, nv), ak4dte(nr, nl), ak4dinb(nr), ak4dipb(nr)
        !gp real(r64) ak5dc(nr, nv), ak5dte(nr, nl), ak5dinb(nr), ak5dipb(nr)
        !gp real(r64) ak6dc(nr, nv), ak6dte(nr, nl), ak6dinb(nr), ak6dipb(nr)
        real(r64) ak2dc(nr, nv, nl), ak2dte(nr, nl), ak2dinb(nr), ak2dipb(nr)
        real(r64) ak3dc(nr, nv, nl), ak3dte(nr, nl), ak3dinb(nr), ak3dipb(nr)
        real(r64) ak4dc(nr, nv, nl), ak4dte(nr, nl), ak4dinb(nr), ak4dipb(nr)
        real(r64) ak5dc(nr, nv, nl), ak5dte(nr, nl), ak5dinb(nr), ak5dipb(nr)
        real(r64) ak6dc(nr, nv, nl), ak6dte(nr, nl), ak6dinb(nr), ak6dipb(nr) !gp end new block

        !call constructors
        temp= integration_(nr)

        !gp 28-oct-04 temp%c(1:nr,1:nv) = intg%c(1:nr,1:nv) + b21 * dt * dc(1:nr,1:nv)
        temp%c(1:nr,1:nv,1:nl) = intg%c(1:nr,1:nv,1:nl) + b21 * dt * dc(1:nr,1:nv,1:nl) !gp
        temp%te(1:nr,1:nl) = intg%te(1:nr,1:nl) + b21* dt * dte(1:nr,1:nl)
        temp%inb(1:nr) = intg%inb(1:nr) + b21* dt * dinb(1:nr)
        temp%ipb(1:nr) = intg%ipb(1:nr) + b21* dt * dipb(1:nr)

        !second step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak2dte, ak2dc, ak2dinb, ak2dipb, t+ a2*dt)
        do i=1, nr

            !gp do j=1, nv
            !gp temp%c(i,j) = intg%c(i,j) + dt* (b31* dc(i,j) + b32 * ak2dc(i,j))
            !gp if (temp%c(i,j)<0) temp%c(i,j)=1.0e-6
            !gp end do
            !gp do k=1, nl
            !gp temp%te(i,k) = intg%te(i,k) + dt* (b31* dte(i,k) + b32 * ak2dte(i,k))
            !gp end do
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt* (b31* dte(i,j) + b32 * ak2dte(i,j))
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt* (b31* dc(i,k,j) + b32 * ak2dc(i,k,j))
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt* (b31* dinb(i) + b32 * ak2dinb(i))
            temp%ipb(i) = intg%ipb(i) + dt* (b31* dinb(i) + b32 * ak2dipb(i))
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !third step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak3dte, ak3dc, ak3dinb, ak3dipb, t + a3*dt)
        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b41* dc(i, j) + b42 * ak2dc(i, j) + &
            !gp b43 * ak3dc(i, j))
            !gp if (temp%c(i,j)<0) temp%c(i,j)=1.0e-6
            !gp end do
            !gp do k=1, nl
            !gp temp%te(i, k) = intg%te(i, k) + dt * (b41 * dte(i, k) + b42 * ak2dte(i, k) + &
            !gp b43 * ak3dte(i, k))
            !gp end do
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt * (b41 * dte(i,j) + b42 * ak2dte(i,j) + &
                    b43 * ak3dte(i, j))
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b41* dc(i,k,j) + b42 * ak2dc(i,k,j) + &
                        b43 * ak3dc(i,k,j))
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt * (b41 * dinb(i) + b42 * ak2dinb(i) + &
                b43 * ak3dinb(i))
            temp%ipb(i) = intg%ipb(i) + dt * (b41 * dipb(i) + b42 * ak2dipb(i) + &
                b43 * ak3dipb(i))
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !fourth step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak4dte, ak4dc, ak4dinb, ak4dipb, t + a4 * dt)

        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b51 * dc(i, j) + b52 * ak2dc(i, j) + &
            !gp b53 * ak3dc(i, j) + b54 * ak4dc(i, j))
            !gp if (temp%c(i,j)<0) temp%c(i,j)=1.0e-6
            !gp end do
            !gp do k=1, nl
            !gp temp%te(i, k) = intg%te(i, k) + dt * (b51 * dte(i, k) + &
            !gp b52 * ak2dte(i, k) + b53 * ak3dte(i, k) + b54 * ak4dte(i, k))
            !gp end do
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt * (b51 * dte(i,j) + &
                    b52 * ak2dte(i,j) + b53 * ak3dte(i,j) + b54 * ak4dte(i,j))
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b51 * dc(i,k,j) + b52 * ak2dc(i,k,j) + &
                        b53 * ak3dc(i,k,j) + b54 * ak4dc(i,k,j))
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt * (b51 * dinb(i) + &
                b52 * ak2dinb(i) + b53 * ak3dinb(i) + b54 * ak4dinb(i))
            temp%ipb(i) = intg%ipb(i) + dt * (b51 * dipb(i) + &
                b52 * ak2dipb(i) + b53 * ak3dipb(i) + b54 * ak4dipb(i))
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !fifth step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak5dte, ak5dc, ak5dinb, ak5dipb, t + a5 * dt)
        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b61 * dc(i, j) + b62 * ak2dc(i, j) + &
            !gp b63 * ak3dc(i, j) + b64 * ak4dc(i, j) + b65 * ak5dc(i, j))
            !gp if (temp%c(i,j)<0) temp%c(i,j)=1.0e-6
            !gp end do
            !gp do k=1, nl
            !gp temp%te(i, k) = intg%te(i, k) + dt * (b61 * dte(i, k) + b62 * ak2dte(i, k) + &
            !gp b63 * ak3dte(i, k) + b64 * ak4dte(i, k) + b65 * ak5dte(i, k))
            !gp end do
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt * (b61 * dte(i,j) + b62 * ak2dte(i,j) + &
                    b63 * ak3dte(i,j) + b64 * ak4dte(i,j) + b65 * ak5dte(i,j))
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b61 * dc(i,k,j) + b62 * ak2dc(i,k,j) + &
                        b63 * ak3dc(i,k,j) + b64 * ak4dc(i,k,j) + b65 * ak5dc(i,k,j))
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt * (b61 * dinb(i) + b62 * ak2dinb(i) + &
                b63 * ak3dinb(i) + b64 * ak4dinb(i) + b65 * ak5dinb(i))
            temp%ipb(i) = intg%ipb(i) + dt * (b61 * dipb(i) + b62 * ak2dipb(i) + &
                b63 * ak3dipb(i) + b64 * ak4dipb(i) + b65 * ak5dipb(i))
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !sixth step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak6dte, ak6dc, ak6dinb, ak6dipb, t + a6 * dt)

        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            !gp outintg%c(i, j) = intg%c(i, j) + dt * (c1 * dc(i, j) + c3 * ak3dc(i, j) + &
            !gp c4 * ak4dc(i, j) + c6 * ak6dc(i, j))
            !gp if (outintg%c(i,j)<0) temp%c(i,j)=1.0e-6
            !gp end do
            !gp do k=1, nl
            !gp outintg%te(i, k) = intg%te(i, k) + dt * (c1 * dte(i, k) + c3 * ak3dte(i, k) + &
            !gp c4 * ak4dte(i, k) + c6 * ak6dte(i, k))
            !gp end do
            do j=1, nl

                !gp 08-jan-10
                !outintg%te(i,j) = intg%te(i,j) + dt * (c1 * dte(i,j) + c3 * ak3dte(i,j) + &
                ! c4 * ak4dte(i,j) + c6 * ak6dte(i,j))
                if (statevariables /= "All except temperature") then
                    outintg%te(i,j) = intg%te(i,j) + dt * (c1 * dte(i,j) + c3 * ak3dte(i,j) + &
                        c4 * ak4dte(i,j) + c6 * ak6dte(i,j))
                end if

                do k=1, nv
                    outintg%c(i,k,j) = intg%c(i,k,j) + dt * (c1 * dc(i,k,j) + c3 * ak3dc(i,k,j) + &
                        c4 * ak4dc(i,k,j) + c6 * ak6dc(i,k,j))
                    if (outintg%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do

            outintg%inb(i) = intg%inb(i) + dt * (c1 * dinb(i) + c3 * ak3dinb(i) + &
                c4 * ak4dinb(i) + c6 * ak6dinb(i))
            outintg%ipb(i) = intg%ipb(i) + dt * (c1 * dipb(i) + c3 * ak3dipb(i) + &
                c4 * ak4dipb(i) + c6 * ak6dipb(i))
            if (outintg%inb(i)<0) temp%inb(i)=1.0e-6
            if (outintg%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do
! outintg%os=intg%os

        !gp 28-oct-04 err%c(1:nr, 1:nv) = dt * (dc1 * dc(1:nr, 1:nv) + dc3 * ak3dc(1:nr, 1:nv) + &
        !gp dc4 * ak4dc(1:nr, 1:nv) + dc5 * ak5dc(1:nr, 1:nv) + &
        !gp dc6 * ak6dc(1:nr, 1:nv))
        err%c(1:nr, 1:nv, 1:nl) = dt * (dc1 * dc(1:nr, 1:nv, 1:nl) + dc3 * ak3dc(1:nr, 1:nv, 1:nl) + &
            dc4 * ak4dc(1:nr, 1:nv, 1:nl) + dc5 * ak5dc(1:nr, 1:nv, 1:nl) + &
            dc6 * ak6dc(1:nr, 1:nv, 1:nl)) !gp
        err%te(1:nr, 1:nl) = dt * (dc1 * dte(1:nr, 1:nl) + dc3 * ak3dte(1:nr, 1:nl) + &
            dc4 * ak4dte(1:nr, 1:nl) + dc5 * ak5dte(1:nr, 1:nl) + &
            dc6 * ak6dte(1:nr, 1:nl))
        err%inb(1:nr) = dt * (dc1 * dinb(1:nr) + dc3 * ak3dinb(1:nr) + dc4 * ak4dinb(1:nr) + &
            dc5 * ak5dinb(1:nr)+ dc6 * ak6dinb(1:nr))
        err%ipb(1:nr) = dt * (dc1 * dipb(1:nr) + dc3 * ak3dipb(1:nr) + dc4 * ak4dipb(1:nr) + &
            dc5 * ak5dipb(1:nr) + dc6 * ak6dipb(1:nr))

    end subroutine

!-------------------------------------------------------------------------------------------
! runge-kutta 4th order integration

    subroutine rk4(nr, meteo, solar, hw, db, hydrau, sys, intg,&
        rates, dte, dc, dinb, dipb, t, dt, statevariables)

        integer(i32), intent(in) :: nr
        type(meteorology_t) :: meteo
        type(solar_type) solar
        type(upstream_boundary_t), intent(in) :: hw
        type(downstream_boundary_t), intent(in) :: db
        type(riverhydraulics_type), intent(in) :: hydrau
        type(systemparams) sys
        type(integral_type), intent(inout) :: intg
        type(rates_t), intent(in) :: rates

        !gp 28-oct-04 real(r64), intent(out):: dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64), intent(out):: dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp

        real(r64), intent(in) :: t, dt !initial time, and step size

        !gp 08-jan-10
        character(len=30), intent(in) :: statevariables !to test for 'all except temperature'

        integer(i32) i,j,k

        type(integral_type) temp
        !runge-kutta derivatives

        real(r64) ak2dc(nr, nv, nl), ak2dte(nr, nl), ak2dinb(nr), ak2dipb(nr)
        real(r64) ak3dc(nr, nv, nl), ak3dte(nr, nl), ak3dinb(nr), ak3dipb(nr)
        real(r64) ak4dc(nr, nv, nl), ak4dte(nr, nl), ak4dinb(nr), ak4dipb(nr) !gp end new block

        !call constructors
        temp= integration_(nr)

        !first step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, intg%te, intg%c, intg%inb, intg%ipb, &
            rates, dte, dc, dinb, dipb, t)
        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt/2 * dte(i,j)
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt/2 * dc(i,k,j)
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt/2 * dinb(i)
            temp%ipb(i) = intg%ipb(i) + dt/2 * dipb(i)
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !second step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak2dte, ak2dc, ak2dinb, ak2dipb, t+ dt/2)

        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt/2 * ak2dte(i,j)
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt/2 * ak2dc(i,k,j)
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt/2 * ak2dinb(i)
            temp%ipb(i) = intg%ipb(i) + dt/2 * ak2dipb(i)
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !third step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak3dte, ak3dc, ak3dinb, ak3dipb, t + dt/2)

        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            do j=1, nl
                temp%te(i,j) = intg%te(i,j) + dt * ak3dte(i,j)
                do k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * ak3dc(i,k,j)
                    if (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            temp%inb(i) = intg%inb(i) + dt * ak3dinb(i)
            temp%ipb(i) = intg%ipb(i) + dt * ak3dipb(i)
            if (temp%inb(i)<0) temp%inb(i)=1.0e-6
            if (temp%ipb(i)<0) temp%ipb(i)=1.0e-6
        end do

        !fourth step
        call derivs(nr, meteo, solar, hw, db, hydrau, sys, temp%te, temp%c, temp%inb, temp%ipb, &
            rates, ak4dte, ak4dc, ak4dinb, ak4dipb, t + dt)

        do i=1, nr

            !gp 28-oct-04 do j=1, nv
            do j=1, nl

                !gp 08-jan-10
                if (statevariables /= "All except temperature") then
                    intg%te(i,j) = intg%te(i,j) + dt/6 * (dte(i,j) + &
                        2 * ak2dte(i,j) + 2 * ak3dte(i,j) + ak4dte(i,j))

                    !gp 12-jan-10
                else
                    if (hydrau%reach(i)%te_ini < 0) then
                        !use diel headwater as diel temp for each reach during intg
                        !call instanteousheadwater(hw, t, intg%te(0,1), intg%c(0,:,1), ph)
                        intg%te(i, j) = intg%te(0, 1)
                    else
                        !use initial temp for each reach as const temp during intg
                        intg%te(i, j) = hydrau%reach(i)%te_ini
                    end if

                end if

                do k=1, nv
                    intg%c(i,k,j) = intg%c(i,k,j) + dt/6 * (dc(i,k,j) + &
                        2 * ak2dc(i,k,j) + 2 * ak3dc(i,k,j) + ak4dc(i,k,j))
                    if (intg%c(i,k,j)<0) intg%c(i,k,j)=1.0e-6
                end do
            end do !gp end new block

            intg%inb(i) = intg%inb(i) + dt/6 * (dinb(i) + &
                2 * ak2dinb(i) + 2 * ak3dinb(i) + ak4dinb(i))
            intg%ipb(i) = intg%ipb(i) + dt/6 * (dipb(i) + &
                2 * ak2dipb(i) + 2 * ak3dipb(i) + ak4dipb(i))
            if (intg%inb(i)<0) intg%inb(i)=1.0e-6
            if (intg%ipb(i)<0) intg%ipb(i)=1.0e-6
        end do
        !outintg%os=intg%os
    end subroutine rk4


    pure function ini_or_prev(ini, prev)
        real(r64), intent(in) :: ini, prev
        real(r64) :: ini_or_prev

        if (ini < 0) then
            ini_or_prev = prev
        else
            ini_or_prev = ini
        end if
    end function ini_or_prev

    subroutine print_report_integration_method(imeth)
        character(len=30), intent(in) :: imeth ! integration method
        select case (imeth)
          case ('Euler')
            write(*,*) 'Integrating: Euler method.'
          case ('Runge-Kutta')
            write(*,*) 'Integrating: Runge-Kutta 4th Order method.'
          case ('Adaptive step')
            write(*,*) 'Integrating: Adaptive-step method.'
          case default
            write(8,*) '** Error: integration method not recognized **' !gp 20-oct-04
            stop 'Error: integration method not recognized'
        end select
    end subroutine print_report_integration_method

end module class_integration
