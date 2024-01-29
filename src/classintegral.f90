!classIntegration

MODULE Class_Integration
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_derivs, only: derivs
    use class_hydraulics, only: riverhydraulics_type
    USE Class_IntegrationData, only: Integral_type, integration_, &
        saveHeatFluxTribs, saveheatfluxadvecdisp, saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, &
        saveheatfluxconv, saveheatfluxevap, saveheatfluxjsed, saveheatfluxjhyporheic, &
        os, phitsave, savedofluxheadwater, saveco2fluxheadwater, savedofluxtribs, savedofluxadvecdisp, &
        saveco2fluxtribs, saveco2fluxadvecdisp, phinsave, phipsave, phicsave, philsave, phitotalsave , &
        savebotalgphoto, savebotalgresp, savebotalgdeath, savebotalgnetgrowth, &
        sodpr, jnh4pr, jno3pr, jch4pr, jsrppr, csodpr, &
        diagfluxdo, diagfluxcbod, diagfluxnh4, diagfluxno3, diagfluxsrp, diagfluxic, savedofluxreaer, &
        savedofluxcbodfast, savedofluxcbodslow, savedofluxnitrif,  savedofluxphytoresp, &
        savedofluxphytophoto, savedofluxbotalgresp, savedofluxbotalgphoto, savedofluxsod, savedofluxcod, &
        saveco2fluxreaer, saveco2fluxcbodfast, saveco2fluxcbodslow, saveco2fluxphytoresp, saveco2fluxphytophoto, &
        saveco2fluxbotalgresp, saveco2fluxbotalgphoto, saveco2fluxsod, hypofluxdo, hypofluxcbod, hypofluxnh4, &
        hypofluxno3, hypofluxsrp, hypofluxic
    USE Class_LightHeat, only: lightheat, heatbalance, lightextinction
    USE Class_Phsolve, ONLY: ChemRates, ModFP2, phsolNewton, pHsolBisect, phsolbrent, ct
    USE Class_SolarCalc, only: solar_type, solarcalc
    USE Class_SourceIn, only: load, sourcescalc
    USE Class_SystemParams, only: SystemParams
    USE m_downstream_boundary, only: downstream_boundary_t, instanteousdownstreamboundary
    USE m_meteorology, only: meteorology_t, instanteousmeteo
    USE m_output, ONLY: outdata_t, output
    use m_oxygen, only: oxygen_inhibition_and_enhancement, oxygen_saturation
    USE m_rates, only: rates_t
    use m_tempadjust, only: temp_adjust
    USE m_upstream_boundary, only: upstream_boundary_t, instanteousheadwater
    USE nrtype, only: nv, nl, adam, bdam, cpw, rhow

    IMPLICIT NONE

CONTAINS
!----------------------------------------------------------------------------------
    SUBROUTINE Integration(sys, Rates, Meteo, Solar, HW, DB, hydrau, pr, nr)


        TYPE(Integral_type) intg !integral data structure
        TYPE(SystemParams) sys
        TYPE(rates_t), INTENT(IN) :: Rates
        TYPE(meteorology_t) Meteo
        TYPE(solar_type) Solar !solar radiation
        TYPE(upstream_boundary_t) HW !headwater
        TYPE(downstream_boundary_t) DB !downstream boundary
        TYPE(RiverHydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
        TYPE(outdata_t), INTENT(OUT) :: pr !print out data structure
        INTEGER(i32), INTENT(IN) :: nr !number of reach
        REAL(r64) t
        !Heat related values

        INTEGER(i32) ip, ic, i, j, k
        REAL(r64) pH
        !gp 27-Oct-04 REAL(r64) dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64) dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp add dim to dc for nl

        !gp 24-Jun-09
        REAL(r64) TSS, TP, TN, DOSat, NH3
        INTEGER(i32) iskip, nskip
        nskip = sys%nc / 32 !for output of dynamic calculations every 45 minutes (32 times per day)
        iskip = 0 !initialize counter for skipping dynamic output

        !gp 24-Jun-09
        !CHARACTER(LEN=30), INTENT(IN) :: IMethpH ! pH solve method

        intg=Integration_(nr)

        !Set Initial Conditions
        dc=0; dINb=0; dIPb=0; dTe=0

        !initial time
        t=0

        !initial temperature and concentrations
        !CALL InstanteousHeadwater(HDboundary, sys%tday, 1, intg%Te(0,1), intg%c(0,:), pH)
        !gp CALL InstanteousHeadwater(HW, t, intg%Te(0,1), intg%c(0,:), pH)
        CALL InstanteousHeadwater(HW, t, intg%Te(0,1), intg%c(0,:,1), pH) !gp

        DO i = 1, nr

            !gp 21-Nov-06
            !DO j = 1, nl
            ! intg%Te(i, j) = intg%Te(0, 1)
            !!gp 27-Oct-04 END DO
            ! DO k = 1, nv - 1
            ! intg%c(i, k, j) = intg%c(0, k, 1)
            ! END DO
            !END DO !gp
            !!gp 11-Jan-05 IF (Rates%typeF == "First-order") intg%c(i, nv) = 1.0
            !intg%c(i, nv, 1) = hydrau%reach(i)%botalg0 * Rates%adc * Rates%aca !'convert mgA/m^2 to gD/m^2
            !IF (Rates%typeF == "First-order" .and. hydrau%reach(i)%botalg0==0) intg%c(i, nv, 1) = 1.0
            DO j = 1, nl
                If (hydrau%reach(i)%Te_ini < 0) Then
                    intg%Te(i, j) = intg%Te(0, 1)
                Else
                    intg%Te(i, j) = hydrau%reach(i)%Te_ini
                End If
                If (hydrau%reach(i)%c01_ini < 0) Then
                    intg%c(i, 1, j) = intg%c(0, 1, 1)
                Else
                    intg%c(i, 1, j) = hydrau%reach(i)%c01_ini
                End If
                If (hydrau%reach(i)%c02_ini < 0) Then
                    intg%c(i, 2, j) = intg%c(0, 2, 1)
                Else
                    intg%c(i, 2, j) = hydrau%reach(i)%c02_ini
                End If
                If (hydrau%reach(i)%c03_ini < 0) Then
                    intg%c(i, 3, j) = intg%c(0, 3, 1)
                Else
                    intg%c(i, 3, j) = hydrau%reach(i)%c03_ini
                End If
                If (hydrau%reach(i)%c04_ini < 0) Then
                    intg%c(i, 4, j) = intg%c(0, 4, 1)
                Else
                    intg%c(i, 4, j) = hydrau%reach(i)%c04_ini
                End If
                If (hydrau%reach(i)%c05_ini < 0) Then
                    intg%c(i, 5, j) = intg%c(0, 5, 1)
                Else
                    intg%c(i, 5, j) = hydrau%reach(i)%c05_ini
                End If
                If (hydrau%reach(i)%c06_ini < 0) Then
                    intg%c(i, 6, j) = intg%c(0, 6, 1)
                Else
                    intg%c(i, 6, j) = hydrau%reach(i)%c06_ini
                End If
                If (hydrau%reach(i)%c07_ini < 0) Then
                    intg%c(i, 7, j) = intg%c(0, 7, 1)
                Else
                    intg%c(i, 7, j) = hydrau%reach(i)%c07_ini
                End If
                If (hydrau%reach(i)%c08_ini < 0) Then
                    intg%c(i, 8, j) = intg%c(0, 8, 1)
                Else
                    intg%c(i, 8, j) = hydrau%reach(i)%c08_ini
                End If
                If (hydrau%reach(i)%c09_ini < 0) Then
                    intg%c(i, 9, j) = intg%c(0, 9, 1)
                Else
                    intg%c(i, 9, j) = hydrau%reach(i)%c09_ini
                End If
                If (hydrau%reach(i)%c10_ini < 0) Then
                    intg%c(i, 10, j) = intg%c(0, 10, 1)
                Else
                    intg%c(i, 10, j) = hydrau%reach(i)%c10_ini
                End If
                If (hydrau%reach(i)%c11_ini < 0) Then
                    intg%c(i, 11, j) = intg%c(0, 11, 1)
                Else
                    intg%c(i, 11, j) = hydrau%reach(i)%c11_ini
                End If
                If (hydrau%reach(i)%c12_ini < 0) Then
                    intg%c(i, 12, j) = intg%c(0, 12, 1)
                Else
                    intg%c(i, 12, j) = hydrau%reach(i)%c12_ini
                End If
                If (hydrau%reach(i)%c13_ini < 0) Then
                    intg%c(i, 13, j) = intg%c(0, 13, 1)
                Else
                    intg%c(i, 13, j) = hydrau%reach(i)%c13_ini
                End If
                If (hydrau%reach(i)%c14_ini < 0) Then
                    intg%c(i, 14, j) = intg%c(0, 14, 1)
                Else
                    intg%c(i, 14, j) = hydrau%reach(i)%c14_ini
                End If
                If (hydrau%reach(i)%c15_ini < 0) Then
                    intg%c(i, nv-2, j) = intg%c(0, nv-2, 1)
                Else
                    intg%c(i, nv-2, j) = hydrau%reach(i)%c15_ini
                End If
                If (hydrau%reach(i)%pH_ini < 0) Then
                    intg%c(i, nv-1, j) = intg%c(0, nv-1, 1)
                Else
                    intg%c(i, nv-1, j) = cT(hydrau%reach(i)%pH_ini, intg%c(i, nv-2, j), intg%Te(i, j), intg%c(i, 1, j))
                End If
            END DO
            If (hydrau%reach(i)%c17_ini < 0) Then !bottom algae gD/m^2
                intg%c(i, nv, 1) = 0
            Else
                intg%c(i, nv, 1) = hydrau%reach(i)%c17_ini
            End If
            IF (Rates%typeF == "First-order" .and. intg%c(i, nv, 1) < 1) intg%c(i, nv, 1) = 1.0

            !gp 02-Nov-07
            !If (hydrau%reach(i)%NINb_ini >= 0) Then !intracellular N mgN/gD
            !intg%INb(i) = hydrau%reach(i)%NINb_ini * intg%c(i, nv, j)
            If (hydrau%reach(i)%NINb_ini >= 0) Then !intracellular N mgN/gD
                intg%INb(i) = hydrau%reach(i)%NINb_ini * intg%c(i, nv, 1)

            End If

            !gp 02-Nov-07
            !If (hydrau%reach(i)%NIPb_ini >= 0) Then !intracellular P mgP/gD
            !intg%IPb(i) = hydrau%reach(i)%NIPb_ini * intg%c(i, nv, j)
            If (hydrau%reach(i)%NIPb_ini >= 0) Then !intracellular P mgP/gD
                intg%IPb(i) = hydrau%reach(i)%NIPb_ini * intg%c(i, nv, 1)

            End If

            IF (sys%simHyporheicWQ == "Level 2" .and. Rates%typeH == "First-order") intg%c(i, nv, 2) = 1.0 !gp 15-Nov-04
        END DO

        SELECT CASE (sys%IMeth)
          CASE ('Euler')
            WRITE(*,*) 'Integrating: Euler method.'
          CASE ('Runge-Kutta')
            WRITE(*,*) 'Integrating: Runge-Kutta 4th Order method.'
          CASE ('Adaptive step')
            WRITE(*,*) 'Integrating: Adaptive-step method.'
          CASE DEFAULT
            WRITE(8,*) '** Error: integration method not recognized **' !gp 20-Oct-04
            STOP 'Error: integration method not recognized'
        END SELECT

        !Integration
        DO ip = 1, sys%np !np


            IF (ip >= sys%np) THEN

                !gp 17-Feb-05
                !CALL Save_init_step(nr, intg, pr, sys%dt, sys%IMethpH, Rates)
                CALL Save_init_step(nr, intg, pr, sys%dt, sys%IMethpH, sys%showDielResults, Rates)

            END IF
            !Euler method
            IF (sys%IMeth == "Euler") THEN
                WRITE(*,*) 'day', ip
                DO ic = 1, sys%nc !Time step in each day
                    CALL derivs(nr, Meteo, Solar,HW, DB, hydrau, sys, intg%Te, intg%c, intg%INb, &
                        intg%IPb, Rates, dTe, dc, dINb, dIPb, t)
                    !finite method: new= old + old * derive
                    DO i = 1, nr
                        DO j = 1, nl

                            !gp 08-Jan-10
                            !intg%Te(i, j) = intg%Te(i, j) + dTe(i, j) * sys%dt
                            If (sys%stateVariables /= "All except temperature") Then
                                intg%Te(i, j) = intg%Te(i, j) + dTe(i, j) * sys%dt

                                !gp 12-Jan-10
                            Else
                                If (hydrau%reach(i)%Te_ini < 0) Then
                                    !use diel headwater as diel temp for each reach during intg
                                    !CALL InstanteousHeadwater(HW, t, intg%Te(0,1), intg%c(0,:,1), pH)
                                    intg%Te(i, j) = intg%Te(0, 1)
                                Else
                                    !use initial temp for each reach as const temp during intg
                                    intg%Te(i, j) = hydrau%reach(i)%Te_ini
                                End If

                            End If

                            !gp 27-Oct-04 END DO
                            ! write (8, *) i, intg%c(i,3), dc(i,3)
                            ! IF (i==1) THEN
                            ! WRITE(8, '(I3, 16F26.18)') i, (intg%c(i, k), k=1, nv)
                            ! WRITE(8, '(I3, 16F26.18)') i, (dc(i, k), k=1, nv)
                            ! WRITE(8, '(I3, 4F26.18)') i, intg%INb(i), dINb(i), intg%IPb(i), dIPb(i)
                            ! END IF
                            DO k = 1, nv
                                intg%c(i, k, j) = intg%c(i, k, j) + dc(i, k, j) * sys%dt
                                IF (intg%c(i, k, j) < 0) intg%c(i, k, j) = 1.0E-6_r64
                            END DO
                        END DO !gp 27-Oct-04
                        intg%INb(i) = intg%INb(i) + dINb(i) * sys%dt
                        intg%IPb(i) = intg%IPb(i) + dIPb(i) * sys%dt
                        IF (intg%INb(i) < 0) intg%INb(i) = 1.0E-6_r64
                        IF (intg%IPb(i) < 0) intg%IPb(i) = 1.0E-6_r64
                    END DO


                    !gp 24-Jun-09 write dynamic output every 45 minutes
                    IF (sys%writeDynamic == "Yes") THEN
                        if (iskip == 0) then
                            DO i=0, nr
                                j = 1
                                TSS = intg%c(i, 11, j) * Rates%ada + intg%c(i, 2, j) + &
                                    intg%c(i, 12, j)
                                TP = intg%c(i, 11, j) * Rates%apa + intg%c(i, 9, j) + &
                                    intg%c(i, 10, j)
                                TN = intg%c(i, 11, j) * Rates%ana + intg%c(i, 6, j) + &
                                    intg%c(i, 7, j) + intg%c(i, 8, j)
                                DOSat = oxygen_saturation(intg%Te(i, j), hydrau%reach(i)%elev)

                                !gp 23-Nov-09
                                !IF (sys%IMethpH == "Newton-Raphson") THEN
                                !CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !ELSE
                                !CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !END IF
                                IF (sys%IMethpH == "Newton-Raphson") THEN
                                    CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                ELSEIF (sys%IMethpH == "Bisection") THEN
                                    CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                ELSE
                                    CALL pHsolBrent(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                END IF

                                NH3 = 1.0_r64/(1 + 10.0_r64 ** (-pH)/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                                    (intg%Te(i, j) + 273.15_r64))) * intg%c(i, 7, j)
                                WRITE(12, '(I13, 41F13.4)') i, t, intg%Te(i, j), &
                                    (intg%c(i, k, j), k=1, nv-2), pH, &
                                    intg%c(i, nv, j), intg%c(i, nv, j)* Rates%mgA / Rates%mgD * 1000, &
                                    TSS, TP, TN , DOSat, NH3, &
                                    intg%INb(i), &
                                    intg%INb(i), &
                                    phitSave(i), philSave(i), &
                                    phinSave(i), phipSave(i), &
                                    phicSave(i), phitotalSave(i), &
                                    saveBotAlgPhoto(i), saveBotAlgResp(i), saveBotAlgDeath(i), saveBotAlgNetGrowth(i), &
                                    saveBotAlgPhoto(i)/Rates%ada, saveBotAlgResp(i)/Rates%ada, &
                                    saveBotAlgDeath(i)/Rates%ada, saveBotAlgNetGrowth(i)/Rates%ada
                            END DO
                        end if
                        iskip = iskip + 1
                        if (iskip == nskip) then
                            iskip = 0
                        end if
                    END IF

                    t = t + sys%dt

                    !save intemediate steps ?
                    IF (ip >= sys%np) THEN

                        !gp 17-Feb-05
                        !CALL Save_a_step(nr, intg,pr, sys%dt, sys%IMethpH, Rates)
                        CALL Save_a_step(nr, intg,pr, sys%dt, sys%IMethpH, sys%showDielResults, Rates)

                    END IF
                END DO
            ELSE IF (sys%IMeth == 'Runge-Kutta') THEN
                WRITE(*,*) 'day', ip
                DO ic = 1, sys%nc !Time step in each day

                    !gp 08-Jan-10
                    !CALL rk4(nr, Meteo, Solar,HW, DB, hydrau, sys, intg, &
                    ! Rates, dTe, dc, dINb, dIPb, t, sys%dt)
                    CALL rk4(nr, Meteo, Solar,HW, DB, hydrau, sys, intg, &
                        Rates, dTe, dc, dINb, dIPb, t, sys%dt, sys%stateVariables)

                    !gp 24-Jun-09 write dynamic output every 45 minutes
                    IF (sys%writeDynamic == "Yes") THEN
                        if (iskip == 0) then
                            DO i=0, nr
                                j = 1
                                TSS = intg%c(i, 11, j) * Rates%ada + intg%c(i, 2, j) + &
                                    intg%c(i, 12, j)
                                TP = intg%c(i, 11, j) * Rates%apa + intg%c(i, 9, j) + &
                                    intg%c(i, 10, j)
                                TN = intg%c(i, 11, j) * Rates%ana + intg%c(i, 6, j) + &
                                    intg%c(i, 7, j) + intg%c(i, 8, j)
                                DOSat = oxygen_saturation(intg%Te(i, j), hydrau%reach(i)%elev)

                                !gp 23-Nov-09
                                !IF (sys%IMethpH == "Newton-Raphson") THEN
                                !CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !ELSE
                                !CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                ! intg%c(i, 1, j))
                                !END IF
                                IF (sys%IMethpH == "Newton-Raphson") THEN
                                    CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                ELSEIF (sys%IMethpH == "Bisection") THEN
                                    CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                ELSE
                                    CALL pHsolBrent(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                                        intg%c(i, 1, j))
                                END IF

                                NH3 = 1.0_r64/(1 + 10.0_r64 ** (-pH)/10.0_r64** -(0.09018_r64 + 2729.92_r64 / &
                                    (intg%Te(i, j) + 273.15_r64))) * intg%c(i, 7, j)
                                WRITE(12, '(I13, 41F13.4)') i, t, intg%Te(i, j), &
                                    (intg%c(i, k, j), k=1, nv-2), pH, &
                                    intg%c(i, nv, j), intg%c(i, nv, j)* Rates%mgA / Rates%mgD * 1000, &
                                    TSS, TP, TN , DOSat, NH3, &
                                    intg%INb(i), &
                                    intg%INb(i), &
                                    phitSave(i), philSave(i), &
                                    phinSave(i), phipSave(i), &
                                    phicSave(i), phitotalSave(i), &
                                    saveBotAlgPhoto(i), saveBotAlgResp(i), saveBotAlgDeath(i), saveBotAlgNetGrowth(i), &
                                    saveBotAlgPhoto(i)/Rates%ada, saveBotAlgResp(i)/Rates%ada, &
                                    saveBotAlgDeath(i)/Rates%ada, saveBotAlgNetGrowth(i)/Rates%ada
                            END DO
                        end if
                        iskip = iskip + 1
                        if (iskip == nskip) then
                            iskip = 0
                        end if
                    END IF

                    t = t + sys%dt

                    !save intemediate steps ?
                    IF (ip >= sys%np) THEN

                        !gp 17-Feb-05
                        !CALL Save_a_step(nr, intg,pr, sys%dt, sys%IMethpH, Rates)
                        CALL Save_a_step(nr, intg,pr, sys%dt, sys%IMethpH, sys%showDielResults, Rates)

                    END IF
                END DO
            ELSE IF (sys%IMeth == "Adaptive step") THEN
                WRITE(*,*) 'day', ip
                IF (ip >= sys%np) THEN

                    !gp 08-Jan-10
                    !CALL odeint(intg, sys, Rates, Meteo, Solar, &
                    ! HW, DB, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .TRUE.)
                    CALL odeint(intg, sys, Rates, Meteo, Solar, &
                        HW, DB, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .TRUE., sys%stateVariables)

                ELSE

                    !gp 08-Jan-10
                    !CALL odeint(intg, sys, Rates, Meteo, Solar, &
                    ! HW, DB, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .FALSE.)
                    CALL odeint(intg, sys, Rates, Meteo, Solar, &
                        HW, DB, hydrau, pr, nr, 0.0_r64, 1.0_r64, sys%dt, .FALSE., sys%stateVariables)

                END IF
            ELSE
                STOP 'Please specify integration method'
            END IF
        END DO

    END SUBROUTINE

!-----------------------------------------------------------------------------------------
    !Save the initial step of the last simulation day
    !called only once

!gp 17-Feb-05
!SUBROUTINE Save_init_step(nr, intg, pr, dt, IMethpH, Rates)
    SUBROUTINE Save_init_step(nr, intg, pr, dt, IMethpH, showDielResults, Rates)

! USE Class_Output, ONLY: outdata_type

        INTEGER(i32), INTENT(IN) :: nr

        REAL(r64), INTENT(IN) :: dt
        CHARACTER(LEN=30), INTENT(IN) :: IMethpH ! pH solve method

        !gp 17-Feb-05
        CHARACTER(LEN=30), INTENT(IN) :: showDielResults !Yes or No (only used in Excel VBA)

        TYPE(rates_t), INTENT(IN) :: Rates
        TYPE(Integral_type), INTENT(IN) :: intg
        TYPE(outdata_t), INTENT(OUT) :: pr

        INTEGER(i32) i, j, k
        !gp 27-Oct-04 REAL(r64) :: pHss(0:nr), Kamm, pH
        REAL(r64) :: pHss(0:nr, 0:nl), Kamm, pH !gp add dim for nl


! !gp debug
! OPEN (unit=9, FILE='debug.out', status='REPLACE', ACTION='WRITE')
! WRITE(9,*) 'Program is now just before TYPE(RiverHydraulics_type) hydrau'
! CLOSE (9)

        !gp 29-Oct-04
        TYPE(RiverHydraulics_type) hydrau !channel dimensions, hydraulics, physical characters

! !gp debug
! OPEN (unit=9, FILE='debug.out', status='REPLACE', ACTION='WRITE')
! WRITE(9,'(A64)') 'Program is now just after TYPE(RiverHydraulics_type) hydrau'
! CLOSE (9)

        pr%tdy(0)=0
        DO i = 0, nr !nr public
            !gp 27-Oct-04 add nl dimension
            !gp pr%Temn(i) = intg%Te(i, 1)
            !gp pr%Temx(i) = intg%Te(i, 1)
            !gp pr%Teav(i) = 0
            !gp pr%osav(i) = 0
            !gp pr%pHsav(i) = 0
            !gp DO k = 1, nv
            !gp pr%cmn(i, k) = intg%c(i, k)
            !gp pr%cmx(i, k) = intg%c(i, k)
            !gp END DO
            !gp IF (IMethpH == "Newton-Raphson") THEN !public IMethpH As String
            !gp CALL phsolNewton(pH, intg%c(i, 15), intg%Te(i, 1), intg%c(i, 14), &
            !gp intg%c(i, 1))
            !gp ELSE
            !gp CALL pHsolBisect(pH, intg%c(i, 15), intg%Te(i, 1), intg%c(i, 14), &
            !gp intg%c(i, 1))
            !gp END IF
            !gp pr%pHav(i) = 0
            !gp pr%pHmn(i) = pH
            !gp pr%pHmx(i) = pH
            !gp pHss(i) = pH
            !gp Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%Te(i, 1) + 273.15_r64)))
            !gp pr%NH3av(i) = 0
            !gp pr%NH3mn(i) = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-pH) / Kamm) * intg%c(i, 7)
            !gp pr%NH3mx(i) = pr%NH3mn(i)
            !gp pr%TPav(i) = 0
            !gp pr%TPmn(i) = intg%c(i, 9) + intg%c(i, 10) + intg%c(i, 11) * Rates%apa
            !gp pr%TPmx(i) = pr%TPmn(i)
            !gp pr%TNav(i) = 0
            !gp pr%TNmn(i) = intg%c(i, 6) + intg%c(i, 7) + intg%c(i, 8) + &
            !gp intg%c(i, 11) * Rates%ana
            !gp pr%TNmx(i) = pr%TNmn(i)
            DO j = 1, nl
                pr%Temn(i, j) = intg%Te(i, j)
                pr%Temx(i, j) = intg%Te(i, j)
                pr%Teav(i, j) = 0
                pr%osav(i, j) = 0
                pr%pHsav(i, j) = 0
                DO k = 1, nv
                    pr%cmn(i, k, j) = intg%c(i, k, j)
                    pr%cmx(i, k, j) = intg%c(i, k, j)
                END DO

                !gp 23-Nov-09
                !IF (IMethpH == "Newton-Raphson") THEN !public IMethpH As String
                ! CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !ELSE
                !CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !END IF
                IF (IMethpH == "Newton-Raphson") THEN
                    CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                ELSEIF (IMethpH == "Bisection") THEN
                    CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                ELSE
                    CALL pHsolBrent(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                END IF

                pr%pHav(i, j) = 0
                pr%pHmn(i, j) = pH
                pr%pHmx(i, j) = pH
                pHss(i, j) = pH
                Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%Te(i, j) + 273.15_r64)))
                pr%NH3av(i, j) = 0
                pr%NH3mn(i, j) = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-pH) / Kamm) * intg%c(i, 7, j)
                pr%NH3mx(i, j) = pr%NH3mn(i, j)
                pr%TPav(i, j) = 0
                pr%TPmn(i, j) = intg%c(i, 9, j) + intg%c(i, 10, j) + intg%c(i, 11, j) * Rates%apa
                pr%TPmx(i, j) = pr%TPmn(i, j)
                pr%TNav(i, j) = 0
                pr%TNmn(i, j) = intg%c(i, 6, j) + intg%c(i, 7, j) + intg%c(i, 8, j) + &
                    intg%c(i, 11, j) * Rates%ana
                pr%TNmx(i, j) = pr%TNmn(i, j) !gp end of new block of code 27-Oct-04
            END DO
            !'gp 15-Nov-04 reach average daily average sediment fluxes
            pr%HypoFluxDOav(i) = 0; pr%HypoFluxCBODav(i) = 0; pr%HypoFluxNH4av(i) = 0
            pr%HypoFluxNO3av(i) = 0; pr%HypoFluxSRPav(i) = 0
            pr%DiagFluxDOav(i) = 0; pr%DiagFluxCBODav(i) = 0; pr%DiagFluxNH4av(i) = 0
            pr%DiagFluxNO3av(i) = 0; pr%DiagFluxSRPav(i) = 0
            !'gp 11-Jan-05 min/max/mean cell quotas mgN/gD and mgP/gD
            If (intg%c(i, nv, 1) > 0) Then
                pr%NINbav(i) = 0
                pr%NINbmn(i) = intg%INb(i) / intg%c(i, nv, 1)
                pr%NINbmx(i) = intg%INb(i) / intg%c(i, nv, 1)
                pr%NIPbav(i) = 0
                pr%NIPbmn(i) = intg%IPb(i) / intg%c(i, nv, 1)
                pr%NIPbmx(i) = intg%IPb(i) / intg%c(i, nv, 1)
            Else
                pr%NINbav(i) = 0; pr%NINbmn(i) = 0; pr%NINbmx(i) = 0
                pr%NIPbav(i) = 0; pr%NIPbmn(i) = 0; pr%NIPbmx(i) = 0
            End If

            !gp 25-Jun-09
            pr%av_BotAlgPhoto(i)=0; pr%av_BotAlgResp(i)=0; pr%av_BotAlgDeath(i)=0; pr%av_BotAlgNetGrowth(i)=0

        END DO

        !gp 17-Feb-05
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        IF (showDielResults == "Yes") THEN
            !x x x x x x x x x x x x x x x x x x x

            pr%nj = 0
            DO i = 0, nr
                DO j = 1, nl
                    pr%Tepr(i, pr%nj, j) = intg%Te(i, j)
                    !gp 27-Oct-04 END DO
                    !gp DO k = 1, nv
                    !gp pr%cpr(i, k, pr%nj) = intg%c(i, k)
                    !gp END DO
                    !gp pr%pHpr(i, pr%nj) = pHss(i)
                    DO k = 1, nv
                        pr%cpr(i, k, pr%nj, j) = intg%c(i, k, j)
                    END DO
                    pr%pHpr(i, pr%nj, j) = pHss(i, j)
                END DO !gp end of new block of code 27-Oct-04
            END DO
            DO i = 1, nr
                !gp 02-Nov-04 pr%NINbpr(i, pr%nj) = intg%INb(i) / intg%c(i, nv) !?
                !gp pr%NIPbpr(i, pr%nj) = intg%IPb(i) / intg%c(i, nv) !?
                pr%NINbpr(i, pr%nj) = intg%INb(i) / intg%c(i, nv, 1)
                pr%NIPbpr(i, pr%nj) = intg%IPb(i) / intg%c(i, nv, 1) !gp 02-Nov-04 end new block
                !gp 20-Oct-04 output growth limitation factors for bottom algae
                pr%phitotalSavepr(i, pr%nj) = phitotalSave(i)
                pr%phitSavepr(i, pr%nj) = phitSave(i)
                pr%philSavepr(i, pr%nj) = philSave(i)
                pr%phinSavepr(i, pr%nj) = phinSave(i)
                pr%phipSavepr(i, pr%nj) = phipSave(i)
                pr%phicSavepr(i, pr%nj) = phicSave(i)
                !gp 28-Oct-04 output hyporheic fluxes of constituents (total surface area)
                pr%HypoFluxDOpr(i, pr%nj) = HypoFluxDO(i) !'dissolved oxygen gO2/m^2/d
                pr%HypoFluxCBODpr(i, pr%nj) = HypoFluxCBOD(i) !'fast CBOD gO2/m^2/d
                pr%HypoFluxNH4pr(i, pr%nj) = HypoFluxNH4(i) !'ammonia mgN/m^2/d
                pr%HypoFluxNO3pr(i, pr%nj) = HypoFluxNO3(i) !'nitrate mgN/m^2/d
                pr%HypoFluxSRPpr(i, pr%nj) = HypoFluxSRP(i) !'SRP mgP/m^2/d
                pr%HypoFluxICpr(i, pr%nj) = HypoFluxIC(i) !'inorganic C gC/m^2/d
                !gp 01-Nov-04 output diagenesis fluxes of constituents (total surface area)
                pr%DiagFluxDOpr(i, pr%nj) = DiagFluxDO(i) !'dissolved oxygen gO2/m^2/d
                pr%DiagFluxCBODpr(i, pr%nj) = DiagFluxCBOD(i) !'fast CBOD gO2/m^2/d
                pr%DiagFluxNH4pr(i, pr%nj) = DiagFluxNH4(i) !'ammonia mgN/m^2/d
                pr%DiagFluxNO3pr(i, pr%nj) = DiagFluxNO3(i) !'nitrate mgN/m^2/d
                pr%DiagFluxSRPpr(i, pr%nj) = DiagFluxSRP(i) !'SRP mgP/m^2/d
                pr%DiagFluxICpr(i, pr%nj) = DiagFluxIC(i) !'inorganic C gC/m^2/d


                !'gp 05-Jul-05
                !'heat fluxes (converted from cal/cm^2/d to W/m^2)
                pr%pr_saveHeatFluxJsnt(i, pr%nj) = saveHeatFluxJsnt(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxLongat(i, pr%nj) = saveHeatFluxLongat(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxBack(i, pr%nj) = saveHeatFluxBack(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxConv(i, pr%nj) = saveHeatFluxConv(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxEvap(i, pr%nj) = saveHeatFluxEvap(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxJsed(i, pr%nj) = saveHeatFluxJsed(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxJhyporheic(i, pr%nj) = saveHeatFluxJhyporheic(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxTribs(i, pr%nj) = saveHeatFluxTribs(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxAdvecDisp(i, pr%nj) = saveHeatFluxAdvecDisp(i) * (4.183076 * 100 * 100 / 86400)
                !'DO fluxes (gO2/m^2/d)
                pr%pr_saveDOfluxReaer(i, pr%nj) = saveDOfluxReaer(i)
                pr%pr_saveDOfluxCBODfast(i, pr%nj) = saveDOfluxCBODfast(i)
                pr%pr_saveDOfluxCBODslow(i, pr%nj) = saveDOfluxCBODslow(i)
                pr%pr_saveDOfluxNitrif(i, pr%nj) = saveDOfluxNitrif(i)
                pr%pr_saveDOfluxPhytoResp(i, pr%nj) = saveDOfluxPhytoResp(i)
                pr%pr_saveDOfluxPhytoPhoto(i, pr%nj) = saveDOfluxPhytoPhoto(i)
                pr%pr_saveDOfluxBotalgResp(i, pr%nj) = saveDOfluxBotalgResp(i)
                pr%pr_saveDOfluxBotalgPhoto(i, pr%nj) = saveDOfluxBotalgPhoto(i)
                pr%pr_saveDOfluxSOD(i, pr%nj) = saveDOfluxSOD(i)
                pr%pr_saveDOfluxCOD(i, pr%nj) = saveDOfluxCOD(i)
                pr%pr_saveDOfluxHyporheic(i, pr%nj) = HypoFluxDO(i)
                pr%pr_saveDOfluxTribs(i, pr%nj) = saveDOfluxTribs(i)
                pr%pr_saveDOfluxAdvecDisp(i, pr%nj) = saveDOfluxAdvecDisp(i)
                !'CO2 fluxes (gC/m^2/d)
                pr%pr_saveCO2fluxReaer(i, pr%nj) = saveCO2fluxReaer(i)
                pr%pr_saveCO2fluxCBODfast(i, pr%nj) = saveCO2fluxCBODfast(i)
                pr%pr_saveCO2fluxCBODslow(i, pr%nj) = saveCO2fluxCBODslow(i)
                pr%pr_saveCO2fluxPhytoResp(i, pr%nj) = saveCO2fluxPhytoResp(i)
                pr%pr_saveCO2fluxPhytoPhoto(i, pr%nj) = saveCO2fluxPhytoPhoto(i)
                pr%pr_saveCO2fluxBotalgResp(i, pr%nj) = saveCO2fluxBotalgResp(i)
                pr%pr_saveCO2fluxBotalgPhoto(i, pr%nj) = saveCO2fluxBotalgPhoto(i)
                pr%pr_saveCO2fluxSOD(i, pr%nj) = saveCO2fluxSOD(i)
                pr%pr_saveCO2fluxHyporheic(i, pr%nj) = HypoFluxIC(i)
                pr%pr_saveCO2fluxTribs(i, pr%nj) = saveCO2fluxTribs(i)
                pr%pr_saveCO2fluxAdvecDisp(i, pr%nj) = saveCO2fluxAdvecDisp(i)


            END DO

            !gp 17-Feb-05
            !x x x x x x x x x x x x x x x x x x x
        END IF
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    END SUBROUTINE Save_init_step

!--------------------------------------------------------------------------------------------
    !save an intermediate step during the last day of the simulation

    !gp 17-Feb-05
    !SUBROUTINE Save_a_step(nr, intg, pr, dt, IMethpH, Rates)
    SUBROUTINE Save_a_step(nr, intg, pr, dt, IMethpH, showDielResults, Rates)


        INTEGER(i32), INTENT(IN) :: nr
        REAL(r64), INTENT(IN) :: dt
        TYPE(rates_t), INTENT(IN) :: Rates
        !CHARACTER(LEN=30) :: IMeth = 'Euler' ! integration method
        CHARACTER(LEN=30), INTENT(IN) :: IMethpH ! pH solve method

        !gp 17-Feb-05
        CHARACTER(LEN=30), INTENT(IN) :: showDielResults !Yes or No (only used in Excel VBA)

        TYPE(Integral_type), INTENT(IN) :: intg
        TYPE(outdata_t), INTENT(INOUT) :: pr
        INTEGER(i32) i, j, k
        REAL(r64) Kamm, NH3, pH, TotP, TotN
        REAL(r64) K1, K2, KW, Kh, CO2sat
        !gp 28-Oct-04 REAL(r64) :: pHss(0:nr)
        REAL(r64) :: pHss(0:nr, nl) !gp

        !gp 29-Oct-04
        TYPE(RiverHydraulics_type) hydrau !channel dimensions, hydraulics, physical characters

        DO i = 0, nr

            !gp 27-Oct-04
            !gp !update maximum and minimum water temperature
            DO j = 1, nl
                IF (intg%Te(i, j) < pr%Temn(i, j)) THEN
                    pr%Temn(i, j) = intg%Te(i, j)
                ELSE IF (intg%Te(i, j) > pr%Temx(i, j)) THEN
                    pr%Temx(i, j) = intg%Te(i, j)
                END IF
                !update average water temperature
                pr%Teav(i, j) = pr%Teav(i, j) + intg%Te(i, j) * dt
                !update average saturation DO
                pr%osav(i, j) = pr%osav(i, j) + os(i, j) * dt
                CALL ChemRates(intg%Te(i, j), K1, K2, KW, Kh, intg%c(i, 1, j))
                CO2sat = Kh * Rates%pco2 !co2 saturation concentration
                CALL ModFP2(4.0_r64, 12.0_r64, pH, CO2sat, intg%Te(i, j), intg%c(i, nv - 2, j), &
                    intg%c(i, 1, j)) !solving pH under co2 saturation condition
                !update average pH under co2 saturation condition
                pr%pHsav(i, j) = pr%pHsav(i, j) + pH * dt
                !
                DO k = 1, nv
                    !update constituents max and min concentrations or values
                    IF (intg%c(i, k, j) < pr%cmn(i, k, j)) THEN
                        pr%cmn(i, k, j) = intg%c(i, k, j)
                    ELSE IF (intg%c(i, k, j) > pr%cmx(i, k, j)) THEN
                        pr%cmx(i, k, j) = intg%c(i, k, j)
                    END IF
                    !update average constituents concentrations or values
                    pr%cav(i, k, j) = pr%cav(i, k, j) + intg%c(i, k, j) * dt
                END DO

                !solve pH value for now

                !gp 23-Nov-09
                !IF (IMethpH == "Newton-Raphson") THEN
                ! CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !ELSE
                ! CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                ! intg%c(i, 1, j))
                !END IF
                IF (IMethpH == "Newton-Raphson") THEN
                    CALL phsolNewton(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                ELSEIF (IMethpH == "Bisection") THEN
                    CALL pHsolBisect(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                ELSE
                    CALL pHsolBrent(pH, intg%c(i, nv - 1, j), intg%Te(i, j), intg%c(i, nv - 2, j), &
                        intg%c(i, 1, j))
                END IF

                pHss(i, j) = pH
                IF (pH < pr%pHmn(i, j)) THEN
                    pr%pHmn(i, j) = pH
                ELSE IF (pH > pr%pHmx(i, j)) THEN
                    pr%pHmx(i, j) = pH
                END IF
                !update average pH value
                pr%pHav(i, j) = pr%pHav(i, j) + pH * dt

                !// NH3 //
                Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (intg%Te(i, j) + 273.15_r64)))
                NH3 = 1.0_r64 / (1.0_r64 + 10.0_r64 ** (-pH) / Kamm) * intg%c(i, 7, j)
                IF (NH3 < pr%NH3mn(i, j)) THEN
                    pr%NH3mn(i, j) = NH3
                ELSE IF (NH3 > pr%NH3mx(i, j)) THEN
                    pr%NH3mx(i, j) = NH3
                END IF
                !average value of NH3
                pr%NH3av(i, j) = pr%NH3av(i, j) + NH3 * dt

                !// Total phosphorus //
                TotP = intg%c(i, 9, j) + intg%c(i, 10, j) + intg%c(i, 11, j) * Rates%apa
                pr%TPav(i, j) = pr%TPav(i, j) + TotP * dt
                IF (TotP < pr%TPmn(i, j)) THEN
                    pr%TPmn(i, j) = TotP
                ELSE IF (TotP > pr%TPmx(i, j)) THEN
                    pr%TPmx(i, j) = TotP
                END IF

                !// Total nitrogen //
                TotN = intg%c(i, 6, j) + intg%c(i, 7, j) + intg%c(i, 8, j) + &
                    intg%c(i, 11, j) * Rates%ana
                pr%TNav(i, j) = pr%TNav(i, j) + TotN * dt
                IF (TotN < pr%TNmn(i, j)) THEN
                    pr%TNmn(i, j) = TotN
                ELSEIF (TotN > pr%TNmx(i, j)) THEN
                    pr%TNmx(i, j) = TotN
                END IF
            END DO !gp 27-Oct-04 end block of new code
            !'gp 15-Nov-04 reach average daily average sediment fluxes
            IF (i > 0) THEN
                pr%HypoFluxDOav(i) = pr%HypoFluxDOav(i) + HypoFluxDO(i) * dt !'dissolved oxygen gO2/m^2/d
                pr%HypoFluxCBODav(i) = pr%HypoFluxCBODav(i) + HypoFluxCBOD(i) * dt !'fast CBOD gO2/m^2/d
                pr%HypoFluxNH4av(i) = pr%HypoFluxNH4av(i) + HypoFluxNH4(i) * dt !'ammonia mgN/m^2/d
                pr%HypoFluxNO3av(i) = pr%HypoFluxNO3av(i) + HypoFluxNO3(i) * dt !'nitrate mgN/m^2/d
                pr%HypoFluxSRPav(i) = pr%HypoFluxSRPav(i) + HypoFluxSRP(i) * dt !'SRP mgP/m^2/d
                pr%DiagFluxDOav(i) = pr%DiagFluxDOav(i) + DiagFluxDO(i) * dt !'dissolved oxygen gO2/m^2/d
                pr%DiagFluxCBODav(i) = pr%DiagFluxCBODav(i) + DiagFluxCBOD(i) * dt !'fast CBOD gO2/m^2/d
                pr%DiagFluxNH4av(i) = pr%DiagFluxNH4av(i) + DiagFluxNH4(i) * dt !'ammonia mgN/m^2/d
                pr%DiagFluxNO3av(i) = pr%DiagFluxNO3av(i) + DiagFluxNO3(i) * dt !'nitrate mgN/m^2/d
                pr%DiagFluxSRPav(i) = pr%DiagFluxSRPav(i) + DiagFluxSRP(i) * dt !'SRP mgP/m^2/d
                !'gp 11-Jan-05 mgN/gD and mgP/gD
                If (intg%c(i, nv, 1) > 0) Then
                    pr%NINbav(i) = pr%NINbav(i) + dt * intg%INb(i) / intg%c(i, nv, 1)
                    If (intg%INb(i) / intg%c(i, nv, 1) < pr%NINbmn(i)) pr%NINbmn(i) = intg%INb(i) / intg%c(i, nv, 1)
                    If (intg%INb(i) / intg%c(i, nv, 1) > pr%NINbmx(i)) pr%NINbmx(i) = intg%INb(i) / intg%c(i, nv, 1)
                    pr%NIPbav(i) = pr%NIPbav(i) + dt * intg%IPb(i) / intg%c(i, nv, 1)
                    If (intg%IPb(i) / intg%c(i, nv, 1) < pr%NIPbmn(i)) pr%NIPbmn(i) = intg%IPb(i) / intg%c(i, nv, 1)
                    If (intg%IPb(i) / intg%c(i, nv, 1) > pr%NIPbmx(i)) pr%NIPbmx(i) = intg%IPb(i) / intg%c(i, nv, 1)
                End If

                !gp 25-Jun-09
                pr%av_BotAlgPhoto(i) = pr%av_BotAlgPhoto(i) + saveBotAlgPhoto(i) * dt !average BotAlgPhoto/Asb gD/m2/d
                pr%av_BotAlgResp(i) = pr%av_BotAlgResp(i) + saveBotAlgResp(i) * dt !average BotAlgReso/Asb gD/m2/d
                pr%av_BotAlgDeath(i) = pr%av_BotAlgDeath(i) + saveBotAlgDeath(i) * dt !average BotAlgDeath/Asb gD/m2/d
                pr%av_BotAlgNetGrowth(i) = pr%av_BotAlgNetGrowth(i) + saveBotAlgNetGrowth(i) * dt !average BotAlg net growth (photo-resp-death) gD/m2/d

            END IF !'gp 15-Nov-04 end new block
        END DO

        !gp 17-Feb-05
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        IF (showDielResults == "Yes") THEN
            !x x x x x x x x x x x x x x x x x x x

            pr%nj = pr%nj + 1
            pr%tdy(pr%nj)= pr%tdy(pr%nj-1) + dt
            DO i = 0, nr
                DO j = 1, nl
                    pr%Tepr(i, pr%nj, j) = intg%Te(i, j)
                    DO k = 1, nv
                        pr%cpr(i, k, pr%nj, j) = intg%c(i, k, j)
                    END DO
                    pr%pHpr(i, pr%nj, j) = pHss(i, j)
                END DO !gp 27-Oct-04
            END DO
            DO i = 1, nr
                !gp 02-Nov-04 pr%NINbpr(i, pr%nj) = intg%INb(i) / intg%c(i, nv)
                !gp pr%NIPbpr(i, pr%nj) = intg%IPb(i) / intg%c(i, nv)
                pr%NINbpr(i, pr%nj) = intg%INb(i) / intg%c(i, nv, 1)
                pr%NIPbpr(i, pr%nj) = intg%IPb(i) / intg%c(i, nv, 1) !gp 02-Nov-04 end new block
                pr%phitotalSavepr(i, pr%nj) = phitotalSave(i) !gp 20-Oct-04
                pr%phitSavepr(i, pr%nj) = phitSave(i) !gp 20-Oct-04
                pr%philSavepr(i, pr%nj) = philSave(i) !gp 20-Oct-04
                pr%phinSavepr(i, pr%nj) = phinSave(i) !gp 20-Oct-04
                pr%phipSavepr(i, pr%nj) = phipSave(i) !gp 20-Oct-04
                pr%phicSavepr(i, pr%nj) = phicSave(i) !gp 20-Oct-04
                !gp 28-Oct-04 output hyporheic fluxes of constituents (total surface area)
                pr%HypoFluxDOpr(i, pr%nj) = HypoFluxDO(i) !'dissolved oxygen gO2/m^2/d
                pr%HypoFluxCBODpr(i, pr%nj) = HypoFluxCBOD(i) !'fast CBOD gO2/m^2/d
                pr%HypoFluxNH4pr(i, pr%nj) = HypoFluxNH4(i) !'ammonia mgN/m^2/d
                pr%HypoFluxNO3pr(i, pr%nj) = HypoFluxNO3(i) !'nitrate mgN/m^2/d
                pr%HypoFluxSRPpr(i, pr%nj) = HypoFluxSRP(i) !'SRP mgP/m^2/d
                pr%HypoFluxICpr(i, pr%nj) = HypoFluxIC(i) !'inorganic C gC/m^2/d
                !gp 01-Nov-04 output diagenesis fluxes of constituents (total surface area)
                pr%DiagFluxDOpr(i, pr%nj) = DiagFluxDO(i) !'dissolved oxygen gO2/m^2/d
                pr%DiagFluxCBODpr(i, pr%nj) = DiagFluxCBOD(i) !'fast CBOD gO2/m^2/d
                pr%DiagFluxNH4pr(i, pr%nj) = DiagFluxNH4(i) !'ammonia mgN/m^2/d
                pr%DiagFluxNO3pr(i, pr%nj) = DiagFluxNO3(i) !'nitrate mgN/m^2/d
                pr%DiagFluxSRPpr(i, pr%nj) = DiagFluxSRP(i) !'SRP mgP/m^2/d
                pr%DiagFluxICpr(i, pr%nj) = DiagFluxIC(i) !'inorganic C gC/m^2/d

                !'gp 05-Jul-05
                !'heat fluxes (converted from cal/cm^2/d to W/m^2)
                pr%pr_saveHeatFluxJsnt(i, pr%nj) = saveHeatFluxJsnt(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxLongat(i, pr%nj) = saveHeatFluxLongat(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxBack(i, pr%nj) = saveHeatFluxBack(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxConv(i, pr%nj) = saveHeatFluxConv(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxEvap(i, pr%nj) = saveHeatFluxEvap(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxJsed(i, pr%nj) = saveHeatFluxJsed(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxJhyporheic(i, pr%nj) = saveHeatFluxJhyporheic(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxTribs(i, pr%nj) = saveHeatFluxTribs(i) * (4.183076 * 100 * 100 / 86400)
                pr%pr_saveHeatFluxAdvecDisp(i, pr%nj) = saveHeatFluxAdvecDisp(i) * (4.183076 * 100 * 100 / 86400)
                !'DO fluxes (gO2/m^2/d)
                pr%pr_saveDOfluxReaer(i, pr%nj) = saveDOfluxReaer(i)
                pr%pr_saveDOfluxCBODfast(i, pr%nj) = saveDOfluxCBODfast(i)
                pr%pr_saveDOfluxCBODslow(i, pr%nj) = saveDOfluxCBODslow(i)
                pr%pr_saveDOfluxNitrif(i, pr%nj) = saveDOfluxNitrif(i)
                pr%pr_saveDOfluxPhytoResp(i, pr%nj) = saveDOfluxPhytoResp(i)
                pr%pr_saveDOfluxPhytoPhoto(i, pr%nj) = saveDOfluxPhytoPhoto(i)
                pr%pr_saveDOfluxBotalgResp(i, pr%nj) = saveDOfluxBotalgResp(i)
                pr%pr_saveDOfluxBotalgPhoto(i, pr%nj) = saveDOfluxBotalgPhoto(i)
                pr%pr_saveDOfluxSOD(i, pr%nj) = saveDOfluxSOD(i)
                pr%pr_saveDOfluxCOD(i, pr%nj) = saveDOfluxCOD(i)
                pr%pr_saveDOfluxHyporheic(i, pr%nj) = HypoFluxDO(i)
                pr%pr_saveDOfluxTribs(i, pr%nj) = saveDOfluxTribs(i)
                pr%pr_saveDOfluxAdvecDisp(i, pr%nj) = saveDOfluxAdvecDisp(i)
                !'CO2 fluxes (gC/m^2/d)
                pr%pr_saveCO2fluxReaer(i, pr%nj) = saveCO2fluxReaer(i)
                pr%pr_saveCO2fluxCBODfast(i, pr%nj) = saveCO2fluxCBODfast(i)
                pr%pr_saveCO2fluxCBODslow(i, pr%nj) = saveCO2fluxCBODslow(i)
                pr%pr_saveCO2fluxPhytoResp(i, pr%nj) = saveCO2fluxPhytoResp(i)
                pr%pr_saveCO2fluxPhytoPhoto(i, pr%nj) = saveCO2fluxPhytoPhoto(i)
                pr%pr_saveCO2fluxBotalgResp(i, pr%nj) = saveCO2fluxBotalgResp(i)
                pr%pr_saveCO2fluxBotalgPhoto(i, pr%nj) = saveCO2fluxBotalgPhoto(i)
                pr%pr_saveCO2fluxSOD(i, pr%nj) = saveCO2fluxSOD(i)
                pr%pr_saveCO2fluxHyporheic(i, pr%nj) = HypoFluxIC(i)
                pr%pr_saveCO2fluxTribs(i, pr%nj) = saveCO2fluxTribs(i)
                pr%pr_saveCO2fluxAdvecDisp(i, pr%nj) = saveCO2fluxAdvecDisp(i)

            END DO

            !gp 17-Feb-05
            !x x x x x x x x x x x x x x x x x x x
        END IF
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    END SUBROUTINE Save_a_step

!----------------------------------------------------------------------------------------------
!make appropriate step
!Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust
!stepsize.

    !gp 08-Jan-10
    !SUBROUTINE odeint (begin, sys, Rates, Meteo, Solar, &
    ! HW, DB, hydrau, pr, nr, t1, t2, dt1, saveSteps)
    SUBROUTINE odeint (begin, sys, Rates, Meteo, Solar, &
        HW, DB, hydrau, pr, nr, t1, t2, dt1, saveSteps, stateVariables)

        !constants
        REAL(r64), PARAMETER :: TINY=1.0E-30
        INTEGER(i32), PARAMETER :: MAXSTP=1000
        !dummy variables
        INTEGER(i32), INTENT(IN) :: nr

        !gp 08-Jan-10
        CHARACTER(LEN=30), INTENT(IN) :: stateVariables !to test for 'All except temperature'

        TYPE(meteorology_t) :: Meteo
        TYPE(solar_type) Solar !solar radiation
! TYPE(HeadwaterDownstream_type), INTENT(IN) :: HDboundary
        TYPE(upstream_boundary_t) HW !headwater
        TYPE(downstream_boundary_t) DB !downstream boundary
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau !channel dimensions, hydraulics, physical characters
        TYPE(SystemParams) sys
        TYPE(Integral_type), INTENT(INOUT) :: begin
        TYPE(rates_t), INTENT(IN) :: Rates
        TYPE(outdata_t), INTENT(OUT) :: pr !output data
        REAL(r64), INTENT(IN) :: t1, t2, dt1
        LOGICAL(4), INTENT(IN) :: saveSteps
        !local variables
        REAL(r64) t, dt !initial time
        INTEGER(i32) nstp, i, j
        REAL(r64) dtnext, dtdid
        TYPE(Integral_type) scal, now
        REAL(r64), DIMENSION(nr,nl) :: Tenow
        REAL(r64), DIMENSION(nr) :: INbnow, IPbnow
        !gp 28-Oct-04 REAL(r64) dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64) dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp
        INTEGER(i32) kmax, kout

        t=t1 !initial time
        dt=dt1

        !constructor
        now = integration_(nr)
        scal= integration_(nr)
        !initialize
        now = begin

        DO nstp = 1, MAXSTP
            ! WRITE(*,*) t
            CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, now%Te, now%c, now%INb, now%IPb, &
                Rates, dTe, dc, dINb, dIPb, t)

            !scaling to monitor accuracy.
            !gp 28-Oct-04 scal%c(1:nr, 1:nv) = ABS(now%c(1:nr, 1:nv)) + ABS(dc(1:nr, 1:nv) * dt) + TINY
            scal%c(1:nr, 1:nv, 1:nl) = ABS(now%c(1:nr, 1:nv, 1:nl)) + ABS(dc(1:nr, 1:nv, 1:nl) * dt) + TINY !gp
            scal%Te(1:nr, 1:nl) = ABS(now%Te(1:nr, 1:nl)) +ABS(dTe(1:nr, 1:nl) *dt) + TINY
            scal%INb(1:nr) = ABS(now%INb(1:nr)) +ABS(dINb(1:nr) *dt) + TINY
            scal%IPb(1:nr) = ABS(now%IPb(1:nr)) +ABS(dIPb(1:nr) *dt) + TINY
            ! if stepsize can overshoot, decrease
            IF (t + dt > t2) dt = t2 - t
            !adaptive step

            !gp 08-Jan-10
            !CALL rkqs(nr, Meteo, Solar, HW, DB, hydrau, sys, now, &
            ! scal, Rates, dTe, dc, dINb, dIPb, t, dt, dtdid,dtnext)
            CALL rkqs(nr, Meteo, Solar, HW, DB, hydrau, sys, now, &
                scal, Rates, dTe, dc, dINb, dIPb, t, dt, dtdid, dtnext, stateVariables)

            IF (saveSteps) THEN

                !gp 17-Feb-05
                !CALL Save_a_step(nr, now, pr, dtdid, sys%IMethpH, Rates)
                CALL Save_a_step(nr, now, pr, dtdid, sys%IMethpH, sys%showDielResults, Rates)

            END IF
            IF (t - t2 >= 0) THEN !are we done?
                begin= now
                RETURN
            END IF
            dt = dtnext
        END DO

    END SUBROUTINE odeint

!-----------------------------------------------------------------------------------------

    !gp 08-Jan-10
    !SUBROUTINE rkqs(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
    ! scal, Rates, dTe, dc, dINb, dIPb, t, dttry, dtdid, dtnext)
    SUBROUTINE rkqs(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
        scal, Rates, dTe, dc, dINb, dIPb, t, dttry, dtdid, dtnext, stateVariables)


        REAL(r64), PARAMETER :: SAFETY=0.9
        REAL(r64), PARAMETER :: PGROW=-0.2
        REAL(r64), PARAMETER :: PSHRNK=-0.25
        REAL(r64), PARAMETER :: ERRCON=0.000189
        REAL(r64), PARAMETER :: EPSTOL=0.0001

        INTEGER(i32), INTENT(IN) :: nr

        !gp 08-Jan-10
        CHARACTER(LEN=30), INTENT(IN) :: stateVariables !to test for 'All except temperature'

        TYPE(meteorology_t) :: Meteo
        TYPE(solar_type) Solar !solar radiation
        TYPE(upstream_boundary_t), INTENT(IN) :: HW !headwater
        TYPE(downstream_boundary_t), INTENT(IN) :: DB !downstream boundary
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau !channel dimensions, hydraulics, physical characters
        TYPE(SystemParams) sys
        TYPE(Integral_type), INTENT(INOUT) :: intg
        TYPE(Integral_type), INTENT(IN) :: scal
        TYPE(rates_t), INTENT(IN) :: Rates
        !gp 28-Oct-04 REAL(r64), INTENT(IN):: dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64), INTENT(IN):: dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp
        REAL(r64), INTENT(INOUT) :: t !initial time, double
        REAL(r64), INTENT(IN) :: dttry
        REAL(r64), INTENT(OUT) :: dtdid, dtnext
        INTEGER(i32) i, j, k
        REAL(r64) errmax, dt, dttemp, tnew

        TYPE(Integral_type) tmp, err !current value and incremented value

        tmp = integration_(nr)
        err = integration_(nr)

        !set the stepsize to the initial trial value
        dt = dttry
        DO WHILE (.TRUE.)
            !take a step

            !gp 08-Jan-10
            !CALL rkck(nr, Meteo, Solar, HW, DB, hydrau, sys, intg, &
            ! Rates, dTe, dc, dINb, dIPb, t, dt, tmp, err)
            CALL rkck(nr, Meteo, Solar, HW, DB, hydrau, sys, intg, &
                Rates, dTe, dc, dINb, dIPb, t, dt, tmp, err, stateVariables)

            !gp 28-Oct-04 errmax = 0.0
            !gp !evaluate accuracy
            !gp DO i = 1 , nr
            !gp DO j = 1 , nv -1
            !gp errmax = MAX(errmax, ABS(err%c(i, j) / scal%c(i, j)))
            !gp ! write(*,*) i, j, errmax
            !gp END DO
            !gp END DO
            !gp
            !gp ! write (*,*) 'temperature'
            !gp DO i = 1 , nr
            !gp DO k=1, nl
            !gp errmax= MAX(errmax, ABS(err%Te(i,k)/scal%Te(i,k)))
            !gp ! write(*,*) i,k, errmax
            !gp END DO
            !gp ! errmax= MAX(errmax, ABS(err%INb(i)/scal%INb(i)))
            !gp ! errmax= MAX(errmax, ABS(err%IPb(i)/scal%IPb(i)))
            !gp END DO
            !evaluate accuracy
            errmax = 0.0
            DO i = 1 , nr
                DO j = 1, nl
                    errmax= MAX(errmax, ABS(err%Te(i, j)/scal%Te(i, j)))
                    DO k = 1 , nv -1
                        errmax = MAX(errmax, ABS(err%c(i, k, j) / scal%c(i, k, j)))
                    END DO
                END DO
            END DO !gp end new block

            !scale relative to required tolerance
            errmax = errmax / EPSTOL
            IF (errmax <= 1.0) THEN
                Exit !Step succeeded, Compute size of next step
            END IF
            dttemp = SAFETY * dt * (errmax ** PSHRNK)

            !Truncation error too large, reduce stepsize
            dt = MAX(dttemp, 0.1 * dt) !no more than a factor of 10.

            tnew = t + dt
            IF (tnew == t) THEN
                ! write(*,*) 'dt = ', dt, 'l = ', l
                STOP "Stepsize underflow in rkqs"
            END IF
        END DO

        ! l=l+1
        IF (errmax > ERRCON) THEN
            dtnext = SAFETY * dt * (errmax ** PGROW)
        ELSE
            dtnext = 5.0 * dt !no more than a factor of 5 increase
        END IF

        dtdid = dt
        t = t + dt
        IF (dt < 0) THEN
            STOP 'Step size underflow in rkqs'
        END IF

! off next 4 lines by hua 09/09/04
! c(1:nr,1:nv) = tmp%c(1:nr,1:nv)
! Te(1:nr,1:nl) = tmp%Te(1:nr,1:nl)
! INb(1:nr) = tmp%INb(1:nr)
! IPb(1:nr) = tmp%IPb(1:nr)
! os=tmp%os
! hua then add
        intg=tmp
! WRITE(*,*) t
    END SUBROUTINE

!---------------------------------------------------------------------------------------------

    !gp 08-Jan-10
    !SUBROUTINE rkck(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
    ! Rates, dTe, dc, dINb, dIPb, t, dt, outIntg, err)
    SUBROUTINE rkck(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
        Rates, dTe, dc, dINb, dIPb, t, dt, outIntg, err, stateVariables)

        INTEGER(i32), INTENT(IN) :: nr

        !gp 08-Jan-10
        CHARACTER(LEN=30), INTENT(IN) :: stateVariables !to test for 'All except temperature'

        TYPE(meteorology_t) :: Meteo
        TYPE(solar_type) Solar !solar radiation
        ! TYPE(HeadwaterDownstream_type), INTENT(IN) :: HDboundary
        TYPE(upstream_boundary_t), INTENT(IN) :: HW !headwater
        TYPE(downstream_boundary_t), INTENT(IN) :: DB !downstream boundary
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau !channel dimensions, hydraulics, physical characters
        TYPE(SystemParams) sys
        TYPE(Integral_type), INTENT(IN) :: intg
        TYPE(Integral_type), INTENT(OUT) :: outIntg, err !current value and incremented value
        TYPE(rates_t), INTENT(IN) :: Rates

        !gp 28-Oct-04 REAL(r64), INTENT(IN):: dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64), INTENT(IN):: dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp
        REAL(r64), INTENT(IN) :: t, dt !initial time, and step size
        INTEGER(i32) i,j,k
        REAL(r64), PARAMETER :: A2 = 0.2, A3 = 0.3, A4 = 0.6, A5 = 1.0, A6 = 0.875, &
            B21 = 0.2, B31 = 3./40., B32 = 9./ 40., B41 = 0.3, &
            b42 = -0.9, b43 = 1.2, b51 = -11./54., b52 = 2.5, &
            b53 = -70./27., b54 = 35./27., b61 = 1631./55296., b62 = 175./512., &
            b63 = 575./13824., b64 = 44275./110592., b65 = 253./4096., c1 = 37./378., &
            c3 = 250./621., c4 = 125./594., c6 = 512./1771., dc5 = -277./14336., &
            dc1 =c1 - 2825./ 27648., dc3 = c3 - 18575./48384., dc4 = c4 - 13525./55296., &
            dc6 = c6 - 0.25
        TYPE(Integral_type) temp
        !Runge-Kutta derivatives
        !gp 28-Oct-04 REAL(r64) ak2dc(nr, nv), ak2dTe(nr, nl), ak2dINb(nr), ak2dIPb(nr)
        !gp REAL(r64) ak3dc(nr, nv), ak3dTe(nr, nl), ak3dINb(nr), ak3dIPb(nr)
        !gp REAL(r64) ak4dc(nr, nv), ak4dTe(nr, nl), ak4dINb(nr), ak4dIPb(nr)
        !gp REAL(r64) ak5dc(nr, nv), ak5dTe(nr, nl), ak5dINb(nr), ak5dIPb(nr)
        !gp REAL(r64) ak6dc(nr, nv), ak6dTe(nr, nl), ak6dINb(nr), ak6dIPb(nr)
        REAL(r64) ak2dc(nr, nv, nl), ak2dTe(nr, nl), ak2dINb(nr), ak2dIPb(nr)
        REAL(r64) ak3dc(nr, nv, nl), ak3dTe(nr, nl), ak3dINb(nr), ak3dIPb(nr)
        REAL(r64) ak4dc(nr, nv, nl), ak4dTe(nr, nl), ak4dINb(nr), ak4dIPb(nr)
        REAL(r64) ak5dc(nr, nv, nl), ak5dTe(nr, nl), ak5dINb(nr), ak5dIPb(nr)
        REAL(r64) ak6dc(nr, nv, nl), ak6dTe(nr, nl), ak6dINb(nr), ak6dIPb(nr) !gp end new block

        !call constructors
        temp= Integration_(nr)

        !gp 28-Oct-04 temp%c(1:nr,1:nv) = intg%c(1:nr,1:nv) + b21 * dt * dc(1:nr,1:nv)
        temp%c(1:nr,1:nv,1:nl) = intg%c(1:nr,1:nv,1:nl) + b21 * dt * dc(1:nr,1:nv,1:nl) !gp
        temp%Te(1:nr,1:nl) = intg%Te(1:nr,1:nl) + b21* dt * dTe(1:nr,1:nl)
        temp%INb(1:nr) = intg%INb(1:nr) + b21* dt * dINb(1:nr)
        temp%IPb(1:nr) = intg%IPb(1:nr) + b21* dt * dIPb(1:nr)

        !Second Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak2dTe, ak2dc, ak2dINb, ak2dIPb, t+ A2*dt)
        DO i=1, nr

            !gp DO j=1, nv
            !gp temp%c(i,j) = intg%c(i,j) + dt* (b31* dc(i,j) + b32 * ak2dc(i,j))
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i,k) = intg%Te(i,k) + dt* (b31* dTe(i,k) + b32 * ak2dTe(i,k))
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt* (b31* dTe(i,j) + b32 * ak2dTe(i,j))
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt* (b31* dc(i,k,j) + b32 * ak2dc(i,k,j))
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt* (b31* dINb(i) + b32 * ak2dINb(i))
            temp%IPb(i) = intg%IPb(i) + dt* (b31* dINb(i) + b32 * ak2dIPb(i))
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Third Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak3dTe, ak3dc, ak3dINb, ak3dIPb, t + A3*dt)
        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b41* dc(i, j) + b42 * ak2dc(i, j) + &
            !gp b43 * ak3dc(i, j))
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt * (b41 * dTe(i, k) + b42 * ak2dTe(i, k) + &
            !gp b43 * ak3dTe(i, k))
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt * (b41 * dTe(i,j) + b42 * ak2dTe(i,j) + &
                    b43 * ak3dTe(i, j))
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b41* dc(i,k,j) + b42 * ak2dc(i,k,j) + &
                        b43 * ak3dc(i,k,j))
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt * (b41 * dINb(i) + b42 * ak2dINb(i) + &
                b43 * ak3dINb(i))
            temp%IPb(i) = intg%IPb(i) + dt * (b41 * dIPb(i) + b42 * ak2dIPb(i) + &
                b43 * ak3dIPb(i))
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Fourth Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak4dTe, ak4dc, ak4dINb, ak4dIPb, t + A4 * dt)

        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b51 * dc(i, j) + b52 * ak2dc(i, j) + &
            !gp b53 * ak3dc(i, j) + b54 * ak4dc(i, j))
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt * (b51 * dTe(i, k) + &
            !gp b52 * ak2dTe(i, k) + b53 * ak3dTe(i, k) + b54 * ak4dTe(i, k))
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt * (b51 * dTe(i,j) + &
                    b52 * ak2dTe(i,j) + b53 * ak3dTe(i,j) + b54 * ak4dTe(i,j))
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b51 * dc(i,k,j) + b52 * ak2dc(i,k,j) + &
                        b53 * ak3dc(i,k,j) + b54 * ak4dc(i,k,j))
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt * (b51 * dINb(i) + &
                b52 * ak2dINb(i) + b53 * ak3dINb(i) + b54 * ak4dINb(i))
            temp%IPb(i) = intg%IPb(i) + dt * (b51 * dIPb(i) + &
                b52 * ak2dIPb(i) + b53 * ak3dIPb(i) + b54 * ak4dIPb(i))
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Fifth Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak5dTe, ak5dc, ak5dINb, ak5dIPb, t + A5 * dt)
        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt * (b61 * dc(i, j) + b62 * ak2dc(i, j) + &
            !gp b63 * ak3dc(i, j) + b64 * ak4dc(i, j) + b65 * ak5dc(i, j))
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt * (b61 * dTe(i, k) + b62 * ak2dTe(i, k) + &
            !gp b63 * ak3dTe(i, k) + b64 * ak4dTe(i, k) + b65 * ak5dTe(i, k))
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt * (b61 * dTe(i,j) + b62 * ak2dTe(i,j) + &
                    b63 * ak3dTe(i,j) + b64 * ak4dTe(i,j) + b65 * ak5dTe(i,j))
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * (b61 * dc(i,k,j) + b62 * ak2dc(i,k,j) + &
                        b63 * ak3dc(i,k,j) + b64 * ak4dc(i,k,j) + b65 * ak5dc(i,k,j))
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt * (b61 * dINb(i) + b62 * ak2dINb(i) + &
                b63 * ak3dINb(i) + b64 * ak4dINb(i) + b65 * ak5dINb(i))
            temp%IPb(i) = intg%IPb(i) + dt * (b61 * dIPb(i) + b62 * ak2dIPb(i) + &
                b63 * ak3dIPb(i) + b64 * ak4dIPb(i) + b65 * ak5dIPb(i))
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Sixth Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak6dTe, ak6dc, ak6dINb, ak6dIPb, t + A6 * dt)

        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp outIntg%c(i, j) = intg%c(i, j) + dt * (c1 * dc(i, j) + c3 * ak3dc(i, j) + &
            !gp c4 * ak4dc(i, j) + c6 * ak6dc(i, j))
            !gp IF (outIntg%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp outIntg%Te(i, k) = intg%Te(i, k) + dt * (c1 * dTe(i, k) + c3 * ak3dTe(i, k) + &
            !gp c4 * ak4dTe(i, k) + c6 * ak6dTe(i, k))
            !gp END DO
            DO j=1, nl

                !gp 08-Jan-10
                !outIntg%Te(i,j) = intg%Te(i,j) + dt * (c1 * dTe(i,j) + c3 * ak3dTe(i,j) + &
                ! c4 * ak4dTe(i,j) + c6 * ak6dTe(i,j))
                If (stateVariables /= "All except temperature") Then
                    outIntg%Te(i,j) = intg%Te(i,j) + dt * (c1 * dTe(i,j) + c3 * ak3dTe(i,j) + &
                        c4 * ak4dTe(i,j) + c6 * ak6dTe(i,j))
                End If

                DO k=1, nv
                    outIntg%c(i,k,j) = intg%c(i,k,j) + dt * (c1 * dc(i,k,j) + c3 * ak3dc(i,k,j) + &
                        c4 * ak4dc(i,k,j) + c6 * ak6dc(i,k,j))
                    IF (outIntg%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO

            outIntg%INb(i) = intg%INb(i) + dt * (c1 * dINb(i) + c3 * ak3dINb(i) + &
                c4 * ak4dINb(i) + c6 * ak6dINb(i))
            outIntg%IPb(i) = intg%IPb(i) + dt * (c1 * dIPb(i) + c3 * ak3dIPb(i) + &
                c4 * ak4dIPb(i) + c6 * ak6dIPb(i))
            IF (outIntg%INb(i)<0) temp%INb(i)=1.0E-6
            IF (outIntg%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO
! outIntg%os=intg%os

        !gp 28-Oct-04 err%c(1:nr, 1:nv) = dt * (DC1 * dc(1:nr, 1:nv) + DC3 * ak3dc(1:nr, 1:nv) + &
        !gp DC4 * ak4dc(1:nr, 1:nv) + DC5 * ak5dc(1:nr, 1:nv) + &
        !gp DC6 * ak6dc(1:nr, 1:nv))
        err%c(1:nr, 1:nv, 1:nl) = dt * (DC1 * dc(1:nr, 1:nv, 1:nl) + DC3 * ak3dc(1:nr, 1:nv, 1:nl) + &
            DC4 * ak4dc(1:nr, 1:nv, 1:nl) + DC5 * ak5dc(1:nr, 1:nv, 1:nl) + &
            DC6 * ak6dc(1:nr, 1:nv, 1:nl)) !gp
        err%Te(1:nr, 1:nl) = dt * (DC1 * dTe(1:nr, 1:nl) + DC3 * ak3dTe(1:nr, 1:nl) + &
            DC4 * ak4dTe(1:nr, 1:nl) + DC5 * ak5dTe(1:nr, 1:nl) + &
            DC6 * ak6dTe(1:nr, 1:nl))
        err%INb(1:nr) = dt * (DC1 * dINb(1:nr) + DC3 * ak3dINb(1:nr) + DC4 * ak4dINb(1:nr) + &
            DC5 * ak5dINb(1:nr)+ DC6 * ak6dINb(1:nr))
        err%IPb(1:nr) = dt * (DC1 * dIPb(1:nr) + DC3 * ak3dIPb(1:nr) + DC4 * ak4dIPb(1:nr) + &
            DC5 * ak5dIPb(1:nr) + DC6 * ak6dIPb(1:nr))

    END SUBROUTINE

!-------------------------------------------------------------------------------------------
! Runge-Kutta 4th order integration

    !gp 08-Jan-10
    !SUBROUTINE rk4(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
    ! Rates, dTe, dc, dINb, dIPb, t, dt)
    SUBROUTINE rk4(nr, Meteo, Solar, HW, DB, hydrau, sys, intg,&
        Rates, dTe, dc, dINb, dIPb, t, dt, stateVariables)

        INTEGER(i32), INTENT(IN) :: nr
        TYPE(meteorology_t) :: Meteo
        TYPE(solar_type) Solar !solar radiation
        !TYPE(HeadwaterDownstream_type), INTENT(IN) :: HDboundary
        TYPE(upstream_boundary_t), INTENT(IN) :: HW !headwater
        TYPE(downstream_boundary_t), INTENT(IN) :: DB !downstream boundary
        TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau !channel dimensions, hydraulics, physical characters
        TYPE(SystemParams) sys
        TYPE(Integral_type), INTENT(INOUT) :: intg
        TYPE(rates_t), INTENT(IN) :: Rates

        !gp 28-Oct-04 REAL(r64), INTENT(OUT):: dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64), INTENT(OUT):: dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp

        REAL(r64), INTENT(IN) :: t, dt !initial time, and step size

        !gp 08-Jan-10
        CHARACTER(LEN=30), INTENT(IN) :: stateVariables !to test for 'All except temperature'

        INTEGER(i32) i,j,k

        TYPE(Integral_type) temp
        !Runge-Kutta derivatives
        !gp 28-Oct-04 REAL(r64) ak2dc(nr, nv), ak2dTe(nr, nl), ak2dINb(nr), ak2dIPb(nr)
        !gp REAL(r64) ak3dc(nr, nv), ak3dTe(nr, nl), ak3dINb(nr), ak3dIPb(nr)
        !gp REAL(r64) ak4dc(nr, nv), ak4dTe(nr, nl), ak4dINb(nr), ak4dIPb(nr)
        REAL(r64) ak2dc(nr, nv, nl), ak2dTe(nr, nl), ak2dINb(nr), ak2dIPb(nr)
        REAL(r64) ak3dc(nr, nv, nl), ak3dTe(nr, nl), ak3dINb(nr), ak3dIPb(nr)
        REAL(r64) ak4dc(nr, nv, nl), ak4dTe(nr, nl), ak4dINb(nr), ak4dIPb(nr) !gp end new block

        !call constructors
        temp= Integration_(nr)

        !First Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, intg%Te, intg%c, intg%INb, intg%IPb, &
            Rates, dTe, dc, dINb, dIPb, t)
        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i, j) = intg%c(i, j) + dt/2 * dc(i, j)
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt/2 * dTe(i, k)
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt/2 * dTe(i,j)
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt/2 * dc(i,k,j)
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt/2 * dINb(i)
            temp%IPb(i) = intg%IPb(i) + dt/2 * dIPb(i)
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Second Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak2dTe, ak2dc, ak2dINb, ak2dIPb, t+ dt/2)

        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i,j) = intg%c(i,j) + dt/2 * ak2dc(i,j)
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt/2 * ak2dTe(i,k)
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt/2 * ak2dTe(i,j)
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt/2 * ak2dc(i,k,j)
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt/2 * ak2dINb(i)
            temp%IPb(i) = intg%IPb(i) + dt/2 * ak2dIPb(i)
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Third Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak3dTe, ak3dc, ak3dINb, ak3dIPb, t + dt/2)

        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp temp%c(i,j) = intg%c(i,j) + dt * ak3dc(i,j)
            !gp IF (temp%c(i,j)<0) temp%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp temp%Te(i, k) = intg%Te(i, k) + dt * ak3dTe(i,k)
            !gp END DO
            DO j=1, nl
                temp%Te(i,j) = intg%Te(i,j) + dt * ak3dTe(i,j)
                DO k=1, nv
                    temp%c(i,k,j) = intg%c(i,k,j) + dt * ak3dc(i,k,j)
                    IF (temp%c(i,k,j)<0) temp%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            temp%INb(i) = intg%INb(i) + dt * ak3dINb(i)
            temp%IPb(i) = intg%IPb(i) + dt * ak3dIPb(i)
            IF (temp%INb(i)<0) temp%INb(i)=1.0E-6
            IF (temp%IPb(i)<0) temp%IPb(i)=1.0E-6
        END DO

        !Fourth Step
        CALL derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, temp%Te, temp%c, temp%INb, temp%IPb, &
            Rates, ak4dTe, ak4dc, ak4dINb, ak4dIPb, t + dt)

        DO i=1, nr

            !gp 28-Oct-04 DO j=1, nv
            !gp intg%c(i, j) = intg%c(i, j) + dt/6 * (dc(i, j) + &
            !gp 2 * ak2dc(i, j) + 2 * ak3dc(i, j) + ak4dc(i, j))
            !gp IF (intg%c(i,j)<0) intg%c(i,j)=1.0E-6
            !gp END DO
            !gp DO k=1, nl
            !gp intg%Te(i, k) = intg%Te(i, k) + dt/6 * (dTe(i, k) + &
            !gp 2 * ak2dTe(i, k) + 2 * ak3dTe(i, k) + ak4dTe(i, k))
            !gp END DO
            DO j=1, nl

                !gp 08-Jan-10
                !intg%Te(i,j) = intg%Te(i,j) + dt/6 * (dTe(i,j) + &
                ! 2 * ak2dTe(i,j) + 2 * ak3dTe(i,j) + ak4dTe(i,j))
                If (stateVariables /= "All except temperature") Then
                    intg%Te(i,j) = intg%Te(i,j) + dt/6 * (dTe(i,j) + &
                        2 * ak2dTe(i,j) + 2 * ak3dTe(i,j) + ak4dTe(i,j))

                    !gp 12-Jan-10
                Else
                    If (hydrau%reach(i)%Te_ini < 0) Then
                        !use diel headwater as diel temp for each reach during intg
                        !CALL InstanteousHeadwater(HW, t, intg%Te(0,1), intg%c(0,:,1), pH)
                        intg%Te(i, j) = intg%Te(0, 1)
                    Else
                        !use initial temp for each reach as const temp during intg
                        intg%Te(i, j) = hydrau%reach(i)%Te_ini
                    End If

                End If

                DO k=1, nv
                    intg%c(i,k,j) = intg%c(i,k,j) + dt/6 * (dc(i,k,j) + &
                        2 * ak2dc(i,k,j) + 2 * ak3dc(i,k,j) + ak4dc(i,k,j))
                    IF (intg%c(i,k,j)<0) intg%c(i,k,j)=1.0E-6
                END DO
            END DO !gp end new block

            intg%INb(i) = intg%INb(i) + dt/6 * (dINb(i) + &
                2 * ak2dINb(i) + 2 * ak3dINb(i) + ak4dINb(i))
            intg%IPb(i) = intg%IPb(i) + dt/6 * (dIPb(i) + &
                2 * ak2dIPb(i) + 2 * ak3dIPb(i) + ak4dIPb(i))
            IF (intg%INb(i)<0) intg%INb(i)=1.0E-6
            IF (intg%IPb(i)<0) intg%IPb(i)=1.0E-6
        END DO
        !outIntg%os=intg%os
    END SUBROUTINE rk4


END MODULE Class_Integration
