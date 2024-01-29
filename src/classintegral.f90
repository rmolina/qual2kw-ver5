!classIntegration

MODULE Class_Integration
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use class_hydraulics, only: riverhydraulics_type
    USE Class_IntegrationData, only: Integral_type, integration_
    USE Class_LightHeat, only: lightheat    
    USE Class_Phsolve, ONLY: ChemRates, ModFP2, phsolNewton, pHsolBisect, phsolbrent, ct
    USE Class_SolarCalc, only: solar_type, solarcalc
    USE Class_SourceIn, only: load, sourcescalc
    USE Class_SystemParams, only: SystemParams
    USE m_downstream_boundary, only: downstream_boundary_t, instanteousdownstreamboundary
    USE m_meteorology, only: meteorology_t, instanteousmeteo
    USE m_output, ONLY: outdata_t, saveheatfluxtribs
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
        INTEGER(i32) nrp
        REAL(r64) TOC, TKN, TSS, TP, TN, BottomAlgae, DOSat, NH3
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
                    CALL Derivs(nr, Meteo, Solar,HW, DB, hydrau, sys, intg%Te, intg%c, intg%INb, &
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
    SUBROUTINE Derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, Te, c, INb, IPb, &
        Rates, dTe, dc, dINb, dIPb, t)

        IMPLICIT NONE
        INTEGER(i32), INTENT(IN) :: nr
        TYPE(meteorology_t) :: Meteo
        TYPE(solar_type) Solar !solar radiation
        TYPE(upstream_boundary_t) HW !headwater
        TYPE(downstream_boundary_t) DB !downstream boundary
        TYPE(RiverHydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
        TYPE(SystemParams) sys
        REAL(r64), DIMENSION(:), POINTER :: INb, IPb
        !REAL(r64), DIMENSION(:), POINTER :: phitotalSave, phitSave, philSave, phinSave, phipSave, phicSave !gp 20-Oct-04
        !gp 27-Oct-04 REAL(r64), DIMENSION(:,:), POINTER :: Te, c
        REAL(r64), DIMENSION(:,:), POINTER :: Te !gp 27-Oct-04
        REAL(r64), DIMENSION(:,:,:), POINTER :: c !gp add dimension for nl
        TYPE(rates_t) Rates
        !gp 27-Oct-04 REAL(r64), INTENT(OUT):: dTe(nr, nl), dc(nr, nv), dINb(nr), dIPb(nr)
        REAL(r64), INTENT(OUT):: dTe(nr, nl), dc(nr, nv, nl), dINb(nr), dIPb(nr) !gp
        REAL(r64), INTENT(IN):: t !time

        INTEGER(i32) i,j, k
        REAL(r64) :: evap, back, flux
        REAL(r64) eair, longat, conv
        !Temperature corrected rates

        !gp 03-Apr-08
        !REAL(r64) khcT(nr), kdcsT(nr), kdcT(nr), kgaFT(nr), kdeaFT(nr), kreaFT(nr), kexaFT(nr)
        REAL(r64) khcT(nr), kdcsT(nr), kdcT(nr), kgaFT(nr), kdeaFT(nr), krea1FT(nr), krea2FT(nr), kexaFT(nr)

        REAL(r64) khnT(nr), khpT(nr), knT(nr), kdtT(nr), kpathT(nr)
        REAL(r64) kiT(nr), kaT(nr), kgaT(nr), kdeaT(nr), kreaT(nr), vdiT(nr), kacT(nr)

        !gp 30-Nov-04
        REAL(r64) kgenT(nr) !gp 30-Nov-04 generic constituent decay

        REAL(r64) :: pH=0
        REAL(r64) fcarb, fnitr, fdenitr, frespp, frespb
        REAL(r64) ke, phip, phint
        REAL(r64) phil, num, den
        REAL(r64) alpha0, alpha1

        REAL(r64) prefam, prefamF
        REAL(r64) SOD, JCH4, JNH4, JNO3, JPO4

        REAL(r64) Jcin, Jnin, Jpin
        REAL(r64) ow, NH3w, NO3w, PO4w, Tw
        REAL(r64) Jamm, Jnitr, Jmeth, Jmethg, Jphos

        REAL(r64) botlight
        REAL(r64) dam, dropft, defa, defb, Iat(nr) !changed by Hua

        !pH
        REAL(r64) K1, K2, KW, Kh, CO2sat

        !gp 01-Nov-07
        !REAL(r64) hh, alp0, alp1, cHCO3CO3
        REAL(r64) hh, alp0, alp1, alp2, cHCO3CO3

        REAL(r64) :: CSOD =0

        REAL(r64) DetrDiss, DetrSettl
        REAL(r64) CBODsHydr, CBODsOxid, CBODfOxid
        REAL(r64) OrgNHydr, NH4Nitrif, Denitr
        REAL(r64) OrgPHydr
        REAL(r64) OrgNSettl, OrgPSettl, InorgPSettl
        REAL(r64) InorgSettl
        REAL(r64) OxReaer
        REAL(r64) PhytoPhoto, PhytoResp, PhytoDeath, PhytoSettl
        REAL(r64) BotAlgPhoto, BotAlgResp, BotAlgDeath, BotAlgExc
        REAL(r64) BotAlgUptakeN, BotAlgUptakeP
        REAL(r64) :: NINb = 0, NIPb =0, FBnb =0, FBpb =0
        REAL(r64) :: phic =0
        REAL(r64) :: PLIM, NLIM

        !gp 15-Nov-04 INTEGER(i32) :: hco3useF

        !pathogens
        REAL(r64) :: ksol=0

        !gp 27-Oct-04 REAL(r64) pHs(0:nr), K1s(0:nr), K2s(0:nr), Khs(0:nr) !New 08/29/04 for less pH solver be called
        REAL(r64) pHs(0:nr, nl), K1s(0:nr, nl), K2s(0:nr, nl), Khs(0:nr, nl) !gp add dim for nl

        !gp 15-Nov-04 move Ehyporheic to hydrau
        !gp 20-Oct-04 sediment and hyporheic heat flux variables
        !gp REAL(r64) Jsed(nr), Jhyporheic(nr), Ehyporheic(nr)
        REAL(r64) Jsed(nr), Jhyporheic(nr)

        !gp 03-Nov-04 hyporheic kinetics variables
        REAL(r64) kgaHT(nr), fcarbH, HeteroGrow

        !gp 15-Nov-04 additional variables for level 2 hyporheic biofilm kinetics
        REAL(r64) kreaHT(nr), kdeaHT(nr), HeteroResp, HeteroDeath, prefamH

        !gp 08-Dec-04 COD oxidation if generic constituent is used as COD
        REAL(r64) CODoxid

        !gp 03-Dec-09 variables for fraction of ionized NH4+ and phosphate speciation
        REAL(r64) Fi, Kamm
        REAL(r64) KPO41, KPO42,KPO43
        REAL(r64) DPO
        REAL(r64) FPO41, FPO42,FPO43
        REAL(r64) K1NH3, KgNH3, HeNH3, vnh3, naus, NH3gas
        REAL(r64) Pcharge, Uw10

        dTe =0; dc =0; dINb =0; dIPb =0
        !gp move kawind to this sub to use interpolated hourly wind speed
        !REAL(r64) kawind

        sys%tday = t - Int(t)

        CALL InstanteousMeteo(nr, t, Meteo)

        !SCC add the following loop to calculate pH and save results 08/29/04
        DO i=0, nr

            !gp 23-Nov-09
            !IF (sys%IMethpH == "Newton-Raphson") THEN
            ! CALL phsolNewton(pH, c(i, nv - 1, 1), Te(i, 1), c(i, nv - 2, 1), c(i, 1, 1))
            !ELSE
            ! CALL phsolBisect(pH, c(i, nv - 1, 1), Te(i, 1), c(i, nv - 2, 1), c(i, 1, 1))
            !END IF
            IF (sys%IMethpH == "Newton-Raphson") THEN
                CALL phsolNewton(pH, c(i, nv - 1, 1), Te(i, 1), c(i, nv - 2, 1), c(i, 1, 1))
            ELSEIF (sys%IMethpH == "Bisection") THEN
                CALL phsolBisect(pH, c(i, nv - 1, 1), Te(i, 1), c(i, nv - 2, 1), c(i, 1, 1))
            ELSE
                CALL phsolBrent(pH, c(i, nv - 1, 1), Te(i, 1), c(i, nv - 2, 1), c(i, 1, 1))
            END IF

            pHs(i, 1) = pH
            pHs(i, 2) = 0 !not used unless hyporheic WQ is simulated
            CALL ChemRates(Te(i, 1), K1, K2, KW, Kh, c(i, 1, 1))
            K1s(i, 1) = K1; K2s(i, 1) = K2; Khs(i, 1) = Kh
            K1s(i, 2) = 0; K2s(i, 2) = 0; Khs(i, 2) = 0 !gp 02-Nov-04 end of new block of code
        END DO

        ! Uw(0) = Uw(1)

        !gp evaluate point source sine functions and distribute loads to reaches for the current time t
        CALL SourcesCalc(t, nr, Hydrau%flag)

        CALL SolarCalc(nr, Solar, Meteo, hydrau, sys) !SolarCalc(tday, i)

        !gp 27-Oct-04 CALL InstanteousHeadwater(HW, t, Te(0,1), c(0,:), pH)
        !gp CALL InstanteousDownstreamboundary(DB, t, Te(nr,1), c(nr,:), Te(nr+1,1), c(nr+1,:))
        CALL InstanteousHeadwater(HW, t, Te(0,1), c(0,:,1), pH)

        !gp 12-Jan-06
        !WRITE(10,*) 'done thru Call InstanteousHeadwater'

        CALL InstanteousDownstreamboundary(DB, t, Te(nr,1), c(nr,:,1), Te(nr+1,1), c(nr+1,:,1)) !gp 27-Oct-04 end new block

        !gp 12-Jan-06
        !WRITE(10,*) 'done thru Call InstanteousDownstreamboundary'
        !WRITE(10,'(32F13.4)') Te(nr,1), (c(nr, k, 1), k=1, nv-2), Te(nr+1,1), (c(nr+1, k, 1), k=1, nv-2)


        !'
        !' ------------------------------------------------------------------------------------
        !' --- Heat transport and flux derivatives for the water column and sediment layers ---
        !' ------------------------------------------------------------------------------------
        !'

        !Heat transport derivatives

        DO i = 1, nr
            ! dTe(i, 1) = hydrau%reach(i-1)%QCMD * Te(i - 1, 1) * rhow * cpw &
            ! - hydrau%reach(i)%QCMD * Te(i, 1) * rhow * cpw &
            ! + hydrau%reach(i)%QptCMD * load(i)%Te * rhow * cpw &
            ! - hydrau%reach(i)%QptaCMD * Te(i, 1) * rhow * cpw &
            ! + hydrau%reach(i-1)%EpCMD * (Te(i - 1, 1) - Te(i, 1)) * rhow * cpw &
            ! + hydrau%reach(i)%EpCMD * (Te(i + 1, 1) - Te(i, 1)) * rhow * cpw
            dTe(i,1) = hydrau%reach(i-1)%QCMD * Te(i - 1, 1) * rhow * cpw
            dTe(i,1) = dTe(i,1) - hydrau%reach(i)%QCMD * Te(i, 1) * rhow * cpw
            dTe(i,1) = dTe(i,1) + hydrau%reach(i)%QptCMD * load(i)%Te * rhow * cpw
            dTe(i,1) = dTe(i,1) - hydrau%reach(i)%QptaCMD * Te(i, 1) * rhow * cpw
            dTe(i,1) = dTe(i,1) + hydrau%reach(i-1)%EpCMD * (Te(i - 1, 1) - Te(i, 1)) * rhow * cpw
            dTe(i,1) = dTe(i,1) + hydrau%reach(i)%EpCMD * (Te(i + 1, 1) - Te(i, 1)) * rhow * cpw

            !'gp 05-Jul-05 (output fluxes in cal/cm^2/d)
            saveHeatFluxTribs(i) = ((hydrau%reach(i)%QptCMD * load(i)%Te &
                - hydrau%reach(i)%QptaCMD * Te(i, 1)) &
                * rhow * cpw * 100 / hydrau%reach(i)%Ast)
            saveHeatFluxAdvecDisp(i) = ((hydrau%reach(i-1)%QCMD * Te(i - 1, 1) &
                - hydrau%reach(i)%QCMD * Te(i, 1) &
                + hydrau%reach(i-1)%EpCMD * (Te(i - 1, 1) - Te(i, 1)) &
                + hydrau%reach(i)%EpCMD * (Te(i + 1, 1) - Te(i, 1))) &
                * rhow * cpw * 100 / hydrau%reach(i)%Ast)

        END DO

        DO i = 1, nr

            !gp 16-Jul-08
            !CALL heatBalance(evap, conv, back, longat, lightheat, Te(i,1), Meteo%Ta(i), &
            ! Meteo%Td(i), Meteo%cc(i), Meteo%Uw(i), hydrau%reach(i)%Ast)
            CALL heatBalance(evap, conv, back, longat, lightheat, Te(i,1), Meteo%Ta(i), &
                Meteo%Td(i), Meteo%cc(i), Meteo%Uw(i), hydrau%reach(i)%Ast, hydrau%reach(i)%SKOP)

            flux = Solar%Jsnt(i) + longat - back - conv - evap !gp changed Jsnt to Jsnt(i)
            dTe(i, 1) = dTe(i, 1) + flux * hydrau%reach(i)%Ast / 100.0_r64
            ! if (i==1) write(8,*) Meteo%Ta(i),Meteo%Td(i) !, Solar%Jsnt(i), longat, back, conv, evap

            !'gp 01-Jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
            saveHeatFluxJsnt(i) = Solar%Jsnt(i)
            saveHeatFluxLongat(i) = longat
            saveHeatFluxBack(i) = -back
            saveHeatFluxConv(i) = -conv
            saveHeatFluxEvap(i) = -evap

        END DO

        DO i = 1 , nr
            !'gp new sediment-water heat flux from sediments into the water in units of cal/cm2/day
            Jsed(i) = (Te(i, 2) - Te(i, 1)) * &
                86400 * 2 * hydrau%reach(i)%sedThermCond / hydrau%reach(i)%HsedCM !'units of (cal/cm2/day), note that sedThermCond is in units of (cal/sec) per (cm deg C)

            !'gp 05-Jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
            saveHeatFluxJsed(i) = Jsed(i)

            dTe(i, 1) = dTe(i, 1) + Jsed(i) * hydrau%reach(i)%Ast/ 100 !'units of (cal/day) * (m3/cm3)
            dTe(i, 2) = -Jsed(i) * hydrau%reach(i)%Ast / 100 !'units of (cal/day) * (m3/cm3)
            SELECT CASE (sys%simHyporheicWQ)
              CASE ('Level 1', 'Level 2') !gp 15-Nov-04
                !'gp new hyporheic exchange flux in units of (cal/day) * (m3/cm3)
                Jhyporheic(i) = hydrau%reach(i)%EhyporheicCMD * &
                    (Te(i, 2) - Te(i, 1)) * rhow * cpw * 100 / hydrau%reach(i)%Ast !'units of (cal/cm2/day)

                !'gp 05-Jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
                saveHeatFluxJhyporheic(i) = Jhyporheic(i)

                dTe(i, 1) = dTe(i, 1) + Jhyporheic(i) * hydrau%reach(i)%Ast / 100 !'units of (cal/day) * (m3/cm3)
                dTe(i, 2) = dTe(i, 2) - Jhyporheic(i) * hydrau%reach(i)%Ast / 100 !'units of (cal/day) * (m3/cm3)
            END SELECT
        END DO
        DO i = 1, nr
            !'gp units for dTe on LHS in the next original line of code is deg C per day
            dTe(i, 1) = dTe(i, 1) / &
                hydrau%reach(i)%Vol / rhow / cpw !'units of deg C per day
            !'gp sediment-water heat flux term dTs in deg C per day
            dTe(i, 2) = dTe(i, 2) / &
                (hydrau%reach(i)%Ast * hydrau%reach(i)%HsedCM / 100) / &
                hydrau%reach(i)%rhoCpSed !'units of deg C per day
        END DO


!gp 03-Feb-05 start of bypass derivs for water quality variables
!unless 'All' state variables are being simulated
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!gp 29-Dec-09
!If (sys%stateVariables == "All") Then
        If (sys%stateVariables /= "Temperature") Then

!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x



            !gp 02-Nov-04 os(0) = oxsat(Te(0, 1), hydrau%reach(0)%elev)
            os(0, 1) = oxygen_saturation(Te(0, 1), hydrau%reach(0)%elev)
            SELECT CASE (sys%simHyporheicWQ)
              CASE ('Level 1', 'Level 2') !gp 15-Nov-04
                os(0, 2) = oxygen_saturation(Te(0, 2), hydrau%reach(0)%elev)
              CASE DEFAULT
                os(0, 2) = 0 !not used in calculation
            END SELECT !gp 02-Nov-04 end new block

            !'
            !' --------------------------------------------------------------------------
            !' --- temperature adjustment of rates and constants for the water column ---
            !' --------------------------------------------------------------------------
            !'

            !gp 03-Apr-08
            !CALL TempAdjust(Rates, hydrau,Meteo, nr, Te, khcT, kdcsT, kdcT, kgaFT, kdeaFT, &
            ! kreaFT, kexaFT, khnT, khpT, knT, kdtT, kpathT, kiT, kaT, kgaT, kdeaT, &
            ! kreaT, vdiT, kacT, kgenT) !gp 30-Nov-04 add kgenT
            CALL temp_adjust(Rates, hydrau,Meteo, nr, Te, khcT, kdcsT, kdcT, kgaFT, kdeaFT, &
                krea1FT, krea2FT, kexaFT, khnT, khpT, knT, kdtT, kpathT, kiT, kaT, kgaT, kdeaT, &
                kreaT, vdiT, kacT, kgenT) !gp 30-Nov-04 add kgenT


            !'gp 20-Oct-04 temp limitation factor for output
            DO i = 1 , nr
                phitSave(i) = Rates%tkgaF ** (Te(i, 1) - 20)
            END DO

            DO i=1 ,nr
                !gp 27-Oct-04 os(i)=oxsat(Te(i,1), hydrau%reach(i)%elev)
                os(i, 1)=oxygen_saturation(Te(i, 1), hydrau%reach(i)%elev) !gp 27-Oct-04
            END DO

            !'
            !' ---------------------------------------------------------------------------
            !' ---------------------------------------------------------------------------
            !' --- ---
            !' --- Water column water quality deriviatives ---
            !' --- Note that units are mass/time in the main i loop ---
            !' --- and at the the derivs are divided by Vol(i) to get mass/volume/time ---
            !' --- or area to get mass/area/time ---
            !' --- ---
            !' ---------------------------------------------------------------------------
            !' ---------------------------------------------------------------------------
            !'

            !transport and small dams
            DO i = 1, nr
                DO k = 1, nv - 1
                    SELECT CASE (k)

                        !gp 31-Mar-05 cT is variable nv-1 (debug Alkalinity deriv)
                      CASE (3,nv-1)

                        dropft = hydrau%reach(i -1)%drop * 3.281_r64 !SCC
                        dam = 1.0_r64 / (1 + 0.116_r64 * adam * bdam * dropft * (1 - 0.034_r64 * dropft) * &
                            (1.0_r64 + 0.046_r64 * Te(i - 1, 1)))
                        IF (k == 3) THEN
                            defa = os(i - 1, 1) - c(i - 1, 3, 1)
                            defb = dam * defa
                            dc(i, k, 1) = hydrau%reach(i-1)%QCMD * (os(i-1, 1) - defb) !gp 27-Oct-04 end new block

                            !'gp 05-Jul-05 (save DO fluxes in gO2/m^2/d)
                            saveDOfluxHeadwater(i) = (hydrau%reach(i-1)%QCMD * (os(i - 1, 1) - defb) &
                                - hydrau%reach(i-1)%QCMD * c(i - 1, k, 1)) / hydrau%reach(i)%Ast

                            !gp 31-Mar-05 cT is variable nv-2
                            !ELSEIF (k == 15) THEN
                        ELSEIF (k == nv-1) THEN

                            CO2sat = Khs(i-1, 1) * Rates%pco2
                            hh = 10.0_r64 ** -pHs(i-1, 1)
                            alp0 = hh * hh / (hh ** 2.0 + K1s(i-1, 1) * hh + K1s(i-1, 1) * K2s(i-1, 1))
                            cHCO3CO3 = (1.0_r64 - alp0) * c(i-1, k, 1)
                            defa = CO2sat - alp0 * c(i-1, k, 1)
                            defb = dam * defa
                            dc(i, k, 1) = hydrau%reach(i-1)%QCMD * (cHCO3CO3 + CO2sat - defb) !gp 27-Oct-04 end new block

                            !'gp 05-Jul-05 (save CO2 fluxes in gC/m^2/d)
                            saveCO2fluxHeadwater(i) = (hydrau%reach(i-1)%QCMD * (cHCO3CO3 + CO2sat - defb) &
                                - hydrau%reach(i-1)%QCMD * c(i - 1, k, 1)) &
                                / hydrau%reach(i)%Ast / Rates%rccc

                        END IF
                      CASE DEFAULT !ELSE
                        dc(i, k, 1) = hydrau%reach(i-1)%QCMD * c(i - 1, k, 1) !gp 27-Oct-04 end new block
                    END SELECT
                    dc(i, k, 1) = dc(i, k, 1) - hydrau%reach(i)%QCMD * c(i, k, 1) &
                        + hydrau%reach(i)%QptCMD * load(i)%c(k) &
                        - hydrau%reach(i)%QptaCMD * c(i, k, 1) &
                        + hydrau%reach(i-1)%EpCMD * (c(i - 1, k, 1) - c(i, k, 1)) &
                        + hydrau%reach(i)%EpCMD * (c(i + 1, k, 1) - c(i, k, 1)) !gp 27-Oct-04 end new block

                    !'gp 05-Jul-05 (save net trib/GW and advec/disp DO fluxes in gO2/m^2/d and CO2 fluxes in gC/m^2/d)
                    Select Case (k)
                      Case (3)
                        saveDOfluxTribs(i) = (hydrau%reach(i)%QptCMD * load(i)%c(k) &
                            - hydrau%reach(i)%QptaCMD * c(i, k, 1)) &
                            / hydrau%reach(i)%Ast
                        saveDOfluxAdvecDisp(i) = saveDOfluxHeadwater(i) &
                            + (-hydrau%reach(i)%QCMD * c(i, k, 1) &
                            + hydrau%reach(i-1)%EpCMD * (c(i - 1, k, 1) - c(i, k, 1)) &
                            + hydrau%reach(i)%EpCMD * (c(i + 1, k, 1) - c(i, k, 1))) &
                            / hydrau%reach(i)%Ast
                      Case (nv - 1)
                        saveCO2fluxTribs(i) = (hydrau%reach(i)%QptCMD * load(i)%c(k) &
                            - hydrau%reach(i)%QptaCMD * c(i, k, 1)) &
                            / hydrau%reach(i)%Ast / Rates%rccc
                        saveCO2fluxAdvecDisp(i) = saveCO2fluxHeadwater(i) &
                            + (-hydrau%reach(i)%QCMD * c(i, k, 1) &
                            + hydrau%reach(i-1)%EpCMD * (c(i - 1, k, 1) - c(i, k, 1)) &
                            + hydrau%reach(i)%EpCMD * (c(i + 1, k, 1) - c(i, k, 1))) &
                            / hydrau%reach(i)%Ast / Rates%rccc
                    End Select

                END DO
            END DO

            ! ----------------
            ! --- kinetics ---
            ! ----------------

            DO i = 1, nr

                ! --- ph dependent variables for later use in deriv calcs in this i loop

                hh = 10.0_r64 ** -pHs(i, 1)
                alp0 = hh * hh / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1))
                alp1 = K1s(i, 1) * hh / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1)) !fraction of cT as HCO3-
                alp2 = K1s(i, 1) * K2s(i, 1) / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1)) !fraction of cT as CO3--

                !'gp 03-Dec-09
                !'fraction of total ammonia that is ionized ammonia NH4+ (Fi)
                Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 1) + 273.15_r64))) !'equilibrium coeff for ammonia dissociation NH4+ = NH3(aq) + H+
                Fi = hh / (hh + Kamm) !'fraction of ionized NH4+ = [NH4+ / (NH4+ + NH3)]
                !'fraction of SRP that is H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                KPO41 = 10.0_r64 ** -2.15_r64 !'= [H+][H2PO4-]/[H3PO4] for eqn H3PO4 = H+ + H2PO4-
                KPO42 = 10.0_r64 ** -7.2_r64 !'= [H+][HPO4--]/[H2PO4-] for eqn H2PO4- = H+ + HPO4--
                KPO43 = 10.0_r64 ** -12.35_r64 !'= [H+][PO4---]/[HPO4--] for eqn HPO4-- = H+ + PO4---
                DPO = 1.0_r64 / (hh ** 3.0_r64 + KPO41 * hh ** 2.0_r64 + KPO41 * KPO42 * hh + KPO41 * KPO42 * KPO43) !'intermediate calc of denominator
                FPO41 = KPO41 * hh ** 2.0_r64 * DPO !'fraction of phosphate as H2PO4-
                FPO42 = KPO41 * KPO42 * hh * DPO !'fraction of phosphate as HPO4--
                FPO43 = KPO41 * KPO42 * KPO43 * DPO !'fraction of phosphate as PO4---
                !'stoichiometric conversion factors for alkalinity (some are pH dependent) (ralkbn & ralkbp calc in classstoch.f90)
                !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqH+/L per mgN/m^3 (1eqH+/moleN / 14gN/moleN / 1000mgN/gN / 1000L/m^3)
                !ralkbp = (FPO41 + 2# * FPO42 + 3# * FPO43) / 30.973762 / 1000# / 1000# !'eqH+/L per mgP/m^3 (eqH+/moleP / 31gP/moleP / 1000mgP/gP / 1000L/m^3)
                Pcharge = (FPO41 + 2.0_r64 * FPO42 + 3.0_r64 * FPO43) !'average charge of P species to multiply by ralkbp
                !'ammonia gas transfer at air/water interface (eqns from Chapra et al QUAL2K manual ver 2.11b8)
                K1NH3 = 1.171_r64 * kaT(i) * hydrau%reach(i)%depth !'liquid film exchange coefficient for NH3 (m/d)
                Uw10 = (10.0_r64 / 7.0_r64) ** 0.15_r64 * Meteo%Uw(i) !'wind speed at 10 m (m/s)
                KgNH3 = 175.287_r64 * Uw10 + 262.9305_r64 !'NH3 mass transfer velocity (m/d)
                HeNH3 = 0.0000136785_r64 * 1.052_r64 ** (Te(i, 1) - 20.0_r64) !'Henry's Law constant for NH3 gas (atm m^3/mole)
                vnh3 = K1NH3 * HeNH3 / (HeNH3 + 0.00008206_r64 * (Te(i, 1) + 273.15_r64) * (K1NH3 / KgNH3)) !'NH3 gas transfer coefficient (m/d)
                naus = 0.000000002_r64 / HeNH3 * 14000.0_r64 !'sat'n conc of NH3 assuming partial press of NH3 in atm = 2e-9 atm (values range from 1-10e-9 in rural and 10-100e-9 in heavily polluted areas)
                NH3gas = vnh3 * hydrau%reach(i)%Ast * (naus - (1.0_r64 - Fi) * c(i, 7, 1)) !'loss or gain of NH3 via gas transfer (mgN/day)

                !determine solar radiation
                Iat(i) = lightheat%PAR * Solar%Jsnt(i)

                CALL oxygen_inhibition_and_enhancement(Rates, c(i, 3, 1), fcarb, fnitr, fdenitr, frespp, frespb) !gp 27-Oct-04

                !light extinction

                !gp 13-Feb-06 include macrophyte extinction if bottom plants are simulated as macrophytes
                IF ((hydrau%reach(i)%NUpWCfrac < 1.0_r64) .OR. (hydrau%reach(i)%PUpWCfrac < 1.0_r64)) THEN
                    !include macrophyte biomass extinction if macrophytes are simulated
                    ke= lightExtinction(lightheat, c(i, 12, 1), c(i, 2, 1), c(i, 11, 1), c(i, nv, 1)/hydrau%reach(i)%depth)
                ELSE
                    !do not include macrophyte biomass extinction of periphyton are simulated
                    ke= lightExtinction(lightheat, c(i, 12, 1), c(i, 2, 1), c(i, 11, 1), 0.0_r64)
                END IF

                !
                ! --- (11) Phytoplankton (mgA/day) ---
                !

                !nutrient limitation

                !gp 09-Dec-09 include Fi
                !IF (hydrau%reach(i)%ksn + c(i, 7, 1) + c(i, 8, 1) <= 0) THEN
                ! phint = 0
                !ELSE
                ! phint = (c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%ksn + c(i, 7, 1) + c(i, 8, 1))
                !END IF
                IF (hydrau%reach(i)%ksn + Fi * c(i, 7, 1) + c(i, 8, 1) <= 0) THEN
                    phint = 0
                ELSE
                    phint = (Fi * c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%ksn + Fi * c(i, 7, 1) + c(i, 8, 1))
                END IF

                IF ((hydrau%reach(i)%ksp + c(i, 10, 1)) <= 0) THEN
                    phip = 0
                ELSE
                    phip = c(i, 10, 1) / (hydrau%reach(i)%ksp + c(i, 10, 1))
                END IF
                IF (phip < phint) phint = phip
                IF (Rates%hco3use == "No") THEN
                    phic = alp0 * c(i, nv-1, 1) / (Rates%ksc + alp0 * c(i, nv-1, 1))
                ELSE
                    phic = (alp0 + alp1) * c(i, nv-1, 1) / (Rates%ksc + (alp0 + alp1) * c(i, nv-1, 1))
                END IF
                IF (phic < phint) phint = phic

                !light limitation

                SELECT CASE (Rates%Ilight)
                  CASE (1) !Half-saturation
                    IF (hydrau%reach(i)%Isat == 0) THEN
                        phil = 1.0
                    ELSE
                        phil = 1.0 / (ke * hydrau%reach(i)%depth) * LOG((hydrau%reach(i)%Isat + Iat(i)) &
                            / (hydrau%reach(i)%Isat + Iat(i) * EXP(-ke * hydrau%reach(i)%depth)))
                    END IF
                  CASE (2) !Smith
                    IF (hydrau%reach(i)%Isat == 0) THEN
                        phil = 1.0
                    ELSE
                        num = Iat(i) / hydrau%reach(i)%Isat + SQRT(1 + (Iat(i) / hydrau%reach(i)%Isat) ** 2)
                        den = Iat(i) * EXP(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%Isat &
                            + SQRT(1 + (Iat(i) * EXP(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%Isat) ** 2)
                        phil = 1.0 / (ke * hydrau%reach(i)%depth) * LOG(num / den)
                    END IF
                  CASE (3) !Steele
                    IF (hydrau%reach(i)%Isat == 0) THEN
                        phil = 1.0
                    ELSE
                        alpha0 = Iat(i) / hydrau%reach(i)%Isat
                        alpha1 = Iat(i) / hydrau%reach(i)%Isat * EXP(-ke * hydrau%reach(i)%depth)
                        phil = EXP(1.0_r64) * (EXP(-alpha1) - EXP(-alpha0)) / (ke * hydrau%reach(i)%depth)
                    END IF
                END SELECT

                PhytoPhoto = phil * phint * kgaT(i) * hydrau%reach(i)%vol * c(i, 11, 1) !mgA/day
                PhytoResp = frespp * kreaT(i) * hydrau%reach(i)%vol * c(i, 11, 1) !mgA/day
                PhytoDeath = kdeaT(i) * hydrau%reach(i)%vol * c(i, 11, 1) !mgA/day
                PhytoSettl = hydrau%reach(i)%va * hydrau%reach(i)%Asd * c(i, 11, 1) !mgA/day
                dc(i, 11, 1) = dc(i, 11, 1) + PhytoPhoto - PhytoResp - PhytoDeath - PhytoSettl !mgA/day

                !
                ! --- (16) Bottom Algae (gD/day) ---
                !

                !luxury uptake calcs

                IF (c(i, nv, 1) > 0) THEN
                    NINb = INb(i) / c(i, nv, 1)
                    NIPb = IPb(i) / c(i, nv, 1)
                ELSE
                    NINb = Rates%ana / Rates%ada
                    NIPb = Rates%apa / Rates%ada
                END IF
                IF (NINb > 0) THEN
                    FBnb = hydrau%reach(i)%NINbmin / NINb
                ELSE
                    FBnb = 1.0_r64
                END IF
                IF (FBnb < 0) FBnb = 0.0_r64
                IF (FBnb > 1.0_r64) FBnb = 1.0_r64
                IF (NIPb > 0.0_r64) THEN
                    FBpb = hydrau%reach(i)%NIPbmin / NIPb
                ELSE
                    FBpb = 1.0_r64
                END IF
                IF (FBpb < 0) FBpb = 0.0_r64
                IF (FBpb > 1.0_r64) FBpb = 1.0_r64

                phint = 1.0_r64 - FBnb
                phinSave(i) = phint !'gp 20-Oct-04
                phip = 1.0_r64 - FBpb
                phipSave(i) = phip !'gp 20-Oct-04
                IF (phip < phint) phint = phip

                IF (Rates%hco3useF == "No") THEN
                    phic = alp0 * c(i, nv-1, 1) / (Rates%kscF + alp0 * c(i, nv-1, 1)) !gp !CO2 is the limiting substrate
                ELSE
                    phic = (alp0 + alp1) * c(i, nv-1, 1) / (Rates%kscF + (alp0 + alp1) * c(i, nv-1, 1)) !gp !HCO3- is the limiting substrate
                END IF
                phicSave(i) = phic !'gp 20-Oct-04
                IF (phic < phint) phint = phic

                !light limitation

                botlight = Iat(i) * EXP(-ke * hydrau%reach(i)%depth)
                IF ((hydrau%reach(i)%NUpWCfrac < 1.0_r64) .OR. (hydrau%reach(i)%PUpWCfrac < 1.0_r64)) THEN
                    !use water column average light for macrophytes
                    SELECT CASE (Rates%IlightF)
                      CASE (1) !Half-saturation
                        IF (hydrau%reach(i)%IsatF == 0) THEN
                            phil = 1.0
                        ELSE
                            phil = 1.0 / (ke * hydrau%reach(i)%depth) * LOG((hydrau%reach(i)%IsatF + Iat(i)) &
                                / (hydrau%reach(i)%IsatF + Iat(i) * EXP(-ke * hydrau%reach(i)%depth)))
                        END IF
                      CASE (2) !Smith
                        IF (hydrau%reach(i)%IsatF == 0) THEN
                            phil = 1.0
                        ELSE
                            num = Iat(i) / hydrau%reach(i)%IsatF + SQRT(1 + (Iat(i) / hydrau%reach(i)%IsatF) ** 2)
                            den = Iat(i) * EXP(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%IsatF &
                                + SQRT(1 + (Iat(i) * EXP(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%IsatF) ** 2)
                            phil = 1.0 / (ke * hydrau%reach(i)%depth) * LOG(num / den)
                        END IF
                      CASE (3) !Steele
                        IF (hydrau%reach(i)%IsatF == 0) THEN
                            phil = 1.0
                        ELSE
                            alpha0 = Iat(i) / hydrau%reach(i)%IsatF
                            alpha1 = Iat(i) / hydrau%reach(i)%IsatF * EXP(-ke * hydrau%reach(i)%depth)
                            phil = EXP(1.0_r64) * (EXP(-alpha1) - EXP(-alpha0)) / (ke * hydrau%reach(i)%depth)
                        END IF
                    END SELECT
                ELSE
                    !use bottom light for periphyton
                    SELECT CASE (Rates%IlightF)
                      CASE (1) !Half-saturation
                        IF (hydrau%reach(i)%IsatF + botlight <= 0) THEN
                            phil = 1.0_r64
                        ELSE
                            phil = botlight / (hydrau%reach(i)%IsatF + botlight)
                        END IF
                      CASE (2) !Smith
                        IF (botlight ** 2 + hydrau%reach(i)%IsatF ** 2 == 0) THEN
                            phil = 1.0_r64
                        ELSE
                            phil = botlight / SQRT(botlight ** 2 + hydrau%reach(i)%IsatF ** 2)
                        END IF
                      CASE (3) !Steele
                        IF (hydrau%reach(i)%IsatF <= 0) THEN
                            phil = 1.0_r64
                        ELSE
                            phil = botlight / hydrau%reach(i)%IsatF * EXP(1 - botlight / hydrau%reach(i)%IsatF)
                        END IF
                    END SELECT
                END IF

                philSave(i) = phil

                IF (Rates%typeF == "Zero-order") THEN
                    BotAlgPhoto = phil * phint * kgaFT(i) * hydrau%reach(i)%Asb !gD/day
                ELSE
                    BotAlgPhoto = phil * phint * kgaFT(i) * hydrau%reach(i)%Asb * c(i, nv, 1) &
                        * (1 - c(i, nv, 1) / hydrau%reach(i)%abmax) !gD/day
                END IF

                phitotalSave(i) = phil * phint * phitSave(i)

                !basal resp
                BotAlgResp = krea1FT(i) * hydrau%reach(i)%Asb * c(i, nv, 1) !gD/day
                !add phot-resp
                If (Rates%typeF == "Zero-order") Then
                    BotAlgResp = BotAlgResp + krea2FT(i) * phil * phint * kgaFT(i) * hydrau%reach(i)%Asb
                Else !'else first-order growth used for photo resp
                    BotAlgResp = BotAlgResp + krea2FT(i) * phil * phint * kgaFT(i) * hydrau%reach(i)%Asb &
                        * c(i, nv, 1) * (1 - c(i, nv, 1) / hydrau%reach(i)%abmax)
                End If
                !'adjust for oxygen attenuation
                BotAlgResp = frespb * BotAlgResp !gD/day

                BotAlgExc = kexaFT(i) * hydrau%reach(i)%Asb * c(i, nv, 1) !gD/day (note that BotAlgExc does not contribute to dc(i,nv,1)
                BotAlgDeath = kdeaFT(i) * hydrau%reach(i)%Asb * c(i, nv, 1) !gD/day
                dc(i, nv, 1) = BotAlgPhoto - BotAlgResp - BotAlgDeath !gD/day

                !gp 25-Jun-09
                saveBotAlgPhoto(i) = BotAlgPhoto / hydrau%reach(i)%Asb
                saveBotAlgResp(i) = BotAlgResp / hydrau%reach(i)%Asb
                saveBotAlgDeath(i) = BotAlgDeath / hydrau%reach(i)%Asb
                saveBotAlgNetGrowth(i) = dc(i, nv, 1) / hydrau%reach(i)%Asb

                !periphyton uptake of N and P

                !gp 03-Dec-09
                !IF ((hydrau%reach(i)%ksnF + c(i, 7, 1) + c(i, 8, 1)) > 0) THEN
                ! NLIM = hydrau%reach(i)%KqN / (hydrau%reach(i)%KqN + NINb - hydrau%reach(i)%NINbmin)
                ! IF (NLIM < 0) THEN
                ! NLIM = 0
                ! ELSE IF (NLIM > 1.0) THEN
                ! NLIM = 1.0_r64
                ! END IF
                ! BotAlgUptakeN = hydrau%reach(i)%Asb * NLIM * c(i, nv, 1) * hydrau%reach(i)%NINbupmax &
                ! * (c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%ksnF + c(i, 7, 1) &
                ! + c(i, 8, 1))
                !ELSE
                ! BotAlgUptakeN = 0
                !END IF
                IF ((hydrau%reach(i)%ksnF + Fi * c(i, 7, 1) + c(i, 8, 1)) > 0) THEN
                    NLIM = hydrau%reach(i)%KqN / (hydrau%reach(i)%KqN + NINb - hydrau%reach(i)%NINbmin)
                    IF (NLIM < 0) THEN
                        NLIM = 0
                    ELSE IF (NLIM > 1.0) THEN
                        NLIM = 1.0_r64
                    END IF
                    BotAlgUptakeN = hydrau%reach(i)%Asb * NLIM * c(i, nv, 1) * hydrau%reach(i)%NINbupmax &
                        * (Fi * c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%ksnF + Fi * c(i, 7, 1) &
                        + c(i, 8, 1))
                ELSE
                    BotAlgUptakeN = 0
                END IF

                IF ((hydrau%reach(i)%kspF + c(i, 10, 1)) > 0) THEN
                    PLIM = hydrau%reach(i)%KqP / (hydrau%reach(i)%KqP + NIPb - hydrau%reach(i)%NIPbmin)
                    IF (PLIM < 0) THEN
                        PLIM = 0
                    ELSE IF (PLIM > 1) THEN
                        PLIM = 1.0_r64
                    END IF
                    BotAlgUptakeP = hydrau%reach(i)%Asb * PLIM * c(i, nv, 1) * hydrau%reach(i)%NIPbupmax &
                        * c(i, 10, 1) / (hydrau%reach(i)%kspF + c(i, 10, 1))
                ELSE
                    BotAlgUptakeP = 0
                END IF

                !gp 03-Dec-09
                !If (BotAlgUptakeN * sys%dt > (c(i, 7, 1) + c(i, 8, 1)) * hydrau%reach(i)%vol) Then
                ! BotAlgUptakeN = (c(i, 7, 1) + c(i, 8, 1)) * hydrau%reach(i)%vol / sys%dt !'mgN/day
                !End If
                !If (BotAlgUptakeP * sys%dt > c(i, 10, 1) * hydrau%reach(i)%vol) Then
                ! BotAlgUptakeP = c(i, 10, 1) * hydrau%reach(i)%vol / sys%dt !'mgP/day
                !End If
                If (BotAlgUptakeN * sys%dt > (Fi * c(i, 7, 1) + c(i, 8, 1)) * hydrau%reach(i)%vol) Then
                    BotAlgUptakeN = (Fi * c(i, 7, 1) + c(i, 8, 1)) * hydrau%reach(i)%vol / sys%dt !'mgN/day
                End If
                If (BotAlgUptakeP * sys%dt > c(i, 10, 1) * hydrau%reach(i)%vol) Then
                    BotAlgUptakeP = c(i, 10, 1) * hydrau%reach(i)%vol / sys%dt !'mgP/day
                End If

                !change in intracellular N and P in periphyton
                dINb(i) = BotAlgUptakeN - NINb * BotAlgDeath - NINb * BotAlgExc !mgN/day
                dIPb(i) = BotAlgUptakeP - NIPb * BotAlgDeath - NIPb * BotAlgExc !mgP/day

                !ammonium preference

                !gp 03-Dec-09
                !prefam = 0
                !IF (c(i, 7, 1) * c(i, 8, 1) > 0) THEN
                ! prefam = c(i, 7, 1) * c(i, 8, 1) / (hydrau%reach(i)%khnx + c(i, 7, 1)) / (hydrau%reach(i)%khnx + c(i, 8, 1)) &
                ! + c(i, 7, 1) * hydrau%reach(i)%khnx / (c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%khnx + c(i, 8, 1))
                !END IF
                !prefamF = 0
                !IF (c(i, 7, 1) + c(i, 8, 1) > 0) THEN
                ! prefamF = c(i, 7, 1) * c(i, 8, 1) / (hydrau%reach(i)%khnxF + c(i, 7, 1)) / (hydrau%reach(i)%khnxF + c(i, 8, 1)) &
                ! + c(i, 7, 1) * hydrau%reach(i)%khnxF / (c(i, 7, 1) + c(i, 8, 1)) / (hydrau%reach(i)%khnxF + c(i, 8, 1))
                !END IF
                prefam = 0
                IF (Fi * c(i, 7, 1) * c(i, 8, 1) > 0) THEN
                    prefam = Fi * c(i, 7, 1) * c(i, 8, 1) / (hydrau%reach(i)%khnx + Fi * c(i, 7, 1)) &
                        / (hydrau%reach(i)%khnx + c(i, 8, 1)) &
                        + Fi * c(i, 7, 1) * hydrau%reach(i)%khnx / (Fi * c(i, 7, 1) + c(i, 8, 1)) &
                        / (hydrau%reach(i)%khnx + c(i, 8, 1))
                END IF
                prefamF = 0
                IF (Fi * c(i, 7, 1) + c(i, 8, 1) > 0) THEN
                    prefamF = Fi * c(i, 7, 1) * c(i, 8, 1) / (hydrau%reach(i)%khnxF + Fi * c(i, 7, 1)) &
                        / (hydrau%reach(i)%khnxF + c(i, 8, 1)) &
                        + Fi * c(i, 7, 1) * hydrau%reach(i)%khnxF / (Fi * c(i, 7, 1) + c(i, 8, 1)) &
                        / (hydrau%reach(i)%khnxF + c(i, 8, 1))
                END IF

                !
                ! --- sediment fluxes of O2, C, N, and P ---
                !

                If (sys%calcSedFlux == "Yes" .OR. sys%calcSedFlux == "Option 1" .OR. sys%calcSedFlux == "Option 2") Then !gp 11-Jan-06

                    Jcin = Rates%roc * Rates%aca * hydrau%reach(i)%va * c(i, 11, 1) &
                        + hydrau%reach(i)%vdt * c(i, 12, 1) / Rates%adc * Rates%roc !gO2/m^2/d
                    Jnin = (Rates%ana * hydrau%reach(i)%va * c(i, 11, 1) + hydrau%reach(i)%von * c(i, 6, 1)) / 1000.0_r64 !gN/m^2/d
                    Jpin = (Rates%apa * hydrau%reach(i)%va * c(i, 11, 1) + hydrau%reach(i)%vop * c(i, 9, 1) &
                        + hydrau%reach(i)%vip * c(i, 10, 1)) / 1000.0_r64 !gP/m^2/d
                    ow = c(i, 3, 1) !mgO2/L = gO2/m^3

                    !gp 09-Dec-09 include Fi
                    !NH3w = c(i, 7, 1) / 1000.0_r64 !gN/m^3
                    NH3w = Fi * c(i, 7, 1) / 1000.0_r64 !gN/m^3

                    NO3w = c(i, 8, 1) / 1000.0_r64 !gN/m^3
                    PO4w = c(i, 10, 1) / 1000.0_r64 !gP/m^3
                    Tw = Te(i, 1) !deg C

                    !note that inputs and outputs of O2, C(O2), N, and P fluxes to/from SedCalcNumNew are g/m^2/day
                    CALL SedCalcNumNew(Jcin, Jnin, Jpin, ow, hydrau%reach(i)%depth, Tw, SOD, Jamm, Jnitr, Jmeth, Jmethg, &
                        Jphos, NH3w, NO3w, PO4w, c(i, 5, 1), CSOD, sys%calcSedFlux)

                    JNH4 = Jamm * 1000.0_r64 !mgN/m^2/d
                    JNO3 = Jnitr * 1000.0_r64 !mgN/m^2/d
                    JCH4 = Jmeth !gO2/m^2/d
                    JPO4 = Jphos * 1000.0_r64 * Rates%kspi / (Rates%kspi + c(i, 3, 1)) !mgP/m^2/d
                END IF

                SODpr(i) = hydrau%reach(i)%SODspec + SOD !gO2/m^2/d
                JNH4pr(i) = hydrau%reach(i)%JNH4spec + JNH4 !mgN/m^2/d
                JNO3pr(i) = JNO3 !mgN/m^2/d
                JCH4pr(i) = hydrau%reach(i)%JCH4spec + JCH4 !gO2/m^2/d
                JSRPpr(i) = hydrau%reach(i)%JSRPspec + JPO4 !mgP/m^2/d
                CSODpr(i) = CSOD !gO2/m^2/d

                !gp 01-Nov-04 output diel diagenesis fluxes of constituents (averaged over total surface area)
                DiagFluxDO(i) = -SODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'dissolved oxygen gO2/m^2/d
                DiagFluxCBOD(i) = JCH4pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'fast CBOD gO2/m^2/d
                DiagFluxNH4(i) = JNH4pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'ammonia mgN/m^2/d
                DiagFluxNO3(i) = JNO3pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'nitrate mgN/m^2/d
                DiagFluxSRP(i) = JSRPpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'SRP mgP/m^2/d
                DiagFluxIC(i) = CSODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast / Rates%roc !'inorganic C gC/m^2/d

                !
                ! --- (12) detritus (gD/day) ---
                !

                DetrDiss = kdtT(i) * hydrau%reach(i)%vol * c(i, 12, 1)
                DetrSettl = hydrau%reach(i)%vdt * hydrau%reach(i)%Asd * c(i, 12, 1)
                dc(i, 12, 1) = dc(i, 12, 1) - DetrDiss
                dc(i, 12, 1) = dc(i, 12, 1) - DetrSettl
                dc(i, 12, 1) = dc(i, 12, 1) + Rates%ada * PhytoDeath
                dc(i, 12, 1) = dc(i, 12, 1) + BotAlgDeath

                !
                ! --- (6) Organic Nitrogen (mgN/day) ---
                !

                OrgNHydr = khnT(i) * hydrau%reach(i)%vol * c(i, 6, 1)
                OrgNSettl = hydrau%reach(i)%von * hydrau%reach(i)%Asd * c(i, 6, 1)
                dc(i, 6, 1) = dc(i, 6, 1) + NINb * BotAlgDeath + Rates%ana * PhytoDeath
                dc(i, 6, 1) = dc(i, 6, 1) - OrgNHydr - OrgNSettl

                !
                ! --- (7) Total Ammonia (NH4+ + NH3) Nitrogen (mgN/day) ---
                !

                !gp 03-Dec-09
                !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%vol * c(i, 7, 1)
                NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%vol * Fi * c(i, 7, 1)

                dc(i, 7, 1) = dc(i, 7, 1) + OrgNHydr
                dc(i, 7, 1) = dc(i, 7, 1) - NH4Nitrif
                dc(i, 7, 1) = dc(i, 7, 1) + Rates%ana * PhytoResp
                dc(i, 7, 1) = dc(i, 7, 1) - prefam * Rates%ana * PhytoPhoto
                dc(i, 7, 1) = dc(i, 7, 1) - prefamF * BotAlgUptakeN * hydrau%reach(i)%NUpWCfrac
                dc(i, 7, 1) = dc(i, 7, 1) + JNH4pr(i) * hydrau%reach(i)%Asd
                dc(i, 7, 1) = dc(i, 7, 1) + NINb * BotAlgExc

                !gp 03-Dec-09
                dc(i, 7, 1) = dc(i, 7, 1) + NH3gas !air/water exchange of NH3 gas

                !
                ! --- (8) Nitrate+Nitrite Nitrogen (mgN/day) ---
                !

                Denitr = c(i, 5, 1) / (0.1_r64 + c(i, 5, 1)) * fdenitr * kiT(i) * hydrau%reach(i)%vol * c(i, 8, 1) !wc denitr
                dc(i, 8, 1) = dc(i, 8, 1) + NH4Nitrif
                dc(i, 8, 1) = dc(i, 8, 1) - Denitr
                dc(i, 8, 1) = dc(i, 8, 1) - (1 - prefam) * Rates%ana * PhytoPhoto
                dc(i, 8, 1) = dc(i, 8, 1) - (1 - prefamF) * BotAlgUptakeN * hydrau%reach(i)%NUpWCfrac
                dc(i, 8, 1) = dc(i, 8, 1) + JNO3pr(i) * hydrau%reach(i)%Asd
                dc(i, 8, 1) = dc(i, 8, 1) - vdiT(i) * hydrau%reach(i)%Asd * c(i, 8, 1) !sed denitr

                !
                ! --- (9) Organic Phosphorus (mgP/day) ---
                !

                OrgPHydr = khpT(i) * hydrau%reach(i)%vol * c(i, 9, 1)
                OrgPSettl = hydrau%reach(i)%vop * hydrau%reach(i)%Asd * c(i, 9, 1)
                dc(i, 9, 1) = dc(i, 9, 1) + NIPb * BotAlgDeath + Rates%apa * PhytoDeath
                dc(i, 9, 1) = dc(i, 9, 1) - OrgPHydr - OrgPSettl

                !
                ! --- (10) Inorganic Soluble Reactive Phosphorus (mgP/day) ---
                !

                dc(i, 10, 1) = dc(i, 10, 1) + OrgPHydr
                dc(i, 10, 1) = dc(i, 10, 1) + Rates%apa * PhytoResp
                dc(i, 10, 1) = dc(i, 10, 1) - Rates%apa * PhytoPhoto
                dc(i, 10, 1) = dc(i, 10, 1) - BotAlgUptakeP * hydrau%reach(i)%PUpWCfrac
                dc(i, 10, 1) = dc(i, 10, 1) + JSRPpr(i) * hydrau%reach(i)%Asd
                InorgPSettl = hydrau%reach(i)%vip * hydrau%reach(i)%Asd * c(i, 10, 1)
                dc(i, 10, 1) = dc(i, 10, 1) - InorgPSettl
                dc(i, 10, 1) = dc(i, 10, 1) + NIPb * BotAlgExc

                !
                ! --- (4) CBOD Slow (gO2/day) ---
                !

                CBODsHydr = khcT(i) * hydrau%reach(i)%vol * c(i, 4, 1)
                dc(i, 4, 1) = dc(i, 4, 1) + Rates%roc * (1.0_r64 / Rates%adc) * DetrDiss
                dc(i, 4, 1) = dc(i, 4, 1) - CBODsHydr
                CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%vol * c(i,4, 1)
                dc(i, 4, 1) = dc(i, 4, 1) - CBODsOxid

                !
                ! --- (5) CBOD Fast (gO2/day) ---
                !

                CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%vol * c(i, 5, 1)
                dc(i, 5, 1) = dc(i, 5, 1) + CBODsHydr
                dc(i, 5, 1) = dc(i, 5, 1) - CBODfOxid
                dc(i, 5, 1) = dc(i, 5, 1) + JCH4pr(i) * hydrau%reach(i)%Asd
                dc(i, 5, 1) = dc(i, 5, 1) - Rates%rondn * Denitr !gp end new block

                !
                ! --- (2) Inorganic Suspended Solids (gD/day) ---
                !

                InorgSettl = hydrau%reach(i)%vss * hydrau%reach(i)%Asd * c(i, 2, 1)
                dc(i, 2, 1) = dc(i, 2, 1) - InorgSettl !gp end new block

                !
                ! --- (14) Generic constituent or COD (user defined units of mass/time, gO2/day if used as COD) ---
                !

                dc(i, 14, 1) = dc(i, 14, 1) - kgenT(i) * hydrau%reach(i)%vol * c(i, 14, 1)
                dc(i, 14, 1) = dc(i, 14, 1) - hydrau%reach(i)%vgen * hydrau%reach(i)%Asd * c(i, 14, 1)
                IF (Rates%useGenericAsCOD == "Yes") THEN
                    CODoxid = kgenT(i) * hydrau%reach(i)%vol * c(i, 14, 1)
                ELSE
                    CODoxid = 0
                END IF

                !
                ! --- (3) Dissolved Oxygen (gO2/day) ---
                !

                OxReaer = kaT(i) * hydrau%reach(i)%vol * (oxygen_saturation(Te(i, 1), hydrau%reach(i)%elev) - c(i, 3, 1))
                dc(i, 3, 1) = dc(i, 3, 1) + OxReaer
                dc(i, 3, 1) = dc(i, 3, 1) - CBODsOxid - CBODfOxid
                dc(i, 3, 1) = dc(i, 3, 1) - Rates%ron * NH4Nitrif
                dc(i, 3, 1) = dc(i, 3, 1) - Rates%roa * PhytoResp
                dc(i, 3, 1) = dc(i, 3, 1) + Rates%roa * PhytoPhoto * prefam
                dc(i, 3, 1) = dc(i, 3, 1) + Rates%roa * PhytoPhoto * (1.0_r64 - prefam) * 138.0_r64 / 107.0_r64
                dc(i, 3, 1) = dc(i, 3, 1) - Rates%roc / Rates%adc * BotAlgResp
                dc(i, 3, 1) = dc(i, 3, 1) + Rates%roc / Rates%adc * BotAlgPhoto * prefamF
                dc(i, 3, 1) = dc(i, 3, 1) + Rates%roc / Rates%adc * BotAlgPhoto * (1 - prefamF) * 138.0_r64 / 107.0_r64
                dc(i, 3, 1) = dc(i, 3, 1) - SODpr(i) * hydrau%reach(i)%Asd
                dc(i, 3, 1) = dc(i, 3, 1) - CODoxid

                !'gp 05-Jul-05 save DO fluxes (gO2/m^2/d)
                saveDOfluxReaer(i) = OxReaer / hydrau%reach(i)%Ast
                saveDOfluxCBODfast(i) = -CBODfOxid / hydrau%reach(i)%Ast
                saveDOfluxCBODslow(i) = -CBODsOxid / hydrau%reach(i)%Ast
                saveDOfluxNitrif(i) = -Rates%ron * NH4Nitrif / hydrau%reach(i)%Ast
                saveDOfluxPhytoResp(i) = -Rates%roa * PhytoResp / hydrau%reach(i)%Ast
                saveDOfluxPhytoPhoto(i) = (Rates%roa * PhytoPhoto * prefam &
                    + Rates%roa * PhytoPhoto * (1 - prefam) * 138.0_r64 / 107.0_r64) &
                    / hydrau%reach(i)%Ast
                saveDOfluxBotalgResp(i) = -Rates%roc / Rates%adc * BotAlgResp / hydrau%reach(i)%Ast
                saveDOfluxBotalgPhoto(i) = (Rates%roc / Rates%adc * BotAlgPhoto * prefamF &
                    + Rates%roc / Rates%adc * BotAlgPhoto * (1 - prefamF) * 138.0_r64 / 107.0_r64) &
                    / hydrau%reach(i)%Ast
                saveDOfluxSOD(i) = -SODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast
                saveDOfluxCOD(i) = -CODoxid / hydrau%reach(i)%Ast

                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) - Rates%rcca * PhytoPhoto
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + Rates%rcca * PhytoResp
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + Rates%rcco * CBODfOxid
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + Rates%rcco * CBODsOxid
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) - Rates%rccd * BotAlgPhoto
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + Rates%rccd * BotAlgResp
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + Rates%rcco * CSOD * hydrau%reach(i)%Asd
                dc(i, nv - 1, 1) = dc(i, nv - 1, 1) + kacT(i) * hydrau%reach(i)%vol * (Khs(i, 1) &
                    * Rates%pco2 - alp0 * c(i, nv - 1, 1)) !gp end new block

                !
                ! --- (nv-2) Alkalinity (gCaCO3/day) ---
                !

                If (sys%simAlk == "Yes") Then

                    !gp 03-Dec-09
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkaa * PhytoPhoto * prefam * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkan * PhytoPhoto * (1.0_r64 - prefam) * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkaa * PhytoResp * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbn * BotAlgUptakeN * prefamF * 50000.0_r64 &
                    ! - Rates%ralkbp * BotAlgUptakeP * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * BotAlgUptakeN * (1.0_r64 - prefamF) * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkn * NH4Nitrif * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkden * Denitr * 50000.0_r64
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * OrgNHydr * 50000.0_r64 + Rates%ralkbp * OrgPHydr * 50000.0_r64 !gp end new block
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkn * JNH4pr(i) * hydrau%reach(i)%Asd * 50000.0_r64 !'sed flux of ammonia
                    !dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbp * JSRPpr(i) * hydrau%reach(i)%Asd * 50000.0_r64 !'sed flux of PO4

                    !gp 03-Dec-09
                    ! alkalinity derivative (factor of 50043.45 converts eqH+/L to mgCaCO3/L or gCaCO3/m^3)
                    ! Fi is already accounted for in nitrification, uptake
                    ! but not in the production of total ammonia (e.g. excretion, OrgN hydrolysis, JNH4pr)
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbn * Fi* Rates%ana * PhytoPhoto * prefam * 50043.45_r64 !'phyto photo uptake of ammonia
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Rates%ana * PhytoPhoto * (1 - prefam) * 50043.45_r64 !'phyto photo uptake of nitrate
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbp * Pcharge * Rates%apa * PhytoPhoto * 50043.45_r64 !'phyto photo uptake of P
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Fi * Rates%ana * PhytoResp * 50043.45_r64 !'phyto resp of N
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbp * Pcharge * Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbn * Fi * BotAlgUptakeN * prefamF &
                        * hydrau%reach(i)%NUpWCfrac * 50043.45_r64 !'periphyton N uptake (ammonia)
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * BotAlgUptakeN * (1 - prefamF) &
                        * hydrau%reach(i)%NUpWCfrac * 50043.45_r64 !'periphyton N uptake (nitrate)
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbp * Pcharge * BotAlgUptakeP &
                        * hydrau%reach(i)%PUpWCfrac * 50043.45_r64 !'periphyton P uptake
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Fi * BotAlgExc * NINb * 50043.45_r64 !'periphyton excretion of N
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbp * Pcharge * BotAlgExc * NIPb * 50043.45_r64 !'periphyton excretion of P
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbn * 2.0_r64 * NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss and nitrate gain
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Denitr * 50043.45_r64 !'water column denitrification
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * vdiT(i) * hydrau%reach(i)%Asd * c(i, 8, 1) * 50043.45_r64 !'sediment denitrification
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Fi * OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbp * Pcharge * OrgPHydr * 50043.45_r64 !'organic P hydrolysis
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbn * Fi * JNH4pr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of ammonia
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbn * JNO3pr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of nitrate
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) - Rates%ralkbp * Pcharge * JSRPpr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of PO4
                    dc(i, nv - 2, 1) = dc(i, nv - 2, 1) + Rates%ralkbp * Pcharge * InorgPSettl * 50043.45_r64 !'settling flux of PO4

                end if

                !'gp 05-Jul-05 save CO2 fluxes (gC/m^2/d)
                saveCO2fluxReaer(i) = (kacT(i) * hydrau%reach(i)%vol &
                    * (Khs(i, 1) * Rates%pco2 - alp0 * c(i, nv - 1, 1))) &
                    / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxCBODfast(i) = Rates%rcco * CBODfOxid / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxCBODslow(i) = Rates%rcco * CBODsOxid / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxPhytoResp(i) = Rates%rcca * PhytoResp / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxPhytoPhoto(i) = -Rates%rcca * PhytoPhoto / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxBotalgResp(i) = Rates%rccd * BotAlgResp / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxBotalgPhoto(i) = -Rates%rccd * BotAlgPhoto / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxSOD(i) = Rates%rcco * CSOD * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast / Rates%rccc

                !
                ! --- Pathogen indicator bacteria (13) ---
                !

                ksol = hydrau%reach(i)%apath * Solar%Jsnt(i) / 24.0_r64 &
                    * (1.0_r64 - EXP(-ke * hydrau%reach(i)%depth)) / (ke * hydrau%reach(i)%depth)
                dc(i, 13, 1) = dc(i, 13, 1) - (kpathT(i) + ksol) * hydrau%reach(i)%vol * c(i, 13, 1)
                dc(i, 13, 1) = dc(i, 13, 1) - hydrau%reach(i)%vpath * hydrau%reach(i)%Asd * c(i, 13, 1)

            END DO


            !
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! --- Hyporheic pore water quality deriviatives ---
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! -------------------------------------------------
            !

            SELECT CASE (sys%simHyporheicWQ)


                !
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                ! --- Level 1 growth kinetics for hyporheic biofilm ---
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                !

              CASE ('Level 1') !gp 15-Nov-04

                !
                ! ------------------------
                ! --- calc and save pH ---
                ! ------------------------
                !

                DO i=0, nr

                    !gp 23-Nov-09
                    !IF (sys%IMethpH == "Newton-Raphson") THEN
                    ! CALL phsolNewton(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    !ELSE
                    ! CALL phsolBisect(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    !END IF
                    IF (sys%IMethpH == "Newton-Raphson") THEN
                        CALL phsolNewton(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    ELSEIF (sys%IMethpH == "Bisection") THEN
                        CALL phsolBisect(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    ELSE
                        CALL phsolBrent(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    END IF

                    pHs(i, 2) = pH
                    CALL ChemRates(Te(i, 2), K1, K2, KW, Kh, c(i, 1, 2))
                    K1s(i, 2) = K1; K2s(i, 2) = K2; Khs(i, 2) = Kh
                END DO

                !
                ! -----------------------------------------------------
                ! --- temperature adjustment of rates and constants ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !gp 08-Feb-06
                    !khcT(i) = Rates%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    !kdcsT(i) = Rates%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    !kdcT(i) = Rates%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    !khnT(i) = Rates%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    !khpT(i) = Rates%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    !knT(i) = Rates%kn * Rates%tkn ** (Te(i, 2) - 20)
                    !kiT(i) = Rates%ki * Rates%tki ** (Te(i, 2) - 20)
                    !os(i, 2)=oxsat(Te(i, 2), hydrau%reach(i)%elev)
                    !kdeaT(i) = Rates%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    !kreaT(i) = Rates%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    !kdtT(i) = Rates%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    !kpathT(i) = Rates%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !enhanced CBOD oxidation rate from hyporheic biofilm
                    !kgenT(i) = Rates%kgen * Rates%tkgen ** (Te(i, 2) - 20) !gp 30-Nov-04
                    khcT(i) = hydrau%reach(i)%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    kdcsT(i) = hydrau%reach(i)%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    kdcT(i) = hydrau%reach(i)%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    khnT(i) = hydrau%reach(i)%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    khpT(i) = hydrau%reach(i)%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    knT(i) = hydrau%reach(i)%kn * Rates%tkn ** (Te(i, 2) - 20)
                    kiT(i) = hydrau%reach(i)%ki * Rates%tki ** (Te(i, 2) - 20)
                    os(i, 2)=oxygen_saturation(Te(i, 2), hydrau%reach(i)%elev)
                    kdeaT(i) = hydrau%reach(i)%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    kreaT(i) = hydrau%reach(i)%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    kdtT(i) = hydrau%reach(i)%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    kpathT(i) = hydrau%reach(i)%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    kgaHT(i) = hydrau%reach(i)%kgaH * Rates%tkgaH ** (Te(i, 2) - 20)
                    kgenT(i) = hydrau%reach(i)%kgen * Rates%tkgen ** (Te(i, 2) - 20)

                END DO

                !
                ! ---------------------------------------------------------------------------------------------
                ! --- mass transfer between water column and hyporheic zone by bulk hyporheic exchange flow ---
                ! ---------------------------------------------------------------------------------------------
                !

                DO i = 1, nr
                    DO k = 1, nv-1
                        !gp un-comment next two lines to activate mass transfer between water column and hyporheic pore water
                        dc(i, k, 1) = dc(i, k, 1) + hydrau%reach(i)%EhyporheicCMD * (c(i, k, 2) - c(i, k, 1))
                        dc(i, k, 2) = hydrau%reach(i)%EhyporheicCMD * (c(i, k, 1) - c(i, k, 2))
                    END DO
                    !store hyporheic fluxes for output (positive flux is source to the water column from the hyporheic zone)
                    HypoFluxDO(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 3, 2) - c(i, 3, 1)) / hydrau%reach(i)%Ast !DO gO2/m^2/d
                    HypoFluxCBOD(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 5, 2) - c(i, 5, 1)) / hydrau%reach(i)%Ast !fast CBOD gO2/m^2/d
                    HypoFluxNH4(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 7, 2) - c(i, 7, 1)) / hydrau%reach(i)%Ast !NH4 mgN/m^2/d
                    HypoFluxNO3(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 8, 2) - c(i, 8, 1)) / hydrau%reach(i)%Ast !NO3 mgN/m^2/d
                    HypoFluxSRP(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 10, 2) - c(i, 10, 1)) / hydrau%reach(i)%Ast !SRP mgP/m^2/d
                    HypoFluxIC(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, nv - 1, 2) - c(i, nv - 1, 1)) &
                        / hydrau%reach(i)%Ast / Rates%rccc !cT gC/m^2/d
                END DO

                !
                ! -----------------------------------------------------
                ! --- mass derivatives for the hyporheic pore water ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !
                    ! --- pH dependent variables for later use in deriv calcs in this i loop
                    !
                    !'gp 03-Dec-09
                    !'fraction of ionized ammonia
                    hh = 10.0_r64 ** (-pHs(i, 2))
                    Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 2) + 273.15_r64)))
                    Fi = hh / (hh + Kamm)
                    !'phosphate fractions of H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                    KPO41 = 10.0_r64 ** (-2.15_r64)
                    KPO42 = 10.0_r64 ** (-7.2_r64)
                    KPO43 = 10.0_r64 ** (-12.35_r64)
                    DPO = 1.0_r64 / (hh ** 3.0_r64 + KPO41 * hh ** 2.0_r64 + KPO41 * KPO42 * hh + KPO41 * KPO42 * KPO43)
                    FPO41 = KPO41 * hh ** 2.0_r64 * DPO !'fraction of phosphate as H2PO4-
                    FPO42 = KPO41 * KPO42 * hh * DPO !'fraction of phosphate as HPO4--
                    FPO43 = KPO41 * KPO42 * KPO43 * DPO !'fraction of phosphate as PO4---
                    !'stoichiometric conversion factors (some are pH dependent)
                    !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqH+/L per mgN/m^3
                    !ralkbp = (FPO41 + 2# * FPO42 + 3# * FPO43) / 30.973762 / 1000# / 1000# !'eqH+/L per mgP/m^3
                    Pcharge = (FPO41 + 2.0_r64 * FPO42 + 3.0_r64 * FPO43) !avg charge of P species for multiplier for ralkbp

                    SELECT CASE (Rates%IkoxC) !'low O2 inhibition of C oxidation by floating heterotrophs
                      CASE (1)
                        fcarb = c(i, 3, 2) / (Rates%Ksocf + c(i, 3, 2))
                      CASE (2)
                        fcarb = 1 - Exp(-Rates%Ksocf * c(i, 3, 2))
                      CASE (3)
                        fcarb = c(i, 3, 2) ** 2 / (Rates%Ksocf + c(i, 3, 2) ** 2)
                    END SELECT

                    !gp 08-Feb-06
                    !SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                    !CASE (1)
                    ! fcarbH = c(i, 3, 2) / (Rates%kinhcH + c(i, 3, 2))
                    !CASE (2)
                    ! fcarbH = 1 - Exp(-Rates%kinhcH * c(i, 3, 2))
                    !CASE (3)
                    ! fcarbH = c(i, 3, 2) ** 2 / (Rates%kinhcH + c(i, 3, 2) ** 2)
                    !END SELECT
                    SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                      CASE (1)
                        fcarbH = c(i, 3, 2) / (hydrau%reach(i)%kinhcH + c(i, 3, 2))
                      CASE (2)
                        fcarbH = 1 - Exp(-hydrau%reach(i)%kinhcH * c(i, 3, 2))
                      CASE (3)
                        fcarbH = c(i, 3, 2) ** 2 / (hydrau%reach(i)%kinhcH + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxN) !'low O2 inhibition of nitrification
                      CASE (1)
                        fnitr = c(i, 3, 2) / (Rates%Ksona + c(i, 3, 2))
                      CASE (2)
                        fnitr = 1 - Exp(-Rates%Ksona * c(i, 3, 2))
                      CASE (3)
                        fnitr = c(i, 3, 2) ** 2 / (Rates%Ksona + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxDN) !'low O2 enhancement of denitrification
                      CASE (1)
                        fdenitr = 1 - c(i, 3, 2) / (Rates%Ksodn + c(i, 3, 2))
                      CASE (2)
                        fdenitr = Exp(-Rates%Ksodn * c(i, 3, 2))
                      CASE (3)
                        fdenitr = 1 - c(i, 3, 2) ** 2 / (Rates%Ksodn + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxP) !'low O2 inhib of phytoplankton respiration
                      CASE (1)
                        frespp = c(i, 3, 2) / (Rates%Ksop + c(i, 3, 2))
                      CASE (2)
                        frespp = 1 - Exp(-Rates%Ksop * c(i, 3, 2))
                      CASE (3)
                        frespp = c(i, 3, 2) ** 2 / (Rates%Ksop + c(i, 3, 2) ** 2)
                    END SELECT

                    !'Phytoplankton (11)
                    PhytoResp = frespp * kreaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 11, 2)
                    PhytoDeath = kdeaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 11, 2)
                    dc(i, 11, 2) = dc(i, 11, 2) - PhytoResp - PhytoDeath

                    !'Enhanced oxidation of fast C in the hyporheic sediment zone
                    !'later versions of Q2K may use c(i, nv, 2) to represent biomass of the biofilm
                    !'this current vesion of Q2K does not simulate the biofilm biomass as a state variable
                    !'limitation from fast DOC and dissolved oxygen:

                    !gp 08-Feb-06
                    !phint = 1
                    !IF ((Rates%kscH + c(i, 5, 2)) <= 0) THEN
                    ! phic = 0
                    !ELSE
                    ! phic = c(i, 5, 2) / (Rates%kscH + c(i, 5, 2))
                    !END IF
                    !IF (phic < phint) THEN
                    ! phint = phic
                    !END IF
                    phint = 1
                    IF ((hydrau%reach(i)%kscH + c(i, 5, 2)) <= 0) THEN
                        phic = 0
                    ELSE
                        phic = c(i, 5, 2) / (hydrau%reach(i)%kscH + c(i, 5, 2))
                    END IF
                    IF (phic < phint) THEN
                        phint = phic
                    END IF

                    !'O2 inhibition of heterotroph growth
                    IF (fcarbH < phint) THEN
                        phint = fcarbH
                    END IF
                    !'convert limited zero-order enhanced hyporheic oxygen flux
                    !'or first-order enhanced hyporheic fast CBOD oxidation
                    !'to equivalent units of gD/d of net growth of heterotrophic bacteria biofilm
                    IF (Rates%typeH == "First-order") THEN
                        !'convert limited first-order d^-1 to gD/d of net growth
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 5, 2)
                    ELSE
                        !'convert limited zero-order gO2/m^2/d to gD/d of net growth
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * hydrau%reach(i)%Ast
                    END IF

                    !'detritus (12)
                    DetrDiss = kdtT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 12, 2)
                    dc(i, 12, 2) = dc(i, 12, 2) - DetrDiss
                    dc(i, 12, 2) = dc(i, 12, 2) + Rates%ada * PhytoDeath

                    !'Organic Nitrogen (6)
                    OrgNHydr = khnT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 6, 2)
                    dc(i, 6, 2) = dc(i, 6, 2) + Rates%ana * PhytoDeath
                    dc(i, 6, 2) = dc(i, 6, 2) - OrgNHydr

                    !'Ammonium Nitrogen (7)

                    !gp 03-Dec-09
                    !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 7, 2)
                    NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * Fi * c(i, 7, 2)

                    dc(i, 7, 2) = dc(i, 7, 2) + OrgNHydr
                    dc(i, 7, 2) = dc(i, 7, 2) - NH4Nitrif
                    dc(i, 7, 2) = dc(i, 7, 2) + Rates%ana * PhytoResp

                    !'Nitrate Nitrogen (8)
                    Denitr = c(i, 5, 2) / (0.1 + c(i, 5, 2)) * fdenitr * kiT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 8, 2)
                    dc(i, 8, 2) = dc(i, 8, 2) + NH4Nitrif
                    dc(i, 8, 2) = dc(i, 8, 2) - Denitr

                    !'Organic Phosphorus (9)
                    OrgPHydr = khpT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 9, 2)
                    dc(i, 9, 2) = dc(i, 9, 2) + Rates%apa * PhytoDeath
                    dc(i, 9, 2) = dc(i, 9, 2) - OrgPHydr

                    !'Inorganic Phosphorus (10)
                    dc(i, 10, 2) = dc(i, 10, 2) + OrgPHydr
                    dc(i, 10, 2) = dc(i, 10, 2) + Rates%apa * PhytoResp

                    !'CBOD Slow (4)
                    CBODsHydr = khcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 4, 2)
                    dc(i, 4, 2) = dc(i, 4, 2) + Rates%roc * (1 / Rates%adc) * DetrDiss
                    dc(i, 4, 2) = dc(i, 4, 2) - CBODsHydr
                    CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 4, 2)
                    dc(i, 4, 2) = dc(i, 4, 2) - CBODsOxid

                    !'CBOD Fast (5)
                    CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 5, 2)
                    CBODfOxid = CBODfOxid + Rates%roc * (1 / Rates%adc) * HeteroGrow !'enhanced CBOD oxidation from hyporheic biofilm growth
                    dc(i, 5, 2) = dc(i, 5, 2) + CBODsHydr
                    dc(i, 5, 2) = dc(i, 5, 2) - CBODfOxid
                    dc(i, 5, 2) = dc(i, 5, 2) - Rates%rondn * Denitr

                    !gp 08-Dec-04
                    !GENERIC CONSTITUENT or COD (14)
                    dc(i, 14, 2) = dc(i, 14, 2) - kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 14, 2)
                    IF (Rates%useGenericAsCOD == "Yes") THEN
                        CODoxid = kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 14, 2)
                    ELSE
                        CODoxid = 0
                    END IF

                    !'Dissolved Oxygen (3)
                    dc(i, 3, 2) = dc(i, 3, 2) - CBODsOxid - CBODfOxid
                    dc(i, 3, 2) = dc(i, 3, 2) - Rates%ron * NH4Nitrif
                    dc(i, 3, 2) = dc(i, 3, 2) - Rates%roa * PhytoResp
                    dc(i, 3, 2) = dc(i, 3, 2) - CODoxid !gp 08-Dec-04

                    !'Alkalinity (nv - 2)

                    If (sys%simAlk == "Yes") Then !gp 26-Oct-07

                        !gp 03-Dec-09
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkaa * PhytoResp * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkn * NH4Nitrif * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkden * Denitr * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * OrgNHydr * 50000 + Rates%ralkbp * OrgPHydr * 50000
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Fi * Rates%ana * PhytoResp * 50043.45_r64 !'phyto resp of N
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbp * Pcharge * Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbn * 2.0_r64 * NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss (Fi) and nitrate gain (+1)
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Denitr * 50043.45_r64 !'denitrification
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Fi * OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbp * Pcharge * OrgPHydr * 50043.45_r64 !'organic P hydrolysis

                    end if !gp 26-Oct-07

                    !'Inorganic Carbon (nv - 1)
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcca * PhytoResp
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcco * CBODfOxid
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcco * CBODsOxid

                    !'PATHOGEN
                    dc(i, 13, 2) = dc(i, 13, 2) - kpathT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 13, 2)

                END DO

                !'
                !' ------------------------------------------------------------------------
                !' --- convert hyporheic pore water deriviatives to concentration units ---
                !' ------------------------------------------------------------------------
                !'

                DO i = 1, nr
                    DO k = 1, nv - 1
                        dc(i, k, 2) = dc(i, k, 2) / hydrau%reach(i)%HypoPoreVol
                    END DO
                END DO


                !gp 15-Nov-04
                !
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                ! --- Level 2 growth kinetics for hyporheic biofilm ---
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                !

              CASE ('Level 2') !gp 15-Nov-04

                !
                ! ------------------------
                ! --- calc and save pH ---
                ! ------------------------
                !

                DO i=0, nr

                    !gp 23-Nov-09
                    !IF (sys%IMethpH == "Newton-Raphson") THEN
                    ! CALL phsolNewton(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    !ELSE
                    ! CALL phsolBisect(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    !END IF
                    IF (sys%IMethpH == "Newton-Raphson") THEN
                        CALL phsolNewton(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    ELSEIF (sys%IMethpH == "Bisection") THEN
                        CALL phsolBisect(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    ELSE
                        CALL phsolBrent(pH, c(i, nv - 1, 2), Te(i, 2), c(i, nv - 2, 2), c(i, 1, 2))
                    END IF

                    pHs(i, 2) = pH
                    CALL ChemRates(Te(i, 2), K1, K2, KW, Kh, c(i, 1, 2))
                    K1s(i, 2) = K1; K2s(i, 2) = K2; Khs(i, 2) = Kh
                END DO

                !
                ! -----------------------------------------------------
                ! --- temperature adjustment of rates and constants ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !gp 08-Feb-06
                    !khcT(i) = Rates%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    !kdcsT(i) = Rates%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    !kdcT(i) = Rates%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    !khnT(i) = Rates%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    !khpT(i) = Rates%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    !knT(i) = Rates%kn * Rates%tkn ** (Te(i, 2) - 20)
                    !kiT(i) = Rates%ki * Rates%tki ** (Te(i, 2) - 20)
                    !os(i, 2)=oxsat(Te(i, 2), hydrau%reach(i)%elev)
                    !kdeaT(i) = Rates%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    !kreaT(i) = Rates%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    !kdtT(i) = Rates%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    !kpathT(i) = Rates%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !enhanced CBOD oxidation rate from hyporheic biofilm
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm growth
                    !kreaHT(i) = Rates%kreaH * Rates%tkreaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm respiration
                    !kdeaHT(i) = Rates%kdeaH * Rates%tkdeaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm death
                    !kgenT(i) = Rates%kgen * Rates%tkgen ** (Te(i, 2) - 20) !gp 30-Nov-04
                    khcT(i) = hydrau%reach(i)%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    kdcsT(i) = hydrau%reach(i)%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    kdcT(i) = hydrau%reach(i)%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    khnT(i) = hydrau%reach(i)%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    khpT(i) = hydrau%reach(i)%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    knT(i) = hydrau%reach(i)%kn * Rates%tkn ** (Te(i, 2) - 20)
                    kiT(i) = hydrau%reach(i)%ki * Rates%tki ** (Te(i, 2) - 20)
                    os(i, 2)=oxygen_saturation(Te(i, 2), hydrau%reach(i)%elev)
                    kdeaT(i) = hydrau%reach(i)%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    kreaT(i) = hydrau%reach(i)%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    kdtT(i) = hydrau%reach(i)%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    kpathT(i) = hydrau%reach(i)%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    kgaHT(i) = hydrau%reach(i)%kgaH * Rates%tkgaH ** (Te(i, 2) - 20)
                    kreaHT(i) = hydrau%reach(i)%kreaH * Rates%tkreaH ** (Te(i, 2) - 20)
                    kdeaHT(i) = hydrau%reach(i)%kdeaH * Rates%tkdeaH ** (Te(i, 2) - 20)
                    kgenT(i) = hydrau%reach(i)%kgen * Rates%tkgen ** (Te(i, 2) - 20)

                END DO

                !
                ! ---------------------------------------------------------------------------------------------
                ! --- mass transfer between water column and hyporheic zone by bulk hyporheic exchange flow ---
                ! ---------------------------------------------------------------------------------------------
                !

                DO i = 1, nr
                    DO k = 1, nv-1
                        !gp un-comment next two lines to activate mass transfer between water column and hyporheic pore water
                        dc(i, k, 1) = dc(i, k, 1) + hydrau%reach(i)%EhyporheicCMD * (c(i, k, 2) - c(i, k, 1))
                        dc(i, k, 2) = hydrau%reach(i)%EhyporheicCMD * (c(i, k, 1) - c(i, k, 2))
                    END DO
                    !store hyporheic fluxes for output (positive flux is source to the water column from the hyporheic zone)
                    HypoFluxDO(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 3, 2) - c(i, 3, 1)) / hydrau%reach(i)%Ast !DO gO2/m^2/d
                    HypoFluxCBOD(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 5, 2) - c(i, 5, 1)) / hydrau%reach(i)%Ast !fast CBOD gO2/m^2/d
                    HypoFluxNH4(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 7, 2) - c(i, 7, 1)) / hydrau%reach(i)%Ast !NH4 mgN/m^2/d
                    HypoFluxNO3(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 8, 2) - c(i, 8, 1)) / hydrau%reach(i)%Ast !NO3 mgN/m^2/d
                    HypoFluxSRP(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, 10, 2) - c(i, 10, 1)) / hydrau%reach(i)%Ast !SRP mgP/m^2/d
                    HypoFluxIC(i) = hydrau%reach(i)%EhyporheicCMD * (c(i, nv - 1, 2) - c(i, nv - 1, 1)) &
                        / hydrau%reach(i)%Ast / Rates%rccc !cT gC/m^2/d
                END DO

                !
                ! -----------------------------------------------------
                ! --- mass derivatives for the hyporheic pore water ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    ! --- pH dependent varialbes for later use in deriv calcs in this i loop

                    !'gp 03-Dec-09
                    !'fraction of ionized ammonia
                    hh = 10.0_r64 ** (-pHs(i, 2))
                    Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 2) + 273.15_r64)))
                    Fi = hh / (hh + Kamm)
                    !'phosphate fractions of H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                    KPO41 = 10.0_r64 ** (-2.15_r64)
                    KPO42 = 10.0_r64 ** (-7.2_r64)
                    KPO43 = 10.0_r64 ** (-12.35_r64)
                    DPO = 1.0_r64 / (hh ** 3.0_r64 + KPO41 * hh ** 2.0_r64 + KPO41 * KPO42 * hh + KPO41 * KPO42 * KPO43)
                    FPO41 = KPO41 * hh ** 2.0_r64 * DPO !'fraction of phosphate as H2PO4-
                    FPO42 = KPO41 * KPO42 * hh * DPO !'fraction of phosphate as HPO4--
                    FPO43 = KPO41 * KPO42 * KPO43 * DPO !'fraction of phosphate as PO4---
                    !'stoichiometric conversion factors (some are pH dependent)
                    !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqH+/L per mgN/m^3
                    !ralkbp = (FPO41 + 2# * FPO42 + 3# * FPO43) / 30.973762 / 1000# / 1000# !'eqH+/L per mgP/m^3
                    Pcharge = (FPO41 + 2.0_r64 * FPO42 + 3.0_r64 * FPO43) !avg charge of P species for multiplier for ralkbp

                    SELECT CASE (Rates%IkoxC) !'low O2 inhibition of C oxidation by floating heterotrophs
                      CASE (1)
                        fcarb = c(i, 3, 2) / (Rates%Ksocf + c(i, 3, 2))
                      CASE (2)
                        fcarb = 1 - Exp(-Rates%Ksocf * c(i, 3, 2))
                      CASE (3)
                        fcarb = c(i, 3, 2) ** 2 / (Rates%Ksocf + c(i, 3, 2) ** 2)
                    END SELECT

                    !gp 08-Feb-06
                    !SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                    !CASE (1)
                    ! fcarbH = c(i, 3, 2) / (Rates%kinhcH + c(i, 3, 2))
                    !CASE (2)
                    ! fcarbH = 1 - Exp(-Rates%kinhcH * c(i, 3, 2))
                    !CASE (3)
                    ! fcarbH = c(i, 3, 2) ** 2 / (Rates%kinhcH + c(i, 3, 2) ** 2)
                    !END SELECT
                    SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                      CASE (1)
                        fcarbH = c(i, 3, 2) / (hydrau%reach(i)%kinhcH + c(i, 3, 2))
                      CASE (2)
                        fcarbH = 1 - Exp(-hydrau%reach(i)%kinhcH * c(i, 3, 2))
                      CASE (3)
                        fcarbH = c(i, 3, 2) ** 2 / (hydrau%reach(i)%kinhcH + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxN) !'low O2 inhibition of nitrification
                      CASE (1)
                        fnitr = c(i, 3, 2) / (Rates%Ksona + c(i, 3, 2))
                      CASE (2)
                        fnitr = 1 - Exp(-Rates%Ksona * c(i, 3, 2))
                      CASE (3)
                        fnitr = c(i, 3, 2) ** 2 / (Rates%Ksona + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxDN) !'low O2 enhancement of denitrification
                      CASE (1)
                        fdenitr = 1 - c(i, 3, 2) / (Rates%Ksodn + c(i, 3, 2))
                      CASE (2)
                        fdenitr = Exp(-Rates%Ksodn * c(i, 3, 2))
                      CASE (3)
                        fdenitr = 1 - c(i, 3, 2) ** 2 / (Rates%Ksodn + c(i, 3, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxP) !'low O2 inhib of phytoplankton respiration
                      CASE (1)
                        frespp = c(i, 3, 2) / (Rates%Ksop + c(i, 3, 2))
                      CASE (2)
                        frespp = 1 - Exp(-Rates%Ksop * c(i, 3, 2))
                      CASE (3)
                        frespp = c(i, 3, 2) ** 2 / (Rates%Ksop + c(i, 3, 2) ** 2)
                    END SELECT

                    !'Phytoplankton (11)
                    PhytoResp = frespp * kreaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 11, 2)
                    PhytoDeath = kdeaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 11, 2)
                    dc(i, 11, 2) = dc(i, 11, 2) - PhytoResp - PhytoDeath

                    !gp 15-Nov-04
                    !Attached heterotrophic bacteria (nv)
                    !nutrient limitation

                    !gp 03-Dec-09
                    !If (hydrau%reach(i)%ksnH + c(i, 7, 2) + c(i, 8, 2) <= 0) Then
                    ! phint = 0
                    !Else
                    ! phint = (c(i, 7, 2) + c(i, 8, 2)) / (hydrau%reach(i)%ksnH + c(i, 7, 2) + c(i, 8, 2))
                    !End If
                    !If (hydrau%reach(i)%kspH + c(i, 10, 2) <= 0) Then
                    ! phip = 0
                    !Else
                    ! phip = c(i, 10, 2) / (hydrau%reach(i)%kspH + c(i, 10, 2))
                    !End If
                    !If (phip < phint) phint = phip
                    If (hydrau%reach(i)%ksnH + Fi * c(i, 7, 2) + c(i, 8, 2) <= 0) Then
                        phint = 0
                    Else
                        phint = (Fi * c(i, 7, 2) + c(i, 8, 2)) / (hydrau%reach(i)%ksnH + Fi * c(i, 7, 2) + c(i, 8, 2))
                    End If
                    If (hydrau%reach(i)%kspH + c(i, 10, 2) <= 0) Then
                        phip = 0
                    Else
                        phip = c(i, 10, 2) / (hydrau%reach(i)%kspH + c(i, 10, 2))
                    End If
                    If (phip < phint) phint = phip

                    !limitation of enhanced oxidation of fast DOC from fast DOC

                    If (hydrau%reach(i)%kscH + c(i, 5, 2) <= 0) Then
                        phic = 0
                    Else
                        phic = c(i, 5, 2) / (hydrau%reach(i)%kscH + c(i, 5, 2))
                    End If

                    If (phic < phint) phint = phic
                    !O2 inhibition of enhanced oxidation of fast DOC
                    If (fcarbH < phint) phint = fcarbH
                    !optional first-order growth kinetics
                    If (Rates%typeH == "First-order") Then

                        !gp 08-Feb-06
                        !HeteroGrow = phint * kgaHT(i) * hydrau%reach(i)%Ast * c(i, nv, 2) * (1 - c(i, nv, 2) / Rates%ahmax) !use first-order growth kinetics
                        HeteroGrow = phint * kgaHT(i) * hydrau%reach(i)%Ast * c(i, nv, 2) &
                            * (1 - c(i, nv, 2) / hydrau%reach(i)%ahmax) !use first-order growth kinetics

                    Else
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * hydrau%reach(i)%Ast !use zero-order, kgaH is gO2/m2/d
                    End If
                    HeteroResp = fcarbH * kreaHT(i) * hydrau%reach(i)%Ast * c(i, nv, 2) !first-order respiration with O2 limitation
                    HeteroDeath = kdeaHT(i) * hydrau%reach(i)%Ast * c(i, nv, 2) !first-order death
                    dc(i, nv, 2) = HeteroGrow - HeteroResp - HeteroDeath !units gD/d deriv for biofilm biomass

                    !ammonium preference

                    !gp 03-Dec-09
                    !prefamH = 0
                    !If (c(i, 7, 2) + c(i, 8, 2) > 0) Then
                    ! prefamH = c(i, 7, 2) * c(i, 8, 2) / (hydrau%reach(i)%khnxH + c(i, 7, 2)) / (hydrau%reach(i)%khnxH + c(i, 8, 2)) &
                    ! + c(i, 7, 2) * hydrau%reach(i)%khnxH / (c(i, 7, 2) + c(i, 8, 2)) / (hydrau%reach(i)%khnxH + c(i, 8, 2))
                    !End If
                    prefamH = 0
                    If (Fi * c(i, 7, 2) + c(i, 8, 2) > 0) Then
                        prefamH = Fi * c(i, 7, 2) * c(i, 8, 2) / (hydrau%reach(i)%khnxH + Fi * c(i, 7, 2)) &
                            / (hydrau%reach(i)%khnxH + c(i, 8, 2)) &
                            + Fi * c(i, 7, 2) * hydrau%reach(i)%khnxH / (Fi * c(i, 7, 2) + c(i, 8, 2)) &
                            / (hydrau%reach(i)%khnxH + c(i, 8, 2))
                    End If

                    !'detritus (12)
                    DetrDiss = kdtT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 12, 2)
                    dc(i, 12, 2) = dc(i, 12, 2) - DetrDiss
                    dc(i, 12, 2) = dc(i, 12, 2) + Rates%ada * PhytoDeath
                    dc(i, 12, 2) = dc(i, 12, 2) + HeteroDeath !gp 22-Nov-04

                    !'Organic Nitrogen (6)
                    OrgNHydr = khnT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 6, 2)
                    dc(i, 6, 2) = dc(i, 6, 2) + Rates%ana * PhytoDeath
                    dc(i, 6, 2) = dc(i, 6, 2) - OrgNHydr
                    dc(i, 6, 2) = dc(i, 6, 2) + HeteroDeath * Rates%ana / Rates%ada !gp 15-Nov-04

                    !'Ammonium Nitrogen (7)

                    !gp 03-Dec-09
                    !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 7, 2)
                    NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * Fi * c(i, 7, 2)

                    dc(i, 7, 2) = dc(i, 7, 2) + OrgNHydr
                    dc(i, 7, 2) = dc(i, 7, 2) - NH4Nitrif
                    dc(i, 7, 2) = dc(i, 7, 2) + Rates%ana * PhytoResp
                    dc(i, 7, 2) = dc(i, 7, 2) + Rates%anc / Rates%adc * HeteroResp !gp 15-Nov-04
                    dc(i, 7, 2) = dc(i, 7, 2) - Rates%anc / Rates%adc * prefamH * HeteroGrow !gp 15-Nov-04

                    !'Nitrate Nitrogen (8)
                    Denitr = c(i, 5, 2) / (0.1 + c(i, 5, 2)) * fdenitr * kiT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 8, 2)
                    dc(i, 8, 2) = dc(i, 8, 2) + NH4Nitrif
                    dc(i, 8, 2) = dc(i, 8, 2) - Denitr
                    dc(i, 8, 2) = dc(i, 8, 2) - Rates%anc / Rates%adc * (1 - prefamH) * HeteroGrow !gp 15-Nov-04

                    !'Organic Phosphorus (9)
                    OrgPHydr = khpT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 9, 2)
                    dc(i, 9, 2) = dc(i, 9, 2) + Rates%apa * PhytoDeath
                    dc(i, 9, 2) = dc(i, 9, 2) - OrgPHydr
                    dc(i, 9, 2) = dc(i, 9, 2) + HeteroDeath * Rates%apa / Rates%ada !gp 15-Nov-04

                    !'Inorganic Phosphorus (10)
                    dc(i, 10, 2) = dc(i, 10, 2) + OrgPHydr
                    dc(i, 10, 2) = dc(i, 10, 2) + Rates%apa * PhytoResp
                    dc(i, 10, 2) = dc(i, 10, 2) + Rates%apc / Rates%adc * HeteroResp !gp 15-Nov-04
                    dc(i, 10, 2) = dc(i, 10, 2) - Rates%apc / Rates%adc * HeteroGrow !gp 15-Nov-04

                    !'CBOD Slow (4)
                    CBODsHydr = khcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 4, 2)
                    dc(i, 4, 2) = dc(i, 4, 2) + Rates%roc * (1 / Rates%adc) * DetrDiss
                    dc(i, 4, 2) = dc(i, 4, 2) - CBODsHydr
                    CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 4, 2)
                    dc(i, 4, 2) = dc(i, 4, 2) - CBODsOxid

                    !'CBOD Fast (5)
                    CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 5, 2) !oxidation by suspended bacteria
                    CBODfOxid = CBODfOxid + Rates%roc * (1 / Rates%adc) * HeteroGrow !enhanced CBOD oxidation from hyporheic biofilm growth
                    dc(i, 5, 2) = dc(i, 5, 2) + CBODsHydr
                    dc(i, 5, 2) = dc(i, 5, 2) - CBODfOxid
                    dc(i, 5, 2) = dc(i, 5, 2) - Rates%rondn * Denitr

                    !gp 08-Dec-04
                    !GENERIC CONSTITUENT or COD (14)
                    dc(i, 14, 2) = dc(i, 14, 2) - kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 14, 2)
                    IF (Rates%useGenericAsCOD == "Yes") THEN
                        CODoxid = kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 14, 2)
                    ELSE
                        CODoxid = 0
                    END IF

                    !'Dissolved Oxygen (3)
                    dc(i, 3, 2) = dc(i, 3, 2) - CBODsOxid - CBODfOxid !CBODfOxid includes oxidation by hyporheic biofilm
                    dc(i, 3, 2) = dc(i, 3, 2) - Rates%ron * NH4Nitrif
                    dc(i, 3, 2) = dc(i, 3, 2) - Rates%roa * PhytoResp
                    dc(i, 3, 2) = dc(i, 3, 2) - Rates%roc / Rates%adc * HeteroResp !gp 15-Nov-04
                    dc(i, 3, 2) = dc(i, 3, 2) - CODoxid !gp 08-Dec-04

                    !'Alkalinity (nv - 2)

                    If (sys%simAlk == "Yes") Then !gp 26-Oct-07

                        !gp 03-Dec-09
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkaa * PhytoResp * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkn * NH4Nitrif * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkden * Denitr * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * OrgNHydr * 50000 + Rates%ralkbp * OrgPHydr * 50000
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkda * HeteroGrow * prefamH * 50000 !gp 15-Nov-04
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkdn * HeteroGrow * (1 - prefamH) * 50000 !gp 15-Nov-04
                        !dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkda * HeteroResp * 50000 !gp 15-Nov-04
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Fi * Rates%ana * PhytoResp * 50043.45_r64 !'phyto resp of N
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbp * Pcharge * Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbn * 2.0_r64 * NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss (Fi) and nitrate gain (+1)
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Denitr * 50043.45_r64 !'denitrification
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Fi * OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbp * Pcharge * OrgPHydr * 50043.45_r64 !'organic P hydrolysis
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbn * Rates%anc / Rates%adc * prefamH &
                            * HeteroGrow * 50043.45_r64 !'hetero uptake of ammonia
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Rates%anc / Rates%adc * (1 - prefamH) &
                            * HeteroGrow * 50043.45_r64 !'hetero uptake of nitrate
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbp * Pcharge * Rates%apc / Rates%adc &
                            * HeteroGrow * 50043.45_r64 !'hetero uptake of SRP
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) + Rates%ralkbn * Fi * Rates%anc / Rates%adc &
                            * HeteroResp * 50043.45_r64 !'hetero resp of ammonia
                        dc(i, nv - 2, 2) = dc(i, nv - 2, 2) - Rates%ralkbp * Pcharge * Rates%apc / Rates%adc &
                            * HeteroResp * 50043.45_r64 !'hetero resp of SRP

                    end if !gp 26-Oct-07

                    !'Inorganic Carbon (nv - 1)
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcca * PhytoResp
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcco * CBODfOxid
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rcco * CBODsOxid
                    dc(i, nv - 1, 2) = dc(i, nv - 1, 2) + Rates%rccd * HeteroResp !gp 15-Nov-04

                    !'PATHOGEN
                    dc(i, 13, 2) = dc(i, 13, 2) - kpathT(i) * hydrau%reach(i)%HypoPoreVol * c(i, 13, 2)

                END DO

                !'
                !' ------------------------------------------------------------------------
                !' --- convert hyporheic pore water deriviatives to concentration units ---
                !' ------------------------------------------------------------------------
                !'

                DO i = 1, nr
                    DO k = 1, nv - 1
                        dc(i, k, 2) = dc(i, k, 2) / hydrau%reach(i)%HypoPoreVol
                    END DO
                    dc(i, nv, 2) = dc(i, nv, 2) / hydrau%reach(i)%Ast !gp 15-Nov-04
                END DO


            END SELECT !'gp end of calculation of hyporheic water quality derivatives


            !'
            !' -------------------------------------------------------------------------
            !' --- convert water column derivatives from mass to concentration units ---
            !' -------------------------------------------------------------------------
            !'

            DO i = 1, nr
                DO k = 1, nv - 1
                    !gp dc(i, k) = dc(i, k) / hydrau%reach(i)%vol
                    dc(i, k, 1) = dc(i, k, 1) / hydrau%reach(i)%vol !gp
                END DO
                !gp dc(i, nv) = dc(i, nv) / hydrau%reach(i)%Asb
                dc(i, nv, 1) = dc(i, nv, 1) / hydrau%reach(i)%Asb !gp
                dINb(i) = dINb(i) / hydrau%reach(i)%Asb
                dIPb(i) = dIPb(i) / hydrau%reach(i)%Asb
            END DO



!gp 03-Feb-05 end of bypass derivs for water quality variables
!unless 'All' state variables are being simulated
!
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

        End If

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




    END SUBROUTINE

!----------------------------------------------------------------------------------------------

    !gp 25-Aug-08
    !SUBROUTINE SedCalcNumNew(Jcin, Jnin, Jpin, O20, depth, Tw, SOD, JNH4, JNO3, JCH4, &
    ! JCH4g, JPO4, NH30, NO30, PO40, CH40, CSOD)
    SUBROUTINE SedCalcNumNew(Jcin, Jnin, Jpin, O20, depth, Tw, SOD, JNH4, JNO3, JCH4, &
        JCH4g, JPO4, NH30, NO30, PO40, CH40, CSOD, calcSedFlux)

        REAL(r64), PARAMETER :: KappaNH3 = 0.131_r64, ThtaNH3 = 1.123_r64, KM_NH3 = 0.728_r64, KM_O2_NH3 = 0.37_r64
        REAL(r64), PARAMETER :: KappaNO3_1 = 0.1_r64, KappaNO3_2 = 0.25_r64, ThtaNO3 = 1.08_r64
        REAL(r64), PARAMETER :: KappaCH4 = 0.7_r64, ThtaCH4 = 1.079_r64
        REAL(r64), PARAMETER :: Dd = 0.001_r64, ThtaDd = 1.08_r64, Dp0 = 0.00012_r64, ThtaDp = 1.117_r64
        REAL(r64), PARAMETER :: KdNH3 = 1.0_r64
        REAL(r64), PARAMETER :: m1 = 0.5_r64, m2 = 0.5_r64
        REAL(r64), PARAMETER :: w2 = 0.000005_r64, H2 = 0.1_r64
        REAL(r64), PARAMETER :: dKDPO41 = 20.0_r64
        REAL(r64), PARAMETER :: KdPO42 = 20.0_r64
        REAL(r64), PARAMETER :: O2critPO4 = 2.0_r64


        REAL(r64), INTENT(IN) :: Jcin, Jnin, Jpin, O20, depth, Tw, CH40, PO40, NH30, NO30

        !gp 25-Aug-08
        CHARACTER(LEN=30), INTENT(IN) :: calcSedFlux

        REAL(r64), INTENT(OUT) :: SOD, JNH4, JNO3, JCH4, JCH4g, JPO4, CSOD
        INTEGER(i32) i, maxit, iter
        REAL(r64) es, ea
        REAL(r64) NSOD, SODold

        REAL(r64) w12, KL12
        REAL(r64) s
        REAL(r64) NH3(2), NH3T(2)
        REAL(r64) CH4(2)
        REAL(r64) NH3toNO3, CH4toCO2
        REAL(r64) JNH3toNO3
        REAL(r64) Denit(2)
        REAL(r64) NO3(2)
        REAL(r64) JDenit(2), JDenitT
        REAL(r64) JO2NO3(2), JO2NO3T, JC_O2equiv
        REAL(r64) CH4SAT, CSODmax
        REAL(r64) SECH_ARG
        REAL(r64) fd1, fp1, fd2, fp2
        REAL(r64) a11, a12, a21, a22, b1, b2

        REAL(r64) fpon(3), fpoc(3), fpop(3)
        REAL(r64) kdiaPOC(3), ThtaPOC(3)
        REAL(r64) kdiaPON(3), ThtaPON(3)
        REAL(r64) kdiaPOP(3), ThtaPOP(3)
        REAL(r64) JPOC(3), JPON(3), JPOP(3)
        REAL(r64) POC2(3), PON2(3), POP2(3)
        REAL(r64) POCT2, PONT2, POPT2
        REAL(r64) Jc, Jn, Jp

        !Phosphorus
        REAL(r64) KdPO41, PO4T(2), PO4(2)

        NH3=0; NH3T=0; CH4=0
        Denit=0; NO3=0; JDenit=0; JO2NO3=0
        KdPO41=0; PO4T=0; PO4=0
        !Compute influxes corrected for settling and refractory
        fpon(1) = 0.65_r64; fpon(2) = 0.25_r64; fpon(3) = 1.0_r64 - fpon(1) - fpon(2)
        fpoc(1) = 0.65_r64; fpoc(2) = 0.2_r64; fpoc(3) = 1.0_r64 - fpoc(1) - fpoc(2)
        fpop(1) = 0.65_r64; fpop(2) = 0.2_r64; fpop(3) = 1.0_r64 - fpop(1) - fpop(2)

        kdiaPON(1) = 0.035_r64; ThtaPON(1) = 1.1_r64
        kdiaPON(2) = 0.0018_r64; ThtaPON(2) = 1.15_r64
        kdiaPON(3) = 0; ThtaPON(3) = 1.17_r64

        kdiaPOC(1) = 0.035_r64; ThtaPOC(1) = 1.1_r64
        kdiaPOC(2) = 0.0018_r64; ThtaPOC(2) = 1.15_r64
        kdiaPOC(3) = 0; ThtaPOC(3) = 1.17_r64

        kdiaPOP(1) = 0.035_r64; ThtaPOP(1) = 1.1_r64
        kdiaPOP(2) = 0.0018_r64; ThtaPOP(2) = 1.15_r64
        kdiaPOP(3) = 0; ThtaPOP(3) = 1.17_r64

        !compute input fluxes
        DO i = 1, 3
            JPOC(i) = Jcin * fpoc(i)
            JPON(i) = Jnin * fpon(i)
            JPOP(i) = Jpin * fpop(i)
        END DO

        !compute particulate organic forms
        DO i = 1, 3
            POC2(i) = JPOC(i) / H2 / (kdiaPOC(i) * ThtaPOC(i) ** (Tw - 20.0_r64) + w2 / H2)
            PON2(i) = JPON(i) / H2 / (kdiaPON(i) * ThtaPON(i) ** (Tw - 20.0_r64) + w2 / H2)
            POP2(i) = JPOP(i) / H2 / (kdiaPOP(i) * ThtaPOP(i) ** (Tw - 20.0_r64) + w2 / H2)
            POCT2 = POCT2 + POC2(i)
            PONT2 = PONT2 + PON2(i)
            POPT2 = POPT2 + POP2(i)
        END DO

        !compute diagenesis fluxes
        Jc = 0; Jn = 0; Jp = 0
        DO i = 1, 3
            Jc = Jc + kdiaPOC(i) * ThtaPOC(i) ** (Tw - 20) * POC2(i) * H2
            Jn = Jn + kdiaPON(i) * ThtaPON(i) ** (Tw - 20) * PON2(i) * H2
            Jp = Jp + kdiaPOP(i) * ThtaPOP(i) ** (Tw - 20) * POP2(i) * H2
        END DO

        maxit = 500
        es = 0.1_r64
        w12 = Dp0 * ThtaDp ** (Tw - 20.0_r64) / (H2 / 2.0_r64)
        KL12 = Dd * ThtaDd ** (Tw - 20.0_r64) / (H2 / 2.0_r64)

        SODold = Jc + 1.714_r64 * Jn
        iter = 0
        ! Saturation conc. of methane in oxygen equivalent units (Equation 10.51)
        CH4SAT = 100.0_r64 * (1.0_r64 + depth / 10.0_r64) * (1.024_r64 ** (20.0_r64 - Tw)) ![gmO*/m3]
        IF ((O20 < 0.001_r64) .OR. (Jc + Jn <= 0)) THEN
            SOD = 0
            JNH4 = Jn
            JCH4 = MIN(SQRT(2.0_r64 * KL12 * CH4SAT * Jc), Jc)
            JCH4g = Jc - JCH4
            JPO4 = Jp
        ELSE
            DO
                s = SODold / O20
                NH3toNO3 = KappaNH3 ** 2 * ThtaNH3 ** (Tw - 20.0_r64) / s * KM_NH3 / (KM_NH3 + NH3(1)) * &
                    O20 / (2.0_r64 * KM_O2_NH3 + O20)
                ! [m/d] = [m2/d2] / [m/d] * [-] * [-]

                ! Calculate dissolved and particulate (sorbed) fractions
                fd1 = (1.0_r64 / (1.0_r64 + m1 * KdNH3))
                fp1 = 1.0_r64 - fd1 != ((m1*KdNH3)/(1 + m1*KdNH3))
                fd2 = (1.0_r64 / (1.0_r64 + m2 * KdNH3))
                fp2 = 1.0_r64 - fd2 != ((m2*KdNH3)/(1 + m2*KdNH3))

                ! Write linear sys of equations around NH3T
                a11 = -fd1 * KL12 - fp1 * w12 - fd1 * NH3toNO3 - fd1 * s - w2
                a12 = fd2 * KL12 + fp2 * w12
                b1 = -s * NH30 ![m/d]*[mg/m3]
                a21 = fd1 * KL12 + fp1 * w12 + w2
                a22 = -fd2 * KL12 - fp2 * w12 - w2
                b2 = -Jn

                CALL Lin_Sys(a11, a12, a21, a22, b1, b2, NH3T(1), NH3T(2))

                ! Dissolved Concentrations
                NH3(1) = fd1 * NH3T(1)
                NH3(2) = fd2 * NH3T(2)

                ! Oxygen Flux due to NH3->NO2, see Equation 23.3 in Chapra (1997)
                JNH3toNO3 = NH3toNO3 * NH3(1)
                NSOD = 2.0_r64 * (32.0_r64 / 14.0_r64) * NH3toNO3 * NH3(1)
                ! [gmO/m2-d] = [mol O2/mol N]*[gm O2/mol O2]/[gm N/mol N]*
                ! [gm/1000mg] * [m/day] * [mgN/m3]

                !:::::::::::::::::::::::::::: BEGIN Nitrate::::::::::::::::::::::::::::

                ! Denitrification in layers 1 and 2 (Equation 4.55)
                Denit(1) = (KappaNO3_1 ** 2 * ThtaNO3 ** (Tw - 20.0_r64) / s)
                Denit(2) = KappaNO3_2 * ThtaNO3 ** (Tw - 20.0_r64)

                ! Layer 1
                a11 = -KL12 - Denit(1) - s - w2
                a12 = KL12
                b1 = -s * NO30 - NH3toNO3 * NH3(1)
                ! Layer 2
                a21 = KL12 + w2
                a22 = -KL12 - Denit(2) - w2
                b2 = 0.0

                CALL Lin_Sys(a11, a12, a21, a22, b1, b2, NO3(1), NO3(2))

                ! Nitrate Flux to water column
                JNO3 = s * (NO3(1) - NO30)

                ! Denitrification Flux [mgN/m2d]
                JDenit(1) = Denit(1) * NO3(1)
                JDenit(2) = Denit(2) * NO3(2)
                JDenitT = JDenit(1) + JDenit(2)

                ! Methane consumption due to denitrification (Equation 9.16)
                ! Layer 1
                JO2NO3(1) = 2.0_r64 * (16.0_r64 / 12.0_r64) * (10.0_r64 / 8.0_r64) * (12.0_r64 / 14.0_r64) * JDenit(1)
                ! [gmO*/m2-d] = [molO/molC]*[(gmO/molO)/(gmC/molC)]*[molC/MolN]*
                ! [(gmC/molC)/(gmN/molN)] * [mg/m2d] * [gm/1000mg]
                ! where 2*16/12 is the ubiquitous 32/12 (= 2.67) for oxidation of carbon
                ! Layer 2
                JO2NO3(2) = (32.0_r64 / 12.0_r64) * (10.0_r64 / 8.0_r64) * (12.0_r64 / 14.0_r64) * JDenit(2)
                ! Sum
                JO2NO3T = JO2NO3(1) + JO2NO3(2)

                ! Calculate methane flux in oxygen equivalent units, adjusted for
                ! the methane consumed in denitrification
                JC_O2equiv = Jc - JO2NO3T
                IF (JC_O2equiv < 0) JC_O2equiv = 0.0

                !***** Methane in O2 equivalents
                !:::::::::::::::::::::::::::: BEGIN methane::::::::::::::::::::::::::::


                ! CSODMAX Equations 10.28 and 10.30
                CSODmax = MIN(SQRT(2.0_r64 * KL12 * CH4SAT * JC_O2equiv), JC_O2equiv) ! [gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])
                !MsgBox CSODmax
                !SECH_ARG = (KappaCH4 * ThtaCH4 ** ((Tw - 20) / 2.0)) / s
                ! CSOD Equation 10.35
                ! The hyperbolic secant is defined as HSec(X) = 2 / (EXP(X) + EXP(-X))
                !IF (SECH_ARG < 400.0) THEN !This is the usual case
                ! CSOD = CSODmax * (1.0 - (2.0 / (EXP(SECH_ARG) + EXP(-SECH_ARG))))
                !ELSE !HSec(SECH_ARG) < 3.8E-174 ~ 0
                ! CSOD = CSODmax
                !END IF
                ! Aqueous methane flux to water column
                !JCH4 = CSODmax - CSOD
                ! Gaseous methane flux to water column
                !JCH4g = JC_O2equiv - JCH4 - CSOD

                CH4toCO2 = (KappaCH4 * KappaCH4 * ThtaCH4 ** ((Tw - 20.0_r64) / 2.0_r64)) / s
                CH4(1) = (CSODmax + s * CH40) / (CH4toCO2 + s)
                CSOD = CH4toCO2 * CH4(1)

                SOD = (SODold + CSOD + NSOD) / 2.0_r64
                iter = iter + 1
                ea = Abs((SOD - SODold) / SOD) * 100.0_r64

                !gp 20-May-09
                !IF (ea <= es) THEN
                IF (ea <= es .or. (calcSedFlux == "Option 2" .AND. iter >= maxit)) THEN

                    EXIT


                    !gp 20-May-09
                    !!gp 25-Aug-08
                    !!ELSEIF (iter >= maxit) THEN
                    !! WRITE(*,*) "SOD iterations exceeded! Turn off sediment diagenesis on QUAL2K sheet."
                    !ELSEIF ((calcSedFlux == "Yes" .or. calcSedFlux == "Option 1") .AND. iter >= maxit) THEN
                    ! WRITE(*,*) "SOD iterations exceeded! Turn off sediment diagenesis on QUAL2K sheet."
                ELSEIF (iter >= maxit) THEN
                    WRITE(*,*) "SOD iterations exceeded! Turn off sediment diagenesis on QUAL2K sheet."


                END IF
                SODold = SOD
            END DO
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ! BEGIN ORTHOPHOSPHATE
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            ! Determine Kd in layer 1 (Kd in layer 2 is a constant)
            IF (O20 > O2critPO4) THEN
                KdPO41 = KdPO42 * dKDPO41
            ELSE
                KdPO41 = KdPO42 * (dKDPO41 ** (O20 / O2critPO4))
            END IF

            ! Calculate dissolved and particulate (sorbed) fractions
            fd1 = (1.0_r64 / (1.0_r64 + m1 * KdPO41))
            fp1 = 1.0_r64 - fd1 != ((m1*KdPO41)/(1 + m1*KdPO41))
            fd2 = (1.0_r64 / (1.0_r64 + m2 * KdPO42))
            fp2 = 1.0_r64 - fd2 != ((m2*KdPO42)/(1 + m2*KdPO42))

            ! Write linear sys of equations around PO4T
            ! Layer 1
            a11 = -fd1 * KL12 - fp1 * w12 - fd1 * s - w2
            a12 = fd2 * KL12 + fp2 * w12
            b1 = -s * PO40
            ! Layer 2
            a21 = fd1 * KL12 + fp1 * w12 + w2
            a22 = -fd2 * KL12 - fp2 * w12 - w2
            b2 = -Jp

            CALL Lin_Sys(a11, a12, a21, a22, b1, b2, PO4T(1), PO4T(2))

            ! Dissolved PO4 concentrations
            PO4(1) = fd1 * PO4T(1)
            PO4(2) = fd2 * PO4T(2)

            ! PO4 flux to water column
            JPO4 = s * (PO4(1) - PO40)
            JNH4 = s * (NH3(1) - NH30)
            JNO3 = s * (NO3(1) - NO30)
            JCH4 = s * (CH4(1) - CH40)
            JCH4g = JC_O2equiv - JCH4 - CSOD
        END IF

    END SUBROUTINE


    SUBROUTINE Lin_Sys(a11, a12, a21, a22, b1, b2, x1, x2)

        REAL(r64), INTENT(IN) :: a11, a12, a21, a22, b1, b2
        REAL(r64), INTENT(OUT) :: x1, x2
        !This subroutine solves a linear sys of 2 equations and 2 unknowns

        x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21)
        x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21)

    END SUBROUTINE


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
