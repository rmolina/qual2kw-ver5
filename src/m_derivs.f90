module m_derivs
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use class_hydraulics, only: riverhydraulics_type
    use class_integrationdata, only: &
        saveheatfluxtribs, saveheatfluxadvecdisp, saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, &
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
    use m_sedcalcnew, only: SedCalcNumNew
    implicit none
    private
    public :: derivs

contains

    SUBROUTINE derivs(nr, Meteo, Solar, HW, DB, hydrau, sys, Te, c, INb, IPb, &
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
                            hh = 10.0_r64 ** (-pHs(i-1, 1))
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

                hh = 10.0_r64 ** (-pHs(i, 1))
                alp0 = hh * hh / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1))
                alp1 = K1s(i, 1) * hh / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1)) !fraction of cT as HCO3-
                alp2 = K1s(i, 1) * K2s(i, 1) / (hh ** 2.0_r64 + K1s(i, 1) * hh + K1s(i, 1) * K2s(i, 1)) !fraction of cT as CO3--

                !'gp 03-Dec-09
                !'fraction of total ammonia that is ionized ammonia NH4+ (Fi)
                Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 1) + 273.15_r64))) !'equilibrium coeff for ammonia dissociation NH4+ = NH3(aq) + H+
                Fi = hh / (hh + Kamm) !'fraction of ionized NH4+ = [NH4+ / (NH4+ + NH3)]
                !'fraction of SRP that is H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                KPO41 = 10.0_r64 ** (-2.15_r64) !'= [H+][H2PO4-]/[H3PO4] for eqn H3PO4 = H+ + H2PO4-
                KPO42 = 10.0_r64 ** (-7.2_r64) !'= [H+][HPO4--]/[H2PO4-] for eqn H2PO4- = H+ + HPO4--
                KPO43 = 10.0_r64 ** (-12.35_r64) !'= [H+][PO4---]/[HPO4--] for eqn HPO4-- = H+ + PO4---
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

end module m_derivs
