! classlightheat.f90
! general variables and functions for calculating light and heat
MODULE Class_LightHeat
    USE m_constants
    IMPLICIT NONE

!	PRIVATE
!	PUBLIC Aa, Light_

    TYPE LightHeat_type
        REAL(DP) :: PAR = 0.47  	!phtosynthetically available radiation, default 0.47
        REAL(DP) :: kep = 0.01		!backgroud light extinction
        REAL(DP) :: kela=0		!linear chlorophyll light extinction
        REAL(DP) :: kenla=0		!nonlinear chlorophyll light extinction
        REAL(DP) :: kess=0		!inorganic SS light extinction
        REAL(DP) :: kepom=0		!detritus light extinction

        !gp 13-Feb-06
        REAL(DP) :: kemac=0		!macrophyte light extinction

        CHARACTER(LEN=30) ::longatMethod = "Brunt"
        !atmospheric longwave emissivity model
        !gp 16-Jul-08
        REAL(DP) :: kbrut=1.24

        !gp 24-Jul-09
        REAL(DP) :: KCL1=0.65
        REAL(DP) :: KCL2=0.17

        CHARACTER(LEN=30) ::fUwMethod = "Brady-Graves-Geyer"
        !evaporation and air convection/conduction
    END TYPE

    TYPE(LightHeat_type) lightheat

CONTAINS

    !gp 13-Feb-06
    !PURE FUNCTION lightExtinction(lightheat, cPom, ss, cAlgae) RESULT(ke)
    !
    !	TYPE(LightHeat_type), INTENT(IN) :: lightheat
    !	REAL(DP), INTENT(IN) :: cPom, ss, cAlgae
    !	REAL(DP) ke
    !
    !	ke = lightheat%kep + lightheat%kepom * cPom + lightheat%kess * ss &
    !				+ lightheat%kela * cAlgae + lightheat%kenla * cAlgae ** (2.0 / 3.0)
    !
    !END FUNCTION lightExtinction
    PURE FUNCTION lightExtinction(lightheat, cPom, ss, cAlgae, cMacrophyte) RESULT(ke)
        TYPE(LightHeat_type), INTENT(IN) :: lightheat
        REAL(DP), INTENT(IN) :: cPom, ss, cAlgae, cMacrophyte	!note: cMacrophyte is gD/m^3
        REAL(DP) ke
        ke = lightheat%kep + lightheat%kepom * cPom + lightheat%kess * ss + lightheat%kemac * cMacrophyte &
            + lightheat%kela * cAlgae + lightheat%kenla * cAlgae ** (2.0 / 3.0)
    END FUNCTION lightExtinction

    !/* public functions */

    !gp 13-Feb-06
    !SUBROUTINE Light_(PAR, kep, kela, kenla, kess, kepom, longatMethod, fUwMethod)

    !gp 24-Jun-09
    !gp 16-Jul-08
    !SUBROUTINE Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, fUwMethod)
    !SUBROUTINE Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, kbrut, fUwMethod)
    SUBROUTINE Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, kbrut, fUwMethod, KCL1, KCL2)

        !gp 13-Feb-06
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom

        !gp 24-Jun-09
        !gp 16-Jul-08
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac, kbrut
        REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac, kbrut, KCL1, KCL2

        CHARACTER(LEN=30) longatMethod, fUwMethod

        lightheat%PAR = PAR
        lightheat%kep = kep
        lightheat%kela = kela
        lightheat%kenla=kenla
        lightheat%kess = kess
        lightheat%kepom = kepom

        !gp 13-Feb-06
        lightheat%kemac = kemac

        lightheat%longatMethod = longatMethod

        !gp 16-Jul-08
        lightheat%kbrut = kbrut

        lightheat%fUwMethod= fUwMethod

        !gp 24-Jun-09
        lightheat%KCL1= KCL1
        lightheat%KCL2= KCL2

    END SUBROUTINE Light_

    Function fUw(fUwMethod, Uw, Ta, Te, Ast, eair, es)

        REAL(DP) fUw
        CHARACTER(LEN=20), INTENT(IN) ::fUwMethod
        REAL(DP),INTENT(IN) :: Uw, Ta, Te, Ast, eair, es
        REAL(DP) wmph, areaa
        REAL(DP) ta2a, tva, tsa, tvs, dtv

        Select Case (fUwMethod)
          Case ("Brady-Graves-Geyer")
            !"Brady, Graves, and Geyer"
            !for this formula the wind speed height is 7 m (see Edinger, et al., 1974)
            fUw = 19.0_dp + 0.95_dp * Uw * Uw     !Chapra eqn 30.22 cal/cm**2/d/mmHg
          Case ("Adams 1")
            !write(*,*) "Adams 1, Wrong!"
            !"East Mesa"
            !from adams, et al., eq. 4.48, p. 4-26
            !for this formula the wind speed height is 2m
            !first convert windspeed from m s-1 to mph
            wmph = Uw * 3600.0_dp / (0.3048_dp * 5280.0_dp)     !Uw(i) is at 7m
            !convert to the formula!s wind speed height using
            !the exponential wind law, Paily et al., 1974
            wmph = wmph * (2.0_dp / 7.0_dp) ** 0.15_dp        !convert wind speed from 7m to 2m height
            !convert surface area in m**2 to area in acres
            areaa = Ast / (0.3048_dp ** 2.0 * 43560.0_dp)
            !next compute virtual temperature difference
            !eagleson, p. s.  1970.  dynamic hydrology.  mcgraw-hill, inc.,
            !new york, new york. p. 56
            ta2a = 9.0_dp / 5.0_dp * Ta + 32.0_dp           !deg F air temp at 2m
            tva = ta2a / (1.0_dp - 0.378_dp * (eair / 760.0_dp))
            tsa = 9.0_dp / 5.0_dp * Te + 32.0_dp !Te(i, 1)
            tvs = tsa / (1.0_dp - 0.378_dp * (es / 760.0_dp))
            dtv = tvs - tva     !original formula in shestt
            If (dtv < 0) dtv = 0
            !next compute fUw in W/m**2/mmHg
            fUw = 0.1313_dp * ((22.4_dp * dtv ** (1.0_dp / 3.0_dp)) ** 2 + (24.2_dp * areaa ** (-0.05_dp) * wmph) ** 2) ** 0.5_dp
            !next convert fUw to cal/cm**2/d/mmHg
            fUw = fUw / (4.183076_dp * 100.0_dp * 100.0_dp / 86400.0_dp)
          Case ("Adams 2")
            !write(*,*) "Adams2, Wrong!"
            !!"East Mesa modified to include Marciano and Harbeck"
            !!from adams, et al., eq. 4.48, p. 4-26
            !for this formula the wind speed height is 2m
            !first convert windspeed from m s-1 to mph
            wmph = Uw * 3600.0_dp / (0.3048_dp * 5280.0_dp)     !Uw(i) is at 7m
            !convert to the formula!s wind speed height using
            !the exponential wind law, Paily et al., 1974
            wmph = wmph * (2.0_dp / 7.0_dp) ** 0.15_dp        !convert wind speed from 7m to 2m height
            !next compute virtual temperature difference
            !eagleson, p. s.  1970.  dynamic hydrology.  mcgraw-hill, inc.,
            !new york, new york. p. 56
            ta2a = 9.0_dp / 5.0_dp * Ta + 32.0_dp           !deg F air temp at 2m
            tva = ta2a / (1.0_dp - 0.378_dp * (eair / 760.0_dp))
            tsa = 9.0_dp / 5.0_dp * Te + 32.0_dp			!Te(i, 1)
            tvs = tsa / (1.0_dp - 0.378_dp * (es / 760.0_dp))
            dtv = tvs - tva     !original formula in shestt
            If (dtv < 0) dtv = 0
            !next compute fUw in W/m**2/mmHg
            fUw = 0.1313_dp * ((22.4_dp * dtv ** (1.0_dp / 3.0_dp)) ** 2 + (17.0_dp * wmph) ** 2) ** 0.5_dp
            !next convert fUw to cal/cm**2/d/mmHg
            fUw = fUw / (4.183076_dp * 100.0_dp * 100.0_dp / 86400.0_dp)
          Case Default
            !write(*,*) "Default, Wrong!"
        End Select

    End Function

!evaperation, convection, back radiation, air longat

    !gp 16-Jul-08
    !SUBROUTINE heatBalance(evap, conv, back, longat, lightheat, Te, Ta, Td, cc, Uw, Ast)
    SUBROUTINE heatBalance(evap, conv, back, longat, lightheat, Te, Ta, Td, cc, Uw, Ast, SKOP)

        TYPE(LightHeat_type), INTENT(IN) :: lightheat

        !gp 16-Jul-08
        !REAL(DP), INTENT(IN) :: Te, Ta, Td, cc, Uw, Ast
        REAL(DP), INTENT(IN) :: Te, Ta, Td, cc, Uw, Ast, SKOP

        REAL(DP), INTENT(OUT) :: evap, conv, back, longat
        REAL(DP) es, eair, emiss, fu

        es = 4.596_dp * Exp(17.27_dp * Te / (237.3_dp + Te))
        eair = 4.596_dp * Exp(17.27_dp * Td / (237.3_dp + Td))
        !gp comment out next 2 lines and replace with new code to use new wind function in module LightAndHeat
        !fUw = 19 + 0.95 * Uw(i) ** 2
        !evap = fUw * (es - eair)
        fu=fUw(lightheat%fUwMethod, Uw, Ta, Te, Ast, eair, es)
        evap = fu * (es - eair)

        !gp 24-Jun-09
        !gp 16-Jul-08
        !emiss=emissivity(lightheat%longatMethod, eair, Ta, cc)
        !emiss=emissivity(lightheat%longatMethod, eair, Ta, cc, lightheat%kbrut)
        emiss=emissivity(lightheat%longatMethod, eair, Ta, cc, lightheat%kbrut, lightheat%KCL1)

        If (lightheat%longatMethod == 'Koberg') Then
            longat = sigma * (Ta + 273.15_dp) ** 4 * emiss * (1 - RL)
        Else

            !gp 24-Jun-09
            !longat = sigma * (Ta + 273.15_dp) ** 4 * emiss * (1 + 0.17_dp * cc ** 2) * (1 - RL)
            longat = sigma * (Ta + 273.15_dp) ** 4 * emiss * (1 + lightheat%KCL2 * cc ** 2) * (1 - RL)

        End If
        !gp comment out next 2 lines and replace with new code to use new wind function in module LightAndHeat
        !conv = Bowen * fUw * (Te(i, 1) - Ta(i))
        conv = Bowen * fu * (Te - Ta)
        back = eps * sigma * (Te + 273.15_dp) ** 4
        !write(*,*) te, es, eair, fu

        !gp 16-Jul=08 adjust longwave for sky opening
        longat = longat * SKOP
        back = back * SKOP

    END SUBROUTINE heatBalance


    !gp 24-Jun-09
    !gp 16-Jul-08
    !PURE FUNCTION Emissivity(longatMethod, eair, Ta, cc)
    !PURE FUNCTION Emissivity(longatMethod, eair, Ta, cc, kbrut)
    PURE FUNCTION Emissivity(longatMethod, eair, Ta, cc, kbrut, KCL1)

        CHARACTER(LEN=20), INTENT(IN) ::longatMethod

        !gp 24-Jun-09
        !gp 16-Jul-08
        !REAL(DP), INTENT(IN) :: eair, Ta, cc
        !REAL(DP), INTENT(IN) :: eair, Ta, cc, kbrut
        REAL(DP), INTENT(IN) :: eair, Ta, cc, kbrut, KCL1

        REAL(DP) Emissivity

        Select Case (longatMethod)

            !gp 16-Jul-08
            !Case ("Brunt")
            !	emissivity = (Acoeff + 0.031_dp * SQRT(eair))
            !Case ("Koberg")
            !	IF (Aa(Ta, (1 - 0.65_dp * cc ** 2.0)) <= 0.735) THEN
            !		emissivity = (Aa(Ta, (1 - 0.65_dp * cc ** 2)) + 0.0263_dp * SQRT(1.333224_dp * eair))
            !	Else
            !		emissivity = (0.735_dp + 0.0263_dp * SQRT(1.333224_dp * eair))
            !	End If
            !Case ("Brutsaert")
            !	emissivity = 1.24_dp * (1.333224_dp * eair / (Ta + 273.15_dp)) ** (1.0_dp / 7.0_dp)   !air vapor pressure is converted to millibars by the factor 1.333224
            !!gp 08-Nov-04 added new options
            !Case ("Satterlund")       !'Satterlund 1979
            !	emissivity = 1.08_dp * (1 - Exp(-(1.333224_dp * eair) ** ((Ta + 273.15_dp) / 2016_dp)))
            !Case ("Idso-Jackson")     !'Idso and Jackson 1969
            !	emissivity = 1_dp - 0.261_dp * Exp(-0.000777_dp * Ta ** 2)
            !Case ("Swinbank")         !'Swinbank 1963
            !	emissivity = 0.0000092_dp * (Ta + 273.15) ** 2		!gp 08-Nov-04 end new block of options
          Case ("Brunt")
            emissivity = (Acoeff + 0.031_dp * SQRT(eair))
          Case ("Koberg")

            !gp 24-Jun-09
            !IF (Aa(Ta, (1_dp - 0.65_dp * cc ** 2.0_dp)) <= 0.735_dp) THEN
            !	emissivity = (Aa(Ta, (1_dp - 0.65_dp * cc ** 2_dp)) + 0.0263_dp * SQRT(1.333224_dp * eair))
            IF (Aa(Ta, (1_dp - KCL1 * cc ** 2.0_dp)) <= 0.735_dp) THEN
                emissivity = (Aa(Ta, (1_dp - KCL1 * cc ** 2_dp)) + 0.0263_dp * SQRT(1.333224_dp * eair))

            Else
                emissivity = (0.735_dp + 0.0263_dp * SQRT(1.333224_dp * eair))
            End If
          Case ("Brutsaert")
            emissivity = kbrut * (1.333224_dp * eair / (Ta + 273.15_dp)) ** (1.0_dp / 7.0_dp)   !air vapor pressure is converted to millibars by the factor 1.333224
          Case ("Satterlund")       !'Satterlund 1979
            emissivity = 1.08_dp * (1_dp - Exp(-(1.333224_dp * eair) ** ((Ta + 273.15_dp) / 2016_dp)))
          Case ("Idso-Jackson")     !'Idso and Jackson 1969
            emissivity = 1_dp - 0.261_dp * Exp(-0.000777_dp * Ta ** 2_dp)
          Case ("Swinbank")         !'Swinbank 1963
            emissivity = 0.0000092_dp * (Ta + 273.15_dp) ** 2_dp
          Case ("Idso 1")
            emissivity = 0.179_dp * ((1.333224_dp * eair) ** (1.0_dp / 7.0_dp)) * Exp(350_dp / (Ta + 273.15_dp)) !'air vapor pressure is converted to millibars by the factor 1.333224
          Case ("Idso 2")
            emissivity = 0.7_dp + 0.0000595_dp * (1.333224_dp * eair) * Exp(1500_dp / (Ta + 273.15_dp)) !'air vapor pressure is converted to millibars by the factor 1.333224

        End Select
    END FUNCTION Emissivity

    ! Aa empirical coefficient (0.5-0.7)
    PURE Function Aa(tempC, clearness)
        !Koberg!s Figure 34 to estimate Brunt!s c coefficient for atmospheric longwave IR
        !inputs: tempC = air temperature deg C
        !        clearness = ratio of estimated measured to clear sky solar radiation
        !output: Koberg!s Brunt!s c coefficient (Aa in Q2K)
        REAL(DP) Aa
        REAL(DP), INTENT(IN) :: tempC, clearness
        REAL(DP) a, b, c

        a = -0.00076437_dp * (clearness ** 3) + 0.00121134_dp * (clearness ** 2) - &
            0.00073087_dp * clearness + 0.0001106_dp
        b = 0.12796842_dp * (clearness ** 3) - 0.2204455_dp * (clearness ** 2) + &
            0.13397992_dp * clearness - 0.02586655_dp
        c = -3.25272249_dp * (clearness ** 3) + 5.65909609_dp * (clearness ** 2) - &
            3.43402413_dp * clearness + 1.43052757_dp
        Aa = a * tempC * tempC + b * tempC + c

    End Function


END MODULE Class_LightHeat
