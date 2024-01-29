! classlightheat.f90
! general variables and functions for calculating light and heat
module class_lightheat
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_constants, only: acoeff, bowen, eps, rl, sigma
    implicit none
    private
    public :: lightheat, light_, heatbalance, lightextinction

    type lightheat_type
        real(r64) :: par = 0.47 !phtosynthetically available radiation, default 0.47
        real(r64) :: kep = 0.01 !backgroud light extinction
        real(r64) :: kela=0 !linear chlorophyll light extinction
        real(r64) :: kenla=0 !nonlinear chlorophyll light extinction
        real(r64) :: kess=0 !inorganic ss light extinction
        real(r64) :: kepom=0 !detritus light extinction

        !gp 13-feb-06
        real(r64) :: kemac=0 !macrophyte light extinction

        character(len=30) ::longatmethod = "Brunt"
        !atmospheric longwave emissivity model
        !gp 16-jul-08
        real(r64) :: kbrut=1.24

        !gp 24-jul-09
        real(r64) :: kcl1=0.65
        real(r64) :: kcl2=0.17

        character(len=30) ::fuwmethod = "Brady-Graves-Geyer"
        !evaporation and air convection/conduction
    end type

    type(lightheat_type) lightheat

contains

    !gp 13-feb-06
    !pure function lightextinction(lightheat, cpom, ss, calgae) result(ke)
    !
    ! type(lightheat_type), intent(in) :: lightheat
    ! real(r64), intent(in) :: cpom, ss, calgae
    ! real(r64) ke
    !
    ! ke = lightheat%kep + lightheat%kepom * cpom + lightheat%kess * ss &
    ! + lightheat%kela * calgae + lightheat%kenla * calgae ** (2.0 / 3.0)
    !
    !end function lightextinction
    pure function lightextinction(lightheat, cpom, ss, calgae, cmacrophyte) result(ke)
        type(lightheat_type), intent(in) :: lightheat
        real(r64), intent(in) :: cpom, ss, calgae, cmacrophyte !note: cmacrophyte is gd/m^3
        real(r64) ke
        ke = lightheat%kep + lightheat%kepom * cpom + lightheat%kess * ss + lightheat%kemac * cmacrophyte &
            + lightheat%kela * calgae + lightheat%kenla * calgae ** (2.0 / 3.0)
    end function lightextinction

    !/* public functions */

    !gp 13-feb-06
    !subroutine light_(par, kep, kela, kenla, kess, kepom, longatmethod, fuwmethod)

    !gp 24-jun-09
    !gp 16-jul-08
    !subroutine light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, fuwmethod)
    !subroutine light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, kbrut, fuwmethod)
    subroutine light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, kbrut, fuwmethod, kcl1, kcl2)

        !gp 13-feb-06
        !real(r64) par, kep, kela, kenla, kess, kepom

        !gp 24-jun-09
        !gp 16-jul-08
        !real(r64) par, kep, kela, kenla, kess, kepom, kemac
        !real(r64) par, kep, kela, kenla, kess, kepom, kemac, kbrut
        real(r64) par, kep, kela, kenla, kess, kepom, kemac, kbrut, kcl1, kcl2

        character(len=30) longatmethod, fuwmethod

        lightheat%par = par
        lightheat%kep = kep
        lightheat%kela = kela
        lightheat%kenla=kenla
        lightheat%kess = kess
        lightheat%kepom = kepom

        !gp 13-feb-06
        lightheat%kemac = kemac

        lightheat%longatmethod = longatmethod

        !gp 16-jul-08
        lightheat%kbrut = kbrut

        lightheat%fuwmethod= fuwmethod

        !gp 24-jun-09
        lightheat%kcl1= kcl1
        lightheat%kcl2= kcl2

    end subroutine light_

    function fuw(fuwmethod, uw, ta, te, ast, eair, es)

        real(r64) fuw
        character(len=20), intent(in) ::fuwmethod
        real(r64),intent(in) :: uw, ta, te, ast, eair, es
        real(r64) wmph, areaa
        real(r64) ta2a, tva, tsa, tvs, dtv

        select case (fuwmethod)
          case ("Brady-Graves-Geyer")
            !"brady, graves, and geyer"
            !for this formula the wind speed height is 7 m (see edinger, et al., 1974)
            fuw = 19.0_r64 + 0.95_r64 * uw * uw !chapra eqn 30.22 cal/cm**2/d/mmhg
          case ("adams 1")
            !write(*,*) "adams 1, wrong!"
            !"east mesa"
            !from adams, et al., eq. 4.48, p. 4-26
            !for this formula the wind speed height is 2m
            !first convert windspeed from m s-1 to mph
            wmph = uw * 3600.0_r64 / (0.3048_r64 * 5280.0_r64) !uw(i) is at 7m
            !convert to the formula!s wind speed height using
            !the exponential wind law, paily et al., 1974
            wmph = wmph * (2.0_r64 / 7.0_r64) ** 0.15_r64 !convert wind speed from 7m to 2m height
            !convert surface area in m**2 to area in acres
            areaa = ast / (0.3048_r64 ** 2.0 * 43560.0_r64)
            !next compute virtual temperature difference
            !eagleson, p. s. 1970. dynamic hydrology. mcgraw-hill, inc.,
            !new york, new york. p. 56
            ta2a = 9.0_r64 / 5.0_r64 * ta + 32.0_r64 !deg f air temp at 2m
            tva = ta2a / (1.0_r64 - 0.378_r64 * (eair / 760.0_r64))
            tsa = 9.0_r64 / 5.0_r64 * te + 32.0_r64 !te(i, 1)
            tvs = tsa / (1.0_r64 - 0.378_r64 * (es / 760.0_r64))
            dtv = tvs - tva !original formula in shestt
            if (dtv < 0) dtv = 0
            !next compute fuw in w/m**2/mmhg
            fuw = 0.1313_r64 * ((22.4_r64 * dtv ** (1.0_r64 / 3.0_r64)) ** 2 &
              + (24.2_r64 * areaa ** (-0.05_r64) * wmph) ** 2) ** 0.5_r64
            !next convert fuw to cal/cm**2/d/mmhg
            fuw = fuw / (4.183076_r64 * 100.0_r64 * 100.0_r64 / 86400.0_r64)
          case ("Adams 2")
            !write(*,*) "adams2, wrong!"
            !!"east mesa modified to include marciano and harbeck"
            !!from adams, et al., eq. 4.48, p. 4-26
            !for this formula the wind speed height is 2m
            !first convert windspeed from m s-1 to mph
            wmph = uw * 3600.0_r64 / (0.3048_r64 * 5280.0_r64) !uw(i) is at 7m
            !convert to the formula!s wind speed height using
            !the exponential wind law, paily et al., 1974
            wmph = wmph * (2.0_r64 / 7.0_r64) ** 0.15_r64 !convert wind speed from 7m to 2m height
            !next compute virtual temperature difference
            !eagleson, p. s. 1970. dynamic hydrology. mcgraw-hill, inc.,
            !new york, new york. p. 56
            ta2a = 9.0_r64 / 5.0_r64 * ta + 32.0_r64 !deg f air temp at 2m
            tva = ta2a / (1.0_r64 - 0.378_r64 * (eair / 760.0_r64))
            tsa = 9.0_r64 / 5.0_r64 * te + 32.0_r64 !te(i, 1)
            tvs = tsa / (1.0_r64 - 0.378_r64 * (es / 760.0_r64))
            dtv = tvs - tva !original formula in shestt
            if (dtv < 0) dtv = 0
            !next compute fuw in w/m**2/mmhg
            fuw = 0.1313_r64 * ((22.4_r64 * dtv ** (1.0_r64 / 3.0_r64)) ** 2 + (17.0_r64 * wmph) ** 2) ** 0.5_r64
            !next convert fuw to cal/cm**2/d/mmhg
            fuw = fuw / (4.183076_r64 * 100.0_r64 * 100.0_r64 / 86400.0_r64)
          case default
            !write(*,*) "default, wrong!"
        end select

    end function

!evaperation, convection, back radiation, air longat

    !gp 16-jul-08
    !subroutine heatbalance(evap, conv, back, longat, lightheat, te, ta, td, cc, uw, ast)
    subroutine heatbalance(evap, conv, back, longat, lightheat, te, ta, td, cc, uw, ast, skop)

        type(lightheat_type), intent(in) :: lightheat

        !gp 16-jul-08
        !real(r64), intent(in) :: te, ta, td, cc, uw, ast
        real(r64), intent(in) :: te, ta, td, cc, uw, ast, skop

        real(r64), intent(out) :: evap, conv, back, longat
        real(r64) es, eair, emiss, fu

        es = 4.596_r64 * exp(17.27_r64 * te / (237.3_r64 + te))
        eair = 4.596_r64 * exp(17.27_r64 * td / (237.3_r64 + td))
        !gp comment out next 2 lines and replace with new code to use new wind function in module lightandheat
        !fuw = 19 + 0.95 * uw(i) ** 2
        !evap = fuw * (es - eair)
        fu=fuw(lightheat%fuwmethod, uw, ta, te, ast, eair, es)
        evap = fu * (es - eair)

        !gp 24-jun-09
        !gp 16-jul-08
        !emiss=emissivity(lightheat%longatmethod, eair, ta, cc)
        !emiss=emissivity(lightheat%longatmethod, eair, ta, cc, lightheat%kbrut)
        emiss=emissivity(lightheat%longatmethod, eair, ta, cc, lightheat%kbrut, lightheat%kcl1)

        if (lightheat%longatmethod == 'Koberg') Then
            longat = sigma * (ta + 273.15_r64) ** 4 * emiss * (1 - rl)
        else

            !gp 24-jun-09
            !longat = sigma * (ta + 273.15_r64) ** 4 * emiss * (1 + 0.17_r64 * cc ** 2) * (1 - rl)
            longat = sigma * (ta + 273.15_r64) ** 4 * emiss * (1 + lightheat%kcl2 * cc ** 2) * (1 - rl)

        end if
        !gp comment out next 2 lines and replace with new code to use new wind function in module lightandheat
        !conv = bowen * fuw * (te(i, 1) - ta(i))
        conv = bowen * fu * (te - ta)
        back = eps * sigma * (te + 273.15_r64) ** 4
        !write(*,*) te, es, eair, fu

        !gp 16-jul=08 adjust longwave for sky opening
        longat = longat * skop
        back = back * skop

    end subroutine heatbalance


    !gp 24-jun-09
    !gp 16-jul-08
    !pure function emissivity(longatmethod, eair, ta, cc)
    !pure function emissivity(longatmethod, eair, ta, cc, kbrut)
    pure function emissivity(longatmethod, eair, ta, cc, kbrut, kcl1)

        character(len=20), intent(in) ::longatmethod

        !gp 24-jun-09
        !gp 16-jul-08
        !real(r64), intent(in) :: eair, ta, cc
        !real(r64), intent(in) :: eair, ta, cc, kbrut
        real(r64), intent(in) :: eair, ta, cc, kbrut, kcl1

        real(r64) emissivity

        select case (longatmethod)

            !gp 16-jul-08
            !case ("brunt")
            ! emissivity = (acoeff + 0.031_r64 * sqrt(eair))
            !case ("koberg")
            ! if (aa(ta, (1 - 0.65_r64 * cc ** 2.0)) <= 0.735) then
            ! emissivity = (aa(ta, (1 - 0.65_r64 * cc ** 2)) + 0.0263_r64 * sqrt(1.333224_r64 * eair))
            ! else
            ! emissivity = (0.735_r64 + 0.0263_r64 * sqrt(1.333224_r64 * eair))
            ! end if
            !case ("brutsaert")
            ! emissivity = 1.24_r64 * (1.333224_r64 * eair / (ta + 273.15_r64)) ** (1.0_r64 / 7.0_r64) !air vapor pressure is converted to millibars by the factor 1.333224
            !!gp 08-nov-04 added new options
            !case ("satterlund") !'satterlund 1979
            ! emissivity = 1.08_r64 * (1 - exp(-(1.333224_r64 * eair) ** ((ta + 273.15_r64) / 2016_r64)))
            !case ("idso-jackson") !'idso and jackson 1969
            ! emissivity = 1_r64 - 0.261_r64 * exp(-0.000777_r64 * ta ** 2)
            !case ("swinbank") !'swinbank 1963
            ! emissivity = 0.0000092_r64 * (ta + 273.15) ** 2 !gp 08-nov-04 end new block of options
          case ("Brunt")
            emissivity = (acoeff + 0.031_r64 * sqrt(eair))
          case ("Koberg")

            !gp 24-jun-09
            !if (aa(ta, (1_r64 - 0.65_r64 * cc ** 2.0_r64)) <= 0.735_r64) then
            ! emissivity = (aa(ta, (1_r64 - 0.65_r64 * cc ** 2_r64)) + 0.0263_r64 * sqrt(1.333224_r64 * eair))
            if (aa(ta, (1_r64 - kcl1 * cc ** 2.0_r64)) <= 0.735_r64) then
                emissivity = (aa(ta, (1_r64 - kcl1 * cc ** 2_r64)) + 0.0263_r64 * sqrt(1.333224_r64 * eair))

            else
                emissivity = (0.735_r64 + 0.0263_r64 * sqrt(1.333224_r64 * eair))
            end if
          case ("Brutsaert")
            emissivity = kbrut * (1.333224_r64 * eair / (ta + 273.15_r64)) ** (1.0_r64 / 7.0_r64) !air vapor pressure is converted to millibars by the factor 1.333224
          case ("Satterlund") !'satterlund 1979
            emissivity = 1.08_r64 * (1_r64 - exp(-(1.333224_r64 * eair) ** ((ta + 273.15_r64) / 2016_r64)))
          case ("Idso-Jackson") !'idso and jackson 1969
            emissivity = 1_r64 - 0.261_r64 * exp(-0.000777_r64 * ta ** 2_r64)
          case ("Swinbank") !'swinbank 1963
            emissivity = 0.0000092_r64 * (ta + 273.15_r64) ** 2_r64
          case ("Idso 1")
            emissivity = 0.179_r64 * ((1.333224_r64 * eair) ** (1.0_r64 / 7.0_r64)) * exp(350_r64 / (ta + 273.15_r64)) !'air vapor pressure is converted to millibars by the factor 1.333224
          case ("Idso 2")
            emissivity = 0.7_r64 + 0.0000595_r64 * (1.333224_r64 * eair) * exp(1500_r64 / (ta + 273.15_r64)) !'air vapor pressure is converted to millibars by the factor 1.333224

        end select
    end function emissivity

    ! aa empirical coefficient (0.5-0.7)
    pure function aa(tempc, clearness)
        !koberg!s figure 34 to estimate brunt!s c coefficient for atmospheric longwave ir
        !inputs: tempc = air temperature deg c
        ! clearness = ratio of estimated measured to clear sky solar radiation
        !output: koberg!s brunt!s c coefficient (aa in q2k)
        real(r64) aa
        real(r64), intent(in) :: tempc, clearness
        real(r64) a, b, c

        a = -0.00076437_r64 * (clearness ** 3) + 0.00121134_r64 * (clearness ** 2) - &
            0.00073087_r64 * clearness + 0.0001106_r64
        b = 0.12796842_r64 * (clearness ** 3) - 0.2204455_r64 * (clearness ** 2) + &
            0.13397992_r64 * clearness - 0.02586655_r64
        c = -3.25272249_r64 * (clearness ** 3) + 5.65909609_r64 * (clearness ** 2) - &
            3.43402413_r64 * clearness + 1.43052757_r64
        aa = a * tempc * tempc + b * tempc + c

    end function


end module class_lightheat
