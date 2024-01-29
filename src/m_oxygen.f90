module m_oxygen
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_rates, only: rates_t
    implicit none
    private
    public oxygen_inhibition_and_enhancement, oxygen_saturation

contains
    pure function oxygen_saturation(temp, elev)

        real(r64) oxygen_saturation
        real(r64), intent(in) :: temp, elev
        real(r64) taa, lnosf

        taa = temp + 273.15_r64
        lnosf = -139.34411_r64 &
            + 157570.1_r64 / taa &
            - 66423080.0_r64 / taa ** 2 &
            + 12438000000.0_r64 / taa ** 3 &
            - 862194900000.0_r64 / taa ** 4
        oxygen_saturation = exp(lnosf) * (1 - 0.0001148_r64 * elev)

    end function oxygen_saturation

    subroutine oxygen_inhibition_and_enhancement(rates, o2, fcarb, fnitr, fdenitr, frespp, frespb)

        real(r64), intent(out) :: fcarb, fnitr, fdenitr, frespp, frespb
        real(r64), intent(in) :: o2
        type(rates_t), intent(in) :: rates

        fcarb = oxygen_inhibition(rates%ikoxc, rates%ksocf, o2) !oxygen inhibition of carbon oxidation
        fnitr = oxygen_inhibition(rates%ikoxn, rates%ksona, o2) !oxygen inhibition of nitrification
        fdenitr = oxygen_enhancement(rates%ikoxdn, rates%ksodn, o2) !oxygen enhancement of denitrification
        frespp = oxygen_inhibition(rates%ikoxp, rates%ksop, o2) !oxygen inhibition of phytoplankton respiration
        frespb = oxygen_inhibition(rates%ikoxb, rates%ksob, o2) !oxygen inhibition of bottom plant respiration

    end subroutine oxygen_inhibition_and_enhancement

    pure function oxygen_inhibition(model, rate, oxygen)
        integer(i32), intent(in) ::  model
        real(r64), intent(in) :: rate, oxygen
        real(r64) :: oxygen_inhibition
        select case (model)
          case (1)
            oxygen_inhibition = half_saturation_attenuation(rate, oxygen)
          case (2)
            oxygen_inhibition = exponential_attenuation(rate, oxygen)
          case (3)
            oxygen_inhibition = second_order_attenuation(rate, oxygen)
        end select
    end function oxygen_inhibition

    pure function oxygen_enhancement(model, rate, oxygen)
        ! oxygen_enhancement = 1 - oxygen_attenuation
        integer(i32), intent(in) ::  model
        real(r64), intent(in) :: rate, oxygen
        real(r64) :: oxygen_enhancement
        oxygen_enhancement = 1 - oxygen_inhibition(model, rate, oxygen)
    end function oxygen_enhancement

    pure function half_saturation_attenuation(rate, oxygen)
        ! oxygen attenuation: half-saturation (EQ122)
        real(r64), intent(in) :: rate, oxygen
        real(r64) :: half_saturation_attenuation
        half_saturation_attenuation = oxygen / (rate + oxygen)
    end function half_saturation_attenuation

    pure function exponential_attenuation(rate, oxygen)
        ! oxygen attenuation: exponential (EQ123)
        real(r64), intent(in) :: rate, oxygen
        real(r64) :: exponential_attenuation
        exponential_attenuation = 1.0_r64 - exp(-rate * oxygen)
    end function exponential_attenuation

    pure function second_order_attenuation(rate, oxygen)
        ! oxygen attenuation: second-order half saturation (EQ124)
        real(r64), intent(in) :: rate, oxygen
        real(r64) :: second_order_attenuation
        second_order_attenuation = oxygen ** 2 / (rate + oxygen ** 2)
    end function second_order_attenuation

end module m_oxygen
