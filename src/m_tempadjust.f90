module m_tempadjust
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_rates, only: rates_t
    use m_oxygen, only: oxygen_saturation
    use class_hydraulics, only: RiverHydraulics_type
    use m_meteorology, only: t_meteorology
    implicit none
    private
    public :: temp_adjust

contains
    !gp 03-apr-08
    !subroutine tempadjust(rates, hydrau, sitemeteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
    ! kreaft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
    ! kreat, vdit, kact, &
    ! kgent) !gp 30-nov-04 add kgent for generic constituent
    subroutine temp_adjust(rates, hydrau, sitemeteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
        krea1ft, krea2ft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
        kreat, vdit, kact, kgent)

        integer(i32), intent(in) :: nr
        type(rates_t) rates
        type(riverhydraulics_type) hydrau
        type(t_meteorology) sitemeteo

        real(r64), dimension(0:,:), intent(in) :: te

        !gp 03-apr-08
        !real(r64), dimension(:), intent(out) :: khct, kdcst, kdct, kgaft, kdeaft, kreaft, kexaft
        real(r64), dimension(:), intent(out) :: khct, kdcst, kdct, kgaft, kdeaft, krea1ft, krea2ft, kexaft

        real(r64), dimension(:), intent(out) :: khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat
        real(r64), dimension(:), intent(out) :: kdeat, kreat, vdit, kact

        !gp 30-nov-04
        real(r64), dimension(:), intent(out) :: kgent !gp 30-nov-04

        integer(i32) i
        real(r64) kawind, uw10
! real(r64), dimension(nr) :: ka

        do i = 1, nr

            !gp 08-feb-06
            !khct(i) = rates%khc * rates%tkhc ** (te(i, 1) - 20.0_r64)
            !kdcst(i) = rates%kdcs * rates%tkdcs **(te(i,1) - 20.0_r64)
            !kdct(i) = rates%kdc * rates%tkdc ** (te(i, 1) - 20.0_r64)
            !khnt(i) = rates%khn * rates%tkhn ** (te(i, 1) - 20.0_r64)
            !khpt(i) = rates%khp * rates%tkhp ** (te(i, 1) - 20.0_r64)
            !knt(i) = rates%kn * rates%tkn ** (te(i, 1) - 20.0_r64)
            !kit(i) = rates%ki * rates%tki ** (te(i, 1) - 20.0_r64)
            !vdit(i) = rates%vdi * rates%tvdi ** (te(i, 1) - 20.0_r64)
            !
            !!use interpolated wind speed to adjust reaeration for selected methods
            !if (hydrau%reach(i)%kaf == "specified") then
            ! hydrau%reach(i)%ka = hydrau%reach(i)%kau
            !else
            ! uw10 = (10.0_r64 / 7.0_r64) ** 0.15_r64 * sitemeteo%uw(i)
            ! select case (rates%kawindmethod)
            ! case ("banks-herrera") !chapra (1997) eqn 20.46
            ! kawind = 0.728_r64 * uw10 ** 0.5_r64 - 0.317_r64 * uw10 + 0.0372_r64 * uw10 ** 2
            ! case ("wanninkhof") !chapra (1997) eqn 20.47
            ! kawind = 0.0986_r64 * uw10 ** 1.64_r64
            ! case default
            ! kawind = 0
            ! end select
            ! hydrau%reach(i)%ka = hydrau%reach(i)%kau + kawind / hydrau%reach(i)%depth
            !end if
            !kat(i) = hydrau%reach(i)%ka * rates%tka ** (te(i, 1) - 20.0_r64)
            !kact(i) = (32.0_r64 / 44.0_r64) ** 0.25_r64 * kat(i)
! !channel(i)%elev = (channel(i)%elev1 + channel(i)%elev2) / 2
            !hydrau%reach(i)%os = oxsat(te(i, 1), hydrau%reach(i)%elev)
            !kgat(i) = rates%kga * rates%tkga ** (te(i, 1) - 20.0_r64)
            !kdeat(i) = rates%kdea * rates%tkdea ** (te(i, 1) - 20.0_r64)
            !kreat(i) = rates%krea * rates%tkrea ** (te(i, 1) - 20.0_r64)
            !kgaft(i) = rates%kgaf * rates%tkgaf ** (te(i, 1) - 20.0_r64)
            !kdeaft(i) = rates%kdeaf * rates%tkdeaf ** (te(i, 1) - 20.0_r64)
            !kreaft(i) = rates%kreaf * rates%tkreaf ** (te(i, 1) - 20.0_r64)
            !kexaft(i) = rates%kexaf * rates%tkexaf ** (te(i, 1) - 20.0_r64)
            !kdtt(i) = rates%kdt * rates%tkdt ** (te(i, 1) - 20.0_r64)
            !kpatht(i) = rates%kpath * rates%tkpath ** (te(i, 1) - 20.0_r64)
            !
            !!gp 30-nov-04
            !kgent(i) = rates%kgen * rates%tkgen ** (te(i, 1) - 20.0_r64)
            khct(i) = hydrau%reach(i)%khc * rates%tkhc ** (te(i, 1) - 20.0_r64)
            kdcst(i) = hydrau%reach(i)%kdcs * rates%tkdcs **(te(i,1) - 20.0_r64)
            kdct(i) = hydrau%reach(i)%kdc * rates%tkdc ** (te(i, 1) - 20.0_r64)
            khnt(i) = hydrau%reach(i)%khn * rates%tkhn ** (te(i, 1) - 20.0_r64)
            khpt(i) = hydrau%reach(i)%khp * rates%tkhp ** (te(i, 1) - 20.0_r64)
            knt(i) = hydrau%reach(i)%kn * rates%tkn ** (te(i, 1) - 20.0_r64)
            kit(i) = hydrau%reach(i)%ki * rates%tki ** (te(i, 1) - 20.0_r64)
            vdit(i) = hydrau%reach(i)%vdi * rates%tvdi ** (te(i, 1) - 20.0_r64)
            !use interpolated wind speed to adjust reaeration for selected methods
            if (hydrau%reach(i)%kaf == "Specified") then
                hydrau%reach(i)%ka = hydrau%reach(i)%kau
            else
                uw10 = (10.0_r64 / 7.0_r64) ** 0.15_r64 * sitemeteo%uw(i)
                select case (rates%kawindmethod)
                  case ("Banks-Herrera") !chapra (1997) eqn 20.46
                    kawind = 0.728_r64 * uw10 ** 0.5_r64 - 0.317_r64 * uw10 + 0.0372_r64 * uw10 ** 2
                  case ("Wanninkhof") !chapra (1997) eqn 20.47
                    kawind = 0.0986_r64 * uw10 ** 1.64_r64
                  case default
                    kawind = 0
                end select
                hydrau%reach(i)%ka = hydrau%reach(i)%kau + kawind / hydrau%reach(i)%depth
            end if
            kat(i) = hydrau%reach(i)%ka * rates%tka ** (te(i, 1) - 20.0_r64)
            kact(i) = (32.0_r64 / 44.0_r64) ** 0.25_r64 * kat(i)
            hydrau%reach(i)%os = oxygen_saturation(te(i, 1), hydrau%reach(i)%elev)
            kgat(i) = hydrau%reach(i)%kga * rates%tkga ** (te(i, 1) - 20.0_r64)
            kdeat(i) = hydrau%reach(i)%kdea * rates%tkdea ** (te(i, 1) - 20.0_r64)
            kreat(i) = hydrau%reach(i)%krea * rates%tkrea ** (te(i, 1) - 20.0_r64)
            kgaft(i) = hydrau%reach(i)%kgaf * rates%tkgaf ** (te(i, 1) - 20.0_r64)
            kdeaft(i) = hydrau%reach(i)%kdeaf * rates%tkdeaf ** (te(i, 1) - 20.0_r64)

            !gp 03-apr-08
            !kreaft(i) = hydrau%reach(i)%kreaf * rates%tkreaf ** (te(i, 1) - 20.0_r64)
            krea1ft(i) = hydrau%reach(i)%krea1f * rates%tkreaf ** (te(i, 1) - 20.0_r64)
            krea2ft(i) = hydrau%reach(i)%krea2f !no temp adj of photo resp because growth rate is temp adjusted

            kexaft(i) = hydrau%reach(i)%kexaf * rates%tkexaf ** (te(i, 1) - 20.0_r64)
            kdtt(i) = hydrau%reach(i)%kdt * rates%tkdt ** (te(i, 1) - 20.0_r64)
            kpatht(i) = hydrau%reach(i)%kpath * rates%tkpath ** (te(i, 1) - 20.0_r64)
            kgent(i) = hydrau%reach(i)%kgen * rates%tkgen ** (te(i, 1) - 20.0_r64)

        end do
    end subroutine temp_adjust

end module m_tempadjust
