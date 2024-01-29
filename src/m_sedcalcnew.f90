module m_sedcalcnew
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_hydraulics, only: riverhydraulics_type
    implicit none
    private
    public :: sedcalcnumnew

contains

    !gp 25-aug-08
    !subroutine sedcalcnumnew(jcin, jnin, jpin, o20, depth, tw, sod, jnh4, jno3, jch4, &
    ! jch4g, jpo4, nh30, no30, po40, ch40, csod)
    subroutine sedcalcnumnew(jcin, jnin, jpin, o20, depth, tw, sod, jnh4, jno3, jch4, &
        jch4g, jpo4, nh30, no30, po40, ch40, csod, calcsedflux)

        real(r64), parameter :: kappanh3 = 0.131_r64, thtanh3 = 1.123_r64, km_nh3 = 0.728_r64, km_o2_nh3 = 0.37_r64
        real(r64), parameter :: kappano3_1 = 0.1_r64, kappano3_2 = 0.25_r64, thtano3 = 1.08_r64
        real(r64), parameter :: kappach4 = 0.7_r64, thtach4 = 1.079_r64
        real(r64), parameter :: dd = 0.001_r64, thtadd = 1.08_r64, dp0 = 0.00012_r64, thtadp = 1.117_r64
        real(r64), parameter :: kdnh3 = 1.0_r64
        real(r64), parameter :: m1 = 0.5_r64, m2 = 0.5_r64
        real(r64), parameter :: w2 = 0.000005_r64, h2 = 0.1_r64
        real(r64), parameter :: dkdpo41 = 20.0_r64
        real(r64), parameter :: kdpo42 = 20.0_r64
        real(r64), parameter :: o2critpo4 = 2.0_r64

        real(r64), intent(in) :: jcin, jnin, jpin, o20, depth, tw, ch40, po40, nh30, no30
        character(len=30), intent(in) :: calcsedflux
        real(r64), intent(out) :: sod, jnh4, jno3, jch4, jch4g, jpo4, csod

        !gp 25-aug-08

        integer(i32) i, maxit, iter
        real(r64) es, ea
        real(r64) nsod, sodold

        real(r64) w12, kl12
        real(r64) s
        real(r64) nh3(2), nh3t(2)
        real(r64) ch4(2)
        real(r64) nh3tono3, ch4toco2
        real(r64) jnh3tono3
        real(r64) denit(2)
        real(r64) no3(2)
        real(r64) jdenit(2), jdenitt
        real(r64) jo2no3(2), jo2no3t, jc_o2equiv
        real(r64) ch4sat, csodmax
        real(r64) sech_arg
        real(r64) fd1, fp1, fd2, fp2
        real(r64) a11, a12, a21, a22, b1, b2

        real(r64) fpon(3), fpoc(3), fpop(3)
        real(r64) kdiapoc(3), thtapoc(3)
        real(r64) kdiapon(3), thtapon(3)
        real(r64) kdiapop(3), thtapop(3)
        real(r64) jpoc(3), jpon(3), jpop(3)
        real(r64) poc2(3), pon2(3), pop2(3)
        real(r64) poct2, pont2, popt2
        real(r64) jc, jn, jp

        !phosphorus
        real(r64) kdpo41, po4t(2), po4(2)

        nh3=0; nh3t=0; ch4=0
        denit=0; no3=0; jdenit=0; jo2no3=0
        kdpo41=0; po4t=0; po4=0
        !compute influxes corrected for settling and refractory
        fpon(1) = 0.65_r64; fpon(2) = 0.25_r64; fpon(3) = 1.0_r64 - fpon(1) - fpon(2)
        fpoc(1) = 0.65_r64; fpoc(2) = 0.2_r64; fpoc(3) = 1.0_r64 - fpoc(1) - fpoc(2)
        fpop(1) = 0.65_r64; fpop(2) = 0.2_r64; fpop(3) = 1.0_r64 - fpop(1) - fpop(2)

        kdiapon(1) = 0.035_r64; thtapon(1) = 1.1_r64
        kdiapon(2) = 0.0018_r64; thtapon(2) = 1.15_r64
        kdiapon(3) = 0; thtapon(3) = 1.17_r64

        kdiapoc(1) = 0.035_r64; thtapoc(1) = 1.1_r64
        kdiapoc(2) = 0.0018_r64; thtapoc(2) = 1.15_r64
        kdiapoc(3) = 0; thtapoc(3) = 1.17_r64

        kdiapop(1) = 0.035_r64; thtapop(1) = 1.1_r64
        kdiapop(2) = 0.0018_r64; thtapop(2) = 1.15_r64
        kdiapop(3) = 0; thtapop(3) = 1.17_r64

        !compute input fluxes
        do i = 1, 3
            jpoc(i) = jcin * fpoc(i)
            jpon(i) = jnin * fpon(i)
            jpop(i) = jpin * fpop(i)
        end do

        !compute particulate organic forms
        do i = 1, 3
            poc2(i) = jpoc(i) / h2 / (kdiapoc(i) * thtapoc(i) ** (tw - 20.0_r64) + w2 / h2)
            pon2(i) = jpon(i) / h2 / (kdiapon(i) * thtapon(i) ** (tw - 20.0_r64) + w2 / h2)
            pop2(i) = jpop(i) / h2 / (kdiapop(i) * thtapop(i) ** (tw - 20.0_r64) + w2 / h2)
            poct2 = poct2 + poc2(i)
            pont2 = pont2 + pon2(i)
            popt2 = popt2 + pop2(i)
        end do

        !compute diagenesis fluxes
        jc = 0; jn = 0; jp = 0
        do i = 1, 3
            jc = jc + kdiapoc(i) * thtapoc(i) ** (tw - 20) * poc2(i) * h2
            jn = jn + kdiapon(i) * thtapon(i) ** (tw - 20) * pon2(i) * h2
            jp = jp + kdiapop(i) * thtapop(i) ** (tw - 20) * pop2(i) * h2
        end do

        maxit = 500
        es = 0.1_r64
        w12 = dp0 * thtadp ** (tw - 20.0_r64) / (h2 / 2.0_r64)
        kl12 = dd * thtadd ** (tw - 20.0_r64) / (h2 / 2.0_r64)

        sodold = jc + 1.714_r64 * jn
        iter = 0
        ! saturation conc. of methane in oxygen equivalent units (equation 10.51)
        ch4sat = 100.0_r64 * (1.0_r64 + depth / 10.0_r64) * (1.024_r64 ** (20.0_r64 - tw)) ![gmo*/m3]
        if ((o20 < 0.001_r64) .or. (jc + jn <= 0)) then
            sod = 0
            jnh4 = jn
            jch4 = min(sqrt(2.0_r64 * kl12 * ch4sat * jc), jc)
            jch4g = jc - jch4
            jpo4 = jp
        else
            do
                s = sodold / o20
                nh3tono3 = kappanh3 ** 2 * thtanh3 ** (tw - 20.0_r64) / s * km_nh3 / (km_nh3 + nh3(1)) * &
                    o20 / (2.0_r64 * km_o2_nh3 + o20)
                ! [m/d] = [m2/d2] / [m/d] * [-] * [-]

                ! calculate dissolved and particulate (sorbed) fractions
                fd1 = (1.0_r64 / (1.0_r64 + m1 * kdnh3))
                fp1 = 1.0_r64 - fd1 != ((m1*kdnh3)/(1 + m1*kdnh3))
                fd2 = (1.0_r64 / (1.0_r64 + m2 * kdnh3))
                fp2 = 1.0_r64 - fd2 != ((m2*kdnh3)/(1 + m2*kdnh3))

                ! write linear sys of equations around nh3t
                a11 = -fd1 * kl12 - fp1 * w12 - fd1 * nh3tono3 - fd1 * s - w2
                a12 = fd2 * kl12 + fp2 * w12
                b1 = -s * nh30 ![m/d]*[mg/m3]
                a21 = fd1 * kl12 + fp1 * w12 + w2
                a22 = -fd2 * kl12 - fp2 * w12 - w2
                b2 = -jn

                call lin_sys(a11, a12, a21, a22, b1, b2, nh3t(1), nh3t(2))

                ! dissolved concentrations
                nh3(1) = fd1 * nh3t(1)
                nh3(2) = fd2 * nh3t(2)

                ! oxygen flux due to nh3->no2, see equation 23.3 in chapra (1997)
                jnh3tono3 = nh3tono3 * nh3(1)
                nsod = 2.0_r64 * (32.0_r64 / 14.0_r64) * nh3tono3 * nh3(1)
                ! [gmo/m2-d] = [mol o2/mol n]*[gm o2/mol o2]/[gm n/mol n]*
                ! [gm/1000mg] * [m/day] * [mgn/m3]

                !:::::::::::::::::::::::::::: begin nitrate::::::::::::::::::::::::::::

                ! denitrification in layers 1 and 2 (equation 4.55)
                denit(1) = (kappano3_1 ** 2 * thtano3 ** (tw - 20.0_r64) / s)
                denit(2) = kappano3_2 * thtano3 ** (tw - 20.0_r64)

                ! layer 1
                a11 = -kl12 - denit(1) - s - w2
                a12 = kl12
                b1 = -s * no30 - nh3tono3 * nh3(1)
                ! layer 2
                a21 = kl12 + w2
                a22 = -kl12 - denit(2) - w2
                b2 = 0.0

                call lin_sys(a11, a12, a21, a22, b1, b2, no3(1), no3(2))

                ! nitrate flux to water column
                jno3 = s * (no3(1) - no30)

                ! denitrification flux [mgn/m2d]
                jdenit(1) = denit(1) * no3(1)
                jdenit(2) = denit(2) * no3(2)
                jdenitt = jdenit(1) + jdenit(2)

                ! methane consumption due to denitrification (equation 9.16)
                ! layer 1
                jo2no3(1) = 2.0_r64 * (16.0_r64 / 12.0_r64) * (10.0_r64 / 8.0_r64) * (12.0_r64 / 14.0_r64) * jdenit(1)
                ! [gmo*/m2-d] = [molo/molc]*[(gmo/molo)/(gmc/molc)]*[molc/moln]*
                ! [(gmc/molc)/(gmn/moln)] * [mg/m2d] * [gm/1000mg]
                ! where 2*16/12 is the ubiquitous 32/12 (= 2.67) for oxidation of carbon
                ! layer 2
                jo2no3(2) = (32.0_r64 / 12.0_r64) * (10.0_r64 / 8.0_r64) * (12.0_r64 / 14.0_r64) * jdenit(2)
                ! sum
                jo2no3t = jo2no3(1) + jo2no3(2)

                ! calculate methane flux in oxygen equivalent units, adjusted for
                ! the methane consumed in denitrification
                jc_o2equiv = jc - jo2no3t
                if (jc_o2equiv < 0) jc_o2equiv = 0.0

                !***** methane in o2 equivalents
                !:::::::::::::::::::::::::::: begin methane::::::::::::::::::::::::::::


                ! csodmax equations 10.28 and 10.30
                csodmax = min(sqrt(2.0_r64 * kl12 * ch4sat * jc_o2equiv), jc_o2equiv) ! [gmo*/m2-d] = sqr([m/d] * [gmo*/m3] * [gmo*/m2-d])
                !msgbox csodmax
                !sech_arg = (kappach4 * thtach4 ** ((tw - 20) / 2.0)) / s
                ! csod equation 10.35
                ! the hyperbolic secant is defined as hsec(x) = 2 / (exp(x) + exp(-x))
                !if (sech_arg < 400.0) then !this is the usual case
                ! csod = csodmax * (1.0 - (2.0 / (exp(sech_arg) + exp(-sech_arg))))
                !else !hsec(sech_arg) < 3.8e-174 ~ 0
                ! csod = csodmax
                !end if
                ! aqueous methane flux to water column
                !jch4 = csodmax - csod
                ! gaseous methane flux to water column
                !jch4g = jc_o2equiv - jch4 - csod

                ch4toco2 = (kappach4 * kappach4 * thtach4 ** ((tw - 20.0_r64) / 2.0_r64)) / s
                ch4(1) = (csodmax + s * ch40) / (ch4toco2 + s)
                csod = ch4toco2 * ch4(1)

                sod = (sodold + csod + nsod) / 2.0_r64
                iter = iter + 1
                ea = abs((sod - sodold) / sod) * 100.0_r64

                !gp 20-may-09
                !if (ea <= es) then
                if (ea <= es .or. (calcsedflux == "Option 2" .and. iter >= maxit)) then

                    exit


                    !gp 20-may-09
                    !!gp 25-aug-08
                    !!elseif (iter >= maxit) then
                    !! write(*,*) "sod iterations exceeded! turn off sediment diagenesis on qual2k sheet."
                    !elseif ((calcsedflux == "yes" .or. calcsedflux == "option 1") .and. iter >= maxit) then
                    ! write(*,*) "sod iterations exceeded! turn off sediment diagenesis on qual2k sheet."
                elseif (iter >= maxit) then
                    write(*,*) "SOD iterations exceeded! Turn off sediment diagenesis on QUAL2K sheet."


                end if
                sodold = sod
            end do
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ! begin orthophosphate
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            ! determine kd in layer 1 (kd in layer 2 is a constant)
            if (o20 > o2critpo4) then
                kdpo41 = kdpo42 * dkdpo41
            else
                kdpo41 = kdpo42 * (dkdpo41 ** (o20 / o2critpo4))
            end if

            ! calculate dissolved and particulate (sorbed) fractions
            fd1 = (1.0_r64 / (1.0_r64 + m1 * kdpo41))
            fp1 = 1.0_r64 - fd1 != ((m1*kdpo41)/(1 + m1*kdpo41))
            fd2 = (1.0_r64 / (1.0_r64 + m2 * kdpo42))
            fp2 = 1.0_r64 - fd2 != ((m2*kdpo42)/(1 + m2*kdpo42))

            ! write linear sys of equations around po4t
            ! layer 1
            a11 = -fd1 * kl12 - fp1 * w12 - fd1 * s - w2
            a12 = fd2 * kl12 + fp2 * w12
            b1 = -s * po40
            ! layer 2
            a21 = fd1 * kl12 + fp1 * w12 + w2
            a22 = -fd2 * kl12 - fp2 * w12 - w2
            b2 = -jp

            call lin_sys(a11, a12, a21, a22, b1, b2, po4t(1), po4t(2))

            ! dissolved po4 concentrations
            po4(1) = fd1 * po4t(1)
            po4(2) = fd2 * po4t(2)

            ! po4 flux to water column
            jpo4 = s * (po4(1) - po40)
            jnh4 = s * (nh3(1) - nh30)
            jno3 = s * (no3(1) - no30)
            jch4 = s * (ch4(1) - ch40)
            jch4g = jc_o2equiv - jch4 - csod
        end if

    end subroutine


    subroutine lin_sys(a11, a12, a21, a22, b1, b2, x1, x2)

        real(r64), intent(in) :: a11, a12, a21, a22, b1, b2
        real(r64), intent(out) :: x1, x2
        !this subroutine solves a linear sys of 2 equations and 2 unknowns

        x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21)
        x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21)

    end subroutine


end module m_sedcalcnew
