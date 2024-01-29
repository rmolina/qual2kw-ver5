module class_phsolve
    use nrtype
    implicit none
    private
    public :: ct, ph_solver, chemrates, modfp2

contains

    !gp 23-nov-09
    !!#####################################################################
    !!12/15/07
    !subroutine phsolver(imethph, ph, ct, te, alk, cond)
    !character(*), intent(in) :: imethph ! ph solve method
    !real(dp), intent(out) :: ph
    !real(dp), intent(in) :: ct, te, alk, cond
    !
    !if (imethph == "newton-raphson") then
    ! call phsolnewton(ph, ct, te, alk, cond)
    !elseif (imethph == "bisection") then !09/23/07
    ! call phsolbisect(ph, ct, te, alk, cond)
    !else !brent method
    ! call phsolfzerotx(ph, ct, te, alk, cond)
    !end if
    !
    !end subroutine
    !!#####################################################################

    !09/23/07 new ph solver, brent method, used as default method

    !gp 23-nov-09
    !subroutine phsolfzerotx(xr, ct, te, alk, cond)
    subroutine phsolbrent(xr, ct, te, alk, cond)

        !fzerotx textbook version of fzero.
        ! x = fzerotx(f,[a,b]) tries to find a zero of f(x) between a and b.
        ! f(a) and f(b) must have opposite signs. fzerotx returns one
        ! end point of a small subinterval of [a,b] where f changes sign.

        !in&out varaibles
        real(dp), intent(out) :: xr
        real(dp), intent(in) :: ct, te, alk, cond

        !local variables
        real(dp) a, b, c, d, e
        real(dp) fa, fb, fc
        real(dp) p, q, r, s
        real(dp) m, tol
        real(dp) alke
        real(dp), parameter:: eps2 = 2.22044604925031e-16

        ! initialize
        a = 0.0
        b = 14.0

        !gp 03-dec-09
        !alke = alk / 50000.0
        alke = alk / 50043.45_dp

        fa = f(a, ct, te, alke, cond)
        fb = f(b, ct, te, alke, cond)
        if (fa * fb >=0) then !if (sgn(fa) == sgn(fb)) then
            !msgbox "function must change sign on the interval"
            !end

            !gp 23-nov-09
            !write(log_file,*) "error: bad ph guesses. try a smaller calculation step!"
            !call abort("bad ph guesses. try a smaller calculation step!")
            stop "Bad pH guesses. Try a smaller calculation step!"

        end if
        c = a
        fc = fa
        d = b - c
        e = d

        ! main loop, exit from middle of the loop
        do while (fb /= 0)
            ! the three current points, a, b, and c, satisfy:
            ! f(x) changes sign between a and b.
            ! abs(f(b)) <= abs(f(a)).
            ! c = previous b, so c might = a.
            ! the next point is chosen from
            ! bisection point, (a+b)/2.
            ! secant point determined by b and c.
            ! inverse quadratic interpolation point determined
            ! by a, b, and c if they are distinct.
            if (fa * fb >=0) then !if (sgn(fa) == sgn(fb)) then
                a = c; fa = fc
                d = b - c; e = d
            end if
            if (abs(fa) < abs(fb)) then
                c = b; b = a; a = c
                fc = fb; fb = fa; fa = fc
            end if
            ! convergence test and possible exit
            m = 0.5 * (a - b)
            tol = 2.0 * eps2 * max(abs(b), 1.0)
            if (abs(m) <= tol .or. fb == 0.0) then
                exit
            end if

            ! choose bisection or interpolation
            if (abs(e) < tol .or. abs(fc) <= abs(fb)) then
                ! bisection
                d = m
                e = m
            else
                ! interpolation
                s = fb / fc
                if (a == c) then
                    ! linear interpolation (secant)
                    p = 2.0 * m * s
                    q = 1.0 - s
                else
                    ! inverse quadratic interpolation
                    q = fc / fa
                    r = fb / fa
                    p = s * (2.0 * m * q * (q - r) - (b - c) * (r - 1.0))
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0)
                end if
                if (p > 0) then
                    q = -q
                else
                    p = -p
                end if
                !is interpolated point acceptable?
                if (2.0 * p < 3.0 * m * q - abs(tol * q) .and. p < abs(0.5 * e * q)) then
                    e = d
                    d = p / q
                else
                    d = m
                    e = m
                end if
            end if
            ! next point
            c = b
            fc = fb
            if (abs(d) > tol) then
                b = b + d
            else
                b= b - sign(tol, b-a) !b = b - sgn(b - a) * tol
            end if
            fb = f(b, ct, te, alke, cond)
        end do

        xr = b

    end subroutine

!#####################################################################

    subroutine phsolbisect(xr, ct, te, alk, cond)
        real(dp), intent(inout) :: xr
        real(dp) ct, te, alk, cond
        real(dp) alke
        real(dp) test, fl, fr
        real(dp) xlo, xup, fu
        integer(i4b) i

        !gp 03-dec-09
        !alke = alk / 50000.0_dp
        alke = alk / 50043.45_dp

        xlo = 3.0_dp
        xup = 13.0_dp
        fl = f(xlo, ct, te, alke, cond)
        fu = f(xup, ct, te, alke, cond)
        if (fl * fu >= 0) then
            !msgbox "bad ph guesses" & chr(13) & "try a smaller calculation step"
            stop "Bad pH guesses. Try a smaller calculation step!"
        end if
        do i = 1, 13
            xr = (xlo + xup) / 2.0_dp
            fr = f(xr, ct, te, alke, cond)
            test = fl * fr
            if (test < 0) then
                xup = xr
            elseif (test > 0) then
                xlo = xr
                fl = fr
            else
                exit
                !msgbox "exact"
            end if
        end do

    end subroutine phsolbisect

    pure function f(ph, ct, te, alk, cond)
        real(dp) f

        real(dp), intent(in) :: ph, ct, te, alk, cond
        real(dp) k1, k2, kw
        real(dp) alp0, alp1, alp2, hh
        real(dp) ta, mu, lam1, lam2

        hh = 10.0_dp ** (-ph)
        ta = te + 273.15_dp
        mu = 0.000016_dp * cond
        lam1 = 10.0_dp ** (-0.5_dp * 1 * sqrt(mu))
        lam2 = 10.0_dp ** (-0.5_dp * 2 ** 2 * sqrt(mu))
        k1 = -356.3094_dp - 0.06091964_dp * ta + 21834.37_dp / ta + 126.8339_dp * log(ta) / log(10.0_dp)
        k1 = k1 - 1684915.0_dp / ta ** 2
        k1 = 10.0_dp ** k1 / lam1 / lam1
        k2 = -107.8871_dp - 0.03252849_dp * ta + 5151.79_dp / ta + 38.92561_dp * log(ta) / log(10.0_dp)
        k2 = k2 - 563713.9_dp / ta ** 2
        k2 = 10.0_dp ** k2 / lam2
        kw = 10.0_dp ** (6.0875_dp - 0.01706_dp * ta - 4470.99_dp / ta) / lam1 / lam1
        alp0 = hh ** 2 / (hh ** 2 + k1 * hh + k1 * k2)
        alp1 = k1 * hh / (hh ** 2 + k1 * hh + k1 * k2)
        alp2 = k1 * k2 / (hh ** 2 + k1 * hh + k1 * k2)
        f = (alp1 + 2.0_dp * alp2) * ct + kw / 10.0_dp ** (-ph) - 10.0_dp ** (-ph) - alk

    end function

!#####################################################################

    subroutine phsolnewton(xr, ct, te, alk, cond)

        ! variables

        real(dp), intent(inout) :: xr
        real(dp) ct, te, alk, cond

        real(dp) alke,lam1,lam2
        integer(i4b) iter
        real(dp) k1, k2, kw
        real(dp) ta !absolute temperature
        real(dp) hh, mu, xnew, ea

        !gp 03-dec-09
        !alke = alk / 50000.0_dp
        alke = alk / 50043.45_dp

        ta = te + 273.15_dp
        xr = 7.0_dp
        mu = 0.000016_dp * cond
        lam1 = 10.0_dp ** (-0.5_dp * 1.0_dp * sqrt(mu))
        lam2 = 10.0_dp ** (-0.5_dp * 2.0_dp ** 2.0_dp * sqrt(mu))
        k1 = 10.0_dp ** (-356.3094_dp - 0.06091964_dp * ta + 21834.37_dp / ta + 126.8339_dp * &
            log(ta) / log(10.0_dp) - 1684915.0_dp / ta ** 2.0_dp) / lam1 / lam1
        k2 = 10.0_dp ** (-107.8871_dp - 0.03252849_dp * ta + 5151.79_dp / ta + 38.92561_dp * &
            log(ta) / log(10.0_dp) - 563713.9_dp / ta ** 2.0) / lam2
        kw = 10.0_dp ** (6.0875_dp - 0.01706_dp * ta - 4470.99_dp / ta) / lam1 / lam1

        iter = 0
        do
            iter = iter + 1
            hh = 10.0_dp ** (-xr)
            xnew = xr - ((k1 * hh + 2.0_dp * k1 * k2) / (hh * hh + k1 * hh + k1 * k2) * &
                ct + kw / hh - hh - alke) / &
                (e * (k1 * hh * ct * (hh * hh + 4.0_dp * k2 * hh + k1 * k2) / &
                (hh * hh + k1 * hh + k1 * k2) ** 2 + kw / hh + hh))
            if (xnew /= 0) then
                ea = abs((xr - xnew) / xnew)
            end if
            xr = xnew
            if ((ea < es) .or. (iter >= imax)) then
                exit
            end if
        end do
    end subroutine phsolnewton

!#####################################################################

    subroutine ph_solver(imethph_0, xr_0, ct_0, te_0, alk_0, cond_0)
        real(dp), intent(inout) :: xr_0
        real(dp), intent(in) :: ct_0, te_0, alk_0, cond_0
        character(len=30), intent(in) :: imethph_0 ! integration method

        select case (imethph_0)
          case ('Newton-Raphson')
            call phsolnewton(xr_0, ct_0, te_0, alk_0, cond_0)
          case ('Bisection')
            call phsolbisect(xr_0, ct_0, te_0, alk_0, cond_0)
          case default ! "Brent"
            call phsolbrent(xr_0, ct_0, te_0, alk_0, cond_0)
        end select

    end subroutine ph_solver


    pure function fnh3(ph, te)
        ! calculate fraction nh3
        real(dp) fnh3
        real(dp), intent(in) :: ph, te

        fnh3 = 1.0_dp / (1.0_dp + 10.0_dp ** (-ph) / &
            (10.0_dp ** ( 0.09018_dp + 2729.92_dp / (te + 273.15_dp))))

    end function fnh3

    !caculate k for tmeperature te
    subroutine chemrates(te, k1, k2, kw, kh, cond)

        real(dp), intent(in) :: te, cond
        real(dp), intent(out) :: k1, k2, kw, kh
        real(dp) ta, mu, lam1, lam2

        ta = te + 273.15_dp
        mu = 0.000016_dp * cond
        lam1 = 10.0_dp ** (-0.5_dp * 1.0_dp * sqrt(mu))
        lam2 = 10.0_dp ** (-0.5_dp * 2.0_dp ** 2 * sqrt(mu))
        k1 = -356.3094_dp - 0.06091964_dp * ta + 21834.37_dp / ta + 126.8339_dp * log(ta) / log(10.0_dp)
        k1 = k1 - 1684915.0_dp / ta ** 2
        k1 = 10.0_dp ** k1 / lam1 / lam1
        k2 = -107.8871_dp - 0.03252849_dp * ta + 5151.79_dp / ta + 38.92561_dp * log(ta) / log(10.0_dp)
        k2 = k2 - 563713.9_dp / ta ** 2
        k2 = 10.0_dp ** k2 / lam2
        kw = 10.0_dp ** (6.0875_dp - 0.01706_dp * ta - 4470.99_dp / ta) / lam1 / lam1
        kh = 10.0_dp ** -(-2385.73_dp / ta - 0.0152642_dp * ta + 14.0184_dp)

    end subroutine chemrates


    function ct(ph, alk, te, cond)

        real(dp) ct
        real(dp), intent(in) :: ph, alk, te, cond
        real(dp) k1, k2, kw, kh
        real(dp) hh, oh, alke, f1, f2

        call chemrates(te, k1, k2, kw, kh, cond)
        hh = 10.0_dp ** -ph
        oh = kw / hh

        !gp 03-dec-09
        !alke = alk / 50000.0_dp
        alke = alk / 50043.45_dp

        f1 = k1 * hh / (hh ** 2 + k1 * hh + k1 * k2)
        f2 = k1 * k2 / (hh ** 2 + k1 * hh + k1 * k2)
        ct = (alke - oh + hh) / (f1 + 2 * f2)

    end function


    subroutine modfp2(xlo1, xup1, xr, co2, te, alk, cond)

        real(dp), intent(in) :: xlo1, xup1
        real(dp) xlo, xup, xr, co2, te, alk, cond

        ! given a function (f) and an initial bracket (xlo, xup),
        ! finds the root (xr) to within a tolerance (tol)
        ! or until a maximum number of iterations (imax) is exceeded.
        ! returns the root along with the actual iterations (iter) and tolerance (tola).
        ! if maximum iterations exceeded, displays error message and terminates.

        integer(i4b) il, iu, iter, imax
        real(dp) xrold, fl, fu, fr
        real(dp) alke, tola, tol


        iter = 0

        !gp 03-dec-09
        !alke = alk / 50000.0_dp
        alke = alk / 50043.45_dp

        imax = 100
        tol = 0.0000001_dp
        tola = 100.0_dp
        xlo=xlo1; xup =xup1

        fl = f2(xlo, co2, te, alke, cond)
        fu = f2(xup, co2, te, alke, cond)
        xr = (xlo + xup) / 2.0_dp

        do
            if (iter > imax) then
                ! msgbox "maximum iterations exceeded"
                stop 'maximum iterations exceeded'
            end if
            if (tola < tol) exit
            xrold = xr
            xr = xup - fu * (xlo - xup) / (fl - fu)
            fr = f2(xr, co2, te, alke, cond)
            iter = iter + 1
            if (xr /= 0) then
                tola = abs(xr - xrold)
            end if
            if (fl * fr < 0) then
                xup = xr
                fu = f2(xup, co2, te, alke, cond)
                iu = 0
                il = il + 1
                if (il >= 1) fl = fl / 2.0_dp
            elseif (fl * fr > 0) then
                xlo = xr
                fl = f2(xlo, co2, te, alke, cond)
                il = 0
                iu = iu + 1
                if (iu >= 1) fu = fu / 2.0_dp
            else
                tola = 0.0
            end if
        end do

    end subroutine

    function f2(ph, co2, te, alk, cond)

        real(dp) f2
        real(dp) ph, co2, te, alk, cond
        real(dp) k1, k2, kw, kh
        real(dp) alp0, alp1, alp2, hh

        hh = 10.0_dp ** -ph
        call chemrates(te, k1, k2, kw, kh, cond)
        alp0 = hh ** 2 / (hh ** 2 + k1 * hh + k1 * k2)
        alp1 = k1 * hh / (hh ** 2 + k1 * hh + k1 * k2)
        alp2 = k1 * k2 / (hh ** 2 + k1 * hh + k1 * k2)
        f2 = (alp1 + 2.0_dp * alp2) * co2 / alp0 + kw / hh - hh - alk

    end function

end module class_phsolve
