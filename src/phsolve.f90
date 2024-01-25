! phsol.f90 
!	PH solver 
! FUNCTIONS/SUBROUTINES exported from phsol.dll:
!	phsoldll      - subroutine 
!
MODULE Class_Phsolve
	USE nrtype
	IMPLICIT NONE
	PRIVATE f

	CONTAINS

!gp 23-Nov-09
!!#####################################################################
!!12/15/07
!SUBROUTINE pHSolver(IMethpH, pH, cT, Te, Alk, Cond)
!CHARACTER(*), INTENT(IN) :: IMethpH       ! pH solve method
!REAL(DP), INTENT(OUT) :: pH
!REAL(DP), INTENT(IN) :: cT, Te, Alk, Cond
!
!IF (IMethpH == "Newton-Raphson") THEN
!    CALL phsolNewton(pH, cT, Te, ALK, COND)      
!ELSEIF (IMethpH == "Bisection") THEN     !09/23/07
!    CALL pHsolBisect(pH, cT, Te, ALK, COND)
!ELSE   !Brent method
!    CALL pHsolFzerotx(pH, cT, Te, ALK, COND) 
!END IF
!
!END SUBROUTINE
!!#####################################################################

!09/23/07 new pH solver, Brent method, used as default method

!GP 23-Nov-09
!SUBROUTINE pHsolFzerotx(xr, cT, Te, Alk, Cond)
SUBROUTINE pHsolBrent(xr, cT, Te, Alk, Cond)

!FZEROTX  Textbook version of FZERO.
!   x = fzerotx(F,[a,b]) tries to find a zero of F(x) between a and b.
!   F(a) and F(b) must have opposite signs.  fzerotx returns one
!   end point of a small subinterval of [a,b] where F changes sign.

!In&Out varaibles
REAL(DP), INTENT(OUT) :: xr
REAL(DP), INTENT(IN) :: cT, Te, Alk, Cond

!Local variables
REAL(DP) a, b, c, d, e
REAL(DP) fa, fb, fc
REAL(DP) p, q, r, s
REAL(DP) m, tol
REAL(DP) Alke
REAL(DP), PARAMETER::  eps2 = 2.22044604925031E-16

! Initialize
a = 0.0
b = 14.0

!gp 03-Dec-09
!Alke = Alk / 50000.0
Alke = Alk / 50043.45_DP

fa = f(a, cT, Te, Alke, Cond)
fb = f(b, cT, Te, Alke, Cond)
IF (fa * fb >=0) THEN      !If (Sgn(fa) == Sgn(fb)) Then
  !MsgBox "Function must change sign on the interval"
  !End
  
  !gp 23-Nov-09
  !WRITE(LOG_FILE,*) "Error: Bad pH guesses. Try a smaller calculation step!"
  !CALL ABORT("Bad pH guesses. Try a smaller calculation step!")
  STOP "Bad pH guesses. Try a smaller calculation step!"

End If
c = a
fc = fa
d = b - c
e = d

! Main loop, exit from middle of the loop
DO While (fb /= 0)
  ! The three current points, a, b, and c, satisfy:
  !    f(x) changes sign between a and b.
  !    abs(f(b)) <= abs(f(a)).
  !    c = previous b, so c might = a.
  ! The next point is chosen from
  !    Bisection point, (a+b)/2.
  !    Secant point determined by b and c.
  !    Inverse quadratic interpolation point determined
  !    by a, b, and c if they are distinct.
  IF (fa * fb >=0) THEN      !If (Sgn(fa) == Sgn(fb)) Then
    a = c;  fa = fc
    d = b - c;  e = d
  End If
  If (Abs(fa) < Abs(fb)) Then
    c = b;    b = a;    a = c
    fc = fb;  fb = fa;  fa = fc
  End If
  ! Convergence test and possible exit
  m = 0.5 * (a - b)
  tol = 2.0 * eps2 * Max(Abs(b), 1.0)
  If (Abs(m) <= tol .OR. fb == 0.0) Then
    Exit
  End If
  
  ! Choose bisection or interpolation
  If (Abs(e) < tol .OR. Abs(fc) <= Abs(fb)) Then
    ! Bisection
    d = m
    e = m
  Else
    ! Interpolation
    s = fb / fc
    If (a == c) Then
      ! Linear interpolation (secant)
      p = 2.0 * m * s
      q = 1.0 - s
    Else
      ! Inverse quadratic interpolation
      q = fc / fa
      r = fb / fa
      p = s * (2.0 * m * q * (q - r) - (b - c) * (r - 1.0))
      q = (q - 1.0) * (r - 1.0) * (s - 1.0)
    End If
    If (p > 0) Then
      q = -q
    Else
      p = -p
    End If
    !Is interpolated point acceptable?
    If (2.0 * p < 3.0 * m * q - Abs(tol * q) .AND. p < Abs(0.5 * e * q)) Then
      e = d
      d = p / q
    Else
      d = m
      e = m
    End If
  End If
  ! Next point
  c = b
  fc = fb
  If (Abs(d) > tol) Then
    b = b + d
  Else
    b= b - sign(tol, b-a)  !b = b - Sgn(b - a) * tol
  End If
  fb = f(b, cT, Te, Alke, Cond)
END DO

xr = b

END SUBROUTINE

!#####################################################################

	SUBROUTINE phsolBisect(xr, cT, Te, Alk, Cond)
		REAL(DP), INTENT(INOUT) :: xr
		REAL(DP) cT, Te, Alk, Cond
		REAL(DP) Alke
		REAL(DP) test, fl, fr
		REAL(DP) xlo, xup, fu
		INTEGER(I4B) i

		!gp 03-Dec-09
		!Alke = Alk / 50000.0_DP
		Alke = Alk / 50043.45_DP

		xlo = 3.0_DP
		xup = 13.0_DP
		fl = f(xlo, cT, Te, Alke, Cond)
		fu = f(xup, cT, Te, Alke, Cond)
		IF (fl * fu >= 0) THEN
			!MsgBox "Bad pH guesses" & Chr(13) & "Try a smaller calculation step"
			STOP "Bad pH guesses. Try a smaller calculation step!"
		END IF
		DO i = 1, 13
			xr = (xlo + xup) / 2.0_DP
			fr = f(xr, cT, Te, Alke, Cond)
			test = fl * fr
			IF (test < 0) THEN
				xup = xr
			ELSEIF (test > 0) THEN
				xlo = xr
				fl = fr
			Else
				EXIT
				!MsgBox "exact"
			END IF
		END DO

	END SUBROUTINE phsolBisect

	PURE FUNCTION f(pH, cT, Te, Alk, Cond)
		REAL(DP) f

		REAL(DP), INTENT(IN) :: pH, cT, Te, Alk, Cond
		REAL(DP) K1, K2, KW
		REAL(DP) alp0, alp1, alp2, hh
		REAL(DP) Ta, mu, lam1, lam2

		hh = 10.0_dp ** -pH
		Ta = Te + 273.15_DP
		mu = 0.000016_DP * Cond
		lam1 = 10.0_DP ** (-0.5_DP * 1 * SQRT(mu))
		lam2 = 10.0_DP ** (-0.5_DP * 2 ** 2 * SQRT(mu))
		K1 = -356.3094_DP - 0.06091964_DP * Ta + 21834.37_DP / Ta + 126.8339_DP * LOG(Ta) / LOG(10.0_DP)
		K1 = K1 - 1684915.0_DP / Ta ** 2
		K1 = 10.0_DP ** K1 / lam1 / lam1
		K2 = -107.8871_DP - 0.03252849_DP * Ta + 5151.79_DP / Ta + 38.92561_DP * LOG(Ta) / LOG(10.0_DP)
		K2 = K2 - 563713.9_DP / Ta ** 2
		K2 = 10.0_DP ** K2 / lam2
		KW = 10.0_DP ** (6.0875_DP - 0.01706_DP * Ta - 4470.99_DP / Ta) / lam1 / lam1
		alp0 = hh ** 2 / (hh ** 2 + K1 * hh + K1 * K2)
		alp1 = K1 * hh / (hh ** 2 + K1 * hh + K1 * K2)
		alp2 = K1 * K2 / (hh ** 2 + K1 * hh + K1 * K2)
		f = (alp1 + 2.0_DP * alp2) * cT + KW / 10.0_DP ** (-pH) - 10.0_DP ** (-pH) - Alk

	END FUNCTION

!#####################################################################

	SUBROUTINE phsolNewton(xr, cT, Te, Alk, Cond)

		! Variables

		REAL(DP), INTENT(INOUT) :: xr
		REAL(DP) cT, Te, Alk, Cond
		
		REAL(DP) Alke,lam1,lam2
		INTEGER(I4B) iter
		REAL(DP) K1, K2, KW
		REAL(DP) Ta																				!absolute temperature
		REAL(DP) hh, mu, xnew, ea
		
		!gp 03-Dec-09
		!Alke = Alk / 50000.0_DP
		Alke = Alk / 50043.45_DP

		Ta = Te + 273.15_DP
		xr = 7.0_DP
		mu = 0.000016_DP * Cond
		lam1 = 10.0_DP ** (-0.5_DP * 1.0_DP * SQRT(mu))
		lam2 = 10.0_DP ** (-0.5_DP * 2.0_DP ** 2.0_DP * SQRT(mu))
		K1 = 10.0_dp ** (-356.3094_DP - 0.06091964_DP * Ta + 21834.37_DP / Ta + 126.8339_DP * &
					LOG(Ta) / LOG(10.0_dp) - 1684915.0_dp / Ta ** 2.0_DP) / lam1 / lam1
		K2 = 10.0_DP ** (-107.8871_dp - 0.03252849_dp * Ta + 5151.79_DP / Ta + 38.92561_DP * & 
					LOG(Ta) / LOG(10.0_dp) - 563713.9_DP / Ta ** 2.0) / lam2
		KW = 10.0_DP ** (6.0875_dp - 0.01706_dp * Ta - 4470.99_DP / Ta) / lam1 / lam1

		iter = 0
		Do
			iter = iter + 1
			hh = 10.0_DP ** (-xr)
			xnew = xr - ((K1 * hh + 2.0_DP * K1 * K2) / (hh * hh + K1 * hh + K1 * K2) * & 
																								cT + KW / hh - hh - Alke) / &
						(e * (K1 * hh * cT * (hh * hh + 4.0_DP * K2 * hh + K1 * K2) / &
													(hh * hh + K1 * hh + K1 * K2) ** 2 + KW / hh + hh))
			IF (xnew /= 0) THEN
				ea = Abs((xr - xnew) / xnew)
			END IF
			xr = xnew
			IF ((ea < es) .Or. (iter >= imax)) THEN
				EXIT
			END IF
		END DO
	END SUBROUTINE phsolNewton

!#####################################################################

	PURE FUNCTION FNH3(pH, Te)
	! Calculate fraction NH3
		REAL(DP) FNH3
		REAL(DP), INTENT(IN) :: pH, Te

		FNH3 = 1.0_DP / (1.0_DP + 10.0_DP ** (-pH) /  &
					(10.0_DP ** ( 0.09018_DP + 2729.92_DP / (Te + 273.15_DP))))

	END FUNCTION FNH3

	!oxygen satuation concentration function
	!temp- temperature
	!	elev- elevation
	PURE FUNCTION oxsat(temp, elev)
		REAL(DP) oxsat 
		REAL(DP), INTENT(IN) :: temp, elev
		REAL(DP) taa, lnosf

		Taa = Temp + 273.15_DP
		lnosf = -139.34411_DP + 157570.1_DP / Taa - 66423080.0_DP / (Taa * Taa)
		lnosf= lnosf + 12438.0e6_DP / (Taa * Taa * Taa) - 8621949.0e5_DP / (Taa * Taa * Taa * Taa)
		oxsat = Exp(lnosf) * (1.0_DP - 0.0001148_DP * elev)
	END FUNCTION

	!CACULATE k FOR TMEPERATURE Te
	SUBROUTINE ChemRates(Te, K1, K2, KW, Kh, Cond)

	REAL(DP), INTENT(IN) :: Te, Cond
	REAL(DP), INTENT(OUT) :: K1, K2, KW, Kh
	REAL(DP) Ta, mu, lam1, lam2

	Ta = Te + 273.15_DP
	mu = 0.000016_DP * Cond
	lam1 = 10.0_DP ** (-0.5_DP * 1.0_DP * SQRT(mu))
	lam2 = 10.0_DP ** (-0.5_DP * 2.0_DP ** 2 * SQRT(mu))
	K1 = -356.3094_DP - 0.06091964_DP * Ta + 21834.37_DP / Ta + 126.8339_DP * LOG(Ta) / LOG(10.0_DP)
	K1 = K1 - 1684915.0_DP / Ta ** 2
	K1 = 10.0_DP ** K1 / lam1 / lam1
	K2 = -107.8871_DP - 0.03252849_DP * Ta + 5151.79_DP / Ta + 38.92561_DP * LOG(Ta) / LOG(10.0_DP)
	K2 = K2 - 563713.9_dp / Ta ** 2
	K2 = 10.0_DP ** K2 / lam2
	KW = 10.0_DP ** (6.0875_DP - 0.01706_DP * Ta - 4470.99_DP / Ta) / lam1 / lam1
	Kh = 10.0_dp ** -(-2385.73_dp / Ta - 0.0152642_dp * Ta + 14.0184_dp)

	END SUBROUTINE


	FUNCTION cT(pH, Alk, Te, Cond)

	REAL(DP) cT
	REAL(DP), INTENT(IN) :: pH, Alk, Te, Cond
	REAL(DP) K1, K2, KW, Kh
	REAL(DP) hh, OH, Alke, f1, f2

	Call ChemRates(Te, K1, K2, KW, Kh, Cond)
	hh = 10.0_DP ** -pH
	OH = KW / hh

	!gp 03-Dec-09
	!Alke = Alk / 50000.0_DP
	Alke = Alk / 50043.45_DP

	f1 = K1 * hh / (hh ** 2 + K1 * hh + K1 * K2)
	f2 = K1 * K2 / (hh ** 2 + K1 * hh + K1 * K2)
	cT = (Alke - OH + hh) / (f1 + 2 * f2)

	END FUNCTION


	SUBROUTINE ModFP2(xlo1, xup1, xr, cO2, Te, Alk, Cond)

	REAL(DP), INTENT(IN) :: xlo1, xup1
	REAL(DP) xlo, xup, xr, cO2, Te, Alk, Cond

	! Given a function (f) and an initial bracket (xlo, xup),
	! finds the root (xr) to within a tolerance (tol)
	! or until a maximum number of iterations (imax) is exceeded.
	! Returns the root along with the actual iterations (iter) and tolerance (tola).
	! If maximum iterations exceeded, displays error message and terminates.

	INTEGER(I4B) il, iu, iter, imax
	REAL(DP) xrold, fl, fu, fr
	REAL(DP) Alke, tola, tol


	iter = 0

	!gp 03-Dec-09
	!Alke = Alk / 50000.0_DP
	Alke = Alk / 50043.45_DP

	imax = 100
	tol = 0.0000001_DP
	tola = 100.0_DP
	xlo=xlo1; xup =xup1

	fl = f2(xlo, cO2, Te, Alke, Cond)
	fu = f2(xup, cO2, Te, Alke, Cond)
	xr = (xlo + xup) / 2.0_DP

	DO
		If (iter > imax) Then
		!	MsgBox "maximum iterations exceeded"
			STOP 'maximum iterations exceeded'
		End If
		If (tola < tol) EXIT
		xrold = xr
		xr = xup - fu * (xlo - xup) / (fl - fu)
		fr = f2(xr, cO2, Te, Alke, Cond)
		iter = iter + 1
		If (xr /= 0) Then
			tola = Abs(xr - xrold)
		End If
		If (fl * fr < 0) Then
			xup = xr
			fu = f2(xup, cO2, Te, Alke, Cond)
			iu = 0
			il = il + 1
			If (il >= 1) fl = fl / 2.0_DP
		ElseIf (fl * fr > 0) Then
			xlo = xr
			fl = f2(xlo, cO2, Te, Alke, Cond)
			il = 0
			iu = iu + 1
			If (iu >= 1) fu = fu / 2.0_DP
		Else
			tola = 0.0
		End If
	END DO

	End SUBROUTINE

	FUNCTION f2(pH, cO2, Te, Alk, Cond)

	REAL(DP) f2
	REAL(DP) pH, cO2, Te, Alk, Cond
	REAL(DP) K1, K2, KW, Kh
	REAL(DP) alp0, alp1, alp2, hh

	hh = 10.0_DP ** -pH
	Call ChemRates(Te, K1, K2, KW, Kh, Cond)
	alp0 = hh ** 2 / (hh ** 2 + K1 * hh + K1 * K2)
	alp1 = K1 * hh / (hh ** 2 + K1 * hh + K1 * K2)
	alp2 = K1 * K2 / (hh ** 2 + K1 * hh + K1 * K2)
	f2 = (alp1 + 2.0_DP * alp2) * cO2 / alp0 + KW / hh - hh - Alk
            
	End FUNCTION

END MODULE Class_Phsolve