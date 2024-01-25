!classhydraulics.f90
!
!hydraulics data structure

MODULE Class_Hydraulics
	USE nrtype
	USE Class_RiverTopo							!only reach number and element number used

	IMPLICIT NONE
	PRIVATE  AllocateHydrauArrays, DepthManning
!	PUBLIC Hydraulics_, MakeHydraulics, Geometry_type, reach
		TYPE Hydraulics_type
			REAL(DP) :: b =0, BB=0, xl =0					!Bottom width
			REAL(DP) :: Hweir =0, Bweir =0					!weir height, width
			REAL(DP) :: alp1=0, bet1=0, alp2=0, bet2 =0			!velocity, stage discharge coefficents
			REAL(DP) :: SS1=0, SS2=0, s =0					!side and channel slopes
			REAL(DP) :: nm=0, kaaa=0, Ediff=0 				!nm - manning's n
			REAL(DP) :: Asb =0, Ast=0, Asd=0
			REAL(DP) :: Frsed=0, Frsod=0
			REAL(DP) :: elev1=0, elev2=0					!elevation
			REAL(DP) :: elev=0						!average elevation
			REAL(DP) :: latr=0, lonr=0					!latitude, longitude radius
			REAL(DP) :: SODspec=0, JCH4spec=0, JNH4spec=0, JSRPspec=0

			REAL(DP) :: Q=0, U=0, trav=0, Qpt=0, Qpta=0			!flow rate, velocity
			REAL(DP) :: depth=0, Ac=0, Vol=0, drop=0
			REAL(DP) :: kau =0, ka=0,  os =0
			CHARACTER(LEN=20) :: kaf='Internal'				!Oxygen reaeration equation
			REAL(DP) :: Ep=0, Epout=0
			REAL(DP) :: QCMD=0, EpCMD=0, QptCMD=0, QptaCMD=0

			!gp 16-Jul-08
			!REAL(DP) :: sedThermCond=0, sedThermDiff=0, HsedCM=0, &
			!			HypoExchFrac=0, porosity=0, rhoCpSed=0, botalg0=0		!gp 11-Jan-05
			REAL(DP) :: sedThermCond=0, sedThermDiff=0, HsedCM=0, &
						HypoExchFrac=0, porosity=0, rhoCpSed=0, botalg0=0, SKOP=0

			REAL(DP) :: HypoPoreVol=0		!gp 03-Nov-04
			REAL(DP) :: EhyporheicCMD=0		!gp 15-Nov-04

			!gp 07-Feb-06
			!Since reach specific rates are not compulsory,
			!these values could be NULL_VAL
			!NULL_VAL means "no value"/ lack of value from users 
			!REAL(DP) :: kaaa= NULL_VAL
			REAL(DP) :: vss =	NULL_VAL				!inorganic suspended solids settling vol
			REAL(DP) :: khc	= NULL_VAL					!slow CBOD hydrolysis rate 
			REAL(DP) :: kdcs = NULL_VAL					!slow CBOD oxidation rate
			REAL(DP) :: kdc	= NULL_VAL					!fast CBOD oxidation rate
			REAL(DP) :: khn	= NULL_VAL					!organic N hydrolysis rate
			REAL(DP) :: von	= NULL_VAL					!Organic N settling velocity
			REAL(DP) :: kn = NULL_VAL					!Ammonium nitrification rate 
			REAL(DP) :: ki = NULL_VAL					!Nitrate denitrification
			REAL(DP) :: vdi	= NULL_VAL					!Nitrate sed denitrification transfer coeff
			REAL(DP) :: khp= NULL_VAL					!organic P hydrolysis
			REAL(DP) :: vop= NULL_VAL					!organic P settling velocity
			REAL(DP) :: vip= NULL_VAL					!inorganic P settling velocity
			REAL(DP) :: kga= NULL_VAL					!Phytoplankton MAX growth rate
			REAL(DP) :: krea= NULL_VAL					!Phytoplankton respiration rate
			REAL(DP) :: kdea= NULL_VAL					!Phytoplankton death rate
			REAL(DP) :: ksn= NULL_VAL					!Phytoplankton N half-sat
			REAL(DP) :: ksp= NULL_VAL					!Phytoplankton P half-sat
			REAL(DP) :: Isat= NULL_VAL					!Phytoplankton light sat
			REAL(DP) :: khnx= NULL_VAL					!Phytoplankton ammonia preference
			REAL(DP) :: va= NULL_VAL					!Phytoplankton settling velocity
			!REAL(DP) :: botalg0= NULL_VAL				!Bottom plant initial bionass
			REAL(DP) :: kgaF= NULL_VAL					!Bottom plant MAX growth rate
			REAL(DP) :: abmax= NULL_VAL					!Bottom plant first-order carrying capacity

			!gp 03-Apr-08
			!REAL(DP) :: kreaF= NULL_VAL					!Bottom plant respiration rate
			REAL(DP) :: krea1F= NULL_VAL					!Bottom plant basal respiration rate
			REAL(DP) :: krea2F= NULL_VAL					!Bottom plant photo respiration rate

			REAL(DP) :: kexaF= NULL_VAL					!Bottom plant excretion rate
			REAL(DP) :: kdeaF= NULL_VAL					!Bottom plant death rate
			REAL(DP) :: ksnF= NULL_VAL					!Bottom plant external N half-sat
			REAL(DP) :: kspF= NULL_VAL					!Bottom plant external P half-sat
			REAL(DP) :: IsatF= NULL_VAL					!Bottom plant light sat
			REAL(DP) :: khnxF= NULL_VAL					!Bottom plant ammonia preference
			REAL(DP) :: NINbmin= NULL_VAL				!Bottom plant subistence quota for N
			REAL(DP) :: NIPbmin= NULL_VAL				!Bottom plant subistence quota for P
			REAL(DP) :: NINbupmax= NULL_VAL				!Bottom plant max uptake rate for N
			REAL(DP) :: NIPbupmax= NULL_VAL				!Bottom plant max uptake rate for P
			REAL(DP) :: KqN= NULL_VAL					!Bottom plant internal N half-sat
			REAL(DP) :: KqP= NULL_VAL					!Bottom plant internal P half-sat
			REAL(DP) :: NUpWCfrac= NULL_VAL				!Bottom plant N uptake fraction from water column
			REAL(DP) :: PUpWCfrac= NULL_VAL				!Bottom plant P uptake fraction from water column
			REAL(DP) :: kdt= NULL_VAL					!POM dissolution rate 
			REAL(DP) :: vdt= NULL_VAL					!POM settling velocity
			REAL(DP) :: kpath= NULL_VAL					!pathogen dieoff rate 
			REAL(DP) :: vpath= NULL_VAL					!pathogen settling velocity
			REAL(DP) :: apath= NULL_VAL					!pathogen light alpha 
			REAL(DP) :: kgaH= NULL_VAL					!Hyporheic heterotrophs MAX growth rate
			REAL(DP) :: kscH= NULL_VAL					!Hyporheic heterotrophs CBOD half-sat
			REAL(DP) :: kinhcH= NULL_VAL				!Hyporheic heterotrophs O2 inhibition
			REAL(DP) :: kreaH= NULL_VAL					!Hyporheic heterotrophs respiration rate
			REAL(DP) :: kdeaH= NULL_VAL					!Hyporheic heterotrophs death rate
			REAL(DP) :: ksnH= NULL_VAL					!Hyporheic heterotrophs N half-sat
			REAL(DP) :: kspH= NULL_VAL					!Hyporheic heterotrophs P half-sat
			REAL(DP) :: khnxH= NULL_VAL					!Hyporheic heterotrophs ammonia preference
			REAL(DP) :: ahmax= NULL_VAL					!Hyporheic heterotrophs first-order carrying capacity
			REAL(DP) :: kgen= NULL_VAL					!generic constituent dissolution rate 
			REAL(DP) :: vgen= NULL_VAL					!generic constituent settling velocity

			!gp 21-Nov-06 initial conditions of state variables
			REAL(DP) :: Te_ini= NULL_VAL				
			REAL(DP) :: c01_ini= NULL_VAL				
			REAL(DP) :: c02_ini= NULL_VAL				
			REAL(DP) :: c03_ini= NULL_VAL				
			REAL(DP) :: c04_ini= NULL_VAL				
			REAL(DP) :: c05_ini= NULL_VAL				
			REAL(DP) :: c06_ini= NULL_VAL				
			REAL(DP) :: c07_ini= NULL_VAL				
			REAL(DP) :: c08_ini= NULL_VAL				
			REAL(DP) :: c09_ini= NULL_VAL				
			REAL(DP) :: c10_ini= NULL_VAL				
			REAL(DP) :: c11_ini= NULL_VAL				
			REAL(DP) :: c12_ini= NULL_VAL				
			REAL(DP) :: c13_ini= NULL_VAL				
			REAL(DP) :: c14_ini= NULL_VAL				
			REAL(DP) :: c15_ini= NULL_VAL				
			REAL(DP) :: pH_ini= NULL_VAL				
			REAL(DP) :: c17_ini= NULL_VAL				
			REAL(DP) :: NINb_ini= NULL_VAL				
			REAL(DP) :: NIPb_ini= NULL_VAL				

		END TYPE Hydraulics_type

		TYPE RiverHydraulics_type
			TYPE(Hydraulics_type), DIMENSION(:), POINTER :: reach
			INTEGER(I4B) :: flag =1
		END TYPE RiverHydraulics_type 		

!		TYPE(RiverHydraulics_type) hydrau

!			TYPE(Geometry_type),	 DIMENSION(:), POINTER :: reach
!			TYPE(Hydraulics_type), DIMENSION(:), POINTER :: hydrau

!			CHARACTER(LEN=20) kai		

!		TYPE(Geometry_type), ALLOCATABLE :: reach(:)
!		TYPE(Hydraulics_type), ALLOCATABLE :: hydrau(:

		CONTAINS
		!/* PUBLIC FUNCTIONS */

		!gp 07-Feb-06
		!FUNCTION Hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, Q, BB, &
		!			SS1, SS2, s, nm, alp1, bet1, alp2, bet2, Ediff, kaaa, &
		!			Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, Bweir, &
		!			sedThermCond, sedThermDiff, HsedCM, &
		!			HypoExchFrac, porosity, rhoCpSed, botalg0) RESULT(hydrau)	!gp 11-Jan-05

		!gp 21-Nov-06
		!FUNCTION Hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, Q, BB, &
		!			SS1, SS2, s, nm, alp1, bet1, alp2, bet2, Ediff, &
		!			Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, Bweir, &
		!			sedThermCond, sedThermDiff, HsedCM, &
		!			HypoExchFrac, porosity, rhoCpSed, &
		!			kaaa, vss, khc, &
		!			kdcs, kdc, khn, &
		!			von, kn , ki , &
		!			vdi, khp, vop, &
		!			vip, kga, krea, &
		!			kdea, ksn, ksp, &
		!			Isat, khnx, va, &
		!			botalg0, kgaF, abmax, &
		!			kreaF, kexaF, kdeaF, &
		!			ksnF, kspF, IsatF, &
		!			khnxF, NINbmin, NIPbmin, &
		!			NINbupmax, NIPbupmax, KqN, &
		!			KqP, NUpWCfrac, PUpWCfrac, &
		!			kdt, vdt, kpath, &
		!			vpath, apath, kgaH, &
		!			kscH, kinhcH, kreaH, &
		!			kdeaH, ksnH, kspH, &
		!			khnxH, ahmax, &
		!			kgen, vgen) RESULT(hydrau)

		!gp 03-Apr-08
		!FUNCTION Hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, Q, BB, &
		!			SS1, SS2, s, nm, alp1, bet1, alp2, bet2, Ediff, &
		!			Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, Bweir, &
		!			sedThermCond, sedThermDiff, HsedCM, &
		!			HypoExchFrac, porosity, rhoCpSed, &
		!			kaaa, vss, khc, &
		!			kdcs, kdc, khn, &
		!			von, kn , ki , &
		!			vdi, khp, vop, &
		!			vip, kga, krea, &
		!			kdea, ksn, ksp, &
		!			Isat, khnx, va, &
		!			kgaF, abmax, &
		!			kreaF, kexaF, kdeaF, &
		!			ksnF, kspF, IsatF, &
		!			khnxF, NINbmin, NIPbmin, &
		!			NINbupmax, NIPbupmax, KqN, &
		!			KqP, NUpWCfrac, PUpWCfrac, &
		!			kdt, vdt, kpath, &
		!			vpath, apath, kgaH, &
		!			kscH, kinhcH, kreaH, &
		!			kdeaH, ksnH, kspH, &
		!			khnxH, ahmax, &
		!			kgen, vgen, &
		!			Te_ini, c01_ini, c02_ini, c03_ini, &
		!			c04_ini, c05_ini, c06_ini, &
		!			c07_ini, c08_ini, c09_ini, &
		!			c10_ini, c11_ini, c12_ini, &
		!			c13_ini, c14_ini, c15_ini, &
		!			pH_ini, c17_ini, NINb_ini, NIPb_ini) RESULT(hydrau)

		!gp 16-Jul-08
		!FUNCTION Hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, Q, BB, &
		!			SS1, SS2, s, nm, alp1, bet1, alp2, bet2, Ediff, &
		!			Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, Bweir, &
		!			sedThermCond, sedThermDiff, HsedCM, &
		!			HypoExchFrac, porosity, rhoCpSed, &
		!			kaaa, vss, khc, &
		!			kdcs, kdc, khn, &
		!			von, kn , ki , &
		!			vdi, khp, vop, &
		!			vip, kga, krea, &
		!			kdea, ksn, ksp, &
		!			Isat, khnx, va, &
		!			kgaF, abmax, &
		!			krea1F, krea2F, kexaF, kdeaF, &
		!			ksnF, kspF, IsatF, &
		!			khnxF, NINbmin, NIPbmin, &
		!			NINbupmax, NIPbupmax, KqN, &
		!			KqP, NUpWCfrac, PUpWCfrac, &
		!			kdt, vdt, kpath, &
		!			vpath, apath, kgaH, &
		!			kscH, kinhcH, kreaH, &
		!			kdeaH, ksnH, kspH, &
		!			khnxH, ahmax, &
		!			kgen, vgen, &
		!			Te_ini, c01_ini, c02_ini, c03_ini, &
		!			c04_ini, c05_ini, c06_ini, &
		!			c07_ini, c08_ini, c09_ini, &
		!			c10_ini, c11_ini, c12_ini, &
		!			c13_ini, c14_ini, c15_ini, &
		!			pH_ini, c17_ini, NINb_ini, NIPb_ini) RESULT(hydrau)
		FUNCTION Hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, Q, BB, &
					SS1, SS2, s, nm, alp1, bet1, alp2, bet2, Ediff, &
					Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, Bweir, &
					sedThermCond, sedThermDiff, HsedCM, &
					HypoExchFrac, porosity, rhoCpSed, SKOP, &
					kaaa, vss, khc, &
					kdcs, kdc, khn, &
					von, kn , ki , &
					vdi, khp, vop, &
					vip, kga, krea, &
					kdea, ksn, ksp, &
					Isat, khnx, va, &
					kgaF, abmax, &
					krea1F, krea2F, kexaF, kdeaF, &
					ksnF, kspF, IsatF, &
					khnxF, NINbmin, NIPbmin, &
					NINbupmax, NIPbupmax, KqN, &
					KqP, NUpWCfrac, PUpWCfrac, &
					kdt, vdt, kpath, &
					vpath, apath, kgaH, &
					kscH, kinhcH, kreaH, &
					kdeaH, ksnH, kspH, &
					khnxH, ahmax, &
					kgen, vgen, &
					Te_ini, c01_ini, c02_ini, c03_ini, &
					c04_ini, c05_ini, c06_ini, &
					c07_ini, c08_ini, c09_ini, &
					c10_ini, c11_ini, c12_ini, &
					c13_ini, c14_ini, c15_ini, &
					pH_ini, c17_ini, NINb_ini, NIPb_ini) RESULT(hydrau)

			TYPE(RiverHydraulics_type) hydrau
			INTEGER(I4B), INTENT(IN) :: nr
			REAL(DP), INTENT(IN) :: xrdn(0:), elev1(0:), elev2(0:), latd(0:), latm(0:)

			!gp 07-Feb-06
			!REAL(DP), INTENT(IN) :: lats(0:), lond(0:), lonm(0:), lons(0:), Q(0:), BB(0:), &
			!			SS1(0:), SS2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
			!			alp2(0:), bet2(0:), Ediff(0:), kaaa(0:), Frsed(0:), &
			!			Frsod(0:), SODspec(0:), JCH4spec(0:), JNH4spec(0:), &
			!			JSRPspec(0:), Hweir(0:), Bweir(0:), &
			!			sedThermCond(0:), sedThermDiff(0:), HsedCM(0:), &
			!			HypoExchFrac(0:), porosity(0:), rhoCpSed(0:), botalg0(0:)  !gp 11-Jan-05

			!gp 21-Nov-06
			!REAL(DP), INTENT(IN) :: lats(0:), lond(0:), lonm(0:), lons(0:), Q(0:), BB(0:), &
			!			SS1(0:), SS2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
			!			alp2(0:), bet2(0:), Ediff(0:), Frsed(0:), &
			!			Frsod(0:), SODspec(0:), JCH4spec(0:), JNH4spec(0:), &
			!			JSRPspec(0:), Hweir(0:), Bweir(0:), &
			!			sedThermCond(0:), sedThermDiff(0:), HsedCM(0:), &
			!			HypoExchFrac(0:), porosity(0:), rhoCpSed(0:), &
			!			kaaa(0:), vss(0:), khc(0:), &
			!			kdcs(0:), kdc(0:), khn(0:), &
			!			von(0:), kn(0:) , ki(0:) , &
			!			vdi(0:), khp(0:), vop(0:), &
			!			vip(0:), kga(0:), krea(0:), &
			!			kdea(0:), ksn(0:), ksp(0:), &
			!			Isat(0:), khnx(0:), va(0:), &
			!			botalg0(0:), kgaF(0:), abmax(0:), &
			!			kreaF(0:), kexaF(0:), kdeaF(0:), &
			!			ksnF(0:), kspF(0:), IsatF(0:), &
			!			khnxF(0:), NINbmin(0:), NIPbmin(0:), &
			!			NINbupmax(0:), NIPbupmax(0:), KqN(0:), &
			!			KqP(0:), NUpWCfrac(0:), PUpWCfrac(0:), &
			!			kdt(0:), vdt(0:), kpath(0:), &
			!			vpath(0:), apath(0:), kgaH(0:), &
			!			kscH(0:), kinhcH(0:), kreaH(0:), &
			!			kdeaH(0:), ksnH(0:), kspH(0:), &
			!			khnxH(0:), ahmax(0:), &
			!			kgen(0:), vgen(0:)

			!gp 03-Apr-08
			!REAL(DP), INTENT(IN) :: lats(0:), lond(0:), lonm(0:), lons(0:), Q(0:), BB(0:), &
			!			SS1(0:), SS2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
			!			alp2(0:), bet2(0:), Ediff(0:), Frsed(0:), &
			!			Frsod(0:), SODspec(0:), JCH4spec(0:), JNH4spec(0:), &
			!			JSRPspec(0:), Hweir(0:), Bweir(0:), &
			!			sedThermCond(0:), sedThermDiff(0:), HsedCM(0:), &
			!			HypoExchFrac(0:), porosity(0:), rhoCpSed(0:), &
			!			kaaa(0:), vss(0:), khc(0:), &
			!			kdcs(0:), kdc(0:), khn(0:), &
			!			von(0:), kn(0:) , ki(0:) , &
			!			vdi(0:), khp(0:), vop(0:), &
			!			vip(0:), kga(0:), krea(0:), &
			!			kdea(0:), ksn(0:), ksp(0:), &
			!			Isat(0:), khnx(0:), va(0:), &
			!			kgaF(0:), abmax(0:), &
			!			kreaF(0:), kexaF(0:), kdeaF(0:), &
			!			ksnF(0:), kspF(0:), IsatF(0:), &
			!			khnxF(0:), NINbmin(0:), NIPbmin(0:), &
			!			NINbupmax(0:), NIPbupmax(0:), KqN(0:), &
			!			KqP(0:), NUpWCfrac(0:), PUpWCfrac(0:), &
			!			kdt(0:), vdt(0:), kpath(0:), &
			!			vpath(0:), apath(0:), kgaH(0:), &
			!			kscH(0:), kinhcH(0:), kreaH(0:), &
			!			kdeaH(0:), ksnH(0:), kspH(0:), &
			!			khnxH(0:), ahmax(0:), &
			!			kgen(0:), vgen(0:), &
			!			Te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
			!			c04_ini(0:), c05_ini(0:), c06_ini(0:), &
			!			c07_ini(0:), c08_ini(0:), c09_ini(0:), &
			!			c10_ini(0:), c11_ini(0:), c12_ini(0:), &
			!			c13_ini(0:), c14_ini(0:), c15_ini(0:), &
			!			pH_ini(0:), c17_ini(0:), NINb_ini(0:), NIPb_ini(0:)

			!gp 16-Jul-08
			!REAL(DP), INTENT(IN) :: lats(0:), lond(0:), lonm(0:), lons(0:), Q(0:), BB(0:), &
			!			SS1(0:), SS2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
			!			alp2(0:), bet2(0:), Ediff(0:), Frsed(0:), &
			!			Frsod(0:), SODspec(0:), JCH4spec(0:), JNH4spec(0:), &
			!			JSRPspec(0:), Hweir(0:), Bweir(0:), &
			!			sedThermCond(0:), sedThermDiff(0:), HsedCM(0:), &
			!			HypoExchFrac(0:), porosity(0:), rhoCpSed(0:), &
			!			kaaa(0:), vss(0:), khc(0:), &
			!			kdcs(0:), kdc(0:), khn(0:), &
			!			von(0:), kn(0:) , ki(0:) , &
			!			vdi(0:), khp(0:), vop(0:), &
			!			vip(0:), kga(0:), krea(0:), &
			!			kdea(0:), ksn(0:), ksp(0:), &
			!			Isat(0:), khnx(0:), va(0:), &
			!			kgaF(0:), abmax(0:), &
			!			krea1F(0:), krea2F(0:), kexaF(0:), kdeaF(0:), &
			!			ksnF(0:), kspF(0:), IsatF(0:), &
			!			khnxF(0:), NINbmin(0:), NIPbmin(0:), &
			!			NINbupmax(0:), NIPbupmax(0:), KqN(0:), &
			!			KqP(0:), NUpWCfrac(0:), PUpWCfrac(0:), &
			!			kdt(0:), vdt(0:), kpath(0:), &
			!			vpath(0:), apath(0:), kgaH(0:), &
			!			kscH(0:), kinhcH(0:), kreaH(0:), &
			!			kdeaH(0:), ksnH(0:), kspH(0:), &
			!			khnxH(0:), ahmax(0:), &
			!			kgen(0:), vgen(0:), &
			!			Te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
			!			c04_ini(0:), c05_ini(0:), c06_ini(0:), &
			!			c07_ini(0:), c08_ini(0:), c09_ini(0:), &
			!			c10_ini(0:), c11_ini(0:), c12_ini(0:), &
			!			c13_ini(0:), c14_ini(0:), c15_ini(0:), &
			!			pH_ini(0:), c17_ini(0:), NINb_ini(0:), NIPb_ini(0:)
			REAL(DP), INTENT(IN) :: lats(0:), lond(0:), lonm(0:), lons(0:), Q(0:), BB(0:), &
						SS1(0:), SS2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
						alp2(0:), bet2(0:), Ediff(0:), Frsed(0:), &
						Frsod(0:), SODspec(0:), JCH4spec(0:), JNH4spec(0:), &
						JSRPspec(0:), Hweir(0:), Bweir(0:), &
						sedThermCond(0:), sedThermDiff(0:), HsedCM(0:), &
						HypoExchFrac(0:), porosity(0:), rhoCpSed(0:), SKOP(0:), &
						kaaa(0:), vss(0:), khc(0:), &
						kdcs(0:), kdc(0:), khn(0:), &
						von(0:), kn(0:) , ki(0:) , &
						vdi(0:), khp(0:), vop(0:), &
						vip(0:), kga(0:), krea(0:), &
						kdea(0:), ksn(0:), ksp(0:), &
						Isat(0:), khnx(0:), va(0:), &
						kgaF(0:), abmax(0:), &
						krea1F(0:), krea2F(0:), kexaF(0:), kdeaF(0:), &
						ksnF(0:), kspF(0:), IsatF(0:), &
						khnxF(0:), NINbmin(0:), NIPbmin(0:), &
						NINbupmax(0:), NIPbupmax(0:), KqN(0:), &
						KqP(0:), NUpWCfrac(0:), PUpWCfrac(0:), &
						kdt(0:), vdt(0:), kpath(0:), &
						vpath(0:), apath(0:), kgaH(0:), &
						kscH(0:), kinhcH(0:), kreaH(0:), &
						kdeaH(0:), ksnH(0:), kspH(0:), &
						khnxH(0:), ahmax(0:), &
						kgen(0:), vgen(0:), &
						Te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
						c04_ini(0:), c05_ini(0:), c06_ini(0:), &
						c07_ini(0:), c08_ini(0:), c09_ini(0:), &
						c10_ini(0:), c11_ini(0:), c12_ini(0:), &
						c13_ini(0:), c14_ini(0:), c15_ini(0:), &
						pH_ini(0:), c17_ini(0:), NINb_ini(0:), NIPb_ini(0:)

			INTEGER(I4B) i, j

			CALL AllocateHydrauArrays(nr, hydrau)
			
			hydrau%reach(0)%Q=Q(0)* 86400.0_DP

			hydrau%reach(0)%nm= nm(0)

			IF (xrdn(2)>xrdn(1)) THEN
				hydrau%flag=1
			ELSE
				hydrau%flag=-1
			END IF

			DO i=0, nr
				!hydrau%reach(i)%Q=Q(i)* 86400.0
				IF (i>0) hydrau%reach(i)%xl = ABS(xrdn(i)-xrdn(i-1))*1000.0_DP
				hydrau%reach(i)%elev1 = elev1(i)
				hydrau%reach(i)%elev2 = elev2(i)
				hydrau%reach(i)%elev = (elev1(i)+elev2(i))/2
				hydrau%reach(i)%SS1   = SS1(i)
				hydrau%reach(i)%SS2   = SS2(i)
				hydrau%reach(i)%s     = s(i)
				IF (Frsod(i) <= 0) THEN				!fraction of bottom SOD
					hydrau%reach(i)%Frsod = 0.00001_DP
				ELSE
					hydrau%reach(i)%Frsod = Frsod(i)			
				END IF
				
				IF (Frsed(i) <= 0) THEN				!fraction of bottom algae coverage
					hydrau%reach(i)%Frsed = 0.00001_DP
				ELSE
					hydrau%reach(i)%Frsed = Frsed(i)			
				END IF
				hydrau%reach(i)%latr  = latd(i) + (latm(i)+lats(i)/60.0_DP)/60.0_DP
				hydrau%reach(i)%lonr  = lond(i) + (lonm(i)+lons(i)/60.0_DP)/60.0_DP
				hydrau%reach(i)%nm		 = nm(i)
				hydrau%reach(i)%Hweir = Hweir(i)
				hydrau%reach(i)%Bweir = Bweir(i)
				IF (((Hweir(i) > 0) .AND. (Bweir(i) <= 0)) .OR. &
					((Bweir(i) > 0) .AND. (Hweir(i) <= 0))) THEN
					STOP 'To include a weir, both height and width must be greater than zero'
			  End IF
				hydrau%reach(i)%BB    = BB(i)
				hydrau%reach(i)%b			= BB(i)

				IF (((Hweir(i) > 0).AND. (Bweir(i) > 0)) .AND. (BB(i) <= 0)) THEN
					STOP 'If a weir is used, the width of the reach upstream on the weir must be entered'
				END IF

				hydrau%reach(i)%alp1  = alp1(i)
				hydrau%reach(i)%alp2  = alp2(i)
				hydrau%reach(i)%bet1	 = bet1(i)
				hydrau%reach(i)%bet2  = bet2 (i)
				hydrau%reach(i)%Ediff = Ediff(i)		!prescribed dispersion rate

				!gp 07-Feb-06
				!hydrau%reach(i)%kaaa  = kaaa(i)		  	!prescribed reaeration rate

				hydrau%reach(i)%SODspec = SODspec(i)  		!prescribed SOD
				hydrau%reach(i)%JCH4spec=JCH4spec(i)
				hydrau%reach(i)%JNH4spec=JNH4spec(i)
				hydrau%reach(i)%JSRPspec=JSRPspec(i)

				hydrau%reach(i)%sedThermCond = sedThermCond(i)	!gp 20-Oct-04
				hydrau%reach(i)%sedThermDiff = sedThermDiff(i)	!gp 20-Oct-04
				hydrau%reach(i)%HsedCM = HsedCM(i)				!gp 20-Oct-04
				hydrau%reach(i)%HypoExchFrac = HypoExchFrac(i)	!gp 20-Oct-04
				hydrau%reach(i)%porosity = porosity(i)			!gp 20-Oct-04
				hydrau%reach(i)%rhoCpSed = rhoCpSed(i)			!gp 20-Oct-04

				!gp 16-Jul-08
				hydrau%reach(i)%SKOP = SKOP(i)

				!gp 07-Feb-06
				!hydrau%reach(i)%botalg0 = botalg0(i)			!gp 11-Jan-05
				hydrau%reach(i)%kaaa  = kaaa(i)		  			!prescribed reaeration rate
				hydrau%reach(i)%vss = vss(i)					!inorganic suspended solids settling vol
				hydrau%reach(i)%khc	= khc(i)					!slow CBOD hydrolysis rate 
				hydrau%reach(i)%kdcs = kdcs(i)					!slow CBOD oxidation rate
				hydrau%reach(i)%kdc	= kdc(i)					!fast CBOD oxidation rate
				hydrau%reach(i)%khn	= khn(i)					!organic N hydrolysis rate
				hydrau%reach(i)%von	= von(i)					!Organic N settling velocity
				hydrau%reach(i)%kn = kn(i)						!Ammonium nitrification rate 
				hydrau%reach(i)%ki = ki(i)						!Nitrate denitrification
				hydrau%reach(i)%vdi	= vdi(i)					!Nitrate sed denitrification transfer coeff
				hydrau%reach(i)%khp = khp(i)					!organic P hydrolysis
				hydrau%reach(i)%vop = vop(i)					!organic P settling velocity
				hydrau%reach(i)%vip = vip(i)					!inorganic P settling velocity
				hydrau%reach(i)%kga = kga(i)					!Phytoplankton MAX growth rate
				hydrau%reach(i)%krea = krea(i)					!Phytoplankton respiration rate
				hydrau%reach(i)%kdea = kdea(i)					!Phytoplankton death rate
				hydrau%reach(i)%ksn = ksn(i)					!Phytoplankton N half-sat
				hydrau%reach(i)%ksp = ksp(i)					!Phytoplankton P half-sat
				hydrau%reach(i)%Isat = Isat(i)					!Phytoplankton light sat
				hydrau%reach(i)%khnx = khnx(i)					!Phytoplankton ammonia preference
				hydrau%reach(i)%va = va(i)						!Phytoplankton settling velocity
				
				!gp 21-Nov-06
				!hydrau%reach(i)%botalg0 = botalg0(i)			!Bottom plant initial bionass
				
				hydrau%reach(i)%kgaF = kgaF(i)					!Bottom plant MAX growth rate
				hydrau%reach(i)%abmax = abmax(i)				!Bottom plant first-order carrying capacity

				!gp 03-Apr-08
				!hydrau%reach(i)%kreaF = kreaF(i)				!Bottom plant respiration rate
				hydrau%reach(i)%krea1F = krea1F(i)				!Bottom plant basal respiration rate
				hydrau%reach(i)%krea2F = krea2F(i)				!Bottom plant photo respiration rate

				hydrau%reach(i)%kexaF = kexaF(i)				!Bottom plant excretion rate
				hydrau%reach(i)%kdeaF = kdeaF(i)				!Bottom plant death rate
				hydrau%reach(i)%ksnF = ksnF(i)					!Bottom plant external N half-sat
				hydrau%reach(i)%kspF = kspF(i)					!Bottom plant external P half-sat
				hydrau%reach(i)%IsatF = IsatF(i)				!Bottom plant light sat
				hydrau%reach(i)%khnxF = khnxF(i)				!Bottom plant ammonia preference
				hydrau%reach(i)%NINbmin = NINbmin(i)			!Bottom plant subistence quota for N
				hydrau%reach(i)%NIPbmin = NIPbmin(i)			!Bottom plant subistence quota for P
				hydrau%reach(i)%NINbupmax = NINbupmax(i)		!Bottom plant max uptake rate for N
				hydrau%reach(i)%NIPbupmax = NIPbupmax(i)		!Bottom plant max uptake rate for P
				hydrau%reach(i)%KqN = KqN(i)					!Bottom plant internal N half-sat
				hydrau%reach(i)%KqP = KqP(i)					!Bottom plant internal P half-sat
				hydrau%reach(i)%NUpWCfrac = NUpWCfrac(i)		!Bottom plant N uptake fraction from water column
				hydrau%reach(i)%PUpWCfrac = PUpWCfrac(i)		!Bottom plant P uptake fraction from water column
				hydrau%reach(i)%kdt = kdt(i)					!POM dissolution rate 
				hydrau%reach(i)%vdt = vdt(i)					!POM settling velocity
				hydrau%reach(i)%kpath = kpath(i)				!pathogen dieoff rate 
				hydrau%reach(i)%vpath = vpath(i)				!pathogen settling velocity
				hydrau%reach(i)%apath = apath(i)				!pathogen light alpha 
				hydrau%reach(i)%kgaH = kgaH(i)					!Hyporheic heterotrophs MAX growth rate
				hydrau%reach(i)%kscH = kscH(i)					!Hyporheic heterotrophs CBOD half-sat
				hydrau%reach(i)%kinhcH = kinhcH(i)				!Hyporheic heterotrophs O2 inhibition
				hydrau%reach(i)%kreaH = kreaH(i)				!Hyporheic heterotrophs respiration rate
				hydrau%reach(i)%kdeaH = kdeaH(i)				!Hyporheic heterotrophs death rate
				hydrau%reach(i)%ksnH = ksnH(i)					!Hyporheic heterotrophs N half-sat
				hydrau%reach(i)%kspH = kspH(i)					!Hyporheic heterotrophs P half-sat
				hydrau%reach(i)%khnxH = khnxH(i)				!Hyporheic heterotrophs ammonia preference
				hydrau%reach(i)%ahmax = ahmax(i)				!Hyporheic heterotrophs first-order carrying capacity
				hydrau%reach(i)%kgen = kgen(i)					!generic constituent dissolution rate 
				hydrau%reach(i)%vgen = vgen(i)					!generic constituent settling velocity

				!gp 21-Nov-06 initial conditions for each reach
				hydrau%reach(i)%Te_ini = Te_ini(i)				!initial temperature
				hydrau%reach(i)%c01_ini = c01_ini(i)			!initial cond
				hydrau%reach(i)%c02_ini = c02_ini(i)			!initial iss
				hydrau%reach(i)%c03_ini = c03_ini(i)			!initial do
				hydrau%reach(i)%c04_ini = c04_ini(i)			!cbod slow
				hydrau%reach(i)%c05_ini = c05_ini(i)			!cbod fast
				hydrau%reach(i)%c06_ini = c06_ini(i)			!org n
				hydrau%reach(i)%c07_ini = c07_ini(i)			!nh4 n
				hydrau%reach(i)%c08_ini = c08_ini(i)			!no23 n
				hydrau%reach(i)%c09_ini = c09_ini(i)			!org p
				hydrau%reach(i)%c10_ini = c10_ini(i)			!srp
				hydrau%reach(i)%c11_ini = c11_ini(i)			!phyto
				hydrau%reach(i)%c12_ini = c12_ini(i)			!detritus
				hydrau%reach(i)%c13_ini = c13_ini(i)			!pathogen
				hydrau%reach(i)%c14_ini = c14_ini(i)			!generic
				hydrau%reach(i)%c15_ini = c15_ini(i)			!alkalinity
				hydrau%reach(i)%pH_ini = pH_ini(i)				!pH
				hydrau%reach(i)%c17_ini = c17_ini(i)			!bottom algae gD/m^2
				hydrau%reach(i)%NINb_ini = NINb_ini(i)			!bottom algae mgN/gD
				hydrau%reach(i)%NIPb_ini = NIPb_ini(i)			!bottom algae mgP/gD

			END DO
			hydrau%reach(0)%elev = elev2(0)
			hydrau%reach(nr+1)%elev1 = elev2(nr)
			hydrau%reach(0)%xl = hydrau%reach(1)%xl
			hydrau%reach(nr+1)%xl = hydrau%reach(nr)%xl


	END FUNCTION Hydraulics_

	
	SUBROUTINE AllocateHydrauArrays(nr, hydrau)

	INTEGER(I4B), INTENT(IN) :: nr
	TYPE(RiverHydraulics_type), INTENT(INOUT) :: hydrau
	INTEGER(I4B) status
	
	IF (.NOT. ASSOCIATED(hydrau%reach)) THEN
		IF (nr > 0) THEN
			!ALLOCATE (reach(0:nr))
			ALLOCATE (hydrau%reach (0:nr+1), STAT= status)
			IF (status==1) stop 'Class_Hydraulics:AllocateHydrauArrays failed. Insufficient Memory!'
		ELSE
			PRINT *, 'ERROR:element number must be great than 0'
			STOP 'Class_Hydraulics:AllocateHydrauArrays failed'
		END IF
	ELSE
		PRINT *, 'Warning: Hydraulics array can only allocate once. Allocation failed!'
	END IF

	END SUBROUTINE AllocateHydrauArrays

	!gp 17-Nov-04 SUBROUTINE MakeHydraulics(nr, hydrau, downstreamBoundary, kai)
	SUBROUTINE MakeHydraulics(nr, hydrau, downstreamBoundary, kai, geoMethod)		!gp 17-Nov-04
	!gp new sub Hydraulics to separate output and save flows in m**3/s and m**3/day
		IMPLICIT NONE
		INTEGER(I4B) ,INTENT(IN) ::nr
		TYPE(RiverHydraulics_type), INTENT(INOUT) :: hydrau

		!gp 17-Nov-04
		CHARACTER(LEN=30), INTENT(IN) :: geoMethod

		LOGICAL(LGT), INTENT(IN) :: downstreamBoundary

		!gp 10-Feb-05
		!CHARACTER(LEN=20) kai					!Reaeration model
		CHARACTER(LEN=30) kai					!Reaeration model

		INTEGER(I4B) i
		REAL(DP) travel, Width1, Width2, Edif, En
		REAL(DP) Eout, dxint, Btop, ctsiv
		REAL(DP) Ushear, FN, Pwet, Rh, Hd

		!Calculate hydraulics
		travel = 0

		DO i = 0, nr !ne()

			IF ((hydrau%reach(i)%Hweir > 0).AND.(hydrau%reach(i)%Bweir > 0)) THEN
				IF (hydrau%reach(i)%b <=0) THEN
					STOP 'Reaches with weirs must have a width entered!'
				END IF
				hydrau%reach(i)%depth = hydrau%reach(i)%Hweir + (hydrau%reach(i)%Q &
																	/ (1.83_DP * hydrau%reach(i)%Bweir)**(2.0_DP/3.0_DP))
				hydrau%reach(i)%Ac = hydrau%reach(i)%depth * hydrau%reach(i)%b
				hydrau%reach(i)%U = hydrau%reach(i)%Q / hydrau%reach(i)%Ac
				Btop = hydrau%reach(i)%b
				Pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
			ELSEIF (hydrau%reach(i)%nm == 0) THEN
				!gp 17-Nov-04
				!gp hydrau%reach(i)%U = hydrau%reach(i)%alp1 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet1
				!gp hydrau%reach(i)%depth = hydrau%reach(i)%alp2 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet2
				!gp hydrau%reach(i)%Ac = hydrau%reach(i)%Q / hydrau%reach(i)%U	!cross section area
				!gp hydrau%reach(i)%b = hydrau%reach(i)%Ac / hydrau%reach(i)%depth	!width
				!gp Btop = hydrau%reach(i)%b
				!gp Pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
				IF (geoMethod == "Depth") THEN
					hydrau%reach(i)%U = hydrau%reach(i)%alp1 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet1
					hydrau%reach(i)%depth = hydrau%reach(i)%alp2 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet2
					hydrau%reach(i)%Ac = hydrau%reach(i)%Q / hydrau%reach(i)%U	!cross section area
					hydrau%reach(i)%b = hydrau%reach(i)%Ac / hydrau%reach(i)%depth	!calc width from depth
					Btop = hydrau%reach(i)%b
					Pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
				ELSE
					hydrau%reach(i)%U = hydrau%reach(i)%alp1 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet1
					hydrau%reach(i)%b = hydrau%reach(i)%alp2 * hydrau%reach(i)%Q ** hydrau%reach(i)%bet2
					hydrau%reach(i)%Ac = hydrau%reach(i)%Q / hydrau%reach(i)%U	!cross section area
					hydrau%reach(i)%depth = hydrau%reach(i)%Ac / hydrau%reach(i)%b	!calc depth from width
					Btop = hydrau%reach(i)%b
					Pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
				END IF		!gp 17-Nov-04 end new block
			ELSE 	!Trapezoidal channel using Manning's equation
				hydrau%reach(i)%depth = DepthManning(hydrau%reach(i)%Q, hydrau%reach(i)%nm, hydrau%reach(i)%b, &
							hydrau%reach(i)%SS1, hydrau%reach(i)%SS2, hydrau%reach(i)%s)
				hydrau%reach(i)%Ac = (hydrau%reach(i)%b + 0.5_dp * (hydrau%reach(i)%SS1 + &
							hydrau%reach(i)%SS2) * hydrau%reach(i)%depth) * hydrau%reach(i)%depth
				Btop = hydrau%reach(i)%b + (hydrau%reach(i)%SS1 + hydrau%reach(i)%SS2) * &
							hydrau%reach(i)%depth
				hydrau%reach(i)%U = hydrau%reach(i)%Q / hydrau%reach(i)%Ac
				Pwet = hydrau%reach(i)%b+ hydrau%reach(i)%depth * SQRT(hydrau%reach(i)%SS1**2 +1) + &
							hydrau%reach(i)%depth * SQRT(hydrau%reach(i)%SS2**2 +1)
			END IF
			Rh = hydrau%reach(i)%Ac / Pwet		!hydraulic radius
  
			IF (i > 0) THEN
				hydrau%reach(i)%Asb = hydrau%reach(i)%Frsed * hydrau%reach(i)%b * hydrau%reach(i)%xl
				hydrau%reach(i)%Ast = Btop * hydrau%reach(i)%xl
				hydrau%reach(i)%Asd = hydrau%reach(i)%Frsod * hydrau%reach(i)%b * hydrau%reach(i)%xl
				hydrau%reach(i)%Vol = hydrau%reach(i)%Ac * hydrau%reach(i)%xl
				!gp 03-Nov-04 hyporheic pore volume in m^3
				hydrau%reach(i)%HypoPoreVol = hydrau%reach(i)%porosity * hydrau%reach(i)%Ast * hydrau%reach(i)%HsedCM * 0.01_DP     
			END IF
  
			IF (hydrau%reach(i)% Ediff > 0) THEN
				Edif = hydrau%reach(i)%Ediff
			ELSE
				IF (hydrau%reach(i)%s >0) THEN
					Edif = 0.011_dp * hydrau%reach(i)%U ** 2 * hydrau%reach(i)%b ** 2 / hydrau%reach(i)%depth &
								/SQRT(grav * hydrau%reach(i)%depth * hydrau%reach(i)%s)

				ELSE IF (hydrau%reach(i)%nm > 0) THEN
					hydrau%reach(i)%s=(hydrau%reach(i)%nm *hydrau%reach(i)%U/Rh**(2.0/3.0))**2
					Edif = 0.011_dp * hydrau%reach(i)%U ** 2 * hydrau%reach(i)%b ** 2 / hydrau%reach(i)%depth &
								/SQRT(grav * hydrau%reach(i)%depth * hydrau%reach(i)%s)
				END IF
			END IF

			En = hydrau%reach(i)%U * hydrau%reach(i)%xl / 2.0_DP
			IF (En <= Edif) THEN
				Edif = Edif - En
				Eout = Edif
			ELSE
				Edif = 0
				Eout = En
			END IF
    
			dxint = (hydrau%reach(i)%xl + hydrau%reach(i+1)%xl) / 2.0_DP
        
			hydrau%reach(i)%Ep = Edif * hydrau%reach(i)%Ac / dxint
			hydrau%reach(i)%Epout = Eout * hydrau%reach(i)%Ac / dxint
  

			SELECT CASE (kai)
			CASE ("Tsivoglou-Neal", "Thackston-Dawson", &
						"USGS(pool-riffle)", "USGS(channel-control)") 
				IF (hydrau%reach(i)%s == 0) THEN
					STOP "kai requires a non-zero channel slope for segment "
				ELSE
					Pwet = hydrau%reach(i)%b + hydrau%reach(i)%depth * &
									SQRT(hydrau%reach(i)%SS1 ** 2 + 1.0) + hydrau%reach(i)%depth * &
									SQRT(hydrau%reach(i)%SS2 ** 2 + 1.0)
					Rh = hydrau%reach(i)%Ac / Pwet     !Hydraulic radius
					Ushear = SQRT(grav * Rh * hydrau%reach(i)%s)
					Hd = hydrau%reach(i)%Ac / Btop     !Hydraulic depth
					FN = hydrau%reach(i)%U / SQRT(grav * Hd)
				END IF
			END SELECT


			IF (hydrau%reach(i)%kaaa > 0) THEN
				hydrau%reach(i)%kau = hydrau%reach(i)%kaaa
				hydrau%reach(i)%kaf = "Specified"
			ELSEIF (kai == "Internal") THEN
				IF (hydrau%reach(i)%depth < 0.61) THEN
					hydrau%reach(i)%kau = 5.32_DP * hydrau%reach(i)%U ** 0.67_DP / hydrau%reach(i)%depth ** 1.85_DP
					hydrau%reach(i)%kaf = "Owens"
				ELSEIF (hydrau%reach(i)%depth > 3.45 * hydrau%reach(i)%U ** 2.5_dp) THEN
					hydrau%reach(i)%kau = 3.93_DP * hydrau%reach(i)%U ** 0.5_DP / hydrau%reach(i)%depth ** 1.5_DP
					hydrau%reach(i)%kaf = "O'Conn"
				ELSE
					hydrau%reach(i)%kau = 5.026_DP * hydrau%reach(i)%U / hydrau%reach(i)%depth ** 1.67_DP
					hydrau%reach(i)%kaf = "Church"
				END IF
			ELSEIF (kai == "O'Connor-Dobbins") THEN
				hydrau%reach(i)%kau = 3.93_DP * hydrau%reach(i)%U ** 0.5_DP / hydrau%reach(i)%depth ** 1.5_DP
				hydrau%reach(i)%kaf = "O'Conn"
			ELSEIF (kai == "Churchill") THEN
				hydrau%reach(i)%kau = 5.026_DP * hydrau%reach(i)%U / hydrau%reach(i)%depth ** 1.67_DP
				hydrau%reach(i)%kaf = "Church"
			ELSEIF (kai == "Owens-Gibbs") THEN
				hydrau%reach(i)%kau = 5.32_DP * hydrau%reach(i)%U ** 0.67_DP / hydrau%reach(i)%depth ** 1.85_DP
				hydrau%reach(i)%kaf = "Owens"
			ELSEIF (kai == "Tsivoglou-Neal") THEN
		    ctsiv = 0.054_DP
				IF (hydrau%reach(i)%Q < (15.0_DP / 3.281_DP ** 3)) ctsiv = 0.11_DP
				hydrau%reach(i)%kau = 86400.0_DP * ctsiv * hydrau%reach(i)%s * hydrau%reach(i)%U * 3.281_DP
				hydrau%reach(i)%kaf = "Tsivoglou"
			ELSEIF (kai == "Thackston-Dawson") THEN
				hydrau%reach(i)%kau = 86400.0_DP * 0.000025_DP * (1.0_DP + 9.0_DP * FN ** 0.25_DP) * Ushear / &
															hydrau%reach(i)%depth
				hydrau%reach(i)%kaf = "Thack-Dawson"
			ELSEIF (kai == "USGS(pool-riffle)") THEN
				IF (hydrau%reach(i)%Q < 0.556) THEN
					hydrau%reach(i)%kau = 517.0_DP * (hydrau%reach(i)%U * hydrau%reach(i)%s) ** 0.524_DP &
															* hydrau%reach(i)%Q ** (-0.242_DP)
				ELSE
					hydrau%reach(i)%kau = 596.0_DP * (hydrau%reach(i)%U * hydrau%reach(i)%s) ** 0.528_DP &
															* hydrau%reach(i)%Q ** (-0.136_DP)
				END IF
				hydrau%reach(i)%kaf = "Pool-riffle"
			ELSEIF (kai == "USGS(channel-control)") THEN
				IF (hydrau%reach(i)%Q < 0.556) THEN
					hydrau%reach(i)%kau = 88.0_DP * (hydrau%reach(i)%U * hydrau%reach(i)%s) ** 0.313_DP &
															* hydrau%reach(i)%depth ** (-0.353_DP)
				ELSE
					hydrau%reach(i)%kau = 142.0_DP * (hydrau%reach(i)%U * hydrau%reach(i)%s) ** 0.333_DP * &
															hydrau%reach(i)%depth ** (-0.66_DP) * Btop ** (-0.243_DP)
				END IF
				hydrau%reach(i)%kaf = "Channel-control"
		  END IF

			travel = travel + hydrau%reach(i)%Vol / hydrau%reach(i)%Q / 86400.0_DP
			hydrau%reach(i)%trav = travel
		END DO

		DO i = 0, nr - 1
			IF ((hydrau%reach(i)%Hweir > 0.01_dp) .OR. &
				(hydrau%reach(i)%elev2 > hydrau%reach(i+1)%elev1 + 0.01)) THEN
				hydrau%reach(i)%Ep = 0
				hydrau%reach(i)%Epout = 0
				hydrau%reach(i)%drop = hydrau%reach(i)%elev2 + hydrau%reach(i)%depth - &
																hydrau%reach(i+1)%elev1 - hydrau%reach(i+1)%depth
			ELSE
				hydrau%reach(i)%drop = 0
			END IF
		END DO

		!gp 17-Jan-06
		IF (.NOT.downstreamBoundary) THEN
		!IF (.NOT.DB%downstreamBound) THEN

		 hydrau%reach(nr)%Ep = 0
		END IF

		!gp 17-Jan-06
		!WRITE(11,*) downstreamBoundary

		!gp convert flows and bulk dispersion to cubic meters per day (CMD) units
		hydrau%reach(0)%QCMD = hydrau%reach(0)%Q * 86400.0_DP
		hydrau%reach(0)%EpCMD = hydrau%reach(0)%Ep * 86400.0_DP
		DO i = 1, nr
			hydrau%reach(i)%QptCMD = hydrau%reach(i)%Qpt * 86400.0_DP
			hydrau%reach(i)%QptaCMD = hydrau%reach(i)%Qpta * 86400.0_DP
			hydrau%reach(i)%QCMD = hydrau%reach(i)%Q * 86400.0_DP
			hydrau%reach(i)%EpCMD = hydrau%reach(i)%Ep * 86400.0_DP
			!gp 15-Nov-04 hyporheic exchange flow in m^3/day
			hydrau%reach(i)%EhyporheicCMD = hydrau%reach(i)%HypoExchFrac * hydrau%reach(i)%QCMD				!'units of (m3/d)

			!gp 17-Jan-06
			!WRITE(11,*) i, hydrau%reach(i)%Ep

		END DO

	END SUBROUTINE MakeHydraulics

	PURE FUNCTION DepthManning(QQ, nmm, BB, sss1, sss2, ss) RESULT(dep)
		IMPLICIT NONE
		REAL(DP),INTENT(IN) :: QQ, nmm, BB, sss1, sss2, ss
		REAL(DP) hold, ea
		REAL(DP) dep

		DO
			hold = dep
			dep = (QQ * nmm) ** 0.6_DP * (BB + dep * &
								SQRT(sss1 * sss1 + 1.0_dp) + dep * &
								SQRT(sss2 * sss2 + 1.0_dp)) ** 0.4_DP/ &
								(BB + 0.5_DP * (sss1 + sss2) * dep) / ss ** 0.3_DP
			ea = ABS((dep - hold) / dep) * 100.0_DP
			IF (ea < 0.001) EXIT
		END DO
	
	END FUNCTION


END MODULE Class_Hydraulics
