MODULE Class_ReadFile

contains
    SUBROUTINE ReadInputfile(system, hydrau, siteMeteo, HW, DB, stochRate, topo, siteSolar)
        USE nrtype
        USE Class_SystemParams    !, ONLY: SystemParams, SystemParams_
        USE Class_RiverTopo, only: rivertopo_type
        USE Class_Hydraulics
        USE m_meteorology, only: t_meteorology

        USE Class_LightHeat
        USE Class_SourceIn
        USE Class_Rates
        USE m_water_quality, ONLY:t_water_quality
        USE Class_Headwater, ONLY: Headwater_, Headwater_type
        USE Class_downstream
        USE Class_SolarCalc

        IMPLICIT NONE

        TYPE(SystemParams), INTENT(OUT) :: system
        TYPE(RiverHydraulics_type), INTENT(OUT) :: hydrau
        TYPE(t_meteorology), INTENT(OUT) :: siteMeteo
        TYPE(Headwater_type) HW
        TYPE(Downstream_type) DB
        TYPE(Rates_type) stochRate				!stoch, reaction, temperature rate
        TYPE(RiverTopo_type) Topo				!river topology
        TYPE(solar_type) :: siteSolar				!solar radiation

        !gp 16-Jul-08
        !gp 21-Nov-06
        !gp 07-Feb-06
        !gp 28-Oct-04 INTEGER(I4B) i, j, k, status(32)
        !INTEGER(I4B) i, j, k, status(39)		!gp 11-Jan-05
        !INTEGER(I4B) i, j, k, status(92)		!gp 11-Jan-05
        !INTEGER(I4B) i, j, k, status(112)		!gp 21-Nov-06
        !INTEGER(I4B) i, j, k, status(113)		!gp 03-Apr-08
        INTEGER(I4B) i, j, k, status(114)		!gp 16-Jul-08

        INTEGER(I4B) nRch, nElem 				!number of reach, num of element
        INTEGER(I4B) nptIn, ndiffIn				!num of point load/src, num of non-point
        !system parameters

        !GP 23-Nov-09
        !CHARACTER(LEN=30) BASINNAME, FILENAME, PATH, TITLE, timezone
        CHARACTER(LEN=30) BASINNAME, FILENAME, PATH, TITLE
        REAL(DP) timezone

        CHARACTER(LEN=30) :: IMeth				!Integration method
        CHARACTER(LEN=30) :: IMethpH				!pH method
        CHARACTER(LEN=30) :: simHyporheicWQ			!gp 28-Oct-04 Yes or No to simulate hyporheic water quality and flux
        CHARACTER(LEN=30) :: showDielResults		!gp 03-Feb-05 Yes or No (only used in Excel VBA)
        CHARACTER(LEN=30) :: stateVariables			!gp 03-Feb-05 All or only Temperature

        !gp 11-Jan-06
        CHARACTER(LEN=30) :: calcSedFlux			!gp 11-Jan-06 Yes or No

        !gp 11-Jan-06
        CHARACTER(LEN=30) :: simAlk					!gp 26-Oct-07 Yes or No

        !gp 24-Jun-09
        CHARACTER(LEN=30) :: writeDynamic			!Yes or No

        REAL(DP) year, month, day
        REAL(DP) dtuser, tf

        !reach data
        CHARACTER(LEN=30) :: geoMethod				!gp 17-Nov-04 Depth or Width for coeff/exponents in col T and U of 'Reach' sheet
        REAL(DP), ALLOCATABLE :: xrdn(:), elev1(:), elev2(:), latd(:), latm(:)

        !gp 07-Feb-07
        !REAL(DP), ALLOCATABLE :: lats(:), lond(:), lonm(:), lons(:), Q(:), BB(:), &
        !				SS1(:), SS2(:), s(:), nm(:) , alp1(:), bet1(:), &
        !				alp2(:), bet2(:), Ediff(:), kaaa(:), Frsed(:), &
        !				Frsod(:), SODspec(:), JCH4spec(:), JNH4spec(:), &
        !				JSRPspec(:), Hweir(:), Bweir(:), &
        !				sedThermCond(:), sedThermDiff(:), HsedCM(:), &
        !				HypoExchFrac(:), porosity(:), rhoCpSed(:), botalg0(:)		!gp 11-Jan-05

        !gp 03-Apr-08
        !REAL(DP), ALLOCATABLE :: lats(:), lond(:), lonm(:), lons(:), Q(:), BB(:), &
        !				SS1(:), SS2(:), s(:), nm(:) , alp1(:), bet1(:), &
        !				alp2(:), bet2(:), Ediff(:), Frsed(:), &
        !				Frsod(:), SODspec(:), JCH4spec(:), JNH4spec(:), &
        !				JSRPspec(:), Hweir(:), Bweir(:), &
        !				sedThermCond(:), sedThermDiff(:), HsedCM(:), &
        !				HypoExchFrac(:), porosity(:), rhoCpSed(:), &
        !				kaaa(:), vss_rch(:), khc_rch(:), &
        !				kdcs_rch(:), kdc_rch(:), khn_rch(:), &
        !				von_rch(:), kn_rch(:) , ki_rch(:) , &
        !				vdi_rch(:), khp_rch(:), vop_rch(:), &
        !				vip_rch(:), kga_rch(:), krea_rch(:), &
        !				kdea_rch(:), ksn_rch(:), ksp_rch(:), &
        !				Isat_rch(:), khnx_rch(:), va_rch(:), &
        !				botalg0(:), kgaF_rch(:), abmax_rch(:), &
        !				kreaF_rch(:), kexaF_rch(:), kdeaF_rch(:), &
        !				ksnF_rch(:), kspF_rch(:), IsatF_rch(:), &
        !				khnxF_rch(:), NINbmin_rch(:), NIPbmin_rch(:), &
        !				NINbupmax_rch(:), NIPbupmax_rch(:), KqN_rch(:), &
        !				KqP_rch(:), NUpWCfrac_rch(:), PUpWCfrac_rch(:), &
        !				kdt_rch(:), vdt_rch(:), kpath_rch(:), &
        !				vpath_rch(:), apath_rch(:), kgaH_rch(:), &
        !				kscH_rch(:), kinhcH_rch(:), kreaH_rch(:), &
        !				kdeaH_rch(:), ksnH_rch(:), kspH_rch(:), &
        !				khnxH_rch(:), ahmax_rch(:), &

        !gp 16-Jul-08
        !REAL(DP), ALLOCATABLE :: lats(:), lond(:), lonm(:), lons(:), Q(:), BB(:), &
        !				SS1(:), SS2(:), s(:), nm(:) , alp1(:), bet1(:), &
        !				alp2(:), bet2(:), Ediff(:), Frsed(:), &
        !				Frsod(:), SODspec(:), JCH4spec(:), JNH4spec(:), &
        !				JSRPspec(:), Hweir(:), Bweir(:), &
        !				sedThermCond(:), sedThermDiff(:), HsedCM(:), &
        !				HypoExchFrac(:), porosity(:), rhoCpSed(:), &
        !				kaaa(:), vss_rch(:), khc_rch(:), &
        !				kdcs_rch(:), kdc_rch(:), khn_rch(:), &
        !				von_rch(:), kn_rch(:) , ki_rch(:) , &
        !				vdi_rch(:), khp_rch(:), vop_rch(:), &
        !				vip_rch(:), kga_rch(:), krea_rch(:), &
        !				kdea_rch(:), ksn_rch(:), ksp_rch(:), &
        !				Isat_rch(:), khnx_rch(:), va_rch(:), &
        !				botalg0(:), kgaF_rch(:), abmax_rch(:), &
        !				krea1F_rch(:), krea2F_rch(:), kexaF_rch(:), kdeaF_rch(:), &
        !				ksnF_rch(:), kspF_rch(:), IsatF_rch(:), &
        !				khnxF_rch(:), NINbmin_rch(:), NIPbmin_rch(:), &
        !				NINbupmax_rch(:), NIPbupmax_rch(:), KqN_rch(:), &
        !				KqP_rch(:), NUpWCfrac_rch(:), PUpWCfrac_rch(:), &
        !				kdt_rch(:), vdt_rch(:), kpath_rch(:), &
        !				vpath_rch(:), apath_rch(:), kgaH_rch(:), &
        !				kscH_rch(:), kinhcH_rch(:), kreaH_rch(:), &
        !				kdeaH_rch(:), ksnH_rch(:), kspH_rch(:), &
        !				khnxH_rch(:), ahmax_rch(:), &
        !
        !				!gp 21-Nov-06
        !				!kgen_rch(:), vgen_rch(:)
        !				kgen_rch(:), vgen_rch(:), &
        !				Te_ini(:), c01_ini(:), c02_ini(:), c03_ini(:), &
        !				c04_ini(:), c05_ini(:), c06_ini(:), &
        !				c07_ini(:), c08_ini(:), c09_ini(:), &
        !				c10_ini(:), c11_ini(:), c12_ini(:), &
        !				c13_ini(:), c14_ini(:), c15_ini(:), &
        !				pH_ini(:), c17_ini(:), NINb_ini(:), NIPb_ini(:)
        REAL(DP), ALLOCATABLE :: lats(:), lond(:), lonm(:), lons(:), Q(:), BB(:), &
            SS1(:), SS2(:), s(:), nm(:) , alp1(:), bet1(:), &
            alp2(:), bet2(:), Ediff(:), Frsed(:), &
            Frsod(:), SODspec(:), JCH4spec(:), JNH4spec(:), &
            JSRPspec(:), Hweir(:), Bweir(:), &
            sedThermCond(:), sedThermDiff(:), HsedCM(:), &
            HypoExchFrac(:), porosity(:), rhoCpSed(:), SKOP(:), &
            kaaa(:), vss_rch(:), khc_rch(:), &
            kdcs_rch(:), kdc_rch(:), khn_rch(:), &
            von_rch(:), kn_rch(:) , ki_rch(:) , &
            vdi_rch(:), khp_rch(:), vop_rch(:), &
            vip_rch(:), kga_rch(:), krea_rch(:), &
            kdea_rch(:), ksn_rch(:), ksp_rch(:), &
            Isat_rch(:), khnx_rch(:), va_rch(:), &
            botalg0(:), kgaF_rch(:), abmax_rch(:), &
            krea1F_rch(:), krea2F_rch(:), kexaF_rch(:), kdeaF_rch(:), &
            ksnF_rch(:), kspF_rch(:), IsatF_rch(:), &
            khnxF_rch(:), NINbmin_rch(:), NIPbmin_rch(:), &
            NINbupmax_rch(:), NIPbupmax_rch(:), KqN_rch(:), &
            KqP_rch(:), NUpWCfrac_rch(:), PUpWCfrac_rch(:), &
            kdt_rch(:), vdt_rch(:), kpath_rch(:), &
            vpath_rch(:), apath_rch(:), kgaH_rch(:), &
            kscH_rch(:), kinhcH_rch(:), kreaH_rch(:), &
            kdeaH_rch(:), ksnH_rch(:), kspH_rch(:), &
            khnxH_rch(:), ahmax_rch(:), &
            kgen_rch(:), vgen_rch(:), &
            Te_ini(:), c01_ini(:), c02_ini(:), c03_ini(:), &
            c04_ini(:), c05_ini(:), c06_ini(:), &
            c07_ini(:), c08_ini(:), c09_ini(:), &
            c10_ini(:), c11_ini(:), c12_ini(:), &
            c13_ini(:), c14_ini(:), c15_ini(:), &
            pH_ini(:), c17_ini(:), NINb_ini(:), NIPb_ini(:)

        CHARACTER(LEN=30), ALLOCATABLE :: rlab1(:), rlab2(:), rname(:)

        !gp 13-Feb-06
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom, nfacBras, atcRyanStolz

        !gp 24-Jun-09
        !gp 16-Jul-08
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac, nfacBras, atcRyanStolz
        !REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac, nfacBras, atcRyanStolz, kbrut
        REAL(DP) PAR, kep, kela, kenla, kess, kepom, kemac, nfacBras, atcRyanStolz, kbrut, KCL1, KCL2

        CHARACTER(LEN=30) solarMethod, longatMethod, fUwMethod

        !point load/source
        CHARACTER(LEN=30), ALLOCATABLE :: PtName(:)
        REAL(DP), ALLOCATABLE :: xptt(:), Qptta(:), Qptt(:), TepttMean(:), TepttAmp(:), TepttMaxTime(:)
        REAL(DP), ALLOCATABLE :: cpttMean(:, :), cpttAmp(:, :), cpttMaxTime(:, :)
        REAL(DP), ALLOCATABLE :: phpttMean(:), phpttAmp(:), phpttMaxTime(:)
        !diffusion load/source
        REAL(DP), ALLOCATABLE :: xdup(:), xddn(:), Qdifa(:), Qdif(:), pHind(:), Tedif(:), cdif(:,:)
        CHARACTER(LEN=30), ALLOCATABLE :: DiffName(:)
        !Rates
        REAL(DP) vss, mgC, mgN, mgP, mgD, mgA
        REAL(DP) tka, roc, ron
        REAL(DP) Ksocf, Ksona, Ksodn, Ksop, Ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von
        REAL(DP) kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        REAL(DP) kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat

        !gp 03-Apr-08
        !REAL(DP) khnx, va, kgaF, tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, abmax
        REAL(DP) khnx, va, kgaF, tkgaF, krea1F, krea2F, tkreaF, kexaF, tkexaF, kdeaF, abmax

        REAL(DP) tkdeaF, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt

        !gp 10-Jan-05 REAL(DP) NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP
        REAL(DP) NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP, KqNratio, KqPratio

        !gp 26-Jan-06
        REAL(DP) NUpWCfrac, PUpWCfrac

        REAL(DP) kpath, tkpath, vpath

        !gp 30-Nov-04
        REAL(DP) apath, kgen, tkgen, vgen		!gp 30-Nov-04 new params for pathogen and generic constituent

        !gp 08-Dec-04
        CHARACTER(LEN=30) useGenericAsCOD

        REAL(DP) pco2
        !org C, N inhition model, Denitrification enhance model, Algae light model, perihyte light model
        CHARACTER(LEN=30) xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, typeF
        CHARACTER(LEN=30) kawindmethod	!Reaeartion wind effect
        !DATA
        INTEGER(I4B) kk, nteda, nhydda, nwqd(3), nrp, nwqdiur
        REAL(DP) junk !Junk data
        CHARACTER(LEN=30) xxx, shtnam(3) !Junk string,

        !GP 23-Nov-09
        !INTEGER(I4B) dlstime
        REAL(DP) dlstime

        LOGICAL(LGT)  downstreamBoundary

        !gp 03-Nov-04 hyporheic biofilm rates
        CHARACTER(LEN=30) typeH, xdum8
        REAL(DP) kgaH, tkgaH, kscH, kinhcH

        !gp 15-Nov-04 HCO3- use by phytoplankton and bottom algae
        CHARACTER(LEN=30) hco3use, hco3useF

        !gp 15-Nov-04 level 2 hyporheic biofilm rates
        REAL(DP) kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax

        TYPE(t_water_quality), DIMENSION(:), ALLOCATABLE :: HwIn	!HEADWATERS
        TYPE(t_water_quality), DIMENSION(:),	POINTER :: DBin		!Downstream bondary

        !METEOROLOGY DATA

        !gp 16-Jul-08
        !REAL(DP), ALLOCATABLE :: shadeHH(:,:), TaHH(:,:), TdHH(:,:), UwHH(:,:), ccHH(:,:)
        REAL(DP), ALLOCATABLE :: shadeHH(:,:), TaHH(:,:), TdHH(:,:), UwHH(:,:), ccHH(:,:), solarHH(:,:)

        !/* Read in input file */
        READ(8,*) BASINNAME, FILENAME, PATH, TITLE
        READ (8,*)  month, day, year
        !gp 28-Oct-04 READ(8,*) timezone, pco2, dtuser, tf, IMeth, IMethpH

        !gp 24-Jun-09
        !READ(8,*) timezone, pco2, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, showDielResults, stateVariables, calcSedFlux, simAlk	!gp 26-Oct-07
        READ(8,*) timezone, pco2, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, showDielResults, &
            stateVariables, calcSedFlux, simAlk, writeDynamic

        !READ(8,*) IMeth, IMethpH

        !gp 24-Jun-09
        !gp 17-Nov-04 system= SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
        !gp 										timezone, dtuser, tf, IMeth, IMethpH)
        !system= SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
        ! 										timezone, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, &
        !										showDielResults, stateVariables, calcSedFlux, simAlk)	!gp 26-Oct-07
        system= SystemParams_(BASINNAME, FILENAME, PATH, TITLE, year, month, day, &
            timezone, dtuser, tf, IMeth, IMethpH, simHyporheicWQ, &
            showDielResults, stateVariables, calcSedFlux, simAlk, writeDynamic)

        !reach data
        !gp 17-Nov-04 READ(8,*) nRch		!total reach number
        READ(8,*) nRch, geoMethod			!gp 17,Nov-04 number of reaches and Depth or Width for col T and U

        ALLOCATE (xrdn(0:nRch+1), STAT=status(1));	ALLOCATE (elev1(0:nRch+1), STAT=status(2))
        ALLOCATE (elev2(0:nRch+1), STAT=status(3));	ALLOCATE (latd(0:nRch+1), STAT=status(4))
        ALLOCATE (latm(0:nRch+1), STAT=status(5));	ALLOCATE (lats(0:nRch+1), STAT=status(6))
        ALLOCATE (lond(0:nRch+1), STAT=status(7));	ALLOCATE (lonm(0:nRch+1), STAT=status(8))
        ALLOCATE (lons(0:nRch+1), STAT=status(9));	ALLOCATE (Q(0:nRch+1), STAT=status(10))
        ALLOCATE (BB(0:nRch+1), STAT=status(11));		ALLOCATE (SS1(0:nRch+1), STAT=status(12))
        ALLOCATE (SS2(0:nRch+1), STAT=status(13));	ALLOCATE (s(0:nRch+1), STAT=status(14))
        ALLOCATE (nm(0:nRch+1), STAT=status(15)); 	ALLOCATE (alp1(0:nRch+1), STAT=status(16))
        ALLOCATE (bet1(0:nRch+1), STAT=status(17));	ALLOCATE (alp2(0:nRch+1), STAT=status(18))
        ALLOCATE (bet2(0:nRch+1), STAT=status(19));	ALLOCATE (Ediff(0:nRch+1), STAT=status(20))
        ALLOCATE (kaaa(0:nRch+1), STAT=status(21));	ALLOCATE (Frsed(0:nRch+1), STAT=status(22))
        ALLOCATE (Frsod(0:nRch+1), STAT=status(23));ALLOCATE (SODspec(0:nRch+1), STAT=status(24))
        ALLOCATE (JCH4spec(0:nRch+1), STAT=status(25));	ALLOCATE (JNH4spec(0:nRch+1), STAT=status(26))
        ALLOCATE (JSRPspec(0:nRch+1), STAT=status(27));	ALLOCATE (Hweir(0:nRch+1), STAT=status(28))
        ALLOCATE (rlab1(0:nRch+1), STAT=status(29));	ALLOCATE (rlab2(0:nRch+1), STAT=status(30))
        ALLOCATE (rname(0:nRch+1), STAT=status(31));ALLOCATE (Bweir(0:nRch+1), STAT=status(32))
        ALLOCATE (sedThermCond(0:nRch+1), STAT=status(33));ALLOCATE (sedThermDiff(0:nRch+1), STAT=status(34))	!gp 20-Oct-04
        ALLOCATE (HsedCM(0:nRch+1), STAT=status(35));ALLOCATE (HypoExchFrac(0:nRch+1), STAT=status(36))			!gp 20-Oct-04
        ALLOCATE (porosity(0:nRch+1), STAT=status(37));ALLOCATE (rhoCpSed(0:nRch+1), STAT=status(38))			!gp 20-Oct-04

        !gp 16-Jul-08
        ALLOCATE (SKOP(0:nRch+1), STAT=status(114))															!gp 11-Jan-05

        ALLOCATE (botalg0(0:nRch+1), STAT=status(39))															!gp 11-Jan-05

        !gp 07-Feb-07
        !ALLOCATE (kaaa(0:nRch+1), STAT=status(39))
        ALLOCATE (vss_rch(0:nRch+1), STAT=status(40))				!inorganic suspended solids settling vol
        ALLOCATE (khc_rch(0:nRch+1), STAT=status(41))					!slow CBOD hydrolysis rate
        ALLOCATE (kdcs_rch(0:nRch+1), STAT=status(42))					!slow CBOD oxidation rate
        ALLOCATE (kdc_rch(0:nRch+1), STAT=status(43))					!fast CBOD oxidation rate
        ALLOCATE (khn_rch(0:nRch+1), STAT=status(44))					!organic N hydrolysis rate
        ALLOCATE (von_rch(0:nRch+1), STAT=status(45))					!Organic N settling velocity
        ALLOCATE (kn_rch(0:nRch+1), STAT=status(46))					!Ammonium nitrification rate
        ALLOCATE (ki_rch(0:nRch+1), STAT=status(47))					!Nitrate denitrification
        ALLOCATE (vdi_rch(0:nRch+1), STAT=status(48))					!Nitrate sed denitrification transfer coeff
        ALLOCATE (khp_rch(0:nRch+1), STAT=status(49))					!organic P hydrolysis
        ALLOCATE (vop_rch(0:nRch+1), STAT=status(50))					!organic P settling velocity
        ALLOCATE (vip_rch(0:nRch+1), STAT=status(51))					!inorganic P settling velocity
        ALLOCATE (kga_rch(0:nRch+1), STAT=status(52))					!Phytoplankton MAX growth rate
        ALLOCATE (krea_rch(0:nRch+1), STAT=status(53))					!Phytoplankton respiration rate
        ALLOCATE (kdea_rch(0:nRch+1), STAT=status(54))					!Phytoplankton death rate
        ALLOCATE (ksn_rch(0:nRch+1), STAT=status(55))					!Phytoplankton N half-sat
        ALLOCATE (ksp_rch(0:nRch+1), STAT=status(56))					!Phytoplankton P half-sat
        ALLOCATE (Isat_rch(0:nRch+1), STAT=status(57))					!Phytoplankton light sat
        ALLOCATE (khnx_rch(0:nRch+1), STAT=status(58))					!Phytoplankton ammonia preference
        ALLOCATE (va_rch(0:nRch+1), STAT=status(59))					!Phytoplankton settling velocity
        !ALLOCATE (botalg0(0:nRch+1), STAT=status(39))				!Bottom plant initial bionass
        ALLOCATE (kgaF_rch(0:nRch+1), STAT=status(60))					!Bottom plant MAX growth rate
        ALLOCATE (abmax_rch(0:nRch+1), STAT=status(61))					!Bottom plant first-order carrying capacity

        !gp 03-Apr-08
        !ALLOCATE (kreaF_rch(0:nRch+1), STAT=status(62))					!Bottom plant respiration rate
        ALLOCATE (krea1F_rch(0:nRch+1), STAT=status(62))					!Bottom plant basal respiration rate
        ALLOCATE (krea2F_rch(0:nRch+1), STAT=status(113))					!Bottom plant photo respiration rate

        ALLOCATE (kexaF_rch(0:nRch+1), STAT=status(63))					!Bottom plant excretion rate
        ALLOCATE (kdeaF_rch(0:nRch+1), STAT=status(64))					!Bottom plant death rate
        ALLOCATE (ksnF_rch(0:nRch+1), STAT=status(65))					!Bottom plant external N half-sat
        ALLOCATE (kspF_rch(0:nRch+1), STAT=status(66))					!Bottom plant external P half-sat
        ALLOCATE (IsatF_rch(0:nRch+1), STAT=status(67))					!Bottom plant light sat
        ALLOCATE (khnxF_rch(0:nRch+1), STAT=status(68))					!Bottom plant ammonia preference
        ALLOCATE (NINbmin_rch(0:nRch+1), STAT=status(69))				!Bottom plant subistence quota for N
        ALLOCATE (NIPbmin_rch(0:nRch+1), STAT=status(70))				!Bottom plant subistence quota for P
        ALLOCATE (NINbupmax_rch(0:nRch+1), STAT=status(71))				!Bottom plant max uptake rate for N
        ALLOCATE (NIPbupmax_rch(0:nRch+1), STAT=status(72))				!Bottom plant max uptake rate for P
        ALLOCATE (KqN_rch(0:nRch+1), STAT=status(73))					!Bottom plant internal N half-sat
        ALLOCATE (KqP_rch(0:nRch+1), STAT=status(74))					!Bottom plant internal P half-sat
        ALLOCATE (NUpWCfrac_rch(0:nRch+1), STAT=status(75))				!Bottom plant N uptake fraction from water column
        ALLOCATE (PUpWCfrac_rch(0:nRch+1), STAT=status(76))				!Bottom plant P uptake fraction from water column
        ALLOCATE (kdt_rch(0:nRch+1), STAT=status(77))					!POM dissolution rate
        ALLOCATE (vdt_rch(0:nRch+1), STAT=status(78))					!POM settling velocity
        ALLOCATE (kpath_rch(0:nRch+1), STAT=status(79))					!pathogen dieoff rate
        ALLOCATE (vpath_rch(0:nRch+1), STAT=status(80))					!pathogen settling velocity
        ALLOCATE (apath_rch(0:nRch+1), STAT=status(81))					!pathogen light alpha
        ALLOCATE (kgaH_rch(0:nRch+1), STAT=status(82))					!Hyporheic heterotrophs MAX growth rate
        ALLOCATE (kscH_rch(0:nRch+1), STAT=status(83))					!Hyporheic heterotrophs CBOD half-sat
        ALLOCATE (kinhcH_rch(0:nRch+1), STAT=status(84))				!Hyporheic heterotrophs O2 inhibition
        ALLOCATE (kreaH_rch(0:nRch+1), STAT=status(85))					!Hyporheic heterotrophs respiration rate
        ALLOCATE (kdeaH_rch(0:nRch+1), STAT=status(86))					!Hyporheic heterotrophs death rate
        ALLOCATE (ksnH_rch(0:nRch+1), STAT=status(87))					!Hyporheic heterotrophs N half-sat
        ALLOCATE (kspH_rch(0:nRch+1), STAT=status(88))					!Hyporheic heterotrophs P half-sat
        ALLOCATE (khnxH_rch(0:nRch+1), STAT=status(89))					!Hyporheic heterotrophs ammonia preference
        ALLOCATE (ahmax_rch(0:nRch+1), STAT=status(90))					!Hyporheic heterotrophs first-order carrying capacity
        ALLOCATE (kgen_rch(0:nRch+1), STAT=status(91))					!generic constituent dissolution rate
        ALLOCATE (vgen_rch(0:nRch+1), STAT=status(92))					!generic constituent settling velocity

        !gp 21-Nov-06 initial conditions of state variables
        ALLOCATE (Te_ini(0:nRch+1), STAT=status(93))
        ALLOCATE (c01_ini(0:nRch+1), STAT=status(94))
        ALLOCATE (c02_ini(0:nRch+1), STAT=status(95))
        ALLOCATE (c03_ini(0:nRch+1), STAT=status(96))
        ALLOCATE (c04_ini(0:nRch+1), STAT=status(97))
        ALLOCATE (c05_ini(0:nRch+1), STAT=status(98))
        ALLOCATE (c06_ini(0:nRch+1), STAT=status(99))
        ALLOCATE (c07_ini(0:nRch+1), STAT=status(100))
        ALLOCATE (c08_ini(0:nRch+1), STAT=status(101))
        ALLOCATE (c09_ini(0:nRch+1), STAT=status(102))
        ALLOCATE (c10_ini(0:nRch+1), STAT=status(103))
        ALLOCATE (c11_ini(0:nRch+1), STAT=status(104))
        ALLOCATE (c12_ini(0:nRch+1), STAT=status(105))
        ALLOCATE (c13_ini(0:nRch+1), STAT=status(106))
        ALLOCATE (c14_ini(0:nRch+1), STAT=status(107))
        ALLOCATE (c15_ini(0:nRch+1), STAT=status(108))
        ALLOCATE (pH_ini(0:nRch+1), STAT=status(109))
        ALLOCATE (c17_ini(0:nRch+1), STAT=status(110))
        ALLOCATE (NINb_ini(0:nRch+1), STAT=status(111))
        ALLOCATE (NIPb_ini(0:nRch+1), STAT=status(112))
        !gp 03-Apr-08 status(113) is used above for krea2F
        !gp 16-Jul-08 status(114) is used above for SKOP

        !classriver module

        !gp 07-Feb-06
        !DO i = 0, nRch + 1
        !	READ(8,*) rlab1(i), rlab2(i), rname(i), xrdn(i), elev1(i), elev2(i), latd(i), latm(i), &
        !        lats(i), lond(i), lonm(i), lons(i), Q(i), BB(i), SS1(i), SS2(i), s(i), nm(i), &
        !        alp1(i), bet1(i), alp2(i), bet2(i), Ediff(i), kaaa(i), Frsed(i), Frsod(i), &
        !        SODspec(i), JCH4spec(i), JNH4spec(i), JSRPspec(i), Hweir(i), Bweir(i), &
        !        sedThermCond(i), sedThermDiff(i), HsedCM(i), &
        !		HypoExchFrac(i), porosity(i), rhoCpSed(i), botalg0(i)			!gp 11-Jan-05
        !
        !END DO

        !gp 16-Jul-06
        !DO i = 0, nRch + 1
        !	READ(8,*) rlab1(i), rlab2(i), rname(i), xrdn(i), elev1(i), elev2(i), latd(i), latm(i), &
        !        lats(i), lond(i), lonm(i), lons(i), Q(i), BB(i), SS1(i), SS2(i), s(i), nm(i), &
        !        alp1(i), bet1(i), alp2(i), bet2(i), Ediff(i), Frsed(i), Frsod(i), &
        !        SODspec(i), JCH4spec(i), JNH4spec(i), JSRPspec(i), Hweir(i), Bweir(i), &
        !        sedThermCond(i), sedThermDiff(i), HsedCM(i), &
        !		HypoExchFrac(i), porosity(i), rhoCpSed(i)
        !END DO
        DO i = 0, nRch + 1
            READ(8,*) rlab1(i), rlab2(i), rname(i), xrdn(i), elev1(i), elev2(i), latd(i), latm(i), &
                lats(i), lond(i), lonm(i), lons(i), Q(i), BB(i), SS1(i), SS2(i), s(i), nm(i), &
                alp1(i), bet1(i), alp2(i), bet2(i), Ediff(i), Frsed(i), Frsod(i), &
                SODspec(i), JCH4spec(i), JNH4spec(i), JSRPspec(i), Hweir(i), Bweir(i), &
                sedThermCond(i), sedThermDiff(i), HsedCM(i), &
                HypoExchFrac(i), porosity(i), rhoCpSed(i), SKOP(i)
        END DO

        !gp 21-Nov-06
        !WRITE(10,*) 'done thru classreadfile ALLOCATE'
        !DO i=1, nRch
        !	READ(8,*) kaaa(i), vss_rch(i), khc_rch(i), &
        !				kdcs_rch(i), kdc_rch(i), khn_rch(i), &
        !				von_rch(i), kn_rch(i) , ki_rch(i) , &
        !				vdi_rch(i), khp_rch(i), vop_rch(i), &
        !				vip_rch(i), kga_rch(i), krea_rch(i), &
        !				kdea_rch(i), ksn_rch(i), ksp_rch(i), &
        !				Isat_rch(i), khnx_rch(i), va_rch(i), &
        !				botalg0(i), kgaF_rch(i), abmax_rch(i), &
        !				kreaF_rch(i), kexaF_rch(i), kdeaF_rch(i), &
        !				ksnF_rch(i), kspF_rch(i), IsatF_rch(i), &
        !				khnxF_rch(i), NINbmin_rch(i), NIPbmin_rch(i), &
        !				NINbupmax_rch(i), NIPbupmax_rch(i), KqN_rch(i), &
        !				KqP_rch(i), NUpWCfrac_rch(i), PUpWCfrac_rch(i), &
        !				kdt_rch(i), vdt_rch(i), kpath_rch(i), &
        !				vpath_rch(i), apath_rch(i), kgaH_rch(i), &
        !				kscH_rch(i), kinhcH_rch(i), kreaH_rch(i), &
        !				kdeaH_rch(i), ksnH_rch(i), kspH_rch(i), &
        !				khnxH_rch(i), ahmax_rch(i), &
        !				kgen_rch(i), vgen_rch(i)
        !END DO
        DO i=1, nRch

            !gp 03-Apr-08
            !READ(8,*) kaaa(i), vss_rch(i), khc_rch(i), &
            !			kdcs_rch(i), kdc_rch(i), khn_rch(i), &
            !			von_rch(i), kn_rch(i) , ki_rch(i) , &
            !			vdi_rch(i), khp_rch(i), vop_rch(i), &
            !			vip_rch(i), kga_rch(i), krea_rch(i), &
            !			kdea_rch(i), ksn_rch(i), ksp_rch(i), &
            !			Isat_rch(i), khnx_rch(i), va_rch(i), &
            !			kgaF_rch(i), abmax_rch(i), &
            !			kreaF_rch(i), kexaF_rch(i), kdeaF_rch(i), &
            !			ksnF_rch(i), kspF_rch(i), IsatF_rch(i), &
            !			khnxF_rch(i), NINbmin_rch(i), NIPbmin_rch(i), &
            !			NINbupmax_rch(i), NIPbupmax_rch(i), KqN_rch(i), &
            !			KqP_rch(i), NUpWCfrac_rch(i), PUpWCfrac_rch(i), &
            !			kdt_rch(i), vdt_rch(i), kpath_rch(i), &
            !			vpath_rch(i), apath_rch(i), kgaH_rch(i), &
            !			kscH_rch(i), kinhcH_rch(i), kreaH_rch(i), &
            !			kdeaH_rch(i), ksnH_rch(i), kspH_rch(i), &
            !			khnxH_rch(i), ahmax_rch(i), &
            !			kgen_rch(i), vgen_rch(i)
            READ(8,*) kaaa(i), vss_rch(i), khc_rch(i), &
                kdcs_rch(i), kdc_rch(i), khn_rch(i), &
                von_rch(i), kn_rch(i) , ki_rch(i) , &
                vdi_rch(i), khp_rch(i), vop_rch(i), &
                vip_rch(i), kga_rch(i), krea_rch(i), &
                kdea_rch(i), ksn_rch(i), ksp_rch(i), &
                Isat_rch(i), khnx_rch(i), va_rch(i), &
                kgaF_rch(i), abmax_rch(i), &
                krea1F_rch(i), krea2F_rch(i), kexaF_rch(i), kdeaF_rch(i), &
                ksnF_rch(i), kspF_rch(i), IsatF_rch(i), &
                khnxF_rch(i), NINbmin_rch(i), NIPbmin_rch(i), &
                NINbupmax_rch(i), NIPbupmax_rch(i), KqN_rch(i), &
                KqP_rch(i), NUpWCfrac_rch(i), PUpWCfrac_rch(i), &
                kdt_rch(i), vdt_rch(i), kpath_rch(i), &
                vpath_rch(i), apath_rch(i), kgaH_rch(i), &
                kscH_rch(i), kinhcH_rch(i), kreaH_rch(i), &
                kdeaH_rch(i), ksnH_rch(i), kspH_rch(i), &
                khnxH_rch(i), ahmax_rch(i), &
                kgen_rch(i), vgen_rch(i)

            !	WRITE(10,'(I5, 4F10.3)') i, kaaa(i), vss_rch(i), kgen_rch(i), vgen_rch(i)
        END DO
        DO i=1, nRch
            READ(8,*) Te_ini(i), pH_ini(i), c17_ini(i), NINb_ini(i), NIPb_ini(i)
            !	WRITE(10,'(I5, 4F10.3)') i, Te_ini(i), c17_ini(i), NINb_ini(i), NIPb_ini(i)
        END DO
        DO i=1, nRch
            READ(8,*) c01_ini(i), c02_ini(i), &
                c03_ini(i), c04_ini(i), c05_ini(i), &
                c06_ini(i), c07_ini(i), c08_ini(i), &
                c09_ini(i), c10_ini(i), c11_ini(i), &
                c12_ini(i), c13_ini(i), c14_ini(i), &
                c15_ini(i)
            !	WRITE(10,'(I5, 4F10.3)') i, c01_ini(i), c02_ini(i), c14_ini(i), c15_ini(i)
        END DO
        !WRITE(10,*) 'done thru classreadfile reading new init cond'

        !gp 20-Oct-04 note units for new inputs:
        !	sedThermCond = (cal/s) per (cm degC)
        !	sedThermDiff = cm^2/sec
        !	HsedCM = cm
        !	HypoExchFrac = unitless fraction of streamflow
        !	porosity = unitless fraction of volume
        !	rhoCpSed = cal/(cm^3 deg C)

        !gp 17-Nov-04 Topo= RiverTopo_(nRch, rlab2, rname, xrdn)
        Topo= rivertopo_type(nRch, rlab2, rname, xrdn, geoMethod)		!gp 17-Nov-04

        !gp 07-Feb-06
        !hydrau= hydraulics_(nRch, xrdn, elev1, elev2, latd, latm, lats, &
        !								 lond, lonm, lons, Q, BB, SS1, SS2, s, nm, alp1, bet1, alp2, bet2, &
        !								 Ediff, kaaa, Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, Hweir, &
        !								 Bweir, &
        !								 sedThermCond, sedThermDiff, HsedCM, &
        !								 HypoExchFrac, porosity, rhoCpSed, botalg0)	!gp 11-Jan-05

        !gp 21-Nov-06
        !hydrau= hydraulics_(nRch, xrdn, elev1, elev2, latd, latm, lats, &
        !								lond, lonm, lons, Q, BB, SS1, SS2, s, nm, alp1, bet1, alp2, bet2, &
        !								Ediff, Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, &
        !								Hweir, Bweir, &
        !								sedThermCond, sedThermDiff, HsedCM, &
        !								HypoExchFrac, porosity, rhoCpSed, &
        !								kaaa, vss_rch, khc_rch, &
        !								kdcs_rch, kdc_rch, khn_rch, &
        !								von_rch, kn_rch, ki_rch, &
        !								vdi_rch, khp_rch, vop_rch, &
        !								vip_rch, kga_rch, krea_rch, &
        !								kdea_rch, ksn_rch, ksp_rch, &
        !								Isat_rch, khnx_rch, va_rch, &
        !								botalg0, kgaF_rch, abmax_rch, &
        !								kreaF_rch, kexaF_rch, kdeaF_rch, &
        !								ksnF_rch, kspF_rch, IsatF_rch, &
        !								khnxF_rch, NINbmin_rch, NIPbmin_rch, &
        !								NINbupmax_rch, NIPbupmax_rch, KqN_rch, &
        !								KqP_rch, NUpWCfrac_rch, PUpWCfrac_rch, &
        !								kdt_rch, vdt_rch, kpath_rch, &
        !								vpath_rch, apath_rch, kgaH_rch, &
        !								kscH_rch, kinhcH_rch, kreaH_rch, &
        !								kdeaH_rch, ksnH_rch, kspH_rch, &
        !								khnxH_rch, ahmax_rch, &
        !								kgen_rch, vgen_rch)

        !gp 03-Apr-08
        !hydrau= hydraulics_(nRch, xrdn, elev1, elev2, latd, latm, lats, &
        !								lond, lonm, lons, Q, BB, SS1, SS2, s, nm, alp1, bet1, alp2, bet2, &
        !								Ediff, Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, &
        !								Hweir, Bweir, &
        !								sedThermCond, sedThermDiff, HsedCM, &
        !								HypoExchFrac, porosity, rhoCpSed, &
        !								kaaa, vss_rch, khc_rch, &
        !								kdcs_rch, kdc_rch, khn_rch, &
        !								von_rch, kn_rch, ki_rch, &
        !								vdi_rch, khp_rch, vop_rch, &
        !								vip_rch, kga_rch, krea_rch, &
        !								kdea_rch, ksn_rch, ksp_rch, &
        !								Isat_rch, khnx_rch, va_rch, &
        !								kgaF_rch, abmax_rch, &
        !								kreaF_rch, kexaF_rch, kdeaF_rch, &
        !								ksnF_rch, kspF_rch, IsatF_rch, &
        !								khnxF_rch, NINbmin_rch, NIPbmin_rch, &
        !								NINbupmax_rch, NIPbupmax_rch, KqN_rch, &
        !								KqP_rch, NUpWCfrac_rch, PUpWCfrac_rch, &
        !								kdt_rch, vdt_rch, kpath_rch, &
        !								vpath_rch, apath_rch, kgaH_rch, &
        !								kscH_rch, kinhcH_rch, kreaH_rch, &
        !								kdeaH_rch, ksnH_rch, kspH_rch, &
        !								khnxH_rch, ahmax_rch, &
        !								kgen_rch, vgen_rch, &
        !								Te_ini, c01_ini, c02_ini, c03_ini, &
        !								c04_ini, c05_ini, c06_ini, &
        !								c07_ini, c08_ini, c09_ini, &
        !								c10_ini, c11_ini, c12_ini, &
        !								c13_ini, c14_ini, c15_ini, &
        !								pH_ini, c17_ini, NINb_ini, NIPb_ini)

        !gp 16-Jul-08
        !hydrau= hydraulics_(nRch, xrdn, elev1, elev2, latd, latm, lats, &
        !								lond, lonm, lons, Q, BB, SS1, SS2, s, nm, alp1, bet1, alp2, bet2, &
        !								Ediff, Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, &
        !								Hweir, Bweir, &
        !								sedThermCond, sedThermDiff, HsedCM, &
        !								HypoExchFrac, porosity, rhoCpSed, &
        !								kaaa, vss_rch, khc_rch, &
        !								kdcs_rch, kdc_rch, khn_rch, &
        !								von_rch, kn_rch, ki_rch, &
        !								vdi_rch, khp_rch, vop_rch, &
        !								vip_rch, kga_rch, krea_rch, &
        !								kdea_rch, ksn_rch, ksp_rch, &
        !								Isat_rch, khnx_rch, va_rch, &
        !								kgaF_rch, abmax_rch, &
        !								krea1F_rch, krea2F_rch, kexaF_rch, kdeaF_rch, &
        !								ksnF_rch, kspF_rch, IsatF_rch, &
        !								khnxF_rch, NINbmin_rch, NIPbmin_rch, &
        !								NINbupmax_rch, NIPbupmax_rch, KqN_rch, &
        !								KqP_rch, NUpWCfrac_rch, PUpWCfrac_rch, &
        !								kdt_rch, vdt_rch, kpath_rch, &
        !								vpath_rch, apath_rch, kgaH_rch, &
        !								kscH_rch, kinhcH_rch, kreaH_rch, &
        !								kdeaH_rch, ksnH_rch, kspH_rch, &
        !								khnxH_rch, ahmax_rch, &
        !								kgen_rch, vgen_rch, &
        !								Te_ini, c01_ini, c02_ini, c03_ini, &
        !								c04_ini, c05_ini, c06_ini, &
        !								c07_ini, c08_ini, c09_ini, &
        !								c10_ini, c11_ini, c12_ini, &
        !								c13_ini, c14_ini, c15_ini, &
        !								pH_ini, c17_ini, NINb_ini, NIPb_ini)
        hydrau= hydraulics_(nRch, xrdn, elev1, elev2, latd, latm, lats, &
            lond, lonm, lons, Q, BB, SS1, SS2, s, nm, alp1, bet1, alp2, bet2, &
            Ediff, Frsed, Frsod, SODspec, JCH4spec, JNH4spec, JSRPspec, &
            Hweir, Bweir, &
            sedThermCond, sedThermDiff, HsedCM, &
            HypoExchFrac, porosity, rhoCpSed, SKOP, &
            kaaa, vss_rch, khc_rch, &
            kdcs_rch, kdc_rch, khn_rch, &
            von_rch, kn_rch, ki_rch, &
            vdi_rch, khp_rch, vop_rch, &
            vip_rch, kga_rch, krea_rch, &
            kdea_rch, ksn_rch, ksp_rch, &
            Isat_rch, khnx_rch, va_rch, &
            kgaF_rch, abmax_rch, &
            krea1F_rch, krea2F_rch, kexaF_rch, kdeaF_rch, &
            ksnF_rch, kspF_rch, IsatF_rch, &
            khnxF_rch, NINbmin_rch, NIPbmin_rch, &
            NINbupmax_rch, NIPbupmax_rch, KqN_rch, &
            KqP_rch, NUpWCfrac_rch, PUpWCfrac_rch, &
            kdt_rch, vdt_rch, kpath_rch, &
            vpath_rch, apath_rch, kgaH_rch, &
            kscH_rch, kinhcH_rch, kreaH_rch, &
            kdeaH_rch, ksnH_rch, kspH_rch, &
            khnxH_rch, ahmax_rch, &
            kgen_rch, vgen_rch, &
            Te_ini, c01_ini, c02_ini, c03_ini, &
            c04_ini, c05_ini, c06_ini, &
            c07_ini, c08_ini, c09_ini, &
            c10_ini, c11_ini, c12_ini, &
            c13_ini, c14_ini, c15_ini, &
            pH_ini, c17_ini, NINb_ini, NIPb_ini)

        !gp 21-Nov-06
        !WRITE(10,*) 'done thru classreadfile hydrau'

        DEALLOCATE (botalg0, STAT=status(39))		!gp 11-Jan-05
        DEALLOCATE (rhoCpSed, STAT=status(38))		!gp 20-Oct-04

        !gp 16-Jul-08
        DEALLOCATE (SKOP, STAT=status(114))		!gp 20-Oct-04

        DEALLOCATE (porosity, STAT=status(37))		!gp 20-Oct-04
        DEALLOCATE (HypoExchFrac, STAT=status(36))	!gp 20-Oct-04
        DEALLOCATE (HsedCM, STAT=status(35))		!gp 20-Oct-04
        DEALLOCATE (sedThermDiff, STAT=status(34))	!gp 20-Oct-04
        DEALLOCATE (sedThermCond, STAT=status(33))	!gp 20-Oct-04
        DEALLOCATE (Bweir, STAT=status(32))			!gp 20-Oct-04 (missing from original code)
        DEALLOCATE (rname, STAT=status(31))
        DEALLOCATE (rlab2, STAT=status(30))
        DEALLOCATE (rlab1, STAT=status(29))
        DEALLOCATE (Hweir, STAT=status(28))
        DEALLOCATE (JSRPspec, STAT=status(27))
        DEALLOCATE (JNH4spec, STAT=status(26))
        DEALLOCATE (JCH4spec, STAT=status(25))
        DEALLOCATE (SODspec, STAT=status(24))
        DEALLOCATE (Frsod, STAT=status(23))
        DEALLOCATE (Frsed, STAT=status(22))
        DEALLOCATE (kaaa, STAT=status(21))
        DEALLOCATE (Ediff, STAT=status(20))
        DEALLOCATE (bet2, STAT=status(19))
        DEALLOCATE (alp2, STAT=status(18))
        DEALLOCATE (bet1, STAT=status(17))
        DEALLOCATE (alp1, STAT=status(16))
        DEALLOCATE (nm, STAT=status(15))
        DEALLOCATE (s, STAT=status(14))
        DEALLOCATE (SS2, STAT=status(13))
        DEALLOCATE (SS1, STAT=status(12))
        DEALLOCATE (BB, STAT=status(11))
        DEALLOCATE (Q, STAT=status(10))
        DEALLOCATE (lons, STAT=status(9))
        DEALLOCATE (lonm, STAT=status(8))
        DEALLOCATE (lond, STAT=status(7))
        DEALLOCATE (lats, STAT=status(6))
        DEALLOCATE (latm, STAT=status(5))
        DEALLOCATE (latd, STAT=status(4))
        DEALLOCATE (elev2, STAT=status(3))
        DEALLOCATE (elev1, STAT=status(2))
        DEALLOCATE (xrdn, STAT=status(1))

        !gp 07-Feb-07
        !DEALLOCATE (kaaa, STAT=status(39))
        DEALLOCATE (vss_rch, STAT=status(40))				!inorganic suspended solids settling vol
        DEALLOCATE (khc_rch, STAT=status(41))					!slow CBOD hydrolysis rate
        DEALLOCATE (kdcs_rch, STAT=status(42))					!slow CBOD oxidation rate
        DEALLOCATE (kdc_rch, STAT=status(43))					!fast CBOD oxidation rate
        DEALLOCATE (khn_rch, STAT=status(44))					!organic N hydrolysis rate
        DEALLOCATE (von_rch, STAT=status(45))					!Organic N settling velocity
        DEALLOCATE (kn_rch, STAT=status(46))					!Ammonium nitrification rate
        DEALLOCATE (ki_rch, STAT=status(47))					!Nitrate denitrification
        DEALLOCATE (vdi_rch, STAT=status(48))					!Nitrate sed denitrification transfer coeff
        DEALLOCATE (khp_rch, STAT=status(49))					!organic P hydrolysis
        DEALLOCATE (vop_rch, STAT=status(50))					!organic P settling velocity
        DEALLOCATE (vip_rch, STAT=status(51))					!inorganic P settling velocity
        DEALLOCATE (kga_rch, STAT=status(52))					!Phytoplankton MAX growth rate
        DEALLOCATE (krea_rch, STAT=status(53))					!Phytoplankton respiration rate
        DEALLOCATE (kdea_rch, STAT=status(54))					!Phytoplankton death rate
        DEALLOCATE (ksn_rch, STAT=status(55))					!Phytoplankton N half-sat
        DEALLOCATE (ksp_rch, STAT=status(56))					!Phytoplankton P half-sat
        DEALLOCATE (Isat_rch, STAT=status(57))					!Phytoplankton light sat
        DEALLOCATE (khnx_rch, STAT=status(58))					!Phytoplankton ammonia preference
        DEALLOCATE (va_rch, STAT=status(59))					!Phytoplankton settling velocity
        !DEALLOCATE (botalg0, STAT=status(39))				!Bottom plant initial bionass
        DEALLOCATE (kgaF_rch, STAT=status(60))					!Bottom plant MAX growth rate
        DEALLOCATE (abmax_rch, STAT=status(61))					!Bottom plant first-order carrying capacity

        !gp 03-Apr-08
        !DEALLOCATE (kreaF_rch, STAT=status(62))					!Bottom plant respiration rate
        DEALLOCATE (krea1F_rch, STAT=status(62))					!Bottom plant basal respiration rate
        DEALLOCATE (krea2F_rch, STAT=status(113))					!Bottom plant photo respiration rate

        DEALLOCATE (kexaF_rch, STAT=status(63))					!Bottom plant excretion rate
        DEALLOCATE (kdeaF_rch, STAT=status(64))					!Bottom plant death rate
        DEALLOCATE (ksnF_rch, STAT=status(65))					!Bottom plant external N half-sat
        DEALLOCATE (kspF_rch, STAT=status(66))					!Bottom plant external P half-sat
        DEALLOCATE (IsatF_rch, STAT=status(67))					!Bottom plant light sat
        DEALLOCATE (khnxF_rch, STAT=status(68))					!Bottom plant ammonia preference
        DEALLOCATE (NINbmin_rch, STAT=status(69))				!Bottom plant subistence quota for N
        DEALLOCATE (NIPbmin_rch, STAT=status(70))				!Bottom plant subistence quota for P
        DEALLOCATE (NINbupmax_rch, STAT=status(71))				!Bottom plant max uptake rate for N
        DEALLOCATE (NIPbupmax_rch, STAT=status(72))				!Bottom plant max uptake rate for P
        DEALLOCATE (KqN_rch, STAT=status(73))					!Bottom plant internal N half-sat
        DEALLOCATE (KqP_rch, STAT=status(74))					!Bottom plant internal P half-sat
        DEALLOCATE (NUpWCfrac_rch, STAT=status(75))				!Bottom plant N uptake fraction from water column
        DEALLOCATE (PUpWCfrac_rch, STAT=status(76))				!Bottom plant P uptake fraction from water column
        DEALLOCATE (kdt_rch, STAT=status(77))					!POM dissolution rate
        DEALLOCATE (vdt_rch, STAT=status(78))					!POM settling velocity
        DEALLOCATE (kpath_rch, STAT=status(79))					!pathogen dieoff rate
        DEALLOCATE (vpath_rch, STAT=status(80))					!pathogen settling velocity
        DEALLOCATE (apath_rch, STAT=status(81))					!pathogen light alpha
        DEALLOCATE (kgaH_rch, STAT=status(82))					!Hyporheic heterotrophs MAX growth rate
        DEALLOCATE (kscH_rch, STAT=status(83))					!Hyporheic heterotrophs CBOD half-sat
        DEALLOCATE (kinhcH_rch, STAT=status(84))				!Hyporheic heterotrophs O2 inhibition
        DEALLOCATE (kreaH_rch, STAT=status(85))					!Hyporheic heterotrophs respiration rate
        DEALLOCATE (kdeaH_rch, STAT=status(86))					!Hyporheic heterotrophs death rate
        DEALLOCATE (ksnH_rch, STAT=status(87))					!Hyporheic heterotrophs N half-sat
        DEALLOCATE (kspH_rch, STAT=status(88))					!Hyporheic heterotrophs P half-sat
        DEALLOCATE (khnxH_rch, STAT=status(89))					!Hyporheic heterotrophs ammonia preference
        DEALLOCATE (ahmax_rch, STAT=status(90))					!Hyporheic heterotrophs first-order carrying capacity
        DEALLOCATE (kgen_rch, STAT=status(91))					!generic constituent dissolution rate
        DEALLOCATE (vgen_rch, STAT=status(92))					!generic constituent settling velocity

        !gp 21-Nov-06 initial conditions of state variables
        DEALLOCATE (Te_ini, STAT=status(93))
        DEALLOCATE (c01_ini, STAT=status(94))
        DEALLOCATE (c02_ini, STAT=status(95))
        DEALLOCATE (c03_ini, STAT=status(96))
        DEALLOCATE (c04_ini, STAT=status(97))
        DEALLOCATE (c05_ini, STAT=status(98))
        DEALLOCATE (c06_ini, STAT=status(99))
        DEALLOCATE (c07_ini, STAT=status(100))
        DEALLOCATE (c08_ini, STAT=status(101))
        DEALLOCATE (c09_ini, STAT=status(102))
        DEALLOCATE (c10_ini, STAT=status(103))
        DEALLOCATE (c11_ini, STAT=status(104))
        DEALLOCATE (c12_ini, STAT=status(105))
        DEALLOCATE (c13_ini, STAT=status(106))
        DEALLOCATE (c14_ini, STAT=status(107))
        DEALLOCATE (c15_ini, STAT=status(108))
        DEALLOCATE (pH_ini, STAT=status(109))
        DEALLOCATE (c17_ini, STAT=status(110))
        DEALLOCATE (NINb_ini, STAT=status(111))
        DEALLOCATE (NIPb_ini, STAT=status(112))
        !gp 03-Apr-08 status(113) is  used abouve for krea2F

        !gp 21-Nov-06
        !WRITE(10,*) 'done thru classreadfile DEALLOCATE'

        !light data

        !gp 13-Feb-06
        !READ(8,*) PAR, kep, kela, kenla, kess, kepom
        READ(8,*) PAR, kep, kela, kenla, kess, kepom, kemac

        !gp 24-Jun-09
        !gp 16-Jul-08
        !READ(8,*) solarMethod, nfacBras, atcRyanStolz, longatMethod, fUwMethod
        !READ(8,*) solarMethod, nfacBras, atcRyanStolz, longatMethod, kbrut, fUwMethod
        READ(8,*) solarMethod, nfacBras, atcRyanStolz, longatMethod, kbrut, fUwMethod, KCL1, KCL2

        !light and heat module

        !gp 13-Feb-06
        !CALL Light_(PAR, kep, kela, kenla, kess, kepom, longatMethod, fUwMethod)

        !gp 24-Jun-09
        !gp 16-Jul-08
        !CALL Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, fUwMethod)
        !CALL Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, kbrut, fUwMethod)
        CALL Light_(PAR, kep, kela, kenla, kess, kepom, kemac, longatMethod, kbrut, fUwMethod, KCL1, KCL2)

        READ(8,*) nptIn

        IF (nptIn>=1) THEN		!if any point load/source
            ALLOCATE (PtName(nptIn), STAT=status(1)); ALLOCATE (xptt(nptIn), STAT=status(2))
            ALLOCATE (Qptta(nptIn), STAT=status(3));	ALLOCATE (Qptt(nptIn), STAT=status(4))
            ALLOCATE (TepttMean(nptIn), STAT=status(5));ALLOCATE (TepttAmp(nptIn), STAT=status(6))
            ALLOCATE (TepttMaxTime(nptIn), STAT=status(7))
            ALLOCATE (cpttMean(nptIn, nv-2), STAT=status(8)); ALLOCATE(cpttAmp(nptIn, nv-2), STAT=status(9))
            ALLOCATE(cpttMaxTime(nptIn, nv-2), STAT=status(10)); ALLOCATE (phpttMean(nptIn), STAT=status(11))
            ALLOCATE (phpttAmp(nptIn), STAT=status(12));ALLOCATE(phpttMaxTime(nptIn), STAT=status(13))
            !Point sources

            DO i=1 , 13

                IF (status(i)==1) THEN
                    STOP 'Class_ReadFile:ReadInputfile. dynamic memory allocation failed!'
                END IF
            END DO
            DO i = 1, nptIn
                READ(8,*) PtName(i), xptt(i), Qptta(i), Qptt(i), TepttMean(i), TepttAmp(i), TepttMaxTime(i)
                DO k = 1, nv - 2
                    READ(8,*) cpttMean(i, k), cpttAmp(i, k), cpttMaxTime(i, k)
                END DO
                READ(8,*) phpttMean(i), phpttAmp(i), phpttMaxTime(i)
            END DO

        END IF


        !Diffuse sources
        READ(8,*) ndiffIn
        IF (ndiffIn>=1) THEN		!if any nonpoint load/source
            ALLOCATE (xdup(ndiffIn), STAT=status(1));	ALLOCATE(xddn(ndiffIn), STAT=status(2))
            ALLOCATE(Qdifa(ndiffIn), STAT=status(3));	ALLOCATE (Qdif(ndiffIn), STAT=status(4))
            ALLOCATE(Tedif(ndiffIn), STAT=status(5));	ALLOCATE(DiffName(ndiffIn), STAT=status(6))
            ALLOCATE(pHind(ndiffIn), STAT=status(7)); ALLOCATE(cdif(ndiffIn, nv-2), STAT=status(8))
            DO i=1, 8
                IF (status(i)==1) THEN
                    STOP 'Class_ReadFile:ReadInputfile. dynamic memory allocation failed!'
                END IF
            END DO
            DO i = 1, ndiffIn
                READ(8,*) DiffName(i), xdup(i), xddn(i), Qdifa(i), Qdif(i), Tedif(i)
                DO k = 1, nv - 2
                    READ(8,*) cdif(i, k)
                END DO
                READ(8,*) pHind(i)
            END DO
        END IF
        CALL sourceIn_(nRch, nptIn, ndiffIn, hydrau%flag, Topo, hydrau, PtName, xptt, Qptta, Qptt, TepttMean, &
            TepttAmp, TepttMaxTime, cpttMean, cpttAmp, cpttMaxTime, phpttMean, &
            phpttAmp, phpttMaxTime, DiffName, xdup, xddn, Qdifa, Qdif, Tedif, cdif, pHind)

        IF (ndiffIn>=1) THEN
            DEALLOCATE(pHind, STAT=status(1));	DEALLOCATE(DiffName, STAT=status(2));
            DEALLOCATE(Tedif, STAT=status(3));	DEALLOCATE (Qdif, STAT=status(4));
            DEALLOCATE(Qdifa, STAT=status(5));	DEALLOCATE(xddn, STAT=status(6))
            DEALLOCATE(xdup, STAT=status(7))
        END IF
        IF (nptIn>=1) THEN
            DEALLOCATE(phpttMaxTime, STAT=status(13));	DEALLOCATE (phpttAmp, STAT=status(12))
            DEALLOCATE (phpttMean, STAT=status(11));		DEALLOCATE(cpttMaxTime, STAT=status(10))
            DEALLOCATE(cpttAmp, STAT=status(9));				DEALLOCATE (cpttMean, STAT=status(8))
            DEALLOCATE (TepttMaxTime, STAT=status(7));	DEALLOCATE (TepttAmp, STAT=status(6))
            DEALLOCATE (TepttMean, STAT=status(5));			DEALLOCATE (Qptt, STAT=status(4))
            DEALLOCATE (Qptta, STAT=status(3));					DEALLOCATE (xptt, STAT=status(2))
            DEALLOCATE (PtName, STAT=status(1))
        END IF

        !Rates
        READ(8,*) vss, mgC, mgN, mgP, mgD, mgA
        READ(8,*) tka, roc, ron
        READ(8,*) Ksocf, Ksona, Ksodn, Ksop, Ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von
        READ(8,*) kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        READ(8,*) kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat

        !gp 03-Apr-08
        !READ(8,*) khnx, va, typeF, kgaF, tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, abmax
        READ(8,*) khnx, va, typeF, kgaF, tkgaF, krea1F, krea2F, tkreaF, kexaF, tkexaF, kdeaF, abmax

        READ(8,*) tkdeaF, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt

        !gp 26-Jan-06
        !gp 10-Jan-05 READ(8,*) NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP
        !READ(8,*) NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqNratio, KqPratio
        READ(8,*) NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqNratio, KqPratio, NUpWCfrac, PUpWCfrac

        KqN = KqNratio * NINbmin
        KqP = KqPratio * NIPbmin

        !gp 30-Nov-04 READ(8,*) kpath, tkpath, vpath
        READ(8,*) kpath, tkpath, vpath, apath, kgen, tkgen, vgen, useGenericAsCOD		!gp 08-Dec-04

        READ(8,*) xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7

        !gp 15-Nov-04 separate hyporheic rates from reaeration and add level 2 rates and HCO3- use
        !gp 03-Nov-04 READ(8,*) kai, kawindmethod	!oxygen reaeration option
        !gp READ(8,*) kai, kawindmethod, typeH, kgaH, tkgaH, kscH, xdum8, kinhcH
        READ(8,*) kai, kawindmethod
        READ(8,*) hco3use, hco3useF
        READ(8,*) typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax

        !gp 17-Jan-06
        !gp 17-Nov-04 CALL MakeHydraulics(nRch, hydrau, downstreamBoundary, kai)
        !CALL MakeHydraulics(nRch, hydrau, downstreamBoundary, kai, Topo%geoMethod)

        !gp 08-Feb-06
        !stochRate= Rates_(mgC,mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, &
        !							Ksodn, Ksop, Ksob, khc, tkhc, kdcs, tkdcs,kdc, tkdc, khn, tkhn, von, &
        !							kn, tkn, ki, tki, vdi, tvdi, khp, tkhp,vop, vip, kspi, kga, &
        !							tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, &
        !							va, typeF, kgaF, tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, &
        !							abmax, tkdeaF, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, &
        !							NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP, &
        !							kpath, tkpath, vpath, pco2, xdum1, xdum2, xdum3, xdum4, xdum5, &
        !							xdum6,xdum7, kai, kawindmethod, &
        !							hco3use, hco3useF, &
        !							typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
        !							kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &	!gp 15-Nov-04
        !							apath, kgen, tkgen, vgen, useGenericAsCOD, &				!gp 08-Dec-04
        !							NUpWCfrac, PUpWCfrac)										!gp 26-Jan-06

        !gp 03-Apr-08
        !stochRate= Rates_(nRch, hydrau, mgC,mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, &
        !							Ksodn, Ksop, Ksob, khc, tkhc, kdcs, tkdcs,kdc, tkdc, khn, tkhn, von, &
        !							kn, tkn, ki, tki, vdi, tvdi, khp, tkhp,vop, vip, kspi, kga, &
        !							tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, &
        !							va, typeF, kgaF, tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, &
        !							abmax, tkdeaF, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, &
        !							NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP, &
        !							kpath, tkpath, vpath, pco2, xdum1, xdum2, xdum3, xdum4, xdum5, &
        !							xdum6,xdum7, kai, kawindmethod, &
        !							hco3use, hco3useF, &
        !							typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
        !							kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &
        !							apath, kgen, tkgen, vgen, useGenericAsCOD, &
        !							NUpWCfrac, PUpWCfrac)
        stochRate= Rates_(nRch, hydrau, mgC,mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, &
            Ksodn, Ksop, Ksob, khc, tkhc, kdcs, tkdcs,kdc, tkdc, khn, tkhn, von, &
            kn, tkn, ki, tki, vdi, tvdi, khp, tkhp,vop, vip, kspi, kga, &
            tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, &
            va, typeF, kgaF, tkgaF, krea1F, krea2F, tkreaF, kexaF, tkexaF, kdeaF, &
            abmax, tkdeaF, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, &
            NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP, &
            kpath, tkpath, vpath, pco2, xdum1, xdum2, xdum3, xdum4, xdum5, &
            xdum6,xdum7, kai, kawindmethod, &
            hco3use, hco3useF, &
            typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
            kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &
            apath, kgen, tkgen, vgen, useGenericAsCOD, &
            NUpWCfrac, PUpWCfrac)

        !day saving time
        READ(8,*) dlstime
        siteSolar= sitesolar_(nRch, timezone, solarMethod, nfacBras, atcRyanStolz, dlstime)
        !calculate sun rise and set
        CALL SunriseSunset(nRch, siteSolar, hydrau, system%today)

        !Boundary data
        READ(8,*) downstreamBoundary	!true/false

        !gp 17-Jan-06
        CALL MakeHydraulics(nRch, hydrau, downstreamBoundary, kai, Topo%geoMethod)

        !gp 12-Jan-06
        !WRITE(10,*) downstreamBoundary
        !WRITE(10,*) system%steadystate

        !HEADWATERS
        IF (system%steadystate) THEN
            ALLOCATE(HwIn(0:HRSDAY-1), STAT=status(1))
            IF (status(1)==1) STOP 'Insufficient memory, dynamic allocation failed!'
            READ(8,*) (HwIn(j)%Te, j=0, HRSDAY-1), xxx
            DO j = 1, nv - 2
                READ(8,*) (HwIn(k)%c(j), k=0, HRSDAY-1), xxx
            END DO
            !hourly pH parameters
            READ(8,*) (HwIn(j)%pH, j=0, HRSDAY-1), xxx
        ELSE
            !DYNAMIC SIMULATION
        END IF

        !gp 11-Jan-06 re-write with downstream boundary read same as headwater
        IF (downstreamBoundary) THEN
            IF (system%steadystate) THEN
                ALLOCATE(DBIn(0:HRSDAY-1), STAT=status(1))
                IF (status(1)==1) STOP 'Insufficient memory, dynamic allocation failed!'
                READ(8,*) (DBIn(j)%Te, j=0, HRSDAY-1), xxx

                !gp 12-Jan-06
                !WRITE(9,*) (DBIn(j)%Te, j=0, HRSDAY-1), xxx

                DO j = 1, nv - 2
                    READ(8,*) (DBIn(k)%c(j), k=0, HRSDAY-1), xxx

                    !gp 12-Jan-06
                    !WRITE(10,*) (DBIn(k)%c(j), k=0, HRSDAY-1), xxx

                END DO
                !hourly pH parameters
                READ(8,*) (DBIn(j)%pH, j=0, HRSDAY-1), xxx

                !gp 12-Jan-06
                !WRITE(10,*) (DBIn(j)%pH, j=0, HRSDAY-1), xxx

            ELSE
                !DYNAMIC SIMULATION
            END IF
            !ELSE !downstream boundary is not specified
        END IF

        HW= Headwater_(HwIn)
        DB= Downstream_(DBin,downstreamBoundary)
        IF (downstreamBoundary) DEALLOCATE(DBin, STAT=status(1))
        DEALLOCATE(HwIn, STAT=status(2))

        !gp 12-Jan-06
        !WRITE(10,*) 'done thru read downstreamBoundary'

        !Met data
        IF (system%steadystate) THEN
            ALLOCATE(shadeHH(0:HRSDAY-1,nRch), STAT=status(1))
            ALLOCATE(TaHH(0:HRSDAY-1,nRch), STAT=status(2))
            ALLOCATE(TdHH(0:HRSDAY-1,nRch), STAT=status(3))
            ALLOCATE(UwHH(0:HRSDAY-1,nRch), STAT=status(4))
            ALLOCATE(ccHH(0:HRSDAY-1,nRch), STAT=status(5))

            !gp 16-Jul-08
            ALLOCATE(solarHH(0:HRSDAY-1,nRch), STAT=status(6))

            DO i=1, 5
                IF (status(i)==1) STOP 'Insufficient memory, dynamic allocation failed!'
            END DO

            DO i = 1, nRch
                !DO j = 0, HRSDAY-1
                READ(8,*) (shadeHH(j, i), j = 0, HRSDAY-1), xxx 	!Reach Hourly shade data
                !END DO
                !READ(8,*) xxx
            END DO
            DO i = 1, nRch
                !DO j = 0, HRSDAY-1						!Reach Hourly Air temperature data
                READ(8,*) (TaHH(j, i), j = 0, HRSDAY-1), xxx
                !END DO
                !READ(8,*) xxx
            END DO
            DO i = 1, nRch
                !DO j = 0, HRSDAY-1
                READ(8,*) (TdHH(j, i), j = 0, HRSDAY-1), xxx		!Reach dew point Air temperature data
                !END DO
                !READ(8,*) xxx
            END DO
            DO i = 1, nRch
                !DO j = 0, HRSDAY-1
                READ(8,*) (UwHH(j, i), j = 0, HRSDAY-1), xxx		!Reach wind speed data
                !END DO
                !READ(8,*) xxx
            END DO
            DO i = 1, nRch
                !DO j = 0, HRSDAY-1
                READ(8,*) (ccHH(j, i), j = 0, HRSDAY-1), xxx		!Reach cloud cover data
                !END DO
                !READ(8,*) xxx
            END DO

            !gp 16-Jul-08
            DO i = 1, nRch
                READ(8,*) (solarHH(j, i), j = 0, HRSDAY-1), xxx		!Reach solar data
            END DO

        ELSE															!dyanmic simulation
        END IF

        !gp 16-Jul-08
        !siteMeteo = MeteoData_(nRch, shadeHH, TaHH, TdHH, UwHH, ccHH)
        siteMeteo = t_meteorology(nRch, shadeHH, TaHH, TdHH, UwHH, ccHH, solarHH)

        !gp 16-Jul-08
        DEALLOCATE(solarHH, STAT=status(6))

        DEALLOCATE(ccHH, STAT=status(5))
        DEALLOCATE(UwHH, STAT=status(4))
        DEALLOCATE(TdHH, STAT=status(3))
        DEALLOCATE(TaHH, STAT=status(2))
        DEALLOCATE(shadeHH, STAT=status(1))

        !Data
!	READ(8,*) nteda
!	DO i = 1, nteda
!		READ(8,*) junk, junk, junk, junk !xteda(i), tedaav(i), tedamn(i), tedamx(i)
!	END DO

!	READ(8,*) nhydda
!	DO i = 1, nhydda
!		READ(8,*) junk, junk, junk, junk, junk ! xhydda(i), Qdata(i), Hdata(i), Udata(i), Travdata(i)
!	END DO

!	DO kk = 1, 3
!		READ(8,*) shtnam(kk), nwqd(kk)
!		DO i = 1, nwqd(kk)
!			READ(8,*) junk !dist(i, kk)
!			READ(8,*) junk, junk, junk, junk, junk, junk, junk, junk, junk, junk
!			READ(8,*) junk, junk, junk, junk, junk, junk, junk, junk, junk, junk
!			READ(8,*) junk, junk, junk, junk, junk, junk, junk, junk

        !cwqdat(i, 1, kk), cwqdat(i, 2, kk), cwqdat(i, 3, kk), cwqdat(i, 4, kk), cwqdat(i, 5, kk), _
        !      cwqdat(i, 6, kk), cwqdat(i, 7, kk), cwqdat(i, 8, kk), cwqdat(i, 9, kk), cwqdat(i, 10, kk)
        !READ(8,*) , cwqdat(i, 11, kk), cwqdat(i, 12, kk), cwqdat(i, 13, kk), cwqdat(i, 14, kk), cwqdat(i, 15, kk), _
        !      cwqdat(i, 16, kk), cwqdat(i, 17, kk), cwqdat(i, 18, kk), cwqdat(i, 19, kk), cwqdat(i, 20, kk)
        !READ(8,*) , cwqdat(i, 21, kk), cwqdat(i, 22, kk), cwqdat(i, 23, kk), cwqdat(i, 24, kk), cwqdat(i, 25, kk), _
        !       cwqdat(i, 26, kk), cwqdat(i, 27, kk), cwqdat(i, 28, kk)
!		END DO
!	END DO

!	READ(8,*) nrp
!	READ(8,*) nwqdiur
!	DO i = 1, nwqdiur
!		READ(8,*) junk !distdiur
!		READ(8,*) junk, junk, junk, junk, junk, junk, junk, junk, junk, junk, junk
!		READ(8,*) junk, junk, junk, junk, junk, junk, junk, junk, junk, junk
        !READ(8,*) , distdiur(i)
        !READ(8,*) , cwqdatdiur(i, 1), cwqdatdiur(i, 2), cwqdatdiur(i, 3), cwqdatdiur(i, 4), cwqdatdiur(i, 5), _
        !       cwqdatdiur(i, 6), cwqdatdiur(i, 7), cwqdatdiur(i, 8), cwqdatdiur(i, 9), cwqdatdiur(i, 10), cwqdatdiur(i, 11)
        !READ(8,*) , cwqdatdiur(i, 12), cwqdatdiur(i, 13), cwqdatdiur(i, 14), cwqdatdiur(i, 15), cwqdatdiur(i, 16), _
        !       cwqdatdiur(i, 17), cwqdatdiur(i, 18), cwqdatdiur(i, 19), cwqdatdiur(i, 20), cwqdatdiur(i, 21)
!	END DO

        !gp 21-Nov-06
        !WRITE(10,*) 'done thru classreadfile SUB ReadInputFile'

    END SUBROUTINE ReadInputfile


END MODULE Class_ReadFile

