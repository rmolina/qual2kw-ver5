!classhydraulics.f90
!
!hydraulics data structure

module class_hydraulics
    use nrtype
    use m_rivertopo !only reach number and element number used

    implicit none
    private allocatehydrauarrays, depthmanning
! public hydraulics_, makehydraulics, geometry_type, reach
    type hydraulics_type
        real(dp) :: b =0, bb=0, xl =0 !bottom width
        real(dp) :: hweir =0, bweir =0 !weir height, width
        real(dp) :: alp1=0, bet1=0, alp2=0, bet2 =0 !velocity, stage discharge coefficents
        real(dp) :: ss1=0, ss2=0, s =0 !side and channel slopes
        real(dp) :: nm=0, kaaa=0, ediff=0 !nm - manning's n
        real(dp) :: asb =0, ast=0, asd=0
        real(dp) :: frsed=0, frsod=0
        real(dp) :: elev1=0, elev2=0 !elevation
        real(dp) :: elev=0 !average elevation
        real(dp) :: latr=0, lonr=0 !latitude, longitude radius
        real(dp) :: sodspec=0, jch4spec=0, jnh4spec=0, jsrpspec=0

        real(dp) :: q=0, u=0, trav=0, qpt=0, qpta=0 !flow rate, velocity
        real(dp) :: depth=0, ac=0, vol=0, drop=0
        real(dp) :: kau =0, ka=0, os =0
        character(len=20) :: kaf='Internal' !oxygen reaeration equation
        real(dp) :: ep=0, epout=0
        real(dp) :: qcmd=0, epcmd=0, qptcmd=0, qptacmd=0

        !gp 16-jul-08
        !real(dp) :: sedthermcond=0, sedthermdiff=0, hsedcm=0, &
        ! hypoexchfrac=0, porosity=0, rhocpsed=0, botalg0=0 !gp 11-jan-05
        real(dp) :: sedthermcond=0, sedthermdiff=0, hsedcm=0, &
            hypoexchfrac=0, porosity=0, rhocpsed=0, botalg0=0, skop=0

        real(dp) :: hypoporevol=0 !gp 03-nov-04
        real(dp) :: ehyporheiccmd=0 !gp 15-nov-04

        !gp 07-feb-06
        !since reach specific rates are not compulsory,
        !these values could be null_val
        !null_val means "no value"/ lack of value from users
        !real(dp) :: kaaa= null_val
        real(dp) :: vss = null_val !inorganic suspended solids settling vol
        real(dp) :: khc = null_val !slow cbod hydrolysis rate
        real(dp) :: kdcs = null_val !slow cbod oxidation rate
        real(dp) :: kdc = null_val !fast cbod oxidation rate
        real(dp) :: khn = null_val !organic n hydrolysis rate
        real(dp) :: von = null_val !organic n settling velocity
        real(dp) :: kn = null_val !ammonium nitrification rate
        real(dp) :: ki = null_val !nitrate denitrification
        real(dp) :: vdi = null_val !nitrate sed denitrification transfer coeff
        real(dp) :: khp= null_val !organic p hydrolysis
        real(dp) :: vop= null_val !organic p settling velocity
        real(dp) :: vip= null_val !inorganic p settling velocity
        real(dp) :: kga= null_val !phytoplankton max growth rate
        real(dp) :: krea= null_val !phytoplankton respiration rate
        real(dp) :: kdea= null_val !phytoplankton death rate
        real(dp) :: ksn= null_val !phytoplankton n half-sat
        real(dp) :: ksp= null_val !phytoplankton p half-sat
        real(dp) :: isat= null_val !phytoplankton light sat
        real(dp) :: khnx= null_val !phytoplankton ammonia preference
        real(dp) :: va= null_val !phytoplankton settling velocity
        !real(dp) :: botalg0= null_val !bottom plant initial bionass
        real(dp) :: kgaf= null_val !bottom plant max growth rate
        real(dp) :: abmax= null_val !bottom plant first-order carrying capacity

        !gp 03-apr-08
        !real(dp) :: kreaf= null_val !bottom plant respiration rate
        real(dp) :: krea1f= null_val !bottom plant basal respiration rate
        real(dp) :: krea2f= null_val !bottom plant photo respiration rate

        real(dp) :: kexaf= null_val !bottom plant excretion rate
        real(dp) :: kdeaf= null_val !bottom plant death rate
        real(dp) :: ksnf= null_val !bottom plant external n half-sat
        real(dp) :: kspf= null_val !bottom plant external p half-sat
        real(dp) :: isatf= null_val !bottom plant light sat
        real(dp) :: khnxf= null_val !bottom plant ammonia preference
        real(dp) :: ninbmin= null_val !bottom plant subistence quota for n
        real(dp) :: nipbmin= null_val !bottom plant subistence quota for p
        real(dp) :: ninbupmax= null_val !bottom plant max uptake rate for n
        real(dp) :: nipbupmax= null_val !bottom plant max uptake rate for p
        real(dp) :: kqn= null_val !bottom plant internal n half-sat
        real(dp) :: kqp= null_val !bottom plant internal p half-sat
        real(dp) :: nupwcfrac= null_val !bottom plant n uptake fraction from water column
        real(dp) :: pupwcfrac= null_val !bottom plant p uptake fraction from water column
        real(dp) :: kdt= null_val !pom dissolution rate
        real(dp) :: vdt= null_val !pom settling velocity
        real(dp) :: kpath= null_val !pathogen dieoff rate
        real(dp) :: vpath= null_val !pathogen settling velocity
        real(dp) :: apath= null_val !pathogen light alpha
        real(dp) :: kgah= null_val !hyporheic heterotrophs max growth rate
        real(dp) :: ksch= null_val !hyporheic heterotrophs cbod half-sat
        real(dp) :: kinhch= null_val !hyporheic heterotrophs o2 inhibition
        real(dp) :: kreah= null_val !hyporheic heterotrophs respiration rate
        real(dp) :: kdeah= null_val !hyporheic heterotrophs death rate
        real(dp) :: ksnh= null_val !hyporheic heterotrophs n half-sat
        real(dp) :: ksph= null_val !hyporheic heterotrophs p half-sat
        real(dp) :: khnxh= null_val !hyporheic heterotrophs ammonia preference
        real(dp) :: ahmax= null_val !hyporheic heterotrophs first-order carrying capacity
        real(dp) :: kgen= null_val !generic constituent dissolution rate
        real(dp) :: vgen= null_val !generic constituent settling velocity

        !gp 21-nov-06 initial conditions of state variables
        real(dp) :: te_ini= null_val
        real(dp) :: c01_ini= null_val
        real(dp) :: c02_ini= null_val
        real(dp) :: c03_ini= null_val
        real(dp) :: c04_ini= null_val
        real(dp) :: c05_ini= null_val
        real(dp) :: c06_ini= null_val
        real(dp) :: c07_ini= null_val
        real(dp) :: c08_ini= null_val
        real(dp) :: c09_ini= null_val
        real(dp) :: c10_ini= null_val
        real(dp) :: c11_ini= null_val
        real(dp) :: c12_ini= null_val
        real(dp) :: c13_ini= null_val
        real(dp) :: c14_ini= null_val
        real(dp) :: c15_ini= null_val
        real(dp) :: ph_ini= null_val
        real(dp) :: c17_ini= null_val
        real(dp) :: ninb_ini= null_val
        real(dp) :: nipb_ini= null_val

    end type hydraulics_type

    type riverhydraulics_type
        type(hydraulics_type), dimension(:), pointer :: reach
        integer(i4b) :: flag =1
    end type riverhydraulics_type

! type(riverhydraulics_type) hydrau

! type(geometry_type), dimension(:), pointer :: reach
! type(hydraulics_type), dimension(:), pointer :: hydrau

! character(len=20) kai

! type(geometry_type), allocatable :: reach(:)
! type(hydraulics_type), allocatable :: hydrau(:

contains
    !/* public functions */

    !gp 07-feb-06
    !function hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, q, bb, &
    ! ss1, ss2, s, nm, alp1, bet1, alp2, bet2, ediff, kaaa, &
    ! frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, hweir, bweir, &
    ! sedthermcond, sedthermdiff, hsedcm, &
    ! hypoexchfrac, porosity, rhocpsed, botalg0) result(hydrau) !gp 11-jan-05

    !gp 21-nov-06
    !function hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, q, bb, &
    ! ss1, ss2, s, nm, alp1, bet1, alp2, bet2, ediff, &
    ! frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, hweir, bweir, &
    ! sedthermcond, sedthermdiff, hsedcm, &
    ! hypoexchfrac, porosity, rhocpsed, &
    ! kaaa, vss, khc, &
    ! kdcs, kdc, khn, &
    ! von, kn , ki , &
    ! vdi, khp, vop, &
    ! vip, kga, krea, &
    ! kdea, ksn, ksp, &
    ! isat, khnx, va, &
    ! botalg0, kgaf, abmax, &
    ! kreaf, kexaf, kdeaf, &
    ! ksnf, kspf, isatf, &
    ! khnxf, ninbmin, nipbmin, &
    ! ninbupmax, nipbupmax, kqn, &
    ! kqp, nupwcfrac, pupwcfrac, &
    ! kdt, vdt, kpath, &
    ! vpath, apath, kgah, &
    ! ksch, kinhch, kreah, &
    ! kdeah, ksnh, ksph, &
    ! khnxh, ahmax, &
    ! kgen, vgen) result(hydrau)

    !gp 03-apr-08
    !function hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, q, bb, &
    ! ss1, ss2, s, nm, alp1, bet1, alp2, bet2, ediff, &
    ! frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, hweir, bweir, &
    ! sedthermcond, sedthermdiff, hsedcm, &
    ! hypoexchfrac, porosity, rhocpsed, &
    ! kaaa, vss, khc, &
    ! kdcs, kdc, khn, &
    ! von, kn , ki , &
    ! vdi, khp, vop, &
    ! vip, kga, krea, &
    ! kdea, ksn, ksp, &
    ! isat, khnx, va, &
    ! kgaf, abmax, &
    ! kreaf, kexaf, kdeaf, &
    ! ksnf, kspf, isatf, &
    ! khnxf, ninbmin, nipbmin, &
    ! ninbupmax, nipbupmax, kqn, &
    ! kqp, nupwcfrac, pupwcfrac, &
    ! kdt, vdt, kpath, &
    ! vpath, apath, kgah, &
    ! ksch, kinhch, kreah, &
    ! kdeah, ksnh, ksph, &
    ! khnxh, ahmax, &
    ! kgen, vgen, &
    ! te_ini, c01_ini, c02_ini, c03_ini, &
    ! c04_ini, c05_ini, c06_ini, &
    ! c07_ini, c08_ini, c09_ini, &
    ! c10_ini, c11_ini, c12_ini, &
    ! c13_ini, c14_ini, c15_ini, &
    ! ph_ini, c17_ini, ninb_ini, nipb_ini) result(hydrau)

    !gp 16-jul-08
    !function hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, q, bb, &
    ! ss1, ss2, s, nm, alp1, bet1, alp2, bet2, ediff, &
    ! frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, hweir, bweir, &
    ! sedthermcond, sedthermdiff, hsedcm, &
    ! hypoexchfrac, porosity, rhocpsed, &
    ! kaaa, vss, khc, &
    ! kdcs, kdc, khn, &
    ! von, kn , ki , &
    ! vdi, khp, vop, &
    ! vip, kga, krea, &
    ! kdea, ksn, ksp, &
    ! isat, khnx, va, &
    ! kgaf, abmax, &
    ! krea1f, krea2f, kexaf, kdeaf, &
    ! ksnf, kspf, isatf, &
    ! khnxf, ninbmin, nipbmin, &
    ! ninbupmax, nipbupmax, kqn, &
    ! kqp, nupwcfrac, pupwcfrac, &
    ! kdt, vdt, kpath, &
    ! vpath, apath, kgah, &
    ! ksch, kinhch, kreah, &
    ! kdeah, ksnh, ksph, &
    ! khnxh, ahmax, &
    ! kgen, vgen, &
    ! te_ini, c01_ini, c02_ini, c03_ini, &
    ! c04_ini, c05_ini, c06_ini, &
    ! c07_ini, c08_ini, c09_ini, &
    ! c10_ini, c11_ini, c12_ini, &
    ! c13_ini, c14_ini, c15_ini, &
    ! ph_ini, c17_ini, ninb_ini, nipb_ini) result(hydrau)
    function hydraulics_(nr, xrdn, elev1, elev2, latd, latm, lats, lond, lonm, lons, q, bb, &
        ss1, ss2, s, nm, alp1, bet1, alp2, bet2, ediff, &
        frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, hweir, bweir, &
        sedthermcond, sedthermdiff, hsedcm, &
        hypoexchfrac, porosity, rhocpsed, skop, &
        kaaa, vss, khc, &
        kdcs, kdc, khn, &
        von, kn , ki , &
        vdi, khp, vop, &
        vip, kga, krea, &
        kdea, ksn, ksp, &
        isat, khnx, va, &
        kgaf, abmax, &
        krea1f, krea2f, kexaf, kdeaf, &
        ksnf, kspf, isatf, &
        khnxf, ninbmin, nipbmin, &
        ninbupmax, nipbupmax, kqn, &
        kqp, nupwcfrac, pupwcfrac, &
        kdt, vdt, kpath, &
        vpath, apath, kgah, &
        ksch, kinhch, kreah, &
        kdeah, ksnh, ksph, &
        khnxh, ahmax, &
        kgen, vgen, &
        te_ini, c01_ini, c02_ini, c03_ini, &
        c04_ini, c05_ini, c06_ini, &
        c07_ini, c08_ini, c09_ini, &
        c10_ini, c11_ini, c12_ini, &
        c13_ini, c14_ini, c15_ini, &
        ph_ini, c17_ini, ninb_ini, nipb_ini) result(hydrau)

        type(riverhydraulics_type) hydrau
        integer(i4b), intent(in) :: nr
        real(dp), intent(in) :: xrdn(0:), elev1(0:), elev2(0:), latd(0:), latm(0:)

        !gp 07-feb-06
        !real(dp), intent(in) :: lats(0:), lond(0:), lonm(0:), lons(0:), q(0:), bb(0:), &
        ! ss1(0:), ss2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
        ! alp2(0:), bet2(0:), ediff(0:), kaaa(0:), frsed(0:), &
        ! frsod(0:), sodspec(0:), jch4spec(0:), jnh4spec(0:), &
        ! jsrpspec(0:), hweir(0:), bweir(0:), &
        ! sedthermcond(0:), sedthermdiff(0:), hsedcm(0:), &
        ! hypoexchfrac(0:), porosity(0:), rhocpsed(0:), botalg0(0:) !gp 11-jan-05

        !gp 21-nov-06
        !real(dp), intent(in) :: lats(0:), lond(0:), lonm(0:), lons(0:), q(0:), bb(0:), &
        ! ss1(0:), ss2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
        ! alp2(0:), bet2(0:), ediff(0:), frsed(0:), &
        ! frsod(0:), sodspec(0:), jch4spec(0:), jnh4spec(0:), &
        ! jsrpspec(0:), hweir(0:), bweir(0:), &
        ! sedthermcond(0:), sedthermdiff(0:), hsedcm(0:), &
        ! hypoexchfrac(0:), porosity(0:), rhocpsed(0:), &
        ! kaaa(0:), vss(0:), khc(0:), &
        ! kdcs(0:), kdc(0:), khn(0:), &
        ! von(0:), kn(0:) , ki(0:) , &
        ! vdi(0:), khp(0:), vop(0:), &
        ! vip(0:), kga(0:), krea(0:), &
        ! kdea(0:), ksn(0:), ksp(0:), &
        ! isat(0:), khnx(0:), va(0:), &
        ! botalg0(0:), kgaf(0:), abmax(0:), &
        ! kreaf(0:), kexaf(0:), kdeaf(0:), &
        ! ksnf(0:), kspf(0:), isatf(0:), &
        ! khnxf(0:), ninbmin(0:), nipbmin(0:), &
        ! ninbupmax(0:), nipbupmax(0:), kqn(0:), &
        ! kqp(0:), nupwcfrac(0:), pupwcfrac(0:), &
        ! kdt(0:), vdt(0:), kpath(0:), &
        ! vpath(0:), apath(0:), kgah(0:), &
        ! ksch(0:), kinhch(0:), kreah(0:), &
        ! kdeah(0:), ksnh(0:), ksph(0:), &
        ! khnxh(0:), ahmax(0:), &
        ! kgen(0:), vgen(0:)

        !gp 03-apr-08
        !real(dp), intent(in) :: lats(0:), lond(0:), lonm(0:), lons(0:), q(0:), bb(0:), &
        ! ss1(0:), ss2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
        ! alp2(0:), bet2(0:), ediff(0:), frsed(0:), &
        ! frsod(0:), sodspec(0:), jch4spec(0:), jnh4spec(0:), &
        ! jsrpspec(0:), hweir(0:), bweir(0:), &
        ! sedthermcond(0:), sedthermdiff(0:), hsedcm(0:), &
        ! hypoexchfrac(0:), porosity(0:), rhocpsed(0:), &
        ! kaaa(0:), vss(0:), khc(0:), &
        ! kdcs(0:), kdc(0:), khn(0:), &
        ! von(0:), kn(0:) , ki(0:) , &
        ! vdi(0:), khp(0:), vop(0:), &
        ! vip(0:), kga(0:), krea(0:), &
        ! kdea(0:), ksn(0:), ksp(0:), &
        ! isat(0:), khnx(0:), va(0:), &
        ! kgaf(0:), abmax(0:), &
        ! kreaf(0:), kexaf(0:), kdeaf(0:), &
        ! ksnf(0:), kspf(0:), isatf(0:), &
        ! khnxf(0:), ninbmin(0:), nipbmin(0:), &
        ! ninbupmax(0:), nipbupmax(0:), kqn(0:), &
        ! kqp(0:), nupwcfrac(0:), pupwcfrac(0:), &
        ! kdt(0:), vdt(0:), kpath(0:), &
        ! vpath(0:), apath(0:), kgah(0:), &
        ! ksch(0:), kinhch(0:), kreah(0:), &
        ! kdeah(0:), ksnh(0:), ksph(0:), &
        ! khnxh(0:), ahmax(0:), &
        ! kgen(0:), vgen(0:), &
        ! te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
        ! c04_ini(0:), c05_ini(0:), c06_ini(0:), &
        ! c07_ini(0:), c08_ini(0:), c09_ini(0:), &
        ! c10_ini(0:), c11_ini(0:), c12_ini(0:), &
        ! c13_ini(0:), c14_ini(0:), c15_ini(0:), &
        ! ph_ini(0:), c17_ini(0:), ninb_ini(0:), nipb_ini(0:)

        !gp 16-jul-08
        !real(dp), intent(in) :: lats(0:), lond(0:), lonm(0:), lons(0:), q(0:), bb(0:), &
        ! ss1(0:), ss2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
        ! alp2(0:), bet2(0:), ediff(0:), frsed(0:), &
        ! frsod(0:), sodspec(0:), jch4spec(0:), jnh4spec(0:), &
        ! jsrpspec(0:), hweir(0:), bweir(0:), &
        ! sedthermcond(0:), sedthermdiff(0:), hsedcm(0:), &
        ! hypoexchfrac(0:), porosity(0:), rhocpsed(0:), &
        ! kaaa(0:), vss(0:), khc(0:), &
        ! kdcs(0:), kdc(0:), khn(0:), &
        ! von(0:), kn(0:) , ki(0:) , &
        ! vdi(0:), khp(0:), vop(0:), &
        ! vip(0:), kga(0:), krea(0:), &
        ! kdea(0:), ksn(0:), ksp(0:), &
        ! isat(0:), khnx(0:), va(0:), &
        ! kgaf(0:), abmax(0:), &
        ! krea1f(0:), krea2f(0:), kexaf(0:), kdeaf(0:), &
        ! ksnf(0:), kspf(0:), isatf(0:), &
        ! khnxf(0:), ninbmin(0:), nipbmin(0:), &
        ! ninbupmax(0:), nipbupmax(0:), kqn(0:), &
        ! kqp(0:), nupwcfrac(0:), pupwcfrac(0:), &
        ! kdt(0:), vdt(0:), kpath(0:), &
        ! vpath(0:), apath(0:), kgah(0:), &
        ! ksch(0:), kinhch(0:), kreah(0:), &
        ! kdeah(0:), ksnh(0:), ksph(0:), &
        ! khnxh(0:), ahmax(0:), &
        ! kgen(0:), vgen(0:), &
        ! te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
        ! c04_ini(0:), c05_ini(0:), c06_ini(0:), &
        ! c07_ini(0:), c08_ini(0:), c09_ini(0:), &
        ! c10_ini(0:), c11_ini(0:), c12_ini(0:), &
        ! c13_ini(0:), c14_ini(0:), c15_ini(0:), &
        ! ph_ini(0:), c17_ini(0:), ninb_ini(0:), nipb_ini(0:)
        real(dp), intent(in) :: lats(0:), lond(0:), lonm(0:), lons(0:), q(0:), bb(0:), &
            ss1(0:), ss2(0:), s(0:), nm(0:) , alp1(0:), bet1(0:), &
            alp2(0:), bet2(0:), ediff(0:), frsed(0:), &
            frsod(0:), sodspec(0:), jch4spec(0:), jnh4spec(0:), &
            jsrpspec(0:), hweir(0:), bweir(0:), &
            sedthermcond(0:), sedthermdiff(0:), hsedcm(0:), &
            hypoexchfrac(0:), porosity(0:), rhocpsed(0:), skop(0:), &
            kaaa(0:), vss(0:), khc(0:), &
            kdcs(0:), kdc(0:), khn(0:), &
            von(0:), kn(0:) , ki(0:) , &
            vdi(0:), khp(0:), vop(0:), &
            vip(0:), kga(0:), krea(0:), &
            kdea(0:), ksn(0:), ksp(0:), &
            isat(0:), khnx(0:), va(0:), &
            kgaf(0:), abmax(0:), &
            krea1f(0:), krea2f(0:), kexaf(0:), kdeaf(0:), &
            ksnf(0:), kspf(0:), isatf(0:), &
            khnxf(0:), ninbmin(0:), nipbmin(0:), &
            ninbupmax(0:), nipbupmax(0:), kqn(0:), &
            kqp(0:), nupwcfrac(0:), pupwcfrac(0:), &
            kdt(0:), vdt(0:), kpath(0:), &
            vpath(0:), apath(0:), kgah(0:), &
            ksch(0:), kinhch(0:), kreah(0:), &
            kdeah(0:), ksnh(0:), ksph(0:), &
            khnxh(0:), ahmax(0:), &
            kgen(0:), vgen(0:), &
            te_ini(0:), c01_ini(0:), c02_ini(0:), c03_ini(0:), &
            c04_ini(0:), c05_ini(0:), c06_ini(0:), &
            c07_ini(0:), c08_ini(0:), c09_ini(0:), &
            c10_ini(0:), c11_ini(0:), c12_ini(0:), &
            c13_ini(0:), c14_ini(0:), c15_ini(0:), &
            ph_ini(0:), c17_ini(0:), ninb_ini(0:), nipb_ini(0:)

        integer(i4b) i, j

        call allocatehydrauarrays(nr, hydrau)

        hydrau%reach(0)%q=q(0)* 86400.0_dp

        hydrau%reach(0)%nm= nm(0)

        if (xrdn(2)>xrdn(1)) then
            hydrau%flag=1
        else
            hydrau%flag=-1
        end if

        do i=0, nr
            !hydrau%reach(i)%q=q(i)* 86400.0
            if (i>0) hydrau%reach(i)%xl = abs(xrdn(i)-xrdn(i-1))*1000.0_dp
            hydrau%reach(i)%elev1 = elev1(i)
            hydrau%reach(i)%elev2 = elev2(i)
            hydrau%reach(i)%elev = (elev1(i)+elev2(i))/2
            hydrau%reach(i)%ss1 = ss1(i)
            hydrau%reach(i)%ss2 = ss2(i)
            hydrau%reach(i)%s = s(i)
            if (frsod(i) <= 0) then !fraction of bottom sod
                hydrau%reach(i)%frsod = 0.00001_dp
            else
                hydrau%reach(i)%frsod = frsod(i)
            end if

            if (frsed(i) <= 0) then !fraction of bottom algae coverage
                hydrau%reach(i)%frsed = 0.00001_dp
            else
                hydrau%reach(i)%frsed = frsed(i)
            end if
            hydrau%reach(i)%latr = latd(i) + (latm(i)+lats(i)/60.0_dp)/60.0_dp
            hydrau%reach(i)%lonr = lond(i) + (lonm(i)+lons(i)/60.0_dp)/60.0_dp
            hydrau%reach(i)%nm = nm(i)
            hydrau%reach(i)%hweir = hweir(i)
            hydrau%reach(i)%bweir = bweir(i)
            if (((hweir(i) > 0) .and. (bweir(i) <= 0)) .or. &
                ((bweir(i) > 0) .and. (hweir(i) <= 0))) then
                stop 'To include a weir, both height and width must be greater than zero'
            end if
            hydrau%reach(i)%bb = bb(i)
            hydrau%reach(i)%b = bb(i)

            if (((hweir(i) > 0).and. (bweir(i) > 0)) .and. (bb(i) <= 0)) then
                stop 'If a weir is used, the width of the reach upstream on the weir must be entered'
            end if

            hydrau%reach(i)%alp1 = alp1(i)
            hydrau%reach(i)%alp2 = alp2(i)
            hydrau%reach(i)%bet1 = bet1(i)
            hydrau%reach(i)%bet2 = bet2 (i)
            hydrau%reach(i)%ediff = ediff(i) !prescribed dispersion rate

            !gp 07-feb-06
            !hydrau%reach(i)%kaaa = kaaa(i) !prescribed reaeration rate

            hydrau%reach(i)%sodspec = sodspec(i) !prescribed sod
            hydrau%reach(i)%jch4spec=jch4spec(i)
            hydrau%reach(i)%jnh4spec=jnh4spec(i)
            hydrau%reach(i)%jsrpspec=jsrpspec(i)

            hydrau%reach(i)%sedthermcond = sedthermcond(i) !gp 20-oct-04
            hydrau%reach(i)%sedthermdiff = sedthermdiff(i) !gp 20-oct-04
            hydrau%reach(i)%hsedcm = hsedcm(i) !gp 20-oct-04
            hydrau%reach(i)%hypoexchfrac = hypoexchfrac(i) !gp 20-oct-04
            hydrau%reach(i)%porosity = porosity(i) !gp 20-oct-04
            hydrau%reach(i)%rhocpsed = rhocpsed(i) !gp 20-oct-04

            !gp 16-jul-08
            hydrau%reach(i)%skop = skop(i)

            !gp 07-feb-06
            !hydrau%reach(i)%botalg0 = botalg0(i) !gp 11-jan-05
            hydrau%reach(i)%kaaa = kaaa(i) !prescribed reaeration rate
            hydrau%reach(i)%vss = vss(i) !inorganic suspended solids settling vol
            hydrau%reach(i)%khc = khc(i) !slow cbod hydrolysis rate
            hydrau%reach(i)%kdcs = kdcs(i) !slow cbod oxidation rate
            hydrau%reach(i)%kdc = kdc(i) !fast cbod oxidation rate
            hydrau%reach(i)%khn = khn(i) !organic n hydrolysis rate
            hydrau%reach(i)%von = von(i) !organic n settling velocity
            hydrau%reach(i)%kn = kn(i) !ammonium nitrification rate
            hydrau%reach(i)%ki = ki(i) !nitrate denitrification
            hydrau%reach(i)%vdi = vdi(i) !nitrate sed denitrification transfer coeff
            hydrau%reach(i)%khp = khp(i) !organic p hydrolysis
            hydrau%reach(i)%vop = vop(i) !organic p settling velocity
            hydrau%reach(i)%vip = vip(i) !inorganic p settling velocity
            hydrau%reach(i)%kga = kga(i) !phytoplankton max growth rate
            hydrau%reach(i)%krea = krea(i) !phytoplankton respiration rate
            hydrau%reach(i)%kdea = kdea(i) !phytoplankton death rate
            hydrau%reach(i)%ksn = ksn(i) !phytoplankton n half-sat
            hydrau%reach(i)%ksp = ksp(i) !phytoplankton p half-sat
            hydrau%reach(i)%isat = isat(i) !phytoplankton light sat
            hydrau%reach(i)%khnx = khnx(i) !phytoplankton ammonia preference
            hydrau%reach(i)%va = va(i) !phytoplankton settling velocity

            !gp 21-nov-06
            !hydrau%reach(i)%botalg0 = botalg0(i) !bottom plant initial bionass

            hydrau%reach(i)%kgaf = kgaf(i) !bottom plant max growth rate
            hydrau%reach(i)%abmax = abmax(i) !bottom plant first-order carrying capacity

            !gp 03-apr-08
            !hydrau%reach(i)%kreaf = kreaf(i) !bottom plant respiration rate
            hydrau%reach(i)%krea1f = krea1f(i) !bottom plant basal respiration rate
            hydrau%reach(i)%krea2f = krea2f(i) !bottom plant photo respiration rate

            hydrau%reach(i)%kexaf = kexaf(i) !bottom plant excretion rate
            hydrau%reach(i)%kdeaf = kdeaf(i) !bottom plant death rate
            hydrau%reach(i)%ksnf = ksnf(i) !bottom plant external n half-sat
            hydrau%reach(i)%kspf = kspf(i) !bottom plant external p half-sat
            hydrau%reach(i)%isatf = isatf(i) !bottom plant light sat
            hydrau%reach(i)%khnxf = khnxf(i) !bottom plant ammonia preference
            hydrau%reach(i)%ninbmin = ninbmin(i) !bottom plant subistence quota for n
            hydrau%reach(i)%nipbmin = nipbmin(i) !bottom plant subistence quota for p
            hydrau%reach(i)%ninbupmax = ninbupmax(i) !bottom plant max uptake rate for n
            hydrau%reach(i)%nipbupmax = nipbupmax(i) !bottom plant max uptake rate for p
            hydrau%reach(i)%kqn = kqn(i) !bottom plant internal n half-sat
            hydrau%reach(i)%kqp = kqp(i) !bottom plant internal p half-sat
            hydrau%reach(i)%nupwcfrac = nupwcfrac(i) !bottom plant n uptake fraction from water column
            hydrau%reach(i)%pupwcfrac = pupwcfrac(i) !bottom plant p uptake fraction from water column
            hydrau%reach(i)%kdt = kdt(i) !pom dissolution rate
            hydrau%reach(i)%vdt = vdt(i) !pom settling velocity
            hydrau%reach(i)%kpath = kpath(i) !pathogen dieoff rate
            hydrau%reach(i)%vpath = vpath(i) !pathogen settling velocity
            hydrau%reach(i)%apath = apath(i) !pathogen light alpha
            hydrau%reach(i)%kgah = kgah(i) !hyporheic heterotrophs max growth rate
            hydrau%reach(i)%ksch = ksch(i) !hyporheic heterotrophs cbod half-sat
            hydrau%reach(i)%kinhch = kinhch(i) !hyporheic heterotrophs o2 inhibition
            hydrau%reach(i)%kreah = kreah(i) !hyporheic heterotrophs respiration rate
            hydrau%reach(i)%kdeah = kdeah(i) !hyporheic heterotrophs death rate
            hydrau%reach(i)%ksnh = ksnh(i) !hyporheic heterotrophs n half-sat
            hydrau%reach(i)%ksph = ksph(i) !hyporheic heterotrophs p half-sat
            hydrau%reach(i)%khnxh = khnxh(i) !hyporheic heterotrophs ammonia preference
            hydrau%reach(i)%ahmax = ahmax(i) !hyporheic heterotrophs first-order carrying capacity
            hydrau%reach(i)%kgen = kgen(i) !generic constituent dissolution rate
            hydrau%reach(i)%vgen = vgen(i) !generic constituent settling velocity

            !gp 21-nov-06 initial conditions for each reach
            hydrau%reach(i)%te_ini = te_ini(i) !initial temperature
            hydrau%reach(i)%c01_ini = c01_ini(i) !initial cond
            hydrau%reach(i)%c02_ini = c02_ini(i) !initial iss
            hydrau%reach(i)%c03_ini = c03_ini(i) !initial do
            hydrau%reach(i)%c04_ini = c04_ini(i) !cbod slow
            hydrau%reach(i)%c05_ini = c05_ini(i) !cbod fast
            hydrau%reach(i)%c06_ini = c06_ini(i) !org n
            hydrau%reach(i)%c07_ini = c07_ini(i) !nh4 n
            hydrau%reach(i)%c08_ini = c08_ini(i) !no23 n
            hydrau%reach(i)%c09_ini = c09_ini(i) !org p
            hydrau%reach(i)%c10_ini = c10_ini(i) !srp
            hydrau%reach(i)%c11_ini = c11_ini(i) !phyto
            hydrau%reach(i)%c12_ini = c12_ini(i) !detritus
            hydrau%reach(i)%c13_ini = c13_ini(i) !pathogen
            hydrau%reach(i)%c14_ini = c14_ini(i) !generic
            hydrau%reach(i)%c15_ini = c15_ini(i) !alkalinity
            hydrau%reach(i)%ph_ini = ph_ini(i) !ph
            hydrau%reach(i)%c17_ini = c17_ini(i) !bottom algae gd/m^2
            hydrau%reach(i)%ninb_ini = ninb_ini(i) !bottom algae mgn/gd
            hydrau%reach(i)%nipb_ini = nipb_ini(i) !bottom algae mgp/gd

        end do
        hydrau%reach(0)%elev = elev2(0)
        hydrau%reach(nr+1)%elev1 = elev2(nr)
        hydrau%reach(0)%xl = hydrau%reach(1)%xl
        hydrau%reach(nr+1)%xl = hydrau%reach(nr)%xl


    end function hydraulics_


    subroutine allocatehydrauarrays(nr, hydrau)

        integer(i4b), intent(in) :: nr
        type(riverhydraulics_type), intent(inout) :: hydrau
        integer(i4b) status

        !if (.not. associated(hydrau%reach)) then
        if (nr > 0) then
            !allocate (reach(0:nr))
            allocate (hydrau%reach (0:nr+1), stat= status)
            if (status==1) stop 'Class_Hydraulics:AllocateHydrauArrays failed. Insufficient Memory!'
        else
            print *, 'ERROR:element number must be great than 0'
            stop 'Class_Hydraulics:AllocateHydrauArrays failed'
        end if
        !else
        !print *, 'warning: hydraulics array can only allocate once. allocation failed!'
        !end if

    end subroutine allocatehydrauarrays

    !gp 17-nov-04 subroutine makehydraulics(nr, hydrau, downstreamboundary, kai)
    subroutine makehydraulics(nr, hydrau, downstreamboundary, kai, geomethod) !gp 17-nov-04
        !gp new sub hydraulics to separate output and save flows in m**3/s and m**3/day
        implicit none
        integer(i4b) ,intent(in) ::nr
        type(riverhydraulics_type), intent(inout) :: hydrau

        !gp 17-nov-04
        character(len=30), intent(in) :: geomethod

        logical(lgt), intent(in) :: downstreamboundary

        !gp 10-feb-05
        !character(len=20) kai !reaeration model
        character(len=30) kai !reaeration model

        integer(i4b) i
        real(dp) travel, width1, width2, edif, en
        real(dp) eout, dxint, btop, ctsiv
        real(dp) ushear, fn, pwet, rh, hd

        !calculate hydraulics
        travel = 0

        do i = 0, nr !ne()

            if ((hydrau%reach(i)%hweir > 0).and.(hydrau%reach(i)%bweir > 0)) then
                if (hydrau%reach(i)%b <=0) then
                    stop 'Reaches with weirs must have a width entered!'
                end if
                hydrau%reach(i)%depth = hydrau%reach(i)%hweir + (hydrau%reach(i)%q &
                    / (1.83_dp * hydrau%reach(i)%bweir)**(2.0_dp/3.0_dp))
                hydrau%reach(i)%ac = hydrau%reach(i)%depth * hydrau%reach(i)%b
                hydrau%reach(i)%u = hydrau%reach(i)%q / hydrau%reach(i)%ac
                btop = hydrau%reach(i)%b
                pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
            elseif (hydrau%reach(i)%nm == 0) then
                !gp 17-nov-04
                !gp hydrau%reach(i)%u = hydrau%reach(i)%alp1 * hydrau%reach(i)%q ** hydrau%reach(i)%bet1
                !gp hydrau%reach(i)%depth = hydrau%reach(i)%alp2 * hydrau%reach(i)%q ** hydrau%reach(i)%bet2
                !gp hydrau%reach(i)%ac = hydrau%reach(i)%q / hydrau%reach(i)%u !cross section area
                !gp hydrau%reach(i)%b = hydrau%reach(i)%ac / hydrau%reach(i)%depth !width
                !gp btop = hydrau%reach(i)%b
                !gp pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
                if (geomethod == "Depth") then
                    hydrau%reach(i)%u = hydrau%reach(i)%alp1 * hydrau%reach(i)%q ** hydrau%reach(i)%bet1
                    hydrau%reach(i)%depth = hydrau%reach(i)%alp2 * hydrau%reach(i)%q ** hydrau%reach(i)%bet2
                    hydrau%reach(i)%ac = hydrau%reach(i)%q / hydrau%reach(i)%u !cross section area
                    hydrau%reach(i)%b = hydrau%reach(i)%ac / hydrau%reach(i)%depth !calc width from depth
                    btop = hydrau%reach(i)%b
                    pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
                else
                    hydrau%reach(i)%u = hydrau%reach(i)%alp1 * hydrau%reach(i)%q ** hydrau%reach(i)%bet1
                    hydrau%reach(i)%b = hydrau%reach(i)%alp2 * hydrau%reach(i)%q ** hydrau%reach(i)%bet2
                    hydrau%reach(i)%ac = hydrau%reach(i)%q / hydrau%reach(i)%u !cross section area
                    hydrau%reach(i)%depth = hydrau%reach(i)%ac / hydrau%reach(i)%b !calc depth from width
                    btop = hydrau%reach(i)%b
                    pwet = hydrau%reach(i)%b + 2 * hydrau%reach(i)%depth
                end if !gp 17-nov-04 end new block
            else !trapezoidal channel using manning's equation
                hydrau%reach(i)%depth = depthmanning(hydrau%reach(i)%q, hydrau%reach(i)%nm, hydrau%reach(i)%b, &
                    hydrau%reach(i)%ss1, hydrau%reach(i)%ss2, hydrau%reach(i)%s)
                hydrau%reach(i)%ac = (hydrau%reach(i)%b + 0.5_dp * (hydrau%reach(i)%ss1 + &
                    hydrau%reach(i)%ss2) * hydrau%reach(i)%depth) * hydrau%reach(i)%depth
                btop = hydrau%reach(i)%b + (hydrau%reach(i)%ss1 + hydrau%reach(i)%ss2) * &
                    hydrau%reach(i)%depth
                hydrau%reach(i)%u = hydrau%reach(i)%q / hydrau%reach(i)%ac
                pwet = hydrau%reach(i)%b+ hydrau%reach(i)%depth * sqrt(hydrau%reach(i)%ss1**2 +1) + &
                    hydrau%reach(i)%depth * sqrt(hydrau%reach(i)%ss2**2 +1)
            end if
            rh = hydrau%reach(i)%ac / pwet !hydraulic radius

            if (i > 0) then
                hydrau%reach(i)%asb = hydrau%reach(i)%frsed * hydrau%reach(i)%b * hydrau%reach(i)%xl
                hydrau%reach(i)%ast = btop * hydrau%reach(i)%xl
                hydrau%reach(i)%asd = hydrau%reach(i)%frsod * hydrau%reach(i)%b * hydrau%reach(i)%xl
                hydrau%reach(i)%vol = hydrau%reach(i)%ac * hydrau%reach(i)%xl
                !gp 03-nov-04 hyporheic pore volume in m^3
                hydrau%reach(i)%hypoporevol = hydrau%reach(i)%porosity * hydrau%reach(i)%ast * hydrau%reach(i)%hsedcm * 0.01_dp
            end if

            if (hydrau%reach(i)% ediff > 0) then
                edif = hydrau%reach(i)%ediff
            else
                if (hydrau%reach(i)%s >0) then
                    edif = 0.011_dp * hydrau%reach(i)%u ** 2 * hydrau%reach(i)%b ** 2 / hydrau%reach(i)%depth &
                        /sqrt(grav * hydrau%reach(i)%depth * hydrau%reach(i)%s)

                else if (hydrau%reach(i)%nm > 0) then
                    hydrau%reach(i)%s=(hydrau%reach(i)%nm *hydrau%reach(i)%u/rh**(2.0/3.0))**2
                    edif = 0.011_dp * hydrau%reach(i)%u ** 2 * hydrau%reach(i)%b ** 2 / hydrau%reach(i)%depth &
                        /sqrt(grav * hydrau%reach(i)%depth * hydrau%reach(i)%s)
                end if
            end if

            en = hydrau%reach(i)%u * hydrau%reach(i)%xl / 2.0_dp
            if (en <= edif) then
                edif = edif - en
                eout = edif
            else
                edif = 0
                eout = en
            end if

            dxint = (hydrau%reach(i)%xl + hydrau%reach(i+1)%xl) / 2.0_dp

            hydrau%reach(i)%ep = edif * hydrau%reach(i)%ac / dxint
            hydrau%reach(i)%epout = eout * hydrau%reach(i)%ac / dxint


            select case (kai)
              case ("Tsivoglou-Neal", "Thackston-Dawson", &
                  "USGS(pool-riffle)", "USGS(channel-control)")
                if (hydrau%reach(i)%s == 0) then
                    stop "kai requires a non-zero channel slope for segment "
                ELSE
                    Pwet = hydrau%reach(i)%b + hydrau%reach(i)%depth * &
                        SQRT(hydrau%reach(i)%SS1 ** 2 + 1.0) + hydrau%reach(i)%depth * &
                        SQRT(hydrau%reach(i)%SS2 ** 2 + 1.0)
                    Rh = hydrau%reach(i)%Ac / Pwet !Hydraulic radius
                    Ushear = SQRT(grav * Rh * hydrau%reach(i)%s)
                    Hd = hydrau%reach(i)%Ac / Btop !Hydraulic depth
                    FN = hydrau%reach(i)%U / SQRT(grav * Hd)
                END IF
            END SELECT


            IF (hydrau%reach(i)%kaaa > 0) THEN
                hydrau%reach(i)%kau = hydrau%reach(i)%kaaa
                hydrau%reach(i)%kaf = "Specified"
            elseif (kai == "Internal") then
                if (hydrau%reach(i)%depth < 0.61) then
                    hydrau%reach(i)%kau = 5.32_dp * hydrau%reach(i)%u ** 0.67_dp / hydrau%reach(i)%depth ** 1.85_dp
                    hydrau%reach(i)%kaf = "Owens"
                elseif (hydrau%reach(i)%depth > 3.45 * hydrau%reach(i)%u ** 2.5_dp) then
                    hydrau%reach(i)%kau = 3.93_dp * hydrau%reach(i)%u ** 0.5_dp / hydrau%reach(i)%depth ** 1.5_dp
                    hydrau%reach(i)%kaf = "O'Conn"
                else
                    hydrau%reach(i)%kau = 5.026_dp * hydrau%reach(i)%u / hydrau%reach(i)%depth ** 1.67_dp
                    hydrau%reach(i)%kaf = "Church"
                end if
            elseif (kai == "O'Connor-Dobbins") then
                hydrau%reach(i)%kau = 3.93_dp * hydrau%reach(i)%u ** 0.5_dp / hydrau%reach(i)%depth ** 1.5_dp
                hydrau%reach(i)%kaf = "O'Conn"
            elseif (kai == "Churchill") then
                hydrau%reach(i)%kau = 5.026_dp * hydrau%reach(i)%u / hydrau%reach(i)%depth ** 1.67_dp
                hydrau%reach(i)%kaf = "Church"
            elseif (kai == "Owens-Gibbs") then
                hydrau%reach(i)%kau = 5.32_dp * hydrau%reach(i)%u ** 0.67_dp / hydrau%reach(i)%depth ** 1.85_dp
                hydrau%reach(i)%kaf = "Owens"
            elseif (kai == "Tsivoglou-Neal") then
                ctsiv = 0.054_dp
                if (hydrau%reach(i)%q < (15.0_dp / 3.281_dp ** 3)) ctsiv = 0.11_dp
                hydrau%reach(i)%kau = 86400.0_dp * ctsiv * hydrau%reach(i)%s * hydrau%reach(i)%u * 3.281_dp
                hydrau%reach(i)%kaf = "Tsivoglou"
            elseif (kai == "Thackston-Dawson") then
                hydrau%reach(i)%kau = 86400.0_dp * 0.000025_dp * (1.0_dp + 9.0_dp * fn ** 0.25_dp) * ushear / &
                    hydrau%reach(i)%depth
                hydrau%reach(i)%kaf = "Thack-Dawson"
            elseif (kai == "USGS(pool-riffle)") then
                if (hydrau%reach(i)%q < 0.556) then
                    hydrau%reach(i)%kau = 517.0_dp * (hydrau%reach(i)%u * hydrau%reach(i)%s) ** 0.524_dp &
                        * hydrau%reach(i)%q ** (-0.242_dp)
                else
                    hydrau%reach(i)%kau = 596.0_dp * (hydrau%reach(i)%u * hydrau%reach(i)%s) ** 0.528_dp &
                        * hydrau%reach(i)%q ** (-0.136_dp)
                end if
                hydrau%reach(i)%kaf = "Pool-riffle"
            elseif (kai == "USGS(channel-control)") then
                if (hydrau%reach(i)%q < 0.556) then
                    hydrau%reach(i)%kau = 88.0_dp * (hydrau%reach(i)%u * hydrau%reach(i)%s) ** 0.313_dp &
                        * hydrau%reach(i)%depth ** (-0.353_dp)
                else
                    hydrau%reach(i)%kau = 142.0_dp * (hydrau%reach(i)%u * hydrau%reach(i)%s) ** 0.333_dp * &
                        hydrau%reach(i)%depth ** (-0.66_dp) * btop ** (-0.243_dp)
                end if
                hydrau%reach(i)%kaf = "Channel-control"
            end if

            travel = travel + hydrau%reach(i)%vol / hydrau%reach(i)%q / 86400.0_dp
            hydrau%reach(i)%trav = travel
        end do

        do i = 0, nr - 1
            if ((hydrau%reach(i)%hweir > 0.01_dp) .or. &
                (hydrau%reach(i)%elev2 > hydrau%reach(i+1)%elev1 + 0.01)) then
                hydrau%reach(i)%ep = 0
                hydrau%reach(i)%epout = 0
                hydrau%reach(i)%drop = hydrau%reach(i)%elev2 + hydrau%reach(i)%depth - &
                    hydrau%reach(i+1)%elev1 - hydrau%reach(i+1)%depth
            else
                hydrau%reach(i)%drop = 0
            end if
        end do

        !gp 17-jan-06
        if (.not.downstreamboundary) then
            !if (.not.db%downstreambound) then

            hydrau%reach(nr)%ep = 0
        end if

        !gp 17-jan-06
        !write(11,*) downstreamboundary

        !gp convert flows and bulk dispersion to cubic meters per day (cmd) units
        hydrau%reach(0)%qcmd = hydrau%reach(0)%q * 86400.0_dp
        hydrau%reach(0)%epcmd = hydrau%reach(0)%ep * 86400.0_dp
        do i = 1, nr
            hydrau%reach(i)%qptcmd = hydrau%reach(i)%qpt * 86400.0_dp
            hydrau%reach(i)%qptacmd = hydrau%reach(i)%qpta * 86400.0_dp
            hydrau%reach(i)%qcmd = hydrau%reach(i)%q * 86400.0_dp
            hydrau%reach(i)%epcmd = hydrau%reach(i)%ep * 86400.0_dp
            !gp 15-nov-04 hyporheic exchange flow in m^3/day
            hydrau%reach(i)%ehyporheiccmd = hydrau%reach(i)%hypoexchfrac * hydrau%reach(i)%qcmd !'units of (m3/d)

            !gp 17-jan-06
            !write(11,*) i, hydrau%reach(i)%ep

        end do

    end subroutine makehydraulics

    pure function depthmanning(qq, nmm, bb, sss1, sss2, ss) result(dep)
        implicit none
        real(dp),intent(in) :: qq, nmm, bb, sss1, sss2, ss
        real(dp) hold, ea
        real(dp) dep

        do
            hold = dep
            dep = (qq * nmm) ** 0.6_dp * (bb + dep * &
                sqrt(sss1 * sss1 + 1.0_dp) + dep * &
                sqrt(sss2 * sss2 + 1.0_dp)) ** 0.4_dp/ &
                (bb + 0.5_dp * (sss1 + sss2) * dep) / ss ** 0.3_dp
            ea = abs((dep - hold) / dep) * 100.0_dp
            if (ea < 0.001) exit
        end do

    end function


end module class_hydraulics
