module m_readfile
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use nrtype, only: lgt, hrsday, nv
    use m_system_params, only: system_params_t
    use m_rivertopo, only: t_rivertopo
    use m_hydraulics, only: riverhydraulics_type, hydraulics_, makehydraulics
    use m_meteorology, only: meteorology_t
    use class_lightheat, only: Light_
    use class_sourcein, only: sourcein_
    use m_rates, only: rates_t
    use m_water_quality, only:water_quality_t
    use m_upstream_boundary, only: upstream_boundary_t
    use m_downstream_boundary, only: downstream_boundary_t
    use class_solarcalc, only: solar_type, sitesolar_, sunrisesunset
    implicit none
    private
    public :: readinputfile
contains
    subroutine readinputfile(system, hydrau, sitemeteo, hw, db, stochrate, topo, sitesolar)
        type(system_params_t), intent(out) :: system
        type(riverhydraulics_type), intent(out) :: hydrau
        type(meteorology_t), intent(out) :: sitemeteo
        type(upstream_boundary_t) hw
        type(downstream_boundary_t) db
        type(rates_t) stochrate !stoch, reaction, temperature rate
        type(t_rivertopo) topo !river topology
        type(solar_type) :: sitesolar !solar radiation

        !gp 16-jul-08
        !gp 21-nov-06
        !gp 07-feb-06
        !gp 28-oct-04 integer(i32) i, j, k, status(32)
        !integer(i32) i, j, k, status(39) !gp 11-jan-05
        !integer(i32) i, j, k, status(92) !gp 11-jan-05
        !integer(i32) i, j, k, status(112) !gp 21-nov-06
        !integer(i32) i, j, k, status(113) !gp 03-apr-08
        integer(i32) i, j, k, status(114) !gp 16-jul-08

        integer(i32) nrch, nelem !number of reach, num of element
        integer(i32) nptin, ndiffin !num of point load/src, num of non-point
        !system parameters

        !gp 23-nov-09
        !character(len=30) basinname, filename, path, title, timezone
        character(len=30) basinname, filename, path, title
        real(r64) timezone

        character(len=30) :: imeth !integration method
        character(len=30) :: imethph !ph method
        character(len=30) :: simhyporheicwq !gp 28-oct-04 yes or no to simulate hyporheic water quality and flux
        character(len=30) :: showdielresults !gp 03-feb-05 yes or no (only used in excel vba)
        character(len=30) :: statevariables !gp 03-feb-05 all or only temperature

        !gp 11-jan-06
        character(len=30) :: calcsedflux !gp 11-jan-06 yes or no

        !gp 11-jan-06
        character(len=30) :: simalk !gp 26-oct-07 yes or no

        !gp 24-jun-09
        character(len=30) :: writedynamic !yes or no

        real(r64) year, month, day
        real(r64) dtuser, tf

        !reach data
        character(len=30) :: geomethod !gp 17-nov-04 depth or width for coeff/exponents in col t and u of 'reach' sheet
        real(r64), allocatable :: xrdn(:), elev1(:), elev2(:), latd(:), latm(:)


        real(r64), allocatable :: lats(:), lond(:), lonm(:), lons(:), q(:), bb(:), &
            ss1(:), ss2(:), s(:), nm(:) , alp1(:), bet1(:), &
            alp2(:), bet2(:), ediff(:), frsed(:), &
            frsod(:), sodspec(:), jch4spec(:), jnh4spec(:), &
            jsrpspec(:), hweir(:), bweir(:), &
            sedthermcond(:), sedthermdiff(:), hsedcm(:), &
            hypoexchfrac(:), porosity(:), rhocpsed(:), skop(:), &
            kaaa(:), vss_rch(:), khc_rch(:), &
            kdcs_rch(:), kdc_rch(:), khn_rch(:), &
            von_rch(:), kn_rch(:) , ki_rch(:) , &
            vdi_rch(:), khp_rch(:), vop_rch(:), &
            vip_rch(:), kga_rch(:), krea_rch(:), &
            kdea_rch(:), ksn_rch(:), ksp_rch(:), &
            isat_rch(:), khnx_rch(:), va_rch(:), &
            botalg0(:), kgaf_rch(:), abmax_rch(:), &
            krea1f_rch(:), krea2f_rch(:), kexaf_rch(:), kdeaf_rch(:), &
            ksnf_rch(:), kspf_rch(:), isatf_rch(:), &
            khnxf_rch(:), ninbmin_rch(:), nipbmin_rch(:), &
            ninbupmax_rch(:), nipbupmax_rch(:), kqn_rch(:), &
            kqp_rch(:), nupwcfrac_rch(:), pupwcfrac_rch(:), &
            kdt_rch(:), vdt_rch(:), kpath_rch(:), &
            vpath_rch(:), apath_rch(:), kgah_rch(:), &
            ksch_rch(:), kinhch_rch(:), kreah_rch(:), &
            kdeah_rch(:), ksnh_rch(:), ksph_rch(:), &
            khnxh_rch(:), ahmax_rch(:), &
            kgen_rch(:), vgen_rch(:), &
            te_ini(:), c01_ini(:), c02_ini(:), c03_ini(:), &
            c04_ini(:), c05_ini(:), c06_ini(:), &
            c07_ini(:), c08_ini(:), c09_ini(:), &
            c10_ini(:), c11_ini(:), c12_ini(:), &
            c13_ini(:), c14_ini(:), c15_ini(:), &
            ph_ini(:), c17_ini(:), ninb_ini(:), nipb_ini(:)

        character(len=30), allocatable :: rlab1(:), rlab2(:), rname(:)

        real(r64) par, kep, kela, kenla, kess, kepom, kemac, nfacbras, atcryanstolz, kbrut, kcl1, kcl2

        character(len=30) solarmethod, longatmethod, fuwmethod

        !point load/source
        character(len=30), allocatable :: ptname(:)
        real(r64), allocatable :: xptt(:), qptta(:), qptt(:), tepttmean(:), tepttamp(:), tepttmaxtime(:)
        real(r64), allocatable :: cpttmean(:, :), cpttamp(:, :), cpttmaxtime(:, :)
        real(r64), allocatable :: phpttmean(:), phpttamp(:), phpttmaxtime(:)
        !diffusion load/source
        real(r64), allocatable :: xdup(:), xddn(:), qdifa(:), qdif(:), phind(:), tedif(:), cdif(:,:)
        character(len=30), allocatable :: diffname(:)
        !rates
        real(r64) vss, mgc, mgn, mgp, mgd, mga
        real(r64) tka, roc, ron
        real(r64) ksocf, ksona, ksodn, ksop, ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von
        real(r64) kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        real(r64) kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, isat

        !gp 03-apr-08
        !real(r64) khnx, va, kgaf, tkgaf, kreaf, tkreaf, kexaf, tkexaf, kdeaf, abmax
        real(r64) khnx, va, kgaf, tkgaf, krea1f, krea2f, tkreaf, kexaf, tkexaf, kdeaf, abmax

        real(r64) tkdeaf, ksnf, kspf, kscf, isatf, khnxf, kdt, tkdt, vdt

        !gp 10-jan-05 real(r64) ninbmin, nipbmin, ninbupmax, nipbupmax, kqn, kqp
        real(r64) ninbmin, nipbmin, ninbupmax, nipbupmax, kqn, kqp, kqnratio, kqpratio

        !gp 26-jan-06
        real(r64) nupwcfrac, pupwcfrac

        real(r64) kpath, tkpath, vpath

        !gp 30-nov-04
        real(r64) apath, kgen, tkgen, vgen !gp 30-nov-04 new params for pathogen and generic constituent

        !gp 08-dec-04
        character(len=30) usegenericascod

        real(r64) pco2
        !org c, n inhition model, denitrification enhance model, algae light model, perihyte light model
        character(len=30) xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, typef
        character(len=30) kawindmethod !reaeartion wind effect
        !data
        integer(i32) kk, nteda, nhydda, nwqd(3), nrp, nwqdiur
        real(r64) junk !junk data
        character(len=30) xxx, shtnam(3) !junk string,

        !gp 23-nov-09
        !integer(i32) dlstime
        real(r64) dlstime

        logical(lgt) downstreamboundary

        !gp 03-nov-04 hyporheic biofilm rates
        character(len=30) typeh, xdum8
        real(r64) kgah, tkgah, ksch, kinhch

        !gp 15-nov-04 hco3- use by phytoplankton and bottom algae
        character(len=30) hco3use, hco3usef

        !gp 15-nov-04 level 2 hyporheic biofilm rates
        real(r64) kreah, tkreah, kdeah, tkdeah, ksnh, ksph, khnxh, ahmax

        type(water_quality_t), dimension(:), allocatable :: hwin !headwaters
        type(water_quality_t), dimension(:), pointer :: dbin !downstream bondary

        !meteorology data

        !gp 16-jul-08
        !real(r64), allocatable :: shadehh(:,:), tahh(:,:), tdhh(:,:), uwhh(:,:), cchh(:,:)
        real(r64), allocatable :: shadehh(:,:), tahh(:,:), tdhh(:,:), uwhh(:,:), cchh(:,:), solarhh(:,:)

        !/* read in input file */
        read(8,*) basinname, filename, path, title
        read (8,*) month, day, year
        !gp 28-oct-04 read(8,*) timezone, pco2, dtuser, tf, imeth, imethph

        !gp 24-jun-09
        !read(8,*) timezone, pco2, dtuser, tf, imeth, imethph, simhyporheicwq, showdielresults, statevariables, calcsedflux, simalk !gp 26-oct-07
        read(8,*) timezone, pco2, dtuser, tf, imeth, imethph, simhyporheicwq, showdielresults, &
            statevariables, calcsedflux, simalk, writedynamic

        !read(8,*) imeth, imethph

        system = system_params_t(basinname, filename, path, title, year, month, day, &
            timezone, dtuser, tf, imeth, imethph, simhyporheicwq, &
            showdielresults, statevariables, calcsedflux, simalk, writedynamic)

        !reach data
        !gp 17-nov-04 read(8,*) nrch !total reach number
        read(8,*) nrch, geomethod !gp 17,nov-04 number of reaches and depth or width for col t and u

        allocate (xrdn(0:nrch+1), stat=status(1)); allocate (elev1(0:nrch+1), stat=status(2))
        allocate (elev2(0:nrch+1), stat=status(3)); allocate (latd(0:nrch+1), stat=status(4))
        allocate (latm(0:nrch+1), stat=status(5)); allocate (lats(0:nrch+1), stat=status(6))
        allocate (lond(0:nrch+1), stat=status(7)); allocate (lonm(0:nrch+1), stat=status(8))
        allocate (lons(0:nrch+1), stat=status(9)); allocate (q(0:nrch+1), stat=status(10))
        allocate (bb(0:nrch+1), stat=status(11)); allocate (ss1(0:nrch+1), stat=status(12))
        allocate (ss2(0:nrch+1), stat=status(13)); allocate (s(0:nrch+1), stat=status(14))
        allocate (nm(0:nrch+1), stat=status(15)); allocate (alp1(0:nrch+1), stat=status(16))
        allocate (bet1(0:nrch+1), stat=status(17)); allocate (alp2(0:nrch+1), stat=status(18))
        allocate (bet2(0:nrch+1), stat=status(19)); allocate (ediff(0:nrch+1), stat=status(20))
        allocate (kaaa(0:nrch+1), stat=status(21)); allocate (frsed(0:nrch+1), stat=status(22))
        allocate (frsod(0:nrch+1), stat=status(23));allocate (sodspec(0:nrch+1), stat=status(24))
        allocate (jch4spec(0:nrch+1), stat=status(25)); allocate (jnh4spec(0:nrch+1), stat=status(26))
        allocate (jsrpspec(0:nrch+1), stat=status(27)); allocate (hweir(0:nrch+1), stat=status(28))
        allocate (rlab1(0:nrch+1), stat=status(29)); allocate (rlab2(0:nrch+1), stat=status(30))
        allocate (rname(0:nrch+1), stat=status(31));allocate (bweir(0:nrch+1), stat=status(32))
        allocate (sedthermcond(0:nrch+1), stat=status(33));allocate (sedthermdiff(0:nrch+1), stat=status(34)) !gp 20-oct-04
        allocate (hsedcm(0:nrch+1), stat=status(35));allocate (hypoexchfrac(0:nrch+1), stat=status(36)) !gp 20-oct-04
        allocate (porosity(0:nrch+1), stat=status(37));allocate (rhocpsed(0:nrch+1), stat=status(38)) !gp 20-oct-04

        !gp 16-jul-08
        allocate (skop(0:nrch+1), stat=status(114)) !gp 11-jan-05

        allocate (botalg0(0:nrch+1), stat=status(39)) !gp 11-jan-05

        !gp 07-feb-07
        !allocate (kaaa(0:nrch+1), stat=status(39))
        allocate (vss_rch(0:nrch+1), stat=status(40)) !inorganic suspended solids settling vol
        allocate (khc_rch(0:nrch+1), stat=status(41)) !slow cbod hydrolysis rate
        allocate (kdcs_rch(0:nrch+1), stat=status(42)) !slow cbod oxidation rate
        allocate (kdc_rch(0:nrch+1), stat=status(43)) !fast cbod oxidation rate
        allocate (khn_rch(0:nrch+1), stat=status(44)) !organic n hydrolysis rate
        allocate (von_rch(0:nrch+1), stat=status(45)) !organic n settling velocity
        allocate (kn_rch(0:nrch+1), stat=status(46)) !ammonium nitrification rate
        allocate (ki_rch(0:nrch+1), stat=status(47)) !nitrate denitrification
        allocate (vdi_rch(0:nrch+1), stat=status(48)) !nitrate sed denitrification transfer coeff
        allocate (khp_rch(0:nrch+1), stat=status(49)) !organic p hydrolysis
        allocate (vop_rch(0:nrch+1), stat=status(50)) !organic p settling velocity
        allocate (vip_rch(0:nrch+1), stat=status(51)) !inorganic p settling velocity
        allocate (kga_rch(0:nrch+1), stat=status(52)) !phytoplankton max growth rate
        allocate (krea_rch(0:nrch+1), stat=status(53)) !phytoplankton respiration rate
        allocate (kdea_rch(0:nrch+1), stat=status(54)) !phytoplankton death rate
        allocate (ksn_rch(0:nrch+1), stat=status(55)) !phytoplankton n half-sat
        allocate (ksp_rch(0:nrch+1), stat=status(56)) !phytoplankton p half-sat
        allocate (isat_rch(0:nrch+1), stat=status(57)) !phytoplankton light sat
        allocate (khnx_rch(0:nrch+1), stat=status(58)) !phytoplankton ammonia preference
        allocate (va_rch(0:nrch+1), stat=status(59)) !phytoplankton settling velocity
        !allocate (botalg0(0:nrch+1), stat=status(39)) !bottom plant initial bionass
        allocate (kgaf_rch(0:nrch+1), stat=status(60)) !bottom plant max growth rate
        allocate (abmax_rch(0:nrch+1), stat=status(61)) !bottom plant first-order carrying capacity

        !gp 03-apr-08
        !allocate (kreaf_rch(0:nrch+1), stat=status(62)) !bottom plant respiration rate
        allocate (krea1f_rch(0:nrch+1), stat=status(62)) !bottom plant basal respiration rate
        allocate (krea2f_rch(0:nrch+1), stat=status(113)) !bottom plant photo respiration rate

        allocate (kexaf_rch(0:nrch+1), stat=status(63)) !bottom plant excretion rate
        allocate (kdeaf_rch(0:nrch+1), stat=status(64)) !bottom plant death rate
        allocate (ksnf_rch(0:nrch+1), stat=status(65)) !bottom plant external n half-sat
        allocate (kspf_rch(0:nrch+1), stat=status(66)) !bottom plant external p half-sat
        allocate (isatf_rch(0:nrch+1), stat=status(67)) !bottom plant light sat
        allocate (khnxf_rch(0:nrch+1), stat=status(68)) !bottom plant ammonia preference
        allocate (ninbmin_rch(0:nrch+1), stat=status(69)) !bottom plant subistence quota for n
        allocate (nipbmin_rch(0:nrch+1), stat=status(70)) !bottom plant subistence quota for p
        allocate (ninbupmax_rch(0:nrch+1), stat=status(71)) !bottom plant max uptake rate for n
        allocate (nipbupmax_rch(0:nrch+1), stat=status(72)) !bottom plant max uptake rate for p
        allocate (kqn_rch(0:nrch+1), stat=status(73)) !bottom plant internal n half-sat
        allocate (kqp_rch(0:nrch+1), stat=status(74)) !bottom plant internal p half-sat
        allocate (nupwcfrac_rch(0:nrch+1), stat=status(75)) !bottom plant n uptake fraction from water column
        allocate (pupwcfrac_rch(0:nrch+1), stat=status(76)) !bottom plant p uptake fraction from water column
        allocate (kdt_rch(0:nrch+1), stat=status(77)) !pom dissolution rate
        allocate (vdt_rch(0:nrch+1), stat=status(78)) !pom settling velocity
        allocate (kpath_rch(0:nrch+1), stat=status(79)) !pathogen dieoff rate
        allocate (vpath_rch(0:nrch+1), stat=status(80)) !pathogen settling velocity
        allocate (apath_rch(0:nrch+1), stat=status(81)) !pathogen light alpha
        allocate (kgah_rch(0:nrch+1), stat=status(82)) !hyporheic heterotrophs max growth rate
        allocate (ksch_rch(0:nrch+1), stat=status(83)) !hyporheic heterotrophs cbod half-sat
        allocate (kinhch_rch(0:nrch+1), stat=status(84)) !hyporheic heterotrophs o2 inhibition
        allocate (kreah_rch(0:nrch+1), stat=status(85)) !hyporheic heterotrophs respiration rate
        allocate (kdeah_rch(0:nrch+1), stat=status(86)) !hyporheic heterotrophs death rate
        allocate (ksnh_rch(0:nrch+1), stat=status(87)) !hyporheic heterotrophs n half-sat
        allocate (ksph_rch(0:nrch+1), stat=status(88)) !hyporheic heterotrophs p half-sat
        allocate (khnxh_rch(0:nrch+1), stat=status(89)) !hyporheic heterotrophs ammonia preference
        allocate (ahmax_rch(0:nrch+1), stat=status(90)) !hyporheic heterotrophs first-order carrying capacity
        allocate (kgen_rch(0:nrch+1), stat=status(91)) !generic constituent dissolution rate
        allocate (vgen_rch(0:nrch+1), stat=status(92)) !generic constituent settling velocity

        !gp 21-nov-06 initial conditions of state variables
        allocate (te_ini(0:nrch+1), stat=status(93))
        allocate (c01_ini(0:nrch+1), stat=status(94))
        allocate (c02_ini(0:nrch+1), stat=status(95))
        allocate (c03_ini(0:nrch+1), stat=status(96))
        allocate (c04_ini(0:nrch+1), stat=status(97))
        allocate (c05_ini(0:nrch+1), stat=status(98))
        allocate (c06_ini(0:nrch+1), stat=status(99))
        allocate (c07_ini(0:nrch+1), stat=status(100))
        allocate (c08_ini(0:nrch+1), stat=status(101))
        allocate (c09_ini(0:nrch+1), stat=status(102))
        allocate (c10_ini(0:nrch+1), stat=status(103))
        allocate (c11_ini(0:nrch+1), stat=status(104))
        allocate (c12_ini(0:nrch+1), stat=status(105))
        allocate (c13_ini(0:nrch+1), stat=status(106))
        allocate (c14_ini(0:nrch+1), stat=status(107))
        allocate (c15_ini(0:nrch+1), stat=status(108))
        allocate (ph_ini(0:nrch+1), stat=status(109))
        allocate (c17_ini(0:nrch+1), stat=status(110))
        allocate (ninb_ini(0:nrch+1), stat=status(111))
        allocate (nipb_ini(0:nrch+1), stat=status(112))
        !gp 03-apr-08 status(113) is used above for krea2f
        !gp 16-jul-08 status(114) is used above for skop

        !classriver module

        do i = 0, nrch + 1
            read(8,*) rlab1(i), rlab2(i), rname(i), xrdn(i), elev1(i), elev2(i), latd(i), latm(i), &
                lats(i), lond(i), lonm(i), lons(i), q(i), bb(i), ss1(i), ss2(i), s(i), nm(i), &
                alp1(i), bet1(i), alp2(i), bet2(i), ediff(i), frsed(i), frsod(i), &
                sodspec(i), jch4spec(i), jnh4spec(i), jsrpspec(i), hweir(i), bweir(i), &
                sedthermcond(i), sedthermdiff(i), hsedcm(i), &
                hypoexchfrac(i), porosity(i), rhocpsed(i), skop(i)
        end do

        do i=1, nrch

            read(8,*) kaaa(i), vss_rch(i), khc_rch(i), &
                kdcs_rch(i), kdc_rch(i), khn_rch(i), &
                von_rch(i), kn_rch(i) , ki_rch(i) , &
                vdi_rch(i), khp_rch(i), vop_rch(i), &
                vip_rch(i), kga_rch(i), krea_rch(i), &
                kdea_rch(i), ksn_rch(i), ksp_rch(i), &
                isat_rch(i), khnx_rch(i), va_rch(i), &
                kgaf_rch(i), abmax_rch(i), &
                krea1f_rch(i), krea2f_rch(i), kexaf_rch(i), kdeaf_rch(i), &
                ksnf_rch(i), kspf_rch(i), isatf_rch(i), &
                khnxf_rch(i), ninbmin_rch(i), nipbmin_rch(i), &
                ninbupmax_rch(i), nipbupmax_rch(i), kqn_rch(i), &
                kqp_rch(i), nupwcfrac_rch(i), pupwcfrac_rch(i), &
                kdt_rch(i), vdt_rch(i), kpath_rch(i), &
                vpath_rch(i), apath_rch(i), kgah_rch(i), &
                ksch_rch(i), kinhch_rch(i), kreah_rch(i), &
                kdeah_rch(i), ksnh_rch(i), ksph_rch(i), &
                khnxh_rch(i), ahmax_rch(i), &
                kgen_rch(i), vgen_rch(i)

            ! write(10,'(i5, 4f10.3)') i, kaaa(i), vss_rch(i), kgen_rch(i), vgen_rch(i)
        end do
        do i=1, nrch
            read(8,*) te_ini(i), ph_ini(i), c17_ini(i), ninb_ini(i), nipb_ini(i)
            ! write(10,'(i5, 4f10.3)') i, te_ini(i), c17_ini(i), ninb_ini(i), nipb_ini(i)
        end do
        do i=1, nrch
            read(8,*) c01_ini(i), c02_ini(i), &
                c03_ini(i), c04_ini(i), c05_ini(i), &
                c06_ini(i), c07_ini(i), c08_ini(i), &
                c09_ini(i), c10_ini(i), c11_ini(i), &
                c12_ini(i), c13_ini(i), c14_ini(i), &
                c15_ini(i)
            ! write(10,'(i5, 4f10.3)') i, c01_ini(i), c02_ini(i), c14_ini(i), c15_ini(i)
        end do
        !write(10,*) 'done thru classreadfile reading new init cond'

        !gp 17-nov-04 topo= rivertopo_(nrch, rlab2, rname, xrdn)
        topo= t_rivertopo(nrch, rlab2, rname, xrdn, geomethod) !gp 17-nov-04

        hydrau= hydraulics_(nrch, xrdn, elev1, elev2, latd, latm, lats, &
            lond, lonm, lons, q, bb, ss1, ss2, s, nm, alp1, bet1, alp2, bet2, &
            ediff, frsed, frsod, sodspec, jch4spec, jnh4spec, jsrpspec, &
            hweir, bweir, &
            sedthermcond, sedthermdiff, hsedcm, &
            hypoexchfrac, porosity, rhocpsed, skop, &
            kaaa, vss_rch, khc_rch, &
            kdcs_rch, kdc_rch, khn_rch, &
            von_rch, kn_rch, ki_rch, &
            vdi_rch, khp_rch, vop_rch, &
            vip_rch, kga_rch, krea_rch, &
            kdea_rch, ksn_rch, ksp_rch, &
            isat_rch, khnx_rch, va_rch, &
            kgaf_rch, abmax_rch, &
            krea1f_rch, krea2f_rch, kexaf_rch, kdeaf_rch, &
            ksnf_rch, kspf_rch, isatf_rch, &
            khnxf_rch, ninbmin_rch, nipbmin_rch, &
            ninbupmax_rch, nipbupmax_rch, kqn_rch, &
            kqp_rch, nupwcfrac_rch, pupwcfrac_rch, &
            kdt_rch, vdt_rch, kpath_rch, &
            vpath_rch, apath_rch, kgah_rch, &
            ksch_rch, kinhch_rch, kreah_rch, &
            kdeah_rch, ksnh_rch, ksph_rch, &
            khnxh_rch, ahmax_rch, &
            kgen_rch, vgen_rch, &
            te_ini, c01_ini, c02_ini, c03_ini, &
            c04_ini, c05_ini, c06_ini, &
            c07_ini, c08_ini, c09_ini, &
            c10_ini, c11_ini, c12_ini, &
            c13_ini, c14_ini, c15_ini, &
            ph_ini, c17_ini, ninb_ini, nipb_ini)

        !gp 21-nov-06
        !write(10,*) 'done thru classreadfile hydrau'

        deallocate (botalg0, stat=status(39)) !gp 11-jan-05
        deallocate (rhocpsed, stat=status(38)) !gp 20-oct-04

        !gp 16-jul-08
        deallocate (skop, stat=status(114)) !gp 20-oct-04

        deallocate (porosity, stat=status(37)) !gp 20-oct-04
        deallocate (hypoexchfrac, stat=status(36)) !gp 20-oct-04
        deallocate (hsedcm, stat=status(35)) !gp 20-oct-04
        deallocate (sedthermdiff, stat=status(34)) !gp 20-oct-04
        deallocate (sedthermcond, stat=status(33)) !gp 20-oct-04
        deallocate (bweir, stat=status(32)) !gp 20-oct-04 (missing from original code)
        deallocate (rname, stat=status(31))
        deallocate (rlab2, stat=status(30))
        deallocate (rlab1, stat=status(29))
        deallocate (hweir, stat=status(28))
        deallocate (jsrpspec, stat=status(27))
        deallocate (jnh4spec, stat=status(26))
        deallocate (jch4spec, stat=status(25))
        deallocate (sodspec, stat=status(24))
        deallocate (frsod, stat=status(23))
        deallocate (frsed, stat=status(22))
        deallocate (kaaa, stat=status(21))
        deallocate (ediff, stat=status(20))
        deallocate (bet2, stat=status(19))
        deallocate (alp2, stat=status(18))
        deallocate (bet1, stat=status(17))
        deallocate (alp1, stat=status(16))
        deallocate (nm, stat=status(15))
        deallocate (s, stat=status(14))
        deallocate (ss2, stat=status(13))
        deallocate (ss1, stat=status(12))
        deallocate (bb, stat=status(11))
        deallocate (q, stat=status(10))
        deallocate (lons, stat=status(9))
        deallocate (lonm, stat=status(8))
        deallocate (lond, stat=status(7))
        deallocate (lats, stat=status(6))
        deallocate (latm, stat=status(5))
        deallocate (latd, stat=status(4))
        deallocate (elev2, stat=status(3))
        deallocate (elev1, stat=status(2))
        deallocate (xrdn, stat=status(1))

        !gp 07-feb-07
        !deallocate (kaaa, stat=status(39))
        deallocate (vss_rch, stat=status(40)) !inorganic suspended solids settling vol
        deallocate (khc_rch, stat=status(41)) !slow cbod hydrolysis rate
        deallocate (kdcs_rch, stat=status(42)) !slow cbod oxidation rate
        deallocate (kdc_rch, stat=status(43)) !fast cbod oxidation rate
        deallocate (khn_rch, stat=status(44)) !organic n hydrolysis rate
        deallocate (von_rch, stat=status(45)) !organic n settling velocity
        deallocate (kn_rch, stat=status(46)) !ammonium nitrification rate
        deallocate (ki_rch, stat=status(47)) !nitrate denitrification
        deallocate (vdi_rch, stat=status(48)) !nitrate sed denitrification transfer coeff
        deallocate (khp_rch, stat=status(49)) !organic p hydrolysis
        deallocate (vop_rch, stat=status(50)) !organic p settling velocity
        deallocate (vip_rch, stat=status(51)) !inorganic p settling velocity
        deallocate (kga_rch, stat=status(52)) !phytoplankton max growth rate
        deallocate (krea_rch, stat=status(53)) !phytoplankton respiration rate
        deallocate (kdea_rch, stat=status(54)) !phytoplankton death rate
        deallocate (ksn_rch, stat=status(55)) !phytoplankton n half-sat
        deallocate (ksp_rch, stat=status(56)) !phytoplankton p half-sat
        deallocate (isat_rch, stat=status(57)) !phytoplankton light sat
        deallocate (khnx_rch, stat=status(58)) !phytoplankton ammonia preference
        deallocate (va_rch, stat=status(59)) !phytoplankton settling velocity
        !deallocate (botalg0, stat=status(39)) !bottom plant initial bionass
        deallocate (kgaf_rch, stat=status(60)) !bottom plant max growth rate
        deallocate (abmax_rch, stat=status(61)) !bottom plant first-order carrying capacity

        !gp 03-apr-08
        !deallocate (kreaf_rch, stat=status(62)) !bottom plant respiration rate
        deallocate (krea1f_rch, stat=status(62)) !bottom plant basal respiration rate
        deallocate (krea2f_rch, stat=status(113)) !bottom plant photo respiration rate

        deallocate (kexaf_rch, stat=status(63)) !bottom plant excretion rate
        deallocate (kdeaf_rch, stat=status(64)) !bottom plant death rate
        deallocate (ksnf_rch, stat=status(65)) !bottom plant external n half-sat
        deallocate (kspf_rch, stat=status(66)) !bottom plant external p half-sat
        deallocate (isatf_rch, stat=status(67)) !bottom plant light sat
        deallocate (khnxf_rch, stat=status(68)) !bottom plant ammonia preference
        deallocate (ninbmin_rch, stat=status(69)) !bottom plant subistence quota for n
        deallocate (nipbmin_rch, stat=status(70)) !bottom plant subistence quota for p
        deallocate (ninbupmax_rch, stat=status(71)) !bottom plant max uptake rate for n
        deallocate (nipbupmax_rch, stat=status(72)) !bottom plant max uptake rate for p
        deallocate (kqn_rch, stat=status(73)) !bottom plant internal n half-sat
        deallocate (kqp_rch, stat=status(74)) !bottom plant internal p half-sat
        deallocate (nupwcfrac_rch, stat=status(75)) !bottom plant n uptake fraction from water column
        deallocate (pupwcfrac_rch, stat=status(76)) !bottom plant p uptake fraction from water column
        deallocate (kdt_rch, stat=status(77)) !pom dissolution rate
        deallocate (vdt_rch, stat=status(78)) !pom settling velocity
        deallocate (kpath_rch, stat=status(79)) !pathogen dieoff rate
        deallocate (vpath_rch, stat=status(80)) !pathogen settling velocity
        deallocate (apath_rch, stat=status(81)) !pathogen light alpha
        deallocate (kgah_rch, stat=status(82)) !hyporheic heterotrophs max growth rate
        deallocate (ksch_rch, stat=status(83)) !hyporheic heterotrophs cbod half-sat
        deallocate (kinhch_rch, stat=status(84)) !hyporheic heterotrophs o2 inhibition
        deallocate (kreah_rch, stat=status(85)) !hyporheic heterotrophs respiration rate
        deallocate (kdeah_rch, stat=status(86)) !hyporheic heterotrophs death rate
        deallocate (ksnh_rch, stat=status(87)) !hyporheic heterotrophs n half-sat
        deallocate (ksph_rch, stat=status(88)) !hyporheic heterotrophs p half-sat
        deallocate (khnxh_rch, stat=status(89)) !hyporheic heterotrophs ammonia preference
        deallocate (ahmax_rch, stat=status(90)) !hyporheic heterotrophs first-order carrying capacity
        deallocate (kgen_rch, stat=status(91)) !generic constituent dissolution rate
        deallocate (vgen_rch, stat=status(92)) !generic constituent settling velocity

        !gp 21-nov-06 initial conditions of state variables
        deallocate (te_ini, stat=status(93))
        deallocate (c01_ini, stat=status(94))
        deallocate (c02_ini, stat=status(95))
        deallocate (c03_ini, stat=status(96))
        deallocate (c04_ini, stat=status(97))
        deallocate (c05_ini, stat=status(98))
        deallocate (c06_ini, stat=status(99))
        deallocate (c07_ini, stat=status(100))
        deallocate (c08_ini, stat=status(101))
        deallocate (c09_ini, stat=status(102))
        deallocate (c10_ini, stat=status(103))
        deallocate (c11_ini, stat=status(104))
        deallocate (c12_ini, stat=status(105))
        deallocate (c13_ini, stat=status(106))
        deallocate (c14_ini, stat=status(107))
        deallocate (c15_ini, stat=status(108))
        deallocate (ph_ini, stat=status(109))
        deallocate (c17_ini, stat=status(110))
        deallocate (ninb_ini, stat=status(111))
        deallocate (nipb_ini, stat=status(112))
        !gp 03-apr-08 status(113) is used abouve for krea2f

        !gp 21-nov-06
        !write(10,*) 'done thru classreadfile deallocate'

        !light data

        !gp 13-feb-06
        !read(8,*) par, kep, kela, kenla, kess, kepom
        read(8,*) par, kep, kela, kenla, kess, kepom, kemac

        !gp 24-jun-09
        !gp 16-jul-08
        !read(8,*) solarmethod, nfacbras, atcryanstolz, longatmethod, fuwmethod
        !read(8,*) solarmethod, nfacbras, atcryanstolz, longatmethod, kbrut, fuwmethod
        read(8,*) solarmethod, nfacbras, atcryanstolz, longatmethod, kbrut, fuwmethod, kcl1, kcl2

        !light and heat module

        !gp 13-feb-06
        !call light_(par, kep, kela, kenla, kess, kepom, longatmethod, fuwmethod)

        !gp 24-jun-09
        !gp 16-jul-08
        !call light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, fuwmethod)
        !call light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, kbrut, fuwmethod)
        call light_(par, kep, kela, kenla, kess, kepom, kemac, longatmethod, kbrut, fuwmethod, kcl1, kcl2)

        read(8,*) nptin

        if (nptin>=1) then !if any point load/source
            allocate (ptname(nptin), stat=status(1)); allocate (xptt(nptin), stat=status(2))
            allocate (qptta(nptin), stat=status(3)); allocate (qptt(nptin), stat=status(4))
            allocate (tepttmean(nptin), stat=status(5));allocate (tepttamp(nptin), stat=status(6))
            allocate (tepttmaxtime(nptin), stat=status(7))
            allocate (cpttmean(nptin, nv-2), stat=status(8)); allocate(cpttamp(nptin, nv-2), stat=status(9))
            allocate(cpttmaxtime(nptin, nv-2), stat=status(10)); allocate (phpttmean(nptin), stat=status(11))
            allocate (phpttamp(nptin), stat=status(12));allocate(phpttmaxtime(nptin), stat=status(13))
            !point sources

            do i=1 , 13

                if (status(i)==1) then
                    stop 'Class_ReadFile:ReadInputfile. dynamic memory allocation failed!'
                end if
            end do
            do i = 1, nptin
                read(8,*) ptname(i), xptt(i), qptta(i), qptt(i), tepttmean(i), tepttamp(i), tepttmaxtime(i)
                do k = 1, nv - 2
                    read(8,*) cpttmean(i, k), cpttamp(i, k), cpttmaxtime(i, k)
                end do
                read(8,*) phpttmean(i), phpttamp(i), phpttmaxtime(i)
            end do

        end if


        !diffuse sources
        read(8,*) ndiffin
        if (ndiffin>=1) then !if any nonpoint load/source
            allocate (xdup(ndiffin), stat=status(1)); allocate(xddn(ndiffin), stat=status(2))
            allocate(qdifa(ndiffin), stat=status(3)); allocate (qdif(ndiffin), stat=status(4))
            allocate(tedif(ndiffin), stat=status(5)); allocate(diffname(ndiffin), stat=status(6))
            allocate(phind(ndiffin), stat=status(7)); allocate(cdif(ndiffin, nv-2), stat=status(8))
            do i=1, 8
                if (status(i)==1) then
                    stop 'Class_ReadFile:ReadInputfile. dynamic memory allocation failed!'
                end if
            end do
            do i = 1, ndiffin
                read(8,*) diffname(i), xdup(i), xddn(i), qdifa(i), qdif(i), tedif(i)
                do k = 1, nv - 2
                    read(8,*) cdif(i, k)
                end do
                read(8,*) phind(i)
            end do
        end if
        call sourcein_(nrch, nptin, ndiffin, hydrau%flag, topo, hydrau, ptname, xptt, qptta, qptt, tepttmean, &
            tepttamp, tepttmaxtime, cpttmean, cpttamp, cpttmaxtime, phpttmean, &
            phpttamp, phpttmaxtime, diffname, xdup, xddn, qdifa, qdif, tedif, cdif, phind)

        if (ndiffin>=1) then
            deallocate(phind, stat=status(1)); deallocate(diffname, stat=status(2));
            deallocate(tedif, stat=status(3)); deallocate (qdif, stat=status(4));
            deallocate(qdifa, stat=status(5)); deallocate(xddn, stat=status(6))
            deallocate(xdup, stat=status(7))
        end if
        if (nptin>=1) then
            deallocate(phpttmaxtime, stat=status(13)); deallocate (phpttamp, stat=status(12))
            deallocate (phpttmean, stat=status(11)); deallocate(cpttmaxtime, stat=status(10))
            deallocate(cpttamp, stat=status(9)); deallocate (cpttmean, stat=status(8))
            deallocate (tepttmaxtime, stat=status(7)); deallocate (tepttamp, stat=status(6))
            deallocate (tepttmean, stat=status(5)); deallocate (qptt, stat=status(4))
            deallocate (qptta, stat=status(3)); deallocate (xptt, stat=status(2))
            deallocate (ptname, stat=status(1))
        end if

        !rates
        read(8,*) vss, mgc, mgn, mgp, mgd, mga
        read(8,*) tka, roc, ron
        read(8,*) ksocf, ksona, ksodn, ksop, ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von
        read(8,*) kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        read(8,*) kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, isat

        !gp 03-apr-08
        !read(8,*) khnx, va, typef, kgaf, tkgaf, kreaf, tkreaf, kexaf, tkexaf, kdeaf, abmax
        read(8,*) khnx, va, typef, kgaf, tkgaf, krea1f, krea2f, tkreaf, kexaf, tkexaf, kdeaf, abmax

        read(8,*) tkdeaf, ksnf, kspf, kscf, isatf, khnxf, kdt, tkdt, vdt

        !gp 26-jan-06
        !gp 10-jan-05 read(8,*) ninbmin, nipbmin, ninbupmax, nipbupmax, kqn, kqp
        !read(8,*) ninbmin, nipbmin, ninbupmax, nipbupmax, kqnratio, kqpratio
        read(8,*) ninbmin, nipbmin, ninbupmax, nipbupmax, kqnratio, kqpratio, nupwcfrac, pupwcfrac

        kqn = kqnratio * ninbmin
        kqp = kqpratio * nipbmin

        !gp 30-nov-04 read(8,*) kpath, tkpath, vpath
        read(8,*) kpath, tkpath, vpath, apath, kgen, tkgen, vgen, usegenericascod !gp 08-dec-04

        read(8,*) xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7

        !gp 15-nov-04 separate hyporheic rates from reaeration and add level 2 rates and hco3- use
        !gp 03-nov-04 read(8,*) kai, kawindmethod !oxygen reaeration option
        !gp read(8,*) kai, kawindmethod, typeh, kgah, tkgah, ksch, xdum8, kinhch
        read(8,*) kai, kawindmethod
        read(8,*) hco3use, hco3usef
        read(8,*) typeh, kgah, tkgah, ksch, xdum8, kinhch, kreah, tkreah, kdeah, tkdeah, ksnh, ksph, khnxh, ahmax

        !gp 17-jan-06
        !gp 17-nov-04 call makehydraulics(nrch, hydrau, downstreamboundary, kai)
        !call makehydraulics(nrch, hydrau, downstreamboundary, kai, topo%geomethod)

        stochrate= rates_t(nrch, hydrau, mgc,mgn, mgp, mgd, mga, vss,tka, roc, ron, ksocf, ksona, &
            ksodn, ksop, ksob, khc, tkhc, kdcs, tkdcs,kdc, tkdc, khn, tkhn, von, &
            kn, tkn, ki, tki, vdi, tvdi, khp, tkhp,vop, vip, kspi, kga, &
            tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, isat, khnx, &
            va, typef, kgaf, tkgaf, krea1f, krea2f, tkreaf, kexaf, tkexaf, kdeaf, &
            abmax, tkdeaf, ksnf, kspf, kscf, isatf, khnxf, kdt, tkdt, vdt, &
            ninbmin, nipbmin, ninbupmax, nipbupmax, kqn, kqp, &
            kpath, tkpath, vpath, pco2, xdum1, xdum2, xdum3, xdum4, xdum5, &
            xdum6,xdum7, kai, kawindmethod, &
            hco3use, hco3usef, &
            typeh, kgah, tkgah, ksch, xdum8, kinhch, &
            kreah, tkreah, kdeah, tkdeah, ksnh, ksph, khnxh, ahmax, &
            apath, kgen, tkgen, vgen, usegenericascod, &
            nupwcfrac, pupwcfrac)

        !day saving time
        read(8,*) dlstime
        sitesolar= sitesolar_(nrch, timezone, solarmethod, nfacbras, atcryanstolz, dlstime)
        !calculate sun rise and set
        call sunrisesunset(nrch, sitesolar, hydrau, system%today)

        !boundary data
        read(8,*) downstreamboundary !true/false

        !gp 17-jan-06
        call makehydraulics(nrch, hydrau, downstreamboundary, kai, topo%geomethod)

        !gp 12-jan-06
        !write(10,*) downstreamboundary
        !write(10,*) system%steadystate

        !headwaters
        if (system%steadystate) then
            allocate(hwin(0:hrsday-1), stat=status(1))
            if (status(1)==1) stop 'Insufficient memory, dynamic allocation failed!'
            read(8,*) (hwin(j)%te, j=0, hrsday-1), xxx
            do j = 1, nv - 2
                read(8,*) (hwin(k)%c(j), k=0, hrsday-1), xxx
            end do
            !hourly ph parameters
            read(8,*) (hwin(j)%ph, j=0, hrsday-1), xxx
        else
            !dynamic simulation
        end if

        !gp 11-jan-06 re-write with downstream boundary read same as headwater
        if (downstreamboundary) then
            if (system%steadystate) then
                allocate(dbin(0:hrsday-1), stat=status(1))
                if (status(1)==1) stop 'Insufficient memory, dynamic allocation failed!'
                read(8,*) (dbin(j)%te, j=0, hrsday-1), xxx

                !gp 12-jan-06
                !write(9,*) (dbin(j)%te, j=0, hrsday-1), xxx

                do j = 1, nv - 2
                    read(8,*) (dbin(k)%c(j), k=0, hrsday-1), xxx

                    !gp 12-jan-06
                    !write(10,*) (dbin(k)%c(j), k=0, hrsday-1), xxx

                end do
                !hourly ph parameters
                read(8,*) (dbin(j)%ph, j=0, hrsday-1), xxx

                !gp 12-jan-06
                !write(10,*) (dbin(j)%ph, j=0, hrsday-1), xxx

            else
                !dynamic simulation
            end if
            !else !downstream boundary is not specified
        end if

        hw= upstream_boundary_t(hwin)
        db= downstream_boundary_t(dbin,downstreamboundary)
        if (downstreamboundary) deallocate(dbin, stat=status(1))
        deallocate(hwin, stat=status(2))

        !gp 12-jan-06
        !write(10,*) 'done thru read downstreamboundary'

        !met data
        if (system%steadystate) then
            allocate(shadehh(0:hrsday-1,nrch), stat=status(1))
            allocate(tahh(0:hrsday-1,nrch), stat=status(2))
            allocate(tdhh(0:hrsday-1,nrch), stat=status(3))
            allocate(uwhh(0:hrsday-1,nrch), stat=status(4))
            allocate(cchh(0:hrsday-1,nrch), stat=status(5))

            !gp 16-jul-08
            allocate(solarhh(0:hrsday-1,nrch), stat=status(6))

            do i=1, 5
                if (status(i)==1) stop 'Insufficient memory, dynamic allocation failed!'
            end do

            do i = 1, nrch
                !do j = 0, hrsday-1
                read(8,*) (shadehh(j, i), j = 0, hrsday-1), xxx !reach hourly shade data
                !end do
                !read(8,*) xxx
            end do
            do i = 1, nrch
                !do j = 0, hrsday-1 !reach hourly air temperature data
                read(8,*) (tahh(j, i), j = 0, hrsday-1), xxx
                !end do
                !read(8,*) xxx
            end do
            do i = 1, nrch
                !do j = 0, hrsday-1
                read(8,*) (tdhh(j, i), j = 0, hrsday-1), xxx !reach dew point air temperature data
                !end do
                !read(8,*) xxx
            end do
            do i = 1, nrch
                !do j = 0, hrsday-1
                read(8,*) (uwhh(j, i), j = 0, hrsday-1), xxx !reach wind speed data
                !end do
                !read(8,*) xxx
            end do
            do i = 1, nrch
                !do j = 0, hrsday-1
                read(8,*) (cchh(j, i), j = 0, hrsday-1), xxx !reach cloud cover data
                !end do
                !read(8,*) xxx
            end do

            !gp 16-jul-08
            do i = 1, nrch
                read(8,*) (solarhh(j, i), j = 0, hrsday-1), xxx !reach solar data
            end do

        else !dyanmic simulation
        end if

        !gp 16-jul-08
        !sitemeteo = meteodata_(nrch, shadehh, tahh, tdhh, uwhh, cchh)
        sitemeteo = meteorology_t(nrch, shadehh, tahh, tdhh, uwhh, cchh, solarhh)

        !gp 16-jul-08
        deallocate(solarhh, stat=status(6))

        deallocate(cchh, stat=status(5))
        deallocate(uwhh, stat=status(4))
        deallocate(tdhh, stat=status(3))
        deallocate(tahh, stat=status(2))
        deallocate(shadehh, stat=status(1))

    end subroutine readinputfile


end module m_readfile

