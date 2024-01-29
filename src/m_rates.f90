module m_rates
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_hydraulics, only: riverhydraulics_type
    implicit none
    private
    public rates_t

    type rates_t
        real(r64) mga, mgd !scc, for v1_3
        real(r64) vss !inorganic suspended solids settling vol
        real(r64) anc, apc, adc
        real(r64) ana, apa, ada, aca
        real(r64) ronp, ano
        real(r64) :: roc =2.67 !roc - o2 for carbon oxidation
        real(r64) :: ron =4.57 !ron - o2 for nh4 nitrification
        real(r64) roa !roa - o2 for chlorophyll
        real(r64) :: tka = 1.024 !temp correction for reaeration

        real(r64) ralkda, ralkdn
        real(r64) racc
        real(r64) acca, acco, accd, accc

        character(len=30) kai !reaeration model
        character(len=30) kawindmethod !reaeartion wind effect
        integer(i32) ikoxc !oxygen inhibition cbod oxidation model
        integer(i32) ikoxn !oxygen inhibition nitrification model
        integer(i32) ikoxdn !oxygen enhancement of denitrification model
        integer(i32) ikoxp !oxygen inhibition of phytoplankton respiration
        integer(i32) ikoxb !
        integer(i32) ilight !light model

        real(r64) ksocf !oxygen inhibition cbod oxidation parameter
        real(r64) ksona !oxygen inhibition nitrification parameter
        real(r64) ksodn !oxygen enhancement of denitrification parameter
        real(r64) ksop
        real(r64) ksob

        !slow cbod
        real(r64) khc !slow cbod hydrolysis rate
        real(r64) :: tkhc =1.05 !slow cbod hydrolysis temp correction
        real(r64) kdcs !slow cbod oxidation rate
        real(r64) :: tkdcs =1.05 !slow cbod oxidation temp correction

        !fast cbod
        real(r64) kdc !fast cbod oxidation rate
        real(r64) :: tkdc=1.05 !fast cbod oxidation temp correction

        !nitrogen
        real(r64) khn !organic n hydrolysis rate
        real(r64) :: tkhn=1.05 !organic n hydrolysis temp correction
        real(r64) von !organic n settling velocity
        real(r64) kn !ammonium nitrification rate
        real(r64) :: tkn=1.05 !ammonium nitrification temp correction
        real(r64) ki !nitrate denitrification
        real(r64) :: tki=1.05 !nitrate denitrification rate
        real(r64) vdi !nitrate sed denitrification transfer coeff
        real(r64) :: tvdi=1.05 !nitrate sed denitrification transfer coeff temp correction

        !phosphorus
        real(r64) khp !organic p hydrolysis
        real(r64) :: tkhp=1.05 !organic p hydrolysis temp correction
        real(r64) vop !organic p settling velocity
        real(r64) vip !inorganic p settling velocity
        real(r64) kspi !sed p oxygen attenuation half sat constant

        !phytoplanton
        real(r64) kga !phytoplankton max growth rate
        real(r64) :: tkga=1.066 !phytoplankton max growth temp correction
        real(r64) krea !phytoplankton respiration rate
        real(r64) :: tkrea=1.05 !phytoplankton respiration temp correction
        real(r64) kdea !phytoplankton death rate
        real(r64) :: tkdea=1.05 !phytoplankton death temp correction
        real(r64) ksn !nitrogen half sat constant
        real(r64) ksp !phosphorus half sat constant
        real(r64) ksc !inorganic carbon half sat
        real(r64) isat !light constant
        real(r64) :: khnx=15.0 !ammonia preference
        real(r64) va !settling velocity of phytoplankton

        !bottom algae
        real(r64) kgaf !max growth rate
        real(r64) :: tkgaf=1.066 !max growth temp correction

        !gp 03-apr-08
        !real(r64) kreaf !botalg respiration rate
        real(r64) krea1f !botalg basal respiration rate
        real(r64) krea2f !botalg photo respiration rate

        real(r64) :: tkreaf=1.05 !phytoplankton respiration temp correction
        real(r64) kdeaf !phytoplankton death rate
        real(r64) :: tkdeaf=1.05 !phytoplankton death temp correction
        real(r64) ksnf !nitrogen half sat constant
        real(r64) kspf !phosphorus half sat constant
        real(r64) kscf !inorganic carbon half sat
        real(r64) abmax !first-order model carrying capacity
        real(r64) kexaf !excretion rate
        real(r64) :: tkexaf =1.05 !excretion rate temp correction

        integer(i32) ilightf !light model
        real(r64) isatf !light constant
        real(r64) :: khnxf=15.0 !ammonia preference
        character(len=30) :: typef ='Zero-order' !bottom algae growth model,zero-first order

        !pom
        real(r64) kdt !pom dissolution rate
        real(r64) :: tkdt=1.05 !pom dissolution temp correction
        real(r64) vdt !pom settling velocity
        !luxury uptake
        real(r64) ninbmin
        real(r64) nipbmin
        real(r64) ninbupmax
        real(r64) nipbupmax
        real(r64) kqn
        real(r64) kqp

        !gp 26-jan-06
        real(r64) nupwcfrac
        real(r64) pupwcfrac

        !pathogens
        real(r64) kpath !decay
        real(r64) :: tkpath=1.07 !decay temp correction
        real(r64) vpath !settling velocity

        !gp 30-nov-04 new rates
        real(r64) apath !alpha constant for light mortality of pathogen indicator
        real(r64) kgen !decay of generic constituent
        real(r64) :: tkgen=1.07 !decay temp correction for generic constituent
        real(r64) vgen !settling velocity of generic constituent

        !gp 08-dec-04
        character(len=30) usegenericascod

        !ph
        real(r64) pco2 !partial pressure of carbon dioxide
        real(r64) ralkaa, ralkan, ralkbn, ralkbp
        real(r64) ralkden, rondn, ralkn
        real(r64) rcca, rcco, rccd, rccc

        !gp 03-nov-04
        character(len=30) typeh, xdum8
        real(r64) :: tkgah=1.047
        real(r64) kgah, ksch, kinhch
        integer(i32) ikoxch

        !gp 15-nov-04
        character(len=30) hco3use, hco3usef

        !gp 15-nov-04
        real(r64) kreah, kdeah, ksnh, ksph, khnxh, ahmax
        real(r64) :: tkreah=1.07, tkdeah=1.07

    end type rates_t

    interface rates_t
        procedure :: rates_ctor
    end interface rates_t

contains

    function rates_ctor(nr, hydrau, mgc, mgn, mgp, mgd, mga, vss,tka, roc, ron, ksocf, ksona, ksodn, ksop, &
        ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von, kn, &
        tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi, kga, tkga, krea, &
        tkrea, kdea, tkdea, ksn, ksp, ksc, isat, khnx, va, typef, kgaf, &
        tkgaf, krea1f, krea2f, tkreaf, kexaf, tkexaf, kdeaf, abmax, tkdeaf, ksnf, &
        kspf, kscf, isatf, khnxf, kdt, tkdt, vdt, ninbmin, nipbmin, &
        ninbupmax, nipbupmax, kqn, kqp, kpath, tkpath, vpath, pco2, &
        xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, &
        kawindmethod, hco3use, hco3usef, typeh, kgah, tkgah, ksch, xdum8, kinhch, &
        kreah, tkreah, kdeah, tkdeah, ksnh, ksph, khnxh, ahmax, &
        apath, kgen, tkgen, vgen, usegenericascod, &
        nupwcfrac, pupwcfrac) result(rates)


        !gp 08-feb-06
        !this function also assigns global rates to the unspecified reach-specific rates in hydrau
        type(riverhydraulics_type), intent(inout) :: hydrau !assign reach-specific rates
        integer(i32), intent(in) :: nr !number of reach

        type(rates_t) rates

        !stoichiometry
        real(r64) :: mgc, mgn, mgp, mgd, mga
        real(r64), intent(in) :: vss, tka, roc, ron
        real(r64), intent(in) :: ksocf, ksona, ksodn, ksop, ksob, khc, kdcs, tkdcs, tkhc, kdc, tkdc, khn
        real(r64), intent(in) :: tkhn, von, kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        real(r64), intent(in) :: kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, isat

        !gp 03-apr-08
        !real(r64), intent(in) :: khnx, va, kgaf, tkgaf, kreaf, tkreaf, kexaf, tkexaf, kdeaf
        real(r64), intent(in) :: khnx, va, kgaf, tkgaf, krea1f, krea2f, tkreaf, kexaf, tkexaf, kdeaf

        real(r64), intent(in) :: tkdeaf, abmax, ksnf, kspf, kscf, isatf, khnxf, kdt, tkdt, vdt
        real(r64), intent(in) :: ninbmin, nipbmin, ninbupmax, nipbupmax, kqn, kqp

        !gp 26-jan-06
        real(r64), intent(in) :: nupwcfrac, pupwcfrac

        real(r64), intent(in) :: kpath, tkpath, vpath, pco2

        !gp 30-nov-04
        real(r64), intent(in) :: apath, kgen, tkgen, vgen !gp 30-nov-04 new paramters for pathogen and generic constituent

        !gp 08-dec-04
        character(len=30), intent(in) :: usegenericascod

        !org c, n inhition model, denitrification enhance model, algae light model, perihyte light model
        character(len=30), intent(in) :: xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, kawindmethod

        !gp 03-nov-04
        character(len=30), intent(in) :: typeh, xdum8
        real(r64), intent(in) :: kgah, tkgah, ksch, kinhch

        !gp 15-nov-04
        character(len=30), intent(in) :: hco3use, hco3usef

        !gp 15-nov-04
        real(r64), intent(in) :: kreah, tkreah, kdeah, tkdeah, ksnh, ksph, khnxh, ahmax


        character(len=30) :: typef !bottom algae growth model,zero-first order

        real(r64) gc, gd

        if (mgc <= 0) mgc = 40.0_r64
        if (mgn <= 0) mgn = 7.2_r64
        if (mgp <= 0) mgp = 1.0_r64
        if (mgd <= 0) mgd = 100.0_r64
        if (mga <= 0) mga = 1.0_r64

        rates%mga= mga; rates%mgd=mgd !scc, for v1_3
        !suspended solids
        rates%vss= vss

        !oxygen
        select case (kai)
            !case 'internal'
          case ("O'Connor-Dobbins", &
              "Churchill", &
              "Owens-Gibbs", &
              "Thackston-Dawson", &
              "Tsivoglou-Neal", &
              "USGS(pool-riffle)", &
              "USGS(channel-control)")
            rates%kai = kai
          case default
            rates%kai = "internal"
        end select

        rates%kawindmethod = kawindmethod
        if (tka>0) rates%tka = tka
        if (roc>0) rates%roc = roc
        if (ron>0) rates%ron = ron

        !oxygen inhibition of carbon oxidation
        rates%ikoxc= ikox(xdum1)
        rates%ksocf = ksocf

        !oxygen inhibition of nitrification
        rates%ikoxn = ikox(xdum2)
        rates%ksona = ksona

        !oxygen enhancement of denitrification
        rates%ikoxdn = ikox(xdum3)
        rates%ksodn = ksodn

        !oxygen inhibition of phytoplankton respiration
        rates%ikoxp = ikox(xdum4)
        rates%ksop= ksop

        !oxygen inhibition of bottom plant respiration
        rates%ikoxb = ikox(xdum5)
        rates%ksob= ksob

        !slow cbod
        rates%khc = khc
        if (tkhc>0) rates%tkhc = tkhc
        rates%kdcs= kdcs
        if (tkdcs>0) rates%tkdcs = tkdcs

        !fast cbod
        rates%kdc = kdc
        if (tkdc>0) rates%tkdc = tkdc

        !organic n
        rates%khn = khn
        if (tkhn>0) rates%tkhn = tkhn
        !organic n settling vel
        rates%von =von

        !ammonium
        rates%kn = kn
        if (tkn>0) rates%tkn = tkn

        !nitrate
        rates%ki = ki
        if (tki > 0) rates%tki=tki
        rates%vdi = vdi
        if (tvdi>0) rates%tvdi = tvdi

        !organic p
        rates%khp = khp
        if (tkhp>0) rates%tkhp = tkhp

        rates%vop = vop
        rates%vip = vip
        rates%kspi = kspi

        !phytoplankton
        rates%kga = kga
        if (tkga>0) rates%tkga =tkga

        rates%krea = krea
        if (tkrea>0) rates%tkrea = tkrea

        rates%kdea = kdea
        if (tkdea>0) rates%tkdea = tkdea

        rates%ksn = ksn
        rates%ksp = ksp
        rates%ksc = ksc

        rates%ilight = ilight(xdum6)

        rates%isat = isat
        if (khnx > 0) rates%khnx = khnx
        rates%va = va

        !bottom algae
        rates%kgaf = kgaf
        if (tkgaf>0) rates%tkgaf = tkgaf

        !gp 03-apr-08
        !rates%kreaf = kreaf
        rates%krea1f = krea1f
        rates%krea2f = krea2f

        if (tkreaf>0) rates%tkreaf = tkreaf

        rates%kexaf = kexaf !excretion rate
        if (tkexaf > 0) rates%tkexaf = tkexaf

        rates%kdeaf = kdeaf
        if (tkdeaf>0) rates%tkdeaf= tkdeaf
        rates%ksnf = ksnf
        rates%kspf = kspf
        rates%kscf = kscf

        !rates%abmax= abmax !scc 08/09/2004

        if (typef/="First-order") then
            typef="Zero-order"
        end if
        rates%typef = typef

        gc = mgc / 1000.0_r64
        gd = mgd / 1000.0_r64


        !gp 03-apr-08 starting with version b42a02, these inputs are now per unit dry weight instead of chl a
        ! !convert bottom algae rates to a mga basis
        ! if (typef == "zero-order") rates%kgaf = kgaf * gd / mga !scc 08/09/2004
        ! rates%abmax = abmax * gd / mga !scc 08/09/2004
        ! rates%ninbmin = ninbmin * mga / gd !scc 08/09/2004
        ! rates%nipbmin = nipbmin * mga / gd !scc 08/09/2004
        ! rates%ninbupmax = ninbupmax * mga / gd !scc 08/09/2004
        ! rates%nipbupmax = nipbupmax * mga / gd !scc 08/09/2004
        ! rates%kqn = kqn * mga / gd !scc 08/09/2004
        ! rates%kqp = kqp * mga / gd !scc 08/09/2004
        if (typef == "Zero-order") rates%kgaf = kgaf !gd/m^2/day
        rates%abmax = abmax !gd/m^2
        rates%ninbmin = ninbmin !mgn/gd
        rates%nipbmin = nipbmin !mgp/gd
        rates%ninbupmax = ninbupmax !mgn/gd/day
        rates%nipbupmax = nipbupmax !mgp/gd/day
        rates%kqn = kqn !mgn/gd
        rates%kqp = kqp !mgp/gd

        ! select case (xdum7)
        !   case ("Smith")
        !     rates%ilightf = 2 !"langleys/d"
        !   case ("Steele")
        !     rates%ilightf = 3 !"langleys/d"
        !   case default !"half saturation"
        !     rates%ilightf = 1 ! "langleys/d"
        ! end select
        rates%ilightf = ilight(xdum7)

        rates%isatf = isatf
        if (khnxf > 0) rates%khnxf = khnxf

        !rates%ninbmin= ninbmin !scc 08/09/2004
        !rates%nipbmin= nipbmin !scc 08/09/2004
        !rates%ninbupmax =ninbupmax !scc 08/09/2004
        !rates%nipbupmax =nipbupmax !scc 08/09/2004
        !rates%kqn = kqn !scc 08/09/2004
        !rates%kqp =kqp !scc 08/09/2004

        !gp 26-jan-06
        rates%nupwcfrac = nupwcfrac
        rates%pupwcfrac = pupwcfrac

        !detritus (pom)
        rates%kdt = kdt
        if (tkdt>0) rates%tkdt = tkdt
        rates%vdt = vdt

        !pathogens
        rates%kpath = kpath
        if (tkpath>0) rates%tkpath = tkpath
        rates%vpath = vpath

        !gp 30-nov-04 new param for pathogens and generic constituent
        rates%apath = apath
        rates%kgen = kgen
        if (tkgen>0) rates%tkgen = tkgen
        rates%vgen = vgen

        !gp 08-dec-04
        rates%usegenericascod = usegenericascod

        !ph
        rates%pco2 = pco2 ! / 1.0e6

        !gp 15-nov-04 hco3- use by phytoplankton and bottom algae
        rates%hco3use = hco3use
        rates%hco3usef = hco3usef

        !gp 03-nov-04 parameters for hyporheic heterotrophic biofilm oxidation of fast cbod
        rates%typeh = typeh
        rates%kgah = kgah
        if (tkgah>0) rates%tkgah = tkgah
        rates%ksch = ksch
        rates%ikoxch= ikox(xdum8)
        rates%kinhch = kinhch

        !gp 15-nov-04 rates for level 2 hyporheic biofilm growth
        rates%kreah = kreah
        if (tkreah>0) rates%tkreah = tkreah
        rates%kdeah = kdeah
        if (tkdeah>0) rates%tkdeah = tkdeah
        rates%ksnh = ksnh
        rates%ksph = ksph
        rates%khnxh = khnxh
        rates%ahmax = ahmax

        rates%anc = mgn / gc
        rates%apc = mgp / gc
        rates%ron = ron / 1000.0_r64
        rates%adc = mgd / mgc
        rates%roa = roc * gc / mga
        rates%ana = mgn / mga
        rates%apa = mgp / mga
        rates%aca = gc / mga
        rates%ada = gd / mga

        !gp 03-dec-09
        !!ph/inorganic carbon stoichiometry
        !rates%ralkaa = 14.0_r64 / 106.0_r64 / 12.0_r64 * rates%aca / 1000.0_r64
        !rates%ralkan = 18.0_r64 / 106.0_r64 / 12.0_r64 * rates%aca / 1000.0_r64
        !rates%ralkda = 14.0_r64 / 106.0_r64 / 12.0_r64 / rates%adc / 1000.0_r64
        !rates%ralkdn = 18.0_r64 / 106.0_r64 / 12.0_r64 / rates%adc / 1000.0_r64
        !rates%ralkn = 2.0_r64 / 14.0_r64 / 1000.0_r64 / 1000.0_r64
        !rates%ralkden = 4.0_r64 / (4.0_r64 * 14.0_r64) / 1000.0_r64 / 1000.0_r64
        !rates%racc = 14.0_r64 / 106.0_r64 / 12.0_r64 / 1000.0_r64
        !rates%rondn = roc * 5.0_r64 * 12.0_r64 / (4.0_r64 * 14.0_r64) / 1000.0_r64
        !rates%acca = gc / mga / 12.0_r64 / 1000.0_r64
        !rates%acco = 1.0_r64 / 12.0_r64 / roc / 1000.0_r64
        !rates%accd = gc / gd / 12.0_r64 / 1000.0_r64
        !rates%accc = 1.0_r64 / 12.0_r64 / 1000.0_r64
        !
        !!ph/inorganic carbon stoichiometry
        !rates%ralkaa = 14.0_r64 / 106.0_r64 / 12.0_r64 * rates%aca / 1000.0_r64
        !rates%ralkan = 18.0_r64 / 106.0_r64 / 12.0_r64 * rates%aca / 1000.0_r64
        !rates%ralkbn = 16.0_r64 / 16.0_r64 / 14.0_r64 / 1000.0_r64 / 1000.0_r64
        !rates%ralkbp = 2.0_r64 / 1.0_r64 / 31.0_r64 / 1000.0_r64 / 1000.0_r64
        !rates%ralkn = 2.0_r64 / 14.0_r64 / 1000.0_r64 / 1000.0_r64
        !rates%ralkden = 4.0_r64 / (4.0_r64 * 14.0_r64) / 1000.0_r64 / 1000.0_r64
        !rates%rondn = rates%roc * 5.0_r64 * 12.0_r64 / (4.0_r64 * 14.0_r64) / 1000.0_r64
        !rates%rcca = gc / mga / 12.0_r64 / 1000.0_r64
        !rates%rcco = 1.0_r64 / 12.0_r64 / rates%roc / 1000.0_r64
        !rates%rccd = gc / gd / 12.0_r64 / 1000.0_r64
        !rates%rccc = 1.0_r64 / 12.0_r64 / 1000.0_r64
        rates%ralkbn = 1.0_r64 / 14.0067_r64 / 1000.0_r64 / 1000.0_r64
        rates%ralkbp = 1.0_r64 / 30.973762_r64 / 1000.0_r64 / 1000.0_r64 !multiplied later by speciation of phosphate in sub derivs
        rates%rondn = rates%roc * 5.0_r64 * 12.0107_r64 / (4.0_r64 * 14.0067_r64) / 1000.0_r64
        rates%rcca = gc / mga / 12.0107_r64 / 1000.0_r64
        rates%rcco = 1.0_r64 / 12.0107_r64 / rates%roc / 1000.0_r64
        rates%rccd = gc / gd / 12.0107_r64 / 1000.0_r64
        rates%rccc = 1.0_r64 / 12.0107_r64 / 1000.0_r64

        !if reach specific rates provided by input file are invalid,
        !then use default general rates provided by input file
        call validate_user_rates(nr, hydrau, rates)

    end function rates_ctor

    subroutine validate_user_rates(nr, hydrau, rates)
        integer(i32), intent(in) :: nr !number of reach
        type(rates_t), intent(in) :: rates
        type(riverhydraulics_type), intent(inout) :: hydrau

        integer(i32) i

        do i=1, nr
            !1.reaeration coefficient
            !leave the value as it is
            !unless it is equal/great than zero, it will be calculated by reaeration model

            !inorganic suspended solids
            call validate_user_rate(hydrau%reach(i)%vss, rates%vss) !settling velocity

            !slow cbod
            call validate_user_rate(hydrau%reach(i)%khc, rates%khc) !hydrolysis rate
            call validate_user_rate(hydrau%reach(i)%kdcs, rates%kdcs) !oxidation rate

            !fast cbod 
            call validate_user_rate(hydrau%reach(i)%kdc, rates%kdc) !oxidation rate (kdc)

            !organic n
            call validate_user_rate(hydrau%reach(i)%khn, rates%khn) !hydrolysis rate (khn)
            call validate_user_rate(hydrau%reach(i)%von, rates%von) !settling velocity (von)

            !ammonium 
            call validate_user_rate(hydrau%reach(i)%kn, rates%kn) !nitrification rate (kn)

            !nitrate
            call validate_user_rate(hydrau%reach(i)%ki, rates%ki) !denitrification rate (ki)
            call validate_user_rate(hydrau%reach(i)%vdi, rates%vdi) !sed denitrification transfer coeff (vdi)

            !organic p 
            call validate_user_rate(hydrau%reach(i)%khp, rates%khp) !hydrolysis (khp)
            call validate_user_rate(hydrau%reach(i)%vop, rates%vop) !settling velocity (vop)

            !inorganic p 
            call validate_user_rate(hydrau%reach(i)%vip, rates%vip) !settling veloctiy (vip)

            !phytoplankton
            call validate_user_rate(hydrau%reach(i)%kga, rates%kga) !max growth rate (kga)
            call validate_user_rate(hydrau%reach(i)%krea, rates%krea) !respiration rate (krea)
            call validate_user_rate(hydrau%reach(i)%kdea, rates%kdea) !death rate (kdea)
            call validate_user_rate(hydrau%reach(i)%ksn, rates%ksn) !n half-sat (ksn)
            call validate_user_rate(hydrau%reach(i)%ksp, rates%ksp) !p half-sat (ksp)
            call validate_user_rate(hydrau%reach(i)%isat, rates%isat) !light sat (isat)
            call validate_user_rate(hydrau%reach(i)%khnx, rates%khnx) !ammonia pref (khnx)
            call validate_user_rate(hydrau%reach(i)%va, rates%va) !settling (va)

            !bottom plant
            !initial biomass (botalg0) ? !leave the value as it is
            call validate_user_rate(hydrau%reach(i)%kgaf, rates%kgaf) !max growth rate (kgaf)
            call validate_user_rate(hydrau%reach(i)%abmax, rates%abmax) !first order carrying capacity (abmax)
            call validate_user_rate(hydrau%reach(i)%krea1f, rates%krea1f) !respiration rate 1 (krea1f)
            call validate_user_rate(hydrau%reach(i)%krea2f, rates%krea2f) !respiration rate 2 (krea2f)
            call validate_user_rate(hydrau%reach(i)%kexaf, rates%kexaf) !excretion rate (kexaf)
            call validate_user_rate(hydrau%reach(i)%kdeaf, rates%kdeaf) !death rate (kdeaf)
            call validate_user_rate(hydrau%reach(i)%ksnf, rates%ksnf) !external n half-sat (ksnf)
            call validate_user_rate(hydrau%reach(i)%kspf, rates%kspf) !external p half-sat (kspf)
            call validate_user_rate(hydrau%reach(i)%isatf, rates%isatf) !light sat (isatf)
            call validate_user_rate(hydrau%reach(i)%khnxf, rates%khnxf) !nh4 pref (khnxf)
            call validate_user_rate(hydrau%reach(i)%ninbmin, rates%ninbmin) !subsistence quota for n (ninbmin)
            call validate_user_rate(hydrau%reach(i)%nipbmin, rates%nipbmin) !subsistence quota for p (nipbmin)
            call validate_user_rate(hydrau%reach(i)%ninbupmax, rates%ninbupmax) !max uptake of n (ninbupmax)
            call validate_user_rate(hydrau%reach(i)%nipbupmax, rates%nipbupmax) !max uptake of p (nipbupmax)
            call validate_user_rate(hydrau%reach(i)%kqn, rates%kqn) !internal n half-sat (kqn)
            call validate_user_rate(hydrau%reach(i)%kqp, rates%kqp) !internal p half-sat (kqp)
            call validate_user_rate(hydrau%reach(i)%nupwcfrac, rates%nupwcfrac) !fraction of n uptake from water column (nupwcfrac)
            call validate_user_rate(hydrau%reach(i)%pupwcfrac, rates%pupwcfrac) !fraction of p uptake from water column (pupwcfrac)

            !detritus(pom)
            call validate_user_rate(hydrau%reach(i)%kdt, rates%kdt) !dissolution rate (kdt)
            call validate_user_rate(hydrau%reach(i)%vdt, rates%vdt) !settling velocity (vdt)

            !pathogen
            call validate_user_rate(hydrau%reach(i)%kpath, rates%kpath) !dieoff rate (kpath)
            call validate_user_rate(hydrau%reach(i)%vpath, rates%vpath) !settling rate (vpath)
            call validate_user_rate(hydrau%reach(i)%apath, rates%apath) !light alpha (apath)

            !hyporheic heterotroph
            call validate_user_rate(hydrau%reach(i)%kgah, rates%kgah) !max growth rate (kgah)
            call validate_user_rate(hydrau%reach(i)%ksch, rates%ksch) !half-sat for cbod (ksch)
            call validate_user_rate(hydrau%reach(i)%kinhch, rates%kinhch) !o2 inhibition parameter (kinhch)
            call validate_user_rate(hydrau%reach(i)%kreah, rates%kreah) !respiration rate (kreah)
            call validate_user_rate(hydrau%reach(i)%khnxh, rates%khnxh) !death rate (kdeah)
            call validate_user_rate(hydrau%reach(i)%khnxh, rates%khnxh) !n half-sat (ksnh)
            call validate_user_rate(hydrau%reach(i)%khnxh, rates%khnxh) !p half-sat (ksph)
            call validate_user_rate(hydrau%reach(i)%khnxh, rates%khnxh) !nh4 preference (khnxh)
            call validate_user_rate(hydrau%reach(i)%ahmax, rates%ahmax) !first-order carrying capacity (ahmax)

            !generic constituent
            call validate_user_rate(hydrau%reach(i)%kgen, rates%kgen) ! first-order loss rate constant (kgen)
            call validate_user_rate(hydrau%reach(i)%vgen, rates%vgen) ! settling rate (vgen)
        end do

    end subroutine validate_user_rates


    pure subroutine validate_user_rate(user_rate, default_rate)
        real(r64), intent(in) :: default_rate
        real(r64), intent(inout) :: user_rate

        if (user_rate < 0) then
            user_rate = default_rate
        end if
    end subroutine validate_user_rate


    pure function ikox(xdum)
        character(len=30), intent(in) :: xdum
        integer(i32) ikox

        select case (xdum)
          case ("Half saturation")
            ikox = 1 !"mgO2/L"
          case ("2nd order")
            ikox = 3 !"mgO2/L"
          case default !"Exponential"
            ikox = 2 !"L/mgO2"
        end select
    end function ikox


    pure function ilight(xdum)
        character(len=30), intent(in) :: xdum
        integer(i32) ilight

        select case (xdum)
          case ("Smith")
            ilight = 2 !"langleys/d"
          case ("Steele")
            ilight = 3 !"langleys/d"
          case default !"half saturation"
            ilight = 1 !"langleys/d"
        end select
    end function ilight


end module m_rates
