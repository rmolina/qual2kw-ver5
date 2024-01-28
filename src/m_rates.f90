module m_rates
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    implicit none
    private
    public rates_t, setoxygeninhibenhance, tempadjust

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
        use class_hydraulics
        type(riverhydraulics_type), intent(inout) :: hydrau !assign reach-specific rates
        integer(i32), intent(in) :: nr !number of reach
        integer(i32) i

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
        rates%ikoxn= ikox(xdum2)
        rates%ksona = ksona

        !oxygen enhancement of denitrification
        rates%ikoxdn= ikox(xdum3)
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

        select case (xdum6)
          case ("Smith")
            rates%ilight = 2 !"langleys/d"
          case ("Steele")
            rates%ilight = 3 !"langleys/d"
          case default !"half saturation"
            rates%ilight = 1 !"langleys/d"
        end select

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

        select case (xdum7)
          case ("Smith")
            rates%ilightf = 2 !"langleys/d"
          case ("Steele")
            rates%ilightf = 3 !"langleys/d"
          case default !"half saturation"
            rates%ilightf = 1 ! "langleys/d"
        end select

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

        !gp 08-feb-06
        !if reach specific rates provided by input file are invalid, then use default general
        !rates provided by input file
        do i=1, nr
            !1.reaeration coefficient
            !leave the value as it is
            !unless it is equal/great than zero, it will be calculated by reaeration model
            !2.inorganic suspended solids settling velocity
            if (hydrau%reach(i)%vss <0) then
                hydrau%reach(i)%vss = rates%vss
            end if
            !3.slow cbod hydrolysis rate
            if (hydrau%reach(i)%khc <0) then
                hydrau%reach(i)%khc = rates%khc
            end if
            !4.slow cbod oxidation rate
            if (hydrau%reach(i)%kdcs <0) then
                hydrau%reach(i)%kdcs = rates%kdcs
            end if
            !5.fast cbod oxidation rate
            if (hydrau%reach(i)%kdc < 0) then
                hydrau%reach(i)%kdc = rates%kdc
            end if
            !6. organic n hydrolysis rate (khn)
            if (hydrau%reach(i)%khn < 0) then
                hydrau%reach(i)%khn = rates%khn
            end if
            !7. organic n settling velocity (von)
            if (hydrau%reach(i)%von < 0) then
                hydrau%reach(i)%von = rates%von
            end if
            !8. ammonium nitrification rate (kn)
            if (hydrau%reach(i)%kn < 0) then
                hydrau%reach(i)%kn = rates%kn
            end if
            !9. nitrate denitrification rate (ki)
            if (hydrau%reach(i)%ki < 0) then
                hydrau%reach(i)%ki = rates%ki
            end if
            !10. nitrate sed denitrification transfer coeff (vdi)
            if (hydrau%reach(i)%vdi < 0) then
                hydrau%reach(i)%vdi = rates%vdi
            end if
            !11. organic p hydrolysis (khp)
            if (hydrau%reach(i)%khp < 0) then
                hydrau%reach(i)%khp = rates%khp
            end if
            !12. organic p settling velocity (vop)
            if (hydrau%reach(i)%vop < 0) then
                hydrau%reach(i)%vop = rates%vop
            end if
            !13. inorganic p settling veloctiy (vip)
            if (hydrau%reach(i)%vip < 0) then
                hydrau%reach(i)%vip = rates%vip
            end if
            !14. phytoplankton max growth rate (kga)
            if (hydrau%reach(i)%kga < 0) then
                hydrau%reach(i)%kga = rates%kga
            end if
            !phytoplankton respiration rate (krea)
            if (hydrau%reach(i)%krea < 0) then
                hydrau%reach(i)%krea = rates%krea
            end if
            !phytoplankton death rate (kdea)
            if (hydrau%reach(i)%kdea < 0) then
                hydrau%reach(i)%kdea = rates%kdea
            end if
            !phytoplankton n half-sat (ksn)
            if (hydrau%reach(i)%ksn < 0) then
                hydrau%reach(i)%ksn = rates%ksn
            end if
            !phytoplankton p half-sat (ksp)
            if (hydrau%reach(i)%ksp < 0) then
                hydrau%reach(i)%ksp = rates%ksp
            end if
            !phytoplankton light sat (isat)
            if (hydrau%reach(i)%isat < 0) then
                hydrau%reach(i)%isat = rates%isat
            end if
            !phytoplankton ammonia pref (khnx)
            if (hydrau%reach(i)%khnx < 0) then
                hydrau%reach(i)%khnx = rates%khnx
            end if
            !phytoplankton settling (va)
            if (hydrau%reach(i)%va < 0) then
                hydrau%reach(i)%va = rates%va
            end if

            !bottom plant initial biomass (botalg0)
            !leave the value as it is
            !bottom plant max growth rate (kgaf)
            if (hydrau%reach(i)%kgaf < 0) then
                hydrau%reach(i)%kgaf = rates%kgaf

                !gp 03-apr-08 already in gd basis
                !else
                ! !convert to mga basis
                ! if (rates%typef == "zero-order") hydrau%reach(i)%kgaf = hydrau%reach(i)%kgaf * gd / mga

            end if

            !bottom plant first order carrying capacity (abmax)
            if (hydrau%reach(i)%abmax < 0) then
                hydrau%reach(i)%abmax = rates%abmax

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%abmax = hydrau%reach(i)%abmax * gd / mga !convert to mga basis

            end if

            !bottom plant respiration rate (kreaf)

            !gp 03-apr-08
            !if (hydrau%reach(i)%kreaf < 0) then
            ! hydrau%reach(i)%kreaf = rates%kreaf
            !end if
            if (hydrau%reach(i)%krea1f < 0) then
                hydrau%reach(i)%krea1f = rates%krea1f
            end if
            if (hydrau%reach(i)%krea2f < 0) then
                hydrau%reach(i)%krea2f = rates%krea2f
            end if

            !bottom plant excretion rate (kexaf)
            if (hydrau%reach(i)%kexaf < 0) then
                hydrau%reach(i)%kexaf = rates%kexaf
            end if
            !bottom plant death rate (kdeaf)
            if (hydrau%reach(i)%kdeaf < 0) then
                hydrau%reach(i)%kdeaf = rates%kdeaf
            end if
            !bottom plant external n half-sat (ksnf)
            if (hydrau%reach(i)%ksnf < 0) then
                hydrau%reach(i)%ksnf = rates%ksnf
            end if
            !bottom plant external p half-sat (kspf)
            if (hydrau%reach(i)%kspf < 0) then
                hydrau%reach(i)%kspf = rates%kspf
            end if
            !bottom plant light sat (isatf)
            if (hydrau%reach(i)%isatf < 0) then
                hydrau%reach(i)%isatf = rates%isatf
            end if
            !bottom plant nh4 pref (khnxf)
            if (hydrau%reach(i)%khnxf < 0) then
                hydrau%reach(i)%khnxf = rates%khnxf
            end if
            !bottom plant subsistence quota for n (ninbmin)
            if (hydrau%reach(i)%ninbmin < 0) then
                hydrau%reach(i)%ninbmin = rates%ninbmin

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%ninbmin = hydrau%reach(i)%ninbmin * mga / gd !convert to mga basis

            end if
            !bottom plant subsistence quota for p (nipbmin)
            if (hydrau%reach(i)%nipbmin < 0) then
                hydrau%reach(i)%nipbmin = rates%nipbmin

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%nipbmin = hydrau%reach(i)%nipbmin * mga / gd !convert to mga basis

            end if
            !bottom plant max uptake of n (ninbupmax)
            if (hydrau%reach(i)%ninbupmax < 0) then
                hydrau%reach(i)%ninbupmax = rates%ninbupmax

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%ninbupmax = hydrau%reach(i)%ninbupmax * mga / gd !convert to mga basis

            end if
            !bottom plant max uptake of p (nipbupmax)
            if (hydrau%reach(i)%nipbupmax < 0) then
                hydrau%reach(i)%nipbupmax = rates%nipbupmax

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%nipbupmax = hydrau%reach(i)%nipbupmax * mga / gd !convert to mga basis

            end if
            !bottom plant internal n half-sat (kqn)
            if (hydrau%reach(i)%kqn < 0) then
                hydrau%reach(i)%kqn = rates%kqn

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%kqn = hydrau%reach(i)%kqn * mga / gd !convert to mga basis

            end if
            !bottom plant internal p half-sat (kqp)
            if (hydrau%reach(i)%kqp < 0) then
                hydrau%reach(i)%kqp = rates%kqp

                !gp 03-apr-08 already in gd basis
                !else
                ! hydrau%reach(i)%kqp = hydrau%reach(i)%kqp * mga / gd !convert to mga basis

            end if
            !bottom plant fraction of n uptake from water column (nupwcfrac)
            if (hydrau%reach(i)%nupwcfrac < 0) then
                hydrau%reach(i)%nupwcfrac = rates%nupwcfrac
            end if
            !bottom plant fraction of p uptake from water column (pupwcfrac)
            if (hydrau%reach(i)%pupwcfrac < 0) then
                hydrau%reach(i)%pupwcfrac = rates%pupwcfrac
            end if
            !17. detritus(pom) dissolution rate (kdt)
            if (hydrau%reach(i)%kdt < 0) then
                hydrau%reach(i)%kdt = rates%kdt
            end if
            !18. detritus(pom) settling velocity (vdt)
            if (hydrau%reach(i)%vdt < 0) then
                hydrau%reach(i)%vdt = rates%vdt
            end if
            !pathogen dieoff rate (kpath)
            if (hydrau%reach(i)%kpath < 0) then
                hydrau%reach(i)%kpath = rates%kpath
            end if
            !pathogen settling rate (vpath)
            if (hydrau%reach(i)%vpath < 0) then
                hydrau%reach(i)%vpath = rates%vpath
            end if
            !pathogen light alpha (apath)
            if (hydrau%reach(i)%apath < 0) then
                hydrau%reach(i)%apath = rates%apath
            end if
            !hyporheic heterotroph max growth rate (kgah)
            if (hydrau%reach(i)%kgah < 0) then
                hydrau%reach(i)%kgah = rates%kgah
            end if
            !hyporheic heterotroph half-sat for cbod (ksch)
            if (hydrau%reach(i)%ksch < 0) then
                hydrau%reach(i)%ksch = rates%ksch
            end if
            !hyporheic heterotroph o2 inhibition parameter (kinhch)
            if (hydrau%reach(i)%kinhch < 0) then
                hydrau%reach(i)%kinhch = rates%kinhch
            end if
            !hyporheic heterotroph respiration rate (kreah)
            if (hydrau%reach(i)%kreah < 0) then
                hydrau%reach(i)%kreah = rates%kreah
            end if
            !hyporheic heterotroph death rate (kdeah)
            if (hydrau%reach(i)%kdeah < 0) then
                hydrau%reach(i)%kdeah = rates%kdeah
            end if
            !hyporheic heterotroph n half-sat (ksnh)
            if (hydrau%reach(i)%ksnh < 0) then
                hydrau%reach(i)%ksnh = rates%ksnh
            end if
            !hyporheic heterotroph p half-sat (ksph)
            if (hydrau%reach(i)%ksph < 0) then
                hydrau%reach(i)%ksph = rates%ksph
            end if
            !hyporheic heterotroph nh4 preference (khnxh)
            if (hydrau%reach(i)%khnxh < 0) then
                hydrau%reach(i)%khnxh = rates%khnxh
            end if
            !hyporheic heterotroph first-order carrying capacity (ahmax)
            if (hydrau%reach(i)%ahmax < 0) then
                hydrau%reach(i)%ahmax = rates%ahmax
            end if
            !generic constituent first-order loss rate constant (kgen)
            if (hydrau%reach(i)%kgen < 0) then
                hydrau%reach(i)%kgen = rates%kgen
            end if
            !generic constituent settling rate (vgen)
            if (hydrau%reach(i)%vgen < 0) then
                hydrau%reach(i)%vgen = rates%vgen
            end if
        end do

    end function

    !gp 03-apr-08
    !subroutine tempadjust(rates, hydrau, sitemeteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
    ! kreaft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
    ! kreat, vdit, kact, &
    ! kgent) !gp 30-nov-04 add kgent for generic constituent
    subroutine tempadjust(rates, hydrau, sitemeteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
        krea1ft, krea2ft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
        kreat, vdit, kact, &
        kgent)
        use class_hydraulics
        use m_meteorology

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
            hydrau%reach(i)%os = oxsat(te(i, 1), hydrau%reach(i)%elev)
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
    end subroutine tempadjust

    pure function oxsat(temp, elev)

        real(r64) oxsat
        real(r64), intent(in) :: temp, elev
        real(r64) taa, lnosf

        taa = temp + 273.15_r64
        lnosf = -139.34411_r64 + 157570.1_r64 / taa - 66423080.0_r64 / taa ** 2 &
            + 12438000000.0_r64 / taa ** 3 - 862194900000.0_r64 / taa ** 4
        oxsat = exp(lnosf) * (1 - 0.0001148_r64 * elev)

    end function oxsat

    subroutine setoxygeninhibenhance(rates, o2, fcarb, fnitr, fdenitr, frespp, frespb)

        real(r64), intent(out) :: fcarb, fnitr, fdenitr, frespp, frespb
        type(rates_t), intent(in) :: rates
        real(r64), intent(in) :: o2

        fcarb = setoxygeninhibenhance_(rates%ikoxc, rates%ksocf, o2)
        fnitr = setoxygeninhibenhance_(rates%ikoxn, rates%ksona, o2)
        fdenitr = 1 - setoxygeninhibenhance_(rates%ikoxdn, rates%ksodn, o2)
        frespp = setoxygeninhibenhance_(rates%ikoxp, rates%ksop, o2)
        frespb = setoxygeninhibenhance_(rates%ikoxb, rates%ksob, o2)

    end subroutine setoxygeninhibenhance

    real(r64) function setoxygeninhibenhance_(model, rate, oxygen)
        integer(i32), intent(in) ::  model
        real(r64), intent(in) :: rate, oxygen
        !real(r64), intent(out) :: res_frespb
        select case (model)
          case (1)
            setoxygeninhibenhance_ = oxygen / (rate + oxygen)
          case (2)
            setoxygeninhibenhance_ = 1.0_r64 - exp(-rate * oxygen)
          case (3)
            setoxygeninhibenhance_ = oxygen ** 2 / (rate + oxygen ** 2)
        end select
    end function setoxygeninhibenhance_

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
    end function

end module m_rates
