! classstoch.f90

!stoichiometry
!INCLUDE "NRTYPE.F90"
MODULE Class_Rates
    USE nrtype
    IMPLICIT NONE

    PRIVATE
    PUBLIC  Rates_, Rates_type, SetOxygenInhibEnhance, TempAdjust
    TYPE Rates_type
        REAL(DP) mgA, mgD							!SCC, for v1_3
        REAL(DP) vss									!inorganic suspended solids settling vol
        REAL(DP) anc, apc, adc
        REAL(DP) ana, apa, ada, aca
        REAL(DP) ronp, ano
        REAL(DP) :: roc =2.67					!roc - O2 for carbon oxidation
        REAL(DP) :: ron =4.57					!ron - O2 for NH4 nitrification
        REAL(DP) roa									!roa - O2 for Chlorophyll
        REAL(DP) :: tka = 1.024				!temp correction for reaeration

        REAL(DP) ralkda, ralkdn
        REAL(DP) racc
        REAL(DP) acca, acco, accd, 	accc

        CHARACTER(LEN=30) kai					!Reaeration model
        CHARACTER(LEN=30) kawindmethod	!Reaeartion wind effect
        INTEGER(I4B) IkoxC						!oxygen inhibition CBOD oxidation model
        INTEGER(I4B) IkoxN						!oxygen inhibition nitrification model
        INTEGER(I4B) IkoxDN						!Oxygen enhancement of denitrification model
        INTEGER(I4B) IkoxP						!oxygen inhibition of Phytoplankton respiration
        INTEGER(I4B) IkoxB						!
        INTEGER(I4B) Ilight						!light model

        REAL(DP) Ksocf								!oxygen inhibition CBOD oxidation parameter
        REAL(DP) Ksona								!oxygen inhibition nitrification parameter
        REAL(DP) Ksodn								!Oxygen enhancement of denitrification parameter
        REAL(DP) Ksop
        REAL(DP) Ksob


        !slow CBOD
        REAL(DP) khc									!slow CBOD hydrolysis rate
        REAL(DP) :: tkhc =1.05				!slow CBOD hydrolysis temp correction
        REAL(DP) kdcs									!slow CBOD oxidation rate
        REAL(DP) :: tkdcs =1.05				!slow CBOD oxidation temp correction
        !Fast CBOD
        REAL(DP) kdc									!fast CBOD oxidation rate
        REAL(DP) :: tkdc=1.05					!fast CBOD oxidation temp correction
        !Nitrogen
        REAL(DP) khn									!organic N hydrolysis rate
        REAL(DP) :: tkhn=1.05					!organic N hydrolysis temp correction
        REAL(DP) von									!Organic N settling velocity
        REAL(DP) kn										!Ammonium nitrification rate
        REAL(DP) :: tkn=1.05					!Ammonium nitrification temp correction
        REAL(DP) ki										!Nitrate denitrification
        REAL(DP) :: tki=1.05					!Nitrate denitrification rate
        REAL(DP) vdi									!Nitrate sed denitrification transfer coeff
        REAL(DP) :: tvdi=1.05					!Nitrate sed denitrification transfer coeff temp correction
        !phosphorus
        REAL(DP) khp									!organic P hydrolysis
        REAL(DP) :: tkhp=1.05					!organic P hydrolysis temp correction
        REAL(DP) vop									!organic P settling velocity
        REAL(DP) vip									!inorganic P settling velocity
        REAL(DP) kspi									!Sed P oxygen attenuation half sat constant

        !Phytoplanton
        REAL(DP) kga									!Phytoplankton max growth rate
        REAL(DP) :: tkga=1.066				!phytoplankton max growth temp correction
        REAL(DP) krea									!phytoplankton respiration rate
        REAL(DP) :: tkrea=1.05				!phytoplankton respiration temp correction
        REAL(DP) kdea									!phytoplankton death rate
        REAL(DP) :: tkdea=1.05				!phytoplankton death temp correction
        REAL(DP) ksn									!nitrogen half sat constant
        REAL(DP) ksp									!phosphorus half sat constant
        REAL(DP) ksc									!inorganic carbon half sat
        REAL(DP) Isat									!light constant
        REAL(DP) :: khnx=15.0					!Ammonia preference
        REAL(DP) va										!settling velocity of phytoplankton

        !bottom algae
        REAL(DP) kgaF									!max growth rate
        REAL(DP) :: tkgaF=1.066				!max growth temp correction

        !gp 03-Apr-08
        !REAL(DP) kreaF								!botalg respiration rate
        REAL(DP) krea1F								!botalg basal respiration rate
        REAL(DP) krea2F								!botalg photo respiration rate

        REAL(DP) :: tkreaF=1.05				!phytoplankton respiration temp correction
        REAL(DP) kdeaF								!phytoplankton death rate
        REAL(DP) :: tkdeaF=1.05				!phytoplankton death temp correction
        REAL(DP) ksnF									!nitrogen half sat constant
        REAL(DP) kspF									!phosphorus half sat constant
        REAL(DP) kscF									!inorganic carbon half sat
        REAL(DP) abmax								!First-order model carrying capacity
        REAL(DP) kexaF								!Excretion rate
        REAL(DP) :: tkexaF =1.05			!Excretion rate temp correction

        INTEGER(I4B) IlightF					!light model
        REAL(DP) IsatF								!light constant
        REAL(DP) :: khnxF=15.0				!Ammonia preference
        CHARACTER(LEN=30)	:: typeF ='Zero-order'       !Bottom algae growth model,zero-first order
        !POM
        REAL(DP) kdt									!POM dissolution rate
        REAL(DP) :: tkdt=1.05					!POM dissolution temp correction
        REAL(DP) vdt									!POM settling velocity
        !Luxury uptake
        REAL(DP) NINbmin
        REAL(DP) NIPbmin
        REAL(DP) NINbupmax
        REAL(DP) NIPbupmax
        REAL(DP) KqN
        REAL(DP) KqP

        !gp 26-Jan-06
        REAL(DP) NUpWCfrac
        REAL(DP) PUpWCfrac

        !Pathogens
        REAL(DP) kpath								!decay
        REAL(DP) :: tkpath=1.07				!decay temp correction
        REAL(DP) vpath								!settling velocity

        !gp 30-Nov-04 new rates
        REAL(DP) apath								!alpha constant for light mortality of pathogen indicator
        REAL(DP) kgen								!decay of generic constituent
        REAL(DP) :: tkgen=1.07						!decay temp correction for generic constituent
        REAL(DP) vgen								!settling velocity of generic constituent

        !gp 08-Dec-04
        CHARACTER(LEN=30) useGenericAsCOD

        !pH
        REAL(DP) pco2									!partial pressure of carbon dioxide
        REAL(DP) ralkaa, ralkan, ralkbn, ralkbp
        REAL(DP) ralkden, rondn, ralkn
        REAL(DP) rcca, rcco, rccd, rccc

        !gp 03-Nov-04
        CHARACTER(LEN=30) typeH, xdum8
        REAL(DP) :: tkgaH=1.047
        REAL(DP) kgaH, kscH, kinhcH
        INTEGER(I4B) IkoxCH

        !gp 15-Nov-04
        CHARACTER(LEN=30) hco3use, hco3useF

        !gp 15-Nov-04
        REAL(DP) kreaH, kdeaH, ksnH, kspH, khnxH, ahmax
        REAL(DP) :: tkreaH=1.07, tkdeaH=1.07


    END TYPE
!	TYPE(Rates_type) Rates

CONTAINS

    !gp 08-Feb-06
    !FUNCTION Rates_(mgC, mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, Ksodn, Ksop, &
    !				Ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von, kn, &
    !				tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi, kga, tkga, krea, &
    !				tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, va, typeF, kgaF, &
    !				tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, abmax, tkdeaF, ksnF, &
    !				kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, NINbmin, NIPbmin, &
    !				NINbupmax, NIPbupmax, KqN, KqP, kpath, tkpath, vpath, pco2, &
    !				xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, &
    !				kawindmethod, hco3use, hco3useF, typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
    !				kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &		!gp 15-Nov-04
    !				apath, kgen, tkgen, vgen, useGenericAsCOD, &					!gp 08-Dec-04
    !				NUpWCfrac, PUpWCfrac) RESULT(Rates)								!gp 26-Jan-06

    !gp 03-Apr-08
    !FUNCTION Rates_(nr, hydrau, mgC, mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, Ksodn, Ksop, &
    !				Ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von, kn, &
    !				tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi, kga, tkga, krea, &
    !				tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, va, typeF, kgaF, &
    !				tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF, abmax, tkdeaF, ksnF, &
    !				kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, NINbmin, NIPbmin, &
    !				NINbupmax, NIPbupmax, KqN, KqP, kpath, tkpath, vpath, pco2, &
    !				xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, &
    !				kawindmethod, hco3use, hco3useF, typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
    !				kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &
    !				apath, kgen, tkgen, vgen, useGenericAsCOD, &
    !				NUpWCfrac, PUpWCfrac) RESULT(Rates)
    FUNCTION Rates_(nr, hydrau, mgC, mgN, mgP, mgD, mgA, vss,tka, roc, ron, Ksocf, Ksona, Ksodn, Ksop, &
        Ksob, khc, tkhc, kdcs, tkdcs, kdc, tkdc, khn, tkhn, von, kn, &
        tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi, kga, tkga, krea, &
        tkrea, kdea, tkdea, ksn, ksp, ksc, Isat, khnx, va, typeF, kgaF, &
        tkgaF, krea1F, krea2F, tkreaF, kexaF, tkexaF, kdeaF, abmax, tkdeaF, ksnF, &
        kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt, NINbmin, NIPbmin, &
        NINbupmax, NIPbupmax, KqN, KqP, kpath, tkpath, vpath, pco2, &
        xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, &
        kawindmethod, hco3use, hco3useF, typeH, kgaH, tkgaH, kscH, xdum8, kinhcH, &
        kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax, &
        apath, kgen, tkgen, vgen, useGenericAsCOD, &
        NUpWCfrac, PUpWCfrac) RESULT(Rates)


        !gp 08-Feb-06
        !this function also assigns global rates to the unspecified reach-specific rates in hydrau
        USE Class_Hydraulics
        TYPE(RiverHydraulics_type), INTENT(INOUT) :: hydrau			!assign reach-specific rates
        INTEGER(I4B), INTENT(IN) :: nr								!number of reach
        INTEGER(I4B) i

        TYPE(Rates_type) Rates

        !stoichiometry
        REAL(DP) :: mgC, mgN, mgP, mgD, mgA
        REAL(DP), INTENT(IN) :: vss, tka, roc, ron
        REAL(DP), INTENT(IN) :: Ksocf, Ksona, Ksodn, Ksop, Ksob, khc, kdcs, tkdcs, tkhc, kdc, tkdc, khn
        REAL(DP), INTENT(IN) :: tkhn, von, kn, tkn, ki, tki, vdi, tvdi, khp, tkhp, vop, vip, kspi
        REAL(DP), INTENT(IN) :: kga, tkga, krea, tkrea, kdea, tkdea, ksn, ksp, ksc, Isat

        !gp 03-Apr-08
        !REAL(DP), INTENT(IN) :: khnx, va, kgaF, tkgaF, kreaF, tkreaF, kexaF, tkexaF, kdeaF
        REAL(DP), INTENT(IN) :: khnx, va, kgaF, tkgaF, krea1F, krea2F, tkreaF, kexaF, tkexaF, kdeaF

        REAL(DP), INTENT(IN) :: tkdeaF, abmax, ksnF, kspF, kscF, Isatf, khnxF, kdt, tkdt, vdt
        REAL(DP), INTENT(IN) :: NINbmin, NIPbmin, NINbupmax, NIPbupmax, KqN, KqP

        !gp 26-Jan-06
        REAL(DP), INTENT(IN) :: NUpWCfrac, PUpWCfrac

        REAL(DP), INTENT(IN) :: kpath, tkpath, vpath, pco2

        !gp 30-Nov-04
        REAL(DP), INTENT(IN) :: apath, kgen, tkgen, vgen		!gp 30-Nov-04 new paramters for pathogen and generic constituent

        !gp 08-Dec-04
        CHARACTER(LEN=30), INTENT(IN) :: useGenericAsCOD

        !org C, N inhition model, Denitrification enhance model, Algae light model, perihyte light model
        CHARACTER(LEN=30), INTENT(IN) :: xdum1, xdum2, xdum3, xdum4, xdum5, xdum6, xdum7, kai, kawindmethod

        !gp 03-Nov-04
        CHARACTER(LEN=30), INTENT(IN) :: typeH, xdum8
        REAL(DP), INTENT(IN) :: kgaH, tkgaH, kscH, kinhcH

        !gp 15-Nov-04
        CHARACTER(LEN=30), INTENT(IN) :: hco3use, hco3useF

        !gp 15-Nov-04
        REAL(DP), INTENT(IN) :: kreaH, tkreaH, kdeaH, tkdeaH, ksnH, kspH, khnxH, ahmax


        CHARACTER(LEN=30) ::	typeF       !Bottom algae growth model,zero-first order

        REAL(DP) gC, gD

        If (mgC <= 0) mgC = 40.0_dp
        If (mgN <= 0) mgN = 7.2_dp
        If (mgP <= 0) mgP = 1.0_dp
        If (mgD <= 0) mgD = 100.0_dp
        If (mgA <= 0) mgA = 1.0_dp

        Rates%mgA= mgA; Rates%mgD=mgD							!SCC, for v1_3
        !suspended solids
        Rates%vss= vss

        !oxygen
        SELECT CASE (kai)
            !CASE 'Internal'
          CASE ("O'Connor-Dobbins", &
              "Churchill", &
              "Owens-Gibbs", &
              "Thackston-Dawson", &
              "Tsivoglou-Neal", &
              "USGS(pool-riffle)", &
              "USGS(channel-control)")
            Rates%kai = kai
          CASE DEFAULT
            Rates%kai = "Internal"
        END SELECT

        Rates%kawindmethod = kawindmethod
        IF (tka>0) Rates%tka = tka
        IF (roc>0) Rates%roc = roc
        IF (ron>0) Rates%ron = ron

        !Oxygen inhibition of carbon oxidation
        Rates%IkoxC= Ikox(xdum1)
        Rates%Ksocf = Ksocf

        !Oxygen inhibition of nitrification
        Rates%IkoxN= Ikox(xdum2)
        Rates%Ksona = Ksona

        !Oxygen enhancement of denitrification
        Rates%IkoxDN= Ikox(xdum3)
        Rates%Ksodn = Ksodn

        !Oxygen inhibition of phytoplankton respiration
        Rates%IkoxP = Ikox(xdum4)
        Rates%Ksop= Ksop

        !Oxygen inhibition of bottom plant respiration
        Rates%IkoxB = Ikox(xdum5)
        Rates%Ksob= Ksob

        !slow CBOD
        Rates%khc =	khc
        IF (tkhc>0) Rates%tkhc = tkhc
        Rates%kdcs= kdcs
        IF (tkdcs>0) Rates%tkdcs = tkdcs

        !fast CBOD
        Rates%kdc = kdc
        IF (tkdc>0) Rates%tkdc = tkdc

        !organic N
        Rates%khn = khn
        IF (tkhn>0) Rates%tkhn = tkhn
        !organic N settling vel
        Rates%von =von

        !Ammonium
        Rates%kn = kn
        IF (tkn>0) 	Rates%tkn = tkn

        !Nitrate
        Rates%ki = ki
        If (tki > 0) Rates%tki=tki
        Rates%vdi = vdi
        if (tvdi>0)	Rates%tvdi = tvdi

        !organic P
        Rates%khp = khp
        IF (tkhp>0) Rates%tkhp = tkhp

        Rates%vop = vop
        Rates%vip = vip
        Rates%kspi = kspi

        !phytoplankton
        Rates%kga = kga
        IF (tkga>0) Rates%tkga =tkga

        Rates%krea = krea
        IF (tkrea>0) Rates%tkrea = tkrea

        Rates%kdea = kdea
        IF (tkdea>0) Rates%tkdea = tkdea

        Rates%ksn = ksn
        Rates%ksp = ksp
        Rates%ksc = ksc

        SELECT CASE (xdum6)
          CASE ("Smith")
            Rates%Ilight = 2			!"langleys/d"
          CASE ("Steele")
            Rates%Ilight = 3			!"langleys/d"
          CASE DEFAULT								!"Half saturation"
            Rates%Ilight = 1			!"langleys/d"
        END SELECT

        Rates%Isat = Isat
        IF (khnx > 0) Rates%khnx = khnx
        Rates%va = va

        !bottom algae
        Rates%kgaF = kgaF
        IF (tkgaF>0) Rates%tkgaF = tkgaF

        !gp 03-Apr-08
        !Rates%kreaF = kreaF
        Rates%krea1F = krea1F
        Rates%krea2F = krea2F

        IF (tkreaF>0) Rates%tkreaF = tkreaF

        Rates%kexaF = kexaF										!excretion rate
        IF (tkexaF > 0) Rates%tkexaF = tkexaF

        Rates%kdeaF = kdeaF
        IF (tkdeaF>0)  Rates%tkdeaF= tkdeaF
        Rates%ksnF = ksnF
        Rates%kspF = kspF
        Rates%kscF = kscF

        !Rates%abmax= abmax  !SCC 08/09/2004

        IF (typeF/="First-order") THEN
            typeF="Zero-order"
        END IF
        Rates%typeF = typeF

        gC = mgC / 1000.0_dp
        gD = mgD / 1000.0_dp


!gp 03-Apr-08 Starting with version b42a02, these inputs are now per unit dry weight instead of chl a
!	!convert bottom algae rates to a mgA basis
!	If (typeF == "Zero-order") Rates%kgaF = kgaF * gD / mgA     !SCC 08/09/2004
!	Rates%abmax = abmax * gD / mgA										!SCC 08/09/2004
!	Rates%NINbmin = NINbmin * mgA / gD								!SCC 08/09/2004
!	Rates%NIPbmin = NIPbmin * mgA / gD								!SCC 08/09/2004
!	Rates%NINbupmax = NINbupmax * mgA / gD						!SCC 08/09/2004
!	Rates%NIPbupmax = NIPbupmax * mgA / gD						!SCC 08/09/2004
!	Rates%KqN = KqN * mgA / gD												!SCC 08/09/2004
!	Rates%KqP = KqP * mgA / gD												!SCC 08/09/2004
        If (typeF == "Zero-order") Rates%kgaF = kgaF	!gD/m^2/day
        Rates%abmax = abmax								!gD/m^2
        Rates%NINbmin = NINbmin							!mgN/gD
        Rates%NIPbmin = NIPbmin							!mgP/gD
        Rates%NINbupmax = NINbupmax						!mgN/gD/day
        Rates%NIPbupmax = NIPbupmax						!mgP/gD/day
        Rates%KqN = KqN									!mgN/gD
        Rates%KqP = KqP									!mgP/gD

        SELECT CASE (xdum7)
          CASE ("Smith")
            Rates%IlightF = 2			!"langleys/d"
          CASE ("Steele")
            Rates%IlightF = 3			!"langleys/d"
          CASE DEFAULT								!"Half saturation"
            Rates%IlightF = 1			! "langleys/d"
        END SELECT

        Rates%Isatf = Isatf
        If (khnxF > 0) Rates%khnxF = khnxF

        !Rates%NINbmin= NINbmin				!SCC 08/09/2004
        !Rates%NIPbmin= NIPbmin				!SCC 08/09/2004
        !Rates%NINbupmax =NINbupmax		!SCC 08/09/2004
        !Rates%NIPbupmax =NIPbupmax		!SCC 08/09/2004
        !Rates%KqN = KqN							!SCC 08/09/2004
        !Rates%KqP =KqP								!SCC 08/09/2004

        !gp 26-Jan-06
        Rates%NUpWCfrac = NUpWCfrac
        Rates%PUpWCfrac = PUpWCfrac

        !detritus (POM)
        Rates%kdt = kdt
        IF (tkdt>0) Rates%tkdt = tkdt
        Rates%vdt = vdt

        !Pathogens
        Rates%kpath = kpath
        IF (tkpath>0) Rates%tkpath = tkpath
        Rates%vpath = vpath

        !gp 30-Nov-04 new param for pathogens and generic constituent
        Rates%apath = apath
        Rates%kgen = kgen
        IF (tkgen>0) Rates%tkgen = tkgen
        Rates%vgen = vgen

        !gp 08-Dec-04
        Rates%useGenericAsCOD = useGenericAsCOD

        !pH
        Rates%pco2 = pco2				! / 1.0e6

        !gp 15-Nov-04 HCO3- use by phytoplankton and bottom algae
        Rates%hco3use = hco3use
        Rates%hco3useF = hco3useF

        !gp 03-Nov-04 parameters for hyporheic heterotrophic biofilm oxidation of fast CBOD
        Rates%typeH = typeH
        Rates%kgaH = kgaH
        IF (tkgaH>0) Rates%tkgaH = tkgaH
        Rates%kscH = kscH
        Rates%IkoxCH= Ikox(xdum8)
        Rates%kinhcH = kinhcH

        !gp 15-Nov-04 rates for level 2 hyporheic biofilm growth
        Rates%kreaH = kreaH
        IF (tkreaH>0) Rates%tkreaH = tkreaH
        Rates%kdeaH = kdeaH
        IF (tkdeaH>0) Rates%tkdeaH = tkdeaH
        Rates%ksnH = ksnH
        Rates%kspH = kspH
        Rates%khnxH = khnxH
        Rates%ahmax = ahmax

        Rates%anc = mgN / gC
        Rates%apc = mgP / gC
        Rates%ron = ron / 1000.0_dp
        Rates%adc = mgD / mgC
        Rates%roa = roc * gC / mgA
        Rates%ana = mgN / mgA
        Rates%apa = mgP / mgA
        Rates%aca = gC / mgA
        Rates%ada = gD / mgA

        !gp 03-Dec-09
        !!pH/Inorganic carbon stoichiometry
        !Rates%ralkaa = 14.0_dp / 106.0_dp / 12.0_dp * Rates%aca / 1000.0_dp
        !Rates%ralkan = 18.0_dp / 106.0_dp / 12.0_dp * Rates%aca / 1000.0_dp
        !Rates%ralkda = 14.0_dp / 106.0_dp / 12.0_dp / Rates%adc / 1000.0_dp
        !Rates%ralkdn = 18.0_dp / 106.0_dp / 12.0_dp / Rates%adc / 1000.0_dp
        !Rates%ralkn = 2.0_dp / 14.0_dp / 1000.0_dp / 1000.0_dp
        !Rates%ralkden = 4.0_dp / (4.0_dp * 14.0_dp) / 1000.0_dp / 1000.0_dp
        !Rates%racc = 14.0_dp / 106.0_dp / 12.0_dp / 1000.0_dp
        !Rates%rondn = roc * 5.0_dp * 12.0_dp / (4.0_dp * 14.0_dp) / 1000.0_dp
        !Rates%acca = gC / mgA / 12.0_dp / 1000.0_dp
        !Rates%acco = 1.0_dp / 12.0_dp / roc / 1000.0_dp
        !Rates%accd = gC / gD / 12.0_dp / 1000.0_dp
        !Rates%accc = 1.0_dp / 12.0_dp / 1000.0_dp
        !
        !!pH/Inorganic carbon stoichiometry
        !Rates%ralkaa = 14.0_dp / 106.0_dp / 12.0_dp * Rates%aca / 1000.0_dp
        !Rates%ralkan = 18.0_dp / 106.0_dp / 12.0_dp * Rates%aca / 1000.0_dp
        !Rates%ralkbn = 16.0_dp / 16.0_dp / 14.0_dp / 1000.0_dp / 1000.0_dp
        !Rates%ralkbp = 2.0_dp / 1.0_dp / 31.0_dp / 1000.0_dp / 1000.0_dp
        !Rates%ralkn = 2.0_dp / 14.0_dp / 1000.0_dp / 1000.0_dp
        !Rates%ralkden = 4.0_dp / (4.0_dp * 14.0_dp) / 1000.0_dp / 1000.0_dp
        !Rates%rondn = Rates%roc * 5.0_dp * 12.0_dp / (4.0_dp * 14.0_dp) / 1000.0_dp
        !Rates%rcca = gC / mgA / 12.0_dp / 1000.0_dp
        !Rates%rcco = 1.0_dp / 12.0_dp / Rates%roc / 1000.0_dp
        !Rates%rccd = gC / gD / 12.0_dp / 1000.0_dp
        !Rates%rccc = 1.0_dp / 12.0_dp / 1000.0_dp
        Rates%ralkbn = 1.0_dp / 14.0067_dp / 1000.0_dp / 1000.0_dp
        Rates%ralkbp = 1.0_dp / 30.973762_dp / 1000.0_dp / 1000.0_dp	!multiplied later by speciation of phosphate in sub derivs
        Rates%rondn = Rates%roc * 5.0_dp * 12.0107_dp / (4.0_dp * 14.0067_dp) / 1000.0_dp
        Rates%rcca = gC / mgA / 12.0107_dp / 1000.0_dp
        Rates%rcco = 1.0_dp / 12.0107_dp / Rates%roc / 1000.0_dp
        Rates%rccd = gC / gD / 12.0107_dp / 1000.0_dp
        Rates%rccc = 1.0_dp / 12.0107_dp / 1000.0_dp

        !gp 08-Feb-06
        !If reach specific rates provided by input file are invalid, then use default general
        !rates provided by input file
        DO i=1, nr
            !1.reaeration coefficient
            !leave the value as it is
            !unless it is equal/great than zero, it will be calculated by reaeration model
            !2.Inorganic suspended solids settling velocity
            IF (hydrau%reach(i)%vss <0) THEN
                hydrau%reach(i)%vss = Rates%vss
            END IF
            !3.Slow CBOD Hydrolysis rate
            IF (hydrau%reach(i)%khc <0) THEN
                hydrau%reach(i)%khc = Rates%khc
            END IF
            !4.Slow CBOD oxidation rate
            IF (hydrau%reach(i)%kdcs <0) THEN
                hydrau%reach(i)%kdcs = Rates%kdcs
            END IF
            !5.Fast CBOD oxidation rate
            IF (hydrau%reach(i)%kdc < 0) THEN
                hydrau%reach(i)%kdc = Rates%kdc
            END IF
            !6. Organic N hydrolysis rate (khn)
            IF (hydrau%reach(i)%khn < 0) THEN
                hydrau%reach(i)%khn = Rates%khn
            END IF
            !7. Organic N Settling velocity (von)
            IF (hydrau%reach(i)%von < 0) THEN
                hydrau%reach(i)%von = Rates%von
            END IF
            !8. Ammonium Nitrification rate (kn)
            IF (hydrau%reach(i)%kn < 0) THEN
                hydrau%reach(i)%kn = Rates%kn
            END IF
            !9. Nitrate Denitrification rate (ki)
            IF (hydrau%reach(i)%ki < 0) THEN
                hydrau%reach(i)%ki = Rates%ki
            END IF
            !10. Nitrate Sed denitrification transfer coeff (vdi)
            IF (hydrau%reach(i)%vdi < 0) THEN
                hydrau%reach(i)%vdi = Rates%vdi
            END IF
            !11. Organic P Hydrolysis (khp)
            IF (hydrau%reach(i)%khp < 0) THEN
                hydrau%reach(i)%khp = Rates%khp
            END IF
            !12. Organic P Settling velocity (vop)
            IF (hydrau%reach(i)%vop < 0) THEN
                hydrau%reach(i)%vop = Rates%vop
            END IF
            !13. Inorganic P settling veloctiy (vip)
            IF (hydrau%reach(i)%vip < 0) THEN
                hydrau%reach(i)%vip = Rates%vip
            END IF
            !14. Phytoplankton Max growth rate (kga)
            IF (hydrau%reach(i)%kga < 0) THEN
                hydrau%reach(i)%kga = Rates%kga
            END IF
            !Phytoplankton respiration rate (krea)
            IF (hydrau%reach(i)%krea < 0) THEN
                hydrau%reach(i)%krea = Rates%krea
            END IF
            !Phytoplankton death rate (kdea)
            IF (hydrau%reach(i)%kdea < 0) THEN
                hydrau%reach(i)%kdea = Rates%kdea
            END IF
            !Phytoplankton N half-sat (ksn)
            IF (hydrau%reach(i)%ksn < 0) THEN
                hydrau%reach(i)%ksn = Rates%ksn
            END IF
            !Phytoplankton P half-sat (ksp)
            IF (hydrau%reach(i)%ksp < 0) THEN
                hydrau%reach(i)%ksp = Rates%ksp
            END IF
            !Phytoplankton light sat (Isat)
            IF (hydrau%reach(i)%Isat < 0) THEN
                hydrau%reach(i)%Isat = Rates%Isat
            END IF
            !Phytoplankton ammonia pref (khnx)
            IF (hydrau%reach(i)%khnx < 0) THEN
                hydrau%reach(i)%khnx = Rates%khnx
            END IF
            !Phytoplankton settling (va)
            IF (hydrau%reach(i)%va < 0) THEN
                hydrau%reach(i)%va = Rates%va
            END IF

            !Bottom plant initial biomass (botalg0)
            !leave the value as it is
            !Bottom plant Max growth rate (kgaF)
            IF (hydrau%reach(i)%kgaF < 0) THEN
                hydrau%reach(i)%kgaF = Rates%kgaF

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	!convert to mgA basis
                !	If (Rates%typeF == "Zero-order") hydrau%reach(i)%kgaF = hydrau%reach(i)%kgaF * gD / mgA

            END IF

            !Bottom plant first order carrying capacity (abmax)
            IF (hydrau%reach(i)%abmax < 0) THEN
                hydrau%reach(i)%abmax = Rates%abmax

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%abmax = hydrau%reach(i)%abmax * gD / mgA		!convert to mgA basis

            END IF

            !Bottom plant respiration rate (kreaF)

            !gp 03-Apr-08
            !IF (hydrau%reach(i)%kreaF < 0) THEN
            !	hydrau%reach(i)%kreaF = Rates%kreaF
            !END IF
            IF (hydrau%reach(i)%krea1F < 0) THEN
                hydrau%reach(i)%krea1F = Rates%krea1F
            END IF
            IF (hydrau%reach(i)%krea2F < 0) THEN
                hydrau%reach(i)%krea2F = Rates%krea2F
            END IF

            !Bottom plant excretion rate (kexaF)
            IF (hydrau%reach(i)%kexaF < 0) THEN
                hydrau%reach(i)%kexaF = Rates%kexaF
            END IF
            !Bottom plant death rate (kdeaF)
            IF (hydrau%reach(i)%kdeaF < 0) THEN
                hydrau%reach(i)%kdeaF = Rates%kdeaF
            END IF
            !Bottom plant external N half-sat (ksnF)
            IF (hydrau%reach(i)%ksnF < 0) THEN
                hydrau%reach(i)%ksnF = Rates%ksnF
            END IF
            !Bottom plant external P half-sat (kspF)
            IF (hydrau%reach(i)%kspF < 0) THEN
                hydrau%reach(i)%kspF = Rates%kspF
            END IF
            !Bottom plant light sat (IsatF)
            IF (hydrau%reach(i)%IsatF < 0) THEN
                hydrau%reach(i)%IsatF = Rates%IsatF
            END IF
            !Bottom plant NH4 pref (khnxF)
            IF (hydrau%reach(i)%khnxF < 0) THEN
                hydrau%reach(i)%khnxF = Rates%khnxF
            END IF
            !Bottom plant subsistence quota for N (NINbmin)
            IF (hydrau%reach(i)%NINbmin < 0) THEN
                hydrau%reach(i)%NINbmin = Rates%NINbmin

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%NINbmin = hydrau%reach(i)%NINbmin * mgA / gD		!convert to mgA basis

            END IF
            !Bottom plant subsistence quota for P (NIPbmin)
            IF (hydrau%reach(i)%NIPbmin < 0) THEN
                hydrau%reach(i)%NIPbmin = Rates%NIPbmin

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%NIPbmin = hydrau%reach(i)%NIPbmin * mgA / gD		!convert to mgA basis

            END IF
            !Bottom plant max uptake of N (NINbupmax)
            IF (hydrau%reach(i)%NINbupmax < 0) THEN
                hydrau%reach(i)%NINbupmax = Rates%NINbupmax

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%NINbupmax = hydrau%reach(i)%NINbupmax * mgA / gD	!convert to mgA basis

            END IF
            !Bottom plant max uptake of P (NIPbupmax)
            IF (hydrau%reach(i)%NIPbupmax < 0) THEN
                hydrau%reach(i)%NIPbupmax = Rates%NIPbupmax

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%NIPbupmax = hydrau%reach(i)%NIPbupmax * mgA / gD	!convert to mgA basis

            END IF
            !Bottom plant internal N half-sat (KqN)
            IF (hydrau%reach(i)%KqN < 0) THEN
                hydrau%reach(i)%KqN = Rates%KqN

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%KqN = hydrau%reach(i)%KqN * mgA / gD				!convert to mgA basis

            END IF
            !Bottom plant internal P half-sat (KqP)
            IF (hydrau%reach(i)%KqP < 0) THEN
                hydrau%reach(i)%KqP = Rates%KqP

                !gp 03-Apr-08 already in gD basis
                !ELSE
                !	hydrau%reach(i)%KqP = hydrau%reach(i)%KqP * mgA / gD				!convert to mgA basis

            END IF
            !Bottom plant fraction of N uptake from water column (NUpWCfrac)
            IF (hydrau%reach(i)%NUpWCfrac < 0) THEN
                hydrau%reach(i)%NUpWCfrac = Rates%NUpWCfrac
            END IF
            !Bottom plant fraction of P uptake from water column (PUpWCfrac)
            IF (hydrau%reach(i)%PUpWCfrac < 0) THEN
                hydrau%reach(i)%PUpWCfrac = Rates%PUpWCfrac
            END IF
            !17. Detritus(POM) dissolution rate (kdt)
            IF (hydrau%reach(i)%kdt < 0) THEN
                hydrau%reach(i)%kdt = Rates%kdt
            END IF
            !18. Detritus(POM) Settling velocity (vdt)
            IF (hydrau%reach(i)%vdt < 0) THEN
                hydrau%reach(i)%vdt = Rates%vdt
            END IF
            !Pathogen dieoff rate (kpath)
            IF (hydrau%reach(i)%kpath < 0) THEN
                hydrau%reach(i)%kpath = Rates%kpath
            END IF
            !Pathogen settling rate (vpath)
            IF (hydrau%reach(i)%vpath < 0) THEN
                hydrau%reach(i)%vpath = Rates%vpath
            END IF
            !Pathogen light alpha (apath)
            IF (hydrau%reach(i)%apath < 0) THEN
                hydrau%reach(i)%apath = Rates%apath
            END IF
            !Hyporheic heterotroph max growth rate (kgaH)
            IF (hydrau%reach(i)%kgaH < 0) THEN
                hydrau%reach(i)%kgaH = Rates%kgaH
            END IF
            !Hyporheic heterotroph half-sat for CBOD (kscH)
            IF (hydrau%reach(i)%kscH < 0) THEN
                hydrau%reach(i)%kscH = Rates%kscH
            END IF
            !Hyporheic heterotroph O2 inhibition parameter (kinhcH)
            IF (hydrau%reach(i)%kinhcH < 0) THEN
                hydrau%reach(i)%kinhcH = Rates%kinhcH
            END IF
            !Hyporheic heterotroph respiration rate (kreaH)
            IF (hydrau%reach(i)%kreaH < 0) THEN
                hydrau%reach(i)%kreaH = Rates%kreaH
            END IF
            !Hyporheic heterotroph death rate (kdeaH)
            IF (hydrau%reach(i)%kdeaH < 0) THEN
                hydrau%reach(i)%kdeaH = Rates%kdeaH
            END IF
            !Hyporheic heterotroph N half-sat (ksnH)
            IF (hydrau%reach(i)%ksnH < 0) THEN
                hydrau%reach(i)%ksnH = Rates%ksnH
            END IF
            !Hyporheic heterotroph P half-sat (kspH)
            IF (hydrau%reach(i)%kspH < 0) THEN
                hydrau%reach(i)%kspH = Rates%kspH
            END IF
            !Hyporheic heterotroph NH4 preference (khnxH)
            IF (hydrau%reach(i)%khnxH < 0) THEN
                hydrau%reach(i)%khnxH = Rates%khnxH
            END IF
            !Hyporheic heterotroph first-order carrying capacity (ahmax)
            IF (hydrau%reach(i)%ahmax < 0) THEN
                hydrau%reach(i)%ahmax = Rates%ahmax
            END IF
            !Generic constituent first-order loss rate constant (kgen)
            IF (hydrau%reach(i)%kgen < 0) THEN
                hydrau%reach(i)%kgen = Rates%kgen
            END IF
            !Generic constituent settling rate (vgen)
            IF (hydrau%reach(i)%vgen < 0) THEN
                hydrau%reach(i)%vgen = Rates%vgen
            END IF
        END DO

    END FUNCTION

    !gp 03-Apr-08
    !SUBROUTINE TempAdjust(Rates, hydrau, siteMeteo, nr, Te, khcT, kdcsT, kdcT, kgaFT, kdeaFT, &
    !										kreaFT, kexaFT, khnT, khpT, knT, kdtT, kpathT, kiT, kaT, kgaT, kdeaT, &
    !										kreaT, vdiT, kacT, &
    !										kgenT)		!gp 30-Nov-04 add kgenT for generic constituent
    SUBROUTINE TempAdjust(Rates, hydrau, siteMeteo, nr, Te, khcT, kdcsT, kdcT, kgaFT, kdeaFT, &
        krea1FT, krea2FT, kexaFT, khnT, khpT, knT, kdtT, kpathT, kiT, kaT, kgaT, kdeaT, &
        kreaT, vdiT, kacT, &
        kgenT)
        USE Class_Hydraulics
        USE m_meteorology

        INTEGER(I4B), INTENT(IN) :: nr
        TYPE(Rates_type) Rates
        TYPE(RiverHydraulics_type) hydrau
        TYPE(t_meteorology) siteMeteo

        REAL(DP), DIMENSION(0:,:), INTENT(IN) :: Te

        !gp 03-Apr-08
        !REAL(DP), DIMENSION(:), INTENT(OUT) :: khcT, kdcsT, kdcT, kgaFT, kdeaFT, kreaFT, kexaFT
        REAL(DP), DIMENSION(:), INTENT(OUT) :: khcT, kdcsT, kdcT, kgaFT, kdeaFT, krea1FT, krea2FT, kexaFT

        REAL(DP), DIMENSION(:), INTENT(OUT) :: khnT, khpT, knT, kdtT, kpathT, kiT, kaT, kgaT
        REAL(DP), DIMENSION(:), INTENT(OUT) :: kdeaT, kreaT, vdiT, kacT

        !gp 30-Nov-04
        REAL(DP), DIMENSION(:), INTENT(OUT) :: kgenT	!gp 30-Nov-04

        INTEGER(I4B) i
        REAL(DP) kawind, Uw10
!		REAL(DP), DIMENSION(nr) :: ka

        DO i = 1, nr

            !gp 08-Feb-06
            !khcT(i) = Rates%khc * Rates%tkhc ** (Te(i, 1) - 20.0_dp)
            !kdcsT(i) = Rates%kdcs * Rates%tkdcs **(Te(i,1) - 20.0_dp)
            !kdcT(i) = Rates%kdc * Rates%tkdc ** (Te(i, 1) - 20.0_dp)
            !khnT(i) = Rates%khn * Rates%tkhn ** (Te(i, 1) - 20.0_dp)
            !khpT(i) = Rates%khp * Rates%tkhp ** (Te(i, 1) - 20.0_dp)
            !knT(i) = Rates%kn * Rates%tkn ** (Te(i, 1) - 20.0_dp)
            !kiT(i) = Rates%ki * Rates%tki ** (Te(i, 1) - 20.0_dp)
            !vdiT(i) = Rates%vdi * Rates%tvdi ** (Te(i, 1) - 20.0_dp)
            !
            !!use interpolated wind speed to adjust reaeration for selected methods
            !If (hydrau%reach(i)%kaf == "Specified") Then
            !	hydrau%reach(i)%ka = hydrau%reach(i)%kau
            !Else
            !	Uw10 = (10.0_dp / 7.0_dp) ** 0.15_dp * siteMeteo%Uw(i)
            !	Select Case (Rates%kawindmethod)
            !		Case ("Banks-Herrera")     !Chapra (1997) eqn 20.46
            !			kawind = 0.728_dp * Uw10 ** 0.5_dp - 0.317_dp * Uw10 + 0.0372_dp * Uw10 ** 2
            !		Case ("Wanninkhof")        !Chapra (1997) eqn 20.47
            !			kawind = 0.0986_dp * Uw10 ** 1.64_dp
            !		Case DEFAULT
            !			kawind = 0
            !	End Select
            !	hydrau%reach(i)%ka = hydrau%reach(i)%kau + kawind / hydrau%reach(i)%depth
            !End If
            !kaT(i) = hydrau%reach(i)%ka * Rates%tka ** (Te(i, 1) - 20.0_dp)
            !kacT(i) = (32.0_dp / 44.0_dp) ** 0.25_dp * kaT(i)
!			!channel(i)%elev = (channel(i)%elev1 + channel(i)%elev2) / 2
            !hydrau%reach(i)%os = oxsat(Te(i, 1), hydrau%reach(i)%elev)
            !kgaT(i) = Rates%kga * Rates%tkga ** (Te(i, 1) - 20.0_dp)
            !kdeaT(i) = Rates%kdea * Rates%tkdea ** (Te(i, 1) - 20.0_dp)
            !kreaT(i) = Rates%krea * Rates%tkrea ** (Te(i, 1) - 20.0_dp)
            !kgaFT(i) = Rates%kgaF * Rates%tkgaF ** (Te(i, 1) - 20.0_dp)
            !kdeaFT(i) = Rates%kdeaF * Rates%tkdeaF ** (Te(i, 1) - 20.0_dp)
            !kreaFT(i) = Rates%kreaF * Rates%tkreaF ** (Te(i, 1) - 20.0_dp)
            !kexaFT(i) = Rates%kexaF * Rates%tkexaF ** (Te(i, 1) - 20.0_dp)
            !kdtT(i) = Rates%kdt * Rates%tkdt ** (Te(i, 1) - 20.0_dp)
            !kpathT(i) = Rates%kpath * Rates%tkpath ** (Te(i, 1) - 20.0_dp)
            !
            !!GP 30-nOV-04
            !kgenT(i) = Rates%kgen * Rates%tkgen ** (Te(i, 1) - 20.0_dp)
            khcT(i) = hydrau%reach(i)%khc * Rates%tkhc ** (Te(i, 1) - 20.0_dp)
            kdcsT(i) = hydrau%reach(i)%kdcs * Rates%tkdcs **(Te(i,1) - 20.0_dp)
            kdcT(i) = hydrau%reach(i)%kdc * Rates%tkdc ** (Te(i, 1) - 20.0_dp)
            khnT(i) = hydrau%reach(i)%khn * Rates%tkhn ** (Te(i, 1) - 20.0_dp)
            khpT(i) = hydrau%reach(i)%khp * Rates%tkhp ** (Te(i, 1) - 20.0_dp)
            knT(i) = hydrau%reach(i)%kn * Rates%tkn ** (Te(i, 1) - 20.0_dp)
            kiT(i) = hydrau%reach(i)%ki * Rates%tki ** (Te(i, 1) - 20.0_dp)
            vdiT(i) = hydrau%reach(i)%vdi * Rates%tvdi ** (Te(i, 1) - 20.0_dp)
            !use interpolated wind speed to adjust reaeration for selected methods
            If (hydrau%reach(i)%kaf == "Specified") Then
                hydrau%reach(i)%ka = hydrau%reach(i)%kau
            Else
                Uw10 = (10.0_dp / 7.0_dp) ** 0.15_dp * siteMeteo%Uw(i)
                Select Case (Rates%kawindmethod)
                  Case ("Banks-Herrera")     !Chapra (1997) eqn 20.46
                    kawind = 0.728_dp * Uw10 ** 0.5_dp - 0.317_dp * Uw10 + 0.0372_dp * Uw10 ** 2
                  Case ("Wanninkhof")        !Chapra (1997) eqn 20.47
                    kawind = 0.0986_dp * Uw10 ** 1.64_dp
                  Case DEFAULT
                    kawind = 0
                End Select
                hydrau%reach(i)%ka = hydrau%reach(i)%kau + kawind / hydrau%reach(i)%depth
            End If
            kaT(i) = hydrau%reach(i)%ka * Rates%tka ** (Te(i, 1) - 20.0_dp)
            kacT(i) = (32.0_dp / 44.0_dp) ** 0.25_dp * kaT(i)
            hydrau%reach(i)%os = oxsat(Te(i, 1), hydrau%reach(i)%elev)
            kgaT(i) = hydrau%reach(i)%kga * Rates%tkga ** (Te(i, 1) - 20.0_dp)
            kdeaT(i) = hydrau%reach(i)%kdea * Rates%tkdea ** (Te(i, 1) - 20.0_dp)
            kreaT(i) = hydrau%reach(i)%krea * Rates%tkrea ** (Te(i, 1) - 20.0_dp)
            kgaFT(i) = hydrau%reach(i)%kgaF * Rates%tkgaF ** (Te(i, 1) - 20.0_dp)
            kdeaFT(i) = hydrau%reach(i)%kdeaF * Rates%tkdeaF ** (Te(i, 1) - 20.0_dp)

            !gp 03-Apr-08
            !kreaFT(i) = hydrau%reach(i)%kreaF * Rates%tkreaF ** (Te(i, 1) - 20.0_dp)
            krea1FT(i) = hydrau%reach(i)%krea1F * Rates%tkreaF ** (Te(i, 1) - 20.0_dp)
            krea2FT(i) = hydrau%reach(i)%krea2F		!no temp adj of photo resp because growth rate is temp adjusted

            kexaFT(i) = hydrau%reach(i)%kexaF * Rates%tkexaF ** (Te(i, 1) - 20.0_dp)
            kdtT(i) = hydrau%reach(i)%kdt * Rates%tkdt ** (Te(i, 1) - 20.0_dp)
            kpathT(i) = hydrau%reach(i)%kpath * Rates%tkpath ** (Te(i, 1) - 20.0_dp)
            kgenT(i) = hydrau%reach(i)%kgen * Rates%tkgen ** (Te(i, 1) - 20.0_dp)

        END DO
    END SUBROUTINE TempAdjust

    PURE FUNCTION oxsat(Temp, elev)

        REAL(DP) oxsat
        REAL(DP), INTENT(IN) :: Temp, elev
        REAL(DP) Taa, lnosf

        Taa = Temp + 273.15_dp
        lnosf = -139.34411_dp + 157570.1_dp / Taa - 66423080.0_dp / Taa ** 2 &
            + 12438000000.0_dp / Taa ** 3 - 862194900000.0_dp / Taa ** 4
        oxsat = Exp(lnosf) * (1 - 0.0001148_dp * elev)

    END FUNCTION oxsat

    SUBROUTINE SetOxygenInhibEnhance(Rates, O2, fcarb, fnitr, fdenitr, frespp, frespb)

        REAL(DP), INTENT(OUT) :: fcarb, fnitr, fdenitr, frespp, frespb
        TYPE(Rates_type), INTENT(IN) :: Rates
        REAL(DP), INTENT(IN) :: O2

        Select Case (Rates%IkoxC)
          Case (1)
            fcarb = O2 / (Rates%Ksocf + O2)
          Case (2)
            fcarb = 1.0_dp - Exp(-Rates%Ksocf * O2)
          Case (3)
            fcarb = O2 ** 2 / (Rates%Ksocf + O2 ** 2)
        End Select

        Select Case (Rates%IkoxN)
          Case (1)
            fnitr = O2 / (Rates%Ksona + O2)
          Case (2)
            fnitr = 1.0_dp - Exp(-Rates%Ksona * O2)
          Case (3)
            fnitr = O2 ** 2 / (Rates%Ksona + O2 ** 2)
        End Select

        Select Case (Rates%IkoxDN)
          Case (1)
            fdenitr = 1.0_dp - O2 / (Rates%Ksodn + O2)
          Case (2)
            fdenitr = Exp(-Rates%Ksodn * O2)
          Case (3)
            fdenitr = 1.0_dp - O2 ** 2 / (Rates%Ksodn + O2 ** 2)
        End Select

        Select Case (Rates%IkoxP)
          Case (1)
            frespp = O2 / (Rates%Ksop + O2)
          Case (2)
            frespp = 1.0_dp - Exp(-Rates%Ksop * O2)
          Case (3)
            frespp = O2 ** 2 / (Rates%Ksop + O2 ** 2)
        End Select

        Select Case (Rates%IkoxB)
          Case (1)
            frespb = O2 / (Rates%Ksob + O2)
          Case (2)
            frespb = 1.0_dp - Exp(-Rates%Ksob * O2)
          Case (3)
            frespb = O2 ** 2 / (Rates%Ksob + O2 ** 2)
        End Select

    END SUBROUTINE SetOxygenInhibEnhance

    PURE FUNCTION Ikox(xdum)
        CHARACTER(LEN=30), INTENT(IN) :: xdum
        INTEGER(I4B) Ikox

        SELECT CASE (xdum)
          CASE ("Half saturation")
            Ikox = 1									!"mgO2/L"
          CASE ("2nd order")
            Ikox = 3									!"mgO2/L"
          CASE DEFAULT								!"Exponential"
            Ikox = 2									!"L/mgO2"  End
        END SELECT
    END FUNCTION
END MODULE Class_Rates
