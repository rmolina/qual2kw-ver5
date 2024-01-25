MODULE Class_IntegrationData

	USE nrtype

	!gp 01-Nov-04
	IMPLICIT NONE

	TYPE Integral_type
!		INTEGER(I4B) np, nc													!days, step in each day
		!dependent variables
		REAL(DP), DIMENSION(:), POINTER :: INb, IPb
		!gp 27-Oct-04 REAL(DP), DIMENSION(:,:), POINTER :: Te, c
		REAL(DP), DIMENSION(:,:), POINTER :: Te		!gp
		REAL(DP), DIMENSION(:,:,:), POINTER :: c	!gp add dimension for nl
	END TYPE	Integral_type

	!sediment nutrient flux
	REAL(DP), DIMENSION(:), POINTER :: SODpr, JNH4pr, JNO3pr, JCH4pr, JSRPpr

	!gp REAL(DP), DIMENSION(:), POINTER :: os
	REAL(DP), DIMENSION(:,:), POINTER :: os		!gp
	
	!gp 20-Oct-04 bottom algae limitation factors (fraction of maximum potential growth)
	REAL(DP), DIMENSION(:), POINTER :: phitotalSave, phitSave, philSave, phinSave, phipSave, phicSave	!gp 20-Oct-04

	!gp 28-Oct-04 diagenesis and hyporheic sediment/water fluxes (positive is source to water from sediment)
	REAL(DP), DIMENSION(:), POINTER :: DiagFluxDO, DiagFluxCBOD, DiagFluxNH4, DiagFluxNO3, DiagFluxSRP, DiagFluxIC	!gp 28-Oct-04
	REAL(DP), DIMENSION(:), POINTER :: HypoFluxDO, HypoFluxCBOD, HypoFluxNH4, HypoFluxNO3, HypoFluxSRP, HypoFluxIC	!gp 28-Oct-04

	!gp 29-Oct-04 
	REAL(DP), DIMENSION(:), POINTER :: CSODpr

	!gp 05-Jul-05 heat/DO/CO2 fluxes
	REAL(DP), DIMENSION(:), POINTER :: saveHeatFluxJsnt, saveHeatFluxLongat, saveHeatFluxBack, saveHeatFluxConv
	REAL(DP), DIMENSION(:), POINTER :: saveHeatFluxEvap, saveHeatFluxJsed, saveHeatFluxJhyporheic, saveHeatFluxTribs
	REAL(DP), DIMENSION(:), POINTER :: saveHeatFluxAdvecDisp
	REAL(DP), DIMENSION(:), POINTER :: saveDOfluxReaer, saveDOfluxCBODfast, saveDOfluxCBODslow, saveDOfluxNitrif
	REAL(DP), DIMENSION(:), POINTER :: saveDOfluxPhytoResp, saveDOfluxPhytoPhoto, saveDOfluxBotalgResp, saveDOfluxBotalgPhoto
	REAL(DP), DIMENSION(:), POINTER :: saveDOfluxSOD, saveDOfluxCOD, saveDOfluxHyporheic, saveDOfluxTribs
	REAL(DP), DIMENSION(:), POINTER :: saveDOfluxAdvecDisp, saveDOfluxHeadwater
	REAL(DP), DIMENSION(:), POINTER :: saveCO2fluxReaer, saveCO2fluxCBODfast, saveCO2fluxCBODslow
	REAL(DP), DIMENSION(:), POINTER :: saveCO2fluxPhytoResp, saveCO2fluxPhytoPhoto, saveCO2fluxBotalgResp, saveCO2fluxBotalgPhoto
	REAL(DP), DIMENSION(:), POINTER :: saveCO2fluxSOD, saveCO2fluxHyporheic, saveCO2fluxTribs
	REAL(DP), DIMENSION(:), POINTER :: saveCO2fluxAdvecDisp, saveCO2fluxHeadwater

	!gp 25-Jun-09
	REAL(DP), DIMENSION(:), POINTER :: saveBotAlgPhoto, saveBotAlgResp, saveBotAlgDeath, saveBotAlgNetGrowth
	!REAL(DP), DIMENSION(:), POINTER :: av_BotAlgPhoto, av_BotAlgResp, av_BotAlgDeath, av_BotAlgNetGrowth

	CONTAINS

	!user-defined type 'Integral_type' constructor
	FUNCTION Integration_(nr) RESULT(integral)

!		USE Class_RiverTopo
		TYPE(Integral_type) integral
		INTEGER(I4B), INTENT(IN) :: nr

		!gp 05-Jul-05 INTEGER(I4B) i, status(29)	!gp 28-Oct-04
		!gp 25-Jun-09 INTEGER(I4B) i, status(64)	!gp 05-Jul-05
		INTEGER(I4B) i, status(68)	!gp 25-Jun-09

		!gp ALLOCATE(os(0:nr), STAT=status(1))
		ALLOCATE(os(0:nr, nl), STAT=status(1))	!gp
		
		ALLOCATE(integral%INb(0:nr), STAT=status(2))
		ALLOCATE(integral%IPb(0:nr), STAT=status(3))

		!gp ALLOCATE(integral%c(0:nr+1, nv), STAT=status(4))
		ALLOCATE(integral%c(0:nr+1, nv, nl), STAT=status(4))

		ALLOCATE(integral%Te(0:nr+1, nl), STAT=status(5))
		ALLOCATE(JNH4pr(0:nr), STAT=status(6))
		ALLOCATE(JNO3pr(0:nr), STAT=status(7))
		ALLOCATE(JCH4pr(0:nr), STAT=status(8))
		ALLOCATE(JSRPpr(0:nr), STAT=status(9))
		ALLOCATE(SODpr(0:nr), STAT=status(10))

		!gp 20-Oct-04 growth limitation factors for bottom algae
		ALLOCATE(phitotalSave(0:nr), STAT=status(11))	
		ALLOCATE(phitSave(0:nr), STAT=status(12))
		ALLOCATE(philSave(0:nr), STAT=status(13))	
		ALLOCATE(phinSave(0:nr), STAT=status(14))
		ALLOCATE(phipSave(0:nr), STAT=status(15))
		ALLOCATE(phicSave(0:nr), STAT=status(16))	!gp end new block

		!gp 28-Oct-04 Diagenesis fluxes between sediment/water
		ALLOCATE(DiagFluxDO(0:nr), STAT=status(17))
		ALLOCATE(DiagFluxCBOD(0:nr), STAT=status(18))
		ALLOCATE(DiagFluxNH4(0:nr), STAT=status(19))
		ALLOCATE(DiagFluxNO3(0:nr), STAT=status(20))
		ALLOCATE(DiagFluxSRP(0:nr), STAT=status(21))
		ALLOCATE(DiagFluxIC(0:nr), STAT=status(22))	!gp end new block
		
		!gp 28-Oct-04 Hyporheic exchange fluxes between sediment/water
		ALLOCATE(HypoFluxDO(0:nr), STAT=status(23))
		ALLOCATE(HypoFluxCBOD(0:nr), STAT=status(24))
		ALLOCATE(HypoFluxNH4(0:nr), STAT=status(25))
		ALLOCATE(HypoFluxNO3(0:nr), STAT=status(26))
		ALLOCATE(HypoFluxSRP(0:nr), STAT=status(27))
		ALLOCATE(HypoFluxIC(0:nr), STAT=status(28))	!gp end new block

		!gp 29-Oct-04
		ALLOCATE(CSODpr(0:nr), STAT=status(29))

		!gp 05-Jul-05 heat/DO/CO2 fluxes
		ALLOCATE(saveHeatFluxJsnt(0:nr), STAT=status(30))	
		ALLOCATE(saveHeatFluxLongat(0:nr), STAT=status(31))	
		ALLOCATE(saveHeatFluxBack(0:nr), STAT=status(32))	
		ALLOCATE(saveHeatFluxConv(0:nr), STAT=status(33))	
		ALLOCATE(saveHeatFluxEvap(0:nr), STAT=status(34))	
		ALLOCATE(saveHeatFluxJsed(0:nr), STAT=status(35))	
		ALLOCATE(saveHeatFluxJhyporheic(0:nr), STAT=status(36))	
		ALLOCATE(saveHeatFluxTribs(0:nr), STAT=status(37))	
		ALLOCATE(saveHeatFluxAdvecDisp(0:nr), STAT=status(38))	
		ALLOCATE(saveDOfluxReaer(0:nr), STAT=status(39))	
		ALLOCATE(saveDOfluxCBODfast(0:nr), STAT=status(40))	
		ALLOCATE(saveDOfluxCBODslow(0:nr), STAT=status(41))	
		ALLOCATE(saveDOfluxCOD(0:nr), STAT=status(42))	
		ALLOCATE(saveDOfluxNitrif(0:nr), STAT=status(43))	
		ALLOCATE(saveDOfluxPhytoResp(0:nr), STAT=status(44))	
		ALLOCATE(saveDOfluxPhytoPhoto(0:nr), STAT=status(45))	
		ALLOCATE(saveDOfluxBotalgResp(0:nr), STAT=status(46))	
		ALLOCATE(saveDOfluxBotalgPhoto(0:nr), STAT=status(47))	
		ALLOCATE(saveDOfluxSOD(0:nr), STAT=status(48))	
		ALLOCATE(saveDOfluxHyporheic(0:nr), STAT=status(49))	
		ALLOCATE(saveDOfluxTribs(0:nr), STAT=status(50))	
		ALLOCATE(saveDOfluxAdvecDisp(0:nr), STAT=status(51))	
		ALLOCATE(saveDOfluxHeadwater(0:nr), STAT=status(52))	
		ALLOCATE(saveCO2fluxReaer(0:nr), STAT=status(53))	
		ALLOCATE(saveCO2fluxCBODfast(0:nr), STAT=status(54))	
		ALLOCATE(saveCO2fluxCBODslow(0:nr), STAT=status(55))	
		ALLOCATE(saveCO2fluxPhytoResp(0:nr), STAT=status(56))	
		ALLOCATE(saveCO2fluxPhytoPhoto(0:nr), STAT=status(57))	
		ALLOCATE(saveCO2fluxBotalgResp(0:nr), STAT=status(58))	
		ALLOCATE(saveCO2fluxBotalgPhoto(0:nr), STAT=status(59))	
		ALLOCATE(saveCO2fluxSOD(0:nr), STAT=status(60))	
		ALLOCATE(saveCO2fluxHyporheic(0:nr), STAT=status(61))	
		ALLOCATE(saveCO2fluxTribs(0:nr), STAT=status(62))	
		ALLOCATE(saveCO2fluxAdvecDisp(0:nr), STAT=status(63))	
		ALLOCATE(saveCO2fluxHeadwater(0:nr), STAT=status(64))	

		!gp 25-Jun-09
		ALLOCATE(saveBotAlgPhoto(0:nr), STAT=status(65))
		ALLOCATE(saveBotAlgResp(0:nr), STAT=status(66))
		ALLOCATE(saveBotAlgDeath(0:nr), STAT=status(67))
		ALLOCATE(saveBotAlgNetGrowth(0:nr), STAT=status(68))
		!ALLOCATE(av_BotAlgPhoto(0:nr), STAT=status(69))
		!ALLOCATE(av_BotAlgResp(0:nr), STAT=status(70))
		!ALLOCATE(av_BotAlgDeath(0:nr), STAT=status(71))
		!ALLOCATE(av_BotAlgNetGrowth(0:nr), STAT=status(72))

		!initialize
		os=0;			
		integral%INb=0;	integral%IPb=0;		integral%c=0;		integral%Te=0
		SODpr=0;	JNH4pr=0;	JNO3pr=0;	JCH4pr=0;	JSRPpr=0

		!gp 20-Oct-04
		phitotalSave=0; phitSave=0; philSave=0; phinSave=0; phipSave=0; phicSave=0	!gp end new block

		!gp 25-Jun-09
		saveBotAlgPhoto=0; saveBotAlgResp=0; saveBotAlgDeath=0; saveBotAlgNetGrowth=0
		!av_BotAlgPhoto=0; av_BotAlgResp=0; av_BotAlgDeath=0; av_BotAlgNetGrowth=0

		!28-Oct-04 
		DiagFluxDO=0; DiagFluxCBOD=0; DiagFluxNH4=0; DiagFluxNO3=0; DiagFluxSRP=0; DiagFluxIC=0
		HypoFluxDO=0; HypoFluxCBOD=0; HypoFluxNH4=0; HypoFluxNO3=0; HypoFluxSRP=0; HypoFluxIC=0	!gp end new block

		!gp 29-Oct-04
		CSODpr=0

		!gp 05-Jul-05
		saveHeatFluxJsnt=0; saveHeatFluxLongat=0; saveHeatFluxBack=0; saveHeatFluxConv=0
		saveHeatFluxEvap=0; saveHeatFluxJsed=0; saveHeatFluxJhyporheic=0; saveHeatFluxTribs=0
		saveHeatFluxAdvecDisp=0
		saveDOfluxReaer=0; saveDOfluxCBODfast=0; saveDOfluxCBODslow=0; saveDOfluxNitrif=0
		saveDOfluxPhytoResp=0; saveDOfluxPhytoPhoto=0; saveDOfluxBotalgResp=0; saveDOfluxBotalgPhoto=0
		saveDOfluxSOD=0; saveDOfluxCOD=0; saveDOfluxHyporheic=0; saveDOfluxTribs=0
		saveDOfluxAdvecDisp=0; saveDOfluxHeadwater=0
		saveCO2fluxReaer=0; saveCO2fluxCBODfast=0; saveCO2fluxCBODslow=0
		saveCO2fluxPhytoResp=0; saveCO2fluxPhytoPhoto=0; saveCO2fluxBotalgResp=0; saveCO2fluxBotalgPhoto=0
		saveCO2fluxSOD=0; saveCO2fluxHyporheic=0; saveCO2fluxTribs=0
		saveCO2fluxAdvecDisp=0; saveCO2fluxHeadwater=0

		!gp 20-Oct-04 DO i=2, 10	
		!gp 05-Jul-05 DO i=2, 29
		DO i=2, 64 
			IF (status(i)==1) THEN 
					WRITE(8,*) '** Class_Integration:Integration_ failed. Insufficient Memory **'
					STOP !Class_Integration:Integration_ failed. Insufficient Memory!'
			END IF
		END DO

	END FUNCTION Integration_

END MODULE Class_IntegrationData