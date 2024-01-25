
!classoutput.f90

MODULE Class_Output
USE nrtype
USE Class_SourceIn
USE Class_IntegrationData

IMPLICIT NONE

TYPE Outdata_type
	!gp 24-Oct-04 INTEGER(I4B) nj
	INTEGER(I4B) nj		!gp need long integer for nj
	REAL(DP), DIMENSION(:), POINTER :: tdy
	!gp 27-Oct-04 add dimension for nl
	!gp REAL(DP), DIMENSION(:,:), POINTER :: cmn, cmx, cav		!constituents concentrations
	!gp REAL(DP), DIMENSION(:), POINTER :: Temn, Temx, Teav, osav	!temperature, saturated DO
	!gp REAL(DP), DIMENSION(:), POINTER :: pHmn, pHmx, pHav, pHsav	!pH
	!gp REAL(DP), DIMENSION(:), POINTER :: TNmn, TNmx, TNav		!Total nitrogen
	!gp REAL(DP), DIMENSION(:), POINTER :: TPmn, TPmx, TPav		!Total Phosphorus
	!gp REAL(DP), DIMENSION(:), POINTER :: NH3mn, NH3mx, NH3av		!NH3
	!gp REAL(DP), DIMENSION(:,:), POINTER :: pHpr, NINbpr, NIPbpr 	!bottom algae luxury update 
	!gp REAL(DP), DIMENSION(:,:), POINTER :: phitotalSavepr, phitSavepr, philSavepr, &
	!gp 						phinSavepr, phipSavepr, phicSavepr 	!gp 20-Oct-04 bottom algae growth limitation factors 
	!gp REAL(DP), DIMENSION(:,:,:), POINTER :: cpr, Tepr		!Print out
	REAL(DP), DIMENSION(:,:,:), POINTER :: cmn, cmx, cav		!constituents concentrations	(nr, nv, nl)
	REAL(DP), DIMENSION(:,:), POINTER :: Temn, Temx, Teav, osav	!temperature, saturated DO		(nr, nl)
	REAL(DP), DIMENSION(:,:), POINTER :: pHmn, pHmx, pHav, pHsav	!pH							(nr, nl)
	REAL(DP), DIMENSION(:,:), POINTER :: TNmn, TNmx, TNav		!Total nitrogen					(nr, nl)
	REAL(DP), DIMENSION(:,:), POINTER :: TPmn, TPmx, TPav		!Total Phosphorus				(nr, nl)
	REAL(DP), DIMENSION(:,:), POINTER :: NH3mn, NH3mx, NH3av		!NH3						(nr, nl)
	REAL(DP), DIMENSION(:,:,:), POINTER :: pHpr 	!bottom algae luxury update					(nr, nj, nl) 
	REAL(DP), DIMENSION(:,:), POINTER :: NINbpr, NIPbpr 	!bottom algae luxury update			(nr, nj) 

	!gp 20-Oct-04 growth limitation factors for bottom algae (nr, nj)
	REAL(DP), DIMENSION(:,:), POINTER :: phitotalSavepr, phitSavepr, philSavepr, phinSavepr, phipSavepr, phicSavepr  

	!gp 28-Oct-04 diagenesis flux between sediment/water (nr, nj)
	REAL(DP), DIMENSION(:,:), POINTER :: DiagFluxDOpr, DiagFluxCBODpr, DiagFluxNH4pr, DiagFluxNO3pr, DiagFluxSRPpr, DiagFluxICpr  

	!gp 28-Oct-04 hyporheic exchange flux between sediment/water (nr, nj)
	REAL(DP), DIMENSION(:,:), POINTER :: HypoFluxDOpr, HypoFluxCBODpr, HypoFluxNH4pr, HypoFluxNO3pr, HypoFluxSRPpr, HypoFluxICpr  

	!gp 15-Nov-04 reach-average daily-average flux between sediment/water (nr)
	REAL(DP), DIMENSION(:), POINTER :: DiagFluxDOav, DiagFluxCBODav, DiagFluxNH4av, DiagFluxNO3av, DiagFluxSRPav  
	REAL(DP), DIMENSION(:), POINTER :: HypoFluxDOav, HypoFluxCBODav, HypoFluxNH4av, HypoFluxNO3av, HypoFluxSRPav  

	!gp 11-Jan-05 reach-min/max/mean cell quota mgN/gD and mgP/gD (nr)
	REAL(DP), DIMENSION(:), POINTER :: NINbmn, NINbmx, NINbav, NIPbmn, NIPbmx, NIPbav  

	REAL(DP), DIMENSION(:,:,:), POINTER :: Tepr		!Print out									(nr, nj, nl)
	REAL(DP), DIMENSION(:,:,:,:), POINTER :: cpr		!Print out								(nr, nv, nj, nl)

	!gp 05-Jul-05 heat/DO/CO2 fluxes
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveHeatFluxJsnt, pr_saveHeatFluxLongat, pr_saveHeatFluxBack, pr_saveHeatFluxConv
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveHeatFluxEvap, pr_saveHeatFluxJsed, pr_saveHeatFluxJhyporheic, pr_saveHeatFluxTribs
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveHeatFluxAdvecDisp
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveDOfluxReaer, pr_saveDOfluxCBODfast, pr_saveDOfluxCBODslow, pr_saveDOfluxNitrif
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveDOfluxPhytoResp, pr_saveDOfluxPhytoPhoto, pr_saveDOfluxBotalgResp, pr_saveDOfluxBotalgPhoto
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveDOfluxSOD, pr_saveDOfluxCOD, pr_saveDOfluxHyporheic, pr_saveDOfluxTribs
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveDOfluxAdvecDisp
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveCO2fluxReaer, pr_saveCO2fluxCBODfast, pr_saveCO2fluxCBODslow
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveCO2fluxPhytoResp, pr_saveCO2fluxPhytoPhoto, pr_saveCO2fluxBotalgResp, pr_saveCO2fluxBotalgPhoto
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveCO2fluxSOD, pr_saveCO2fluxHyporheic, pr_saveCO2fluxTribs
	REAL(DP), DIMENSION(:,:), POINTER :: pr_saveCO2fluxAdvecDisp

	!gp 25-Jun-09
	REAL(DP), DIMENSION(:), POINTER :: av_BotAlgPhoto, av_BotAlgResp, av_BotAlgDeath, av_BotAlgNetGrowth  

END TYPE

CONTAINS

!gp 17-Nov-04 FUNCTION	outData_(nr) RESULT(pr) 
FUNCTION	outData_(nr, sys) RESULT(pr)		!gp 17-Nov-04 pass sys to minimize size of dynamic diel arrays

	!gp 17-Nov-04
	USE Class_SystemParams

	INTEGER(I4B), INTENT(IN) :: nr
	TYPE(OutData_type) pr

	!gp 05-Jul-05 INTEGER(I4B) i,status(0:60)		!gp 11-Jan-05
	!gp 25-Jun-09 INTEGER(I4B) i,status(0:93)		!gp 05-Jul-05
	INTEGER(I4B) i,status(0:97)		!gp 25-Jun-09

	!gp 17-Nov-04
	TYPE(SystemParams) sys
	INTEGER(I4B) nsteps
	IF (sys%Imeth == "Adaptive step") THEN
		nsteps = 2400
	ELSE
		nsteps = sys%nc		!minimizes array sizes for Euler and RK4 integration 
	END IF

	status=0

	!gp 17-Nov-04 ALLOCATE(pr%tdy(0:2400), STAT=status(0))
	ALLOCATE(pr%tdy(0:nsteps), STAT=status(0))		!gp 17-Nov-04 replace 2400 with nsteps for all dynamic diel output arrays below

	!gp ALLOCATE(pr%cpr(0:nr, nv, 0:2400), STAT=status(1)) 
	!gp ALLOCATE(pr%Tepr(0:nr, 0:2400, 2), STAT=status(2))
	!gp ALLOCATE(pr%pHpr(0:nr, 0:2400), STAT=status(3))
	ALLOCATE(pr%cpr(0:nr, nv, 0:nsteps, nl), STAT=status(1)) 
	ALLOCATE(pr%Tepr(0:nr, 0:nsteps, nl), STAT=status(2))
	ALLOCATE(pr%pHpr(0:nr, 0:nsteps, nl), STAT=status(3))		!gp end new block

	ALLOCATE(pr%NINbpr(0:nr, 0:nsteps),  STAT=status(4))
	ALLOCATE(pr%NIPbpr(0:nr, 0:nsteps),  STAT=status(5))

	!gp ALLOCATE(pr%Temn(0:nr), STAT=status(7))
	!gp ALLOCATE(pr%Temx(0:nr), STAT=status(8))
	!gp ALLOCATE(pr%Teav(0:nr), STAT=status(9))
	!gp ALLOCATE(pr%osav(0:nr), STAT=status(10))
	!gp ALLOCATE(pr%pHsav(0:nr), STAT=status(11))
	!gp ALLOCATE(pr%cmn(0:nr,nv), STAT=status(12))
	!gp ALLOCATE(pr%cmx(0:nr,nv), STAT=status(13))
	!gp ALLOCATE(pr%cav(0:nr,nv), STAT=status(14))
	!gp ALLOCATE(pr%pHmn(0:nr), STAT=status(15))
	!gp ALLOCATE(pr%pHmx(0:nr), STAT=status(16))
	!gp ALLOCATE(pr%pHav(0:nr), STAT=status(17))
	!gp ALLOCATE(pr%TNmn(0:nr), STAT=status(18))
	!gp ALLOCATE(pr%TNmx(0:nr), STAT=status(19))
	!gp ALLOCATE(pr%TNav(0:nr), STAT=status(20))
	!gp ALLOCATE(pr%TPmn(0:nr), STAT=status(21))
	!gp ALLOCATE(pr%TPmx(0:nr), STAT=status(22))
	!gp ALLOCATE(pr%TPav(0:nr), STAT=status(23))
	!gp ALLOCATE(pr%NH3mn(0:nr), STAT=status(24))
	!gp ALLOCATE(pr%NH3mx(0:nr), STAT=status(25))
	!gp ALLOCATE(pr%NH3av(0:nr), STAT=status(26))
	ALLOCATE(pr%Temn(0:nr, nl), STAT=status(7))
	ALLOCATE(pr%Temx(0:nr, nl), STAT=status(8))
	ALLOCATE(pr%Teav(0:nr, nl), STAT=status(9))
	ALLOCATE(pr%osav(0:nr, nl), STAT=status(10))
	ALLOCATE(pr%pHsav(0:nr, nl), STAT=status(11))
	ALLOCATE(pr%cmn(0:nr,nv, nl), STAT=status(12))
	ALLOCATE(pr%cmx(0:nr,nv, nl), STAT=status(13))
	ALLOCATE(pr%cav(0:nr,nv, nl), STAT=status(14))
	ALLOCATE(pr%pHmn(0:nr, nl), STAT=status(15))
	ALLOCATE(pr%pHmx(0:nr, nl), STAT=status(16))
	ALLOCATE(pr%pHav(0:nr, nl), STAT=status(17))
	ALLOCATE(pr%TNmn(0:nr, nl), STAT=status(18))
	ALLOCATE(pr%TNmx(0:nr, nl), STAT=status(19))
	ALLOCATE(pr%TNav(0:nr, nl), STAT=status(20))
	ALLOCATE(pr%TPmn(0:nr, nl), STAT=status(21))
	ALLOCATE(pr%TPmx(0:nr, nl), STAT=status(22))
	ALLOCATE(pr%TPav(0:nr, nl), STAT=status(23))
	ALLOCATE(pr%NH3mn(0:nr, nl), STAT=status(24))
	ALLOCATE(pr%NH3mx(0:nr, nl), STAT=status(25))
	ALLOCATE(pr%NH3av(0:nr, nl), STAT=status(26))		!gp end new block

	!gp 20-Oct-04 growth limitation factors for bottom algae
	ALLOCATE(pr%phitotalSavepr(0:nr, 0:nsteps), STAT=status(27))	
	ALLOCATE(pr%phitSavepr(0:nr, 0:nsteps), STAT=status(28))		
	ALLOCATE(pr%philSavepr(0:nr, 0:nsteps), STAT=status(29))	
	ALLOCATE(pr%phinSavepr(0:nr, 0:nsteps), STAT=status(30))		
	ALLOCATE(pr%phipSavepr(0:nr, 0:nsteps), STAT=status(31))		
	ALLOCATE(pr%phicSavepr(0:nr, 0:nsteps), STAT=status(32))		!gp 20-Oct-04 end new block

	!gp 28-Oct-04 diagenesis flux between sediment/water
	ALLOCATE(pr%DiagFluxDOpr(0:nr, 0:nsteps), STAT=status(33))	
	ALLOCATE(pr%DiagFluxCBODpr(0:nr, 0:nsteps), STAT=status(34))	
	ALLOCATE(pr%DiagFluxNH4pr(0:nr, 0:nsteps), STAT=status(35))	
	ALLOCATE(pr%DiagFluxNO3pr(0:nr, 0:nsteps), STAT=status(36))	
	ALLOCATE(pr%DiagFluxSRPpr(0:nr, 0:nsteps), STAT=status(37))	
	ALLOCATE(pr%DiagFluxICpr(0:nr, 0:nsteps), STAT=status(38))	!gp end new block	

	!gp 28-Oct-04 diagenesis flux between sediment/water
	ALLOCATE(pr%HypoFluxDOpr(0:nr, 0:nsteps), STAT=status(39))	
	ALLOCATE(pr%HypoFluxCBODpr(0:nr, 0:nsteps), STAT=status(40))	
	ALLOCATE(pr%HypoFluxNH4pr(0:nr, 0:nsteps), STAT=status(41))	
	ALLOCATE(pr%HypoFluxNO3pr(0:nr, 0:nsteps), STAT=status(42))	
	ALLOCATE(pr%HypoFluxSRPpr(0:nr, 0:nsteps), STAT=status(43))	
	ALLOCATE(pr%HypoFluxICpr(0:nr, 0:nsteps), STAT=status(44))	!gp end new block	

	!gp 15-Nov-04 reach-average daily-average flux between sediment/water
	ALLOCATE(pr%DiagFluxDOav(0:nr), STAT=status(45))	
	ALLOCATE(pr%DiagFluxCBODav(0:nr), STAT=status(46))	
	ALLOCATE(pr%DiagFluxNH4av(0:nr), STAT=status(47))	
	ALLOCATE(pr%DiagFluxNO3av(0:nr), STAT=status(48))	
	ALLOCATE(pr%DiagFluxSRPav(0:nr), STAT=status(49))	
	ALLOCATE(pr%HypoFluxDOav(0:nr), STAT=status(50))	
	ALLOCATE(pr%HypoFluxCBODav(0:nr), STAT=status(51))	
	ALLOCATE(pr%HypoFluxNH4av(0:nr), STAT=status(52))	
	ALLOCATE(pr%HypoFluxNO3av(0:nr), STAT=status(53))	
	ALLOCATE(pr%HypoFluxSRPav(0:nr), STAT=status(54))	

	!gp 11-Jan-05 cell quota mgN/gD and mgP/gD
	ALLOCATE(pr%NINbmn(0:nr), STAT=status(55))	
	ALLOCATE(pr%NINbmx(0:nr), STAT=status(56))	
	ALLOCATE(pr%NINbav(0:nr), STAT=status(57))	
	ALLOCATE(pr%NIPbmn(0:nr), STAT=status(58))	
	ALLOCATE(pr%NIPbmx(0:nr), STAT=status(59))	
	ALLOCATE(pr%NIPbav(0:nr), STAT=status(60))	

	!gp 05-Jul-05 heat/DO/CO2 fluxes
	ALLOCATE(pr%pr_saveHeatFluxJsnt(0:nr, 0:nsteps), STAT=status(61))	
	ALLOCATE(pr%pr_saveHeatFluxLongat(0:nr, 0:nsteps), STAT=status(62))	
	ALLOCATE(pr%pr_saveHeatFluxBack(0:nr, 0:nsteps), STAT=status(63))	
	ALLOCATE(pr%pr_saveHeatFluxConv(0:nr, 0:nsteps), STAT=status(64))	
	ALLOCATE(pr%pr_saveHeatFluxEvap(0:nr, 0:nsteps), STAT=status(65))	
	ALLOCATE(pr%pr_saveHeatFluxJsed(0:nr, 0:nsteps), STAT=status(66))	
	ALLOCATE(pr%pr_saveHeatFluxJhyporheic(0:nr, 0:nsteps), STAT=status(67))	
	ALLOCATE(pr%pr_saveHeatFluxTribs(0:nr, 0:nsteps), STAT=status(68))	
	ALLOCATE(pr%pr_saveHeatFluxAdvecDisp(0:nr, 0:nsteps), STAT=status(69))	
	ALLOCATE(pr%pr_saveDOfluxReaer(0:nr, 0:nsteps), STAT=status(70))	
	ALLOCATE(pr%pr_saveDOfluxCBODfast(0:nr, 0:nsteps), STAT=status(71))	
	ALLOCATE(pr%pr_saveDOfluxCBODslow(0:nr, 0:nsteps), STAT=status(72))	
	ALLOCATE(pr%pr_saveDOfluxCOD(0:nr, 0:nsteps), STAT=status(73))	
	ALLOCATE(pr%pr_saveDOfluxNitrif(0:nr, 0:nsteps), STAT=status(74))	
	ALLOCATE(pr%pr_saveDOfluxPhytoResp(0:nr, 0:nsteps), STAT=status(75))	
	ALLOCATE(pr%pr_saveDOfluxPhytoPhoto(0:nr, 0:nsteps), STAT=status(76))	
	ALLOCATE(pr%pr_saveDOfluxBotalgResp(0:nr, 0:nsteps), STAT=status(77))	
	ALLOCATE(pr%pr_saveDOfluxBotalgPhoto(0:nr, 0:nsteps), STAT=status(78))	
	ALLOCATE(pr%pr_saveDOfluxSOD(0:nr, 0:nsteps), STAT=status(79))	
	ALLOCATE(pr%pr_saveDOfluxHyporheic(0:nr, 0:nsteps), STAT=status(80))	
	ALLOCATE(pr%pr_saveDOfluxTribs(0:nr, 0:nsteps), STAT=status(81))	
	ALLOCATE(pr%pr_saveDOfluxAdvecDisp(0:nr, 0:nsteps), STAT=status(82))	
	ALLOCATE(pr%pr_saveCO2fluxReaer(0:nr, 0:nsteps), STAT=status(83))	
	ALLOCATE(pr%pr_saveCO2fluxCBODfast(0:nr, 0:nsteps), STAT=status(84))	
	ALLOCATE(pr%pr_saveCO2fluxCBODslow(0:nr, 0:nsteps), STAT=status(85))	
	ALLOCATE(pr%pr_saveCO2fluxPhytoResp(0:nr, 0:nsteps), STAT=status(86))	
	ALLOCATE(pr%pr_saveCO2fluxPhytoPhoto(0:nr, 0:nsteps), STAT=status(87))	
	ALLOCATE(pr%pr_saveCO2fluxBotalgResp(0:nr, 0:nsteps), STAT=status(88))	
	ALLOCATE(pr%pr_saveCO2fluxBotalgPhoto(0:nr, 0:nsteps), STAT=status(89))	
	ALLOCATE(pr%pr_saveCO2fluxSOD(0:nr, 0:nsteps), STAT=status(90))	
	ALLOCATE(pr%pr_saveCO2fluxHyporheic(0:nr, 0:nsteps), STAT=status(91))	
	ALLOCATE(pr%pr_saveCO2fluxTribs(0:nr, 0:nsteps), STAT=status(92))	
	ALLOCATE(pr%pr_saveCO2fluxAdvecDisp(0:nr, 0:nsteps), STAT=status(93))	

	!gp 25-Jun-09
	ALLOCATE(pr%av_BotAlgPhoto(0:nr), STAT=status(94))	
	ALLOCATE(pr%av_BotAlgResp(0:nr), STAT=status(95))	
	ALLOCATE(pr%av_BotAlgDeath(0:nr), STAT=status(96))	
	ALLOCATE(pr%av_BotAlgNetGrowth(0:nr), STAT=status(97))	

	pr%Temn=0;		pr%Temx=0;		pr%Teav=0
	pr%osav=0;		pr%pHsav=0;	pr%cmn=0
	pr%cmx=0;		pr%cav=0;		pr%pHmn=0
	pr%pHmx=0;		pr%pHav=0;		pr%TNmn=0
	pr%TNmx=0;		pr%TNav=0;   pr%TPmn=0
	pr%TPmx=0;		pr%TPav=0;		pr%NH3mn=0
	pr%NH3mx=0;		pr%NH3av=0;
	pr%cpr=0;		pr%Tepr=0;		pr%pHpr=0
	pr%NINbpr=0;	pr%NIPbpr=0

	!gp 20-Oct-04 
	pr%phitotalSavepr=0; pr%phitSavepr=0; pr%philSavepr=0; pr%phinSavepr=0; pr%phipSavepr=0; pr%phicSavepr=0

	!gp 28-Oct-04 sed fluxes at each calc step
	pr%DiagFluxDOpr=0; pr%DiagFluxCBODpr=0; pr%DiagFluxNH4pr=0; pr%DiagFluxNO3pr=0; pr%DiagFluxSRPpr=0; pr%DiagFluxICpr=0 
	pr%HypoFluxDOpr=0; pr%HypoFluxCBODpr=0; pr%HypoFluxNH4pr=0; pr%HypoFluxNO3pr=0; pr%HypoFluxSRPpr=0; pr%HypoFluxICpr=0 

	!gp 15-Nov-04 reach-average daily average sed fluxes
	pr%DiagFluxDOav=0; pr%DiagFluxCBODav=0; pr%DiagFluxNH4av=0; pr%DiagFluxNO3av=0; pr%DiagFluxSRPav=0 
	pr%HypoFluxDOav=0; pr%HypoFluxCBODav=0; pr%HypoFluxNH4av=0; pr%HypoFluxNO3av=0; pr%HypoFluxSRPav=0 

	!gp 11-Jan-05 cell quota mgN/gD and mgP/gD
	pr%NINbmn=0; pr%NINbmx=0; pr%NINbav=0; pr%NIPbmn=0; pr%NIPbmx=0; pr%NIPbav=0 


	!gp 05-Jul-05
	pr%pr_saveHeatFluxJsnt=0; pr%pr_saveHeatFluxLongat=0; pr%pr_saveHeatFluxBack=0; pr%pr_saveHeatFluxConv=0
	pr%pr_saveHeatFluxEvap=0; pr%pr_saveHeatFluxJsed=0; pr%pr_saveHeatFluxJhyporheic=0; pr%pr_saveHeatFluxTribs=0
	pr%pr_saveHeatFluxAdvecDisp=0
	pr%pr_saveDOfluxReaer=0; pr%pr_saveDOfluxCBODfast=0; pr%pr_saveDOfluxCBODslow=0; pr%pr_saveDOfluxNitrif=0
	pr%pr_saveDOfluxPhytoResp=0; pr%pr_saveDOfluxPhytoPhoto=0; pr%pr_saveDOfluxBotalgResp=0; pr%pr_saveDOfluxBotalgPhoto=0
	pr%pr_saveDOfluxSOD=0; pr%pr_saveDOfluxCOD=0; pr%pr_saveDOfluxHyporheic=0; pr%pr_saveDOfluxTribs=0
	pr%pr_saveDOfluxAdvecDisp=0
	pr%pr_saveCO2fluxReaer=0; pr%pr_saveCO2fluxCBODfast=0; pr%pr_saveCO2fluxCBODslow=0
	pr%pr_saveCO2fluxPhytoResp=0; pr%pr_saveCO2fluxPhytoPhoto=0; pr%pr_saveCO2fluxBotalgResp=0; pr%pr_saveCO2fluxBotalgPhoto=0
	pr%pr_saveCO2fluxSOD=0; pr%pr_saveCO2fluxHyporheic=0; pr%pr_saveCO2fluxTribs=0
	pr%pr_saveCO2fluxAdvecDisp=0

	!gp 25-Jun-09
	pr%av_BotAlgPhoto=0; pr%av_BotAlgResp=0; pr%av_BotAlgDeath=0; pr%av_BotAlgNetGrowth=0

	!gp 15-Nov-04 DO i=0, 26
	!gp 05-Jul-05 DO i=0, 60
	DO i=0, 93
		IF (status(i)==1) THEN 
				WRITE(8,*) '** Class_Output:outData_ failed. Insufficient memory for dynamic diel output arrays. **'
				CLOSE (8)		!gp 17-Nov-04
				STOP !Class_Integration:Integration_ failed. Insufficient Memory!'
		END IF
			
		!gp debug
		!OPEN (unit=9, FILE='debug2.out', status='REPLACE', ACTION='WRITE')
		!WRITE(9,*) 'status(', i, ') =', status(i)
		!CLOSE (9)

	END DO

END FUNCTION outData_

SUBROUTINE Output(pr, nr, topo, hydrau, Rates, system)
	USE Class_Hydraulics
	USE Class_Rates
	USE Class_RiverTopo
	USE Class_SystemParams
	USE Class_Phsolve
!	USE Class_Integration

	TYPE(Outdata_type), INTENT(IN) :: pr
	TYPE(RiverTopo_type) Topo
	TYPE(RiverHydraulics_type), INTENT(IN) :: hydrau
	TYPE(Rates_type) Rates
	TYPE(SystemParams) system

	INTEGER(I4B), INTENT(IN) :: nr	
	!gp INTEGER(I4B) i, j, nrp
	INTEGER(I4B) i, j, nrp, k	!gp
	REAL(DP) TOC, TKN, TSS, TP, TN, BottomAlgae, DOSat, NH3
	REAL(DP) kawind, ka
	CHARACTER(LEN=30) reaFormular
	INTEGER(I4B) ihour
	REAL(DP) t, pH, CBODu, tdy
	!gp output of hydraulics
	!Sheets("Hydraulics Summary")
	
	!calculate ka(0) for output, not used in calculations
	!kawind = 0.728_dp * Uw(0) ** 0.5_dp - 0.317_dp * Uw(0) + 0.0372_dp * Uw(0) ** 2
	!ka(0) = kau(0) + kawind / depth(0)
	WRITE(8,*) '** Hydraulics Summary **'
	WRITE (8,'(A9, 9A12, A25)')'Downstream', 'Hydraulics', "E'", 'H', 'B', 'Ac', 'U', 'trav time', &
					'slope', 'Reaeration', 'Reaeration formulas'
	WRITE (8,'(A9, 9A12, A25)') 'distance', 'Q,m3/s', 'm3/s', 'm', 'm', 'm2', 'mps', 'd', '', &
					'ka,20,/d', 'water/wind'
	DO i = 0, nr
		IF (i==0) THEN
			reaFormular =''
		ELSE
			IF (Rates%kawindmethod=='None') THEN
				reaFormular = TRIM(hydrau%reach(i)%kaf) // '/No wind'
			ELSE
				reaFormular = TRIM(hydrau%reach(i)%kaf) // Rates%kawindmethod
			END IF
		END IF
		WRITE(8,'(F9.4, 9F12.5, 1A35)') topo%reach(i)%xrdn, hydrau%reach(i)%Q, hydrau%reach(i)%Epout, &
					hydrau%reach(i)%depth, hydrau%reach(i)%b, hydrau%reach(i)%Ac, &
					hydrau%reach(i)%U, hydrau%reach(i)%trav, hydrau%reach(i)%s, &
					hydrau%reach(i)%ka, reaFormular
	END DO

!output of Hourly summary of loads
!sheets("Source Summary")
	WRITE(8,*)
	WRITE(8,*) '** Source summary **'
	!gp 30-Nov-04
	!gp WRITE(8,'(A6, A2, 2A15, 20A12)') 'Time', '', 'Reach', 'Downstream', 'UpDist', 'Down Dist', 'Abstraction', 'Inflow', &
	!gp 				'Temp', 'Cond', 'Iss', 'Oxygen', 'CBODs', 'CBODf', 'No', 'NH4', &
	!gp 													 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogens', 'Alk', 'pH'
	!gp WRITE(8,'(A8, 2A15, 20A12)') '', 'Label', 'Label', 'x(km)', 'x(km)', 'cms', 'cms', 'C', 'umhos', 'mgD/L', &
	!gp 				'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', 'ugN/L', 'ugN/L', &
	!gp 				'ugN/L', 'ugN/L', 'ugN/L', 'mgD/L', 'cfu/100mL', 'mgCaCO3/L'
	WRITE(8,'(A6, A2, 2A15, 21A12)') 'Time', '', 'Reach', 'Downstream', 'UpDist', 'Down Dist', 'Abstraction', 'Inflow', &
					'Temp', 'Cond', 'Iss', 'Oxygen', 'CBODs', 'CBODf', 'No', 'NH4', &
														 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogens', 'Generic', 'Alk', 'pH'
	WRITE(8,'(A8, 2A15, 21A12)') '', 'Label', 'Label', 'x(km)', 'x(km)', 'cms', 'cms', 'C', 'umhos', 'mgD/L', &
					'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', 'ugN/L', 'ugN/L', &
					'ugN/L', 'ugN/L', 'ugN/L', 'mgD/L', 'cfu/100mL', 'user define', 'mgCaCO3/L'	!gp 30-Nov-04 end new block
	IF (npt /= 0 .OR. ndiff /= 0) THEN
	!	rlab2(0) = rlab1(1)
		DO ihour = 0, 23						!loop through output of hourly sources
			t = ihour / 24.0_dp       !current output time in days
			!evaluate sine functions and distribute loads to reaches at time t
			CALL SourcesCalc(t, nr, Hydrau%flag)
			
			!output of simulation date+time
			!ActiveCell.Value = DateSerial(xyear, xmon, xday) + t
			
			DO i = 1, nr
				
				IF (load(i)%c(nv - 2) > 0 .AND. load(i)%c(nv - 1) > 0) THEN											!gp 03-Dec-04

						!gp 23-Nov-09
						!IF (system%IMethpH == "Newton-Raphson") THEN
						!	CALL phsolNewton(pH, load(i)%c(nv - 1), load(i)%Te, load(i)%c(nv - 2), load(i)%c(1))	!gp 03-Dec-04
						!ELSE
						!	CALL phsolBisect(pH, load(i)%c(nv - 1), load(i)%Te, load(i)%c(nv - 2), load(i)%c(1))	!gp 03-Dec-04
						!END IF
						IF (system%IMethpH == "Newton-Raphson") THEN
							CALL phsolNewton(pH, load(i)%c(nv - 1), load(i)%Te, load(i)%c(nv - 2), load(i)%c(1))
						ELSEIF (system%IMethpH == "Bisection") THEN
							CALL phsolBisect(pH, load(i)%c(nv - 1), load(i)%Te, load(i)%c(nv - 2), load(i)%c(1))
						ELSE
							CALL phsolBrent(pH, load(i)%c(nv - 1), load(i)%Te, load(i)%c(nv - 2), load(i)%c(1))	
						END IF

				ELSE
						pH=0.0
				END IF
				!gp 30-Nov-04
				!gp WRITE(8,'(F6.5, A2, 2A15, 20F12.4)') t,'', Topo%reach(i)%rname, Topo%reach(i)%rlab, Topo%reach(i-1)%xrdn, &
				!gp 		Topo%reach(i)%xrdn, hydrau%reach(i)%Qpta,  hydrau%reach(i)%Qpt, &
				!gp 		load(i)%Te, (load(i)%c(j), j=1, nv-2), pH
				WRITE(8,'(F6.5, A2, 2A15, 21F12.4)') t,'', Topo%reach(i)%rname, Topo%reach(i)%rlab, Topo%reach(i-1)%xrdn, &
						Topo%reach(i)%xrdn, hydrau%reach(i)%Qpta,  hydrau%reach(i)%Qpt, &
						load(i)%Te, (load(i)%c(j), j=1, nv-2), pH		!gp 30-Nov-04 add generic constituent
			END DO
		END DO
	END IF

! Output temperature for the water column
! Sheets("Temperature Output")
	WRITE(8,*)
	WRITE(8,*) '** Temperature summary (water column temperature) **'
	WRITE(8,'(A5, A10, 4A10)') 'Reach', '', 'Distance','Temp(C)', 'Temp(C)', 'Temp(C)'
	WRITE(8,'(A5, A10, 4A10)') 'Label', '', 'x(km)','Average', 'Minimum', 'Maximum'
	j = 1	!gp 27-Oct-04 water column is layer 1
	DO i=0, nr
		!gp 27-Oct-04 WRITE (8,'(A15, 4F10.4)') Topo%reach(i)%rname, Topo%reach(i)%xpm,  pr%Teav(i), &
		!gp				pr%Temn(i), pr%Temx(i)
		WRITE (8,'(A15, 4F10.4)') Topo%reach(i)%rname, Topo%reach(i)%xpm,  pr%Teav(i, j), &
						pr%Temn(i, j), pr%Temx(i, j)	!gp add nl dimension
	END DO

! Output concentrations for the water column
	WRITE(8,*)
	WRITE(8,*) '** Daily average water quality summary (water column constituents) **'
	!gp 30-Nov-04
	!gp WRITE(8,'(A5, A10, 27A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
	!gp 					'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Alk', 'pH', &
	!gp 					'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
	!gpWRITE(8,'(A5, A10, 27A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
	!gp						'ugN/L', 'ugN/L', &
	!gp					'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'Alk', 'pH', &
	!gp					'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', &
	!gp					'', '', ''

	!gp 25-Jun-09
	!WRITE(8,'(A5, A10, 33A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
	!					'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
	!					'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
	!											'NH3', 'DO sat', 'pH sat', 'Hypo biofilm', &
	!											'NINbav', 'NIPbav', 'NINbav', 'NIPbav' 
	!WRITE(8,'(A5, A10, 33A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
	!					'ugN/L', 'ugN/L', &
	!					'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
	!					'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', &
	!					'', '', '', 'gD/m2', 'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD'		!gp 11-Jan-05 end new block
	WRITE(8,'(A5, A10, 41A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
						'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
												'NH3', 'DO sat', 'pH sat', 'Hypo biofilm', &
												'NINbav', 'NIPbav', 'NINbav', 'NIPbav', &
												'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow', &
												'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow' 
	WRITE(8,'(A5, A10, 41A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
						'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', &
						'', '', '', 'gD/m2', 'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD', &
						'gD/m2/d', 'gD/m2/d', 'gD/m2/d', 'gD/m2/d', &
						'mgA/m2/d', 'mgA/m2/d', 'mgA/m2/d', 'mgA/m2/d'	

	DO i = 0, nr
		!gp 27-Oct-04 add nl dimension
		!gpTOC = (pr%cav(i, 4) + pr%cav(i, 5)) / Rates%roc + &
		!gp			Rates%aca * pr%cav(i, 11) + Rates%aca / Rates%ada * pr%cav(i, 12)
		!gpTKN = pr%cav(i, 6) + pr%cav(i, 7) + Rates%ana * pr%cav(i, 11)
		!gpTSS = Rates%ada * pr%cav(i, 11) + pr%cav(i, 2) + pr%cav(i, 12)
		!gpCBODu = pr%cav(i, 4) + pr%cav(i, 5) + Rates%roa * pr%cav(i, 11) + &
		!gp				Rates%roc * Rates%aca / Rates%ada * pr%cav(i, 12)
		!gp!Bottom Algae as Chl a
		!gpBottomAlgae= pr%cav(i, 16) / (Rates%adc * Rates%aca)
		!gpWRITE(8,'(A15, 27F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, j), j=1, nv-2), pr%pHav(i) , &
		!gp				pr%cav(i, nv), TOC, pr%TNav(i), pr%TPav(i), &
		!gp				TKN, TSS, CBODu, BottomAlgae, pr%NH3av(i), pr%osav(i), &
		!gp				pr%pHsav(i)
		j = 1	!gp water column is layer 1
		TOC = (pr%cav(i, 4, j) + pr%cav(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cav(i, 11, j) + Rates%aca / Rates%ada * pr%cav(i, 12, j)
		TKN = pr%cav(i, 6, j) + pr%cav(i, 7, j) + Rates%ana * pr%cav(i, 11, j)
		TSS = Rates%ada * pr%cav(i, 11, j) + pr%cav(i, 2, j) + pr%cav(i, 12, j)
		CBODu = pr%cav(i, 4, j) + pr%cav(i, 5, j) + Rates%roa * pr%cav(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cav(i, 12, j)
		!Bottom Algae as Chl a
		BottomAlgae= pr%cav(i, nv, j) / (Rates%adc * Rates%aca)
		!gp 30-Nov-04
		!gp WRITE(8,'(A15, 27F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
		!gp 				pr%cav(i, nv, j), TOC, pr%TNav(i, j), pr%TPav(i, j), &
		!gp 				TKN, TSS, CBODu, BottomAlgae, pr%NH3av(i, j), pr%osav(i, j), &
		!gp 				pr%pHsav(i, j)

		!gp 25-Jun-09
		!IF (i == 0) THEN		!make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
		!	WRITE(8,'(A15, 33F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
		!				pr%cav(1, nv, j), TOC, pr%TNav(i, j), pr%TPav(i, j), &
		!				TKN, TSS, CBODu, pr%cav(1, nv, j) / (Rates%adc * Rates%aca), pr%NH3av(i, j), pr%osav(i, j), &
		!				pr%pHsav(i, j), pr%cav(1, nv, 2), &
		!				pr%NINbav(1) * Rates%mgD / Rates%mgA / 1000, pr%NIPbav(1) * Rates%mgD / Rates%mgA / 1000, &
		!				pr%NINbav(1), pr%NIPbav(1)	!gp 11-Jan-05 end new block
		!ELSE
		!	WRITE(8,'(A15, 33F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
		!				pr%cav(i, nv, j), TOC, pr%TNav(i, j), pr%TPav(i, j), &
		!				TKN, TSS, CBODu, BottomAlgae, pr%NH3av(i, j), pr%osav(i, j), &
		!				pr%pHsav(i, j), pr%cav(i, nv, 2), &
		!				pr%NINbav(i) * Rates%mgD / Rates%mgA / 1000, pr%NIPbav(i) * Rates%mgD / Rates%mgA / 1000, &
		!				pr%NINbav(i), pr%NIPbav(i)	!gp 11-Jan-05 end new block
		!END IF
		IF (i == 0) THEN		!make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
			WRITE(8,'(A15, 41F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
						pr%cav(1, nv, j), TOC, pr%TNav(i, j), pr%TPav(i, j), &
						TKN, TSS, CBODu, pr%cav(1, nv, j) / (Rates%adc * Rates%aca), pr%NH3av(i, j), pr%osav(i, j), &
						pr%pHsav(i, j), pr%cav(1, nv, 2), &
						pr%NINbav(1) * Rates%mgD / Rates%mgA / 1000, pr%NIPbav(1) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbav(1), pr%NIPbav(1), &
						pr%av_BotAlgPhoto(1), pr%av_BotAlgResp(1), pr%av_BotAlgDeath(1), pr%av_BotAlgNetGrowth(1), &
						pr%av_BotAlgPhoto(1)/Rates%ada, pr%av_BotAlgResp(1)/Rates%ada, &
						pr%av_BotAlgDeath(1)/Rates%ada, pr%av_BotAlgNetGrowth(1)/Rates%ada	
		ELSE
			WRITE(8,'(A15, 41F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
						pr%cav(i, nv, j), TOC, pr%TNav(i, j), pr%TPav(i, j), &
						TKN, TSS, CBODu, BottomAlgae, pr%NH3av(i, j), pr%osav(i, j), &
						pr%pHsav(i, j), pr%cav(i, nv, 2), &
						pr%NINbav(i) * Rates%mgD / Rates%mgA / 1000, pr%NIPbav(i) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbav(i), pr%NIPbav(i), &
						pr%av_BotAlgPhoto(i), pr%av_BotAlgResp(i), pr%av_BotAlgDeath(i), pr%av_BotAlgNetGrowth(i), &
						pr%av_BotAlgPhoto(i)/Rates%ada, pr%av_BotAlgResp(i)/Rates%ada, &
						pr%av_BotAlgDeath(i)/Rates%ada, pr%av_BotAlgNetGrowth(i)/Rates%ada	
		END IF

	END DO

! Output min concentrations
	WRITE(8,*)
	WRITE(8,*) '** Daily minimum water quality summary (water column constituents) **'
	!gp 15-Nov-04
	!gp WRITE(8,'(A5, A10, 27A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
	!gp 					'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Alk', 'pH', &
	!gp 					'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
	!gp 					'NH3'
	!gp WRITE(8,'(A5, A10, 27A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
	!gp 					'ugN/L', 'ugN/L', &
	!gp 					'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'Alk', 'pH', &
	!gp 					'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', ''
	!WRITE(8,'(A5, A10, 31A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
	WRITE(8,'(A5, A10, 31A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
						'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
						'NH3', 'Hypo biofilm', 'NINbmn', 'NIPbmn', 'NINbmn', 'NIPbmn'
	!WRITE(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
	WRITE(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'Alk', 'pH', &
						'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', '', 'gD/m2', &
						'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD'		!30-Nov-04 add hypo biofilm and generic const
	DO i = 0, nr
		!gp 27-Oct-04 add dimension for nl
		!gp TOC = (pr%cmn(i, 4) + pr%cmn(i, 5)) / Rates%roc + &
		!gp 			Rates%aca * pr%cmn(i, 11) + Rates%aca / Rates%ada * pr%cmn(i, 12)
		!gp TKN = pr%cmn(i, 6) + pr%cmn(i, 7) + Rates%ana * pr%cmn(i, 11)
		!gp TSS = Rates%ada * pr%cmn(i, 11) + pr%cmn(i, 2) + pr%cmn(i, 12)
		!gp CBODu = pr%cmn(i, 4) + pr%cmn(i, 5) + Rates%roa * pr%cmn(i, 11) + &
		!gp 				Rates%roc * Rates%aca / Rates%ada * pr%cmn(i, 12)
		!gp !Bottom Algae as Chl a
		!gp BottomAlgae= pr%cmn(i, 16) / (Rates%adc * Rates%aca)
		!gp WRITE(8,'(A15, 25F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmn(i, j), j=1, nv-2), pr%pHmn(i) , &
		!gp 				pr%cmn(i, nv), TOC, pr%TNmn(i), pr%TPmn(i), &
		!gp 				TKN, TSS, CBODu, BottomAlgae, pr%NH3mn(i)
		j = 1	!gp water column is layer 1
		TOC = (pr%cmn(i, 4, j) + pr%cmn(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cmn(i, 11, j) + Rates%aca / Rates%ada * pr%cmn(i, 12, j)
		TKN = pr%cmn(i, 6, j) + pr%cmn(i, 7, j) + Rates%ana * pr%cmn(i, 11, j)
		TSS = Rates%ada * pr%cmn(i, 11, j) + pr%cmn(i, 2, j) + pr%cmn(i, 12, j)
		CBODu = pr%cmn(i, 4, j) + pr%cmn(i, 5, j) + Rates%roa * pr%cmn(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cmn(i, 12, j)
		!Bottom Algae as Chl a
		BottomAlgae= pr%cmn(i, nv, j) / (Rates%adc * Rates%aca)
		!gp 30-Nov-04
		!gp WRITE(8,'(A15, 25F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%pHmn(i, j) , &
		!gp 				pr%cmn(i, nv, j), TOC, pr%TNmn(i, j), pr%TPmn(i, j), &
		!gp 				TKN, TSS, CBODu, BottomAlgae, pr%NH3mn(i, j)
		IF (i == 0) THEN		!make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
			WRITE(8,'(A15, 31F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%pHmn(i, j) , &
						pr%cmn(1, nv, j), TOC, pr%TNmn(i, j), pr%TPmn(i, j), &
						TKN, TSS, CBODu, pr%cmn(1, nv, j) / (Rates%adc * Rates%aca), pr%NH3mn(i, j), pr%cmn(1, nv, 2), &		!gp 15-Nov-04 end new block
						pr%NINbmn(1) * Rates%mgD / Rates%mgA / 1000, pr%NIPbmn(1) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbmn(1), pr%NIPbmn(1)								!gp 11-Jan-05 end new block
		ELSE
			WRITE(8,'(A15, 31F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%pHmn(i, j) , &
						pr%cmn(i, nv, j), TOC, pr%TNmn(i, j), pr%TPmn(i, j), &
						TKN, TSS, CBODu, BottomAlgae, pr%NH3mn(i, j), pr%cmn(i, nv, 2), &
						pr%NINbmn(i) * Rates%mgD / Rates%mgA / 1000, pr%NIPbmn(i) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbmn(i), pr%NIPbmn(i)								!gp 11-Jan-05 end new block
		END IF
	END DO

! Output max concentrations
	WRITE(8,*)
	WRITE(8,*) '** Daily maximum water quality summary (water column constituents) **'
	!gp 30-Nov-04
	!gp WRITE(8,'(A5, A10, 27A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
	!gp 					'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Alk', 'pH', &
	!gp 					'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
	!gp 					'NH3'
	!gp WRITE(8,'(A5, A10, 27A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
	!gp 					'ugN/L', 'ugN/L', &
	!gp 					'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'Alk', 'pH', &
	!gp 					'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', ''
	WRITE(8,'(A5, A10, 31A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', 'NH4', 'NO3', &
						'Porg', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'Bot Alg', 'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', 'Bot Alg', &
						'NH3', 'Hypo biofilm', 'NINbmx', 'NIPbmx', 'NINbmx', 'NIPbmx'
	WRITE(8,'(A5, A10, 31A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
						'gD/m2', '', '', '', '', 'mgD/L', '', 'mgA/m2', '', 'gD/m2', &
						'mgN/mgA', 'mgP/mgA', 'mgN/gD', 'mgP/gD'		!gp 11-Jan-05 end new block
	DO i = 0, nr
		!gp 27-Oct-04 add dimension for nl
		!gp TOC = (pr%cmx(i, 4) + pr%cmx(i, 5)) / Rates%roc + &
		!gp 			Rates%aca * pr%cmx(i, 11) + Rates%aca / Rates%ada * pr%cmx(i, 12)
		!gp TKN = pr%cmx(i, 6) + pr%cmx(i, 7) + Rates%ana * pr%cmx(i, 11)
		!gp TSS = Rates%ada * pr%cmx(i, 11) + pr%cmx(i, 2) + pr%cmx(i, 12)
		!gp CBODu = pr%cmx(i, 4) + pr%cmx(i, 5) + Rates%roa * pr%cmx(i, 11) + &
		!gp 				Rates%roc * Rates%aca / Rates%ada * pr%cmx(i, 12)
		!gp !Bottom Algae as Chl a
		!gp BottomAlgae= pr%cmx(i, 16) / (Rates%adc * Rates%aca)
		!gp WRITE(8,'(A15, 25F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmx(i, j), j=1, nv-2), pr%pHmx(i) , &
		!gp 				pr%cmx(i, nv), TOC, pr%TNmx(i), pr%TPmx(i), &
		!gp 				TKN, TSS, CBODu, BottomAlgae, pr%NH3mx(i)
		j = 1	!gp water column is layer 1
		TOC = (pr%cmx(i, 4, j) + pr%cmx(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cmx(i, 11, j) + Rates%aca / Rates%ada * pr%cmx(i, 12, j)
		TKN = pr%cmx(i, 6, j) + pr%cmx(i, 7, j) + Rates%ana * pr%cmx(i, 11, j)
		TSS = Rates%ada * pr%cmx(i, 11, j) + pr%cmx(i, 2, j) + pr%cmx(i, 12, j)
		CBODu = pr%cmx(i, 4, j) + pr%cmx(i, 5, j) + Rates%roa * pr%cmx(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cmx(i, 12, j)
		!Bottom Algae as Chl a
		BottomAlgae= pr%cmx(i, nv, j) / (Rates%adc * Rates%aca)
		!gp 30-Nov-04
		!gp WRITE(8,'(A15, 25F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%pHmx(i, j) , &
		!gp 				pr%cmx(i, nv, j), TOC, pr%TNmx(i, j), pr%TPmx(i, j), &
		!gp 				TKN, TSS, CBODu, BottomAlgae, pr%NH3mx(i, j)
		IF (i == 0) THEN		!make the reach 0 bottom algae and biofilm equal to reach 1 for output charts to look good
			WRITE(8,'(A15, 31F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%pHmx(i, j) , &
						pr%cmx(1, nv, j), TOC, pr%TNmx(i, j), pr%TPmx(i, j), &
						TKN, TSS, CBODu, pr%cmx(1, nv, j) / (Rates%adc * Rates%aca), pr%NH3mx(i, j), pr%cmx(1, nv, 2), &
						pr%NINbmx(1) * Rates%mgD / Rates%mgA / 1000, pr%NIPbmx(1) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbmx(1), pr%NIPbmx(1)		!gp 11-Jan-05 end new block
		ELSE
			WRITE(8,'(A15, 31F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%pHmx(i, j), &
						pr%cmx(i, nv, j), TOC, pr%TNmx(i, j), pr%TPmx(i, j), &
						TKN, TSS, CBODu, BottomAlgae, pr%NH3mx(i, j), pr%cmx(i, nv, 2), &
						pr%NINbmx(i) * Rates%mgD / Rates%mgA / 1000, pr%NIPbmx(i) * Rates%mgD / Rates%mgA / 1000, &
						pr%NINbmx(i), pr%NIPbmx(i)		!gp 11-Jan-05 end new block
		END IF
	END DO

! Output sediment fluxes !gp 15-Nov-04 reach-averaged and daily-averaged and include hyporheic and total flux
	WRITE(8,*)
	WRITE(8,*) '** Sediment fluxes (reach-average daily-average) **'
	!gp 15-Nov-04
	!gp WRITE(8,'(A5, A10, 6A13)') 'Reach', '', 'Distance', 'SOD', 'Flux CH4', 'Flux NH4', &
	!gp 					'Flux InorgP', 'Flux NO3'
	!gp WRITE(8,'(A5, A10, 6A13)') 'Label', '', 'x(km)', 'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', &
	!gp 					'mgP/m2/d', 'mgN/m2/d'
	WRITE(8,'(A5, A10, 16A13)') 'Reach', '', 'Distance', &
								'DiagFluxDO', 'DiagFluxCBOD', 'DiagFluxNH4', 'DiagFluxSRP', 'DiagFluxNO3', &
								'HypoFluxDO', 'HypoFluxCBOD', 'HypoFluxNH4', 'HypoFluxSRP', 'HypoFluxNO3', &
								'TotFluxDO', 'TotFluxCBOD', 'TotFluxNH4', 'TotFluxSRP', 'TotFluxNO3'
	WRITE(8,'(A5, A10, 16A13)') 'Label', '', 'x(km)', &
								'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d', &
								'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d', &
								'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d', 'mgP/m2/d', 'mgN/m2/d'		!gp 15-Nov-04 end new block
	DO i=1, nr
		!gp 15-Nov-04 reach-average daily-average sediment fluxes
		!gp WRITE(8,'(A15, 6F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, SODpr(i), &
		!gp				JCH4pr(i), JNH4pr(i), JSRPpr(i), JNO3pr(i)
		WRITE(8,'(A15, 16F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, &
					pr%DiagFluxDOav(i), pr%DiagFluxCBODav(i), pr%DiagFluxNH4av(i), pr%DiagFluxSRPav(i), pr%DiagFluxNO3av(i), &
					pr%HypoFluxDOav(i), pr%HypoFluxCBODav(i), pr%HypoFluxNH4av(i), pr%HypoFluxSRPav(i), pr%HypoFluxNO3av(i), &
					pr%DiagFluxDOav(i) + pr%HypoFluxDOav(i), &
					pr%DiagFluxCBODav(i) + pr%HypoFluxCBODav(i), &
					pr%DiagFluxNH4av(i) + pr%HypoFluxNH4av(i), &
					pr%DiagFluxSRPav(i) + pr%HypoFluxSRPav(i), &
					pr%DiagFluxNO3av(i) + pr%HypoFluxNO3av(i)		!gp 15-Nov-04 end new block
	END DO


!gp 17-Feb-05
!only output diel if showDielResults = "Yes"
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

If (system%showDielResults == "Yes") Then

!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

!diel output result for water column
	WRITE(8,*)
	WRITE(8,*) '** Diel water quality in the water column **'
	
	!gp 01-Nov-04 WRITE(8,'(32A13)') 't', 'Tempw', ' Temps', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', &
	WRITE(8,'(32A13)') 't', 'Tempw', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', &		!gp 30-Nov-04
						'NH4', 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &	!gp 30-Nov-04
									'Bot Algae', 'TSS', 'TP', 'TN', 'DOsat', 'NH3', 'IntN', 'Int P', &
									'phiTemp', 'phiLight', 'phiNitr', 'phiPhos', 'phiCarb', 'phiTotal'		!gp 20-Oct-04

	!gp 01-Nov-04 WRITE(8,'(32A13)') 'hr', 'c', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', &
	WRITE(8,'(32A13)') 'hr', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', &	!gp 32-Nov-04
						'ugN/L', 'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', '', '', '', '', 'gD/m2', 'mgD/L', &	!gp 30-Nov-04
						'ugP/L', 'ugN/L', 'mg/L', 'ugN/L', 'mgN/mgA', 'mgP/mgA', &
						'frac', 'frac', 'frac', 'frac', 'frac', 'frac'					!gp 20-Oct-04 

	!gp WRITE(8,*) pr%nj
	WRITE(8,'(I13)') pr%nj
	DO nrp=0, nr	
		!tdy=0
		DO i=0, pr%nj
			!gp 27-Oct-04 add dimension for nl
			!gp TSS = pr%cpr(nrp, 11, i) * Rates%ada + pr%cpr(nrp, 2, i) + &
			!gp 							pr%cpr(nrp, 12, i)
			!gp TP = pr%cpr(nrp, 11, i) * Rates%apa + pr%cpr(nrp, 9, i) + &
			!gp 							pr%cpr(nrp, 10, i)
			!gp TN = pr%cpr(nrp, 11, i) * Rates%ana + pr%cpr(nrp, 6, i) + &
			!gp 							pr%cpr(nrp, 7, i) + pr%cpr(nrp, 8, i)
			!gp DOSat = oxsat(pr%Tepr(nrp, i, 1), hydrau%reach(nrp)%elev)
			!gp NH3 = 1.0_dp/(1 + 10.0_dp ** (-pr%pHpr(nrp, i))/10.0_dp** -(0.09018_dp + 2729.92_dp / &
			!gp 		(pr%Tepr(nrp, i, 1) + 273.15_dp))) * pr%cpr(nrp, 7, i)
			!gp WRITE(8, '(32F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, 1),  pr%Tepr(nrp, i, 2), &
			!gp 			(pr%cpr(nrp, j, i), j=1, nv-2), pr%pHpr(nrp, i), &
			!gp 			 pr%cpr(nrp, nv, i)* Rates%mgA / Rates%mgD * 1000, & !SCC 08/09/2004
			!gp 			 TSS, TP, TN , DOSat, NH3, &
			!gp 			 pr%NINbpr(nrp, i)* Rates%mgD / Rates%mgA / 1000, & !SCC 08/09/2004
			!gp 			pr%NIPbpr(nrp, i)* Rates%mgD / Rates%mgA / 1000, &	 !SCC 08/09/2004										
			!gp 			pr%phitSavepr(nrp, i), 	pr%philSavepr(nrp, i), &	!gp 20-Oct-04										
			!gp 			pr%phinSavepr(nrp, i), 	pr%phipSavepr(nrp, i), &		!gp 20-Oct-04										
			!gp 			pr%phicSavepr(nrp, i), 	pr%phitotalSavepr(nrp, i)			!gp 20-Oct-04										
			!gp !	tdy=tdy+ system%dt * 24
			j = 1	!gp water column is layer 1
			TSS = pr%cpr(nrp, 11, i, j) * Rates%ada + pr%cpr(nrp, 2, i, j) + &
										pr%cpr(nrp, 12, i, j)
			TP = pr%cpr(nrp, 11, i, j) * Rates%apa + pr%cpr(nrp, 9, i, j) + &
										pr%cpr(nrp, 10, i, j)
			TN = pr%cpr(nrp, 11, i, j) * Rates%ana + pr%cpr(nrp, 6, i, j) + &
										pr%cpr(nrp, 7, i, j) + pr%cpr(nrp, 8, i, j)
			DOSat = oxsat(pr%Tepr(nrp, i, j), hydrau%reach(nrp)%elev)
			NH3 = 1.0_dp/(1 + 10.0_dp ** (-pr%pHpr(nrp, i, j))/10.0_dp** -(0.09018_dp + 2729.92_dp / &
					(pr%Tepr(nrp, i, j) + 273.15_dp))) * pr%cpr(nrp, 7, i, j)

			!gp 01-Nov-04 WRITE(8, '(32F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, 1),  pr%Tepr(nrp, i, 2), &
			!gp 05-Jul-05 WRITE(8, '(33F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, j),  &		
			WRITE(8, '(32F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, j),  &		!gp 05-Jul-05
						(pr%cpr(nrp, k, i, j), k=1, nv-2), pr%pHpr(nrp, i, j), &
						 pr%cpr(nrp, nv, i, j)* Rates%mgA / Rates%mgD * 1000, & !SCC 08/09/2004
						 TSS, TP, TN , DOSat, NH3, &
						 pr%NINbpr(nrp, i)* Rates%mgD / Rates%mgA / 1000, & !SCC 08/09/2004
						pr%NIPbpr(nrp, i)* Rates%mgD / Rates%mgA / 1000, &	 !SCC 08/09/2004										
						pr%phitSavepr(nrp, i), pr%philSavepr(nrp, i), &		!gp 20-Oct-04										
			 			pr%phinSavepr(nrp, i), pr%phipSavepr(nrp, i), &		!gp 20-Oct-04										
						pr%phicSavepr(nrp, i), pr%phitotalSavepr(nrp, i)	!gp 20-Oct-04										

		END DO
	END DO
 
	!gp 27-Oct-04 (all code below here is new)
	!
	! ---------------------------------------------------
	! --- output of hyporheic pore water constituents --- 
	! ---------------------------------------------------
	!

	!gp diel output result for hyporheic pore water
	WRITE(8,*)
	WRITE(8,*) '** Diel hyporheic pore water and sediment flux **'
	
	WRITE(8,'(41A13)') 't', ' Temps', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', &
						'NH4', 'NO3', 'Po', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
									'TSS', 'TP', 'TN', 'DOsat', 'NH3', &
									'DiagFluxDO', 'DiagFluxCBOD', 'DiagFluxNH4', 'DiagFluxNO3', 'DiagFluxSRP', 'DiagFluxIC', &	
									'HypoFluxDO', 'HypoFluxCBOD', 'HypoFluxNH4', 'HypoFluxNO3', 'HypoFluxSRP', 'DiagFluxIC', &	
									'TotFluxDO', 'TotFluxCBOD', 'TotFluxNH4', 'TotFluxNO3', 'TotFluxSRP', 'TotFluxIC'			
	WRITE(8,'(41A13)') 'hr', 'c', 'umhos', 'mg/L', 'mg/L', 'mgO2/L', 'ugN/L', 'ugN/L', &
						'ugN/L', 'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', '', '', '', '', 'mgD/L', &
						'ugP/L', 'ugN/L', 'mg/L', 'ugN/L', &
						'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d', &	 
						'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d', & 	 
						'gO2/m2/d', 'gO2/m2/d', 'mgN/m2/d','mgN/m2/d','mgP/m2/d','gC/m2/d'		
	WRITE(8,*) pr%nj
	DO nrp=0, nr	
		DO i=0, pr%nj
			SELECT CASE (system%simHyporheicWQ)
			CASE ('Level 1', 'Level 2')		
 				j = 2	!gp hyporheic pore water is layer 2
				TSS = pr%cpr(nrp, 11, i, j) * Rates%ada + pr%cpr(nrp, 2, i, j) + &
												pr%cpr(nrp, 12, i, j)
				TP = pr%cpr(nrp, 11, i, j) * Rates%apa + pr%cpr(nrp, 9, i, j) + &
											pr%cpr(nrp, 10, i, j)
				TN = pr%cpr(nrp, 11, i, j) * Rates%ana + pr%cpr(nrp, 6, i, j) + &
											pr%cpr(nrp, 7, i, j) + pr%cpr(nrp, 8, i, j)
				NH3 = 1.0_dp/(1 + 10.0_dp ** (-pr%pHpr(nrp, i, j))/10.0_dp** -(0.09018_dp + 2729.92_dp / &
						(pr%Tepr(nrp, i, j) + 273.15_dp))) * pr%cpr(nrp, 7, i, j)	
 				DOSat = oxsat(pr%Tepr(nrp, i, j), hydrau%reach(nrp)%elev)
				WRITE(8, '(41F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, j), &
							(pr%cpr(nrp, k, i, j), k=1, nv-2), pr%pHpr(nrp, i, j), &
							 TSS, TP, TN , DOSat, NH3, &
							pr%DiagFluxDOpr(nrp, i), pr%DiagFluxCBODpr(nrp, i), pr%DiagFluxNH4pr(nrp, i), &	
							pr%DiagFluxNO3pr(nrp, i), pr%DiagFluxSRPpr(nrp, i), pr%DiagFluxICpr(nrp, i), &	
							pr%HypoFluxDOpr(nrp, i), pr%HypoFluxCBODpr(nrp, i), pr%HypoFluxNH4pr(nrp, i), &	
							pr%HypoFluxNO3pr(nrp, i), pr%HypoFluxSRPpr(nrp, i), pr%HypoFluxICpr(nrp, i), &	
							pr%DiagFluxDOpr(nrp, i) + pr%HypoFluxDOpr(nrp, i), &		
							pr%DiagFluxCBODpr(nrp, i) + pr%HypoFluxCBODpr(nrp, i), &	
							pr%DiagFluxNH4pr(nrp, i) + pr%HypoFluxNH4pr(nrp, i), &	
							pr%DiagFluxNO3pr(nrp, i) + pr%HypoFluxNO3pr(nrp, i), &		
							pr%DiagFluxSRPpr(nrp, i) + pr%HypoFluxSRPpr(nrp, i), &		
							pr%DiagFluxICpr(nrp, i) + pr%HypoFluxICpr(nrp, i)			
			CASE DEFAULT	!gp only write sediment temperatures and diagenesis fluxes if hyporheic wq is not being simulated		
				WRITE(8, '(41F13.4)') pr%tdy(i)*24, pr%Tepr(nrp, i, 2), (0*k,k=3,23), &
							pr%DiagFluxDOpr(nrp, i), pr%DiagFluxCBODpr(nrp, i), pr%DiagFluxNH4pr(nrp, i), &	
							pr%DiagFluxNO3pr(nrp, i), pr%DiagFluxSRPpr(nrp, i), pr%DiagFluxICpr(nrp, i), &	
							0,0,0,0,0,0, &
							pr%DiagFluxDOpr(nrp, i), pr%DiagFluxCBODpr(nrp, i), pr%DiagFluxNH4pr(nrp, i), &	
							pr%DiagFluxNO3pr(nrp, i), pr%DiagFluxSRPpr(nrp, i), pr%DiagFluxICpr(nrp, i)	
			END SELECT
		END DO
	END DO


!gp 05-Jul-05 diel fluxes for heat/DO/CO2
	WRITE(8,*)
	WRITE(8,*) '** Diel fluxes of heat (W/m^2), DO (gO2/m^2/d), and CO2 (gC/m^2/d) **'
	
	WRITE(8,'(34A13)') 't', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', 'Heat', &	
							'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', 'DO', &	
							'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2', 'CO2'

	WRITE(8,'(34A13)') 'hr', 'solar', 'longat', 'back', 'air conv', 'evap', 'sed cond', 'hyporheic', 'tribs/GW', 'advec/disp', &	
						'reaeration', 'fast CBOD', 'slow CBOD', 'COD', 'nitrif', &
						'phyto resp', 'phyto photo', 'botalg resp', 'botalg photo', &
						'SOD', 'hyporheic', 'tribs/GW', 'advec/disp', &
						'reaeration', 'fast CBOD', 'slow CBOD', &
						'phyto resp', 'phyto photo', 'botalg resp', 'botalg photo', &
						'SOD', 'hyporheic', 'tribs/GW', 'advec/disp'

	WRITE(8,'(I13)') pr%nj
	DO nrp=0, nr	
		DO i=0, pr%nj
			WRITE(8, '(34F13.4)') pr%tdy(i)*24, &		
			pr%pr_saveHeatFluxJsnt(nrp, i), &				!'heat fluxes
			pr%pr_saveHeatFluxLongat(nrp, i), &
			pr%pr_saveHeatFluxBack(nrp, i), &
			pr%pr_saveHeatFluxConv(nrp, i), &
			pr%pr_saveHeatFluxEvap(nrp, i), &
			pr%pr_saveHeatFluxJsed(nrp, i), &
			pr%pr_saveHeatFluxJhyporheic(nrp, i), &
			pr%pr_saveHeatFluxTribs(nrp, i), &
			pr%pr_saveHeatFluxAdvecDisp(nrp, i), &
			pr%pr_saveDOfluxReaer(nrp, i), &				!'DO fluxes
			pr%pr_saveDOfluxCBODfast(nrp, i), &
			pr%pr_saveDOfluxCBODslow(nrp, i), &
			pr%pr_saveDOfluxCOD(nrp, i), &
			pr%pr_saveDOfluxNitrif(nrp, i), &
			pr%pr_saveDOfluxPhytoResp(nrp, i), &
			pr%pr_saveDOfluxPhytoPhoto(nrp, i), &
			pr%pr_saveDOfluxBotalgResp(nrp, i), &
			pr%pr_saveDOfluxBotalgPhoto(nrp, i), &
			pr%pr_saveDOfluxSOD(nrp, i), &
			pr%pr_saveDOfluxHyporheic(nrp, i), &
			pr%pr_saveDOfluxTribs(nrp, i), &
			pr%pr_saveDOfluxAdvecDisp(nrp, i), &
			pr%pr_saveCO2fluxReaer(nrp, i), &				!'CO2 fluxes
			pr%pr_saveCO2fluxCBODfast(nrp, i), &
			pr%pr_saveCO2fluxCBODslow(nrp, i), &
			pr%pr_saveCO2fluxPhytoResp(nrp, i), &
			pr%pr_saveCO2fluxPhytoPhoto(nrp, i), &
			pr%pr_saveCO2fluxBotalgResp(nrp, i), &
			pr%pr_saveCO2fluxBotalgPhoto(nrp, i), &
			pr%pr_saveCO2fluxSOD(nrp, i), &
			pr%pr_saveCO2fluxHyporheic(nrp, i), &
			pr%pr_saveCO2fluxTribs(nrp, i), &
			pr%pr_saveCO2fluxAdvecDisp(nrp, i)

		END DO
	END DO


!gp 17-Feb-05
!only output diel if showDielResults = "Yes"
!
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

End If

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	! Output temperature for the hyporheic pore water
	WRITE(8,*)
	WRITE(8,*) '** Temperature summary (hyporheic pore water temperature) **'
	WRITE(8,'(A5, A10, 4A10)') 'Reach', '', 'Distance','Temp(C)', 'Temp(C)', 'Temp(C)'
	WRITE(8,'(A5, A10, 4A10)') 'Label', '', 'x(km)','Average', 'Minimum', 'Maximum'
	j = 2	!gp 27-Oct-04 hypoprheic pore water is layer 2
	DO i=0, nr
		WRITE (8,'(A15, 4F10.4)') Topo%reach(i)%rname, Topo%reach(i)%xpm,  pr%Teav(i, j), &
						pr%Temn(i, j), pr%Temx(i, j)	!gp add nl dimension
	END DO

	SELECT CASE (system%simHyporheicWQ)
	CASE ('Level 1', 'Level 2')		

	! Output concentrations for the hyporheic pore water constituents
		WRITE(8,*)
		WRITE(8,*) '** Daily average water quality summary (hyporheic pore water constituents) **'
		WRITE(8,'(A5, A10, 26A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
						'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
												'NH3', 'DO sat', 'pH sat'
		WRITE(8,'(A5, A10, 26A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', '', 'user defined', 'mgCaCO3/L', 's.u.', &
						'', '', '', '', 'mgD/L', '', &
						'', '', ''
		DO i = 0, nr
			j = 2	!gp hyporheic pore water is layer 2
			TOC = (pr%cav(i, 4, j) + pr%cav(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cav(i, 11, j) + Rates%aca / Rates%ada * pr%cav(i, 12, j)
			TKN = pr%cav(i, 6, j) + pr%cav(i, 7, j) + Rates%ana * pr%cav(i, 11, j)
			TSS = Rates%ada * pr%cav(i, 11, j) + pr%cav(i, 2, j) + pr%cav(i, 12, j)
			CBODu = pr%cav(i, 4, j) + pr%cav(i, 5, j) + Rates%roa * pr%cav(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cav(i, 12, j)
			WRITE(8,'(A15, 26F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cav(i, k, j), k=1, nv-2), pr%pHav(i, j) , &
						TOC, pr%TNav(i, j), pr%TPav(i, j), &
						TKN, TSS, CBODu, pr%NH3av(i, j), pr%osav(i, j), &
						pr%pHsav(i, j)
		END DO

	! Output min concentrations
		WRITE(8,*)
		WRITE(8,*) '** Daily minimum water quality summary (hyporheic pore water constituents) **'
		WRITE(8,'(A5, A10, 24A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
						'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
						'NH3'
		WRITE(8,'(A5, A10, 24A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', 'cfu/100mL', 'user defined', 'mgCaCO3/L', 's.u.', &
						'', '', '', '', 'mgD/L', '', ''
		DO i = 0, nr
			j = 2	!gp hyporheic pore water is layer 2
			TOC = (pr%cmn(i, 4, j) + pr%cmn(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cmn(i, 11, j) + Rates%aca / Rates%ada * pr%cmn(i, 12, j)
			TKN = pr%cmn(i, 6, j) + pr%cmn(i, 7, j) + Rates%ana * pr%cmn(i, 11, j)
			TSS = Rates%ada * pr%cmn(i, 11, j) + pr%cmn(i, 2, j) + pr%cmn(i, 12, j)
			CBODu = pr%cmn(i, 4, j) + pr%cmn(i, 5, j) + Rates%roa * pr%cmn(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cmn(i, 12, j)
			WRITE(8,'(A15, 24F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmn(i, k, j), k=1, nv-2), pr%pHmn(i, j) , &
						TOC, pr%TNmn(i, j), pr%TPmn(i, j), &
						TKN, TSS, CBODu, pr%NH3mn(i, j)
		END DO

	! Output max concentrations
		WRITE(8,*)
		WRITE(8,*) '** Daily maximum water quality summary (hyporheic pore water constituents) **'
		WRITE(8,'(A5, A10, 24A13)') 'Reach', '', 'x', 'cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'No', 'NH4', 'NO3', &
						'PO', 'InorgP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'TOC', 'TN', 'TP', 'TKN', 'TSS', 'CBODu', &
						'NH3'
		WRITE(8,'(A5, A10, 24A13)') 'Label', '', 'km', 'umhos', 'mgD/L', 'mgO2/L', 'mgO2/L', 'mgO2/L', 'ugN/L', &
						'ugN/L', 'ugN/L', &
						'ugP/L', 'ugP/L', 'ugA/L', 'mgD/L', 'cfu/100mL', 'user defined', 'mgCaCO3/L', 's.u.', &
						'', '', '', '', 'mgD/L', '', ''
		DO i = 0, nr
			j = 2	!gp hyporheic pore water is layer 2
			TOC = (pr%cmx(i, 4, j) + pr%cmx(i, 5, j)) / Rates%roc + &
					Rates%aca * pr%cmx(i, 11, j) + Rates%aca / Rates%ada * pr%cmx(i, 12, j)
			TKN = pr%cmx(i, 6, j) + pr%cmx(i, 7, j) + Rates%ana * pr%cmx(i, 11, j)
			TSS = Rates%ada * pr%cmx(i, 11, j) + pr%cmx(i, 2, j) + pr%cmx(i, 12, j)
			CBODu = pr%cmx(i, 4, j) + pr%cmx(i, 5, j) + Rates%roa * pr%cmx(i, 11, j) + &
						Rates%roc * Rates%aca / Rates%ada * pr%cmx(i, 12, j)
			WRITE(8,'(A15, 24F13.6)')Topo%reach(i)%rname, Topo%reach(i)%xpm, (pr%cmx(i, k, j), k=1, nv-2), pr%pHmx(i, j) , &
						TOC, pr%TNmx(i, j), pr%TPmx(i, j), &
						TKN, TSS, CBODu, pr%NH3mx(i, j)
		END DO

	END SELECT

END SUBROUTINE
END MODULE
