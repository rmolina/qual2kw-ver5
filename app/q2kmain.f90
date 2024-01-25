!  Q2KMAIN.f90 
!
!****************************************************************************
!
!  PROGRAM: Q2KFORTRANTEST
!
!  PURPOSE:  Entry point for the console application.
!
! The following water quality state variables are simulated in QUAL2Kw ver5.x
! (i=reach, j=layer [1=water column, 2=sediment temperature/hyporheic pore water])
!
! Variable       Variable name                             Units
! ----------     -------------                             -----
! c(i, 1, j)     Conductivity                              umho/cm
! c(i, 2, j)     Inorganic Suspended Solids                mgD/L (D=dry weight)
! c(i, 3, j)     Dissolved Oxygen                          mgO/L (O=oxygen)
! c(i, 4, j)     Slow CBOD                                 mgO/L
! c(i, 5, j)     Fast CBOD                                 mgO/L
! c(i, 6, j)     Organic N (dissolved and particulate)     ugN/L
! c(i, 7, j)     Ammonia N                                 ugN/L
! c(i, 8, j)     Nitrate + Nitrite N                       ugN/L
! c(i, 9, j)     Organic P (dissolved and particulate)     ugP/L
! c(i, 10, j)    Soluble Reactive P                        ugP/L
! c(i, 11, j)    Phytoplankton                             ugA/L (A=chlorophyll a)
! c(i, 12, j)    Particulate Organic Mater (POM)           mgD/L
! c(i, 13, j)    Pathogen indicator organism               number/volume
! c(i, 14, j)    Generic constituent                       user defined
! c(i, nv-2, j)  Alkalinity                                mgCaCO3/L 
! c(i, nv-1, j)  Total Inorganic C                         moles/L  
! c(i, nv, 1)    bottom algae                              gD/m^2
! c(i, nv, 2)    hyporheic biofilm                         gD/m^2
! Te(i, j)       temperature                               deg C
! INb(i)         Internal N in bottom algae cells          mgN/m^2
! IPb(i)         Internal P in bottom algae cells          mgP/m^2
!
!****************************************************************************

	program Q2KMain

	! Variables

	USE nrtype
	USE Class_SystemParams
	USE Class_Hydraulics
	USE Class_Meteo
	USE Class_Headwater
	USE Class_Downstream
	USE Class_Rates
	USE Class_SolarCalc
!	USE Class_RiverTopo
	USE Class_ReadFile
	USE Class_Integration
!	USE Class_Output


!gp 04-Feb-05
!	USE DFPORT		!11/16/04


	IMPLICIT NONE




  ! Variables
	CHARACTER(LEN=260) :: infile, outfile	!input & output file name 
											!gp long file names are limited to 255 characters (260 for full paths)
	TYPE(RiverHydraulics_type) hydrau	!channel dimensions, hydraulics, physical characters
	TYPE(Meteo_type) Meteo			!meteology information
	TYPE(Headwater_type) HW			!headwater
	TYPE(Downstream_type) DB		!downstream boundary
	TYPE(Rates_type) Rates			!stoch, reaction, temperature and all other rate
	TYPE(solar_type) :: Solar		!solar radiation 
	TYPE(SystemParams) sys			!declare the system parameter variables 
	TYPE(RiverTopo_type) Topo		!river topology
	TYPE(Outdata_type) prOut
	INTEGER(I4B) i, j, k
	INTEGER(I4B) begintime, endtime


!gp 04-Feb-05
!	CHARACTER(LEN=260) ::msgFile, dirname	!11/16/04
!	INTEGER(I4B) istat							!11/16/04
	CHARACTER(LEN=260) ::msgFile 				!gp long file names are limited to 255 characters (260 for full paths)


	CALL SYSTEM_CLOCK(begintime)

	WRITE(*,*)
	WRITE(*,'(50A)') ' QUAL2Kw version 5.1'
	WRITE(*,'(50A)') ' Department of Ecology and Tufts University'
	WRITE(*,*)
	WRITE(*,'(50A)') ' G.J. Pelletier, S.C. Chapra, and Hua Tao'
	!WRITE(*,'(50A)') '        modified from QUAL2K ver 1.4 by Hua Tao and S.C. Chapra'
	WRITE(*,*)
	WRITE(*,'(50A)') ' Program is running, please wait...'
	WRITE(*,*)
	

	msgFile = 'message.dat'


	open (unit=8, File=msgFile, status='OLD', ACTION='READ')
	READ(8,*) infile, outfile

	CLOSE (8)


!gp debug 04-Feb-05
!open (unit=11, File='debug.out', status='replace', ACTION='READWRITE')
!write(11,*) sDir
!write(11,*) msgfile
!write(11,*) infile, outfile
!close (11)


	!infile='C:\research\qual2k\q2k\input\BC092187.q2k'
	!infile='C:\research\qual2k\input\BC092187ThackP1.q2k'
!	infile='C:\research\qual2k\input\BC092187v1_2Fortran.q2k'
!	infile='C:\research\qual2k\input\SRdummyFortran.q2k'
	OPEN (unit=8, FILE=infile, status='OLD', ACTION='READ')
	!read in file

!gp 12-Jan-06, 21-Nov-06
!	OPEN (unit=10, FILE='debug.txt', status='REPLACE', ACTION='WRITE')
!	OPEN (unit=11, FILE='debug2.txt', status='REPLACE', ACTION='WRITE')

	call ReadInputfile(sys, hydrau, Meteo, HW, DB, Rates, Topo, Solar)


	!finishing READING
	CLOSE (8)

	!do simulation	
	!gp 17-Nov-04 prOut= Outdata_(topo%nr)
	prOut= Outdata_(topo%nr, sys)		!gp 17-Nov-04 pass sys to allocate dynamic diel arrays

!	outfile= 'C:\research\qual2k\input\BC092187v1_2Fortran.out'
!	outfile= 'C:\research\qual2k\input\SRdummyFortran.out'	
	OPEN (unit=8, FILE=outfile, status='REPLACE', ACTION='WRITE')

	!gp 23-Jun-09
	IF (sys%writeDynamic == "Yes") THEN
		open (unit=12, File='dynamic.txt', status='replace', ACTION='WRITE')
		WRITE(12,'(44A13)') 'reach', 't', 'Tempw', 'Cond', 'ISS', 'DO', 'CBODs', 'CBODf', 'Norg', &
						'NH4N', 'NO23N', 'Porg', 'SRP', 'Phyto', 'Detritus', 'Pathogen', 'Generic', 'Alk', 'pH', &
						'Bot Algae', 'Bot Algae', 'TSS', 'TP', 'TN', 'DOsat', 'NH3ui', 'Int N', 'Int P', &
						'phiTemp', 'phiLight', 'phiNitr', 'phiPhos', 'phiCarb', 'phiTotal', &
						'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow', &
						'BotAlgPhoto', 'BotAlgResp', 'BotAlgDeath', 'BotAlgGrow'		
	END IF


	CALL Integration(sys, Rates, Meteo, Solar, HW, DB, &
												hydrau, prOut, topo%nr)

	!gp 23-Jun-09
	IF (sys%writeDynamic == "Yes") THEN
		close (12)
	END IF

	CALL Output(prOut, topo%nr, topo, hydrau, Rates, sys)

	CALL SYSTEM_CLOCK(endtime)
	
	WRITE(*,*) 'elapsed time: ', (endtime-begintime)/1000.0, ' seconds'
	CLOSE (8)

!	CALL Output(Temn, Temx, Teav, osav, pHsav, cmn, cmx, cav, &
!								pHmn, pHmx, pHav, TNmn, TNmx, TNav, &
!								TPmn, TPmx, TPav, NH3mn, NH3mx, NH3av, time1, nj)

	! Body of Q2KFORTRANTEST
	
!gp 12-Jan-06, 21-Nov-06
!	CLOSE (10)
!	CLOSE (11)

	end program Q2KMain




