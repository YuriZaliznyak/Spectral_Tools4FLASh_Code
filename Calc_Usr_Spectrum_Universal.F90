! Routine calculates the Fourier transform (sign=-1, scale=1./(NxNyNz))
! on the variable with the name "VarString"
! and returns real and imaginary parts of the result in the global variables
! with the names "ReSpString" and "ImSpString" (these variables should be
! defined in the Config file and compiled into the global datastructure!).
! Universal: works for any n(xyz)b at any number of CPUs, employes MPI.
! Limitations: leaf nodes should belong to the same refinement level, works only in
! three dimensions, uses library "BLACS" and native IBM's libraries ("ESSL", "pESSL").
!
! Written by Yu. Zaliznyak for FLASH hydro/mhd code,
! tested with FLASH2.4 at IBM p690 in IPP Garching.
!
! 21 Dec. 2004, Garching bei Munchen.
! Yuriy Zaliznyak, mailto: zalik AT kinr.kiev.ua
!
! Updated 15 Jun. 2005. Bug in accessor definitions.

SUBROUTINE Calc_Usr_Spectrum_Universal (VarString, ReSpString, ImSpString)

USE perfmon
USE dBase, ONLY: &
       ndim, nxb, nyb, nzb, nguard, maxblocks, &
       k2d, k3d, &
       dBaseKeyNumber,&
       dBaseGetCoords, &
       dBasePropertyInteger, &
       GC,                             &
       dBaseGetDataPtrSingleBlock, &
       dBaseReleaseDataPtrSingleBlock, &
       dBaseNodeType
USE runtime_parameters, ONLY: get_parm_from_context, GLOBAL_PARM_CONTEXT

IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER _REINDEX3DTO1D_, _HowMuchOf_A_IsNotIn_B_
EXTERNAL _REINDEX3DTO1D_, _HowMuchOf_A_IsNotIn_B_

Character*4 :: VarString, ReSpString, ImSpString

INTEGER, PARAMETER :: q = max(nxb+2*nguard, nyb+2*nguard*k2d, nzb+2*nguard*k3d)
REAL, DIMENSION(q) :: XCoordArr, YCoordArr, ZCoordArr

INTEGER       :: IDMasterCPU, IDCurrCPU, iloop, &
		 iLocNumBlocks, iLocNumLeafNodes, iTotNumLeafNodes, &
		 iNumBlocksXDir, iNumBlocksYDir, iNumBlocksZDir, ix, iy, iz, &
		 ix1, iy1, iz1, iLocNumPlanes, iNumBricks4Fourier, &
		 iNumNeedBlocks, iTotNumBlocks2Recieve, iCPU1, iCPU2, &
		 iNumBlocks4CurrSend, iCurrCopiedBlockNumber, iNF1, iNF2, &
		 iBlockNumber

INTEGER, ALLOCATABLE :: iLocalBlockIndex(:), iGlobalBlockIndex1D(:), &
                        iRecVCounts(:), iDispls(:), iGlobalBlockIndex(:,:), &
			iPlaneIndexArr(:), iHaveBlocks(:), iNeedBlocks(:), &
			iTempArr1(:), iTempArr2(:), iTempArr3(:), &
			ALLRecieveIndex(:), iGlobalPlaneDistributionArr(:,:), &
			iCPUPlaneDistributionArr(:,:)

REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: ALLRecieveData, CurrTransmitData
REAL, DIMENSION(:,:,:,:), POINTER :: ALLBlockData

COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: CIn, COut

INTEGER, SAVE :: iXcoord, iYcoord, iZcoord, izn, &
                 NumCPUs
INTEGER :: iVar, iSpectrumRE, iSpectrumIM

INTEGER :: Spectr_Timer
COMMON /SpectrTimerData/ Spectr_Timer

REAL,    SAVE :: XMin, XMax, YMin, YMax, ZMin, ZMax, SizeX, SizeY, SizeZ

LOGICAL, SAVE :: FirstCall = .TRUE.

INTEGER  :: iTemp1, iTemp2, iTemp3

INTEGER  :: IERROR, MPI_STATUS(MPI_STATUS_SIZE), MessTAG1, MessTAG2

REAL      :: Scale4FT
INTEGER   :: ICONTEXT, IP(40)

CHARACTER*4       :: CPUNumString
CHARACTER*15      :: Filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MessTAG1 = 11; MessTAG2=19
iVar        = dBaseKeyNumber(VarString)
iSpectrumRE = dBaseKeyNumber(ReSpString)
iSpectrumIM = dBaseKeyNumber(ImSpString)
IDMasterCPU = dBasePropertyInteger("MasterProcessor")
IDCurrCPU   = dBasePropertyInteger("MyProcessor")
CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

IF (FirstCall) THEN
   FirstCall = .FALSE.
   IF (IDMasterCPU.NE.0) THEN !!? check :: is it important or not?
     Print *, 'FLASH aborted in Calc_Usr_Spectrum_Universal: '
     Print *, 'ID of Master CPU is not equal to zero! '
     CALL Abort_FLASH()
   EndIF
   CALL Timer_Create("SpectralCalc", Spectr_Timer)
   iXcoord     =  dBaseKeyNumber("xCoord")
   iYcoord     =  dBaseKeyNumber("yCoord")
   iZcoord     =  dBaseKeyNumber("zCoord")
   izn         =  dBaseKeyNumber('zn')
   NumCPUs = dBasePropertyInteger("NumberOfProcessors")
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'xmin', XMin)
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'xmax', XMax)
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'ymin', YMin)
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'ymax', YMax)
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'zmin', ZMin)
   CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'zmax', ZMax)
   SizeX = ABS(XMax-XMin); SizeY = ABS(YMax-YMin); SizeZ = ABS(ZMax-ZMin)
EndIF !! FirstCall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL Timer_START (Spectr_Timer) ! Spectral calculations begin !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!############################################################################
!!!###!!! Define numbers of leaf-blocks in x, y and z directions
iTemp1 = -1
iloop  =  0
 Do While (iTemp1.EQ.-1)
    iloop = iloop + 1
    IF (dBaseNodeType(iloop) .EQ. 1) THEN !! WARNING: all leaf nodes MUST belong !!
      iTemp1 = 0                          !! to the same refinement level        !!
      CALL dBaseGetCoords(izn, iXcoord, iloop, XCoordArr)
      CALL dBaseGetCoords(izn, iYcoord, iloop, YCoordArr)
      CALL dBaseGetCoords(izn, iZcoord, iloop, ZCoordArr)
      iNumBlocksXDir = NINT( (SizeX)/ABS(XCoordArr(nxb + nguard) - XCoordArr(nguard)) )
      iNumBlocksYDir = NINT( (SizeY)/ABS(YCoordArr(nyb + nguard) - YCoordArr(nguard)) )
      iNumBlocksZDir = NINT( (SizeZ)/ABS(ZCoordArr(nzb + nguard) - ZCoordArr(nguard)) )
    EndIF
 EndDo
!!!###!!!
!############################################################################

iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")

!============================================================================
!!!===!!! Calculate the amount of leaf-nodes at current process
iLocNumLeafNodes = 0
 Do iloop = 1, iLocNumBlocks
     IF (dBaseNodeType(iloop).EQ.1) &  !! the leaf node
	    iLocNumLeafNodes = iLocNumLeafNodes + 1
 EndDo
CALL MPI_REDUCE(iLocNumLeafNodes, iTotNumLeafNodes, 1, MPI_INTEGER, &
 	                MPI_SUM, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(iTotNumLeafNodes, 1, MPI_INTEGER, IDMasterCPU, &
	                MPI_COMM_WORLD, IERROR)
!!!===!!!
!============================================================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%!!! Create index for local leaf-blocks and enumerate iHaveBlocks...
ALLOCATE(iLocalBlockIndex(iLocNumLeafNodes*3)) ! 3:: 1D_index, CPUNum, LocalBlockNumber
ALLOCATE(iHaveBlocks(iLocNumLeafNodes)) ! 1D block numbers for local leaf nodes
iTemp1 = -2 ! counter for "iLocalBlockIndex"
iTemp2 =  0 ! counter for "iHaveBlocks"
 Do iloop = 1, iLocNumBlocks
     IF (dBaseNodeType(iloop).EQ.1) THEN  !! the leaf node
            iTemp1 = iTemp1 + 3
	    iTemp2 = iTemp2 + 1
	    CALL dBaseGetCoords(izn, iXcoord, iloop, XCoordArr)
            CALL dBaseGetCoords(izn, iYcoord, iloop, YCoordArr)
            CALL dBaseGetCoords(izn, iZcoord, iloop, ZCoordArr)
            ix = (XCoordArr(q/2)*iNumBlocksXDir/SizeX) + 1
            iy = (YCoordArr(q/2)*iNumBlocksYDir/SizeY) + 1
            iz = (ZCoordArr(q/2)*iNumBlocksZDir/SizeZ) + 1
            iLocalBlockIndex(iTemp1) = _REINDEX3DTO1D_(iNumBlocksXDir, &
	                                iNumBlocksYDir, ix, iy, iz)
	    iLocalBlockIndex(iTemp1 + 1) = IDCurrCPU
	    iLocalBlockIndex(iTemp1 + 2) = iloop
	    iHaveBlocks(iTemp2) = iLocalBlockIndex(iTemp1)
     EndIF
 EndDo
 CALL ISORT(iHaveBlocks, +1, iLocNumLeafNodes)

 ALLOCATE(iRecVCounts(NumCPUs), iDispls(NumCPUs))
 ALLOCATE(iGlobalBlockIndex1D(iTotNumLeafNodes*3)) ! 3:: 1DNumber, CPUNumber, LocalNumber
 CALL MPI_GATHER  (iLocNumLeafNodes, 1, MPI_INTEGER, iRecVCounts, 1, MPI_INTEGER, &
                  IDMasterCPU, MPI_COMM_WORLD, IERROR)
 CALL MPI_BCAST(iRecVCounts, NumCPUs, MPI_INTEGER, IDMasterCPU, MPI_COMM_WORLD, IERROR)
!for each portion of indices form different CPU, the displacements for MPI_GATHERV are calculated
 iDispls(1) = 0
  Do iloop = 2, NumCPUs
     iDispls(iloop) = iDispls(iloop - 1) + 3*iRecVCounts(iloop - 1)
  EndDo
 CALL MPI_GATHERV (iLocalBlockIndex, iLocNumLeafNodes*3, MPI_INTEGER, iGlobalBlockIndex1D, &
                  iRecVCounts*3, iDispls, MPI_INTEGER, IDMasterCPU, MPI_COMM_WORLD, IERROR)
 DEALLOCATE(iDispls, iRecVCounts, iLocalBlockIndex)
! Now, MasterCPU has 1D array of global block distribution, which is to be rearranged and sorted
! Each process allocates 2D global block distribution array...
 ALLOCATE(iGlobalBlockIndex(iTotNumLeafNodes, 3))
! Mater CPU makes arrangements
 IF (IDCurrCPU.EQ.IDMasterCPU) THEN
   iTemp1 = 1
   Do iloop = 1, iTotNumLeafNodes
      iGlobalBlockIndex(iloop, 1) = iGlobalBlockIndex1D(iTemp1)
      iGlobalBlockIndex(iloop, 2) = iGlobalBlockIndex1D(iTemp1 + 1)
      iGlobalBlockIndex(iloop, 3) = iGlobalBlockIndex1D(iTemp1 + 2)
      iTemp1 = iTemp1 + 3
   EndDo
   DEALLOCATE(iGlobalBlockIndex1D)
   ALLOCATE(iTempArr1(iTotNumLeafNodes), iTempArr2(iTotNumLeafNodes), iTempArr3(iTotNumLeafNodes))
   iTempArr2(:) = iGlobalBlockIndex(:, 2)
   iTempArr3(:) = iGlobalBlockIndex(:, 3)

   CALL ISORTX(iGlobalBlockIndex(:,1), +1, iTotNumLeafNodes, iTempArr1)

   Do iloop = 1, iTotNumLeafNodes
    iGlobalBlockIndex(iloop, 2) = iTempArr2(iTempArr1(iloop))
    iGlobalBlockIndex(iloop, 3) = iTempArr3(iTempArr1(iloop))
   EndDo
   DEALLOCATE(iTempArr1, iTempArr2, iTempArr3)
 EndIF ! master CPU finished sorting of global block distribution array
! and MasterCPU then broadcasts global block distribution map to all CPUs
 CALL MPI_BCAST(iGlobalBlockIndex, 3*iTotNumLeafNodes, MPI_INTEGER, IDMasterCPU, &
                  MPI_COMM_WORLD, IERROR)
!!!%%%!!! index for leaf-blocks was created
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!^^^!!! Calculate amount of planes per CPU and indexate needed blocks
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ALLOCATE(iTempArr1(NumCPUs))
iTempArr1(:) = 0
iTemp1 = 0
iloop = 0
iTemp2 = (iNumBlocksZDir*nzb + NumCPUs - 1)/NumCPUs ! number of planes per CPU, pESSL
Do While (iNumBlocksZDir*nzb-iTemp1.GE.iTemp2)
  iloop = iloop + 1
  iTempArr1(iloop) = iTemp2
  iTemp1 = iTemp1 + iTemp2
EndDo
iTempArr1(iloop + 1) = iNumBlocksZDir*nzb - iTemp1
iNumBricks4Fourier = iloop + 1

iLocNumPlanes = iTempArr1(IDCurrCPU + 1)

!-------
ALLOCATE(iPlaneIndexArr(iLocNumPlanes))
iTemp1 = 0 ! for each process calculate starting plane number as a sum of previous...
Do iloop = 0, IDCUrrCPU - 1
   iTemp1 = iTemp1 + iTempArr1(iloop + 1)
EndDo
Do iloop = 1, iLocNumPlanes ! and put the plane numbers into array "iPlaneIndexArr"
   iPlaneIndexArr(iloop) = iTemp1 + iloop
EndDo
!-------

ALLOCATE(iGlobalPlaneDistributionArr(iNumBlocksZDir*nzb, 2), &
         iCPUPlaneDistributionArr(NumCPUs, 2))
 ! maps data planes distribution between processes ...
iGlobalPlaneDistributionArr(:,:) = 0
iloop = 1 ; iTemp1 = 1 ; iTemp2 = 1
Do While (iloop.LE.iNumBlocksZDir*nzb)
 iGlobalPlaneDistributionArr(iloop, 1) = iTemp2 - 1
 iGlobalPlaneDistributionArr(iloop, 2) = iTemp1
 iloop = iloop + 1
 iTemp1 = iTemp1 + 1
 IF (iTemp1.GT.iTempArr1(iTemp2)) THEN
      iTemp1 = 1
      iTemp2 = iTemp2 + 1
 EndIF
EndDo

Do iloop = 1, NumCPUs
 iTemp2 = 0
   Do iTemp1 = 1, iloop - 1
     iTemp2 = iTemp2 + iTempArr1(iTemp1)
   EndDo
 iCPUPlaneDistributionArr(iloop, 1) = iTemp2 + 1
 iCPUPlaneDistributionArr(iloop, 2) = iTempArr1(iloop)
EndDo

DEALLOCATE(iTempArr1)! iTempArr1 contained numbers of planes per CPU, first element -- CPU0, etc.

IF (iLocNumPlanes.NE.0) THEN
  iTemp1 = iPlaneIndexArr(1)/(nzb + 1e-9) + 1
  iTemp2 = iPlaneIndexArr(iLocNumPlanes)/(nzb + 1e-9) + 1
  iNumNeedBlocks = iNumBlocksXDir*iNumBlocksYDir*(iTemp2 - iTemp1 + 1)
                        ELSE !! iLocNumPlanes = 0
  iTemp1 = 1; iTemp2 = 0; iNumNeedBlocks = 0
EndIF

ALLOCATE(iNeedBlocks(iNumNeedBlocks))
iloop = 0
Do iz = iTemp1, iTemp2
  Do iy = 1, iNumBlocksYDir
    Do ix = 1, iNumBlocksXDir
      iloop = iloop + 1
      iNeedBlocks(iloop) = _REINDEX3DTO1D_(iNumBlocksXDir, &
	                                iNumBlocksYDir, ix, iy, iz)
    EndDo
  EndDo
EndDo
CALL ISORT(iNeedBlocks, +1, iNumNeedBlocks)! not really needed, but ... who knows ?

!!!^^^!!!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!Communication part of the routine -- processes send blocks to each other $
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! Each process allocates space for transmitted data
iTotNumBlocks2Recieve = _HowMuchOf_A_IsNotIn_B_(IDCurrCPU, iNumNeedBlocks, iLocNumLeafNodes, &
                                              iNeedBlocks, iHaveBlocks)

ALLOCATE(ALLRecieveData(iTotNumBlocks2Recieve, nxb, nyb, nzb), &
          ALLRecieveIndex(iTotNumBlocks2Recieve))
iCurrCopiedBlockNumber = 0 ! Number of current block, which will be copied
                           ! into the global data array "ALLRecieveData"
!Processors communicate by pairs
Do iCPU1 = 0, NumCPUs - 2
 Do iCPU2 = iCPU1 + 1, NumCPUs-1
   IF (IDCurrCPU.EQ.iCPU1.OR.IDCurrCPU.EQ.iCPU2) THEN !CPUs iCPU1 and iCPU2 communicate!
!!`````````````````````````````````````````````````````````````````````````````````````````
!! iCPU1 WORKS
  IF (IDCurrCPU.EQ.iCPU1) THEN
  ! iCPU1 counts the number of blocks which will be recieved from iCPU2
     iNumBlocks4CurrSend = 0
     Do iloop = 1, iNumNeedBlocks
        IF (iGlobalBlockIndex(iNeedBlocks(iloop),2).EQ.iCPU2) &
	             iNumBlocks4CurrSend = iNumBlocks4CurrSend + 1
     EndDo
  ! iCPU1 after that calls MPI_SEND to inform iCPU2 about iNumBlocks4CurrSend
     CALL MPI_SEND(iNumBlocks4CurrSend, 1, MPI_INTEGER, iCPU2, MessTAG1, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU1 counted number of blocks to be sent from iCPU2

!! iCPU2 WORKS
  IF (IDCurrCPU.EQ.iCPU2) &
     CALL MPI_RECV(iNumBlocks4CurrSend, 1, MPI_INTEGER, iCPU1, MessTAG1, MPI_COMM_WORLD, &
                   MPI_STATUS, IERROR)

 ! all other stuff is done only if there is something to send...
  IF (iNumBlocks4CurrSend.NE.0) THEN

!! iCPU1 and iCPU2 WORK
  ALLOCATE(iTempArr1(iNumBlocks4CurrSend), &          ! index and data arrays are allocated at
     CurrTransmitData(iNumBlocks4CurrSend, nxb, nyb, nzb)) ! both CPUs (iCPU1 and iCPU2)

!! iCPU1 WORKS
  IF (IDCurrCPU.EQ.iCPU1) THEN
    ! iCPU1 indexates blocks which will be recieved from iCPU2
     iTemp1 = 0
     Do iloop = 1, iNumNeedBlocks
        IF (iGlobalBlockIndex(iNeedBlocks(iloop),2).EQ.iCPU2) THEN
	    iTemp1 = iTemp1 + 1
	    iTempArr1(iTemp1) = iGlobalBlockIndex(iNeedBlocks(iloop), 1)
	EndIF
     EndDo
  ! iCPU1 then MPI_SEND index array "iTempArr1" to iCPU2
     CALL MPI_SEND(iTempArr1, iNumBlocks4CurrSend, MPI_INTEGER, &
                   iCPU2, MessTAG2, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU1 indexated blocks and sent index array "iTempArr1" to iCPU2

!! iCPU2 WORKS
  IF (IDCurrCPU.EQ.iCPU2) THEN
    ! iCPU2 recieves index array from iCPU1
      CALL MPI_RECV(iTempArr1, iNumBlocks4CurrSend, MPI_INTEGER, iCPU1, &
                    MessTAG2, MPI_COMM_WORLD, MPI_STATUS, IERROR)

   ! iCPU2 buffers the blocks with 1D numbers from "iTempArr1" into "CurrTransmitData" array
      Do iloop = 1, iNumBlocks4CurrSend
        iTemp2 = iGlobalBlockIndex(iTempArr1(iloop), 3) ! local number of block
        ALLBlockData => dBaseGetDataPtrSingleBlock(iTemp2, GC)
        CurrTransmitData(iloop, :, :, :) = ALLBlockData(iVar, &
	           nguard+1:nguard+nxb,nguard+1:nguard+nyb,nguard+1:nguard+nzb)
	CALL dBaseReleaseDataPtrSingleBlock(iTemp2, AllBlockData)
      EndDo

   ! iCPU2 sends buffered block to iCPU1
      CALL MPI_SEND(CurrTransmitData, iNumBlocks4CurrSend*nxb*nyb*nzb, MPI_DOUBLE_PRECISION, &
                      iCPU1, MessTAG1, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU2 buffered the data into "CurrTransmitData" array and sent it to iCPU1

!! iCPU1 WORKS
  IF (IDCurrCPU.EQ.iCPU1) THEN
    ! iCPU1 recieves data array "CurrTransmitData" from iCPU2
      CALL MPI_RECV(CurrTransmitData, iNumBlocks4CurrSend*nxb*nyb*nzb, MPI_DOUBLE_PRECISION, &
                      iCPU2, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)

    ! and stores recieved data in the array "ALLRecieveData"
      Do iloop = 1, iNumBlocks4CurrSend
        iCurrCopiedBlockNumber = iCurrCopiedBlockNumber + 1
	ALLRecieveData(iCurrCopiedBlockNumber, :, :, :) = CurrTransmitData(iloop, :, :, :)
	ALLRecieveIndex(iCurrCopiedBlockNumber) = iTempArr1(iloop)
      EndDo
  EndIF

  DEALLOCATE(iTempArr1, CurrTransmitData) ! iTempArr1 contained transmit block
                                              ! index for current communication iCPU1--iCPU2
 EndIF ! There was something to send...
!!`````````````````````````````````````````````````````````````````````````````````````

!!!===========================================================================
! The "mirrored iCPU1 <--> iCPU2" communication: now iCPU1 sends needed blocks
! needed blocks _to_ iCPU2
!!!===========================================================================
!!`````````````````````````````````````````````````````````````````````````````````````
!! iCPU2 WORKS
  IF (IDCurrCPU.EQ.iCPU2) THEN
  ! iCPU2 counts the number of blocks which will be recieved from iCPU1
     iNumBlocks4CurrSend = 0
     Do iloop = 1, iNumNeedBlocks
        IF (iGlobalBlockIndex(iNeedBlocks(iloop),2).EQ.iCPU1) &
	             iNumBlocks4CurrSend = iNumBlocks4CurrSend + 1
     EndDo
  ! iCPU2 after that calls MPI_SEND to inform iCPU2 about iNumBlocks4CurrSend
     CALL MPI_SEND(iNumBlocks4CurrSend, 1, MPI_INTEGER, iCPU1, MessTAG1, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU2 counted number of blocks to be sent from iCPU1

!! iCPU1 WORKS
  IF (IDCurrCPU.EQ.iCPU1) &
     CALL MPI_RECV(iNumBlocks4CurrSend, 1, MPI_INTEGER, iCPU2, MessTAG1, MPI_COMM_WORLD, &
                   MPI_STATUS, IERROR)

! all other stuff is done only if there is something to send...
  IF (iNumBlocks4CurrSend.NE.0) THEN

!! iCPU1 and iCPU2 WORK
  ALLOCATE(iTempArr1(iNumBlocks4CurrSend), &          ! index and data arrays are allocated at
     CurrTransmitData(iNumBlocks4CurrSend, nxb, nyb, nzb)) ! both CPUs (iCPU1 and iCPU2)

!! iCPU2 WORKS
  IF (IDCurrCPU.EQ.iCPU2) THEN
    ! iCPU2 indexates blocks which will be recieved from iCPU1
     iTemp1 = 0
     Do iloop = 1, iNumNeedBlocks
        IF (iGlobalBlockIndex(iNeedBlocks(iloop),2).EQ.iCPU1) THEN
	    iTemp1 = iTemp1 + 1
	    iTempArr1(iTemp1) = iGlobalBlockIndex(iNeedBlocks(iloop), 1)
	EndIF
     EndDo
  ! iCPU2 then MPI_SEND index array "iTempArr1" to iCPU1
     CALL MPI_SEND(iTempArr1, iNumBlocks4CurrSend, MPI_INTEGER, &
                   iCPU1, MessTAG2, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU2 indexated blocks and sent index array "iTempArr1" to iCPU1

!! iCPU1 WORKS
  IF (IDCurrCPU.EQ.iCPU1) THEN
    ! iCPU1 recieves index array from iCPU2
      CALL MPI_RECV(iTempArr1, iNumBlocks4CurrSend, MPI_INTEGER, iCPU2, &
                    MessTAG2, MPI_COMM_WORLD, MPI_STATUS, IERROR)

   ! iCPU1 buffers the blocks with 1D numbers from "iTempArr1" into "CurrTransmitData" array
      Do iloop = 1, iNumBlocks4CurrSend
        iTemp2 = iGlobalBlockIndex(iTempArr1(iloop), 3) ! local number of block
        ALLBlockData => dBaseGetDataPtrSingleBlock(iTemp2, GC)
        CurrTransmitData(iloop, :, :, :) = ALLBlockData(iVar, &
	           nguard+1:nguard+nxb,nguard+1:nguard+nyb,nguard+1:nguard+nzb)
	CALL dBaseReleaseDataPtrSingleBlock(iTemp2, AllBlockData)
      EndDo

   ! iCPU1 sends buffered block to iCPU2
      CALL MPI_SEND(CurrTransmitData, iNumBlocks4CurrSend*nxb*nyb*nzb, MPI_DOUBLE_PRECISION, &
                      iCPU2, MessTAG1, MPI_COMM_WORLD, IERROR)
  EndIF ! iCPU1 buffered the data into "CurrTransmitData" array and sent it to iCPU2

!! iCPU2 WORKS
  IF (IDCurrCPU.EQ.iCPU2) THEN
    ! iCPU2 recieves data array "CurrTransmitData" from iCPU1
      CALL MPI_RECV(CurrTransmitData, iNumBlocks4CurrSend*nxb*nyb*nzb, MPI_DOUBLE_PRECISION, &
                      iCPU1, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)

    ! and stores recieved data in the array "ALLRecieveData"
      Do iloop = 1, iNumBlocks4CurrSend
        iCurrCopiedBlockNumber = iCurrCopiedBlockNumber + 1
	ALLRecieveData(iCurrCopiedBlockNumber, :, :, :) = CurrTransmitData(iloop, :, :, :)
	ALLRecieveIndex(iCurrCopiedBlockNumber) = iTempArr1(iloop)
      EndDo
  EndIF

  DEALLOCATE(iTempArr1, CurrTransmitData) ! iTempArr1 contained transmit block
                                              ! index for current communication iCPU1--iCPU2

 EndIF ! There was something to send...

 EndIF ! iCPU1 and iCPU2 finished communications
 CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR) ! other CPUs are waiting
                                          ! until iCPU1 and iCPU2 finish talking
 EndDo ! iCPU2
EndDo ! iCPU1
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!! At this point, communication are finished, all processes have all
!!! data needed _to_ create local arrays for Fourier transform
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
iNF1 = iNumBlocksXDir*nxb
iNF2 = iNumBlocksXDir*nyb

ALLOCATE(CIn(iNF1, iNF2, iLocNumPlanes))

Do iloop = 1, iNumNeedBlocks
  iTemp1 = iNeedBlocks(iloop) ! 1D number of current block
  IF (iGlobalBlockIndex(iTemp1, 2).EQ.IDCurrCPU) THEN ! block is avaliable locally

     iTemp2 = iGlobalBlockIndex(iTemp1, 3) ! local number of current block
     ALLBlockData => dBaseGetDataPtrSingleBlock(iTemp2, GC)
     Do ix = 1, nxb
      Do iy = 1, nyb
       Do iz = 1, nzb ! can be improved by checking once _for_ current "iz"...
	 CALL _REINDEX2FOURIER_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
	                     iLocNumPlanes, &
                             iPlaneIndexArr, iTemp1, ix, iy, iz, &
			     iTemp3, ix1, iy1, iz1)
	 IF (ix1.GT.0) & ! reindexator returns ix=-1, _if_ point is NOT used in brick
	   CIn(ix1, iy1, iz1) = CMPLX(ALLBlockData(iVar, ix+nguard, iy+nguard, iz+nguard), 0.)
       EndDo
      EndDo
     EndDo
     CALL dBaseReleaseDataPtrSingleBlock(iTemp2, AllBlockData)

						ELSE ! NOT avaliable locally
     iTemp2 = -1
     iTemp3 = 0
     Do While (iTemp2.EQ.-1)
       iTemp3 = iTemp3 + 1
       	IF (ALLRecieveIndex(iTemp3).EQ.iTemp1) &
	       iTemp2 = iTemp3 ! iTemp2 now contains index of needed block
	                       ! in "ALLRecieveData" array
     EndDo

     Do ix = 1, nxb
      Do iy = 1, nyb
       Do iz = 1, nzb ! can be improved by checking once _for_ current "iz"...
	 CALL _REINDEX2FOURIER_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
	                     iLocNumPlanes, &
                             iPlaneIndexArr, iTemp1, ix, iy, iz, &
			     iTemp3, ix1, iy1, iz1)
	 IF (ix1.GT.0) & ! reindexator returns ix=-1, _if_ point is NOT used in brick
	   CIn(ix1, iy1, iz1) = CMPLX(ALLRecieveData(iTemp2, ix, iy, iz), 0.)
       EndDo
      EndDo
     EndDo
  EndIF

EndDo !iloop! All blocks from "iNeedBlocks" array were processed
! Now, array "ALLRecieveData" is not needed anymore, and can be deallocated
DEALLOCATE(ALLRecieveData, ALLRecieveIndex)
DEALLOCATE(iPlaneIndexArr)

CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!!! Finally the Fourier transform routine is called...
ALLOCATE(COut(iNF1, iNF2, iLocNumPlanes))
CALL BLACS_GET(0, 0, ICONTEXT)
CALL BLACS_GRIDINIT(ICONTEXT, 'R', 1, NumCPUs)
Scale4FT = 1.0/(iNF1*iNF2*iNumBlocksZDir*nzb)
IP(:)  = 0; IP(1)  = 1; IP(2)  = 1
CALL PDCFT3(CIn, COut, iNF1, iNF2, iNumBlocksZDir*nzb, +1, Scale4FT, ICONTEXT, IP)
CALL BLACS_GRIDEXIT(ICONTEXT)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

IF (1.eq.10) THEN !! some test output -- may be helpful...
      Write(CPUNumString, 1007) 1000 + IDCurrCPU
1007	FORMAT(i4)
      Filename = 'CPUDATA'//CPUNumString//'.dat'
      Open (311, File = Filename)
      Do iz = 1, iLocNumPlanes
       Do iy = 1, iNF2
        Do ix = 1, iNF1

          Write(311,*)REAL(CIn(ix, iy, iz)), AIMAG(CIn(ix, iy, iz)), &
	              REAL(COut(ix, iy, iz)), AIMAG(COut(ix, iy, iz))
	EndDo
       EndDo
      EndDo
      Close(311)
   Print *, 'CPU::', IDCUrrCPU, ' has finished to write::', Filename, &
            '  num. planes::', iLocNumPlanes
EndIF

DEALLOCATE(CIn)
!!!%%Fourier transform is done and result is stored in local arrays COut
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!------------------------------------------------------------------------
!!!----------------------------------------------------------------------
! As the next step, back data redistribution is done _to_ fill local
! "resp" and "imsp" blocks
!!!----------------------------------------------------------------------
!------------------------------------------------------------------------
!Processors communicate by pairs: iCPU1 <--> iCPU2
Do iCPU1 = 0, NumCPUs - 2
 Do iCPU2 = iCPU1 + 1, NumCPUs-1
   IF (IDCurrCPU.EQ.iCPU1.OR.IDCurrCPU.EQ.iCPU2) THEN !CPUs iCPU1 and iCPU2 communicate!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ALLOCATE(iTempArr3(3*MAXBLOCKS*iNumBlocksZDir*nzb))
 IF (IDCurrCPU.EQ.iCPU1) THEN !! 1:iCPU1 works
   ! iCPU1 checks all local blocks whether they intersect planes from iCPU2
     iTemp3 = -2 !-- current position in send-recieve index array "iTempArr3"
   ! _for_ each local block, all planes from iCPU2 are tested
     Do iloop = 1, iLocNumLeafNodes ! iloop -- current number of local block
      Do iTemp1 = iCPUPlaneDistributionArr(iCPU2 + 1, 1), & ! iTemp1 -- number of plane being tested
          iCPUPlaneDistributionArr(iCPU2 + 1, 1) + iCPUPlaneDistributionArr(iCPU2 + 1, 2) - 1
         CALL _BlockPlaneIntersection_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
	                              iHaveBlocks(iloop), iTemp1, ix, iy)
	IF (ix.GT.0) THEN ! there IS intersection of block "iloop" and plane "iTemp1"
           iTemp3 = iTemp3 + 3
           iTempArr3(iTemp3)     = iTemp1
	   iTempArr3(iTemp3 + 1) = ix
	   iTempArr3(iTemp3 + 2) = iy
	EndIF
      EndDo !iloop
     EndDo !iTemp2
     CALL MPI_SEND(iTemp3, 1, MPI_INTEGER, iCPU2, MessTAG1, MPI_COMM_WORLD, IERROR)
   EndIF !! 1:iCPU1 works

 IF (IDCurrCPU.EQ.iCPU2) & ! iCPU2 recieves SIZE of index array
     CALL MPI_RECV(iTemp3, 1, MPI_INTEGER, iCPU1, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)

!________
 IF (iTemp3.GT.0) THEN
!________All other stuff is done ONLY_IF_THERE_IS something to send from iCPU2 to iCPU1
 ALLOCATE(CurrTransmitData((iTemp3 + 2)/3, nxb, nyb, 2))! at both CPUs

 IF (IDCurrCPU.EQ.iCPU1) &
          CALL MPI_SEND(iTempArr3, iTemp3+2, MPI_INTEGER, iCPU2, MessTAG2, MPI_COMM_WORLD, IERROR)

 IF (IDCurrCPU.EQ.iCPU2) THEN !! 2:iCPU2 works
   ! iCPU2 recieves index array
       CALL MPI_RECV(iTempArr3, iTemp3+2, MPI_INTEGER, iCPU1, MessTAG2, &
                   MPI_COMM_WORLD, MPI_STATUS, IERROR)

   ! iCPU2 buffers data into the array "CurrTransmitData"
     iTemp1 = 0 ! current position in "CurrTransmitData" array
     Do iloop = 1, iTemp3, 3 !! iloop -- current position in "iTempArr3" array
       iTemp1 = iTemp1 + 1
       iz  = iTempArr3(iloop) ! plane number in global units
       iz1 = iGlobalPlaneDistributionArr(iz, 2) ! local plane number in the brick
       ix1 = iTempArr3(iloop + 1)
       iy1 = iTempArr3(iloop + 2)
       Do ix = 1, nxb
        Do iy = 1, nyb
         CurrTransmitData(iTemp1, ix, iy, 1) = REAL  (COut(ix1+ix-1, iy1+iy-1, iz1))
         CurrTransmitData(iTemp1, ix, iy, 2) = AIMAG (COut(ix1+ix-1, iy1+iy-1, iz1))
	EndDo
       EndDo
     EndDo !iloop
     CALL MPI_SEND(CurrTransmitData,  ((iTemp3 + 2)/3)*nxb*nyb*2, MPI_DOUBLE_PRECISION, &
                   iCPU1, MessTAG1, MPI_COMM_WORLD, IERROR)
 EndIF !! 2:iCPU2

 IF (IDCurrCPU.EQ.iCPU1) THEN !! 3:iCPU1 works
     CALL MPI_RECV(CurrTransmitData,  ((iTemp3 + 2)/3)*nxb*nyb*2, MPI_DOUBLE_PRECISION, &
                   iCPU2, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)
     iTemp1 = 0 ! current position in "CurrTransmitData" array
     CALL _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, iTempArr3(2), &
                          iTempArr3(3), iTempArr3(1), iBlockNumber, iz1)
     iBlockNumber = iGlobalBlockIndex(iBlockNumber, 3)
     ALLBlockData => dBaseGetDataPtrSingleBlock(iBlockNumber, GC)
     Do iloop = 1, iTemp3, 3 !! iloop -- current position in "iTempArr3" array
       iTemp1 = iTemp1 + 1 ! current position in "CurrTransmitData" array
       iz = iTempArr3(iloop) ! plane number (iz-coordinate) in global units
       ix = iTempArr3(iloop + 1) ! ix-coordinate of block start corner, in global units
       iy = iTempArr3(iloop + 2) ! iy-  ---||---
       CALL _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, ix, iy, iz, &
                             iTemp2, iz1)
       iTemp2 = iGlobalBlockIndex(iTemp2, 3)
       IF (iTemp2.NE.iBlockNumber) THEN
          CALL dBaseReleaseDataPtrSingleBlock(iBlockNumber, ALLBlockData)
          iBlockNumber = iTemp2
	  iTemp2 = 0
          ALLBlockData => dBaseGetDataPtrSingleBlock(iBlockNumber, GC)
       EndIF

       Do ix = 1, nxb
       Do iy = 1, nyb
         ALLBlockData(iSpectrumRE, nguard+ix,nguard+iy, nguard+iz1) = &
           CurrTransmitData(iTemp1, ix, iy, 1)
         ALLBlockData(iSpectrumIM, nguard+ix,nguard+iy, nguard+iz1) = &
           CurrTransmitData(iTemp1, ix, iy, 2)
       EndDo
       EndDo
     EndDo !iloop
     CALL dBaseReleaseDataPtrSingleBlock(iBlockNumber, ALLBlockData)
 EndIF!! 3:iCPU1 works

 DEALLOCATE(CurrTransmitData)
!________There WAS something to send from iCPU2 to iCPU1
 EndIF
!________
 DEALLOCATE(iTempArr3)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!@@@@@@@@@ iCPU1 recieved the data from iCPU2. Now it's time _for_ iCPU2 @@@
!!!!!@@@@@@@@@ to do the same: recieve data from iCPU1. As usual, almost     @@@
!!!!!@@@@@@@@@ mirrored piece of code.                                       @@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 ALLOCATE(iTempArr3(3*MAXBLOCKS*iNumBlocksZDir*nzb))
 IF (IDCurrCPU.EQ.iCPU2) THEN !! 4:iCPU2 works
    ! iCPU2 checks all local blocks whether they intersect planes from iCPU1
     iTemp3 = -2 !-- current position in send-recieve index array "iTempArr3"
     Do iloop = 1, iLocNumLeafNodes ! iloop -- current number of local block
    ! _for_ each local block, all planes from iCPU2 are tested
      Do iTemp1 = iCPUPlaneDistributionArr(iCPU1 + 1, 1), & ! iTemp1 -- number of plane being tested
          iCPUPlaneDistributionArr(iCPU1 + 1, 1) + iCPUPlaneDistributionArr(iCPU1 + 1, 2) - 1
	CALL _BlockPlaneIntersection_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
	                              iHaveBlocks(iloop), iTemp1, ix, iy)
	IF (ix.GT.0) THEN ! there IS intersection of block "iloop" and plane "iTemp1"
           iTemp3 = iTemp3 + 3
           iTempArr3(iTemp3)     = iTemp1
	   iTempArr3(iTemp3 + 1) = ix
	   iTempArr3(iTemp3 + 2) = iy
	EndIF
      EndDo !iloop
     EndDo !iTemp1
     CALL MPI_SEND(iTemp3, 1, MPI_INTEGER, iCPU1, MessTAG1, MPI_COMM_WORLD, IERROR)
   EndIF !! 4:iCPU2 works

 IF (IDCurrCPU.EQ.iCPU1) & ! iCPU1 recieves SIZE of index array
     CALL MPI_RECV(iTemp3, 1, MPI_INTEGER, iCPU2, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)

!________
 IF (iTemp3.GT.0) THEN
!________All other stuff is done ONLY_IF_THERE_IS something to send from iCPU1 to iCPU2
 ALLOCATE(CurrTransmitData((iTemp3 + 2)/3, nxb, nyb, 2))! at both CPUs

 IF (IDCurrCPU.EQ.iCPU2) &
          CALL MPI_SEND(iTempArr3, iTemp3+2, MPI_INTEGER, iCPU1, MessTAG2, MPI_COMM_WORLD, IERROR)

 IF (IDCurrCPU.EQ.iCPU1) THEN !! 5:iCPU1 works
   ! iCPU1 recieves index array
     CALL MPI_RECV(iTempArr3, iTemp3+2, MPI_INTEGER, iCPU2, MessTAG2, &
                   MPI_COMM_WORLD, MPI_STATUS, IERROR)

   ! iCPU1 buffers data into the array "CurrTransmitData"
     iTemp1 = 0 ! current position in "CurrTransmitData" array
     Do iloop = 1, iTemp3, 3 !! iloop -- current position in "iTempArr3" array

       iTemp1 = iTemp1 + 1
       iz  = iTempArr3(iloop) ! plane number in global units
       iz1 = iGlobalPlaneDistributionArr(iz, 2) ! local plane number in the brick
       ix1 = iTempArr3(iloop + 1)
       iy1 = iTempArr3(iloop + 2)
       Do ix = 1, nxb
        Do iy = 1, nyb
         CurrTransmitData(iTemp1, ix, iy, 1) = REAL  (COut(ix1+ix-1, iy1+iy-1, iz1))
         CurrTransmitData(iTemp1, ix, iy, 2) = AIMAG (COut(ix1+ix-1, iy1+iy-1, iz1))
	EndDo
       EndDo
     EndDo !iloop
     CALL MPI_SEND(CurrTransmitData,  ((iTemp3 + 2)/3)*nxb*nyb*2, MPI_DOUBLE_PRECISION, &
                   iCPU2, MessTAG1, MPI_COMM_WORLD, IERROR)
 EndIF !! 5:iCPU1

 IF (IDCurrCPU.EQ.iCPU2) THEN !! 6:iCPU2 works
     CALL MPI_RECV(CurrTransmitData,  ((iTemp3 + 2)/3)*nxb*nyb*2, MPI_DOUBLE_PRECISION, &
                   iCPU1, MessTAG1, MPI_COMM_WORLD, MPI_STATUS, IERROR)
     iTemp1 = 0 ! current position in "CurrTransmitData" array
     CALL _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, iTempArr3(2), &
                          iTempArr3(3), iTempArr3(1), iBlockNumber, iz1)
     iBlockNumber = iGlobalBlockIndex(iBlockNumber, 3)
     ALLBlockData => dBaseGetDataPtrSingleBlock(iBlockNumber, GC)
     Do iloop = 1, iTemp3, 3 !! iloop -- current position in "iTempArr3" array
       iTemp1 = iTemp1 + 1 ! current position in "CurrTransmitData" array
       iz = iTempArr3(iloop) ! plane number (iz-coordinate) in global units
       ix = iTempArr3(iloop + 1) ! ix-coordinate of block start corner, in global units
       iy = iTempArr3(iloop + 2) ! iy-  ---||---
       CALL _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, ix, iy, iz, &
                             iTemp2, iz1)
       iTemp2 = iGlobalBlockIndex(iTemp2, 3)
       IF (iTemp2.NE.iBlockNumber) THEN
          CALL dBaseReleaseDataPtrSingleBlock(iBlockNumber, ALLBlockData)
          iBlockNumber = iTemp2
          ALLBlockData => dBaseGetDataPtrSingleBlock(iBlockNumber, GC)
       EndIF
       Do ix = 1, nxb
       Do iy = 1, nyb
         ALLBlockData(iSpectrumRE, nguard+ix,nguard+iy, nguard+iz1) = &
           CurrTransmitData(iTemp1, ix, iy, 1)
         ALLBlockData(iSpectrumIM, nguard+ix,nguard+iy, nguard+iz1) = &
           CurrTransmitData(iTemp1, ix, iy, 2)
       EndDo
       EndDo
     EndDo !iloop
     CALL dBaseReleaseDataPtrSingleBlock(iBlockNumber, ALLBlockData)
 EndIF!! 6:iCPU2 works

 DEALLOCATE(CurrTransmitData)
!________There WAS something to send from iCPU1 to iCPU2
 EndIF
!________
 DEALLOCATE(iTempArr3)

 EndIF ! iCPU1 and iCPU2 finished communication

 CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)! other processors wait until those two communicate

 EndDo ! loop for iCPU1
EndDo ! loop for iCPU2

!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! End of communication part of the routine $$$$$$$$$$$$$$$$$$$$$$$$$$
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!--------------------------------------------------------------------------
! Now each CPU will fill the block data from the planes avaliable locally
!--------------------------------------------------------------------------
    ! current CPU checks all local blocks whether they intersect local planes
      Do iloop = 1, iLocNumLeafNodes ! iloop -- current number of local block
   ! _for_ each local block, all local planes from Current CPU are tested
       Do iTemp1 = iCPUPlaneDistributionArr(IDCurrCPU + 1, 1), & ! iTemp1 -- number of plane being tested
          iCPUPlaneDistributionArr(IDCurrCPU + 1, 1) + iCPUPlaneDistributionArr(IDCurrCPU + 1, 2) - 1

	CALL _BlockPlaneIntersection_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
	                              iHaveBlocks(iloop), iTemp1, ix, iy)
	IF (ix.GT.0) THEN ! there IS intersection of block "iloop" and plane "iTemp1"
         iz = iTemp1 ! plane number (iz-coordinate) in global units
	 iTemp3 = iGlobalPlaneDistributionArr(iTemp1, 2) ! local plane number in COut array
         CALL _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, ix, iy, iz, &
                             iTemp2, iz1)
         iBlockNumber = iGlobalBlockIndex(iTemp2, 3)

         ALLBlockData => dBaseGetDataPtrSingleBlock(iBlockNumber, GC)! some over-use of dBase calls...
	 Do ix1 = 1, nxb
	  Do iy1 = 1, nyb
	   ALLBlockData(iSpectrumRE, nguard+ix1, nguard+iy1, nguard+iz1) = &
             REAL(COut(ix+ix1-1, iy+iy1-1, iTemp3))
	   ALLBlockData(iSpectrumIM, nguard+ix1, nguard+iy1, nguard+iz1) = &
             AIMAG(COut(ix+ix1-1, iy+iy1-1, iTemp3))
	  EndDo !iy1
	 EndDo !ix1
         CALL dBaseReleaseDataPtrSingleBlock(iBlockNumber, ALLBlockData)
	EndIF
      EndDo !iTemp1
     EndDo !iloop

DEALLOCATE(COut)
DEALLOCATE(iGlobalBlockIndex)
DEALLOCATE(iHaveBlocks)
DEALLOCATE(iNeedBlocks)
DEALLOCATE(iGlobalPlaneDistributionArr, iCPUPlaneDistributionArr)

CALL Timer_STOP (Spectr_Timer)

CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION _REINDEX3DTO1D_(iNBX, iNBY, ix, iy, iz)
IMPLICIT NONE !iNBZ and iNBY -- numbers of blocks in X and Y directions
INTEGER :: iNBX, iNBY, ix, iy, iz
_REINDEX3DTO1D_ = (iz-1)*iNBX*iNBY + (iy-1)*iNBX + ix
RETURN
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function finds how much of elements of array iA are not presented in array iB.
!Integer arrays iA and iB MUST be sorted accending without duplicate entries!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION _HowMuchOf_A_IsNotIn_B_(IDCUrrCPU, iNA, iNB, iA, iB)
IMPLICIT NONE
INTEGER, INTENT(IN) :: IDCUrrCPU, iNA, iNB
INTEGER, DIMENSION(iNA), INTENT(IN) :: iA
INTEGER, DIMENSION(iNB), INTENT(IN) :: iB
INTEGER :: iPosA, iPosB, iNum

iPosA = 1
iPosB = 1
iNum = 0
Do While (iPosA.LE.iNA.AND.iPosB.LE.iNB+1)
  IF (iA(iPosA).EQ.iB(iPosB)) THEN
       iPosA = iPosA + 1
       iPosB = iPosB + 1
    ELSE IF (iA(iPosA).LT.iB(iPosB)) THEN
       iNum = iNum + 1
       iPosA = iPosA + 1
    ELSE IF (iA(iPosA).GT.iB(iPosB)) THEN
       iPosB = iPosB + 1
  EndIF
EndDo
IF (iPosA.LT.iNA) iNum = iNum + iNA - iPosA + 1

_HowMuchOf_A_IsNotIn_B_ = iNum

RETURN
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE _REINDEX2FOURIER_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, iLocNumPlanes, &
                             iPlaneIndexArr, i1DBlockNum, ix, iy, iz, &
			     iNumPlaneStart, ix1, iy1, iz1)

USE dBase, ONLY: nxb, nyb, nzb

IMPLICIT NONE
INTEGER, INTENT(IN)    :: IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, iLocNumPlanes, &
                          iPlaneIndexArr(iLocNumPlanes), i1DBlockNum, ix, iy, iz
INTEGER, INTENT(OUT)   :: iNumPlaneStart, ix1, iy1, iz1

INTEGER :: NPlanes, NLines, NCells

NPlanes = i1DBlockNum/(iNumBlocksXDir*iNumBlocksYDir + 1e-10)

iNumPlaneStart = NPlanes*nzb + 1
IF (iNumPlaneStart+iz-1.LT.iPlaneIndexArr(1) &
    .OR.iNumPlaneStart+iz-1.GT.iPlaneIndexArr(iLocNumPlanes)) THEN
 ! given data point is not used in the filled brick !
    ix1 = -1; iy1 = -1; iz1 = -1
    RETURN
EndIF

NLines = (i1DBlockNum - iNumBlocksXDir*iNumBlocksYDir*NPlanes)&
          /(iNumBlocksXDir + 1e-10)
NCells = i1DBlockNum - iNumBlocksXDir*iNumBlocksYDir*NPlanes &
         - iNumBlocksXDir*NLines - 1

ix1 = ix + NCells*nxb
iy1 = iy + NLines*nyb
iz1 = iz + NPlanes*nzb - iPlaneIndexArr(1) + 1 ! index _for_ "CIn" array

RETURN
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE _BlockPlaneIntersection_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, &
                                    iBlock, iPlane, ix, iy)
! returns coordinates (ix, iy) of a starting point of block-plane intersection region
USE dBase, ONLY: nxb, nyb, nzb
IMPLICIT NONE
INTEGER, INTENT(IN)  :: IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, iBlock, iPlane
INTEGER, INTENT(OUT) :: ix, iy

INTEGER :: NPlanes, NLines, NCells, iNumPlaneStart

NPlanes = iBlock/(iNumBlocksXDir*iNumBlocksYDir + 1e-10)
iNumPlaneStart = NPlanes*nzb + 1

IF (iNumPlaneStart.GT.iPlane.OR.iNumPlaneStart+nzb-1.LT.iPlane) THEN
 !-- plane and block _do_ not intersect
   ix = -1; iy = -1
   RETURN
EndIF

NLines = (iBlock - iNumBlocksXDir*iNumBlocksYDir*NPlanes)&
          /(iNumBlocksXDir + 1e-10)
NCells = iBlock - iNumBlocksXDir*iNumBlocksYDir*NPlanes &
         - iNumBlocksXDir*NLines - 1

ix = NCells*nxb + 1
iy = NLines*nyb + 1

RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! returns 1D block number and iz1 -- z-index inside that block for given
! global indices ix, iy, iz
SUBROUTINE _i1DBlockNumber_(IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, ix, iy, iz, &
                             iBlockNumber, iz1)
USE dBase, ONLY : nxb, nyb, nzb
IMPLICIT NONE
INTEGER, INTENT(IN)  :: IDCurrCPU, iNumBlocksXDir, iNumBlocksYDir, ix, iy, iz
INTEGER, INTENT(OUT) :: iBlockNumber, iz1

INTEGER :: NPlanes, NLines, NCells

NPlanes = iz/(nzb + 1e-10) ; NLines  = iy/(nyb + 1e-10) ; NCells  = ix/(nxb + 1e-10)

iBlockNumber = NPlanes*iNumBlocksXDir*iNumBlocksYDir + &
               NLines*iNumBlocksXDir + NCells + 1
iz1 = iz - NPlanes*nzb ! Z-index inside the block (remember ghost cells!)

RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of Calc_Usr_Spectrum_Universal routine file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
