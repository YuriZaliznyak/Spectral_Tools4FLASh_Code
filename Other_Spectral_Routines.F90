SUBROUTINE Do_ALL_Spectral_Stuff

USE dBase, ONLY: dBasePropertyInteger, dBaseKeyNumber, GC, &
                 dBaseGetDataPtrSingleBlock, dBaseReleaseDataPtrSingleBlock, &
		 nxb, nyb, nzb, k2d, k3d, nguard, &
		 dBaseNodeType, dBaseGetCoords, dBasePropertyInteger
USE runtime_parameters, ONLY: get_parm_from_context, GLOBAL_PARM_CONTEXT
USE perfmon
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER :: Spectr_Timer
COMMON /SpectrTimerData/ Spectr_Timer !! shares Spectr_Timer with Calc_Usr_Spectrum_Universal...

INTEGER, PARAMETER :: q = max(nxb+2*nguard, nyb+2*nguard*k2d, nzb+2*nguard*k3d)
REAL, DIMENSION(q) :: XCoordArr, YCoordArr, ZCoordArr

INTEGER, SAVE :: iCurrSpectr, iDens, iPres, iVelx, iVely, iVelz, izn, iVar1, iVar2, &
                 iReSp, iImSp, iXCoord, iYCoord, iZCoord, iznl
LOGICAL, SAVE :: FirstCall = .TRUE.

INTEGER :: iSpectrSave, ix, iy, iz, iCurrBlockNum, iLocNumBlocks, IDMasterCPU, IDCurrCPU, &
           iStepNumber

REAL, DIMENSION (:,:,:,:), POINTER :: ALLBlockData
REAL :: CurrVx, CurrVy, CurrVz, CurrABSVel, CurrPress, rTemp1, rTemp2, &
        dVxDx, dVxDy, dVxDz, dVyDx, dVyDy, dVyDz, dVzDx, dVzDy, dVzDz, &
	DeltaX, DeltaY, DeltaZ, conductivity_constant, dPDx, dPDy, dPDz
INTEGER :: IERROR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Print *, 'Lazha-lazha-lazha!!!'


CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

IDMasterCPU = dBasePropertyInteger("MasterProcessor")
IDCurrCPU   = dBasePropertyInteger("MyProcessor")
iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")

IF (FirstCall) THEN
  FirstCall = .FALSE.
  iCurrSpectr = 1
  iVar1 = dBaseKeyNumber("var1")
  iDens = dBaseKeyNumber("dens")
  iPres = dBaseKeyNumber("pres")
  iVelx = dBaseKeyNumber("velx")
  iVely = dBaseKeyNumber("vely")
  iVelz = dBaseKeyNumber("velz")
  izn   = dBaseKeyNumber("zn")
  iznl  = dBaseKeyNumber("znl")
  iXcoord     =  dBaseKeyNumber("xCoord")
  iYcoord     =  dBaseKeyNumber("yCoord")
  iZcoord     =  dBaseKeyNumber("zCoord")

  iReSp = dBaseKeyNumber("resp")
  iImSp = dBaseKeyNumber("imsp")
EndIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writes solenoidal and compressible 1D spectra into files "SpectrumD..."
CALL Decompose_Spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates 1D condensed spectrum for the quantity
! $\rho \cdot \vec{v}^2$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Print *, 'Zdesja::', nguard, nzb+nguard*2, nyb+nguard*2, nxb+nguard*2

Do iCurrBlockNum = 1, iLocNumBlocks
  ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = 1, nzb+nguard*2
      Do iy = 1, nyb+nguard*2
        Do ix = 1, nxb+nguard*2
          CurrVx = ALLBlockData(iVelx, ix, iy, iz)
          CurrVy = ALLBlockData(iVely, ix, iy, iz)
          CurrVz = ALLBlockData(iVelz, ix, iy, iz)
          CurrABSVel = SQRT(CurrVx**2 + CurrVy**2 + CurrVz**2)

	  !ALLBlockData(iVar1,ix,iy,iz) = ALLBlockData(iDens,ix,iy,iz)*CurrABSVel**2

	  ALLBlockData(iVar1,ix,iy,iz) = CurrVx**2

        EndDo
      EndDo
    EndDo
  CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
EndDo


!iStepNumber = dBasePropertyInteger("CurrentStepNumber")
!IF (iStepNumber.EQ.25) THEN
!  CALL checkpoint_wr (1234, 9991.)
!  CALL Calc_Usr_Spectrum_Universal("var1")
!  CALL checkpoint_wr (1235, 9992.)
!EndIF

   CALL Calc_AND_Save_Spectrum_Universal("var1", "ENRG")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RETURN
END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine calculates and saves 1D condensed spectra for solenoidal and !!!
! compressible components of momentum density (\rho \cdot \vec{v}).    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Decompose_Spectrum

USE dBase, ONLY: dBasePropertyInteger, dBaseKeyNumber, GC, &
                 dBaseGetDataPtrSingleBlock, dBaseReleaseDataPtrSingleBlock, &
		 nxb, nyb, nzb, k2d, k3d, nguard, ndim, &
		 dBaseNodeType, dBaseGetCoords
USE runtime_parameters, ONLY: get_parm_from_context, set_parm_in_context, GLOBAL_PARM_CONTEXT

USE perfmon

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: Spectr_Timer
COMMON /SpectrTimerData/ Spectr_Timer !! shares Spectr_Timer with Calc_Usr_Spectrum_Universal...

LOGICAL, SAVE :: FirstCALL = .TRUE.
INTEGER, SAVE :: iCurrSpectr = 0, NumCPUs, iSpectrSave, iReSp, iImSp, iDens, &
                 iVelx, iVely, iVelz, iVar1, izn, iXcoord, iYcoord, iZcoord
INTEGER       :: IDCurrCPU, IDMasterCPU

INTEGER  :: IERROR, MPI_STATUS(MPI_STATUS_SIZE), MessTag

INTEGER, PARAMETER :: q = max(nxb+2*nguard, nyb+2*nguard*k2d, nzb+2*nguard*k3d)
REAL, DIMENSION(q) :: XCoordArr, YCoordArr, ZCoordArr
REAL, SAVE         :: XMin, XMax, YMin, YMax, ZMin, ZMax, SizeX, SizeY, SizeZ
REAL, DIMENSION (:,:,:,:), POINTER :: ALLBlockData

REAL, ALLOCATABLE :: ALLLocalSpectralData(:,:,:,:,:)

INTEGER, ALLOCATABLE  :: iTempArr1(:)
REAL, ALLOCATABLE     :: Local1DSpectrum(:,:), Global1DSpectrum(:,:), TempSpectrArr(:,:), &
                         r1DTempArr1(:), r1DTempArr2(:),r1DTempArr3(:), r2DTempArr1(:,:), &
			 Work(:), Condensed1DSpectrum(:,:)

REAL :: CurrVx, CurrVy, CurrVz, CurrABSVel, DeltaX, DeltaY, DeltaZ, OrtK(NDim), &
        AbsK, ReSpVal, ImSpVal, TempVar1, TempVar2, FactorKx, FactorKy, FactorKz, rTemp1, &
	ReVx_k, ReVy_k, ReVz_k, ImVx_k, ImVy_k, ImVz_k, ScalProdRe, ScalProdIm, &
	rPerpendicularSUM, rParallelSUM

INTEGER :: iCurrBlockNum, iLocNumBlocks, ix, iy, iz, ixStart, iyStart, izStart, &
           iNumBlocksXDir, iNumBlocksYDir, iNumBlocksZDir, &
	   Kx, Ky, Kz, NNx, NNy, NNz, iloop, iCount, iCountALL, iCPU, iTemp1, iTemp2, &
	   iPosOld, iPosNew, iPosition, lWork, iStepNumber, KIntCurr

CHARACTER*18 :: Filename
CHARACTER*5  :: StepNumString


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IDMasterCPU = dBasePropertyInteger("MasterProcessor")
IDCurrCPU   = dBasePropertyInteger("MyProcessor")
MessTag = 17
CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

IF (FirstCALL) THEN
  FirstCALL = .FALSE.
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, "iSpectrSave", iSpectrSave)
  NumCPUs = dBasePropertyInteger("NumberOfProcessors")
  iDens = dBaseKeyNumber("dens")
  iVelx = dBaseKeyNumber("velx")
  iVely = dBaseKeyNumber("vely")
  iVelz = dBaseKeyNumber("velz")
  iReSp = dBaseKeyNumber("resp")
  iImSp = dBaseKeyNumber("imsp")
  iVar1 = dBaseKeyNumber("var1")
  izn   = dBaseKeyNumber("zn")
  iXcoord     =  dBaseKeyNumber("xCoord")
  iYcoord     =  dBaseKeyNumber("yCoord")
  iZcoord     =  dBaseKeyNumber("zCoord")
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'xmin', XMin)
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'xmax', XMax)
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'ymin', YMin)
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'ymax', YMax)
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'zmin', ZMin)
  CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'zmax', ZMax)
  SizeX = ABS(XMax-XMin); SizeY = ABS(YMax-YMin); SizeZ = ABS(ZMax-ZMin)
  !CALL Calc_Usr_Spectrum_Universal ('iDens') ! needed to initialize "Spectr_Timer"
EndIF !FirstCALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL Timer_START (Spectr_Timer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!IF (iCurrSpectr.LT.iSpectrSave) THEN ! it's not a time to save the spectrum
!  iCurrSpectr = iCurrSpectr + 1
!  CALL Timer_STOP (Spectr_Timer)
!  RETURN
!EndIF

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
      NNx = iNumBlocksXDir*nxb
      NNy = iNumBlocksYDir*nyb
      NNz = iNumBlocksZDir*nzb
    EndIF
 EndDo
!!!###!!!
!############################################################################

iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")
! Storage space for spectral data is allocated at each CPU
ALLOCATE(ALLLocalSpectralData(iLocNumBlocks, nxb, nyb, nzb, 4))
!! Fourier transform is called three times, for the
!! quantities $\rho \cdot v_x$, $\rho \cdot v_y$, $\rho \cdot v_z$

! X-component
Do iCurrBlockNum = 1, iLocNumBlocks
  ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
!          ALLBlockData(iVar1, ix, iy, iz) = &
!	      ALLBlockData(iDens, ix, iy, iz)*(ALLBlockData(iVelx, ix, iy, iz))**2


          ALLBlockData(iVar1, ix, iy, iz) = &
	      (ALLBlockData(iVelx, ix, iy, iz))**2

        EndDo
      EndDo
    EndDo
  CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
EndDo
!! Make spectral transformation !!
CALL Timer_STOP (Spectr_Timer)
 CALL Calc_Usr_Spectrum_Universal ("var1")
CALL Timer_START (Spectr_Timer)

! Transform results are saved into local array "ALLLocalSpectralData" (positions 1 and 2)

Do iCurrBlockNum = 1, iLocNumBlocks
 IF (dBaseNodeType(iCurrBlockNum).EQ.1) THEN !! leaf node
   ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = 1, nzb
      Do iy = 1, nyb
        Do ix = 1, nxb
          ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 1) = &
	      ALLBlockData(iReSp, ix+nguard, iy+nguard, iz+nguard)
	  ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 2) = &
	      ALLBlockData(iImSp, ix+nguard, iy+nguard, iz+nguard)
        EndDo
      EndDo
    EndDo
   CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
 EndIF
EndDo

! Y-component
Do iCurrBlockNum = 1, iLocNumBlocks
  ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
!          ALLBlockData(iVar1, ix, iy, iz) = &
!	      ALLBlockData(iDens, ix, iy, iz)*(ALLBlockData(iVely, ix, iy, iz))**2

          ALLBlockData(iVar1, ix, iy, iz) = &
	      (ALLBlockData(iVely, ix, iy, iz))**2

        EndDo
      EndDo
    EndDo
  CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
EndDo
!! Make spectral transformation !!
CALL Timer_STOP (Spectr_Timer)
 CALL Calc_Usr_Spectrum_Universal ("var1")
CALL Timer_START (Spectr_Timer)

! Transform results are saved into local array "ALLLocalSpectralData" (positions 3 and 4)
Do iCurrBlockNum = 1, iLocNumBlocks
 IF (dBaseNodeType(iCurrBlockNum).EQ.1) THEN !! leaf node
   ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = 1, nzb
      Do iy = 1, nyb
        Do ix = 1, nxb
          ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 3) = &
	      ALLBlockData(iReSp, ix+nguard, iy+nguard, iz+nguard)
	  ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 4) = &
	      ALLBlockData(iImSp, ix+nguard, iy+nguard, iz+nguard)
        EndDo
      EndDo
    EndDo
   CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
 EndIF
EndDo

! Z-component:: results will remain in "ReSp" and "ImSp" variables...
Do iCurrBlockNum = 1, iLocNumBlocks
  ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
!          ALLBlockData(iVar1, ix, iy, iz) = &
!	      ALLBlockData(iDens, ix, iy, iz)*(ALLBlockData(iVelz, ix, iy, iz))**2

	   ALLBlockData(iVar1, ix, iy, iz) = &
	      (ALLBlockData(iVelz, ix, iy, iz))**2

        EndDo
      EndDo
    EndDo
  CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)
EndDo
!! Make spectral transformation !!
CALL Timer_STOP (Spectr_Timer)
 CALL Calc_Usr_Spectrum_Universal ("var1")
CALL Timer_START (Spectr_Timer)

iCount = 0 ! counter for the amount of 1D spectral records
ALLOCATE(Local1DSpectrum(NNx*NNy*NNz/4, 3)) !maximum possible size
Do iCurrBlockNum = 1, iLocNumBlocks
  IF (dBaseNodeType(iCurrBlockNum).EQ.1) THEN !! leaf node

    CALL dBaseGetCoords(izn, iXcoord, iCurrBlockNum, XCoordArr)
    CALL dBaseGetCoords(izn, iYcoord, iCurrBlockNum, YCoordArr)
    CALL dBaseGetCoords(izn, iZcoord, iCurrBlockNum, ZCoordArr)
    DeltaX = XCoordArr(2) - XCoordArr(1)
    DeltaY = YCoordArr(2) - YCoordArr(1)
    DeltaZ = ZCoordArr(2) - ZCoordArr(1)
    ixStart = (XCoordArr(nguard + 1) - XMin)/DeltaX + 1
    iyStart = (YCoordArr(nguard + 1) - YMin)/DeltaY + 1
    izStart = (ZCoordArr(nguard + 1) - ZMin)/DeltaZ + 1

    ALLBlockData => dBaseGetDataPtrSingleBlock(iCurrBlockNum, GC)
    Do iz = 1, nzb
      Do iy = 1, nyb
        Do ix = 1, nxb
	  Kx = ixStart + ix - 2
	  Ky = iyStart + iy - 2
	  Kz = izStart + iz - 2
	  AbsK = SQRT(FLOAT(Kx)**2+FLOAT(Ky)**2+FLOAT(Kz)**2)

	  IF (AbsK.GT.0..AND.Kx.LE.NNx/2.-1.AND.Ky.LE.NNy/2.-1.AND.Kz.LE.NNz/2.-1) THEN ! wavenumber is
	                                                          ! in the interesting region
	    iCount = iCount + 1

	! make projection into K-vector
            OrtK(1) = Kx/AbsK; OrtK(2) = Ky/AbsK; OrtK(3) = Kz/AbsK
            ReVx_k = ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 1)
	    ImVx_k = ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 2)
	    ReVy_k = ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 3)
	    ImVy_k = ALLLocalSpectralData(iCurrBlockNum, ix, iy, iz, 4)
	    ReVz_k = ALLBlockData(iReSp, ix+nguard,iy+nguard,iz+nguard)
	    ImVz_k = ALLBlockData(iImSp, ix+nguard,iy+nguard,iz+nguard)
	! Scalar product for real and imagimany parts
	    ScalProdRe = ReVx_k*OrtK(1) + ReVy_k*OrtK(2) + ReVz_k*OrtK(3)
	    ScalProdIm = ImVx_k*OrtK(1) + ImVy_k*OrtK(2) + ImVz_k*OrtK(3)

  ! If the perpendicular component is needed, then the parallel part should be substracted:
            ReVx_k = ReVx_k - OrtK(1)*ScalProdRe
            ReVy_k = ReVy_k - OrtK(2)*ScalProdRe
            ReVz_k = ReVz_k - OrtK(3)*ScalProdRe

            ImVx_k = ImVx_k - OrtK(1)*ScalProdIm
            ImVy_k = ImVy_k - OrtK(2)*ScalProdIm
            ImVz_k = ImVz_k - OrtK(3)*ScalProdIm

            ReSpVal = SQRT(ReVx_k**2 + ReVy_k**2 + ReVz_k**2)
            ImSpVal = SQRT(ImVx_k**2 + ImVy_k**2 + ImVz_k**2)

	    FactorKx = 1.; FactorKy = 1.; FactorKz = 1.
	    IF (Kx.EQ.0) FactorKx = 2.
	    IF (Ky.EQ.0) FactorKy = 2.
	    IF (Kz.EQ.0) FactorKz = 2.

	    TempVar1 = (8./(FactorKx*FactorKy*FactorKz))*(ScalProdRe**2 + ScalProdIm**2) !parallel
            TempVar2 = (8./(FactorKx*FactorKy*FactorKz))*(ReSpVal**2 + ImSpVal**2) !perpendicular

	    Local1DSpectrum(iCount, 1) = AbsK
	    Local1DSpectrum(iCount, 2) = TempVar1
	    Local1DSpectrum(iCount, 3) = TempVar2
	  EndIF
        EndDo ! iz
      EndDo ! iy
    EndDo ! ix
    CALL dBaseReleaseDataPtrSingleBlock(iCurrBlockNum, ALLBlockData)

  EndIF ! leaf node
EndDo ! do for all blocks

DEALLOCATE(ALLLocalSpectralData)

IF (iCount.EQ.0) DEALLOCATE(Local1DSpectrum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
IF (iCount.GT.0) THEN !"Local1DSpectrum" is managed if it is non-empty

! Remove non-used elements from "Local1DSpectrum"
 ALLOCATE(r2DTempArr1(iCount, 3))
 Do iloop = 1, iCount
  r2DTempArr1(iloop, :) = Local1DSpectrum(iloop, :)
 EndDo
 DEALLOCATE(Local1DSpectrum)
 ALLOCATE(Local1DSpectrum(iCount, 3))
 Local1DSpectrum(:,:) = r2DTempArr1(:,:)
 DEALLOCATE(r2DTempArr1)

! Sort Local1DSpectrum array accending
 ALLOCATE(iTempArr1(iCount), r1DTempArr1(iCount), r1DTempArr2(iCount), r1DTempArr3(iCount))
 r1DTempArr3(:) = Local1DSpectrum(:, 3)
 r1DTempArr1(:) = Local1DSpectrum(:, 2)
 r1DTempArr2(:) = Local1DSpectrum(:, 1)
 lWork = iCountALL
 ALLOCATE(Work(8*lWork))
 CALL DSORTS(r1DTempArr2, +1, iCount, iTempArr1, Work, 8*lWork)
 DEALLOCATE(Work)
 Do iloop = 1, iCount
     Local1DSpectrum(iloop, 1) = r1DTempArr2(iloop)
     Local1DSpectrum(iloop, 2) = r1DTempArr1(iTempArr1(iloop))
     Local1DSpectrum(iloop, 3) = r1DTempArr3(iTempArr1(iloop))
 EndDo
 DEALLOCATE(iTempArr1, r1DTempArr1, r1DTempArr2, r1DTempArr3)

! condense the spectrum by calculating the sum for repeating harmonics before communication
 ALLOCATE(r2DTempArr1(iCount, 3))
 iPosOld = 2
 iPosNew = 1
 r2DTempArr1(1, :) = Local1DSpectrum(1, :)
 Do While (iPosOld.LE.iCount)
  IF (r2DTempArr1(iPosNew, 1).EQ.Local1DSpectrum(iPosOld, 1)) THEN
    r2DTempArr1(iPosNew, 2) = r2DTempArr1(iPosNew, 2) + Local1DSpectrum(iPosOld, 2)
    r2DTempArr1(iPosNew, 3) = r2DTempArr1(iPosNew, 3) + Local1DSpectrum(iPosOld, 3)
                                                            ELSE
    iPosNew = iPosNew + 1
    r2DTempArr1(iPosNew, :) = Local1DSpectrum(iPosOld, :)
  EndIF
  iPosOld = iPosOld + 1
 EndDo ! while

 iCount = iPosNew
 DEALLOCATE(Local1DSpectrum)
 ALLOCATE(Local1DSpectrum(iCount, 3))
 Local1DSpectrum(1:iCount,:) = r2DTempArr1(1:iCount,:)
 DEALLOCATE(r2DTempArr1)

EndIF !iCount > 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! At this point, each process has spectral array "Local1DSpectrum" of the size 2*iCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! MasterCPU gathers information on the amount of iCount's from all processes
ALLOCATE(iTempArr1(NumCPUs))
CALL MPI_GATHER(iCount, 1, MPI_INTEGER, iTempArr1, 1, MPI_INTEGER, &
                  IDMasterCPU, MPI_COMM_WORLD, IERROR)
IF (IDCurrCPU.EQ.IDMasterCPU) THEN
 iCountALL = iTempArr1(1) !! iCountALL -- total amount of spectral points at all CPUs
 Do iloop = 2, NumCPUs
  iCountALL = iCountALL + iTempArr1(iloop)
 EndDo
 ALLOCATE(Global1DSpectrum(iCountALL, 3))
EndIF
CALL MPI_BCAST(iTempArr1, NumCPUs, MPI_INTEGER, IDMasterCPU, &
               MPI_COMM_WORLD, IERROR)

!!! First, MasterCPU manages its own part of spectrum
iPosition = 0
IF (IDCurrCPU.EQ.IDMasterCPU) THEN
 IF (iTempArr1(IDMasterCPU + 1).GT.0) THEN
   Do iloop = 1, iTempArr1(IDMasterCPU + 1)
     iPosition = iPosition + 1
     Global1DSpectrum(iPosition, :) = Local1DSpectrum(iloop, :)
   EndDo

 EndIF
   DEALLOCATE (Local1DSpectrum) ! At MasterCPU, the "Local1DSpectrum" array
                                ! was processed and can be deallocated
EndIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! COMMUNICATION: MasterCPU recieves spectral data from other CPUs !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Do iCPU = 0, NumCPUs - 1
 IF (IDCurrCPU.EQ.iCPU.OR.IDCurrCPU.EQ.IDMasterCPU) THEN ! iCPU communicates with MasteCPU
   IF (iTempArr1(iCPU + 1).GT.0.AND.iCPU.NE.IDMasterCPU) THEN ! there IS something to send

     IF (IDCurrCPU.EQ.IDMasterCPU) ALLOCATE(Local1DSpectrum(iTempArr1(iCPU + 1), 3))

       IF (IDCurrCPU.EQ.iCPU) &
          CALL MPI_SEND(Local1DSpectrum, 3*iTempArr1(iCPU + 1), MPI_DOUBLE_PRECISION, &
	              IDMasterCPU, MessTag, MPI_COMM_WORLD, IERROR)

       IF (IDCurrCPU.EQ.IDMasterCPU) THEN !123
          CALL MPI_RECV(Local1DSpectrum, 3*iTempArr1(iCPU + 1), MPI_DOUBLE_PRECISION, &
	              iCPU, MessTag, MPI_COMM_WORLD, MPI_STATUS, IERROR)
        ! MasterCPU stores the data in "Global1DSpectrum" array
          Do iloop = 1, iTempArr1(iCPU + 1)
	     iPosition = iPosition + 1
             Global1DSpectrum(iPosition, :) = Local1DSpectrum(iloop, :)
          EndDo
       EndIF !123

     IF (IDCurrCPU.EQ.IDMasterCPU) DEALLOCATE(Local1DSpectrum)

   EndIF ! there WAS something to send
 EndIF ! iCPU communicates with MasteCPU
EndDo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! COMMUNICATION ENDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DEALLOCATE(iTempArr1)

IF (IDCurrCPU.EQ.IDMasterCPU) THEN ! other stuff is done at MasterCPU only

! Sort "Global1DSpectrum" array accending
  ALLOCATE(iTempArr1(iCountALL), r1DTempArr1(iCountALL), r1DTempArr2(iCountALL), r1DTempArr3(iCountALL))
  r1DTempArr2(:) = Global1DSpectrum(:, 1)
  r1DTempArr1(:) = Global1DSpectrum(:, 2)
  r1DTempArr3(:) = Global1DSpectrum(:, 3)
  lWork = iCountALL
  ALLOCATE(Work(8*lWork))
  CALL DSORTS(r1DTempArr2, +1, iCountALL, iTempArr1, Work, 8*lWork)
  DEALLOCATE(Work)
  Do iloop = 1, iCountALL
      Global1DSpectrum(iloop, 1) = r1DTempArr2(iloop)
      Global1DSpectrum(iloop, 2) = r1DTempArr1(iTempArr1(iloop))
      Global1DSpectrum(iloop, 3) = r1DTempArr3(iTempArr1(iloop))
  EndDo
  DEALLOCATE(iTempArr1, r1DTempArr1, r1DTempArr2, r1DTempArr3)

! Condense "Global1DSpectrum" array
  ALLOCATE(r2DTempArr1(iCountALL, 3))
  iPosOld = 2
  iPosNew = 1
  r2DTempArr1(1, :) = Global1DSpectrum(1, :)
  Do While (iPosOld.LE.iCountALL)
   IF (r2DTempArr1(iPosNew, 1).EQ.Global1DSpectrum(iPosOld, 1)) THEN
     r2DTempArr1(iPosNew, 2) = r2DTempArr1(iPosNew, 2) + Global1DSpectrum(iPosOld, 2)
     r2DTempArr1(iPosNew, 3) = r2DTempArr1(iPosNew, 3) + Global1DSpectrum(iPosOld, 3)
                                                             ELSE
     iPosNew = iPosNew + 1
     r2DTempArr1(iPosNew, :) = Global1DSpectrum(iPosOld, :)
   EndIF
   iPosOld = iPosOld + 1
  EndDo ! while

  iCountALL = iPosNew
  DEALLOCATE(Global1DSpectrum)
  ALLOCATE(Global1DSpectrum(iCountALL, 3))
  Global1DSpectrum(1:iCountALL,:) = r2DTempArr1(1:iCountALL,:)
  DEALLOCATE(r2DTempArr1)

! Condense the 1D spectrum once more, leaving only integer wavenumbers
 ALLOCATE(Condensed1DSpectrum(iCountALL, 3), iTempArr1(iCountALL))

 Condensed1DSpectrum(:,:) = 0.0
 KIntCurr = NINT(Global1DSpectrum(1, 1))
 iTempArr1(1) = KIntCurr
 iPosNew = 1
 iCount = 0
 Do iloop = 1, iCountALL
   iTemp1 = NINT(Global1DSpectrum(iloop, 1))
   IF (iTemp1.NE.KIntCurr) THEN
        Condensed1DSpectrum(iPosNew,:) = &
	     Condensed1DSpectrum(iPosNew,:) / iCount
        iCount = 0
        KIntCurr = iTemp1
	iPosNew = iPosNew + 1
	iTempArr1(iPosNew) = KIntCurr
   EndIF
   Condensed1DSpectrum(iPosNew, :) = &
	     Condensed1DSpectrum(iPosNew, :) + Global1DSpectrum(iloop, :)
   iCount = iCount + 1
 EndDo

 DEALLOCATE(Global1DSpectrum)
 iCountALL = iPosNew

! MasterCPU writes spectrum to the file
  iStepNumber = dBasePropertyInteger("CurrentStepNumber")
  Write(StepNumString, 1019) 10000 + iStepNumber
1019	FORMAT(i5)
  Filename = 'SpectrumD'//StepNumString//'.dat'
  Open (1271, File=Filename)
  Write (1271, *) iCountALL
  Do iloop = 1, iCountALL
    Write(1271, *) iTempArr1(iloop), Condensed1DSpectrum(iloop, 2), Condensed1DSpectrum(iloop, 3)
  EndDo
  Close(1271)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  rParallelSUM = 0.0
!  rPerpendicularSUM = 0.0
!  Do iloop = 1, iCountALL
!    rParallelSUM = rParallelSUM + Condensed1DSpectrum(iloop, 2)
!    rPerpendicularSUM = rPerpendicularSUM + Condensed1DSpectrum(iloop, 3)
!  EndDo
!  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_3", rParallelSUM)
!  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_4", rPerpendicularSUM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DEALLOCATE(Condensed1DSpectrum, iTempArr1)

EndIF ! All other stuff was done by MasterCPU...

CALL Timer_STOP (Spectr_Timer)

CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

RETURN
END

