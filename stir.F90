!*******************************************************************************

!  Routine: stir()

!  Description: Wrapper for stirring routine.  This
!       routine handles the application of the stirer to all leaf-node
!       blocks.

SUBROUTINE stir
!===============================================================================
USE perfmon
USE runtime_parameters
USE dBase, ONLY: nxb, nyb, nzb, nguard,  dBaseGetCellVolume, &
     &                  k1d, k2d, k3d, ndim, maxblocks, GC, &
     &                  iHi_gc, jHi_gc, kHi_gc, &
     &                  dBaseKeyNumber, dBasePropertyReal, &
     &                  dBasePropertyInteger, dBaseBlockCoord, &
     &                  dBaseNodeType, dBaseGetDataPtrSingleBlock, &
                        dBaseReleaseDataPtrSingleBlock
IMPLICIT NONE
INCLUDE 'mpif.h'

REAL E_KIN
EXTERNAL E_KIN

LOGICAL, SAVE :: FirstCall = .TRUE.
INTEGER, SAVE :: ivelx, ively, ivelz, idens, iVar2, iVar3, iVar4, &
                 KStirMin, KStirMax, NModes

INTEGER :: IDMasterCPU, IDCurrCPU, IERROR
INTEGER :: iLocNumBlocks, iBlock, Date_Time(8), iTemp1, iTemp2, &
           iSign, iAlreadyHave, iCurrModeNumber, ix, iy, iz, iStepNumber
CHARACTER (LEN = 12) REAL_CLOCK (3)

REAL, SAVE :: st_stirmin, st_stirmax, st_energy, st_decay, RandSeed

REAL :: dT, CurrTime, rTemp1, rTemp2, rTempArr3(ndim), OrtK(ndim), AbsK, &
        ScalProd_RE, ScalProd_IM, rSum1, rSum2, &
	rLocalSum1, rLocalSum2, rGlobalSum1, rGlobalSum2, &
	rLocalBlockCounter, rGlobalBlockCounter, CurrKinEner, &
	EnergyInputPerCurrTimeStep, ForceFactor

REAL, ALLOCATABLE :: ModesArr(:,:), rMomentum_K_RE(:,:), rMomentum_K_IM(:,:), &
                     rMomentum_R(:,:,:,:)
REAL, DIMENSION (:,:,:,:), POINTER :: ALLBlockData
INTEGER :: ix1, iy1, iz1
!===============================================================================

!!! NO STIR!!!
!!RETURN


!!Print *, 'ENTERED STIR!!!'

 IDMasterCPU   = dBasePropertyInteger("MasterProcessor")
 IDCurrCPU     = dBasePropertyInteger("MyProcessor")
 CurrTime      = dBasePropertyReal("Time")
 dT            = dBasePropertyReal("TimeStep")
 iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")

IF(FirstCall) THEN
  FirstCall = .FALSE.
  iVelx = dBaseKeyNumber('velx')
  iVely = dBaseKeyNumber('vely')
  iVelz = dBaseKeyNumber('velz')
  iDens = dBaseKeyNumber('dens')

  iVar2 = dBaseKeyNumber('var2')
  iVar3 = dBaseKeyNumber('var3')
  iVar4 = dBaseKeyNumber('var4')

  CALL get_parm_from_context(global_parm_context, &
     &                                'st_stirmin', st_stirmin)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_stirmax', st_stirmax)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_energy', st_energy)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_decay', st_decay)

! Initialize random number generator for further "DURAND" calls with system time
  CALL DATE_AND_TIME(REAL_CLOCK(1), REAL_CLOCK(2), &
                     REAL_CLOCK(3), DATE_TIME)
  RandSeed = DBLE(Date_Time(8)+Date_Time(7)+Date_Time(6)+Date_Time(5))
EndIF !! FirstCall

KStirMin = st_stirmin
KStirMax = st_stirmax
NModes = st_decay
!!Print *, 'IN STIR::', KStirMin, KStirMax, NModes

ALLOCATE(ModesArr(ndim, 2*NModes))

!! At each call to "stir" routine, master processor defines set of perturbed
!! modes for the current stirring act. These modes are kept in the
!! "ModesArr(ndim, NModes)" array.
IF (IDCurrCPU.EQ.IDMasterCPU) THEN

! First mode will be always equal to (KStirMin, 0, 0)
  ModesArr(1,1) = KStirMin; ModesArr(2,1) = 0.; ModesArr(3,1) = 0.

  iCurrModeNumber = 2

  Do While (iCurrModeNumber.LE.NModes)
! Next element of modes array "ModesArr(*, iCurrModeNumber+1)" is initialized
! in a random way, with the components from the interval (KStirMin, KStirMax)
    Do iTemp1 = 1, ndim
      CALL DURAND(RandSeed, 1, rTemp1)
      !CALL RANDOM_NUMBER(rTemp1)
      rTempArr3(iTemp1) = FLOAT(KStirMin + NINT(rTemp1*(KStirMax - KStirMin)))
    EndDo

! Next, the mode is compared to other (already assigned) modes in order to
! avoid duplicate entries in the "ModesArr" array and aslo zero-modes (0,0,0)
    iAlreadyHave = 0
Do_iTemp1: Do iTemp1 = 1, iCurrModeNumber
            IF ((ModesArr(1,iTemp1).EQ.rTempArr3(1).AND.ModesArr(2,iTemp1).EQ.&
		 rTempArr3(2).AND.ModesArr(3,iTemp1).EQ.rTempArr3(3)).OR. &
		(rTempArr3(1)**2+rTempArr3(2)**2+rTempArr3(3)**2).EQ.0) THEN
                iAlreadyHave = 1
		EXIT Do_iTemp1
            EndIF
           EndDo Do_iTemp1
    IF (iAlreadyHave.EQ.0) THEN
     ModesArr(:, iCurrModeNumber) = rTempArr3(:)
     iCurrModeNumber = iCurrModeNumber + 1
    EndIF
  EndDo ! While

!! Modes array is then expanded by factor two, and the new part will
!! contain symmetric modes (k --> -k)...
  NModes = NModes*2
  Do iTemp1 = NModes/2+1, NModes
    ModesArr(:, iTemp1) = - ModesArr(:, iTemp1 - NModes/2)
  EndDo

EndIF !MasterCPU assigned random modes, array "ModesArr(ndim, NModes)"

!! MasterPE then broadcasts nmodes and "ModesArr" array to all other processors
CALL MPI_BCAST(NModes, 1, MPI_INTEGER, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(ModesArr, NModes*ndim, MPI_DOUBLE_PRECISION, &
	                  IDMasterCPU, MPI_COMM_WORLD, IERROR)

! Each CPU allocates space for spectral components of momentum perturbation
ALLOCATE(rMomentum_K_RE(ndim, NModes), rMomentum_K_IM(ndim, NModes))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The random values for the spectral components of
!! momentum $\rho \cdot \vec{v}$ are assigned for the current stirring act
!! and kept in the arrays "rMomentum_K_RE(IM)(ndim, NModes)"
IF (IDCurrCPU.EQ.IDMasterCPU) THEN
! Compute random forces
Do iTemp1 = 1, NModes/2 ! for the first half of "mode" array
 Do iTemp2 = 1, ndim ! for all dimensions, 1...3
  CALL DURAND(RandSeed, 1, rTemp1)
  CALL DURAND(RandSeed, 1, rTemp2)
  IF (rTemp2 < 0.5 ) THEN
     iSign=+1
                     ELSE
     iSign=-1
  EndIF
  rMomentum_K_RE(iTemp2, iTemp1) = iSign*rTemp1
 EndDo ! iTemp2
 Do iTemp2 = 1, ndim ! for all dimensions, 1...3
  CALL DURAND(RandSeed, 1, rTemp1)
  CALL DURAND(RandSeed, 1, rTemp2)
  IF (rTemp2 < 0.5 ) THEN
     iSign=+1
                     ELSE
     iSign=-1
  EndIf
  rMomentum_K_IM(iTemp2, iTemp1) = iSign*rTemp1
 EndDo ! iTemp2
EndDo !iTemp1. first half of "mode" array obtained random forces

!! The second part of "rMomentum_K_RE(IM)" array is filled as complex conjugate of
!! corresponding first part members
Do iTemp1 = NModes/2+1, NModes
  Do iTemp2 = 1, ndim
   rMomentum_K_RE(iTemp2, iTemp1) = + rMomentum_K_RE(iTemp2, iTemp1-NModes/2)
   rMomentum_K_IM(iTemp2, iTemp1) = - rMomentum_K_IM(iTemp2, iTemp1-NModes/2)
  EndDo
EndDo


!! Each spectral component of the momentum vector vector
!! (it's real and imaginary parts separately) is then
!! projected into the plane perpendicular to the
!! corresponding wavevector "ModesArr" by substracting the parallel to K component
Do iTemp1 = 1, NModes
! Calculates |K|
  AbsK = 0.0
  Do iTemp2 = 1, ndim
    AbsK = AbsK + ModesArr(iTemp2, iTemp1)**2
  EndDo
  AbsK = SQRT(AbsK)
! Calculates scalar products "rMomentum_K" by (K/|K|)
  ScalProd_RE = 0.0; ScalProd_IM = 0.0
  Do iTemp2 = 1, ndim
    OrtK(iTemp2) = ModesArr(iTemp2, iTemp1)/AbsK
    ScalProd_RE = ScalProd_RE + rMomentum_K_RE(iTemp2, iTemp1)*OrtK(iTemp2)
    ScalProd_IM = ScalProd_IM + rMomentum_K_IM(iTemp2, iTemp1)*OrtK(iTemp2)
  EndDo
! Updates "rMomentum_K_RE(IM)" by substracting parallel to K component
  Do iTemp2 = 1, ndim
    rMomentum_K_RE(iTemp2, iTemp1) = rMomentum_K_RE(iTemp2, iTemp1) &
                                     - OrtK(iTemp2)*ScalProd_RE
    rMomentum_K_IM(iTemp2, iTemp1) = rMomentum_K_IM(iTemp2, iTemp1) &
                                     - OrtK(iTemp2)*ScalProd_IM
  EndDo

EndDo !iTem1, modes


! Normalizes momenta to the total number of stirred modes
rMomentum_K_RE(:,:) = rMomentum_K_RE(:,:)/DBLE(NModes)
rMomentum_K_IM(:,:) = rMomentum_K_IM(:,:)/DBLE(NModes)


!Do iTemp1 = 1, NModes
!    Write(*,1121) iTemp1, ModesArr(:,iTemp1), rMomentum_K_RE(:,iTemp1), &
!                   rMomentum_K_IM(:,iTemp1)
!1121 FORMAT (I3, 3(1x,F5.2), 6(1x, E10.4))
!EndDo


EndIF ! MasterCPU has assigned random momenta "rMomentum_K_RE(IM)" in spectral space

!! MasterCPU broadcastss rMomentum_K_RE(IM) arrays to all other processors
CALL MPI_BCAST(rMomentum_K_RE, NModes*ndim, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	                  MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rMomentum_K_IM, NModes*ndim, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	                  MPI_COMM_WORLD, IERROR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Each CPU transforms spectral momentum arrays "rMomentum_K_RE(IM)" into
!! real array "rMomentum_R" and calculates sums rSum1, rSum2 for the
!! kinetic energy balance
ALLOCATE(rMomentum_R(ndim, nxb+2*nguard*k1d, nyb+2*nguard*k2d, nzb+2*nguard*k3d))
rLocalBlockCounter = 0.0; rLocalSum1 = 0.0; rLocalSum2 = 0.0
Do iBlock = 1, iLocNumBlocks
 IF (dBaseNodeType(iBlock).EQ.1) THEN !! the leaf node
  CALL calc_accel (iBlock, NModes, ModesArr, rMomentum_K_RE, rMomentum_K_IM, &
	                     rMomentum_R, rSum1, rSum2, 1)
  rLocalSum1 = rLocalSum1 + rSum1
  rLocalSum2 = rLocalSum2 + rSum2
  rLocalBlockCounter = rLocalBlockCounter + 1.

!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! The variable rMomentum_R components are saved into the global
!!! variables "var2", "var3", "var4".
    ALLBlockData => dBaseGetDataPtrSingleBlock(iBlock, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
          ALLBlockData(iVar2, ix, iy, iz) = rMomentum_R(1,ix,iy,iz)
          ALLBlockData(iVar3, ix, iy, iz) = rMomentum_R(2,ix,iy,iz)
          ALLBlockData(iVar4, ix, iy, iz) = rMomentum_R(3,ix,iy,iz)
        EndDo
      EndDo
    EndDo
    CALL dBaseReleaseDataPtrSingleBlock(iBlock, ALLBlockData)
!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 EndIF ! leaf node was processed
EndDo


CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
rGlobalSum1 = 0.0; rGlobalSum2 = 0.0
CALL MPI_REDUCE(rLocalSum1, rGlobalSum1, 1, MPI_DOUBLE_PRECISION, &
 	                MPI_SUM, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_REDUCE(rLocalSum2, rGlobalSum2, 1, MPI_DOUBLE_PRECISION, &
 	                MPI_SUM, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_REDUCE(rLocalBlockCounter, rGlobalBlockCounter, 1, MPI_DOUBLE_PRECISION, &
 	                MPI_SUM, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rGlobalSum1, 1, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	               MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rGlobalSum2, 1, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	               MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rGlobalBlockCounter, 1, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	               MPI_COMM_WORLD, IERROR)
IF (IDCurrCPU.EQ.IDMasterCPU) THEN

 !!! IF (rGlobalSum1.LT.1e-5) rGlobalSum1 = 1.

  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_10", rGlobalSum1)
  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_11", rGlobalSum2)
  Print *, 'rGlobalSum1 == ', rGlobalSum1 !!! , ' rGlobalSum2 == ', rGlobalSum2
EndIF

Do iBlock = 1, iLocNumBlocks
 IF (dBaseNodeType(iBlock).EQ.1) THEN !! the leaf node
  ALLBlockData => dBaseGetDataPtrSingleBlock(iBlock, GC)
  ALLBlockData(iVar2, :, :, :) = ALLBlockData(iVar2, :, :, :)/rGlobalSum1
  ALLBlockData(iVar3, :, :, :) = ALLBlockData(iVar3, :, :, :)/rGlobalSum1
  ALLBlockData(iVar4, :, :, :) = ALLBlockData(iVar4, :, :, :)/rGlobalSum1
  CALL dBaseReleaseDataPtrSingleBlock(iBlock, ALLBlockData)
 EndIF
EndDo
!!! For the gravity ideology, only rMomentum_R components are needed...
CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
GoTo 17895



!! Calculate kinetic energy BEFORE stirring
 CurrKinEner =  E_KIN("velx","vely","velz")
IF(IDCurrCPU.EQ.IDMasterCPU) Print *, 'KINETIC ENERGY ___BEFORE___ stirring::', &
                                       CurrKinEner


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Stirring:: velocities are updated
iStepNumber = dBasePropertyInteger("CurrentStepNumber")
IF (iStepNumber.GE.23.AND.iStepNumber.LE.27) THEN

EnergyInputPerCurrTimeStep = st_energy*dt

IF (rGlobalSum1.GE.0.) &
  ForceFactor = (-rGlobalSum1 + &
              SQRT(rGlobalSum1**2 + 4.*EnergyInputPerCurrTimeStep*rGlobalSum2))&
              /(2.*rGlobalSum2)
IF (rGlobalSum1.LT.0.) &
  ForceFactor = (-rGlobalSum1 - &
              SQRT(rGlobalSum1**2 + 4.*EnergyInputPerCurrTimeStep*rGlobalSum2))&
              /(2.*rGlobalSum2)

IF(IDCurrCPU.EQ.IDMasterCPU) Print *, 'ForceFactor==', ForceFactor

! Each processor updates velocities in it's set of blocks...

Do iBlock = 1, iLocNumBlocks
  IF (dBaseNodeType(iBlock).EQ.1) THEN !! the leaf node
    CALL calc_accel (iBlock, NModes, ModesArr, rMomentum_K_RE, rMomentum_K_IM, &
	                     rMomentum_R, rSum1, rSum2, 0)


!CALL Switch_To_Rotor(rMomentum_R)
!CALL Calc_Divergence(rMomentum_R, rTemp1, ix1, iy1, iz1)
!Print *, 'Max divergence for the block no.', iBlock, ' ::', rTemp1, ix1, iy1, iz1

    ALLBlockData => dBaseGetDataPtrSingleBlock(iBlock, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
	  rTemp1 = ForceFactor/ALLBlockData(iDens, ix, iy, iz)
          ALLBlockData(iVelx, ix, iy, iz) = &
		       ALLBlockData(iVelx, ix, iy, iz) + 1e-3!rMomentum_R(1,ix,iy,iz)*rTemp1
          ALLBlockData(iVely, ix, iy, iz) = &
	               ALLBlockData(iVely, ix, iy, iz) ! + rMomentum_R(2,ix,iy,iz)*rTemp1
	  ALLBlockData(iVelz, ix, iy, iz) = &
		       ALLBlockData(iVelz, ix, iy, iz) ! + rMomentum_R(3,ix,iy,iz)*rTemp1
        EndDo
      EndDo
    EndDo
  EndIF !leaf node
  CALL dBaseReleaseDataPtrSingleBlock(iBlock, ALLBlockData)
EndDo


!! Calculate kinetic energy AFTER stirring
rTemp1 = CurrKinEner
 CurrKinEner =  E_KIN("velx","vely","velz")
IF(IDCurrCPU.EQ.IDMasterCPU) Print *, 'KINETIC ENERGY ___AFTER___ stirring::', &
            CurrKinEner, ' dE/dT=', (CurrKinEner - rTemp1)/dT

EndIF ! iStepNumber




!!! If the only rMomentum_R components are needed (gravity ideology), jump here...
17895     CONTINUE

DEALLOCATE(rMomentum_R)
DEALLOCATE(ModesArr, rMomentum_K_RE, rMomentum_K_IM)




!Do iTemp1 = 1, NModes
!    Write(*,1121) iTemp1, ModesArr(:,iTemp1), rMomentum_K_RE(:,iTemp1), &
!                   rMomentum_K_IM(:,iTemp1)
!1121 FORMAT (I3, 3(1x,F5.2), 6(1x, E10.4))
!EndDo


RETURN
END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Build_Gravity
!===============================================================================
USE perfmon
USE runtime_parameters
USE dBase, ONLY: nxb, nyb, nzb, nguard,  dBaseGetCellVolume, &
     &                  k1d, k2d, k3d, ndim, maxblocks, GC, &
     &                  iHi_gc, jHi_gc, kHi_gc, &
     &                  dBaseKeyNumber, dBasePropertyReal, &
     &                  dBasePropertyInteger, dBaseBlockCoord, &
     &                  dBaseNodeType, dBaseGetDataPtrSingleBlock, &
                        dBaseReleaseDataPtrSingleBlock
IMPLICIT NONE
INCLUDE 'mpif.h'

REAL E_KIN
EXTERNAL E_KIN

LOGICAL, SAVE :: FirstCall = .TRUE.
INTEGER, SAVE :: ivelx, ively, ivelz, idens, iVar2, iVar3, iVar4, &
                 KStirMin, KStirMax, NModes

INTEGER :: IDMasterCPU, IDCurrCPU, IERROR
INTEGER :: iLocNumBlocks, iBlock, Date_Time(8), iTemp1, iTemp2, &
           iSign, iAlreadyHave, iCurrModeNumber, ix, iy, iz, iStepNumber
CHARACTER (LEN = 12) REAL_CLOCK (3)

REAL, SAVE :: st_stirmin, st_stirmax, st_energy, st_decay, RandSeed

REAL :: dT, CurrTime, rTemp1, rTemp2, rTempArr3(ndim), OrtK(ndim), AbsK, &
        ScalProd_RE, ScalProd_IM, rSum1, rSum2, &
	rLocalSum1, rLocalSum2, rGlobalSum1, rGlobalSum2, &
	rLocalBlockCounter, rGlobalBlockCounter, CurrKinEner, &
	EnergyInputPerCurrTimeStep, ForceFactor

REAL, ALLOCATABLE :: ModesArr(:,:), rMomentum_K_RE(:,:), rMomentum_K_IM(:,:), &
                     rMomentum_R(:,:,:,:)
REAL, DIMENSION (:,:,:,:), POINTER :: ALLBlockData
INTEGER :: ix1, iy1, iz1
!===============================================================================

 IDMasterCPU   = dBasePropertyInteger("MasterProcessor")
 IDCurrCPU     = dBasePropertyInteger("MyProcessor")
 CurrTime      = dBasePropertyReal("Time")
 dT            = dBasePropertyReal("TimeStep")
 iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")

IF(FirstCall) THEN
  FirstCall = .FALSE.
  iVelx = dBaseKeyNumber('velx')
  iVely = dBaseKeyNumber('vely')
  iVelz = dBaseKeyNumber('velz')
  iDens = dBaseKeyNumber('dens')

  iVar2 = dBaseKeyNumber('var2')
  iVar3 = dBaseKeyNumber('var3')
  iVar4 = dBaseKeyNumber('var4')

  CALL get_parm_from_context(global_parm_context, &
     &                                'st_stirmin', st_stirmin)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_stirmax', st_stirmax)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_energy', st_energy)
  CALL get_parm_from_context(global_parm_context, &
     &                                'st_decay', st_decay)

! Initialize random number generator for further "DURAND" calls with system time
  CALL DATE_AND_TIME(REAL_CLOCK(1), REAL_CLOCK(2), &
                     REAL_CLOCK(3), DATE_TIME)
  RandSeed = DBLE(Date_Time(8)+Date_Time(7)+Date_Time(6)+Date_Time(5))

EndIF !! FirstCall

KStirMin = st_stirmin
KStirMax = st_stirmax
NModes = st_decay
!!Print *, 'IN STIR::', KStirMin, KStirMax, NModes

ALLOCATE(ModesArr(ndim, 2*NModes))

!! At each call to "stir" routine, master processor defines set of perturbed
!! modes for the current stirring act. These modes are kept in the
!! "ModesArr(ndim, NModes)" array.
IF (IDCurrCPU.EQ.IDMasterCPU) THEN

! First mode will be always equal to (KStirMin, 0, 0)
  ModesArr(1,1) = KStirMin; ModesArr(2,1) = 0.; ModesArr(3,1) = 0.

  iCurrModeNumber = 2

  Do While (iCurrModeNumber.LE.NModes)
! Next element of modes array "ModesArr(*, iCurrModeNumber+1)" is initialized
! in a random way, with the components from the interval (KStirMin, KStirMax)
    Do iTemp1 = 1, ndim
     ! CALL DURAND(RandSeed, 1, rTemp1)

   CALL RANDOM_NUMBER(rTemp1)

      rTempArr3(iTemp1) = FLOAT(KStirMin + NINT(rTemp1*(KStirMax - KStirMin)))
    EndDo

! Next, the mode is compared to other (already assigned) modes in order to
! avoid duplicate entries in the "ModesArr" array and aslo zero-modes (0,0,0)
    iAlreadyHave = 0
Do_iTemp1: Do iTemp1 = 1, iCurrModeNumber
            IF ((ModesArr(1,iTemp1).EQ.rTempArr3(1).AND.ModesArr(2,iTemp1).EQ.&
		 rTempArr3(2).AND.ModesArr(3,iTemp1).EQ.rTempArr3(3)).OR. &
		(rTempArr3(1)**2+rTempArr3(2)**2+rTempArr3(3)**2).EQ.0) THEN
                iAlreadyHave = 1
		EXIT Do_iTemp1
            EndIF
           EndDo Do_iTemp1
    IF (iAlreadyHave.EQ.0) THEN
     ModesArr(:, iCurrModeNumber) = rTempArr3(:)
     iCurrModeNumber = iCurrModeNumber + 1
    EndIF
  EndDo ! While

!! Modes array is then expanded by factor two, and the new part will
!! contain symmetric modes (k --> -k)...
  NModes = NModes*2
  Do iTemp1 = NModes/2+1, NModes
    ModesArr(:, iTemp1) = - ModesArr(:, iTemp1 - NModes/2)
  EndDo

!Do iTemp1 = 1, NModes/2
!  Write(*,10201) ModesArr(:, iTemp1)
!10201 FORMAT(100(1x, F4.1))
!EndDo

EndIF !MasterCPU assigned random modes, array "ModesArr(ndim, NModes)"

!! MasterPE then broadcasts nmodes and "ModesArr" array to all other processors
CALL MPI_BCAST(NModes, 1, MPI_INTEGER, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(ModesArr, NModes*ndim, MPI_DOUBLE_PRECISION, &
	                  IDMasterCPU, MPI_COMM_WORLD, IERROR)

! Each CPU allocates space for spectral components of momentum perturbation
ALLOCATE(rMomentum_K_RE(ndim, NModes), rMomentum_K_IM(ndim, NModes))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The random values for the spectral components of
!! momentum $\rho \cdot \vec{v}$ are assigned for the current stirring act
!! and kept in the arrays "rMomentum_K_RE(IM)(ndim, NModes)"
IF (IDCurrCPU.EQ.IDMasterCPU) THEN
! Compute random forces
Do iTemp1 = 1, NModes/2 ! for the first half of "mode" array
 Do iTemp2 = 1, ndim ! for all dimensions, 1...3
  !CALL DURAND(RandSeed, 1, rTemp1)
  !CALL DURAND(RandSeed, 1, rTemp2)

  CALL RANDOM_NUMBER(rTemp1)
  CALL RANDOM_NUMBER(rTemp2)

  IF (rTemp2 < 0.5 ) THEN
     iSign=+1
                     ELSE
     iSign=-1
  EndIF
  rMomentum_K_RE(iTemp2, iTemp1) = iSign*rTemp1
 EndDo ! iTemp2
 Do iTemp2 = 1, ndim ! for all dimensions, 1...3
  !CALL DURAND(RandSeed, 1, rTemp1)
  !CALL DURAND(RandSeed, 1, rTemp2)

  CALL RANDOM_NUMBER(rTemp1)
  CALL RANDOM_NUMBER(rTemp2)

  IF (rTemp2 < 0.5 ) THEN
     iSign=+1
                     ELSE
     iSign=-1
  EndIf
  rMomentum_K_IM(iTemp2, iTemp1) = iSign*rTemp1
 EndDo ! iTemp2
EndDo !iTemp1. first half of "mode" array obtained random forces

!! The second part of "rMomentum_K_RE(IM)" array is filled as complex conjugate of
!! corresponding first part members
Do iTemp1 = NModes/2+1, NModes
  Do iTemp2 = 1, ndim
   rMomentum_K_RE(iTemp2, iTemp1) = + rMomentum_K_RE(iTemp2, iTemp1-NModes/2)
   rMomentum_K_IM(iTemp2, iTemp1) = - rMomentum_K_IM(iTemp2, iTemp1-NModes/2)
  EndDo
EndDo


!! Each spectral component of the momentum vector vector
!! (it's real and imaginary parts separately) is then
!! projected into the plane perpendicular to the
!! corresponding wavevector "ModesArr" by substracting the parallel to K component
Do iTemp1 = 1, NModes
! Calculates |K|
  AbsK = 0.0
  Do iTemp2 = 1, ndim
    AbsK = AbsK + ModesArr(iTemp2, iTemp1)**2
  EndDo
  AbsK = SQRT(AbsK)
! Calculates scalar products "rMomentum_K" by (K/|K|)
  ScalProd_RE = 0.0; ScalProd_IM = 0.0
  Do iTemp2 = 1, ndim
    OrtK(iTemp2) = ModesArr(iTemp2, iTemp1)/AbsK
    ScalProd_RE = ScalProd_RE + rMomentum_K_RE(iTemp2, iTemp1)*OrtK(iTemp2)
    ScalProd_IM = ScalProd_IM + rMomentum_K_IM(iTemp2, iTemp1)*OrtK(iTemp2)
  EndDo
! Updates "rMomentum_K_RE(IM)" by substracting parallel to K component
  Do iTemp2 = 1, ndim
    rMomentum_K_RE(iTemp2, iTemp1) = rMomentum_K_RE(iTemp2, iTemp1) &
                                     - OrtK(iTemp2)*ScalProd_RE
    rMomentum_K_IM(iTemp2, iTemp1) = rMomentum_K_IM(iTemp2, iTemp1) &
                                     - OrtK(iTemp2)*ScalProd_IM
  EndDo

EndDo !iTem1, modes


! Normalizes momenta to the total number of stirred modes
rMomentum_K_RE(:,:) = rMomentum_K_RE(:,:)/DBLE(NModes)
rMomentum_K_IM(:,:) = rMomentum_K_IM(:,:)/DBLE(NModes)

EndIF ! MasterCPU has assigned random momenta "rMomentum_K_RE(IM)" in spectral space

!! MasterCPU broadcastss rMomentum_K_RE(IM) arrays to all other processors
CALL MPI_BCAST(rMomentum_K_RE, NModes*ndim, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	                  MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rMomentum_K_IM, NModes*ndim, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	                  MPI_COMM_WORLD, IERROR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Each CPU transforms spectral momentum arrays "rMomentum_K_RE(IM)" into
!! real array "rMomentum_R" and calculates sums rSum1, rSum2 for the
!! kinetic energy balance
ALLOCATE(rMomentum_R(ndim, nxb+2*nguard*k1d, nyb+2*nguard*k2d, nzb+2*nguard*k3d))
rLocalBlockCounter = 0.0; rLocalSum1 = 0.0; rLocalSum2 = 0.0
Do iBlock = 1, iLocNumBlocks
 IF (dBaseNodeType(iBlock).EQ.1) THEN !! the leaf node
  CALL calc_accel (iBlock, NModes, ModesArr, rMomentum_K_RE, rMomentum_K_IM, &
	                     rMomentum_R, rSum1, rSum2, 1)


!Print*, 'GForce::', MAXVAL(rMomentum_R(:,:,:,:)), MINVAL(rMomentum_R(:,:,:,:))


  rLocalSum1 = rLocalSum1 + rSum1
  rLocalSum2 = rLocalSum2 + rSum2
  rLocalBlockCounter = rLocalBlockCounter + 1.

!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! The variable rMomentum_R components are saved into the global
!!! variables "var2", "var3", "var4".
    ALLBlockData => dBaseGetDataPtrSingleBlock(iBlock, GC)
    Do iz = nguard*k3d+1, nguard*k3d+nzb
      Do iy = nguard*k2d+1, nguard*k2d+nyb
        Do ix = nguard+1, nguard+nxb
          ALLBlockData(iVar2, ix, iy, iz) = rMomentum_R(1,ix,iy,iz)
          ALLBlockData(iVar3, ix, iy, iz) = rMomentum_R(2,ix,iy,iz)
          ALLBlockData(iVar4, ix, iy, iz) = rMomentum_R(3,ix,iy,iz)
        EndDo
      EndDo
    EndDo

!Print*, 'VelX::', MAXVAL(ALLBlockData(iVelx, :, :, :)), MINVAL(ALLBlockData(iVelx, :, :, :))
!Print*, 'VelY::', MAXVAL(ALLBlockData(iVely, :, :, :)), MINVAL(ALLBlockData(iVely, :, :, :))
!Print*, 'VelZ::', MAXVAL(ALLBlockData(iVelz, :, :, :)), MINVAL(ALLBlockData(iVelz, :, :, :))

    CALL dBaseReleaseDataPtrSingleBlock(iBlock, ALLBlockData)
!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 EndIF ! leaf node was processed
EndDo


!Print*, MAXVAL(rMomentum_R(:,:,:,:))


DEALLOCATE(rMomentum_R)
DEALLOCATE(ModesArr, rMomentum_K_RE, rMomentum_K_IM)

CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

RETURN
END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Find_Force_Factor
USE runtime_parameters
USE dBase, ONLY: nxb, nyb, nzb, nguard,  k1d, k2d, k3d, ndim, maxblocks, GC, &
     &                  dBaseKeyNumber, dBasePropertyReal, &
     &                  dBasePropertyInteger, dBaseGetCellVolume, &
     &                  dBaseNodeType, dBaseGetDataPtrSingleBlock, &
                        dBaseReleaseDataPtrSingleBlock
IMPLICIT NONE
INCLUDE 'mpif.h'

LOGICAL, SAVE :: FirstCall = .TRUE.
INTEGER, SAVE :: ivelx, ively, ivelz, idens, iVar2, iVar3, iVar4

INTEGER :: IDMasterCPU, IDCurrCPU, IERROR
INTEGER :: iLocNumBlocks, iStepNumber, iBlock, ix, iy, iz

REAL, SAVE :: st_energy
REAL :: DensVal, VelxVal, VelyVal, VelzVal, &
        GravityX_Val, GravityY_Val, GravityZ_Val, &
	rSum1Local, rSum1Global, dVol
REAL, DIMENSION(:,:,:,:), POINTER :: AllBlockData
LOGICAL :: FlagISNAN, ISNAN
EXTERNAL ISNAN
!===============================================================================

 IDMasterCPU   = dBasePropertyInteger("MasterProcessor")
 IDCurrCPU     = dBasePropertyInteger("MyProcessor")
 iLocNumBlocks = dBasePropertyInteger("LocalNumberOfBlocks")
 iStepNumber   = dBasePropertyInteger("CurrentStepNumber")
IF(FirstCall) THEN
  FirstCall = .FALSE.
  iVelx = dBaseKeyNumber('velx')
  iVely = dBaseKeyNumber('vely')
  iVelz = dBaseKeyNumber('velz')
  iDens = dBaseKeyNumber('dens')

  iVar2 = dBaseKeyNumber('var2')
  iVar3 = dBaseKeyNumber('var3')
  iVar4 = dBaseKeyNumber('var4')

  CALL get_parm_from_context(global_parm_context, &
     &                                'st_energy', st_energy)
EndIF !! FirstCall

rSum1Local = 0.0; rSum1Global = 0.0
Do iBlock = 1, iLocNumBlocks
 IF (dBaseNodeType(iBlock).EQ.1) THEN !! the leaf node
  AllBlockData => dBaseGetDataPtrSingleBlock(iBlock, GC)

  CALL dBaseGetCellVolume(nxb/2, nyb/2, nzb/2, iBlock, dVol)
  Do iz = k3d*nguard+1, k3d*nguard+nzb
   Do iy = k2d*nguard+1, k2d*nguard+nyb
    Do ix = nguard+1, nguard+nxb
     VelxVal = AllBlockData(iVelx,ix,iy,iz)
     VelyVal = AllBlockData(iVely,ix,iy,iz)
     VelzVal = AllBlockData(iVelz,ix,iy,iz)
     DensVal = AllBlockData(iDens,ix,iy,iz)
     GravityX_Val = AllBlockData(iVar2,ix,iy,iz)
     GravityY_Val = AllBlockData(iVar3,ix,iy,iz)
     GravityZ_Val = AllBlockData(iVar4,ix,iy,iz)

     rSum1Local = rSum1Local + DensVal*( ( VelxVal*GravityX_Val + &
            VelyVal*GravityY_Val + VelzVal*GravityZ_Val) )*dVol

    EndDo !x, ix-direction
   EndDo !y, iy-direction
  EndDo !z, iz-direction
  CALL dBaseReleaseDataPtrSingleBlock(iBlock, AllBlockData)
 EndIF ! leaf node
EndDo

CALL MPI_REDUCE(rSum1Local, rSum1Global, 1, MPI_DOUBLE_PRECISION, &
 	                MPI_SUM, IDMasterCPU, MPI_COMM_WORLD, IERROR)
CALL MPI_BCAST(rSum1Global, 1, MPI_DOUBLE_PRECISION, IDMasterCPU, &
	               MPI_COMM_WORLD, IERROR)

IF (IDCurrCPU.EQ.IDMasterCPU) Print *, 'ForceFactor::', rSum1Global

FlagISNAN = ISNAN(rSum1Global)
IF (FlagISNAN) THEN
  Print *, "CPU::", IDCurrCPU, "  Force factor is invalid!"
!  CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  CALL abort_flash()
EndIF


CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_10", rSum1Global)

RETURN
END
