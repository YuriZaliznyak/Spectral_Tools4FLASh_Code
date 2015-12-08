!*******************************************************************************

!  Routine:     evolve()

!  Description: Apply the operators needed to advance a simulation
!               through one timestep.

!===============================================================================

subroutine evolve ()

!  interface: null
!  local variables: first_call
!  external modules: common

  use hydro, ONLY: hydro_3d

  use Gravity, ONLY: GravPotentialAllBlocks

  use runtime_parameters, ONLY: get_parm_from_context, set_parm_in_context, &
                                GLOBAL_PARM_CONTEXT
  use source_termsModule, ONLY: source_terms

  use dBase, ONLY: dBasePropertyReal,    &
                   dBaseSetProperty,     &
                   sweep_order_xyz,      &
                   sweep_order_zyx,      &
                   dbaseNodeType,        &
                   dbasePropertyInteger, &
                   maxblocks, ndim

  use ParticleModule, ONLY: AdvanceParticles, &
                            RunningParticles, ReDistributeParticles

  use RadiationModule, ONLY:  AdvanceRadiationField

  use Cosmology, ONLY: SolveFriedmannEquation, RedshiftHydroQuantities

  use perfmon

  implicit none

  INCLUDE 'mpif.h'

  logical,save  :: first_call = .true.

  integer, save :: itemp_limit
  integer, save :: igrav

  real          :: time, dt

  integer       :: block_count, lnblocks, n
  integer       :: block_list(maxblocks)


  INTEGER, SAVE :: iGravityCounter, iLocalTimestepCounter, KMin, KMax
  INTEGER :: iStepNumber, iSpectrSave, IERROR, IDMasterCPU, IDCurrCPU, iTemp1
  REAL, SAVE :: PreviousForceFactor, CurrentForceFactor
  LOGICAL    :: RandomGravityOK
  INTEGER, PARAMETER    :: NModes = 30, iGravityON = 1, iProjectFlag = 1, &
                           iMagneticFieldClean = 0
  REAL, SAVE :: ModesArr(ndim, NModes), rMomentum_K_RE(ndim, NModes), &
                rMomentum_K_IM(ndim, NModes), PhasesArr(ndim), &
		ModesArr2(ndim, NModes), rMomentum_K_RE2(ndim, NModes), &
                rMomentum_K_IM2(ndim, NModes), rLocalTimePrevious, rLocalTimeCurr, &
		st_energy, ZeroLevelOfInternalEnergy
  REAL, PARAMETER :: PICONST = 3.1415927410125732421875D0, ForceDecayTime = 5.0
  REAL       :: rTemp1, rTemp2, rTemp3, rTemp4, CurrKinEner, TTT, rTemp22

  REAL VAR_INTEGRATE
  EXTERNAL VAR_INTEGRATE
  REAL :: PrevEInt, CurrEInt, PrevETot, CurrETot, CurrEKin, rTemp001
  LOGICAL :: RestartFlag

  LOGICAL, PARAMETER ::  InternalEnergyCorrectorFlag = .FALSE., &
                         ForcingUpdateFlag = .FALSE.

  REAL, SAVE  :: Forcing_Energy_CURR,Forcing_Energy_MIN, Forcing_Energy_MAX, &
                     MachMIN, MachMAX, MachUNDERMIN, MachOVERMAX, MachMiddle, &
		     dEnergy, dEnergyBasic, PreviousMachNumber, &
		     rFactor1, rFactor2

!-------------------------------------------------------------------------------
  IDMasterCPU   = dBasePropertyInteger("MasterProcessor")
  IDCurrCPU     = dBasePropertyInteger("MyProcessor")
  lnblocks      = dBasePropertyInteger("LocalNumberOfBlocks")

!               On the first time through the loop, get any necessary runtime
!               parameters before proceeding.
  IF (first_call) THEN
     first_call = .FALSE.

     !!CALL Rescale_Magnetic_Field(2.0)

     PhasesArr(:) = 0.
     iLocalTimestepCounter = 0
     IF (IDCurrCPU.EQ.IDMasterCPU) &
	       Print *, "InternalEnergyCorrectorFlag::", InternalEnergyCorrectorFlag, &
              "  ForcingUpdateFlag::", ForcingUpdateFlag, &
	      "  iMagneticFieldClean::", iMagneticFieldClean
     CALL get_parm_from_context(global_parm_context, "st_stirmin", rTemp1); KMin = rTemp1
     CALL get_parm_from_context(global_parm_context, "st_stirmax", rTemp1); KMax = rTemp1
     CALL get_parm_from_context(global_parm_context, "Forcing_Energy_CURR", Forcing_Energy_CURR)
     CALL get_parm_from_context(global_parm_context, "Forcing_Energy_MIN", Forcing_Energy_MIN)
     CALL get_parm_from_context(global_parm_context, "Forcing_Energy_MAX", Forcing_Energy_MAX)
     CALL get_parm_from_context(global_parm_context, "MachMIN", MachMIN)
     CALL get_parm_from_context(global_parm_context, "MachMAX", MachMAX)
     CALL get_parm_from_context(global_parm_context, "MachMIN", MachMIN)
     CALL get_parm_from_context(global_parm_context, "MachUNDERMIN", MachUNDERMIN)
     CALL get_parm_from_context(global_parm_context, "MachOVERMAX", MachOVERMAX)
     CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_13", Forcing_Energy_CURR)
     dEnergyBasic = (Forcing_Energy_MAX - Forcing_Energy_MIN)/25.
     dEnergy = dEnergyBasic
     MachMiddle = (MachMIN + MachMAX)/2.
     call get_parm_from_context(global_parm_context, "itemp_limit", itemp_limit)
     call get_parm_from_context(global_parm_context, "igrav", igrav)

     dt   = dBasePropertyReal("TimeStep")
     call dBaseSetProperty("OldTimeStep", dt)

     iGravityCounter = 0

     CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, 'Restart_HDF_File_Num', iTemp1)
     IF (iTemp1.NE.-1) THEN
        IF (IDCurrCPU.EQ.IDMasterCPU) &
	          Print *, 'Restart_HDF_File_Num ==', iTemp1, &
		           ' --> HDF regridder was applied...'
	!!CALL Load_Modes_And_Forces(NModes, ModesArr, rMomentum_K_RE, rMomentum_K_IM)
     EndIF

     CALL get_parm_from_context(global_parm_context, "restart", RestartFlag)
     IF (RestartFlag) THEN
        IF (IDCurrCPU.EQ.IDMasterCPU) &
	          Print *, 'RestartFlag == .TRUE., simulation restarted...'
		      ELSE
      	IF (IDCurrCPU.EQ.IDMasterCPU) &
	          Print *, 'RestartFlag == .FALSE., simulation started from scratch...'
     EndIF

        CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr)
	CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr2)
        CALL Build_Spectral_Forces_RE(NModes, ModesArr, rMomentum_K_RE, iProjectFlag)
        CALL Build_Spectral_Forces_IM(NModes, ModesArr, rMomentum_K_IM, iProjectFlag)
	CALL Build_Spectral_Forces_RE(NModes, ModesArr2, rMomentum_K_RE2, iProjectFlag)
        CALL Build_Spectral_Forces_IM(NModes, ModesArr2, rMomentum_K_IM2, iProjectFlag)
        CALL Build_Real_Forces(NModes, ModesArr, &
                      rMomentum_K_RE, rMomentum_K_IM, PhasesArr, 1., 1)
        CALL Build_Real_Forces(NModes, ModesArr2, &
                      rMomentum_K_RE2, rMomentum_K_IM2, PhasesArr, 1., 0)
        CALL Find_Force_Factor

     CALL Calc_Usr_Integrals
     ! needed to initialize "Spectr_Timer"
     CALL Calc_Usr_Spectrum_Universal ("velx", "resp", "imsp")

     rFactor1 = 1.0; rFactor2 = 0.0
     rLocalTimePrevious = 0.0; rLocalTimeCurr = 0.0

     CALL get_parm_from_context(global_parm_context, &
     &                                "Forcing_Energy_CURR", st_energy)
     ZeroLevelOfInternalEnergy = VAR_INTEGRATE("eint")


  endif !! FirstCALL

iLocalTimestepCounter = iLocalTimestepCounter + 1
rTemp4 = st_energy*TANH(0.15*iLocalTimestepCounter)
CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Forcing_Energy_CURR", rTemp4)
CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_13", rTemp4)

!IF (IDCurrCPU.EQ.IDMasterCPU) Print *, 'Current Forcing Energy::', rTemp4


! Calculates user-supplied integrals to be written into "flash.dat"
CALL Calc_Usr_Integrals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  block_count = 0
  do n = 1, lnblocks
     if (dBaseNodeType(n) == 1) then      ! keep only leaf_blocks
        block_count = block_count + 1
        block_list(block_count) = n
     endif
  enddo
  cu_block_count = cu_block_count+block_count

!-------------------------------------------------------------------------------

!               Use second-order operator splitting (Strang 1968).
!               First set of operator calls...

! Time update

  time = dBasePropertyReal("Time")
  dt   = dBasePropertyReal("TimeStep")

  rLocalTimeCurr = rLocalTimeCurr + dt
  iGravityCounter = iGravityCounter + 1
  IF (iGravityCounter.EQ.iGravityON) THEN
    IF (rLocalTimeCurr.GE.ForceDecayTime.AND.rLocalTimePrevious.LT.ForceDecayTime) THEN
!      IF (IDCurrCPU.EQ.IDMasterCPU) &
!         Print *, 'Time::', time, rLocalTimeCurr, ' ... Modes0, Forces0 recalculated'
      CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr)
      CALL Build_Spectral_Forces_RE(NModes, ModesArr, rMomentum_K_RE, iProjectFlag)
      CALL Build_Spectral_Forces_IM(NModes, ModesArr, rMomentum_K_IM, iProjectFlag)
    EndIF
    IF (rLocalTimeCurr.GE.2.*ForceDecayTime.AND.rLocalTimePrevious.LT.2.*ForceDecayTime) THEN
      rLocalTimeCurr = 0.
      CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr2)
      CALL Build_Spectral_Forces_RE(NModes, ModesArr2, rMomentum_K_RE2, iProjectFlag)
      CALL Build_Spectral_Forces_IM(NModes, ModesArr2, rMomentum_K_IM2, iProjectFlag)
    EndIF
    rFactor1 = ABS(SIN(0.5*PICONST*(1.+rLocalTimeCurr/ForceDecayTime)))
    rFactor2 = ABS(SIN(0.5*PICONST*rLocalTimeCurr/ForceDecayTime))

!IF (IDCurrCPU.EQ.IDMasterCPU) &
!  Print *, 'rFactor1::', rFactor1, '  rFactor2::', rFactor2

    CALL Build_Real_Forces(NModes, ModesArr, &
                    rMomentum_K_RE, rMomentum_K_IM, PhasesArr, rFactor1, 1)
    CALL Build_Real_Forces(NModes, ModesArr2, &
                    rMomentum_K_RE2, rMomentum_K_IM2, PhasesArr, rFactor2, 0)
    CALL Find_Force_Factor
    iGravityCounter = 2571
  EndIF

!IF (IDCurrCPU.EQ.IDMasterCPU) &
!         Print *, 'Time::', time, '  rLocalTimeCurr::', rLocalTimeCurr



!!  PrevEInt = VAR_INTEGRATE('eint')
  PrevEInt = ZeroLevelOfInternalEnergy

  call SolveFriedmannEquation(time, dt)
  call dBaseSetProperty("Time", time+dt)

  IF (iMagneticFieldClean.EQ.1) CALL Clean_Magnetic_Field()
  call hydro_3d (sweep_order_xyz)               ! Hydrodynamics
  IF (iMagneticFieldClean.EQ.1) CALL Clean_Magnetic_Field()

!  CALL Find_MinMax("var4", rTemp1, rTemp2)
!  Print *, "CPU::", IDCurrCPU, ' Min::', rTemp1, ' Max::', rTemp2


  call source_terms(block_list, block_count, .false., dt )
                                                ! Apply source terms
                                                ! (e.g. nuclear burning)
  call RedshiftHydroQuantities                  ! Cosmological redshift
  call AdvanceRadiationField                    ! Radiation
  call AdvanceParticles                         ! Particles

  if (RunningParticles) call ReDistributeParticles(.false.)

  CurrEInt = VAR_INTEGRATE("eint")

  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_14", PrevEInt)
  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_15", CurrEInt)

  IF(InternalEnergyCorrectorFlag) THEN
     IF(PrevEInt/CurrEInt.LT.1.) THEN
       !!CALL Rescale_Variable("eint", PrevEInt/CurrEInt)
       CALL Rescale_Variable_2("eint", PrevEInt, CurrEInt)
     EndIF
  EndIF

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!               Second set of operator calls...

! Time update

  time = dBasePropertyReal("Time")
  dt   = dBasePropertyReal("TimeStep")

  rLocalTimePrevious = rLocalTimeCurr
  rLocalTimeCurr = rLocalTimeCurr + dt

  IF (iGravityCounter.EQ.2571) THEN
    IF (rLocalTimeCurr.GE.ForceDecayTime.AND.rLocalTimePrevious.LT.ForceDecayTime) THEN
      CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr)
      CALL Build_Spectral_Forces_RE(NModes, ModesArr, rMomentum_K_RE, iProjectFlag)
      CALL Build_Spectral_Forces_IM(NModes, ModesArr, rMomentum_K_IM, iProjectFlag)

    EndIF
    IF (rLocalTimeCurr.GE.2.*ForceDecayTime.AND.rLocalTimePrevious.LT.2.*ForceDecayTime) THEN
      rLocalTimeCurr = 0.
      CALL Build_Modes_Array(KMin, KMax, NModes, ModesArr2)
      CALL Build_Spectral_Forces_RE(NModes, ModesArr2, rMomentum_K_RE2, iProjectFlag)
      CALL Build_Spectral_Forces_IM(NModes, ModesArr2, rMomentum_K_IM2, iProjectFlag)
    EndIF
    rFactor1 = ABS(SIN(0.5*PICONST*(1.+rLocalTimeCurr/ForceDecayTime)))
    rFactor2 = ABS(SIN(0.5*PICONST*rLocalTimeCurr/ForceDecayTime))
!IF (IDCurrCPU.EQ.IDMasterCPU) &
!  Print *, 'rFactor1::', rFactor1, '  rFactor2::', rFactor2

    CALL Build_Real_Forces(NModes, ModesArr, &
                    rMomentum_K_RE, rMomentum_K_IM, PhasesArr, rFactor1, 1)
    CALL Build_Real_Forces(NModes, ModesArr2, &
                    rMomentum_K_RE2, rMomentum_K_IM2, PhasesArr, rFactor2, 0)
    CALL Find_Force_Factor
    iGravityCounter = 0
  EndIF

!IF (IDCurrCPU.EQ.IDMasterCPU) &
!         Print *, 'Time::', time, '  rLocalTimeCurr::', rLocalTimeCurr


!!  PrevEInt = VAR_INTEGRATE("eint")
  PrevEInt = ZeroLevelOfInternalEnergy

  call SolveFriedmannEquation(time, dt)
  call dBaseSetProperty("Time", time+dt)
! need to set the old timestep for consistency
  call dBaseSetProperty("OldTimeStep", dt)

  IF (iMagneticFieldClean.EQ.1) CALL Clean_Magnetic_Field()
  call hydro_3d (sweep_order_zyx)               ! Hydrodynamics
  IF (iMagneticFieldClean.EQ.1) CALL Clean_Magnetic_Field()
  call source_terms(block_list, block_count, .false., dt )
                                                ! Apply source terms
                                                ! (e.g. nuclear burning)
  call RedshiftHydroQuantities                  ! Cosmological redshift
  call AdvanceRadiationField                    ! Radiation
  call AdvanceParticles                         ! Particles
  if (RunningParticles) call ReDistributeParticles(.false.)

  CurrEInt = VAR_INTEGRATE("eint")
  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_14", PrevEInt)
  CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_15", CurrEInt)


  IF(InternalEnergyCorrectorFlag) THEN
     IF(PrevEInt/CurrEInt.LT.1.) THEN
       !!CALL Rescale_Variable("eint", PrevEInt/CurrEInt)
       CALL Rescale_Variable_2('eint', PrevEInt, CurrEInt)
     EndIF
  EndIF

!-------------------------------------------------------------------------------

! precompute dT/T using the old temperature -- this needs to be done before
! we refine
  if (itemp_limit .EQ. 1) call tstep_temperature_precompute()

!-------------------------------------------------------------------------------
!               Refine/derefine the AMR block structure as necessary.

  call mesh_update_grid_refinement()

!===============================================================================

rLocalTimePrevious = rLocalTimeCurr

!!! Added by Yu. Zaliznyak, 27JAN2005
CALL get_parm_from_context(GLOBAL_PARM_CONTEXT, "iSpectrSave", iSpectrSave)
iStepNumber = dBasePropertyInteger("CurrentStepNumber")
IF (MOD(iStepNumber, iSpectrSave).EQ.0) CALL Do_ALL_Spectral_Stuff


!! update the forcing if necessary
 IF (ForcingUpdateFlag) THEN
    CALL get_parm_from_context(global_parm_context, "Integral_1", rTemp1)
    CALL get_parm_from_context(global_parm_context, "Forcing_Energy_CURR", rTemp2)
    !!rTemp2 = 0.5*(Forcing_Energy_MIN +Forcing_Energy_MAX)
    rTemp3 = rTemp2

    IF (rTemp1.GT.MachMIN.AND.rTemp1.LT.MachMAX) &
       dEnergy = dEnergyBasic*COS(3.1415927*(rTemp1-MachMin)/(MachMAX-MachMIN))
    IF (rTemp1.GE.MachMAX) dEnergy = - dEnergyBasic
    IF (rTemp1.LE.MachMIN) dEnergy =   dEnergyBasic

    IF(rTemp1.GT.MachMiddle.AND.PreviousMachNumber.GT.rTemp1) dEnergy = 0.0
    IF(rTemp1.LT.MachMiddle.AND.PreviousMachNumber.LT.rTemp1) dEnergy = 0.0

    rTemp3 = rTemp2 + dEnergy

    IF (rTemp1.GT.MachOVERMAX)  rTemp3 = Forcing_Energy_MIN
    IF (rTemp1.LT.MachUNDERMIN) rTemp3 = Forcing_Energy_MAX

    IF (rTemp3.GT.Forcing_Energy_MAX) rTemp3 = Forcing_Energy_MAX
    IF (rTemp3.LT.Forcing_Energy_MIN) rTemp3 = Forcing_Energy_MIN

    CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Forcing_Energy_CURR", rTemp3)
    CALL set_parm_in_context(GLOBAL_PARM_CONTEXT, "Integral_13", rTemp3)

    IF (IDCurrCPU.EQ.IDMasterCPU) Print *, 'CPU::', IDCurrCPU, &
                ' MachFromDatabase::', rTemp1, ' CurrForcingEnergy::', rTemp3

 EndIF

  PreviousMachNumber = rTemp1

  CALL Conserve_Momentum()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return
  end

