!!****f* source/driver/time_dep/flash
!!
!! NAME
!!
!!  flash
!!
!!
!! SYNOPSIS
!!
!!  main
!!
!!
!! DESCRIPTION
!!
!!  This is the main FLASH driver for time-dependent problems.  
!!  Initiializations are performed, and then the main time step
!!  loop is begun.  After each timestep, output is called to
!!  dump out any checkpoint or plot files, and to write out
!!  global sums to flash.dat.  
!!
!!
!! NOTES
!!
!!  Please see flash.F90 in the main driver directory for a 
!!  description of FLASH and the author list
!!
!!***

program flash

  use runtime_parameters 
 
  use perfmon

  use logfile, ONLY: stamp_logfile

  use dBase, ONLY: dBasePropertyInteger, dBasePropertyReal
  
  use dBaseDeclarations, ONLY: nstep

  implicit none

#include "mpif.h"

  integer :: myPE, masterPE, nbegin, nend, nstep_
  real    :: time, tmax, dt, redshift, scale, zfinal

  real :: wall_clock_time_start
  real :: wall_clock_time_limit

  logical :: endtime
  integer :: ierr

! for logfile output
  character(len=80), dimension(3,2) :: str_buff
  character(len=80), dimension(2,2) :: str_buff2
  character(len=15) :: num_to_str
!===============================================================================

! Initialize the simulation.

  call init_flash()

! Initialize the time that the simulation started -- we will use this
! for a wall clock limit on the simulation, so we can fit it in a queue

  wall_clock_time_start = MPI_Wtime()

  myPE     = dBasePropertyInteger("MyProcessor")
  masterPE = dBasePropertyInteger("MasterProcessor")
  nbegin   = dBasePropertyInteger("BeginStepNumber")

! Note that these calls must come after first call to init_flash

  call get_parm_from_context("tmax", tmax)
  call get_parm_from_context("zfinal", zfinal)
  call get_parm_from_context("nend", nend)
  call get_parm_from_context("wall_clock_time_limit", wall_clock_time_limit)

!===============================================================================

! Main timestep iteration.

  call stamp_logfile ("[FLASH] Enter evolution loop..." )

  call timer_start ("evolution")

  do nstep_ = nbegin, nbegin+nend-1

     ! make the step number independent on the way loop terminates: exit works fine
     ! but passing through do-loop end increases the do-loop counter by extra 1

     nstep    = nstep_
     time     = dBasePropertyReal("Time")
     dt       = dBasePropertyReal("TimeStep")
     redshift = dBasePropertyReal("Redshift")
     scale    = dBasePropertyReal("ScaleFactor")

     if (MyPE == MasterPE) then
       write (num_to_str(1:), '(I10)') nstep
       write (str_buff(1,1), "(A)") "n"
       write (str_buff(1,2), "(A)") trim(adjustl(num_to_str))

       write (num_to_str(1:), "(1PE12.6)") time
       write (str_buff(2,1), "(A)") "t"
       write (str_buff(2,2), "(A)") trim(adjustl(num_to_str))

       write (num_to_str(1:), "(1PE12.6)") dt
       write (str_buff(3,1), "(A)") "dt"
       write (str_buff(3,2), "(A)") trim(adjustl(Num_to_str))

       call stamp_logfile (str_buff, 3, 2, "step")
     end if

     if ((MyPE == MasterPE) .and. (redshift /= 0.)) then
       write (num_to_str(1:), '(ES14.6)') redshift
       write (str_buff2(1,1),'(A)') 'z'
       write (str_buff2(1,2),'(A)') trim(adjustl(num_to_str))

       write (num_to_str(1:), '(ES14.6)') scale
       write (str_buff2(2,1),'(A)') 'a'
       write (str_buff2(2,2),'(A)') trim(adjustl(num_to_str))

       call stamp_logfile (str_buff2, 2, 2, 'info', 'type="list" count="0"')
     endif

     call evolve()            ! Evolve system one timestep

     ! The time was just updated, so grab the new one

     time     = dBasePropertyReal("Time")
     redshift = dBasePropertyReal("Redshift")

     call timestep()          ! Compute new timestep

     ! dt was just updated, so get the new value

     dt       = dBasePropertyReal("TimeStep")

     call output(time, dt, nstep+1, nbegin) ! Output if needed

     ! Call visualization routines

     call visualization()

     ! Termination criteria: the simulation time is greater than the
     ! maximum time (tmax), or if the simulation redshift is less than the
     ! stopping redshift (zfinal), or if more than wall_clock_time_limit seconds
     ! have elapsed since the start of the simulation.

     if ( time >= tmax       ) exit
     if ( redshift <= zfinal ) exit

     endtime = .false.

     if (myPE == masterPE) then
         if (MPI_Wtime() - wall_clock_time_start > wall_clock_time_limit) endtime = .true.
     endif

     call mpi_bcast(endtime, 1, MPI_LOGICAL, masterpe, MPI_COMM_WORLD, ierr)

     if (endtime) exit


  enddo

  call timer_stop ("evolution")

  call stamp_logfile ("[FLASH] Exit evolution loop..." )

!===============================================================================

! End of time iteration.  Clean up and terminate.

  call end_flash (nstep-nbegin+1)

!===============================================================================

end program flash
