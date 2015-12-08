This repository contains some examples of code written 
by Yuri Zaliznyak for FLASH hydro/mhd code developed in the University if Chicago
(for details and current state of the code see webpage http://flash.uchicago.edu/site/flashcode/).

All routines are in fortran90 and use MPI approach for parallel computations. Implemented, 
tested and employed for production runs with FLASH core version 2.4 on IBM p390 mainframe machine 
at Garching bei Munchen supercomputer centre in 2003-2005. Successfully scaled up to 512 cpus.

Main routine is Calc_Usr_Spectrum_Universal.F90. It contains spectral transformator, which first aggregate 
PARAMESH adaptive blocks structure into contiguous planes at each cpu, then applies parallel Fourier 
transform, and, finally, builds a spectrum of prescribed initial mhd variable and distribute it over 
cpus in a manner suitable for further processing.

Routine Other_Spectral_Routines.F90 contains supplementary subroutines called from Calc_Usr_Spectrum_Universal.

Routine stir.f90 implements consistent energy and momentum input into mhd system in a way similar to stirring tea
by spoon in a cup. Routine was used to study spectral properties of mhd forced stationary turbulence.