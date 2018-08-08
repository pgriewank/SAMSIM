!>
!! @mainpage SAMSIM Semi-Adaptive Multi-phase Sea-Ice Model V1.0
!! 
!! This model was developed from scratch by Philipp Griewank during and after his PhD at  Max Planck Institute of Meteorology from 2010-2014.
!! Most elements of the model are described in my PhD thesis (A 1D model study of brine dynamics in sea ice, 2013, Hamburg) and in the publication "Insights into brine dynamics and sea ice desalination from a 1-D model study of gravity drainage" by Griewank & Notz 2013 JGR: Oceans. 
!!
!! Important changes from the version 1.0 to the previous version use for the my PhD thesis are: 
!! - the switch to a harmonic mean permeability for gravity drainage
!! - A few important bug fixes, which have a large affect on flushing
!!  
!! SAMSIM.f90 is the root program of the SAMSIM, the 1D thermodynamic Semi-Adaptive Multi-phase Sea-Ice Model.
!! The code is intended to be understandable and subroutines, modules, functions, parameters, and global variables all have doxygen compatible descriptions. 
!! However, in SAMSIM.f90 only the testcase and description thread are specified, which are then passed on to mo_grotz, which is where most of the actual work is done, including timestepping. 
!!
!! WARNING: SAMSIM was developed and used for scientific purposes. It surely contains at least some undetected bugs, can easily be crashed by using non-logical input settings, and some of the descriptions and comments may be outdated.  Always check the plausibility of the model results!
!!
!! Getting started:
!! - A number of testcases are implemented in SAMSIM. 
!! Testcases 1, 2, 3, and 4 are intended as standard testcases which should give a first time user a feel for the model capabilities and serve as a basis to set up custom testcases. 
!! To familiarize yourself with the model I suggest running testcases 1-3 and plotting the output with the python plotting scripts provided. 
!! The details of each testcase are commented in mo_init.f90, and each plot script begins with a  list of steps required. 
!! 
!! Running SAMSIM the first times.
!! - Make sure that all .f90 files are located in the same folder with the makefile.
!! - Open the makefile with your editor of choice and choose the compiler and flags of choice.
!! - Open SAMSIM.f90, set a testcase from 1-3, and edit the description string to fit your purpose. 
!! - Use make to compile the code, which produces the executable samsim.x .
!! - Make sure a folder "output" is located in the folder with samsim.x .
!! - Execute SAMSIM by running samsim.x .
!! - Go into output folder
!! - Copy the plot script from plotscripts to output 
!! - Follow the directions written in the plotscripts to plot the output.
!!
!! Running testcase 4.
!! In contrast to testcase 1-3, testcase four requires input files. 
!! Input data for testcase is provided in the input folder. 
!! Choose one of the subfolders from  input/ERA-interim/, copy the *.input files into the folder with the code, and run the executable .samsim.x .
!!
!! Following modules have a good documentation (both in the code and refman.pdf)
!! - mo_heat_fluxes.f90
!! - mo_layer_dynamics.f90
!! - mo_init.f90
!!
!! Biogeochemical tracers can be activated with  bgc_flag=2. 
!! - Warning! This feature is still relatively new, and has not been standardized yet. 
!! - The model will track Nbgc number of individual tracer. 
!! - Especially if you are interested in dissolved gases, you should first make yourself familiar with the bgc_advection subroutine in mo_mass.f90. 
!!
!! Know issues/Tips and Tricks:
!! - If code changes have no effect, run "make clean" and then "make", for unknown reasons this is often needed when making changes to mo_parameters.f90
!! - When bug hunting increase thick_0 and dt, this way the model runs faster, and the output is easier to sort through. 
!! - Use debug_flag= 2 to output data of each layer at each timestep. Be careful, the output size can become very large! 
!! - Check dat_settings to keep track of runs, and use the description variable to keep track of experiments.
!! - Contact me :)
!! 
!! @author Philipp Griewank
!!
!!
!!
!!  COPYRIGHT
!!
!! This file is part of SAMSIM.
!!
!!  SAMSIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!  SAMSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License along with SAMSIM. If not, see <http://www.gnu.org/licenses/>.
!!
!! @par Revision History
!! Started by Philipp Griewank 2014-05-05
!! nothing changed here by Niels Fuchs, MPIMET (2017-03-01)
!! License changed by Philipp Griewank 2018-05-22
!!
PROGRAM SAMSIM
  USE mo_grotz 



  IMPLICIT NONE
  INTEGER:: testcase
  CHARACTER*12000                       :: description   !< String to describe simulation which is outputted into dat_settings
  
  !##########################################################################################
  !Initialization
  !##########################################################################################
  testcase    = 1                  !sets the testcase
  description = 'ttest_v2'  !is written to dat_settings 
  !##########################################################################################
    
  PRINT*,'SAMSIM is getting ready'
  CALL grotz (testcase,description)
  PRINT*,'SAMSIM is finished'



END PROGRAM SAMSIM


