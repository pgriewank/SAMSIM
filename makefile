############################################################################################################
#Three compiler options
#-Nagfor is usefull for finding bugs, but can't compile minpack which is necesarry for optimization runs
#-Ifort is fast
#-Gfortran is free!
############################################################################################################
#Choose your compiler by commenting it out with #
############################################################################################################

#objects = mo_kind.o mo_parameters.o mo_data.o mo_functions.o mo_init.o mo_thermo_functions.o mo_mass.o mo_grav_drain.o mo_output.o mo_layer_dynamics.o mo_flush.o mo_snow.o mo_flood.o mo_heat_fluxes.o mo_testcase_specifics.o mo_grotz.o
objects = mo_parameters.o mo_data.o mo_functions.o mo_init.o mo_thermo_functions.o mo_mass.o mo_grav_drain.o mo_output.o mo_layer_dynamics.o mo_flush.o mo_snow.o mo_flood.o mo_heat_fluxes.o mo_testcase_specifics.o mo_grotz.o

#Nagfor
#Comp		  = nagfor
#FLAGS	  = -C=all -gline -nan -g


#ifort
#Comp		  = ifort
#FLAGS		 = -fast 

#gfortran
Comp		  = gfortran
FLAGS		 = -Wall -Wextra -fbounds-check



all: samsim.x

samsim.x : SAMSIM.f90 $(objects)
	$(Comp) ${FLAGS} -o $@   $(objects) SAMSIM.f90


%.o : %.f90
	$(Comp) ${FLAGS} -c  $<


#remove temporary files
clean:
	rm  *o *mod 


