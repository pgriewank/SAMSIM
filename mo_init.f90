!>
!! Allocates Arrays and sets initial data for a given testcase for SAMSIM.
!!
!!
!! @author Philipp Griewank
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
!
!!
!!
!! @par Revision History
!! first version created to deal with first multi-layer tests. by Philipp Griewank, IMPRS (2010-07-22)
!! Add Testcases: 101-105 are simulations of master theses experiments 1-5, 111 can be used to compare SAMSIM with salinity harps field data by Niels Fuchs, MPIMET (2017-03-01)
!!
MODULE mo_init

  USE mo_parameters
  USE mo_data

  IMPLICIT NONE

CONTAINS
  !>
  !! Sets initial conditions according to which testcase is chosen
  !!
  !! For different initial conditions the Arrays are allocated and the initial values are set.
  !! Following must always be:
  !! 1. Nlayer = N_top+N_middle+N_bottom
  !! 2. N_active is set correctly, N_active <= Nlayer
  !! 3. fl_q_bottom >= 0
  !! 4. T_bottom > freezing point of for S_bu_bottom 
  !! 5. A too high dt for a too small thick_0 leads to numerical thermodynamic instability. For a conservative guess dt [s] should be smaller than 250000 * (dz [m])**2 
  !!
  !! Testcase 1
  !! - Testcase 1 is a replication of lab experiments conducted in tanks cooled from above by a cooling plate using the boundflux_flag 1. 
  !! - In this testcase the cooling plate Temperature T_top changes every 12 hours to imitate the experiments Dirk Notz conducted in his PhD.
  !! - This testcase was used to optimize the free parameters of the gravity drainage parametrization (see Griewank Notz 2013/14).
  !! - Can also be run with bgc tracers.
  !!
  !!
  !! Testcase 2
  !! - Testcase is an example of how to simulate ice growth and melt in cooling chambers.
  !! - Boundflux_flag 3 is used, which uses T2m as the air temperature in the cooling chamber.
  !! - The surface flux heat flux is proportional to the ice-air temperature difference (T_top-T2m).
  !! - When reproducing cooling chamber experiments the alpha flux parameters need to be tuned, and a module in mo_testcase_specifics is needed to set/ T2m over time.
  !! - The heat flux in the water from below (fl_q_bottom) for such experiments can be very hard to reproduce if the heat input is not carefully measured from all pumps or similar devices used. 
  !!
  !! Testcase 3
  !! - Uses interpolated climate mean forcing from Notz and a constant oceanic heat flux (fl_q_bottom) to grow idealized arctic sea ice. 
  !! - Is generally intended as a numerically cheap testcase to check for effects of code changes. 
  !! - Is also useful when runs over many years are needed. 
  !! - The amount of liquid and solid precipitation is set in sub_test3 of mo_testcase specifics. 
  !!
  !! Testcase 4
  !! - Uses three hourly reanalysis forcing over 4.5 years. 
  !! - Is set up to start in July.
  !! - Prescribes annual cycle of oceanic heat flux.
  !! - Requires the proper input data to be copied into the executable folder (see sub_input). 
  !! - Is more computer intensive
  !! - Was used a lot for Griewank & Notz 2013/2014 
  !!  
  !!
  !! @par Revision History
  !! First set up by Philipp Griewank, IMPRS (2010-07-22>)
  SUBROUTINE init (testcase)


    INTEGER, INTENT(in)     :: testcase 
    INTEGER :: i,jj

    !##########################################################################################
    !#######################    DEFAULT SETTINGS   ############################################
    !##########################################################################################
    !________________________top heat flux____________
    boundflux_flag  = 1
    atmoflux_flag   = 1
    albedo_flag     = 2
    !________________________brine_dynamics____________
    grav_heat_flag  = 1
    flush_heat_flag = 1
    flood_flag      = 2
    flush_flag      = 5
    grav_flag       = 2
    harmonic_flag   = 2
    !________________________Salinity____________
    prescribe_flag  = 1
    salt_flag       = 1
    !________________________bottom setting______________________
    turb_flag       = 2
    bottom_flag     = 1
    tank_flag       = 1
    !________________________snow______________________
    precip_flag     = 0
    freeboard_snow_flag = 0    !< Niels, 2017
    snow_flush_flag = 0    !< Niels, 2017
    snow_precip_flag = 1   !< Niels, 2017
    !________________________debugging_____________________
    debug_flag      = 1   !set to 2 for output of all ice layers each timestep
    !________________________bgc_______________________
    bgc_flag        = 1
    N_bgc           = 1


    
  !##########################################################################################
  !Here I just set a lot of things to zero, arrays are set to zero after they are allocated
  !##########################################################################################
    thick_snow    = 0.0_wp
    m_snow        = 0.0_wp
    psi_g_snow    = 0.0_wp
    psi_l_snow    = 0.0_wp
    psi_s_snow    = 0.0_wp
    H_abs_snow    = 0.0_wp
    S_abs_snow    = 0.0_wp
    phi_s         = 0.0_wp
    T_snow        = 0.0_wp
    liquid_precip = 0.0_wp
    solid_precip  = 0.0_wp
    fl_sw         = 0.0_wp
    fl_rest       = 0.0_wp
    albedo        = 0.0_wp
    T_top         = 0.0_wp
    T2m           = 0.0_wp
    fl_q_bottom   = 0.0_wp


    !################################################################################
    !TESTCASES
    !I consider testcases 1-4 the default testcases. I advise leaving those
    !intact to check how code changes affect them, and making your own for all
    !other purposes.
    !################################################################################
    IF (testcase==111) THEN !##########################################################
       !!< Niels, 2017 add: Testcase 111 can be used to compare SAMSIM with salinity harp field measurements. The uppermost temperature sensor of
       ! the harps provides the temperature T_Top. Was used for Greenland data in 2017.
       
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer         = 100
       N_active       = 1
       N_top          = 10
       N_bottom       = 10
       N_middle       = Nlayer-N_top-N_bottom
       length_input_lab = 860333
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       turb_flag      = 1 
       boundflux_flag = 1
       grav_heat_flag = 1
       flush_flag     = 1
       salt_flag      = 2


       !*************************************************************************************************************************
       !Salinity, temperature, and bottom heatflux settings
       !*************************************************************************************************************************
       T_top       =-2.0_wp
       T_bottom    =-1.67_wp
       S_bu_bottom = 33.4079_wp
       fl_q_bottom = 0._wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0     = 0.01_wp
       dt          = 3.0_wp !0.5_wp
       time       = 0.0_wp
       time_out   = 3600._wp*2._wp
       time_total = 2580996._wp

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m(1)     = thick(1)*rho_l
       S_abs(1) = S_bu_bottom*m(1)
       H_abs(1) = m(1)*(T_bottom)*c_l
       
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************
       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          !*************************************************************************************************************************
          !Setting number of tracers
          !*************************************************************************************************************************
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)


          !*************************************************************************************************************************
          !Setting bottom concentrations
          !*************************************************************************************************************************
          bgc_bottom(1) = 400._wp
          bgc_bottom(2) = 500._wp

          !*************************************************************************************************************************
          !setting the initial values of the top and only layer
          !*************************************************************************************************************************
          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
        end if

    
       
   ELSE IF (testcase==101) THEN !#########################################################################
       !Testcase by Niels Fuchs, Master Thesis Experiment 1
       !Simulation setup used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 0._wp           !Heat flux from below
       alpha_flux_instable = 22._wp !30._wp !22._wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18 std:22.0
       alpha_flux_stable   = 21._wp !29._wp !21._wp    !Chosen at random std=15
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 200
       N_bottom = 10
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       length_input_lab = 1628263
       
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 5 !1 !5
       flood_flag = 2 !2
       grav_flag = 2
       lab_snow_flag = 1
       freeboard_snow_flag = 1
       snow_flush_flag = 1
       flush_heat_flag = 2
       snow_precip_flag = 1

       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = 0._wp
       T_top       = 0._wp
       T_bottom    = -1.3_wp !-1.8_wp!0.0
       S_bu_bottom = 25.6664555556_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 60._wp * 60._wp
       time_total = 1625000._wp
       dt         = 1._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
   
   ELSE IF (testcase==102) THEN !#########################################################################
       !Testcase by Niels Fuchs, Master Thesis Experiment 2
       !Simulation setup used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 
	
       fl_q_bottom         = 0._wp           !Heat flux from below
       alpha_flux_instable = 22._wp !22._wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18 std:22.0
       alpha_flux_stable   = 21._wp !21._wp    !Chosen at random std=15
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 200
       N_bottom = 10
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       length_input_lab = 1124187
       
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 5 !1 !5
       flood_flag = 2 !2
       grav_flag = 2
       lab_snow_flag = 1
       freeboard_snow_flag = 1
       snow_flush_flag = 1
       flush_heat_flag = 2
       snow_precip_flag = 1

       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = 0._wp
       T_top       = 0._wp
       T_bottom    = -1.3_wp !-1.8_wp!0.0
       S_bu_bottom = 26.1336777778_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 60._wp * 60._wp
       time_total = 1124000._wp
       dt         = 1._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
   
   
   
   ELSE IF (testcase==103) THEN !#########################################################################
       !Testcase by Niels Fuchs, Master Thesis Experiment 3
       !Simulation setup used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 0._wp           !Heat flux from below
       alpha_flux_instable = 22._wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18 std:22.0
       alpha_flux_stable   = 21._wp    !Chosen at random std=15
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 200
       N_bottom = 10
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       length_input_lab = 1283092
       
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 5 !1 !5
       flood_flag = 2 !2
       grav_flag = 2
       lab_snow_flag = 1
       freeboard_snow_flag = 1
       snow_flush_flag = 1
       flush_heat_flag = 2
       snow_precip_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = 0._wp
       T_top       = 0._wp
       T_bottom    = -1.3_wp !-1.8_wp!0.0
       S_bu_bottom = 26.0335888889_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 60._wp * 60._wp
       time_total = 1283000._wp
       dt         = 1._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
       
   ELSE IF (testcase==104) THEN !#########################################################################
       !Testcase by Niels Fuchs, Master Thesis Experiment 4
       !Simulation setup used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 0._wp           !Heat flux from below
       alpha_flux_instable = 22._wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18 std:22.0
       alpha_flux_stable   = 21._wp    !Chosen at random std=15
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 200
       N_bottom = 10
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       
       length_input_lab = 2439729
       
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 5 !1 !5
       flood_flag = 2 !2
       grav_flag = 2
       lab_snow_flag = 1
       freeboard_snow_flag = 1
       snow_flush_flag = 1
       flush_heat_flag = 2
       snow_precip_flag = 1

       

       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = 0._wp
       T_top       = 0._wp
       T_bottom    = -1.3_wp !-1.8_wp!0.0
       S_bu_bottom = 27.0363_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 60._wp * 60._wp
       time_total = 2439000._wp
       dt         = 1._wp
       
       
       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
       
   ELSE IF (testcase==105) THEN !#########################################################################
       !Testcase by Niels Fuchs, Master Thesis Experiment 5
       !Simulation setup used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 0._wp           !Heat flux from below
       alpha_flux_instable = 22._wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18 std:22.0
       alpha_flux_stable   = 21._wp    !Chosen at random std=15
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 200
       N_bottom = 10
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       
       length_input_lab = 1549323
       
       CALL sub_allocate(Nlayer,length_input_lab)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 5 !1 !5
       flood_flag = 2 !2
       grav_flag = 2
       lab_snow_flag = 1
       freeboard_snow_flag = 1
       snow_flush_flag = 1
       flush_heat_flag = 2
       snow_precip_flag = 1

       

       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = 0._wp
       T_top       = 0._wp
       T_bottom    = -1.3_wp !-1.8_wp!0.0
       S_bu_bottom = 31.5625333333_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 60._wp * 60._wp
       time_total = 1549000._wp
       dt         = 1._wp
       
       
       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
       
   ELSE IF (testcase==99) THEN !#########################################################################
       !Testcase by Niels Fuchs, Testcase for snow on ice in chamber experiments, simple parametrization without
       !changes in the code so far, 2016/05/17

       fl_q_bottom         = 5._wp           !Heat flux from below
       alpha_flux_instable = 22.0_wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 15._wp    !Chosen at random 
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 20
       N_bottom = 5
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       precip_flag =  0
       grav_heat_flag = 1
       flush_flag = 1
       flood_flag = 1
       grav_flag = 2


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -5._wp
       T_top       = -2._wp
       T_bottom    = -1.8_wp !-1.8_wp!0.0
       S_bu_bottom = 34_wp
       
       !*************************************************************************************************************************
       !Initial Snow settings
       !*************************************************************************************************************************
       
       
       

       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.05_wp
       time_out   = 60._wp * 10._wp
       time_total = 3600._wp*24._wp*7._wp
       dt         = 10._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
    
    ELSE IF (testcase==1) THEN !##########################################################
       !Testcase 1 is a replication of lab experiments conducted in tanks
       !cooled from above by a cooling plate using the boundflux_flag 1. 
       !In this testcase the cooling plate Temperature T_top changes every 12 hours to imitate the experiments Dirk Notz conducted in his PhD.
       !This testcase was used to optimize the free parameters of the gravity drainage parametrization (see Griewank Notz 2013/14).
       !Can also be run with bgc tracers
       
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer         = 90
       N_active       = 1
       N_top          = 5
       N_bottom       = 5
       N_middle       = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       turb_flag      = 1 
       boundflux_flag = 1
       grav_heat_flag = 1
       flush_flag     = 1
       salt_flag      = 2


       !*************************************************************************************************************************
       !Salinity, temperature, and bottom heatflux settings
       !*************************************************************************************************************************
       T_top       =-5.0
       T_bottom    =-1.!735
       S_bu_bottom = 34.
       fl_q_bottom = 0._wp          


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0     = 0.002_wp
       dt          = 1.0_wp !0.5_wp
       time       = 0.0_wp
       time_out   = 3600._wp
       time_total = time_out*72._wp


       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m(1)     = thick(1)*rho_l
       S_abs(1) = S_bu_bottom*m(1)
       H_abs(1) = m(1)*(T_bottom)*c_l
       
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************
       bgc_flag    = 2

       IF (bgc_flag==2) THEN
          !*************************************************************************************************************************
          !Setting number of tracers
          !*************************************************************************************************************************
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)


          !*************************************************************************************************************************
          !Setting bottom concentrations
          !*************************************************************************************************************************
          bgc_bottom(1) = 400._wp
          bgc_bottom(2) = 500._wp

          !*************************************************************************************************************************
          !setting the initial values of the top and only layer
          !*************************************************************************************************************************
          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
        end if


    ELSE IF (testcase==2) THEN !#########################################################################
       !Simulation setup mostly used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 10._wp           !Heat flux from below
       alpha_flux_instable = 22.0_wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 15._wp    !Chosen at random 
       tank_depth          = 1._wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 10
       N_top    = 3
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       grav_heat_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -20._wp
       T_top       = -18._wp
       T_bottom    = 0.0_wp!-1.5
       S_bu_bottom = 31.2_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time_out   = 3600._wp*6._wp
       time_total = time_out*4._wp*30._wp
       dt         = 30._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 2

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF



    ELSE IF (testcase==3) THEN !###############################################################################################
       !Uses interpolated climate mean forcing from Notz and a constant oceanic
       !heat flux
       !Is generally intended as a way to take a quick look at possible changes
       !or for runs over many many years where repeating and smooth forcing is
       !wanted. 
       
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 20
       N_bottom = 5
       N_top    = 5
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       atmoflux_flag  = 1
       precip_flag    = 0
       boundflux_flag = 2
       !flush_flag     = 4
       !grav_flag      = 1
       !grav_heat_flag = 2
       !prescribe_flag = 2

       !*************************************************************************************************************************
       !Salinity, temperature, and oceanic heatflux settings
       !*************************************************************************************************************************
       fl_q_bottom  = 8._wp
       T_bottom     = -1.0
       S_bu_bottom  = 34._wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.03_wp
       time       = 0.0_wp
       time_out   = 86400._wp*3.5_wp
       time_total = time_out*54._wp*2._wp*2._wp!0._wp
       dt         = 60._wp


       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = 0._wp
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************
       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          !*************************************************************************************************************************
          !Setting number of tracers
          !*************************************************************************************************************************
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          !*************************************************************************************************************************
          !Setting bottom concentrations
          !*************************************************************************************************************************
          bgc_bottom(1) = 400._wp
          bgc_bottom(2) = 500._wp

          !*************************************************************************************************************************
          !setting the initial values of the top and only layer
          !*************************************************************************************************************************
          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
        end if



    ELSE IF (testcase==4) THEN !###############################################################################################
       !Imports lw and sw fluxes, T2m, and precip information form reanalysis data. 
       !Data has to be prepared as ascii data in the right format.  
       !The run depends on which input files are in the code folder. 
       !This test is the setup used to produce "realistic" ice in Griewank and Notz 2013,2014
       !Runs begin in July and runs for 4.5 years.  
       !A simple annual heat flux is used for the oceanic heat flux. 
       
       
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 20
       N_top    = 20
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       atmoflux_flag  = 2
       precip_flag    = 1
       boundflux_flag = 2
       snow_flush_flag = 0
       flush_heat_flag = 2
       snow_precip_flag = 1


       !*************************************************************************************************************************
       !Salinity, temperature, and oceanic heatflux settings
       !*************************************************************************************************************************
       T_bottom    = -1.0
       S_bu_bottom = 34._wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 86400._wp
       time_total = time_out*365._wp*4.5_wp
       dt         = 10._wp



       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = 0._wp
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************
       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          !*************************************************************************************************************************
          !Setting number of tracers
          !*************************************************************************************************************************
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          !*************************************************************************************************************************
          !Setting bottom concentrations
          !*************************************************************************************************************************
          bgc_bottom(1) = 400._wp
          bgc_bottom(2) = 500._wp

          !*************************************************************************************************************************
          !setting the initial values of the top and only layer
          !*************************************************************************************************************************
          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
        end if


    ELSE IF (testcase==5) THEN !#####################################################################
       !Set to test top melt
       !Melts a two meter thick block of ice with a fixed salinity
       !Important for paper 2!
       !bottom_flag    = 1!2
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_active = Nlayer
       N_bottom = 10!3
       N_top    = 20!3
       N_middle = Nlayer-N_top-N_bottom

       CALL sub_allocate(Nlayer)


       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       boundflux_flag = 2
       atmoflux_flag  = 3
       flush_heat_flag = 2

       !Complex
       flush_flag = 5
       grav_flag  = 2
       grav_flag  = 1
       flood_flag = 1
 

       !*************************************************************************************************************************
       !Heat fluxes
       !*************************************************************************************************************************
       fl_sw       = 0._wp
       fl_rest     = 290._wp**4*sigma!287._wp**4 *sigma 
       fl_q_bottom = 15._wp!15._wp
       
       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       S_bu_bottom = 5.00000_wp
       T_bottom = 0.

       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0 = 1._wp/500._wp*2._wp
       thick_0 = 0.01_wp 
       time    = 0.0_wp
       time_out   = 3600._wp*3._wp 
       time_total = time_out*24._wp*10._wp !30 war original
       dt         = 1._wp
       dt         = 10._wp


       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       
       thick   = thick_0
       m       = thick*rho_l
       S_abs   = m*S_bu_bottom
       H_abs   = m*(-90.0)*c_l !80 was the original value used in the paper




    ELSE IF (testcase==6) THEN !#########################################################################
       fl_q_bottom  = 35._wp  !(unknown heat flux from pumps)
       alpha_flux_instable = 22.0_wp   !Is tuned to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 11._wp   
       tank_depth          = 0.159_wp
       
       

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 40
       N_bottom = 3
       N_top    = 3
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       grav_heat_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -18._wp
       T_top       = -18._wp
       T_bottom    = 0.0_wp!-1.5
       S_bu_bottom = 31.2_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0  = 0.0025_wp
       time_out   = 1800._wp/2._wp
       time_total = time_out*39._wp*2._wp*2._wp
       dt         = 0.5

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 2

       IF (bgc_flag==2) THEN
          N_bgc = 1
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          !Oxygen
          bgc_bottom(1) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
       END IF
    ELSE IF (testcase==7) THEN !###############################################################################################
       !Identical to 4, but used for comparisons of thermal properties
       !Which means grav_heat_flag and flush_heat_flag are changed, and
       !albedo_flag is set to 1 to ensure albedo feedbacks don't affect things. 
       
       
       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 20
       N_top    = 20
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       atmoflux_flag  = 2
       precip_flag    = 1
       boundflux_flag = 2
       albedo_flag     = 1
       grav_heat_flag  = 2
       flush_heat_flag = 2

       !prescribe
       !flush_flag = 4
       !grav_flag  = 1
       !flood_flag = 1
       !prescribe_flag = 2

       !simple
       flush_flag = 4
       grav_flag  = 3
       flood_flag = 3


       !*************************************************************************************************************************
       !Salinity, temperature, and oceanic heatflux settings
       !*************************************************************************************************************************
       T_bottom    = -1.0
       S_bu_bottom = 34._wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.01_wp
       time       = 0.0_wp
       time_out   = 86400._wp/2._wp
       time_total = time_out*365._wp*9_wp
       dt         = 10._wp


       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = 0._wp
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************
       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          !*************************************************************************************************************************
          !Setting number of tracers
          !*************************************************************************************************************************
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          !*************************************************************************************************************************
          !Setting bottom concentrations
          !*************************************************************************************************************************
          bgc_bottom(1) = 400._wp
          bgc_bottom(2) = 500._wp

          !*************************************************************************************************************************
          !setting the initial values of the top and only layer
          !*************************************************************************************************************************
          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 
        end if



    ELSE IF (testcase==8) THEN !#########################################################################
       !Imports temperatures from the field data. 
       !T_top is read in from Tinput.txt, which must have been set to fit with dt 
       !WARNING: settings are likely outdated

       Nlayer   = 50
       N_active = 1
       N_bottom = 5
       N_top    = 4
       N_middle = Nlayer-N_top-N_bottom

       T_top       = -5._wp
       T_bottom    = -1.8
       T_test      = -1._wp
       S_bu_bottom = 34._wp

       CALL sub_allocate(Nlayer)

       boundflux_flag = 1
       fl_q_bottom    = 15._wp


       grav_flag  = 2
       flush_flag = 5
       flood_flag = 2


       thick_0  = 0.005_wp
       m        = 0._wp
       S_abs    = 0._wp
       H_abs    = 0._wp
       thick    = 0._wp
       thick(1) = thick_0
       m(1)     = thick(1)*rho_l
       S_abs(1) = S_bu_bottom*m(1)
       H_abs(1) = m(1)*(T_bottom)*c_l

       time       = 0.0_wp
       time_out   = 3600._wp!*12._wp
       time_total = time_out*12._wp*12._wp!138._wp
       dt         = 1.0_wp





    ELSE IF (testcase==50) THEN !#########################################################################
       !See Griewank Notz 2012
       !Used to generate stable initial conditions for convection either due to top warming or/and bottom melt

       fl_sw          = 0.0_wp
       boundflux_flag = 2

       fl_rest        = sigma*(zeroK-20._wp)**4._wp
       fl_q_bottom    = 20.0_wp

       Nlayer   = 70
       N_bottom = 5
       N_top    = 5
       N_middle = Nlayer-N_top-N_bottom

       T_top       = -20._wp
       T_bottom    = -1.72_wp
       T_test      = -1._wp
       S_bu_bottom = 34._wp
       N_active    = 1

       CALL sub_allocate(Nlayer)

       thick_0  = 0.005_wp

       thick(1) = thick_0
       m(1)     = thick(1)*rho_l
       S_abs(1) = S_bu_bottom*m(1)
       H_abs(1) = m(1)*(T_bottom)*c_l


       !__________________________Time Settings__________________________________________________________________________________
       time       = 0.0_wp
       time_out   = 3600._wp*24._wp*30._wp
       dt         = 10.0_wp
       time_total =time_out*12._wp*3._wp

    ELSE IF (testcase==51) THEN !#########################################################################
       !Testcase to generate convection from the initial conditions generated by
       !testcase 50.
       !See Griewank Notz 2012

       !Stability values
       fl_rest=sigma*(zeroK-20._wp)**4._wp
       fl_q_bottom =20.0_wp
       T_top    =-16.7_wp

       fl_sw=0.0_wp
       boundflux_flag=1
       turb_flag=1 

       !For melt tests
       flush_flag=5
       boundflux_flag=2
       grav_flag=2


       fl_rest=sigma*(zeroK+10._wp)**4._wp


       !Bottom warming values
       !fl_q_bottom =100.0_wp

       !Top warming value
       !fl_rest=sigma*(zeroK-10._wp)**4._wp
       !T_top =-5._wp


       Nlayer   = 70
       N_bottom = 5
       N_top    = 5
       N_middle = Nlayer-N_top-N_bottom

       T_bottom    = -1.72_wp
       T_test      = -1._wp
       S_bu_bottom = 34._wp
       N_active    = 70

       CALL sub_allocate(Nlayer)

       thick_0 = 0.01_wp

       thick= (/ 1.000000000000000E-002, 1.000000000000000E-002, 1.000000000000000E-002, &
            1.000000000000000E-002, 1.000000000000000E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 2.433333333333330E-002, &
            2.433333333333330E-002, 2.433333333333330E-002, 1.000000000000000E-002, &
            1.000000000000000E-002, 1.000000000000000E-002, 1.000000000000000E-002, &
            1.000000000000000E-002 /)
       H_abs=  (/ -3067803.68479882,   -3157066.88222312,   -3245320.33966741,&
            -3279773.33899002,   -3287465.23823668,   -7982644.18078968,& 
            -7983295.60211406,   -7974144.08828951,   -7961801.99509219,& 
            -7950114.81389292,   -7939225.05603074,   -7928882.39340588,& 
            -7919122.23011157,   -7909879.66329468,   -7900896.94648847,& 
            -7891878.51063782,   -7882619.36868540,   -7873028.35565817,& 
            -7863088.05602079,   -7852808.88216294,   -7842202.94107025,& 
            -7831274.84587167,   -7820019.82118742,   -7808424.47653668,& 
            -7796469.99130447,   -7784137.98764326,   -7771418.06645274,& 
            -7758314.78941944,   -7744851.50834870,   -7731068.93060507,& 
            -7717017.64975141,   -7702746.00575044,   -7688286.98524375,& 
            -7673649.06388404,   -7658814.62185122,   -7643746.04933834,& 
            -7628396.01512377,   -7612717.00471593,   -7596666.75539792,& 
            -7580209.06470259,   -7563311.51912057,   -7545942.09056554,& 
            -7528065.86048900,   -7509642.38021483,   -7490623.82950644,& 
            -7470953.99330339,   -7450567.83564391,   -7429391.17550865,& 
            -7407339.99686197,   -7384319.35256003,   -7360221.87604925,& 
            -7334923.70692784,   -7308269.51461524,   -7280040.12907148,& 
            -7249953.56216598,   -7217405.54376676,   -7181461.64725477,& 
            -7141061.72081779,   -7094766.84356921,   -7040432.37369522,& 
            -6974656.23455786,   -6891692.47650013,   -6781195.01997001,& 
            -6618939.74338630,   -6375126.79490718,   -2507755.81456431,& 
            -2381863.95043677,   -2128107.88136274,   -537715.435477377,& 
            -58736.8820714995 /)
       S_abs= (/  192.926263213596, 134.117690286622, 76.3530714340286,&   
            53.1728736974625, 47.1042506294751, 120.364319404813,&   
            113.025150825946, 111.866418637160, 112.667470674304,&   
            113.029332437771, 112.874754124485, 112.365709288828,&   
            111.488808633605, 110.290845483411, 108.930490214493,&   
            107.583123838865, 106.365739721918, 105.325934349434,&   
            104.467523396504, 103.778220161281, 103.244881509094,&   
            102.858595517385, 102.615491783745, 102.516041543011,&   
            102.563089234326, 102.758571394302, 103.099576228544,&   
            103.574968190804, 104.163890649505, 104.837081774730,&   
            105.561155807062, 106.304900144283, 107.045560297543,&   
            107.772705388992, 108.488132214214, 109.202107418496,&   
            109.927935520001, 110.677270219149, 111.457638880733,&   
            112.272210797958, 113.120914489105, 114.001923559312,&   
            114.912925252744, 115.851964411225, 116.817825625550,&   
            117.809984400219, 118.828252619582, 119.872330382230,&   
            120.941431746516, 122.033969093113, 123.147256819261,&   
            124.277905641201, 125.425304063697, 126.599300649984,&   
            127.815886165244, 129.165537852646, 130.795740225498,&   
            132.836931453737, 135.453554058069, 138.892330435293,&   
            143.552140461959, 150.115395232649, 159.805508716573,&   
            175.686084197530, 200.163341017208, 93.8743106256185,&   
            108.336464433014, 138.880147601434, 334.426718007957,&   
            349.545260822269    /) 

       m =(/ 9.30703496197182,  9.27463917806153,  9.24262135776094 ,&  
            9.22977109352861,  9.22645249692149,  22.4546108250397 ,&  
            22.4509423829188,  22.4507622333214,  22.4517137893723 ,&  
            22.4524280806754,  22.4528547936790,  22.4530833331591 ,&  
            22.4531014113224,  22.4529324045383,  22.4526681662252 ,&  
            22.4524138356116,  22.4522422072450,  22.4521855355380 ,&  
            22.4522498356824,  22.4524312332869,  22.4527252582962 ,&  
            22.4531301070479,  22.4536473097727,  22.4542814501597 ,&  
            22.4550390299766,  22.4559263825359,  22.4569470004799 ,&  
            22.4580990652460,  22.4593740991186,  22.4607574858411 ,&  
            22.4622311290298,  22.4637777603522,  22.4653855792912 ,&  
            22.4670514872595,  22.4687816329404,  22.4705892332300 ,&  
            22.4724909243966,  22.4745033790201,  22.4766413819690 ,&  
            22.4789175475646,  22.4813431270825,  22.4839292152183 ,&  
            22.4866879089724,  22.4896332391652,  22.4927818177603 ,&  
            22.4961531945452,  22.4997700019047,  22.5036580631821 ,&  
            22.5078466299829,  22.5123687623642,  22.5172618466498 ,&  
            22.5225690254500,  22.5283444745393,  22.5346648097843 ,&  
            22.5416286813881,  22.5494483137879,  22.5584520587677 ,&  
            22.5690120220543,  22.5816332341286,  22.5970686893094 ,&  
            22.6165133803753,  22.6419864066045,  22.6771245216542 ,&  
            22.7304485274531,  22.8124251454444,  9.41332494279381 ,&  
            9.45698473247615,  9.54567091554276,  10.1050233306145 ,&  
            10.2800000000000    /) 


       !__________________________Time Settings__________________________________________________________________________________
       time       = 0.0_wp
       time_out   = 3600._wp
       dt         = 10.0_wp
       time_total = time_out*24._wp*7._wp*10._wp
       
    ELSE IF (testcase==9) THEN !#########################################################################
       !Simulation setup mostly used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 

       fl_q_bottom         = 10._wp           !Heat flux from below
       alpha_flux_instable = 22.0_wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 15._wp    !Chosen at random 
       tank_depth          = 0.8_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 10
       N_top    = 3
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       grav_heat_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -15._wp
       T_top       = -10._wp
       T_bottom    = -0.07_wp !-1.8_wp!0.0
       S_bu_bottom = 34.6_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.005_wp
       time_out   = 3600._wp*2._wp
       time_total = time_out*12._wp*6._wp
       dt         = 10._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
       
    ELSE IF (testcase==33) THEN !#########################################################################
       !Simulation setup mostly used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 
       
       ! Experiment 3, Freshwater, Niels Fuchs

       fl_q_bottom         = 10._wp           !Heat flux from below
       alpha_flux_instable = 22.0_wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 15._wp    !Chosen at random 
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 10
       N_top    = 3
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       grav_heat_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -15._wp
       T_top       = -10._wp
       T_bottom    = 0.5_wp !-1.8_wp!0.0
       S_bu_bottom = 0.13_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.005_wp
       time_out   = 60._wp*5._wp
       time_total = time_out*12._wp*6._wp
       dt         = 10._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF
       
    ELSE IF (testcase==34) THEN !#########################################################################
       !Simulation setup mostly used to imitate ice growth and melt in cooling
       !chambers.
       !Boundflux_flag 3 is used, under which T2m is the air temp in the
       !cooling chamber, and the surface flux is proportional is to the ice air
       !temperature difference.
       !When reproducing cooling chamber experiments the alpha flux parameters
       !need to be tuned, and a module in mo_testcase_specifics is needed to set
       !T2m over time.
       !fl_q_bottom for such setups can be very hard to define if the heat input
       !is not carefully measured and pumps or similar devices are involved. 
       
       ! Experiment 3, Freshwater, Niels Fuchs

       fl_q_bottom         = 10._wp           !Heat flux from below
       alpha_flux_instable = 22.0_wp   !Was chosen to produce 9 cm of ice in 42 hours for T2m = -18
       alpha_flux_stable   = 15._wp    !Chosen at random 
       tank_depth          = 0.94_wp

       !*************************************************************************************************************************
       !Layer settings and allocation
       !*************************************************************************************************************************
       Nlayer   = 100
       N_bottom = 10
       N_top    = 3
       N_active = 1
       N_middle = Nlayer-N_top-N_bottom
       CALL sub_allocate(Nlayer)

       !*************************************************************************************************************************
       !Flags
       !*************************************************************************************************************************
       tank_flag      = 2
       boundflux_flag = 3
       grav_heat_flag = 1


       !*************************************************************************************************************************
       !Salinity and temperature settings
       !*************************************************************************************************************************
       T2m         = -15._wp
       T_top       = -10._wp
       T_bottom    = 0.5_wp !-1.8_wp!0.0
       S_bu_bottom = 34.9_wp


       !*************************************************************************************************************************
       !Initial layer thickness, timestep, output time and simulation time
       !*************************************************************************************************************************
       thick_0    = 0.005_wp
       time_out   = 60._wp*10._wp
       time_total = 86400._wp * 10._wp
       dt         = 10._wp

       !*************************************************************************************************************************
       !total mass and salinity are calculated with the water depth in the tank to account for changing salinity
       !beneath the ice.
       !*************************************************************************************************************************
       m_total       = rho_l            *tank_depth             !Multiply with tank depth in meters
       S_total       = rho_l*S_bu_bottom*tank_depth

       !*************************************************************************************************************************
       !setting the initial values of the top and only layer
       !*************************************************************************************************************************
       thick(1) = thick_0
       m        = thick*rho_l
       S_abs    = S_bu_bottom*m
       H_abs    = m*T_bottom
       
       !*************************************************************************************************************************
       !BGC settings, only relevant if bgc_flag is set to 2 
       !*************************************************************************************************************************

       bgc_flag    = 1

       IF (bgc_flag==2) THEN
          N_bgc = 2
          CALL sub_allocate_bgc(Nlayer,N_bgc)

          
          bgc_bottom(1) = 385._wp
          bgc_bottom(2) = 385._wp

          !total values to calculate the bulk values beneath the ice in the tank
          !*************************************************************************************************************************
          !total bgc amount are calculated with the water depth in the tank to
          !account for changing concentrations in the tank
          !*************************************************************************************************************************
          bgc_total(:)  = bgc_bottom(:)*rho_l*tank_depth 

          bgc_abs(1,:)  = bgc_bottom(:) * m(1) 
          bgc_bu (1,:)  = bgc_bottom(:) 
          bgc_br (1,:)  = bgc_bottom(:) 


       END IF

    ELSE

       PRINT*,'selected testcase does not exist'
       STOP 4321

    END IF


    !Varius default settings
    T          = T_bottom
    S_bu       = S_bu_bottom
    psi_s      = 0._wp
    ray        = 0.0_wp
    phi        = 0.0_wp
    psi_l      = 1.0_wp
    fl_rad     = 0.0_wp
    bulk_salin = 0._wp
    grav_salt  = 0._wp
    grav_drain = 0._wp

    thick_min = thick_0/2._wp
    melt_thick_output(:) = 0._wp !< Niels, 2017

    !Calculates the number of timesteps for given time_total and dt as well as
    !timesteps between output 
    i_time     = INT(time_total/dt)
    i_time_out = INT(time_out  /dt)
    n_time_out = 0



    grav_drain = 0._wp
    grav_salt  = 0._wp
    grav_temp  = 0._wp
    melt_thick = 0
    thickness  = 0._wp
    bulk_salin = SUM(S_abs(1:N_active))/SUM(m(1:N_active))

    IF(N_top<3) THEN
       PRINT*,'Problem occurs when N_top smaller then 3, so just change it to 3 or more'
       STOP 666
    END IF


    IF(bgc_flag==2 .AND. (grav_flag+flush_flag+flood_flag).NE.9) THEN
       PRINT*,'WARNING: Biogeochemistry is on, but none or not all complex brine parametrizations & 
            & are  activated. Make sure this is your intent'
    END IF

    PRINT*,'Initialization of testcase complete, testcase:',testcase


  END SUBROUTINE init

  !>
  !! Allocates Arrays.
  !!
  !! For a given number of layers Nlayers all arrays are allocated
  SUBROUTINE sub_allocate (Nlayer,length_input_lab)

    INTEGER, INTENT(in)    :: Nlayer  !<  number of layers
    INTEGER, INTENT(in), OPTIONAL    :: length_input_lab  !< Niels, 2017 add:  dimension of input arrays

    !allocated Nlayer

    ALLOCATE(H(Nlayer),H_abs(Nlayer))
    ALLOCATE(Q(Nlayer))
    ALLOCATE(T(Nlayer))
    ALLOCATE(S_abs(Nlayer),S_bu(Nlayer),S_br(Nlayer))
    ALLOCATE(thick(Nlayer))
    ALLOCATE(m(Nlayer))
    ALLOCATE(V_s(Nlayer),V_l(Nlayer),V_g(Nlayer))
    ALLOCATE(V_ex(Nlayer))
    ALLOCATE(phi(Nlayer))
    ALLOCATE(perm(Nlayer))
    ALLOCATE(flush_v(Nlayer))   !< Niels, 2017
    ALLOCATE(flush_h(Nlayer))   !< Niels, 2017
    ALLOCATE(flush_v_old(Nlayer))   !< Niels, 2017
    ALLOCATE(flush_h_old(Nlayer))   !< Niels, 2017
    ALLOCATE(psi_s(Nlayer),psi_l(Nlayer),psi_g(Nlayer))
    ALLOCATE(fl_rad(Nlayer))
    IF (present(length_input_lab)) THEN
       ALLOCATE(Tinput(length_input_lab))   !< Niels, 2017
       ALLOCATE(precipinput(length_input_lab))  !< Niels, 2017
       ALLOCATE(ocean_T_input(length_input_lab))    !< Niels, 2017
       ALLOCATE(ocean_flux_input(length_input_lab)) !< Niels, 2017
       ALLOCATE(styropor_input(length_input_lab))   !< Niels, 2017 
       ALLOCATE(Ttop_input(length_input_lab))   !< Niels, 2017
    END IF
  
    !allocated Nlayer+1

    ALLOCATE(fl_Q(Nlayer+1))
    ALLOCATE(fl_m(Nlayer+1))
    ALLOCATE(fl_s(Nlayer+1))


    !Allocate Nlayer-1
    ALLOCATE(ray(Nlayer-1))

       m     = 0._wp
       S_abs = 0._wp
       H_abs = 0._wp
       thick = 0._wp
       
       flush_v(:) = 0._wp   !< Niels, 2017
       flush_h(:) = 0._wp   !< Niels, 2017

  END SUBROUTINE sub_allocate
  
  !>
  !! Allocates BGC Arrays.
  !!
  SUBROUTINE sub_allocate_bgc (Nlayer,N_bgc)

    INTEGER, INTENT(in)    :: Nlayer,N_bgc  

    !Brine flux matrix
    ALLOCATE(fl_brine_bgc(Nlayer+1,Nlayer+1))
    !Chemical matrices
    ALLOCATE(bgc_abs(Nlayer,N_bgc),bgc_bu(Nlayer,N_bgc),bgc_br(Nlayer,N_bgc))
    !Bottom values
    ALLOCATE(bgc_bottom(N_bgc),bgc_total(N_bgc))
   
   bgc_abs      = 0.0_wp
   fl_brine_bgc = 0.0_wp
   bgc_bu       = 0.0_wp
   bgc_br       = 0.0_wp
   bgc_bottom   = 0.0_wp 
   bgc_total    = 0.0_wp 

  END SUBROUTINE sub_allocate_bgc
  
  
  !>
  !! Deallocates Arrays.
  !!
  SUBROUTINE sub_deallocate 
    DEALLOCATE(H,H_abs)
    DEALLOCATE(Q)
    DEALLOCATE(T)
    DEALLOCATE(S_abs,S_bu,S_br)
    DEALLOCATE(thick)
    DEALLOCATE(m)
    DEALLOCATE(V_s,V_l,V_g)
    DEALLOCATE(V_ex)
    DEALLOCATE(phi)
    DEALLOCATE(perm)
    DEALLOCATE(psi_s,psi_l,psi_g)
    DEALLOCATE(fl_rad)

    !allocated Nlayer+1
    DEALLOCATE(fl_Q)
    DEALLOCATE(fl_m)
    DEALLOCATE(fl_s)

    !Allocate Nlayer-1
    DEALLOCATE(ray)

  END SUBROUTINE sub_deallocate

END MODULE mo_init

