!>
!! Computes all heat fluxes 
!!
!! Everything related to heat fluxes happens in sub_heat_fluxes, which is why it is a very crucial part of SAMSIM.
!! 
!!
!! @author Philipp Griewank
!!
!!  COPYRIGHT
!!        Copyright (c) 2014 Max-Planck-Institut fuer Meteorologie, Hamburg,
!!        Germany
!!
!!        Copying and distribution of this file, with or without modification,
!!        are permitted in any medium without royalty provided the copyright
!!        notice and this notice are preserved.  This file is offered as-is,
!!        without any warranty.
!!
!!
!!
!! @par Revision History
!! Copy and pasted into existence by Philipp Griewank (2014-04-02)
MODULE mo_heat_fluxes

  USE mo_parameters
  USE mo_thermo_functions
  USE mo_functions
  IMPLICIT NONE

  PUBLIC
CONTAINS
  !>
  !! Computes surface temperature and heatfluxes. 
  !!
  !! Major subroutine, calculates all atmospheric energy fluxes and applies both atmospheric and oceanic fluxes.
  !! Is one of the only subroutines to directly use mo_data because so many variables are needed. 
  !! 
  !! There are three different ways to calculate atmospheric heat fluxes implemented which are defined using boundflux_flag. 
  !!
  !! - Boundflux_flag: 1 imitates top cooling plate by setting a fixed surface temperature, heat flux is derived from the T gradient from the surface to the top layer
  !! - Boundflux_flag: 2 balances incoming and outgoing radiation to determine the surface temperature, heat flux is then calculated as in boundflux_flag 1. Some of the ice penetrates into the ice as is absorbed according to Beer's law. Optical properties are defined by the parameters emissivity_ice, emissivity_snow, extinct, and penetr. 
  !! - Boundflux_flag: 3 assumes the atmospheric heat flux is proportional to the difference between the top layer temperature and the air temperature. 
  !!
  !! For 1 and 2 the surface temperature in turn determines the atmospheric heat flux into the snow or ice.
  !! Atmoflux_flag is important for boundflux_flag 2, as it determines which atmospheric fluxes are used. 
  !! - Atmoflux_flag: 1 Mean climatology fluxes of Notz are used (see sub_notz)
  !! - Atmoflux_flag: 2 Imported values are used, see sub_input for more info on reading in data. 
  !! - Atmoflux_flag: 3 Prescribed values are used (e.g. testcase 5).
  !! 
  !! Melting occurs when the surface T is above the melting temperature of the top layer
  !! - Boundflux_flag: 1 atmospheric flux is limited by the parameter max_flux_plate which represents the maximum heating capacity of the plate
  !! - Boundflux_flag: 2 the atmospheric heat flux is given by the difference between incoming and outgoing radiation 
  !! - Boundflux_flag: 3 works the same during melt and freezing, but a different proportionality parameter is used (alpha_flux_stable) because the air above the ice is assumed to be stably stratified.
  !!
  !! Boundflux_flag 1 and 3 are not made to work with snow. If you need snow you'll have to implement snow cover yourself.
  !! For a detailed look at what is happening see the source code.
  !!
  !! The snow layer is treated differently based on the snow thickness. 
  !! - If the snow layer is thinner than thick_min/100 it is simply ignored.
  !! - If the snow layer is thinner than thick_min but thicker than thick_min/100 the snow and top ice layer are assumed to have the same temperature and are coupled using snow_coupling.
  !! - If the snow layer is thicker than thick_min it is treated totally separately.
  !!
  !!
  !!
  !! @par Revision History
  !! First version by Philipp Griewank (2014-04-02)
  !! Second version by Niels Fuchs (2017-02-02)

  Subroutine  sub_heat_fluxes()
    USE mo_snow
    USE mo_data
    USE mo_thermo_functions
    USE mo_functions
    IMPLICIT NONE
    REAL(wp) ::T_old, k_1, T_b,emi, pen, temp, temp1, temp2


       IF (boundflux_flag==1) THEN
          T_freeze = func_T_freeze(S_abs(1)/m(1),salt_flag)
          IF (func_S_br(T_top)>S_abs(1)/m(1)) THEN
             CALL sub_fl_Q_0(psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),T_top,-1,fl_Q(1))
          else
             fl_Q(1)   = -max_flux_plate
          end if
       end if
          


       IF (boundflux_flag==2) THEN
          !Calculates surface temperature by balancing incoming fluxes with a linear approximated black body emission.
          !Simplifications: Linearisation of black body emission, snow layer is ignored when smaller then thick_min.
          !Default fluxes are defined by sub_fluxnotz
          albedo = func_albedo(thick_snow,T_snow,psi_l(1),thick_min,albedo_flag) 
          IF (atmoflux_flag == 1) THEN
             CALL sub_notzflux(time+86400._wp*180._wp,fl_sw,fl_rest)
          else if (atmoflux_flag == 2) then
             IF (time==time_input(time_counter)) THEN
                fl_sw   = fl_sw_input(time_counter)
                fl_lw   = fl_lw_input(time_counter)
             ELSE
                temp=(time-time_input(time_counter-1))/(time_input(time_counter)-time_input(time_counter-1))
                fl_sw  = (1._wp-temp)*fl_sw_input(time_counter-1)   +temp*fl_sw_input(time_counter)
                fl_lw  = (1._wp-temp)*fl_lw_input(time_counter-1)   +temp*fl_lw_input(time_counter)
             END IF
             !Sensible and latent heat fluxes are ignored for now
             fl_sen  = 0._wp
             fl_lat  = 0._wp

             fl_rest = fl_lw+fl_sen+fl_lat
          END IF
          

          IF (thick_snow<thick_min) THEN
             T_old = T(1)
          ELSE
             T_old = T_snow
          END IF
          IF (thick_snow<thick_min) THEN
             !k_1 = psi_s(1)*k_s+psi_l(1)*k_l
             !T_b = T(1)!_bottom
             emi = emissivity_ice
             pen = penetr
          ELSE
             !k_1 = 0.3
             !T_b = T_snow
             emi = emissivity_snow
             pen = 0.0_wp
          END IF
          !T_b   = T_b   +zeroK !< Niels, 2017 commented, since it's unused afterwards
          T_old = T_old + zeroK

          !T_top is calculated by balancing incoming and outgoing radiation, and
          !is iterated once for some additional precision
          temp1 = (1._wp-albedo)*(1._wp-pen)*fl_sw+fl_rest
          temp1 = temp1+emi*3._wp*sigma*T_old**4._wp
          temp1 = temp1/(emi*4._wp*sigma*T_old**3._wp)
          !temp2 = (1._wp-albedo)*fl_sw+fl_rest+3._wp*sigma*T_old**4._wp-temp1*4.*sigma*T_old**3._wp !< Niels, 2017 commented, since it's unused afterwards
          temp1 = temp1-zeroK
          
          T_old = temp1 + zeroK 
          temp1 = (1._wp-albedo)*(1._wp-pen)*fl_sw+fl_rest
          temp1 = temp1+emi*3._wp*sigma*T_old**4._wp
          temp1 = temp1/(emi*4._wp*sigma*T_old**3._wp)
          !temp2 = (1._wp-albedo)*fl_sw+fl_rest+3._wp*sigma*T_old**4._wp-temp1*4.*sigma*T_old**3._wp !< Niels, 2017 add: commented, since it isnt used afterwards
          temp1 = temp1-zeroK

          T_top = temp1

          !Absorption of penetrating sw radiation with a depth independent extinction coefficient.
          temp2 = pen*(1._wp-albedo)*fl_sw
          DO k = 1,N_active
             fl_rad(k) = temp2-temp2*EXP(-extinc*thick(k))
             temp2     = temp2      *EXP(-extinc*thick(k))
          END DO
          
          
          IF (thick_snow>=thick_min/100._wp) THEN
             T_freeze = 0._wp
          ELSE 
             T_freeze = func_T_freeze(S_abs(1)/m(1),salt_flag)
          END IF
          
           


          IF (T_top>T_freeze.AND. N_active>1) THEN 
          !If T_top greater then T_freeze and ice is present (n_active>1) T_top is set to T_freeze and the top flux 
          !has to be calculated from the radiation imbalance.
             temp1 = emi*sigma*(T_freeze+zeroK)**4._wp-(1._wp-albedo)*(1._wp-pen)*fl_sw-fl_rest
             IF (thick_snow>=thick_min) THEN 
                fl_Q_snow = temp1
                CALL sub_fl_Q_snow(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),fl_Q(1))
             else IF (thick_snow>=thick_min/100._wp) THEN
                fl_Q_snow = temp1
                fl_Q(1)   = 0._wp
             ELSE 
                fl_Q(1)   = temp1
             END IF
             T_top = T_freeze
          
          ELSE
!           !If T_top below freezing
                    !Boundary fluxes
             IF (thick_snow>=thick_min) THEN 
                CALL sub_fl_Q_snow(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),fl_Q(1))
                CALL sub_fl_Q_0_snow(m_snow,thick_snow,T_snow,T_top,fl_Q_snow)
             ELSE IF (thick_snow>thick_min/100._wp .AND. thick_snow<thick_min) THEN
                fl_Q(1) = 0.0_wp
                CALL sub_fl_Q_0_snow_thin(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),T_top,fl_Q_snow)
             ELSE
                CALL sub_fl_Q_0(psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),T_top,-1,fl_Q(1))
             END IF
          end if
       END IF

       
       
       !Meant to reproduce lab freezing
       !Relies on two parameters, alpha_flux_stable < alpha_flux_instable, which
       !vary from experiment to experiment, set values as needed in mo_init.f90
       IF (boundflux_flag==3) THEN
       
          !< Niels, 2017 add: old lab setup or snow free
          
          IF (lab_snow_flag==0 .OR. thick_snow<=thick_min/100._wp) THEN
          
             T_freeze = MIN(func_T_freeze(S_abs(N_active)/m(N_active),salt_flag),0._wp)
             T_top    = T(1)
             fl_Q(1)  = alpha_flux_instable*(T_top-T2m)

             IF(fl_Q(1)<0._wp) THEN
                T_top   = MAX(T_freeze,T(1))
                fl_Q(1) = alpha_flux_stable*(T_top-T2m)  
             END IF
             !< Niels, 2017 add: if styropor simulates snow heat fluxes in the lab
             IF (thick_snow == 0._wp .AND. lab_snow_flag==1 .AND. styropor_flag == 1) THEN
                CALL sub_fl_Q_styropor(k_styropor,fl_Q(1))
             END IF
                          
                                       
          !< Niels, 2017 add:  newly implemented snowphysics for the lab based on boundflux_flag 2
          
          ELSE IF (lab_snow_flag==1) THEN
          
             T_freeze = func_T_freeze(S_abs_snow/m_snow,salt_flag)  !< Niels, 2017 add: T_freeze of snow cover

             T_top = T_snow !< Niels, 2017 add: T_top is T_snow
             
             temp1 = alpha_flux_instable*(T_top-T2m) !< cooling
             
             IF (temp1 >= 0.0_wp) THEN
                IF (thick_snow>=thick_min) THEN 
                   fl_Q_snow = temp1
                   CALL sub_fl_Q_snow(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),fl_Q(1))
                else IF (thick_snow>=thick_min/100._wp) THEN
                   CALL sub_fl_Q_0_snow_thin(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),&
                   &(T2m+T_top)/2._wp,fl_Q_snow) !< Niels, 2017 add: (T2m+T_top)/2. not an optimal solution, but better than T_top here
                   fl_Q(1)   = 0.0_wp
                END IF

             ELSE
             
                temp1 = alpha_flux_stable*(T_top-T2m) !< warming
               
                IF (thick_snow>=thick_min) THEN 
                   fl_Q_snow = temp1
                   CALL sub_fl_Q_snow(m_snow,thick_snow,T_snow,psi_s(1),psi_l(1),psi_g(1),thick(1),T(1),fl_Q(1))
                else IF (thick_snow>=thick_min/100._wp) THEN
                   fl_Q_snow = temp1
                   fl_Q(1)   = 0._wp
                END IF
                
             END IF
             

          END IF 
       END IF


       !Heatflux from below
       fl_q(N_Active+1) = fl_q_bottom


       !Built in energy check, makes sure that energy is conserved when
       !heatfluxes are applied, temp1 is the total enthalpy plus the fluxes
       !and temp2 is the total energy after heat fluxes are applied
       !Temp1 should equal temp2 down to machine precision
       temp1 = sum(H_abs)+H_abs_snow

       !Intermediate fluxes
       DO k = 2,N_active
          CALL sub_fl_Q(psi_s(k-1),psi_l(k-1),psi_g(k-1),thick(k-1),T(k-1),psi_s(k),psi_l(k),psi_g(k),thick(k),T(k),fl_Q(k))
       END DO

       !New Enthalpies due to heat flux (explicit)
       DO k = 1,N_active
          H_abs(k) = H_abs(k)+(fl_Q(k+1)-fl_Q(k))*dt
       END DO

       !New Enthalpies due to radiation absorption
       DO k = 1,N_active
          H_abs(k) = H_abs(k)+fl_rad(N_active)*dt
          temp1 = temp1+fl_rad(N_active)*dt
       END DO


       !##########################################################################################
       !Snow treatment
       !##########################################################################################
       IF(thick_snow>=thick_min/100._wp .AND. thick_snow<thick_min)  THEN
          H_abs_snow = H_abs_snow-fl_Q_snow*dt
          CALL snow_coupling (H_abs_snow,phi_s,T_snow,H_abs(1),H(1),phi(1),T(1),m_snow,S_abs_snow,m(1),S_bu(1))
          
          temp1=temp1+fl_q_bottom*dt-fl_Q_snow*dt
       ELSE IF (thick_snow>=thick_min) THEN
          H_abs_snow = H_abs_snow+(fl_Q(1)-fl_q_snow)*dt
          
          temp1=temp1+fl_q_bottom*dt-fl_Q_snow*dt
       ELSE
          
          temp1=temp1+fl_q_bottom*dt-fl_Q(1)*dt
       END IF

       temp2 = sum(H_abs)+H_abs_snow

       IF (abs((temp1-temp2)/dt) > 0.00001_wp) then
          print*,'energie discrepencies during heat flux calculations higher than 0.00001 Joules per second',(temp1-temp2)/dt
          STOP 431
       end if

  END subroutine sub_heat_fluxes

end module mo_heat_fluxes
