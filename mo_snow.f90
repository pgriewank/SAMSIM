
!>
!!  Module contains all things directly related to snow. 
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
!  
!!
!! @par Revision History
!! Provided for by Philipp Griewank 2010-12-13
!!
!!
!!
MODULE mo_snow

  USE mo_parameters
  USE mo_thermo_functions

  IMPLICIT NONE

  PUBLIC :: snow_thermo
  PUBLIC :: snow_precip_0
  PUBLIC :: snow_precip
  PUBLIC :: sub_fl_q_snow
  PUBLIC :: sub_fl_q_0_snow
  PUBLIC :: sub_fl_Q_0_snow_thin
  PUBLIC :: snow_coupling
  PUBLIC :: func_k_snow
  PUBLIC :: snow_thermo_meltwater   !< Niels, 2017
  PRIVATE


CONTAINS 
  
  !>
  !! Subroutine to couple a thin snow layer to the upper ice layer.
  !!
  !! Subroutine is activated when thick_snow<thick_min.
  !! The enthalpies of the two layers are adjusted until both layers have the same temperatures.
  !! The following approach is used. 
  !! 1. The enthalpies are adjusted so T_snow=0,  and phi_s=1.
  !! 2. The temperatures are calculated.
  !! 3. If the ice temperature is greater 0 the balanced enthalpies are calculated directly.
  !!    ELSE they are calculated iteratively.
  !!
  !! @par Revision History
  !! Written by Philipp Griewank,  IMPRS (2011-01-20)
  !!
  SUBROUTINE snow_coupling(H_abs_snow, phi_s, T_snow, H_abs, H, phi, T, m_snow, S_abs_snow, m, S_bu)
    REAL(wp), INTENT(inout) :: H_abs_snow, phi_s, T_snow
    REAL(wp), INTENT(inout) :: H_abs, H, phi, T
    REAL(wp), INTENT(in)    :: m_snow, S_abs_snow
    REAL(wp), INTENT(in)    :: m, S_bu

    INTEGER                 :: jj

    H_abs      = H_abs+m_snow*latent_heat+H_abs_snow
    H_abs_snow =-m_snow*latent_heat
    H          = H_abs/m

    CALL getT(H_abs_snow/m_snow, S_abs_snow/m_snow,  T_snow, T_snow, phi_s, 5701) 
    CALL getT(H, S_bu, T, T, phi, 5702)

    IF (T>0 .AND. H_abs <= -H_abs_snow) THEN
       H_abs_snow = H_abs_snow+H_abs
       H_abs      = 0._wp
       CALL getT(H_abs_snow/m_snow, S_abs_snow/m_snow, T_snow, T_snow, phi_s, 5701) 
       CALL getT(H, S_bu, T, T, phi, 5702)
    ELSE IF (T>0. .AND. H_abs>-H_abs_snow) THEN
       H_abs      = (H_abs+H_abs_snow)*m/m_snow/(1._wp+m/m_snow)
       H_abs_snow = H_abs*m_snow/m
       CALL getT(H_abs_snow/m_snow, S_abs_snow/m_snow, T_snow, T_snow, phi_s, 5701) 
       CALL getT(H, S_bu, T, T, phi, 5702)
    ELSE
       jj = 0
       DO WHILE (ABS(T-T_snow)>0.1 .AND. jj<201)
          H_abs_snow = H_abs_snow -SIGN(MAX(ABS(T_snow-(T_snow+T)/2._wp), 0.1_wp), (T_snow-(T_snow+T)/2._wp)) *c_s*m_snow
          H_abs      = H_abs +SIGN(MAX(ABS(T_snow-(T_snow+T)/2._wp), 0.1_wp), (T_snow-(T_snow+T)/2._wp)) *c_s*m_snow
          jj = jj+1
          H  = H_abs/m
          IF(jj>100) THEN
             PRINT*, 'T, T_snow', T, T_snow, jj
          END IF
          CALL getT(H_abs_snow/m_snow, S_abs_snow/m_snow, T_snow, T_snow, phi_s, 5701) 
          CALL getT(H, S_bu, T, T, phi, 5702)
       END DO
       IF(jj>200 .AND. ABS(T-T_snow)>1._wp) THEN
          PRINT*, 'coupling of thin snow and top layer crashed'
          STOP 16
       END IF
    END IF
  END SUBROUTINE snow_coupling
  
  
  
  !>
  !! Subroutine for calculating precipitation on an existing snow cover.
  !!
  !! Can optionally deal with separate solid and liquid precipitation or a single liquid input.
  !! The 2 meter temperature  determines the temperature of the precipitation.
  !! In case of single input the 2 meter temperature determines if snow or rain falls.
  !! Snow makes the thickness grow according to the density of new snow(rho_snow),  while rain falls into the snow without increasing snow depth.
  !! It is necessary to calculate the new psi_s_snow to ensure proper melting in snow_thermo.
  !!
  !!
  !! @par Revision History
  !! Sired by Philipp Griewank,  IMPRS (2010-12-14)
  !!
  SUBROUTINE snow_precip (m_snow, H_abs_snow, thick_snow, psi_s_snow, dt, liquid_precip_in, T2m, solid_precip_in)
    REAL(wp), INTENT(inout)        :: m_snow, H_abs_snow, thick_snow
    REAL(wp), INTENT(inout)        :: psi_s_snow 
    REAL(wp), INTENT(in)           :: liquid_precip_in, T2m, dt
    REAL(wp), INTENT(in), OPTIONAL :: solid_precip_in
    REAL(wp)                       :: solid_precip, liquid_precip, d_thick

    IF(PRESENT(solid_precip_in)) THEN
       solid_precip  = solid_precip_in
       liquid_precip = liquid_precip_in
    ELSE
       IF (T2m>0.0_wp) THEN
          solid_precip  = 0._wp
          liquid_precip = liquid_precip_in
       ELSE
          solid_precip  = liquid_precip_in
          liquid_precip = 0._wp
       END IF
    END IF
    d_thick    =             dt*solid_precip*rho_l/rho_snow
    m_snow     = m_snow     +dt*rho_l*(liquid_precip+solid_precip)
    thick_snow = thick_snow +d_thick
    H_abs_snow = H_abs_snow +dt*T2m*liquid_precip*rho_l*c_l
    H_abs_snow = H_abs_snow +dt*-1._wp*solid_precip*rho_l*c_s
    H_abs_snow = H_abs_snow -dt*solid_precip*rho_l*latent_heat

  END SUBROUTINE snow_precip
  
  
  !>
  !! Subroutine for calculating precipitation into the ocean.
  !!
  !! Can optionally deal with separate solid and liquid precipitation or a single liquid input.
  !! The 2 meter temperature  determines the temperature of the precipitation.
  !! In case of single input the 2 meter temperature determines if snow or rain falls.
  !! It is important,  that the mass, energy and salt leaving the upper layer must be outputted.
  !! This is not the case. Temp!
  !!
  !!
  !!
  !! @par Revision History
  !! Copy and Pasted by Philipp Griewank,  IMPRS (2011-01-10)
  !!
  SUBROUTINE snow_precip_0 (H_abs, S_abs, m, T, dt, liquid_precip_in, T2m, solid_precip_in)
    REAL(wp), INTENT(inout)        :: H_abs, S_abs !>Values of water layer
    REAL(wp), INTENT(in)           :: m, T !>Value of water layer
    REAL(wp), INTENT(in)           :: liquid_precip_in, T2m, dt
    REAL(wp), INTENT(in), OPTIONAL :: solid_precip_in
    REAL(wp)                       :: solid_precip, liquid_precip, d_thick

    IF(PRESENT(solid_precip_in)) THEN
       solid_precip  = solid_precip_in
       liquid_precip = liquid_precip_in
    ELSE
       IF (T2m>0.0_wp) THEN
          solid_precip  = 0._wp
          liquid_precip = liquid_precip_in
       ELSE
          solid_precip  = liquid_precip_in
          liquid_precip = 0._wp
       END IF
    END IF

    H_abs = H_abs+(liquid_precip+solid_precip)*(T2m-T)*dt 
    H_abs = H_abs-solid_precip*latent_heat*dt
 
    S_abs = S_abs-(liquid_precip+solid_precip)*S_abs/m*dt 

  END SUBROUTINE snow_precip_0

  !>
  !! Subroutine for calculating snow thermodynamics
  !!
  !! Behaves similar to mushy layer sea ice. 
  !! Important differences are:
  !! 1. no expulsion,  thick_snow is raised if the volume expands.
  !! 2. The liquid fraction is limited.
  !! 3. When the liquid fraction exceeds it's limit the thickness of the snow layer is reduced.
  !! This is done as follows:
  !! Only applies if the fluid fraction is above the irreducible water content as defined in Coleuo-Lasaffre 98.
  !! thick_snow=thick_snow*(1._wp-(psi_s_old-psi_s_snow)/psi_s_old)
  !! Warning: the formula for liquid water content in Coleuo-Lasaffre contains 2 typos
  !! When the water exceeds the limit water runs down to the bottom of the snow layer.
  !! The saturated lower layer is added to the top ice layer.
  !!
  !! @par Revision History
  !! Fabricated by Philipp Griewank,  IMPRS (2010-12-14)
  !! Major redo,  water saturated bottom snow added to top ice layer by Philipp Griewank (2010-12-14) 
  SUBROUTINE snow_thermo (psi_l_snow, psi_s_snow, psi_g_snow, thick_snow, S_abs_snow, H_abs_snow, m_snow, T_snow, m, thick, H_abs)
    REAL(wp), INTENT(inout) :: psi_l_snow, psi_s_snow, psi_g_snow, T_snow
    REAL(wp), INTENT(inout) :: S_abs_snow, H_abs_snow, m_snow
    REAL(wp), INTENT(inout) :: thick_snow
    REAL(wp), INTENT(inout) :: m, thick, H_abs                               !<Top ice layer variables
    REAL(wp)                :: sat_snow                                      !<height of saturated layer in snow [m]
    REAL(wp)                :: psi_s_old, phi_snow, T_in, H_snow, S_bu_snow
    REAL(wp)                :: max_lwc                                       !<max liquid water content,  in mass percent ,  see Coleuo-Lasaffre 98.
    REAL(wp)                :: max_lwc_v                                     !<max liquid water content,  in volume percent

    H_snow    = H_abs_snow/m_snow
    S_bu_snow = S_abs_snow/m_snow
    psi_s_old = psi_s_snow


    T_in = T_snow
    CALL getT(H_snow, S_bu_snow, T_in, T_snow, phi_snow, 5700)

    psi_s_snow = m_snow*phi_snow/rho_s/thick_snow
    psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
    IF (psi_s_snow+psi_l_snow>1._wp) THEN
       thick_snow = m_snow*(phi_snow/rho_s+(1._wp-phi_snow)/rho_l)
       psi_s_snow = m_snow*phi_snow/rho_s/thick_snow
       psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
       IF (ABS(psi_s_snow+psi_l_snow -1._wp)>0.0000001_wp) THEN
          PRINT*, 'fault in recalculating height in snow'
          PRINT*, psi_s_snow, psi_l_snow
          STOP 345
       END IF
    END IF


    psi_g_snow = 1._wp-psi_s_snow-psi_l_snow
    IF (psi_s_snow>0._wp) THEN
       max_lwc = 0.057_wp*(1._wp-psi_s_snow)/(psi_s_snow) + 0.017_wp !Coleuo-Lasaffre have 2 serious typos in the  equation.
    ELSE
       max_lwc = 0.0_wp
    END IF

    IF (psi_s_old>psi_s_snow .AND. psi_s_snow>0.0_wp) THEN
       IF ((1._wp-phi_snow)>max_lwc) THEN
          thick_snow = thick_snow*(1._wp-(psi_s_old-psi_s_snow)/psi_s_old)
       END IF
       IF (thick_snow<(phi_snow*m_snow/rho_s+(1._wp-phi_snow)*m_snow/rho_l)) THEN
          thick_snow = (phi_snow*m_snow/rho_s+(1._wp-phi_snow)*m_snow/rho_l)
       END IF
       psi_s_snow = m_snow*phi_snow/rho_s/thick_snow
       psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
       psi_g_snow = 1._wp-psi_s_snow-psi_l_snow
       psi_g_snow = ABS(psi_g_snow)
    ELSE IF (psi_s_snow<0.000001_wp) THEN
       thick_snow = m_snow/rho_l
       psi_s_snow = 0.0_wp
       psi_g_snow = 0.0_wp
       psi_l_snow = 1.0_wp
    END IF


    !Calculate thickness of saturated layer at bottom of snow,  again see Coleuo-Lasaffre 98.
    !Saturated layer is then added to the top ice layer
    !The parameter gas_snow_ice2 controls the amount of gas which is entrapped in the saturated layer.
    !Added some security because I had problems with very slightly negative
    !psi_g.
    !If psi_g<0 then the snow layer instantly completely saturated.
    
    IF ((1._wp-phi_snow)>max_lwc .and. psi_g_snow>0.0_wp) THEN
       max_lwc_v = max_lwc*m_snow/(rho_l*thick_snow)
       sat_snow  = thick_snow*(psi_l_snow-max_lwc_v)
       sat_snow  = sat_snow/(1._wp-psi_s_snow-max_lwc_v-MIN(gas_snow_ice2, psi_g_snow))

       !Ok,  LAYER transformation!
       IF(sat_snow<0.0) THEN
          PRINT*, 'wtf', sat_snow, phi_snow, psi_l_snow, max_lwc_v
       END IF
       thick_snow = thick_snow -sat_snow
       thick      = thick +sat_snow
       m_snow     = m_snow -sat_snow*(psi_s_snow*rho_s+(1._wp-psi_s_snow-gas_snow_ice2)*rho_l)
       m          = m +sat_snow*(psi_s_snow*rho_s+(1._wp-psi_s_snow-gas_snow_ice2)*rho_l)

       H_abs_snow = H_abs_snow -sat_snow*psi_s_snow*rho_s*c_s*T_snow
       H_abs      = H_abs +sat_snow*psi_s_snow*rho_s*c_s*T_snow
       H_abs_snow = H_abs_snow +sat_snow*psi_s_snow*rho_s*latent_heat
       H_abs      = H_abs -sat_snow*psi_s_snow*rho_s*latent_heat
       H_abs_snow = H_abs_snow -sat_snow*(1._wp-psi_s_snow)*rho_l*c_l*T_snow
       H_abs      = H_abs +sat_snow*(1._wp-psi_s_snow)*rho_l*c_l*T_snow
       
    Else IF ( psi_g_snow.le.0.0_wp) THEN
       sat_snow   = thick_snow
       
       H_abs      = H_abs+H_abs_snow
       m          = m+m_snow
       thick      = thick+thick_snow
       H_abs_snow = 0.0_wp
       m_snow     = 0.0_wp
       thick_snow = 0.0_wp
       psi_g_snow = 0.0_wp
       psi_s_snow = 0.0_wp
       psi_l_snow = 0.0_wp
    ELSE
       sat_snow = 0.0_wp
    END IF


   IF (psi_g_snow<0._wp) then
     print*,'negative snow in snow_thermo',psi_g_snow,psi_s_snow,psi_l_snow,thick_snow
     stop 09876
   end if

  END SUBROUTINE snow_thermo

  !> Niels, 2017 add: 
  !! Subroutine for calculating snow thermodynamics
  !! most of the physics are taken from snow_thermo()
  !! based on lab observations: parts of the snow meltwater percolate directly into the ice
  !! @par Revision History
  !! introduced by Niels Fuchs (2016-10-13)
  
  SUBROUTINE snow_thermo_meltwater (psi_l_snow, psi_s_snow, psi_g_snow, thick_snow, S_abs_snow, H_abs_snow, m_snow, T_snow, m, &
  &thick, H_abs, melt_thick_snow)
    REAL(wp), INTENT(inout) :: psi_l_snow, psi_s_snow, psi_g_snow, T_snow
    REAL(wp), INTENT(inout) :: S_abs_snow, H_abs_snow, m_snow
    REAL(wp), INTENT(inout) :: thick_snow, melt_thick_snow
    REAL(wp), INTENT(inout) :: m, thick, H_abs                               !<Top ice layer variables
    REAL(wp)                :: sat_snow                                      !<height of saturated layer in snow [m]
    REAL(wp)                :: psi_s_old, phi_snow, T_in, H_snow, S_bu_snow
    REAL(wp)                :: max_lwc                                       !<max liquid water content,  in mass percent ,  see Coleuo-Lasaffre 98.
    REAL(wp)                :: max_lwc_v                                     !<max liquid water content,  in volume percent
    REAL(wp)                :: psi_l_snow_slush, psi_l_snow_flush            !<volume fraction of excess water that goes into flushing resp. slush(-ing)

    H_snow    = H_abs_snow/m_snow
    S_bu_snow = S_abs_snow/m_snow
    psi_s_old = psi_s_snow


    T_in = T_snow
    CALL getT(H_snow, S_bu_snow, T_in, T_snow, phi_snow, 5700)	! new T & phi
    psi_s_snow = m_snow*phi_snow/rho_s/thick_snow  ! Calculate psi_s & psi_l from new phi value
    psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
    IF (psi_s_snow+psi_l_snow>1._wp) THEN
       thick_snow = m_snow*(phi_snow/rho_s+(1._wp-phi_snow)/rho_l)  ! recalculate thick_snow and reduce psi_s and psi_l so that the sum is 1
       psi_s_snow = m_snow*phi_snow/rho_s/thick_snow
       psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
       IF (ABS(psi_s_snow+psi_l_snow -1._wp)>0.0000001_wp) THEN
          PRINT*, 'fault in recalculating height in snow'
          PRINT*, psi_s_snow, psi_l_snow
          STOP 345
       END IF
    END IF


    psi_g_snow = 1._wp-psi_s_snow-psi_l_snow
    IF (psi_s_snow>0._wp) THEN
       max_lwc = 0.057_wp*(1._wp-psi_s_snow)/(psi_s_snow) + 0.017_wp !Coleuo-Lasaffre have 2 serious typos in the  equation.
    ELSE
       max_lwc = 0.0_wp
    END IF

    IF (psi_s_old>psi_s_snow .AND. psi_s_snow>0.0_wp) THEN
       IF ((1._wp-phi_snow)>max_lwc) THEN
          thick_snow = thick_snow*(1._wp-(psi_s_old-psi_s_snow)/psi_s_old)	! Sum A in 2015 paper
       END IF
       IF (thick_snow<(phi_snow*m_snow/rho_s+(1._wp-phi_snow)*m_snow/rho_l)) THEN
          thick_snow = (phi_snow*m_snow/rho_s+(1._wp-phi_snow)*m_snow/rho_l)
       END IF
       psi_s_snow = m_snow*phi_snow/rho_s/thick_snow
       psi_l_snow = m_snow*(1._wp-phi_snow)/rho_l/thick_snow
       psi_g_snow = 1._wp-psi_s_snow-psi_l_snow
       psi_g_snow = ABS(psi_g_snow)
    ELSE IF (psi_s_snow<0.000001_wp) THEN
       thick_snow = m_snow/rho_l
       psi_s_snow = 0.0_wp
       psi_g_snow = 0.0_wp
       psi_l_snow = 1.0_wp
    END IF


    !Calculate thickness of saturated layer at bottom of snow,  again see Coleuo-Lasaffre 98.
    !Saturated layer is then added to the top ice layer
    !The parameter gas_snow_ice2 controls the amount of gas which is entrapped in the saturated layer.
    !Added some security because I had problems with very slightly negative
    !psi_g.
    !If psi_g<0 then the snow layer instantly completely saturated.
    
    !IF ((1._wp-phi_snow)>max_lwc .and. psi_g_snow>0.0_wp) THEN
    IF ((1._wp-phi_snow)>max_lwc .and. psi_l_snow>0.0_wp .and. psi_g_snow>0.0_wp) THEN
       max_lwc_v = max_lwc*m_snow/(rho_l*thick_snow)
       
       !!! newly included separation of excess water into slush and flush source by Niels Fuchs 2016
       
       psi_l_snow_slush = (psi_l_snow-max_lwc_v)*(1._wp - k_snow_flush)
       psi_l_snow_flush = (psi_l_snow-max_lwc_v)*k_snow_flush
       
       melt_thick_snow = thick_snow * psi_l_snow_flush
	
       sat_snow  = thick_snow*(psi_l_snow_slush)	! equation 5 in 2015 paper, value B
       sat_snow  = sat_snow/(1._wp-psi_s_snow-max_lwc_v-MIN(gas_snow_ice2, psi_g_snow))
       
       
       !Ok,  LAYER transformation!
       IF(sat_snow<0.0) THEN
          PRINT*, 'wtf', sat_snow, phi_snow, psi_l_snow, max_lwc_v
       END IF
       

       thick_snow = thick_snow - sat_snow - melt_thick_snow
       thick      = thick + sat_snow
       m_snow     = m_snow -sat_snow*(psi_s_snow*rho_s+(1._wp-psi_s_snow-MIN(gas_snow_ice2, psi_g_snow))*rho_l) - &
       & melt_thick_snow*rho_l
       m          = m +sat_snow*(psi_s_snow*rho_s+(1._wp-psi_s_snow-MIN(gas_snow_ice2, psi_g_snow))*rho_l)

       H_abs_snow = H_abs_snow -sat_snow*psi_s_snow*rho_s*c_s*T_snow 
       H_abs      = H_abs +sat_snow*psi_s_snow*rho_s*c_s*T_snow
       H_abs_snow = H_abs_snow +sat_snow*psi_s_snow*rho_s*latent_heat
       H_abs      = H_abs -sat_snow*psi_s_snow*rho_s*latent_heat
       H_abs_snow = H_abs_snow - sat_snow*(1._wp-psi_s_snow-MIN(gas_snow_ice2, psi_g_snow))*rho_l*c_l*T_snow - &
       &melt_thick_snow*rho_l*c_l*T_snow
       H_abs      = H_abs +sat_snow*(1._wp-psi_s_snow-MIN(gas_snow_ice2, psi_g_snow))*rho_l*c_l*T_snow

    Else IF ( psi_g_snow.le.0.0_wp) THEN
       sat_snow   = thick_snow
       
       H_abs      = H_abs+H_abs_snow
       m          = m+m_snow
       thick      = thick+thick_snow
       H_abs_snow = 0.0_wp
       m_snow     = 0.0_wp
       thick_snow = 0.0_wp
       psi_g_snow = 0.0_wp
       psi_s_snow = 0.0_wp
       psi_l_snow = 0.0_wp
    ELSE
       sat_snow = 0.0_wp
    END IF


   IF (psi_g_snow<0._wp) then
     print*,'negative snow in snow_thermo',psi_g_snow,psi_s_snow,psi_l_snow,thick_snow
     stop 09876
   end if

  END SUBROUTINE snow_thermo_meltwater
  !>
  !! Determines conductive Heat flux for combined top ice and snow layer.
  !!
  !! When thick_snow<thick_min.
  !! 
  !!
  !!
  !! @par Revision History
  !! first version by Philipp Griewank (2011-01-19)
  !!
  SUBROUTINE sub_fl_Q_0_snow_thin(m_snow, thick_snow, T_snow, psi_s, psi_l, psi_g, thick,  T_bound, fl_Q_snow)


    REAL(wp),  INTENT(in)  :: m_snow
    REAL(wp),  INTENT(in)  :: psi_s, psi_l, psi_g
    REAL(wp),  INTENT(in)  :: thick_snow, thick
    REAL(wp),  INTENT(in)  :: T_snow,  T_bound

    REAL(wp),  INTENT(out) :: fl_Q_snow

    REAL(wp)               :: R         !< Thermal resistance
    REAL(wp)               :: k_snow, k !< Thermal conductivity

    k_snow = func_k_snow(m_snow, thick_snow)
    k = psi_s*k_s+psi_l*k_l+psi_g*0._wp
    k = thick_snow/(thick_snow+thick)*k_snow+thick/(thick_snow+thick)*k

    R = (thick_snow+thick)/(2._wp*k)

    fl_Q_snow = (T_snow-T_bound)/R

  END SUBROUTINE sub_fl_Q_0_snow_thin


  !>
  !! Determines conductive Heat flux between Snow and top ice layer.
  !!
  !! Standard approach.
  !!
  !! @par Revision History
  !! first version by Philipp Griewank (2010-12-15)
  !!
  SUBROUTINE sub_fl_Q_snow(m_snow, thick_snow, T_snow, psi_s_2, psi_l_2, psi_g_2, thick_2, T_2, fl_Q)


    REAL(wp),  INTENT(in)  :: m_snow
    REAL(wp),  INTENT(in)  :: psi_s_2, psi_l_2, psi_g_2
    REAL(wp),  INTENT(in)  :: thick_snow, thick_2
    REAL(wp),  INTENT(in)  :: T_snow, T_2
    REAL(wp),  INTENT(out) :: fl_Q

    REAL(wp) :: R           !< Thermal resistance
    REAL(wp) :: k_snow, k_2 !< Thermal conductivity


    k_snow = func_k_snow(m_snow, thick_snow)
    k_2    = psi_s_2*k_s+psi_l_2*k_l

    R      = thick_snow/(2._wp*k_snow)+thick_2/(2._wp*k_2)

    fl_Q   = (T_2-T_snow)/R

  END SUBROUTINE sub_fl_Q_snow

  !>
  !! Determines conductive Heat between snow layer and upper boundary layer.
  !! A limiting factor is added to increase stability of layers thinner then thick_min.
  !!
  !! @par Revision History
  !! first version by Philipp Griewank (2010-12-15)
  !! Artificial limitation introduced by Philipp Griewank (2011-01-17)
  SUBROUTINE sub_fl_Q_0_snow(m_snow, thick_snow, T_snow, T_bound, fl_Q)


    REAL(wp),  INTENT(in)  :: m_snow
    REAL(wp),  INTENT(in)  :: thick_snow
    REAL(wp),  INTENT(in)  :: T_snow, T_bound !< T_bound temperature of boundary layer
    REAL(wp),  INTENT(out) :: fl_Q

    REAL(wp)               :: R !< Thermal resistance
    REAL(wp)               :: k !< Thermal conductivity

    k    = func_k_snow(m_snow, thick_snow)

    R    = thick_snow/(2._wp*k)

    fl_Q = (T_snow-T_bound)/R

  END SUBROUTINE sub_fl_Q_0_snow



  !>
  !! Calculates the thermal conductivity of the snow layer as a function of the density.
  !!
  !! Based on the Sturm et al 1997 data fit for densities greater then 0.156 g/cm**3.
  !! Warning,  Sturm et al use g/cm**3,  I use kg/m**3
  !! Snow density probability functions can be included lated to raise the effective conductivity. 
  !! Warning!: added 0.15 to the thermal conductivity. 
  !! 
  !! @par Revision History
  !! Forged by Philipp Griewank (2010-12-13)
  !!
  FUNCTION func_k_snow(m_snow, thick_snow) RESULT (k_snow)
    IMPLICIT NONE

    REAL(wp),  PARAMETER     :: c0 = 0.138_wp
    REAL(wp),  PARAMETER     :: c1 = -1.01_wp/1000._wp
    REAL(wp),  PARAMETER     :: c2 = 3.233_wp/1000000._wp

    REAL(wp),  INTENT(in)    :: m_snow, thick_snow
    REAL(wp)                 :: k_snow 

    k_snow = c0 + c1*m_snow/thick_snow + c2*(m_snow/thick_snow)**2._wp
    k_snow = k_snow+0.15 !WARNING

  END FUNCTION func_k_snow


END MODULE mo_snow
