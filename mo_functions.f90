!>
!!  Module houses functions which have no home :(.
!!
!!  Created because I wanted to calculate the freeboard separately and didn't know where to put it.
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
!! Ribbon cut by Philipp Griewank 2011-01-07
!!
!!
!!
!!
!!
!!
MODULE mo_functions

  USE mo_parameters

  IMPLICIT NONE
  PUBLIC

CONTAINS
  !>
  !! Calculates the physical density for given S and T.
  !!
  !! 
  !!
  !! Although the model treats Salinity as a massless tracer, sometimes it is necessary to determine the exact density for specific purposes. 
  !! First implemented to calculate simple turbulence between liquid layer and ocean.
  !! Uses following simplification of Frank J. Millero and Alain Poisson 1981:
  !! Density = density_0 +A*S+B*S**1.5
  !! 
  !!
  !! @par Revision History
  !! Started by Philipp Griewank (2011-02-24)
  FUNCTION func_density(T,S) RESULT (density)

    REAL(wp), INTENT(in) :: T,S
    REAL(wp)             :: density
    REAL(wp)             :: A,B,density_0

    density_0 = 999.842594_wp+6.8_wp/100._wp*T
    A         = 0.825_wp
    B         =-5.7_wp/1000._wp
    density   = density_0+A*S+B*MAX(S,0._wp)**1.5_wp

  END FUNCTION func_density

  !>
  !! Calculates the freeboard of the 1d ice column.
  !!
  !! 
  !!
  !! The freeboard is calculated by first finding out which layer is at water level, and then finding out how deep the layer is submerged. For the correct freeboard the mass above water equals the buoyancy of the submerged part. Since the density of each layer is constant, step two can be calculated explicitly.
  !! The freeboard is the distance from the top of the ice to the water level. 
  !! If snow pushes the ice underwater the freeboard becomes negative
  !!
  !! @par Revision History
  !! Built to spill by Philipp Griewank (2011-01-07)
  !! Negative freeboard included by Philipp Griewank (2011-01-09)
  !! Patched bug by Philipp Griewank (2011-03-10)
  !! Add freeboard_snow_flag calculation of snow mass, check the code for further explanations by Niels Fuchs, MPIMET (2017-03-91)
  
  FUNCTION func_freeboard(N_active,Nlayer,psi_s,psi_g,m,thick,m_snow,freeboard_snow_flag) RESULT (freeboard)

    INTEGER,  INTENT(in)                   :: N_active,Nlayer,freeboard_snow_flag

    REAL(wp), INTENT(IN),DIMENSION(Nlayer) :: psi_s,psi_g,m,thick

    REAL(wp), INTENT(in)                   :: m_snow
    REAL(wp)                               :: freeboard, snowmass   !< Niels, 2017 add: snowmass
    REAL(wp)                               :: test1,test2
    INTEGER                                :: k
    
    !< Niels, 2017 add: if loop, important for lab studies with large snow thickness only on a small area of the ice
    IF (freeboard_snow_flag==0) THEN
       snowmass = m_snow
    ELSE
       snowmass = 0._wp
    END IF
    
    !step 0.5, checking if the freeboard is negative (mass_snow greater then buoyancy of ice)
    IF(snowmass>SUM(psi_s(1:N_active)*thick(1:N_active))*(rho_l-rho_s)+SUM(psi_g(1:N_active)*thick(1:N_active))*rho_l)    THEN
       !snow underwater
       test2=SUM(psi_s(1:N_active)*thick(1:N_active))*(rho_l-rho_s)+SUM(psi_g(1:N_active)*thick(1:N_active))*rho_l
       freeboard=test2-snowmass
       freeboard=freeboard/rho_l



    ELSE
       !snow above water
       !Step 1, finding the layer___________________________________________________
       test1 =0._wp
       test2 =1._wp


       k=0
       DO WHILE(test1<test2)
          k=k+1
          test2=SUM(psi_s(k+1:N_active)*thick(k+1:N_active))*(rho_l-rho_s)+SUM(psi_g(k+1:N_active)*thick(k+1:N_active))*rho_l
          test1=SUM(m(1:k))+snowmass
       END DO

       !Step 2, find out where in the layer the waterline is________________________
       test1=SUM(m(1:k-1))+snowmass
       test2=test2 !stays the same

       freeboard=test2-test1+(rho_l-m(k)/thick(k))*thick(k)
       freeboard=freeboard/rho_l
       freeboard=freeboard+SUM(thick(1:k-1))

    END IF

  END FUNCTION func_freeboard


  !>
  !! Calculates the albedo.
  !!
  !! 
  !!
  !! Calculates the albedo according to top conditions.
  !! This is not a good albedo scheme! It is only a quick approach.
  !! Non-continuous switching between wet and dry ice. 
  !! Linear change from wet ice to water.
  !! Linear change from ice_dry snow for snow thinner than 30cm. 
  !!
  !! psi_l(1)> 0.75 water
  !! psi_l(1)> 0.6 linear change from wet ice to water
  !! psi_l(1)> 0.2  wet ice
  !! psi_l(1)< 0.2 -> dry ice
  !! T_snow  = 0 -> wet snow
  !! T_snow  < 0 -> dry snow
  !!
  !!
  !!
  !!
  !! @par Revision History
  !! Built to spill by Philipp Griewank (2011-02-12)

  FUNCTION func_albedo(thick_snow,T_snow,psi_l,thick_min,albedo_flag) RESULT (albedo)
    REAL(wp), INTENT(in) :: thick_snow,T_snow,psi_l,thick_min
    INTEGER,  INTENT(in) :: albedo_flag
    REAL(wp)             :: albedo
    REAL(wp)             :: ice_dry,ice_wet,snow_dry,snow_wet,water

    ice_dry  = 0.75
    ice_wet  = 0.6
    snow_dry = 0.85
    snow_wet = 0.75
    water    = 0.2

    IF (thick_snow>thick_min) THEN
       IF (T_snow<-0.01) THEN
          albedo = snow_dry
       ELSE
          albedo = snow_wet
       END IF

       !snow depth approach
       albedo = ice_dry+(albedo-ice_dry)*MIN(1._wp,thick_snow/0.3_wp)

    ELSE
       IF (psi_l>0.9_wp) THEN
          albedo = water
       ELSE IF (psi_l>0.6_wp ) THEN
          albedo = ice_wet + (water-ice_wet)*((psi_l-0.6_wp)/0.3_wp)
       ELSE IF (psi_l>0.2_wp ) THEN
          albedo = ice_wet
       ELSE
          albedo = ice_dry
       END IF
    END IF

    !Simple version
    If (albedo_flag==1) then
      IF (thick_snow>thick_min) THEN
         IF (T_snow<-0.01) THEN
            albedo = snow_dry
         ELSE
            albedo = snow_wet
         END IF
      ELSE  
         IF (psi_l<0.8) THEN 
            albedo = ice_dry
         ELSE
            albedo = water
         END IF
      END IF
    END IF

  END FUNCTION func_albedo

  !>
  !! Calculates the oxygen saturation as a function of salinity and temperature.
  !!
  !! Calculates the  concentration of oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere
  !! The value should be umol/kg.
  !! I switched to the solubility of nitrogen, oxygen and argon in water and sea wate from Weiss R.F. 1970 because I couldn't get the other one to work out
  !! @par Revision History
  !! Written by Dr. Philipp Griewank (2014-02-25)

  FUNCTION func_sat_O2(T,S_bu) RESULT (sat_O2)
    REAL(wp), INTENT(in) :: T,S_bu
    REAL(wp)             :: sat_O2, TT

    TT = T+273.16


    sat_O2 = 1.42905_wp*EXP(-173.4292_wp   +24963.39_wp/TT  +143.3483*LOG(TT/100._wp) -0.218492_wp*TT)
    sat_O2 = sat_O2*EXP(S_bu*(-0.033096_wp +0.00014259_wp*TT - 0.0017_wp*TT**2/10000._wp)    )
    sat_O2 = sat_O2/0.032
  END FUNCTION func_sat_O2


  !>
  !! Calculates the freezing temperature.
  !! Salt_flag determines if either ocean salt or NAcl is used.
  !!
  !! @par Revision History
  !! Written to procrastinate by Philipp Griewank (2011-05-05)

  FUNCTION func_T_freeze(S_bu,salt_flag) RESULT (T_freeze)
    REAL(wp), INTENT(in) :: S_bu
    INTEGER,  INTENT(in) :: salt_flag
    REAL(wp)             :: T_freeze


    IF (salt_flag==2) THEN
       T_freeze = -0.0592_wp*S_bu -9.37*S_bu**2.0 -5.33*10.0**(-7.0)*S_bu**3.0 
    ELSE IF (salt_flag==1) THEN
       T_freeze = -0.0575_wp*S_bu +1.710523*1e-3*S_bu**1.5 -2.154996*1e-4*S_bu**2.0
    END IF
  END FUNCTION func_T_freeze


  !>
  !! Calculates the incoming shortwave and other fluxes according to p. 193-194 PhD Notz.
  !!
  !! 
  !!
  !! Simplified version of the Untersteiner Fluxes.
  !! Returns only two fluxes as a function of time. 
  !! Simplified Year, 12 months of 30 days.
  !! fl_sw is set to zero for November till February
  !! Returns fluxes for day with day zero being 1. Jan. 
  !! Depending on when the run starts the time should be modified when calling 
  !!
  !!
  !!
  !! @par Revision History
  !! Ripped from Dirk by Philipp Griewank (2011-02-13)

  SUBROUTINE sub_notzflux(time,fl_sw,fl_rest)
    REAL(wp), INTENT(in) :: time
    REAL(wp) :: day
    REAL(wp), INTENT(out) :: fl_sw,fl_rest

    day    = time/86400._wp
    day    = day
    DO WHILE (day>360)
       day = day-360
    END DO


    fl_sw   = 314._wp*EXP(-0.5_wp*((day-164._wp)/47.9)**2._wp)
    fl_rest = 118._wp*EXP(-0.5_wp*((day-206._wp)/53.1)**2._wp)+179._wp

    IF (day<60. .OR. day>300.) THEN
       fl_sw = 0.0_wp
    END IF

  END SUBROUTINE sub_notzflux

  !>
  !! Reads in data for atmoflux_flag ==2.
  !!
  !! Standard setup used for testcase 4 and all Griewank & Notz 2013/14 reanalysis forced runs is 4.5 years of three hourly values of shortwave incoming, longwave incoming, two meter T, and total precipitation. 
  !! Data is read from ascii files and stored in long 1D arrays. 
  !! ERA-interim derived input files in the standard length for various Arctic locations are located under /input/ERA/
  !! Latent and sensible heat fluxes are not included, but could be added if needed. 
  !!
  !!
  !!
  !! @par Revision History
  !! Moved here from mo_grotz by Philipp Griewank (2014-04-20)

  SUBROUTINE  sub_input(length_input,fl_sw_input,fl_lw_input,T2m_input,precip_input,time_input)
    INTEGER,  INTENT(in)               :: length_input
    INTEGER                            :: k
    REAL(wp), DIMENSION(:), allocatable, INTENT(out) :: fl_sw_input,fl_lw_input,T2m_input,precip_input,time_input
     ALLOCATE(fl_sw_input(Length_Input),fl_lw_input(Length_Input),T2m_input(Length_Input),precip_input(Length_Input),&
            &time_input(Length_Input))

       OPEN(1234,file='flux_lw.txt.input',status='old')
       READ(1234,*)fl_lw_input
       CLOSE(1234)
       OPEN(1234,file='flux_sw.txt.input',status='old')
       READ(1234,*)fl_sw_input
       CLOSE(1234)
       OPEN(1234,file='T2m.txt.input'    ,status='old')
       READ(1234,*)T2m_input
       CLOSE(1234)
       OPEN(1234,file='precip.txt.input' ,status='old')
       READ(1234,*)precip_input
       CLOSE(1234)
       DO k=1,Length_Input
          time_input(k) = (REAL(k)-1._wp)*3600._wp*3._wp
       END DO

  END SUBROUTINE sub_input


  !>
  !! Calculates salt and tracer mixing between lowest layer and underlying water.
  !!
  !! 
  !!Very simple turbulence assumption which mixes the lowest layer with the underlying water.
  !!Based on assumption that there is a constant amount of turbulence A.
  !!This turbulence is amplified when the lowest layer is denser then the ocean mixed layer.
  !!And also dampened when the lowest layer is less dense then the mixed layer.
  !!Assumption; turb=A*exp(B(density_layer-density_ocean))
  !!A and B set in parameters.i
  !!A = turb_A , B = turb_B
  !!
  !!
  !!
  !! @par Revision History
  !! Moved from grotz by Philipp Griewank (2014-04-2)

  SUBROUTINE sub_turb_flux(T_bottom,S_bu_bottom,T,S_abs,m,dt,N_bgc,bgc_bottom,bgc_abs)
    REAL(wp), INTENT(in)                                :: T_bottom, S_bu_bottom, T, m, dt 
    INTEGER,  INTENT(in)                                :: N_bgc
    REAL(wp), INTENT(inout)                             :: S_abs
    REAL(wp), INTENT(inout), DIMENSION(N_bgc), OPTIONAL :: bgc_abs
    REAL(wp), INTENT(in),    DIMENSION(N_bgc), OPTIONAL :: bgc_bottom
    REAL(wp)                                            :: turb !< mass flow between lowest layer and underlying water over the timestep


    turb  = turb_A*EXP(turb_B*(-func_density(T_bottom,S_bu_bottom)+func_density(T,S_abs/m)))*dt
    S_abs = S_abs-turb*(S_abs/m-S_bu_bottom)
    IF( PRESENT( bgc_bottom ) ) THEN
       bgc_abs(:) = bgc_abs(:)-turb*(bgc_abs(:)/m-bgc_bottom(:))
    END IF


  END SUBROUTINE sub_turb_flux



  !>
  !! Calculates the thickness of the meltwater film.
  !!
  !! 
  !!
  !! If the top ice layer is being melted (T_top>T_freeze) it is assumed that a thin meltwater film appears at the top.
  !! The thickness of this film is determined by the amount of incoming heat and diffusive transport.
  !! The incoming heat is an input (fl_q(1)) and the diffusive heat is (T(1)-T_freeze)/R. 
  !! See the thermodynamics section for R.
  !! The thickness of the meltlayer is determined by dividing the heat intake of the meltwater film by the amount of latent heat needed to melt the solid fraction of the top layer. 
  !! If the solid fractions sinks below a given threshold (psi_s_top_min) a different approach is used.
  !! The melt thickness is then calculated by assuming that the ice below the meltwater film has a solid fraction of psi_s_top_min.
  !! Although the thickness can be reduced, variations of mass, salinity and enthalpy are calculated in the flushing subroutine.
  !!
  !!
  !!
  !! @par Revision History
  !! Introduced by Philipp Griewank (2011-05-09)

  SUBROUTINE sub_melt_thick(psi_l,psi_s,psi_g,T,T_freeze,T_top,fl_Q,thick_snow,dt,melt_thick,thick,thick_min)
    REAL(wp), INTENT(in)    :: psi_l,psi_s,psi_g,T,T_freeze,T_top,fl_Q,thick_snow,dt,thick_min
    REAL(wp), INTENT(out)   :: melt_thick
    REAL(wp), INTENT(inout) :: thick

    melt_thick = 0.0_wp

    !Assumption Q flux from film to layer fl_mid = 2._wp*(psi_l(1)*k_l+psi_s(1)*k_s)/thick*(T_freeze-T(1)).
    !The energy different between fl_Q(1) and fl_mid is the amount of energy that goes to melt the meltwater film.
    !This only applies when the snow thickness is thin
    IF (thick_snow<thick_min .AND. T_top>=T_freeze ) THEN

       melt_thick = -fl_Q-2._wp*(psi_l*k_l+psi_s*k_s)/thick*(T_freeze-T)
       melt_thick = melt_thick*dt/MAX((latent_heat*rho_s*psi_s),0.000000000000001_wp)
       melt_thick = MIN(psi_l*thick,melt_thick)

    END IF


    !The other way that top meltwater is formed is when the solid fraction of the top layer is less then
    !a prescribed value (psi_s_top_min).
    !The melt layer is assigned so that the rest of the layer achieves a solid fraction that equals
    !the prescribed (psi_s_top_min)
    !When snow is present meltwater floods the snow layer..
    !At the moment this melt_thick assignment has higher priority then the T_top method. 

    IF (psi_s<psi_s_top_min ) THEN
       melt_thick = thick*(1._wp-psi_s/psi_s_top_min)
    END IF


    !Making sure that air percentage is kept to gas_snow_ice2
    IF(melt_thick>0.0_wp .AND. psi_g>gas_snow_ice2) THEN
       IF(melt_thick>(psi_g-gas_snow_ice2)*thick) THEN
          melt_thick = melt_thick-(psi_g-gas_snow_ice2)*thick
          thick      = thick*(1._wp-(psi_g-gas_snow_ice2))
       ELSE
          thick      = thick-melt_thick
          melt_thick = 0._wp
       END IF
    END IF

  END SUBROUTINE sub_melt_thick
 
  !>
  !! Calculates how the meltwater film interacts with snow.
  !!
  !! 
  !! Is activated when a thin snow layer (thinner then thick_min) is on top of meltwater.
  !! The snow is flooded and turned into ice.
  !! 
  !!
  !!
  !!
  !! @par Revision History
  !! Put together by Philipp Griewank (2011-10-17)

  SUBROUTINE sub_melt_snow(melt_thick,thick,thick_snow,H_abs,H_abs_snow,m,m_snow,psi_g_snow)
    REAL(wp), INTENT(inout) :: melt_thick,thick,thick_snow,H_abs,H_abs_snow,m,m_snow,psi_g_snow
    REAL(wp)                :: shift

    shift = 1._wp/MAX(psi_g_snow,0.01_wp)*melt_thick

    IF (shift>=thick_snow) THEN
       melt_thick = melt_thick-thick_snow*psi_g_snow
       H_abs      = H_abs+H_abs_snow
       m          = m+m_snow
       thick      = thick+(1._wp-psi_g_snow)*thick_snow

       thick_snow = 0.0_wp
       m_snow     = 0.0_wp
       H_abs_snow = 0.0_wp
    ELSE

       H_abs      = H_abs +shift/thick_snow*H_abs_snow
       H_abs_snow = H_abs_snow -shift/thick_snow*H_abs_snow

       m      = m +shift/thick_snow*m_snow
       m_snow = m_snow -shift/thick_snow*m_snow

       thick      = thick +shift-melt_thick
       thick_snow = thick_snow -shift

       melt_thick = 0.0_wp

    END IF


  END SUBROUTINE sub_melt_snow

end module mo_functions
