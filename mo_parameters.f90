
!>
!!  Module determines physical constants to be used by the SAMSIM Seaice model.
!!
!!  Many values are taken from Notz 2005, Table 5.2.
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
!!
!! @par Revision History
!! Started by Philipp Griewank 2010-07-08 \n
!! add parameters: heat conductivity of styropor under special sea ice lab conditions and ratio of penetrating melt water by Niels Fuchs, MPIMET (2017-03-01)
!!
!!
MODULE mo_parameters


  IMPLICIT NONE
  PUBLIC
  
  
  
  INTEGER, PARAMETER :: wp =  SELECTED_REAL_KIND(12,307) !< set working precision _wp
  



  REAL, PARAMETER:: pi   = 3.1415_wp
  REAL, PARAMETER:: grav = 9.8061_wp !< gravitational constant [m/s^2]


  ! ! Physical constants
  !
  ! ! Almost none of the parameters are truly constant, but will be assumed as such based on sensitivity studies detailed in Notz 2005 subsection 5.7
  !--------------------------------------------------------------------
  REAL(wp), PARAMETER:: k_s  =  2.2_wp    !< solid heat conductivity [J / m s K] 2.2
  REAL(wp), PARAMETER:: k_l  =  0.523_wp  !< liquid heat conductivity [J / m s K] 0.523
  
  REAL(wp), PARAMETER:: c_s          = 2020.0_wp       !< solid heat capacity [J/ kg K]
  REAL(wp), PARAMETER:: c_s_beta     = 7.6973_wp       !< linear solid heat capacity approximation [J/ kg K^2] c_s = c_s+c_s_beta*T
  REAL(wp), PARAMETER:: c_l          = 3400._wp        !< liquid heat capacity [J/ kg K]
  REAL(wp), PARAMETER:: rho_s        = 920._wp         !< density of solid [kg / m^3]
  REAL(wp), PARAMETER:: rho_l        = 1028.0_wp       !< density of liquid [kg / m^3]
  REAL(wp), PARAMETER:: latent_heat  = 333500._wp      !< latent heat release [J/kg]
  REAL(wp), PARAMETER:: zeroK        = 273.15_wp       !< Zero degrees Celsius in Kelvin [K]
  REAL(wp), PARAMETER:: bbeta        = 0.8_wp*1e-3     !< concentration expansion coefficient  [kg / (m^3 ppt)]
  REAL(wp), PARAMETER:: mu           = 2.55_wp*1e-3    !< dynamic viscosity [kg /m s]
  REAL(wp), PARAMETER:: kappa_l      = k_l/rho_l/c_l   !< heat diffusivity of water
  REAL(wp), PARAMETER:: sigma        = 5.6704_wp*1e-8  !< Stefan Boltzmann constant [W/(m^2*K^4)]


  ! ! Model constants
  !
  ! ! Parameters used by the model which may be based on physical values, but are primarily intended to 
  ! ! Keep the model running well. Physical accuracy is of second order priority.


  !Layer Dynamics
  REAL(wp), PARAMETER:: psi_s_min  =  0.05_wp !<The amount of ice that the lowest layer can have before it counts as an ice layer 
  REAL(wp), PARAMETER:: neg_free   = -0.05_wp !<The distance the freeboard can be below 0 before water starts flooding through cracks. 


  !Gravity Drainage, alternative values included, see Griewank & Notz 2014 for details and 
  REAL(wp), PARAMETER:: x_grav     = 0.000584_wp !0.00051_wp !0.000681_wp  !< links Rayleigh number to grav_drain
  REAL(wp), PARAMETER:: ray_crit   = 4.89_wp      !7.10_wp    !3.23_wp     !< critical Rayleigh number 

  !Flushing and snow cover, serves to stabilize thermal fluxes.
  !Flush_flag =  = 5
  REAL(wp), PARAMETER:: para_flush_horiz = 1.0_wp       !<determines relationship of horizontal flow distance in during flushing (guess 1)
  !Flush_flag =  = 6
  REAL(wp), PARAMETER:: para_flush_gamma = 0.9_wp       !<Strength of desalination  per timestep (guess)

  REAL(wp), PARAMETER:: psi_s_top_min    = 0.40_wp      !<if psi_s is below this value meltwater forms  (guess) 0.4
  !Flooding
  REAL(wp), PARAMETER:: ratio_flood      = 1.50_wp      !<Ratio of flooded to dissolve snow, plays an important role in subroutine flood
  !Freshwater calculation
  REAL(wp), PARAMETER:: ref_salinity     = 34._wp       !<Reference salinity [g/kg] used to calculate freshwater column
  

  !Snow
  REAL(wp), PARAMETER:: rho_snow      = 330._wp !<density of new snow [kg/m**3], !< Niels, 2017 add: can be adjusted to lab values if they are measured
  REAL(wp), PARAMETER:: gas_snow_ice  = 0.10_wp !<volume of gas percentage in new snow ice due to flooding, no longer used
  REAL(wp), PARAMETER:: gas_snow_ice2 = 0.20_wp !<volume of gas percentage in new snow ice due to snow melting (Eicken 95)

  !Radiation        all values token from Notz PhD if not otherwise mentioned.
  REAL(wp), PARAMETER:: emissivity_ice  = 0.95_wp    !< Emissivity of water and ice
  REAL(wp), PARAMETER:: emissivity_snow = 1.00_wp    !< Emissivity of Snow
  REAL(wp), PARAMETER:: penetr          = 0.30_wp    !< Amount of penetrating sw radiation
  REAL(wp), PARAMETER:: extinc          = 2.00_wp    !< Extinction coefficient of ice

  !Bottom turbulence
  REAL(wp), PARAMETER:: Turb_A  = 0.1_wp*0.05_wp*rho_l/86400._wp !< Standard turbulence [kg/s] WARNING no source, just set so that 5cm of water are overturned each day
  REAL(wp), PARAMETER:: Turb_B  = 0.05_wp                        !< Exponential turbulence slope [m**3/kg] WARNING no source, simple guess


  !Limitations
  REAL(wp) :: max_flux_plate = 50.0 !< Maximal heating rate of a heating plate

  !Snow melt process
  REAL(wp) :: k_snow_flush = 0.75_wp !< Niels, 2017 add:  Percentage of excess liquid water content in the snow that is used for flushing instead of forming slush

  REAL(wp) :: k_styropor = 0.8_wp !< Niels, 2017 add:  heat conduction of styropor (empirical value to fit measurement data)

END MODULE mo_parameters

