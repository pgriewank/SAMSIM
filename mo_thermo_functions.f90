!>
!!  Contains subroutines and functions related to multi-phase thermodynamics.
!!
!!  See the subroutine and function descriptions for details.
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
!! Started by Philipp Griewank 2010-07-08
!! Add function for styropor cover by Niels Fuchs, MPIMET (2017-01-03)
!! Modified salinity functions by Philipp Griewank, Uni K (2018-08-01)

MODULE mo_thermo_functions

  USE mo_parameters


  IMPLICIT NONE

  PUBLIC :: getT
  PUBLIC :: Expulsion
  PUBLIC :: sub_fl_Q
  PUBLIC :: sub_fl_Q_0
  PUBLIC :: func_S_br
  PUBLIC :: func_ddT_S_br
  PUBLIC :: sub_fl_Q_styropor   !> Niels, 2017

  PRIVATE


CONTAINS

  !>
  !! Determines equilibrium Temperature of a layer for given S_bu and H as well as solid fraction.
  !!
  !! The temperature of a fully liquid layer is used to see if the resulting brine salinity is lower than the bulk salinity.
  !! After checking if the layer is a fluid or a mushy layer the temperature is calculated
  !! by solving f(T) = 0 using the Newton method.
  !! f(T) = -latent_heat-H+latent_heat*S_bu/S_br(T) + c_s*T+c_s_beta*T^2/2
  !! f'(T) = c_s+c_s_beta*T-latent_heat*S_bu*S_br'(T)/S_br^2
  !! Described in Notz2005, subsubsection 5.6.1.
  !! See func_S_br(T) and func_ddT_S_br(T).
  !! First guess T_0 must be given, low first guess lead to overshooting which would lead to very high Temperatures.
  !! To avoid this, an if loop sets T to freezing T when T>0.
  !! Freezing T is also calculated at the beginning using the Newton-Method.
  !! If S_bu<0.001 then it is treated as pure ice.
  !!
  !! @par Revision History
  !! first version by Philipp Griewank (2010-07-13)
  !! Freezing temperature is calculated and introduced if T goes above 0 by Philipp Griewank (2010-07-13) 
  !! Added if loops to deal with saltless ice by Philipp Griewank (2010-11-27)
  SUBROUTINE getT(H,S_bu,T_in,T,phi,k)


    REAL(wp), INTENT(in)           :: H      !<  Enthalpy [J/kg]
    REAL(wp), INTENT(in)           :: S_bu   !<  Bulk Salinity [g/kg]
    REAL(wp), INTENT(in)           :: T_in   !<  input Temperature for T_0 [C]

    REAL(wp), INTENT(out)          :: T      !<  Temperature [C]
    REAL(wp), INTENT(out)          :: phi    !<  solid fraction

    REAL(wp)                       :: T_0    !<  First guess and iterative Temperature
    REAL(wp)                       :: f,ddT_f!<  f(t) and f'(T) as described above
    REAL(wp)                       :: T_fr   !<  Freezing temperature
    INTEGER,  INTENT(in), OPTIONAL :: k
    INTEGER                        :: i



    T = H/c_l

    IF(func_S_br(T,S_bu)>S_bu .AND. S_bu>0.001_wp) THEN 

       !After it has been checked that ice is present, the freezing temp is calculated
       T_fr = -1._wp

       DO WHILE(ABS(func_S_br(T_fr)/S_bu-1.0)>0.0001) 
          T_0   = T_fr
          f     = func_S_br(T_0)-S_bu
          ddT_f = func_ddT_S_br(T_0)
          T_fr  = T_0-f/ddT_f
       END DO

       T_0    = T_in
       f      = -latent_heat-H+latent_heat*S_bu/MAX(func_S_br(T_0),0.000000001_wp) + c_s*T_0+c_s_beta*T_0*T_0/2._wp
       ddT_f  =  c_s+c_s_beta*T_0-latent_heat*S_bu*func_ddT_S_br(T_0)/MAX(func_S_br(T_0)**2._wp,0.0000000001_wp)
       T      =  T_0-f/ddT_f
       i = 0
       DO WHILE (ABS(f)>1_wp)
          T_0 = T
          IF(T_0>0. .OR. T_0<-200.) THEN 
             T_0 = T_fr 
          END IF
          f      = -latent_heat-H+latent_heat*S_bu/MAX(func_S_br(T_0),0.0000000001_wp) + c_s*T_0+c_s_beta*T_0*T_0/2._wp
          ddT_f  = c_s+c_s_beta*T_0-latent_heat*S_bu*func_ddT_S_br(T_0)/MAX(func_S_br(T_0)**2,0.0000000001_wp)
          T      = T_0-f/ddT_f
          i      = i+1


          IF (i>200) THEN
             PRINT*,'converging problem',T,T_0,f/ddt_f,i
          END IF

          IF (i==260) THEN
             OPEN(99,FILE = 'errorlog.txt',STATUS = 'replace')
             WRITE(99,*)'In layer',k,'getT does not converge for H,S_bu = ',H,S_bu
             WRITE(99,*)'Such problems have occured when Temperatures are well below -100C or have bad intial values'
             WRITE(99,*)'Negative salinity values or very very low Salinity can also be a problem.'
             WRITE(99,*)'Often the problem can be traced to a too big time step.'
             WRITE(*,*)'better luck next time, check errorlog.txt'
             CLOSE(99)
             STOP 99
          END IF
       END DO
       phi = 1._wp-S_bu/func_S_br(T,S_bu)

    ELSE IF (S_bu<0.001_wp) THEN
       IF (H>0.0_wp) THEN
          phi = 0._wp
          T   = H/c_l
       ELSE IF (H<=-latent_heat) THEN
          phi = 1._wp
          T   = (H+latent_heat)/c_s
       ELSE IF (H<=0.0_wp .AND. -latent_heat<H) THEN 
          T   = 0._wp
          phi = -H/latent_heat
       END IF
    ELSE
       phi = 0.0_wp
    END IF


  END SUBROUTINE getT


  !>
  !! Determines Brine flux expelled from out of a layer due to freezing.
  !!
  !! If the volume of ice and brine exceed the Volume of the layer brine is expelled.
  !! The volume of the ejected brine is calculated and exported.
  !! The volume fractions are also calculated.
  !!
  !! @par Revision History
  !! first version by Philipp Griewank, (2010-07-19)
  !! changes to mass, Enthalpy and Salinity are now computed in subroutine mass_transfer by Philipp Griewank, (2010-08-24)
  !!
  SUBROUTINE Expulsion(phi,thick,m,psi_s,psi_l,psi_g,V_ex)


    REAL(wp), INTENT(in)    :: phi
    REAL(wp), INTENT(in)    :: thick

    REAL(wp), INTENT(inout) :: m
    REAL(wp), INTENT(out)   :: psi_s,psi_l,psi_g,V_ex
    REAL(wp) :: V_s,V_l

    V_s = m*phi/rho_s
    V_l = m*(1._wp-phi)/rho_l

    IF(V_s+V_l>thick) THEN
       V_ex = V_l+V_s-thick
    ELSE
       V_ex = 0._wp
    END IF

    psi_s  = V_s/thick
    psi_l  = (V_l-V_ex)/thick
    psi_g  = (thick-V_l-V_s+V_ex)/thick

    IF (psi_l<0._wp) THEN
       psi_l = 0.0_wp
    END IF
    IF (psi_g<0) THEN
       psi_g = 0.0_wp
    END IF

  END SUBROUTINE Expulsion

  !>
  !! Determines conductive heat flux between two layers.
  !!
  !! Details can be found in Notz 2005, especially equation 5.7.
  !! The gas volume is assumed to have no thermal properties at all.
  !! First the thermal resistance R is calculated using the approximated thermal conductivity of the mushy layer (see Notz 2005 eq. 3.41.).
  !! Then the heat flux Q is simply (T_1-T_2)/R
  !! "_1" denotes the upper layer and "_2" the lower layer. A positive heat flux is from lower to upper layer.
  !!
  !! @par Revision History
  !! First version by Philipp Griewank (2010-07-21)
  !!
  SUBROUTINE sub_fl_Q(psi_s_1,psi_l_1,psi_g_1,thick_1,T_1,psi_s_2,psi_l_2,psi_g_2,thick_2,T_2,fl_Q)


    REAL(wp), INTENT(in)   :: psi_s_1,psi_l_1,psi_g_1
    REAL(wp), INTENT(in)   :: psi_s_2,psi_l_2,psi_g_2
    REAL(wp), INTENT(in)   :: thick_1,thick_2
    REAL(wp), INTENT(in)   :: T_1,T_2

    REAL(wp), INTENT(out)  :: fl_Q

    REAL(wp)               :: R       !< Thermal resistance
    REAL(wp)               :: k_1,k_2 !< Thermal conductivity


    k_1 = psi_s_1*k_s+psi_l_1*k_l+psi_g_1*0._wp
    k_2 = psi_s_2*k_s+psi_l_2*k_l+psi_g_2*0._wp


    R   = thick_1/(2._wp*k_1)+thick_2/(2._wp*k_2)

    fl_Q = (T_2-T_1)/R

  END SUBROUTINE sub_fl_Q


  !>
  !! Determines conductive Heat flux between layer and boundary temperatures.
  !!
  !! Details can be found in Notz 2005, especially equation 5.10 and 5.11.
  !! The gas volume is assumed to have no thermal properties.
  !! direct_flag denotes if the boundary layer is above or below the layer.
  !! 1 : =  layer above boundary
  !! -1: =  layer below boundary
  !!
  !! @par Revision History
  !! first version by Philipp Griewank (2010-07-21)
  !!
  SUBROUTINE sub_fl_Q_0(psi_s,psi_l,psi_g,thick,T,T_bound,direct_flag,fl_Q)


    REAL(wp), INTENT(in)    :: psi_s,psi_l,psi_g
    REAL(wp), INTENT(in)    :: thick
    REAL(wp), INTENT(in)    :: T,T_bound !< T_bound temperature of boundary layer
    INTEGER                 :: direct_flag
    REAL(wp), INTENT(out)   :: fl_Q

    REAL(wp)   :: R !< Thermal resistance
    REAL(wp)   :: k !< Thermal conductivity

    k = psi_s*k_s+psi_l*k_l+psi_g*0._wp
    R = thick/(2._wp*k)


    IF (direct_flag==1) THEN
       fl_Q = (T_bound-T)/R
    ELSE IF (direct_flag==-1) THEN
       fl_Q = (T-T_bound)/R
    ELSE IF (direct_flag .NE. -1 .AND. direct_flag .NE. 1 ) THEN
       OPEN(99,FILE = 'errorlog.txt',STATUS = 'replace')
       WRITE(99,*)'sub_fl_Q_0 is given unavailable direct_flag'
       CLOSE(99)
       STOP 98

    END IF

  END SUBROUTINE sub_fl_Q_0
  
  !> Niels, 2017 add: 
  !! Determines conductive Heat flux below styropor cover
  !!
  !! Standard approach.
  !!
  !! @par Revision History
  !! first version by Niels Fuchs, MPIMET (2017-01-03)
  !!
  SUBROUTINE sub_fl_Q_styropor(k_styropor,fl_Q)

    REAL(wp),  INTENT(in)  :: k_styropor
    REAL(wp),  INTENT(inout) :: fl_Q
    REAL(wp)               :: R


    R      = k_styropor
    fl_Q   = fl_Q*R

  END SUBROUTINE sub_fl_Q_styropor


  !>
  !! Computes salinity of brine pockets for given temperature in Celsius of mushy layer
  !!
  !! Subroutine computes one S_br for one given T in Celsius by third-order polynomial.
  !! NaCl solutions and seawater produce slight variations. Which solution is used is determined by salt_flag.  
  !! S_br = c1+c2*T+c3*T^2+c4*T^3
  !! Based on Notz 2005 p. 36
  !! 
  !! For T below -20 we simply use a linear extension based on "Composition of sea ice and its tensile strength".
  !! The actual salinity function is much more complicated and depends on the salt composition, but the linear fit is far better than using the polynomial fit.
  !!
  !! The NaCL precipitates at -22, leading to ice/salt kristall mix below -22. Ideally the whole code would be modified to take the non-continous transition at -22 but given that there is currently little interest I can't be bothered to put in the effor.
  !! 
  !! @par Revision History
  !! First version by Philipp Griewank (2010-07-12)
  !! Changed to go through 0 by Philipp Griewank (2014-05-07)
  !! Added linear bit by Philipp Griewank (2018-07-22)

  FUNCTION func_S_br(T,S_bu) RESULT (S_br)
    USE mo_data, ONLY:salt_flag
    IMPLICIT NONE


    REAL(wp)                      :: c1,c2,c3,c4,T_crit,S_0,ddT_S_0
    REAL(wp), INTENT(in)          :: T  !<  Temperature in Celsius
    REAL(wp), INTENT(in),OPTIONAL :: S_bu

    REAL(wp)                      :: S_br !<  Brine salinity



    T_crit = -20._wp
    IF (salt_flag==1) THEN
       !Saltwater coefficients
       c1  =  0.0_wp
       c2  =  -21.4_wp
       c3  =  -0.886_wp
       c4  =  -0.0170_wp
    ELSE IF (salt_flag==2) THEN
       !NaCl coefficients
       c1  =  0.0_wp
       c2  =  -17.6_wp
       c3  =  -0.389_wp
       c4  =  -0.00362_wp

    END IF

    S_br = c1 + c2*T + c3*T**2._wp + c4*T**3._wp
    
    IF (T.lt.T_crit) THEN

        S_0     = c1 + c2*T_crit + c3*T_crit**2._wp + c4*T_crit**3._wp
        ddT_S_0 = c2 + 2._wp*c3*T_crit + 3._wp*c4*T_crit**2._wp
        S_br    = S_0+ddT_S_0*(T-T_crit)

    END IF


    IF(PRESENT(S_bu))THEN
       IF( S_br<S_bu) THEN
          S_br = S_bu
       END IF
    ELSE
    END IF

  END FUNCTION func_S_br

  !>
  !! Computes temperature derivative of brine pocket salinity for given temperature in Celsius of mushy layer
  !!
  !! Subroutine computes one ddT_S_br for one given T
  !! NaCl solutions and seawater produce slight variations. Which solution is used is specified by salt_flag.
  !! Based on Notz 2005 p. 36
  !! ddT_S_br = c2+2*c3*T+3*c4*T^2
  !!
  !! For T below -20 we simply use a linear extension based on "Composition of sea ice and its tensile strength".
  !! The actual salinity function is much more complicated and depends on the salt composition, but the linear fit is far better than using the polynomial fit.
  !!
  !!
  !! The NaCL precipitates at -22, leading to ice/salt kristall mix below -22. Ideally the whole code would be modified to take the non-continous transition at -22 but given that there is currently little interest I can't be bothered to put in the effor.
  !!
  !! @par Revision History
  !! First version by Philipp Griewank (2010-07-13)
  !! Added linear bit by Philipp Griewank (2018-07-22)

  FUNCTION func_ddT_S_br(T) RESULT (ddT_S_br)
    USE mo_data, ONLY:salt_flag
    IMPLICIT NONE


    REAL(wp)             :: c2,c3,c4,T_crit

    REAL(wp), INTENT(in) :: T  !<  Temperature in Celsius

    REAL(wp)             :: ddT_S_br !<  derivative of Brine salinity


    T_crit = -20._wp
    IF (salt_flag==1) THEN
       !Saltwater coefficients
       c2  =  -21.4_wp
       c3  =  -0.886_wp
       c4  =  -0.0170_wp
    ELSE IF (salt_flag==2) THEN
       !NaCl coefficients
       c2  =  -17.6_wp
       c3  =  -0.389_wp
       c4  =  -0.00362_wp

    END IF

    ddT_S_br = c2 + 2._wp*c3*T + 3._wp*c4*T**2._wp
    
    IF (T.lt.T_crit) THEN

        ddT_S_br = c2 + 2._wp*c3*T_crit + 3._wp*c4*T_crit**2._wp

    END IF

  END FUNCTION func_ddT_S_br
  
END MODULE mo_thermo_functions
