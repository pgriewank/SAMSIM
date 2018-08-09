!>
!! Contains various subroutines for flushing.
!!
!! Which subroutine is called is determined by flush_flag.
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
!!
!! @par Revision History
!! Sang into existence for very 1D column by Philipp Griewank, IMPRS (2010-10-20) \n
!! First stable release by Philipp Griewank, IMPRS (2010-11-27) \n
!! Freeboard calculation outsourced to mo_functions by Philipp Griewank, IMPRS (2010-11-27) \n
!! Drainage through cracks is added by Philipp Griewank, IMPRS (2011-02-24) \n
!! Changes in subroutine flush3 by Niels Fuchs, MPIMET (2017-03-01) \n

MODULE mo_flush

  USE mo_parameters
  USE mo_thermo_functions, ONLY: func_S_br
  USE mo_functions,        ONLY: func_density
  USE mo_output,           ONLY: output_raw_lay
  USE mo_mass,             ONLY: mass_transfer
  IMPLICIT NONE

  PRIVATE
  PUBLIC::flush3,flush4


CONTAINS


  !>
  !! Subroutine for complex flushing.
  !!
  !! 
  !! Each layer splits the flushing brine into a fraction that moves downward, and a fraction that leaves the ice.
  !! A fraction of the top layer is considered melt water. 
  !! This approach uses hydraulic resistivity R = mu*thick/perm .
  !! The hydraulic head is assumed to be the freeboard. 
  !! The vertical resistance R_v of each layer is a determined by its viscosity * thickness divided by it's permeability.
  !! Additionally, each layer is given horizontal resistivity R_h.
  !! It is assumed that there is an average length horizontally which brine needs to flow to reach a drainage feature in the ice.
  !! We assume this length is a linear function of the ice thickness.
  !! The only tuning parameter is para_flush_horiz.
  !! The total resistance of layer i to the bottom is R.
  !! 
  !! For flush_heat_flag==2 the amount of heat which leaves by dynamics from the lowest layer is added to the lowest layer to keep results comparable to the other approaches.
  !! See PhD Griewank for details
  !!
  !! @par Revision History
  !! Invented by Philipp Griewank, IMPRS (2012-06-15) \n
  !! Trying to add brine fluxes by Philipp Griewank, IMPRS (2014-02-01) \n
  !! Changed: Permeability calculation (only for snow_flush_flag==1), hydraulic head and output data by Niels Fuchs, MPIMET (2017-03-01) \n
  !! 
  !!
  SUBROUTINE flush3 (freeboard,psi_l,thick,thick_0,S_abs,H_abs,m,T,dt,Nlayer,N_active,&
       & T_bottom,S_bu_bottom,melt_thick,debug_flag,flush_heat_flag,melt_err,perm,flush_v,flush_h,&
       & psi_g,thick_snow,rho_l,snow_flush_flag,fl_brine_bgc)
    INTEGER,                      INTENT(in)    :: Nlayer,debug_flag,flush_heat_flag,snow_flush_flag !< Niels, 2017 add: snow_flush_flag
    INTEGER,                      INTENT(inout) :: N_active
    REAL(wp), DIMENSION(Nlayer),  INTENT(in)    :: T   !< Niels, 2017 add: moved psi_l -> INTENT(inout)
    REAL(wp),                     INTENT(in)    :: dt,T_bottom,S_bu_bottom,freeboard,thick_0,thick_snow,rho_l
    REAL(wp), DIMENSION(Nlayer),  INTENT(inout) :: psi_l,S_abs,H_abs,m,thick,psi_g !< Niels, 2017 add: psi_l, psi_g
    REAL(wp),                     INTENT(inout) :: melt_thick, melt_err
    REAL(wp), DIMENSION(N_active), INTENT(inout)  :: flush_v         !< mass of vertically flushed brine of each layer [kg] !< Niels, 2017 add: inout
    REAL(wp), DIMENSION(N_active), INTENT(inout)  :: flush_h         !< mass of brine which leaves the ice of each layer [kg] !< Niels, 2017 add: inout
    REAL(wp), DIMENSION(N_active)               :: R_h,R_v,R       !< hydraulic resistivity, see subroutine description []
    REAL(wp), DIMENSION(Nlayer), INTENT(out)    :: perm !< Niels, 2017 add: out
    REAL(wp), DIMENSION(Nlayer)                 :: S_bu
    REAL(wp), DIMENSION(Nlayer+1)               :: fl_m
    REAL(wp)                                    :: const           !< horizontal flow, determined by thickness and para_flush_horiz [m]
    REAL(wp)                                    :: flush_total     !< total amount of brine in motion [kg]
    REAL(wp)                                    :: loss_S_abs,loss_H_abs
    INTEGER                                     :: k
    REAL(wp), DIMENSION(Nlayer+1,Nlayer+1), INTENT(inout),OPTIONAL :: fl_brine_bgc 




    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flu3.1')
    END IF

    S_bu     = 0._wp
    fl_m     = 0._wp
    flush_v  = 0._wp
    flush_h  = 0._wp
    S_bu(1:N_active) = S_abs(1:N_active)/m(1:N_active)

    !determine length of horizontal flow
    const = SUM(thick(1:N_active))*para_flush_horiz

    !First stabilization,melt_thick can not be larger then the fluid fraction

    melt_thick = MIN(melt_thick,psi_l(1)*thick(1))
    !second stabilization, melt_thick may not be larger then thick_0/3
    melt_thick = MIN(melt_thick,thick_0/3._wp)
    
    if (snow_flush_flag == 1) THEN
       !< Niels, 2017 add: if loop, enhanced the permeability, revise
       !First step calculate delta z_mp
       perm             = 0._wp
       
       !< Niels, 2017 add: psi_g to permeability calculation, improved the results but must be checked
       perm(1:N_active) = 10.0_wp**(-17.0_wp) * (1000.0_wp*ABS(psi_l(1:N_active)+2.*psi_g(1:N_active)))**3.10_wp    
       DO k=1,N_active
          IF (perm(k) == 0._wp) THEN
             perm(k) = 1._wp
          END IF
       END DO
    ELSE IF (snow_flush_flag == 0) THEN
       !First step calculate delta z_mp
       perm             = 1._wp
       perm(1:N_active) = 10.0_wp**(-17.0_wp) * (1000.0_wp*ABS(psi_l(1:N_active)))**3.10_wp
    END IF
    
    !Hydraulic resitivities are comfl_md
    DO k = 1,N_active
       R_v(k)=mu*thick(k)        / MAX(perm(k),0.00000000000000000000001_wp)
       !Really important to include thickness!
       R_h(k)=mu*const/(thick(k) * MAX(perm(k),0.00000000000000000000001_wp))
    END DO
    R(N_active)   = 0._wp
    R(N_active-1) = R_v(N_active-1)

    IF (N_active > 2) THEN
       DO k = N_active-2,1,-1
          R(k) = R(k+1)+R_v(k)
          R(k) = ((R(k))*R_h(k))/(R(k)+R_h(k))
       END DO
    END IF



    !Total amount of brine to move, flush_total
    !flush_total = freeboard/R(1)*grav*dt*func_density(T(1),func_S_br(T(1)))*rho_l !< Niels, 2017 add: commented
    flush_total = (freeboard+melt_thick)/R(1)*grav*dt*func_density(T(1),func_S_br(T(1)))*rho_l  !< Niels, 2017 add: melt thich is on top of the ice and therefore also part of the hydraulic head

    !Stabilization, flush_total can be no larger than melt_thick*rho_l
    flush_total = MIN(flush_total,melt_thick*rho_l)
    melt_err = melt_err + melt_thick-MIN(flush_total/rho_l,melt_thick) !< Niels, 2017 add: melt_err, check how much meltwater vanishes in the line above

    !Total flush is divided into the vertical and horizontal fluxes of each layer
    flush_h(1) = flush_total*(R(2)+R_v(1)) / (R(2)+R_v(1)+R_h(1))
    flush_v(1) = flush_total*R_h(1)        / (R(2)+R_v(1)+R_h(1))
    DO k = 2,N_active-1
       flush_h(k) = flush_v(k-1)* (R(k+1)+R_v(k)) / (R(k+1)+R_v(k)+R_h(k))
       flush_v(k) = flush_v(k-1)* R_h(k)          / (R(k+1)+R_v(k)+R_h(k))
    END DO
    flush_v(N_active) = flush_v(N_active-1)
    flush_h(N_active) = 0._wp

    !If bgc tracers are active fl_brine_bgc records the fluxes to calculate bgc advection
    IF( PRESENT( fl_brine_bgc ) ) THEN
       fl_brine_bgc(1:(N_active-1),N_active) = fl_brine_bgc(1:(N_active-1),N_active) + flush_h(1:(N_active-1))
       fl_brine_bgc(N_active,N_active+1)     = fl_brine_bgc(N_active,N_active+1)     + SUM(flush_h(:))
       DO k = 1,N_active 
          fl_brine_bgc(k,k+1) = fl_brine_bgc(k,k+1) + flush_v(k)
       END DO

    END IF

    !Vertical fluxes are dealt with mass transfer subroutine
    fl_m(1)            = 0.0_wp
    fl_m(2:N_active+1) = -flush_v(1:N_active)

    CALL mass_transfer (Nlayer,N_active,T,H_abs,S_abs,S_bu,T_bottom,S_bu_bottom,fl_m)

    !to balance heat loss
    IF (flush_heat_flag==2) THEN
       H_abs(N_active) = H_abs(N_active) - fl_m(N_active+1)*T(N_active)*c_l
    END IF

    !Mass and thickness of the upper layer is reduced by the brine flushing away.
    m(1)     = m(1)     -flush_total
    thick(1) = thick(1) -flush_total/rho_l



    !Now we deal with horizontal flushing, which transports water directly to the lowest layer
    DO k = 1,N_active-1
       loss_S_abs = flush_h(k)*func_S_br(T(k),S_abs(k)/m(k))
       loss_H_abs = flush_h(k)*T(k)*c_l
       S_abs(k)   = S_abs(k)-loss_S_abs
       H_abs(k)   = H_abs(k)-loss_H_abs

       H_abs(N_active) = H_abs(N_active)+loss_H_abs
       S_abs(N_active) = S_abs(N_active)+loss_S_abs


    END DO
    loss_S_abs     = SUM(flush_h)*S_bu(N_active)
    loss_H_abs     = SUM(flush_h)*T(N_active)*c_l

    !Balances heat loss from horizontal flushing
    IF (flush_heat_flag==2) THEN
       H_abs(N_active) = H_abs(N_active)-loss_H_abs
    END IF
    S_abs(N_active) = S_abs(N_active)-loss_S_abs



    IF (MINVAL(S_abs)<-0.00000000000000000000000001_wp) THEN
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flu3.2')
       END IF
       S_abs(1:N_active) = MAX(S_abs(1:N_active),0._wp)
    END IF

    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flu3.3')
    END IF


    IF (ABS(m(1))<0.000001_wp) THEN
       PRINT*,'negative mass after flushing 3, aborted',m
       STOP 9876
    END IF




  END SUBROUTINE flush3


  !>
  !! An alternative subroutine for calculating flushing.
  !!
  !! Simplified approach.
  !! Melt_thick of top layer is simply removed with brine salinity. 
  !! Salinity of a layer is reduced if the solid fraction is lower than that of the layer above it.
  !! Flushing stops as soon as a layer has a higher solid fraction than the layer below it.
  !! 
  !! 
  !!
  !! @par Revision History
  !! Invented by Philipp Griewank, IMPRS (2012-07-9)
  !!
  SUBROUTINE flush4 (psi_l,thick,T,thick_0,S_abs,H_abs,m,dt,Nlayer,N_active,N_top,N_middle,N_bottom,melt_thick,debug_flag)
    INTEGER,                     INTENT(in)    :: Nlayer,N_middle,N_top,N_bottom,debug_flag,N_active
    REAL(wp), DIMENSION(Nlayer), INTENT(in)    :: psi_l,T
    REAL(wp),                    INTENT(in)    :: dt,thick_0
    REAL(wp), DIMENSION(Nlayer), INTENT(inout) :: S_abs,H_abs,m,thick
    REAL(wp),                    INTENT(inout) :: melt_thick
    REAL(wp), DIMENSION(Nlayer)                :: S_bu
    REAL(wp)                                   :: loss_S_abs,loss_H_abs
    INTEGER                                    :: k

    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flu4.1')
    END IF

    !Setup
    S_bu             = 0._wp
    S_bu(1:N_active) = S_abs(1:N_active)/m(1:N_active)


    !Salt and enthalpy loss of the melt water
    H_abs(1)      = H_abs(1)  -melt_thick*rho_l*c_l*T(1)
    S_abs(1)      = S_abs(1)  -melt_thick*rho_l*func_S_br(T(1),S_bu(1))
    thick(1)      = thick(1)-melt_thick  
    m(1)          = m(1)-melt_thick*rho_l
    melt_thick    = 0.0_wp

    k = 2
    DO WHILE (psi_l(k)>psi_l(k-1) )
       S_abs(k) = para_flush_gamma*S_abs(k)
       k        = k+1
    END DO

    S_abs(1)      = MAX(S_abs(1),0.00_wp)

    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flu4.2')
    END IF


    IF(MINVAL(S_abs)<0.00000_wp) THEN
       PRINT*,'negative salinity after flushing 4, aborted',S_abs
       STOP 9876
    END IF
  END SUBROUTINE flush4
END MODULE mo_flush
