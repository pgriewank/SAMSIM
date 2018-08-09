!>
!! Computes the Salt fluxes caused by gravity drainage.
!!
!!
!!
!!
!! @author Philipp Griewank
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
!
!!
!!
!! @par Revision History
!! Injected with life by Philipp Griewank, IMPRS (<2010-08-27>)
!! 
!!

MODULE mo_grav_drain

  USE mo_parameters
  USE mo_thermo_functions, ONLY:func_S_br
  USE mo_mass
  IMPLICIT NONE

  PUBLIC :: fl_grav_drain,fl_grav_drain_simple

CONTAINS


  !>
  !! Calculates fluxes caused by gravity drainage.
  !!
  !! If the Rayleigh number of a layer is higher then the critical value, brine leaves the layer by a theoretical brine channel.
  !! The discharged brine flows downward through all layers directly into the underlying ocean.
  !! To preserve mass the same amount of water flows upwards through all lower layers. 
  !! In contrast to the downward flux the upward flux is assumed to be in thermal equilibrium thus moving salt and heat to each layer.
  !! The upward flux is a standard upwind advection.
  !! The downward flux of a layer over the timestep is = x*(Ray-Ray_crit)*dt*thick.
  !! _____________________________________
  !! |            ->__|__<-
  !! |______________| v |__________________
  !! |  ^         ->|   |<-         ^
  !! |______________|   |__________________
  !! |  ^ ^       ->|||||<-       ^ ^ 
  !! |______________|vvv|__________________
  !!    ^ ^ ^                   ^ ^ ^
  !! This superb ascii art is supposed to show how the assumed fluxes flow downward through a brine channel and upwards through the layer.
  !! x and r are passed on to enable easy optimization. 
  !! The effect of the upward moving brine is calculated in mass_transfer.
  !!  
  !! IMPORTANT: The height assumptions are special. The bottom of the ice edge is assumed to be at psi_s(N_active)/psi_s_min *thick_0 
  !! 
  !! The first approach assumed that brine drainage occurred between two layers but performed poorly.
  !!
  !! If grav_heat_flag is set to 2 the amount of heat transported out of the ice will be compensated in the lowest layer
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-08-27) \n
  !! Completely revised to assume brine channels by Philipp Griewank , IMPRS (2010-11-05) \n
  !! Mass_transfer is used to advect H and S by Philipp Griewank, IMPRS (2010-11-05) \n
  !! Added condition S_br(k)>S_br(k+1) by Philipp Griewank. IMPRS (2011-04-29) \n
  !! Added harmonic mean for permeability by Philipp Griewank (2014-01-05)

  SUBROUTINE fl_grav_drain (S_br,S_bu,psi_l,psi_s,psi_g,thick,S_abs,H_abs,T,m,dt,Nlayer,N_active,ray, &
       T_bottom,S_bu_bottom,grav_drain,grav_temp,grav_salt,grav_heat_flag,harmonic_flag, fl_brine_bgc)
    INTEGER,                       INTENT(in)    :: Nlayer, N_active,grav_heat_flag,harmonic_flag
    REAL(wp), DIMENSION(Nlayer),   INTENT(in)    :: S_br,S_bu,thick,T,psi_l,psi_s,psi_g
    REAL(wp),                      INTENT(in)    :: dt,S_bu_bottom,T_bottom
    REAL(wp),                      INTENT(inout) :: grav_drain,grav_temp,grav_salt
    REAL(wp), DIMENSION(Nlayer),   INTENT(inout) :: S_abs,H_abs,m
    REAL(wp), DIMENSION(Nlayer-1), INTENT(out)   :: ray                   !< Rayleigh number
    REAL(wp), DIMENSION(N_active)                :: fl_up                 !< Upward brine flux [kg] 1 is the flux from 2 to 1, N_active flux from ocean to N_active
    REAL(wp), DIMENSION(N_active)                :: fl_down               !< Downward brine flux [kg] 1 is the flux from 1 to N_active, N_active is from N_active to the ocean
    REAL(wp), DIMENSION(Nlayer)                  :: perm                  !< Permeability
    REAL(wp), DIMENSION(Nlayer)                  :: harmonic_perm         !< Harmonic mean permeability
    REAL(wp), DIMENSION(Nlayer+1)                :: fl_m  
    REAL(wp)                                     :: flux                  !< Downward brine flux [kg]
    REAL(wp)                                     :: test1
    REAL(wp)                                     :: d_S_br,height,ray_mini
    REAL(wp)                                     :: heat_loss             !< Amount of heat transported from the ice [J]
    INTEGER                                      :: k,kk
    REAL(wp), DIMENSION(Nlayer+1,Nlayer+1), INTENT(inout),OPTIONAL :: fl_brine_bgc 

    ray_mini       = ray_crit!Arises from the optimization of the r
    heat_loss      = 0._wp
    perm           = 0.0_wp
    perm(N_active) = 9999999._wp !lowest layer is considered fully liquid
    ray            = 0._wp
    fl_up          = 0._wp
    fl_down        = 0._wp
    harmonic_perm  = 0.0_wp

    !Permeability is determined
    DO k=1,N_active
       perm(k)=10.0_wp**(-17.0_wp) * (1000.0_wp*ABS(psi_l(k)))**3.10_wp
    END DO

    !Harmonic mean permeability is determined
    IF (harmonic_flag == 2) THEN
       DO k=1,N_active-1
          test1 = minval(perm(k:N_active-1))
          IF (test1 .lt. 10._wp**(-14._wp)) THEN
             harmonic_perm(k)=0.0_wp
          ELSE
             DO kk=k,N_active-1
                harmonic_perm(k)=harmonic_perm(k)+thick(kk)/perm(kk)
             END DO
             !bottom layer is included in a linear fashion
             harmonic_perm(k)=harmonic_perm(k)+(thick(N_active)*psi_s(N_active)/psi_s_min)/perm(N_active)
             harmonic_perm(k)=(SUM(thick(k:N_active-1))+thick(N_active)*psi_s(N_active)/psi_s_min)/harmonic_perm(k)
          END IF
       END DO
    END IF

    !___________Raleigh number is calculated (see Notz 2005)__________________________
    DO k=1,N_active-1
       d_S_br = S_br(k)-S_br(N_active)
       height = SUM(thick((k+1):N_active-1))+thick(N_active)*psi_s(N_active)/psi_s_min !Height adjustment
       IF (harmonic_flag == 1) THEN
          ray(k) = grav*rho_l*bbeta*d_S_br*height* MINVAL(perm(k:N_active))
       ELSE IF (harmonic_flag == 2) THEN
          ray(k) = grav*rho_l*bbeta*d_S_br*height* harmonic_perm(k)
       END IF
       ray(k) = ray(k)/(kappa_l*mu)
       ray(k) = MAX(ray(k),0._wp)
    END DO




    grav_salt = grav_salt + SUM(S_abs(:))


    DO k = 1,N_active-1
       IF (ray(k)>ray_mini .AND. psi_s(k)>0.001_wp .AND. S_abs(k)/m(k)>0.1_wp .AND. S_br(k)>S_br(k+1)) THEN  !if psi_s equal zero then flushing is assumed
          flux     = x_grav*(ray(k)-ray_mini) *dt*thick(k)
          flux     = MIN(flux,psi_l(k)*rho_l*thick(k)) !should limit flux to brine volume of layer
          S_abs(k) = S_abs(k) -flux*S_br(k)
          IF (S_abs(k)<0.0_wp) THEN
             PRINT*,'Negative Salinity due to gravity drainige overdrive in layer',k
             PRINT*,S_abs(k),S_abs(k)+flux*S_br(k),flux,psi_l(k)*rho_l*thick(k)
             STOP 21234
          END IF

          grav_temp  = grav_temp + flux*T(k)

          H_abs(k)   = H_abs(k) -flux*c_l*T(k)
          heat_loss  = heat_loss +flux*c_l*T(k)


          fl_down(k) = flux
          FORALL(kk=k:N_active)
             fl_up(kk) = fl_up(kk)+flux 
          END FORALL

          fl_up(k) = MIN(fl_up(k),psi_l(k)*rho_l*thick(k)) ! should limit flux to brine volume of layer


       END IF
    END DO

    grav_salt  = grav_salt - SUM(S_abs(:))

    !Influence of summed up upward transport resulting from the downward drainage
    fl_m(1)=0.0_wp
    fl_m(2:N_active+1)=fl_up(1:N_active)

    IF( PRESENT( fl_brine_bgc ) ) THEN
       fl_brine_bgc(1:(N_active-1),N_active+1) = fl_brine_bgc(1:(N_active-1),N_active) + fl_down(1:(N_active-1))
       !fl_brine_bgc(N_active,N_active+1)     = fl_brine_bgc(N_active,N_active+1)     + sum(fl_down(:))
       DO k = 1,N_active 
          fl_brine_bgc(k+1,k) = fl_brine_bgc(k+1,k) + fl_up(k)
       END DO

    END IF

    CALL mass_transfer (Nlayer,N_active,T,H_abs,S_abs,S_bu,T_bottom,S_bu_bottom,fl_m)

    grav_drain = grav_drain+fl_m(N_active+1)

    !Heatflux correction
    IF (grav_heat_flag==2) THEN
       H_abs(N_active) = H_abs(N_active) +heat_loss-fl_up(N_active)*c_l*T_bottom       
    END IF


    IF(MINVAL(S_abs)<0.0_wp) THEN
       PRINT*,'negative salinity after gravity drainige, aborted',S_abs
       STOP 1337
    END IF
  END SUBROUTINE fl_grav_drain





  !>
  !! Calculates salinity to imitate the effects gravity drainage.
  !! 
  !! Based on the assumption that super critical Rayleigh numbers are quickly reduced below the critical Rayleigh number.
  !! Proposed as a very simplified parametrisation of gravity drainage.
  !! Includes no fluxes of any kind, instead bulk salinity is simply reduced when ever the Rayleigh number is above the critical values.
  !! The parametrization begins from the bottom layers and moves upward.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2012-01-01)

  SUBROUTINE fl_grav_drain_simple (psi_s,psi_l,thick,S_abs,S_br,Nlayer,N_active,ray, &
       grav_drain,harmonic_flag)
    INTEGER,                       INTENT(in)    :: Nlayer, N_active, harmonic_flag
    REAL(wp), DIMENSION(Nlayer),   INTENT(in)    :: thick,psi_l,psi_s,S_br
    REAL(wp),                      INTENT(inout) :: grav_drain
    REAL(wp), DIMENSION(Nlayer),   INTENT(inout) :: S_abs
    REAL(wp), DIMENSION(Nlayer-1), INTENT(out)   :: ray   !< Rayleigh number
    REAL(wp), DIMENSION(Nlayer)                  :: perm  !< Permeability
    REAL(wp), DIMENSION(Nlayer)                  :: harmonic_perm  !< Harmonic permeability
    REAL(wp)                                     :: d_S_br,height,ray_mini,temp
    INTEGER                                      :: k,kk

    ray_mini       = ray_crit 
    perm           = 0.0_wp
    perm(N_active) = 9999999._wp
    ray            = 0._wp

    !Permeability is determined
    DO k = 1,N_active
       perm(k)=10.0_wp**(-17.0_wp) * (1000.0_wp*ABS(psi_l(k)))**3.10_wp
    END DO

    !Harmonic mean permeability is determined
    IF (harmonic_flag == 2) THEN
       DO k=1,N_active-1
          temp = minval(perm(k:N_active-1))
          IF (temp .lt. 10._wp**(-14._wp)) THEN
             harmonic_perm(k)=0.0_wp
          ELSE
             DO kk=k,N_active-1
                harmonic_perm(k)=harmonic_perm(k)+thick(kk)/perm(kk)
             END DO
             !bottom layer is included in a linear fashion
             harmonic_perm(k)=harmonic_perm(k)+(thick(N_active)*psi_s(N_active)/psi_s_min)/perm(N_active)
             harmonic_perm(k)=(SUM(thick(k:N_active-1))+thick(N_active)*psi_s(N_active)/psi_s_min)/harmonic_perm(k)
          END IF
       END DO
    END IF

    !___________Raleigh number is calculated (see Notz 2005)__________________________
    DO k=1,N_active-1
       d_S_br = S_br(k)-S_br(N_active)
       height = SUM(thick((k+1):N_active-1))+thick(N_active)*psi_s(N_active)/psi_s_min !Height adjustment
       IF (harmonic_flag == 1) THEN
          ray(k) = grav*rho_l*bbeta*d_S_br*height* MINVAL(perm(k:N_active))
       ELSE IF (harmonic_flag == 2) THEN
          ray(k) = grav*rho_l*bbeta*d_S_br*height* harmonic_perm(k)
       END IF
       ray(k) = ray(k)/(kappa_l*mu)
       ray(k) = MAX(ray(k),0._wp)
    END DO


    !##########Where the magic happens#################################
    DO k = N_active-1,1,-1
       IF (ray(k) > ray_mini ) THEN  
          S_abs(k)      = S_abs(k)*0.99
       END IF
    END DO
    grav_drain = 0._wp
  END SUBROUTINE fl_grav_drain_simple

END MODULE mo_grav_drain

