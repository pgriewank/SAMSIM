!>
!! Computes the fluxes caused by liquid flooding the snow layer.
!!
!! Water floods the snow layer instantly transforming it to ice which is added to the top layer.
!! As long as the negative freeboard is smaller then a certain parameter (neg_free) the flood strength is limited by the harmonic mean permeability of the whole ice layer driven by the freeboard.  
!! When this parameter is exceed, instant flooding is assumed. 
!! Based on Ted Maksyms work, brine is moved from the ocean to the snow without interacting with the ice in between.
!! Very little of the process is well understood, so this parametrisation is
!! ID mostly speculation.
!! Ratio_flood is a very important parameter, as it regulates how much wicking into the snow layer occurs during melting which dilutes the flooded snow. 
!! Ratio of two should lead  to the snow pack being reduced twice as much as the top layer grows.
!! 
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
!!
!!
!! @par Revision History
!! Copy and pasted into existence by Philipp Griewank, IMPRS (2011-01-21)
MODULE mo_flood

  USE mo_parameters
  USE mo_thermo_functions, ONLY: func_S_br
  USE mo_functions
  USE mo_output,           ONLY: output_raw_lay
  USE mo_mass,             ONLY: mass_transfer
  IMPLICIT NONE

  PRIVATE
  PUBLIC::flood, flood_simple


CONTAINS

  !>
  !! Subroutine for calculating flooding.
  !!
  !! Details explained in module description.
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, IMPRS (2011-01-21)
  !! Cleaned and commented by Philipp Griewank, (2014-04-19)
  !!
  SUBROUTINE flood (freeboard,psi_s,psi_l,S_abs,H_abs,m,T,thick,dt,Nlayer,N_active,T_bottom,S_bu_bottom,&
       &H_abs_snow,m_snow,thick_snow,psi_g_snow,debug_flag,fl_brine_bgc)
    REAL(wp),                    INTENT(inout) :: H_abs_snow,m_snow,thick_snow
    INTEGER,                     INTENT(in)    :: N_active,Nlayer,debug_flag
    REAL(wp), DIMENSION(Nlayer), INTENT(in)    :: psi_s,psi_l,T
    REAL(wp),                    INTENT(in)    :: dt,T_bottom,S_bu_bottom,freeboard,psi_g_snow
    REAL(wp), DIMENSION(Nlayer), INTENT(inout) :: S_abs,H_abs,m,thick
    REAL(wp)                                   :: flood_brine    !< mass of flooded brine [kg], is the same for all layers, is directed upwards
    REAL(wp)                                   :: shift_ice      !< distance that the top ice layer grows
    REAL(wp)                                   :: shift_snow     !< distance that the snow layer shrinks
    REAL(wp), DIMENSION(Nlayer)                :: perm,S_bu
    REAL(wp)                                   :: shift
    REAL(wp)                                   :: harmonic_perm  !< Permeability of total ice layer calculated using the harmonic mean 
    
    REAL(wp), DIMENSION(Nlayer+1,Nlayer+1), INTENT(inout),OPTIONAL :: fl_brine_bgc 
    INTEGER                                    :: k


    perm(1:N_active) =10.0_wp**(-17.0_wp) * (1000.0_wp*psi_l(1:N_active))**3.10_wp
    harmonic_perm    =0._wp


       DO k=1,N_active-1
            harmonic_perm=harmonic_perm+thick(k)/perm(k)
       END DO
       !bottom layer is included in a linear fashion
       harmonic_perm=harmonic_perm+(thick(N_active)*psi_s(N_active)/psi_s_min)/perm(N_active)
       harmonic_perm=(SUM(thick(1:N_active-1))+thick(N_active)*psi_s(N_active)/psi_s_min)/harmonic_perm


    flood_brine = -dt*grav*rho_l*rho_l*harmonic_perm*(freeboard)/(mu*SUM(thick(1:N_active)))



    shift_ice  = flood_brine/(rho_l*psi_g_snow/ratio_flood)
    shift_snow = shift_ice*(1+psi_g_snow/(1._wp-psi_g_snow)*(1._wp-1._wp/ratio_flood))


    DO k = 1,N_active
       S_bu(k) = S_abs(k)/m(k)
    END DO

    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flood1')
    END IF


    S_abs(1)    = S_abs(1) +flood_brine*S_bu(N_active)
    H_abs(1)    = H_abs(1) +flood_brine*H_abs(N_active)/m(N_active)
    m(1)        = m(1)     +flood_brine

    !Top layer grows with flood water by amount shift
    thick(1)   = thick(1)   +shift_ice
    H_abs(1)   = H_abs(1)   +shift_snow/thick_snow*H_abs_snow
    H_abs_snow = H_abs_snow -shift_snow/thick_snow*H_abs_snow
    m(1)       = m(1)       +shift_snow/thick_snow*m_snow
    m_snow     = m_snow     -shift_snow/thick_snow*m_snow
    thick_snow = thick_snow -shift_snow



    !If after the shift the freeboard is still greater then neg_free, additional water from the lowest layer is transported to the snow layer making the top ice layer grow. 
    IF (freeboard+shift_ice<neg_free ) THEN
       shift       = neg_free-(freeboard+shift_ice)
       flood_brine = shift*(psi_g_snow)*rho_l

       !first the flood_brine effects on the top and lowest layer is calculated
       S_abs(N_active) = S_abs(N_active)+(S_bu_bottom-S_bu(N_active))*flood_brine
       H_abs(N_active) = H_abs(N_active)+(T_bottom-T(N_active))*c_l*flood_brine

       S_abs(1) = S_abs(1) +S_bu(N_active)*flood_brine
       H_abs(1) = H_abs(1) +T(N_active)*c_l*flood_brine
       m(1)     = m(1)     +flood_brine


       !Second we have the transformation of snow to snow-ice
       thick(1)   = thick(1)   +shift
       H_abs(1)   = H_abs(1)   +shift/thick_snow*H_abs_snow
       H_abs_snow = H_abs_snow -shift/thick_snow*H_abs_snow
       m(1)       = m(1)       +shift/thick_snow*m_snow
       m_snow     = m_snow     -shift/thick_snow*m_snow
       thick_snow = thick_snow -shift

    END IF

    IF( PRESENT( fl_brine_bgc ) ) THEN
       fl_brine_bgc(N_active,1)              = fl_brine_bgc(N_active,1)     + flood_brine
       fl_brine_bgc(N_active+1,N_active)     = fl_brine_bgc(N_active+1,N_active)     + flood_brine

    END IF


    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flood2')
    END IF

  END SUBROUTINE flood





  !>
  !! Subroutine for calculating flooding.
  !!
  !! Simplified version of flood.
  !! Flooding occurs instantly to fill the negative freeboard until it reaches neg_free with underlying ocean water.
  !!
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, IMPRS (2012-07-16)
  !! Added neg_free limitation. 
  SUBROUTINE flood_simple (freeboard,S_abs,H_abs,m,thick,T_bottom,S_bu_bottom,H_abs_snow,&
       &m_snow,thick_snow,psi_g_snow,Nlayer,N_active,debug_flag)
    REAL(wp),                    INTENT(inout) :: H_abs_snow,m_snow,thick_snow
    INTEGER,                     INTENT(in)    :: Nlayer,debug_flag,N_active
    REAL(wp),                    INTENT(in)    :: T_bottom,S_bu_bottom,freeboard,psi_g_snow
    REAL(wp), DIMENSION(Nlayer), INTENT(inout) :: S_abs,H_abs,m,thick
    REAL(wp)                                   :: flood_brine      !< mass of flooded brine [kg]
    REAL(wp)                                   :: shift            !< thickness of flooded snow [m]



    shift = freeboard-neg_free

    flood_brine = -shift*psi_g_snow*rho_l



    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flooS1')
    END IF



    !Top layer grows with flood water by freeboard
    thick(1) = thick(1) -shift
    S_abs(1) = S_abs(1) +S_bu_bottom*flood_brine
    H_abs(1) = H_abs(1) -shift/thick_snow*H_abs_snow
    H_abs(1) = H_abs(1) +T_bottom*c_l*flood_brine
    m(1)     = m(1)     -shift/thick_snow*m_snow
    m(1)     = m(1)     +flood_brine


    H_abs_snow = H_abs_snow +shift/thick_snow*H_abs_snow
    m_snow     = m_snow     +shift/thick_snow*m_snow
    thick_snow = thick_snow +shift




    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'flooS2')
    END IF

  END SUBROUTINE flood_simple

END MODULE mo_flood
