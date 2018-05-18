!>
!! Regulates mass transfers and their results.
!!
!! Ultimately all processes which involve a mass flux should be stored here.
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
!! Begin implementing Expulsion by Philipp Griewank, IMPRS (2010-08-24)
!!
MODULE mo_mass

  USE mo_parameters,       ONLY:wp,rho_l,c_l
  USE mo_thermo_functions, ONLY:func_S_br
  USE mo_functions,        ONLY:func_sat_O2
  IMPLICIT NONE

  PUBLIC :: mass_transfer,expulsion_flux,bgc_advection


CONTAINS

  !>
  !! Calculates the effects of mass transfers on H_abs and S_abs.
  !!
  !! The effects of brine displaced by expulsion, flushing or drainage expansion lead to changes in mass, salt ans enthalpy.
  !! This subroutine calculates the effects on S_abs and H_abs.
  !! A very simple upwind strategy is employed, Brine from below has T and S_br of the lower layer, and brine from above T and S_br of the upper layer.
  !! To avoid negative salinity, the maximum amount of advective salt is the total salt content of the layer.
  !! The amount of mass transfered is calculated in other subroutines.
  !! 
  !! This subroutine was started as a quick and dirty way to simulate the bottom freezing experiment described in Notz 2005 p. 85
  !! IMPORTANT: Before this subroutine expelled brine was removed from the system and its effects were determined in subroutine expulsion.
  !! S_bu must be up to date!
  !!
  !! @par Revision History
  !! Brought to life by Philipp Griewank, IMPRS (2010-08-24)
  !! Modified to work with all processes by Philipp Griewank, IMPRS (2010-11-27)
  SUBROUTINE mass_transfer (Nlayer,N_active,T,H_abs,S_abs,S_bu,T_bottom,S_bu_bottom,fl_m)

    IMPLICIT NONE

    INTEGER,                       INTENT(in)    :: Nlayer, N_active
    REAL(wp),                      INTENT(in)    :: T_bottom,S_bu_bottom
    REAL(wp), DIMENSION(Nlayer),   INTENT(in)    :: T,S_bu
    REAL(wp), DIMENSION(Nlayer),   INTENT(inout) :: H_abs,S_abs
    REAL(wp), DIMENSION(Nlayer+1), INTENT(in)    :: fl_m
    REAL(wp), DIMENSION(Nlayer+1)                :: TT,SS_bu,SS_abs !Same as T and S_bu but expanded to include bottom values
    INTEGER                                      :: k


    TT(1:N_active)      =  T(1:N_active)
    SS_bu(1:N_active)   =  S_bu(1:N_active)
    SS_abs(1:N_active)  =  S_abs(1:N_active)

    TT(N_active+1)      =  T_bottom
    SS_bu(N_active+1)   =  S_bu_bottom
    SS_abs(N_active+1)  =  S_bu_bottom*2000._wp



    DO k = 1,N_active

       IF (fl_m(k+1)>0.) THEN
          H_abs(k) =  H_abs(k) +fl_m(k+1)*TT(k+1)*c_l
          S_abs(k) =  S_abs(k) +MIN(fl_m(k+1)*func_S_br(TT(k+1),SS_bu(k+1)), SS_abs(k+1))
       ELSE IF (fl_m(k+1)<0.) THEN
          H_abs(k) =  H_abs(k) +fl_m(k+1)*TT(k)*c_l
          S_abs(k) =  S_abs(k) +MAX(fl_m(k+1)*func_S_br(TT(k),SS_bu(k)), -S_abs(k))
       END IF

       IF (fl_m(k)>0.) THEN
          H_abs(k) =  H_abs(k) -fl_m(k)*TT(k)*c_l
          S_abs(k) =  S_abs(k) -MIN(fl_m(k)*func_S_br(TT(k),SS_bu(k)), S_abs(k))
       ELSE IF (fl_m(k)<0) THEN
          H_abs(k) =  H_abs(k) -fl_m(k)*TT(k-1)*c_l
          S_abs(k) =  S_abs(k) -MAX(fl_m(k)*func_S_br( TT(k-1),SS_bu(k-1)), -S_abs(k-1))
       END IF


    END DO
  END SUBROUTINE mass_transfer



  !>
  !! Generates the fluxes caused by expulsion.
  !!
  !! Brine displaced by expansion of a freezing mushy layer lead to a mass, enthalpy and salt flux.
  !! This subroutine calculates the amount of brine which moves between the layers caused by V_ex and how the mass in the layers changes.
  !! Vary basic assumptions are made. Brine always moves downward (negative), no horizontal movement are allowed and gas pockets can be filled.
  !! The upper boundary layer is not permeable but the bottom one is.
  !! This subroutine was started as a quick and dirty way to simulate the bottom freezing experiment described in Notz 2005 p. 85
  !!
  !! @par Revision History
  !! Brought to life by Philipp Griewank, IMPRS (2010-08-24)
  !! Simplified by Philipp Griewank, IMPRS (2010-11-27)
  SUBROUTINE expulsion_flux (thick,V_ex,Nlayer,N_active,psi_g,fl_m,m)

    INTEGER,                        INTENT(in)    :: Nlayer, N_active
    REAL(wp), DIMENSION(Nlayer),    INTENT(in)    :: V_ex,thick
    REAL(wp), DIMENSION(Nlayer),    INTENT(inout) :: psi_g,m
    REAL(wp), DIMENSION(Nlayer+1),  INTENT(out)   :: fl_m

    INTEGER::k

    fl_m(1:Nlayer+1)  =  0._wp
    fl_m(2)           = -V_ex(1)*rho_l
    DO k = 2,N_active
       IF (psi_g(k)<0.001) THEN
          fl_m(k+1)   = -V_ex(k)*rho_l+fl_m(k)
       ELSE
          fl_m(k+1)   = -MAX((V_ex(k)-psi_g(k)*thick(k))*rho_l    ,0.0_wp)
          psi_g(k)    =  MAX((psi_g(k)*thick(k)-V_ex(k))/thick(k) ,0.0_wp)
       END IF
    END DO

    DO k = 1,N_active
       m(k) =  m(k) +fl_m(k+1)-fl_m(k)
    END DO

  END SUBROUTINE expulsion_flux

  !>
  !! Calculates how the brine fluxes stored in fl_brine_bgc advect bgc tracers
  !!
  !! A very simple upwind strategy is employed.
  !! To avoid negative tracer densities, the maximum amount of advection is restricted to the current tracer content in a layer divided by three.
  !! Three is chosen as a limit as currently each layer can have a maximum of three flows leaving the layer (to the layer above, the layer below, and the lowest layer).
  !! The advection scheme is likely overly diffusive, but given the limitations we are working with (e.g. changing brine volumes) nothing more sophisticated can be applied easily.
  !! 
  !! For gases it might make sense to limit the brine density to saturation value in advecting brine, to take bubble formation into account. This needs to be specified in bgc_advection, and is a first attempt (both scientifically and code wise) which should be used with caution!
  !!
  !! @par Revision History
  !! Brought to life by Philipp Griewank, IMPRS (2014-02-10)
  SUBROUTINE bgc_advection (Nlayer,N_active,N_bgc,fl_brine_bgc,bgc_abs,psi_l,T,S_abs,m,thick,bgc_bottom)

    IMPLICIT NONE

    INTEGER,                                INTENT(in)    :: Nlayer,N_active,N_bgc
    REAL(wp), DIMENSION(Nlayer),            INTENT(in)    :: psi_l,thick,T,S_abs,m
    REAL(wp), DIMENSION(N_bgc),             INTENT(in)    :: bgc_bottom
    REAL(wp), DIMENSION(Nlayer+1,Nlayer+1), INTENT(in)    :: fl_brine_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc),      INTENT(inout) :: bgc_abs
    REAL(wp), DIMENSION(Nlayer,N_bgc)                     :: bgc_temp
    REAL(wp), DIMENSION(Nlayer,N_bgc)                     :: bgc_br
    REAL(wp), DIMENSION(N_bgc)                            :: flux
    REAL(wp)                                              :: sat_O2
    INTEGER                                      ::k,i,j

    bgc_temp = bgc_abs

    !Calculating brine concentration of advecting brine, saturation limitations must be included manually!
    DO k = 1,N_active
       bgc_br(k,:)  =  bgc_abs(k,:)/(MAX(psi_l(k)*thick(k)*rho_l,0.000000000000001_wp))
       
       !How to limit the first tracer to the oxygen saturation limit
       !sat_O2 = func_sat_O2(T(k),S_abs(k)/m(k))
       !bgc_br(k,1) = MIN(bgc_br(k,1),1.25*sat_O2)
    END DO

    !Internal flows
    DO i = 1,N_active
       DO j = 1,N_active
          DO k = 1,N_bgc
             IF (fl_brine_bgc(i,j)*bgc_br(i,k) > bgc_abs(i,k)/3._wp) THEN
             END IF
             flux(k)      = MIN(fl_brine_bgc(i,j)*bgc_br(i,k),bgc_abs(i,k)/3._wp)
          END DO
          bgc_temp(i,:) = bgc_temp(i,:) - flux(:)
          bgc_temp(j,:) = bgc_temp(j,:) + flux(:)

       END DO
    END DO

    !Flows which leave the domain
    DO i = 1,N_active
       j = N_active+1
       DO k = 1,N_bgc
          flux(k)      = MIN(fl_brine_bgc(i,j)*bgc_br(i,k),bgc_abs(i,k)/3._wp)
       END DO
       bgc_temp(i,:) = bgc_temp(i,:) - flux(:)
    END DO

    !Flows which enter the domain
    i = N_active+1
    DO j = 1,N_active
       DO k = 1,N_bgc
          flux(k)      = fl_brine_bgc(i,j)*bgc_bottom(k)
       END DO
       bgc_temp(j,:) = bgc_temp(j,:) + flux(:)
    END DO

    bgc_abs = bgc_temp 
  END SUBROUTINE bgc_advection



END MODULE mo_mass
