!>
!! Mo_layer_dynamics contains all subroutines for the growth and shrinking of layer thickness.
!!
!! The middle layers have flexible thickness in contrast to the lower and upper layers which have static thickness.
!! The details are provided in the separate subroutines.
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
!! @par Revision History
!! Shrinking and growth at the bottom are started by Philipp Griewank, IMPRS (2010-07-28) \n
!! add melt_thick_output by Niels Fuchs, MPIMET (2017-03-01)
!!
MODULE mo_layer_dynamics

  USE mo_parameters
  USE mo_output,   ONLY:output_raw_lay


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: layer_dynamics,top_melt,top_grow



CONTAINS


  !>
  !! Organizes the Semi-Adaptive grid SAMSIM uses.
  !!
  !! Modifies the grid and all core variables due to growth or melt.
  !! Calls the different subroutines according to current conditions.
  !! All subroutines can be called with or without biogeochemical tracers active, which is triggered by providing bgc_abs when calling the subroutine.
  !! See Griewank PhD thesis for a full description of the grid.
  !!
  !! Conditions under which following layer dynamics subroutines are called:
  !! - bottom_melt:          lowest layer is ice free, second lowest layer has a solid fraction smaller than phi_s_min/2, and all Nlayer layers are active.
  !! - bottom_melt_simple:   lowest layer is ice free, second lowest layer has a solid fraction smaller than phi_s_min/2, and not all Nlayer layers are active.
  !! - bottom_melt_simple:   lowest layer is ice free, second lowest layer has a solid fraction smaller than phi_s_min/2, all Nlayer layers are active, and the thickness of the middle layers equals thick_0
  !! - bottom_growth_simple: lowest layer has a solid fraction higher than psi_s_min, and not all Nlayer layers are active
  !! - bottom_growth:        lowest layer has a solid fraction higher than psi_s_min, and all Nlayer layers are active
  !! - top_grow:             top layer thicker than 3/2 * thick_0 
  !! - top_melt:             top layer thinner than 1/2 * thick_0 
  !!
  !! If debug_flag is set to 2 the layer values will be written into the debug output (thermoXX.dat) before and after layer dynamics with a string to identify which subroutine was called 
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-07-29) \n
  !! first complete and hopefully stable version by Philipp Griewank, IMPRS (2010-08-10)

  SUBROUTINE layer_dynamics (phi,N_active,Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,thick_0,T_bottom,S_bu_bottom,&
       & bottom_flag,debug_flag,melt_thick_output,N_bgc,bgc_abs,bgc_bottom)

    INTEGER,                                     INTENT(in)    :: Nlayer,N_bottom,N_middle,N_top,bottom_flag,debug_flag
    INTEGER,                                     INTENT(inout) :: N_active
    REAL(wp), DIMENSION(Nlayer),                 INTENT(in)    :: phi
    REAL(wp),                                    INTENT(in)    :: T_bottom,S_bu_bottom,thick_0
    REAL(wp), DIMENSION(Nlayer),                 INTENT(inout) :: m,S_abs,H_abs,thick
    INTEGER,                                     INTENT(in)    :: N_bgc
    REAL(wp),                                    INTENT(inout) :: melt_thick_output    !< Niels, 2017 add: melt_thick_output !OBS: only 3rd element in standard melt_thick_output vector!
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    REAL(wp), DIMENSION(N_bgc),        INTENT(in),    OPTIONAL :: bgc_bottom


    IF (debug_flag==2) THEN
       CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'layer_')
    END IF

    !_____________________________________________________________________________________________________________
    !___Bottom melt_______________________________________________________________________________________________
    !_____________________________________________________________________________________________________________
    IF (phi(Nlayer-1)<=psi_s_min/2._wp .AND. phi(N_active)<0.00001_wp  .AND. N_active==Nlayer & 
         &                              .AND. thick(N_top+1)/thick_0>1.000001_wp .AND. bottom_flag==1) THEN
       IF ( PRESENT(bgc_abs)) THEN
          CALL bottom_melt(Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,N_bgc,bgc_abs)
       ELSE
          CALL bottom_melt(Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'BoMelt')
       END IF
    ELSE IF ( N_active > 1  .AND. N_active<Nlayer .AND. phi(N_active)<0.00001_wp &
         &                   .AND. phi(MAX(N_active-1,1))<=psi_s_min/2._wp .AND. bottom_flag==1) THEN
       IF ( PRESENT(bgc_abs)) THEN
          CALL bottom_melt_simple(N_active,Nlayer,m,S_abs,H_abs,thick,N_bgc,bgc_abs)
       ELSE
          CALL bottom_melt_simple(N_active,Nlayer,m,S_abs,H_abs,thick,N_bgc)
       END IF

       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'BoMeS1')
       END IF
    ELSE IF ( N_active > 1  .AND. phi(N_active)<0.00001_wp  .AND. phi(MAX(N_active-1,1))<=psi_s_min/2._wp &
         .AND. (thick(N_top+1)/thick_0)<1.01_wp .AND. bottom_flag==1 ) THEN
       IF ( PRESENT(bgc_abs)) THEN
          CALL bottom_melt_simple(N_active,Nlayer,m,S_abs,H_abs,thick,N_bgc,bgc_abs)
       ELSE
          CALL bottom_melt_simple(N_active,Nlayer,m,S_abs,H_abs,thick,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'BoMeS2')
       END IF


    !_____________________________________________________________________________________________________________
    !___Bottom growth_____________________________________________________________________________________________
    !_____________________________________________________________________________________________________________

    ELSE IF (phi(N_active)>psi_s_min .AND. N_active<Nlayer .AND. bottom_flag==1) THEN
       IF ( PRESENT(bgc_abs)) THEN
          CALL bottom_growth_simple(N_active,Nlayer,N_top,m,S_abs,H_abs,thick,thick_0,T_bottom,S_bu_bottom,N_bgc,bgc_abs,bgc_bottom)
       ELSE
          CALL bottom_growth_simple(N_active,Nlayer,N_top,m,S_abs,H_abs,thick,thick_0,T_bottom,S_bu_bottom,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'BoGrSi')
       END IF

    ELSE IF (phi(Nlayer)>psi_s_min .AND. bottom_flag==1)  THEN
       IF ( PRESENT(bgc_abs)) THEN
          CALL bottom_growth (Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,T_bottom,S_bu_bottom,N_bgc,bgc_abs,bgc_bottom)
       ELSE
          CALL bottom_growth (Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,T_bottom,S_bu_bottom,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'BoGrow')
       END IF

    !_____________________________________________________________________________________________________________
    !___Top growth________________________________________________________________________________________________
    !_____________________________________________________________________________________________________________
    ELSE IF (thick(1)>1.5_wp*thick_0) THEN
       melt_thick_output = melt_thick_output - thick(1) !< Niels, 2017 add: subtract top growth from melt thick output
       IF ( PRESENT(bgc_abs)) THEN
          CALL top_grow(Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc,bgc_abs)
       ELSE
          CALL top_grow(Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'ToGro2')
       END IF
       melt_thick_output = melt_thick_output + thick(1)

    !_____________________________________________________________________________________________________________
    !___Top melt__________________________________________________________________________________________________
    !_____________________________________________________________________________________________________________
    ELSE IF (thick(1)<0.5_wp*thick_0) THEN
       melt_thick_output = melt_thick_output - thick(1)
       IF ( PRESENT(bgc_abs)) THEN
          CALL top_melt(Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc,bgc_abs)
       ELSE
          CALL top_melt(Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc)
       END IF
       IF (debug_flag==2) THEN
          CALL output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,'ToMel2')
       END IF
       melt_thick_output = melt_thick_output + thick(1)


    END IF

  END SUBROUTINE layer_dynamics


  !!>
  !! Controls melting of the surface layer.
  !!
  !! Subroutine is called when thick(1)<0.5*thick_0(1).
  !! If N_active<Nlayer then N_active = N_active-1.
  !! If N_active=Nlayer and middle layer thickness = thick_0 than N_active = N_active-1.
  !! Otherwise if N_active=Nlayer the middle layers are shrunk.
  !!
  !! @par Revision History
  !! Pasted by Philipp Griewank, IMPRS (2011-05-10)
  !!


  SUBROUTINE top_melt (Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc,bgc_abs)

    INTEGER,                           INTENT(in)              :: Nlayer,N_bottom,N_middle,N_top
    INTEGER,                           INTENT(inout)           :: N_active
    REAL(wp),                          INTENT(in)              :: thick_0
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    INTEGER                                                    :: k
    REAL(wp)                                                   :: loss_m,loss_S_abs,loss_H_abs
    REAL(wp)                                                   :: shift     !> distance the border between layer moves
    REAL(wp), DIMENSION(Nlayer)                                :: rho       !> density of layer
    REAL(wp), DIMENSION(Nlayer)                                :: H,S_bu
    REAL(wp), DIMENSION(N_bgc)                                 :: loss_bgc      
    REAL(wp), DIMENSION(Nlayer,N_bgc)                          :: bgc_temp, bgc_bulk

    
    IF (PRESENT(bgc_abs)) THEN
       bgc_temp = bgc_abs
    ELSE
       bgc_temp = 0._wp
    END IF

    !Calculating rho, H and S_bu 
    DO k = 1,N_active
       rho(k)        = m(k)/thick(k)
       S_bu(k)       = S_abs(k)/m(k)
       H(k)          = H_abs(k)/m(k)
       bgc_bulk(k,:) = bgc_temp(k,:)/m(k)
    END DO

    !Adjusting the top layer from thick(1)<0.5*thick_0(1) to thick(1)=thick(1)+thick_0(1)
    loss_m      = thick_0    *rho(1)
    loss_S_abs  = loss_m     *S_bu(1)
    loss_H_abs  = loss_m     *H(1)
    loss_bgc(:) = loss_m     *bgc_bulk(1,:)

    m(1)          = m(1)     +m(2)
    S_abs(1)      = S_abs(1) +S_abs(2)
    H_abs(1)      = H_abs(1) +H_abs(2)
    bgc_temp(1,:) = bgc_temp(1,:) +bgc_temp(2,:)
    thick(1)      = thick(1) +thick(2)

    !Adjusting top layers
    DO k = 2,MIN(N_top-1,N_active-1)
       m(k)          = rho (k+1)          *thick_0
       S_abs(k)      = S_bu(k+1)          *rho(k+1)*thick_0
       H_abs(k)      = H (k+1)            *rho(k+1)*thick_0
       bgc_temp(k,:) = bgc_bulk (k+1,:)   *rho(k+1)*thick_0
    END DO


    !Removing bottom layer if N_active<N_top
    IF (N_active<=N_top) THEN
       m(N_active)          = 0.0_wp 
       S_abs(N_active)      = 0.0_wp  
       H_abs(N_active)      = 0.0_wp 
       thick(N_active)      = 0.0_wp
       bgc_temp(N_active,:) = 0.0_wp 

       N_active = N_active-1

    ELSE IF (N_active>N_top .AND. N_active<=Nlayer .AND. thick(N_top+1)/thick_0<1.00001_wp) THEN
       !Removing bottom layer and moving all other layers if N_active<Nlayer .and. N_active>N_top.
       !Is also activated if thick(N_top+1)=thick_0
       DO k = N_top,N_active-1
          m(k)          = rho (k+1)        *thick_0
          S_abs(k)      = S_bu(k+1)        *rho(k+1)*thick_0
          H_abs(k)      = H (k+1)          *rho(k+1)*thick_0
          bgc_temp(k,:) = bgc_bulk (k+1,:) *rho(k+1)*thick_0
       END DO
       m(N_active)          = 0.0_wp 
       S_abs(N_active)      = 0.0_wp  
       H_abs(N_active)      = 0.0_wp 
       bgc_temp(N_active,:) = 0.0_wp 
       thick(N_active)      = 0.0_wp
       

       N_active = N_active-1
    END IF

    IF (N_active==Nlayer .AND. thick(N_top+1)-thick_0>=0.000001_wp ) THEN
       !Middle layers are adjusted__________________________
       loss_m     = thick_0*rho(N_top+1)
       loss_S_abs = loss_m*S_bu(N_top+1)
       loss_H_abs = loss_m*H(N_top+1)
       loss_bgc   = loss_m*bgc_bulk(N_top+1,:)

       !layer N_top needs new values________________________
       m(N_top)          = loss_m
       S_abs(N_top)      = loss_S_abs
       H_abs(N_top)      = loss_H_abs
       bgc_temp(N_top,:) = loss_bgc(:)


       DO k=N_top+1,N_middle+N_top
          m(k)     = m(k)     -loss_m
          H_abs(k) = H_abs(k) -loss_H_abs
          S_abs(k) = S_abs(k) -loss_S_abs
          bgc_temp(k,:) = bgc_temp(k,:) -loss_bgc(:)

          shift = thick_0*REAL(N_middle-k+N_top)/REAL(N_middle) !distance the border moves

          loss_m     = shift  *rho(k+1)
          loss_S_abs = loss_m *S_bu(k+1)
          loss_H_abs = loss_m *H(k+1)
          loss_bgc   = loss_m *bgc_bulk(k+1,:)

          m(k)     = m(k)     +loss_m
          H_abs(k) = H_abs(k) +loss_H_abs
          S_abs(k) = S_abs(k) +loss_S_abs
          bgc_temp(k,:) = bgc_temp(k,:) +loss_bgc(:)
       END DO





       !Middle layer thickness is adjusted________________________
       DO k = N_top+1,N_top+N_middle
          thick(k) = thick(k)-thick_0/REAL(N_middle)
       END DO
    END IF

 
    !Checks if current thick vector is possible for current N_active and thick_0
    !Uses 0.501 instead of 0.500 to give a bit of wriggle room for numerical noise
    IF(thick_0*(N_active+0.501_wp)<=SUM(thick) .AND. N_active<Nlayer) THEN
       PRINT*,'wtf layer issue',thick_0,SUM(thick),N_active,Nlayer
       STOP 7889
    END IF
    
    IF (PRESENT(bgc_abs)) THEN
       bgc_abs = bgc_temp
    END IF
  END SUBROUTINE top_melt


  !>
  !! Controls melting at the bottom.
  !!
  !! Subroutine is called if the second lowest layer is fully liquid and the middle layers are thicker then the bottom layers.
  !! When this occurs the thickness of the middle layers is reduced by 1/N_middle*thick_0.
  !! The lower layers are adjusted upwards, layer(k)=layer(k-1).
  !!
  !!
  !! @par Revision History
  !! Started by Philipp Griewank, IMPRS (2010-07-28) \n
  !! Linear profiles replaced with simple upstream by Philipp Griewank, IMPRS (2010-08-25)
  !!
  SUBROUTINE bottom_melt (Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,N_bgc,bgc_abs)

    INTEGER,                           INTENT(in)              :: Nlayer,N_bottom,N_middle,N_top
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    INTEGER                                                    :: k
    REAL(wp)                                                   :: loss_m,loss_S_abs,loss_H_abs,shift
    REAL(wp), DIMENSION(Nlayer)                                :: rho       !>density of layer
    REAL(wp), DIMENSION(Nlayer)                                :: H,S_bu
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    REAL(wp), DIMENSION(N_bgc)                                 :: loss_bgc      
    REAL(wp), DIMENSION(Nlayer,N_bgc)                          :: bgc_temp, bgc_bulk


    IF (PRESENT(bgc_abs)) THEN
       bgc_temp = bgc_abs
    ELSE
       bgc_temp = 0._wp
    END IF
    

    !Calculating rho, H and S_bu 
    DO k = N_top+1,Nlayer
       rho(k)        = m(k)/thick(k)
       S_bu(k)       = S_abs(k)/m(k)
       H(k)          = H_abs(k)/m(k)
       bgc_bulk(k,:) = bgc_temp(k,:)/m(k)
 
    END DO

    loss_m     = 0._wp
    loss_S_abs = 0._wp
    loss_H_abs = 0._wp
    loss_bgc   = 0._wp

    !_________________________Middle layers are adjusted__________________________
    DO k = N_top+1,N_top+N_middle

       !each layer receives all that the layer above loses
       m(k)          = m(k)          +loss_m
       H_abs(k)      = H_abs(k)      +loss_H_abs
       S_abs(k)      = S_abs(k)      +loss_S_abs
       bgc_temp(k,:) = bgc_temp(k,:) +loss_bgc(:)

       

       shift    = thick(Nlayer)*(k-N_top)/ REAL(N_middle) !distance the border moves

       loss_m         = shift*rho(k)
       loss_H_abs     = loss_m*H(k)
       loss_S_abs     = loss_m*S_bu(k)
       loss_bgc(:)    = loss_m*bgc_bulk(k,:)


       m(k)          = m(k)          -loss_m
       H_abs(k)      = H_abs(k)      -loss_H_abs
       S_abs(k)      = S_abs(k)      -loss_S_abs
       bgc_temp(k,:) = bgc_temp(k,:) -loss_bgc(:)


    END DO


    !________________________Middle layer thickness is adjusted

    DO k = N_top+1,N_top+N_middle
       thick(k) = thick(k)-thick(Nlayer)/REAL(N_middle)
    END DO
    !_________________________Bottom layers are adjusted__________________________

    DO k = N_top+N_middle+1,Nlayer
       H_abs(k)      = rho(k-1)*thick(k)*H(k-1)
       S_abs(k)      = rho(k-1)*thick(k)*S_bu(k-1)
       m(k)          = rho(k-1)*thick(k)
       bgc_temp(k,:) = rho(k-1)*thick(k)*bgc_bulk(k-1,:)
    END DO

    IF (PRESENT(bgc_abs)) THEN
       bgc_abs = bgc_temp
    END IF

  END SUBROUTINE bottom_melt



  !>
  !! Controls bottom growth.
  !!
  !! Subroutine is called if the lowest layer is not fully liquid and the middle layers are thicker then thick_0.
  !! When this occurs the thickness of the middle layers is increased by 1/N_middle*thick(bottom).
  !! The lower layers are adjusted downwards, layer(k)=layer(k+1).
  !!
  !!
  !! @par Revision History
  !! Started by Philipp Griewank, IMPRS (2010-07-29) \n
  !! Linear profiles removed, replaced with simple upwind by Philipp Griewank, IMPRS (2010-08-25)
  !!
  SUBROUTINE bottom_growth (Nlayer,N_bottom,N_middle,N_top,m,S_abs,H_abs,thick,T_bottom,S_bu_bottom,N_bgc,bgc_abs,bgc_bottom)

    INTEGER,                           INTENT(in)              :: Nlayer,N_bottom,N_middle,N_top
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    REAL(wp),                          INTENT(in)              :: T_bottom,S_bu_bottom
    INTEGER                                                    :: k
    REAL(wp)                                                   :: gain_m,gain_S_abs,gain_H_abs,shift
    REAL(wp), DIMENSION(Nlayer)                                :: rho       !>density of layer
    REAL(wp), DIMENSION(Nlayer)                                :: H,S_bu
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    REAL(wp), DIMENSION(N_bgc),        INTENT(in),    OPTIONAL :: bgc_bottom      
    REAL(wp), DIMENSION(N_bgc)                                 :: gain_bgc      
    REAL(wp), DIMENSION(Nlayer,N_bgc)                          :: bgc_temp, bgc_bulk


    IF (PRESENT(bgc_abs)) THEN
       bgc_temp = bgc_abs
    ELSE
       bgc_temp = 0._wp
    END IF


    !Calculating rho, H and S_bu 
    DO k = N_top+1,N_top+N_middle+1
       rho(k)        = m(k)/thick(k)
       S_bu(k)       = S_abs(k)/m(k)
       H(k)          = H_abs(k)/m(k)
       bgc_bulk(k,:) = bgc_temp(k,:)/m(k)
    END DO

    gain_m     = 0._wp
    gain_S_abs = 0._wp
    gain_H_abs = 0._wp
    gain_bgc   = 0._wp

    !_________________________Middle layers are adjusted__________________________
    DO k = N_top+1,N_top+N_middle

       !each layer loses all that the upper layer gained
       m(k)          = m(k)          -gain_m
       H_abs(k)      = H_abs(k)      -gain_H_abs
       S_abs(k)      = S_abs(k)      -gain_S_abs
       bgc_temp(k,:) = bgc_temp(k,:) -gain_bgc(:)

       shift = thick(Nlayer)*(k-N_top)/ REAL(N_middle) !distance the border moves

       gain_m     = shift*rho(k+1)
       gain_H_abs = gain_m*H(k+1)
       gain_S_abs = gain_m*S_bu(k+1)
       gain_bgc(:)= gain_m*bgc_bulk(k+1,:)

       m(k)          = m(k)          +gain_m
       H_abs(k)      = H_abs(k)      +gain_H_abs
       S_abs(k)      = S_abs(k)      +gain_S_abs
       bgc_temp(k,:) = bgc_temp(k,:) +gain_bgc(:)

    END DO

    !________________________Middle layer thickness is adjusted
    DO k=N_top+1,N_top+N_middle
       thick(k) = thick(k)+thick(Nlayer)/REAL(N_middle)
    END DO

    !_________________________Bottom layers are adjusted__________________________
    DO k = Nlayer-N_bottom+1,Nlayer-1
       H_abs(k)       = H_abs(k+1)
       S_abs(k)       = S_abs(k+1)
       m(k)           = m(k+1)
       bgc_temp(k,:)  = bgc_temp(k+1,:)
    END DO

    !_________________________Lowest layer is set using T_bottom and S_bu_bottom
    m(Nlayer)     = thick(Nlayer)*rho_l
    H_abs(Nlayer) = m(Nlayer)*T_bottom*c_l
    S_abs(Nlayer) = m(Nlayer)*S_bu_bottom

    IF (PRESENT(bgc_abs)) THEN
       bgc_temp(Nlayer,:) = m(Nlayer)*bgc_bottom(:)
       bgc_abs = bgc_temp
    END IF

  END SUBROUTINE bottom_growth




  !>
  !! Controls growth if the middle layers are not all active.
  !!
  !! Subroutine is called if the lowest layer is not fully liquid and the middle layers are not all activated.
  !! A new layer is activated.
  !! The new lowest layer is assigned standard values of T_bottom and S_bu_bottom.
  !!
  !!
  !! @par Revision History
  !! Started by Philipp Griewank, IMPRS (2010-07-29) \n
  !! Expanded to deal with not full top layers by Philipp Griewank, IMPRS (2011-01-13) \n
  !! Expansion removed and bgc added by Philipp Griewank, IMPRS (2014-02-04) \n
  SUBROUTINE bottom_growth_simple (N_active,Nlayer,N_top,m,S_abs,H_abs,thick,thick_0,T_bottom,S_bu_bottom,N_bgc,bgc_abs,bgc_bottom)

    INTEGER,                           INTENT(in)              :: Nlayer,N_top
    INTEGER,                           INTENT(inout)           :: N_active
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    REAL(wp),                          INTENT(in)              :: thick_0
    REAL(wp),                          INTENT(in)              :: T_bottom,S_bu_bottom
    REAL(wp)                                                   :: shift
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    REAL(wp), DIMENSION(N_bgc),        INTENT(in),    OPTIONAL :: bgc_bottom


    N_active        = N_active+1
    thick(N_active) = thick_0

    m(N_active)     = thick(N_active)*rho_l
    H_abs(N_active) = m(N_active)    *T_bottom*c_l
    S_abs(N_active) = m(N_active)    *S_bu_bottom

    IF (PRESENT(bgc_abs)) THEN
       bgc_abs(N_active,:) = bgc_bottom(:)*m(N_active)
    END IF

  END SUBROUTINE bottom_growth_simple

  !>
  !! Controls bottom melting if the middle layers are not all active.
  !!
  !! Subroutine is called if the second lowest active layer is fully liquid and the middle layers are not all activated, or if the middle layers are thick_0 thick.
  !! The lowest layer is simply deactivated.
  !!
  !!
  !! @par Revision History
  !! Started by Philipp Griewank, IMPRS (2010-07-30>)
  !!
  SUBROUTINE bottom_melt_simple(N_active,Nlayer,m,S_abs,H_abs,thick,N_bgc,bgc_abs)

    INTEGER,                           INTENT(in)              :: Nlayer
    INTEGER,                           INTENT(inout)           :: N_active
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs

    thick(N_active) = 0._wp
    m(N_active)     = 0._wp
    S_abs(N_active) = 0._wp
    H_abs(N_active) = 0._wp
    IF(PRESENT(bgc_abs)) THEN
       bgc_abs(N_active,:) = 0._wp
    END IF
    N_active = N_active-1

  END SUBROUTINE bottom_melt_simple



  !>
  !! Top grow subroutine.
  !!
  !! Should be called when the top layer is thicker then 1.5 *thick_0. 
  !! If N_active=Nlayer middle layers are expanded by thick_0/N_middle and top layers are moved one down.
  !! IF N_active<Nlayer then N_active=N_active+1 and all layers are shifted downwards.
  !!
  !!
  !! @par Revision History
  !! Started by Philipp Griewank, IMPRS (2011-05-10>)
  !!


  SUBROUTINE top_grow (Nlayer,N_active,N_bottom,N_middle,N_top,thick_0,m,S_abs,H_abs,thick,N_bgc,bgc_abs)

    INTEGER,                           INTENT(in)              :: Nlayer,N_bottom,N_middle,N_top
    INTEGER,                           INTENT(inout)           :: N_active
    REAL(wp),                          INTENT(in)              :: thick_0
    REAL(wp), DIMENSION(Nlayer),       INTENT(inout)           :: m,S_abs,H_abs,thick
    INTEGER                                                    :: k
    REAL(wp)                                                   :: loss_m,loss_S_abs,loss_H_abs,shift
    REAL(wp), DIMENSION(Nlayer)                                :: rho       !>density of layer
    REAL(wp), DIMENSION(Nlayer)                                :: H,S_bu
    INTEGER,                           INTENT(in)              :: N_bgc
    REAL(wp), DIMENSION(Nlayer,N_bgc), INTENT(inout), OPTIONAL :: bgc_abs
    REAL(wp), DIMENSION(N_bgc)                                 :: loss_bgc      
    REAL(wp), DIMENSION(Nlayer,N_bgc)                          :: bgc_temp, bgc_bulk
    
    IF (PRESENT(bgc_abs)) THEN
       bgc_temp = bgc_abs
    ELSE
       bgc_temp = 0._wp
    END IF
    !Calculating rho, H and S_bu 
    DO k = 1,N_active
       rho(k)        = m(k)/thick(k)
       S_bu(k)       = S_abs(k)/m(k)
       H(k)          = H_abs(k)/m(k)
       bgc_bulk(k,:) = bgc_temp(k,:)/m(k)
    END DO

    !Adjusting the top layer from thick(1)>1.5*thick_0(1) to thick(1)=thick(1)-thick_0(1)
    loss_m      = thick_0*rho(1)
    loss_S_abs  = loss_m*S_bu(1)
    loss_H_abs  = loss_m*H(1)
    loss_bgc(:) = loss_m*bgc_bulk(1,:)

    m(1)          = m(1)          -loss_m
    S_abs(1)      = S_abs(1)      -loss_S_abs
    H_abs(1)      = H_abs(1)      -loss_H_abs
    bgc_temp(1,:) = bgc_temp(1,:) -loss_bgc(:)
    thick(1)      = thick(1)      -thick_0

    !Adjusting top layers
    DO k = 2,MIN(N_top,N_active)
       m(k)          = rho (k-1)         *thick_0 
       S_abs(k)      = S_bu(k-1)         *rho(k-1)*thick_0 
       H_abs(k)      = H (k-1)           *rho(k-1)*thick_0
       bgc_temp(k,:) = bgc_bulk(k-1,:)   *rho(k-1)*thick_0
    END DO

    !Adding new bottom layer if N_active<=N_top+1
    IF (N_active<=N_top) THEN
       N_active            = N_active+1
       m(N_active)         = rho (N_active-1)         *thick_0
       S_abs(N_active)     = S_bu(N_active-1)         *thick_0*rho(N_active-1) 
       H_abs(N_active)     = H (N_active-1)           *thick_0*rho(N_active-1)
       bgc_temp(N_active,:) = bgc_bulk(N_active-1,:)  *thick_0*rho(N_active-1)
       thick(N_active)     = thick_0

       !Adding new bottom layer and moving all other layers if N_active<Nlayer .and. N_active>N_top
    ELSE IF (N_active>N_top .AND. N_active<Nlayer) THEN
       DO k = N_top+1,N_active
          m(k)         = rho (k-1)         *thick_0
          S_abs(k)     = S_bu(k-1)         *rho(k-1)*thick_0 
          H_abs(k)     = H (k-1)           *rho(k-1)*thick_0 
          bgc_temp(k,:) = bgc_bulk(k-1,:)  *rho(k-1)*thick_0
       END DO
       N_active = N_active+1
       m(N_active)         = rho (N_active-1)         *thick_0 
       S_abs(N_active)     = S_bu(N_active-1)         *thick_0 *rho(N_active-1)
       H_abs(N_active)     = H (N_active-1)           *thick_0 *rho(N_active-1)
       bgc_temp(N_active,:) = bgc_bulk(N_active-1,:)  *thick_0 *rho(N_active-1)
       thick(N_active)     = thick_0

    ELSE IF (N_active==Nlayer) THEN
       !Middle layers are adjusted__________________________
       loss_m      = thick_0*rho(N_top)
       loss_S_abs  = loss_m*S_bu(N_top)
       loss_H_abs  = loss_m*H(N_top)
       loss_bgc(:) = loss_m*bgc_bulk(N_top,:)

       DO k = N_top+1,N_middle+N_top
          m(k)          = m(k)          +loss_m
          H_abs(k)      = H_abs(k)      +loss_H_abs
          S_abs(k)      = S_abs(k)      +loss_S_abs
          bgc_temp(k,:) = bgc_temp(k,:) +loss_bgc(:)

          shift = thick_0*REAL(N_middle-k+N_top)/REAL(N_middle) !distance the border moves

          loss_m      = shift  *rho(k)
          loss_S_abs  = loss_m *S_bu(k)
          loss_H_abs  = loss_m *H(k)
          loss_bgc(:) = loss_m *bgc_bulk(k,:)

          m(k)          = m(k)          -loss_m
          H_abs(k)      = H_abs(k)      -loss_H_abs
          S_abs(k)      = S_abs(k)      -loss_S_abs
          bgc_temp(k,:) = bgc_temp(k,:) -loss_bgc(:)
       END DO

       !Middle layer thickness is adjusted
       DO k = N_top+1,N_top+N_middle
          thick(k) = thick(k)+thick_0/REAL(N_middle)
       END DO
    END IF


    IF (PRESENT(bgc_abs)) THEN
       bgc_abs = bgc_temp
    END IF

  END SUBROUTINE top_grow


END MODULE mo_layer_dynamics
