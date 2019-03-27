!>
!! All things output.
!!
!! Used to clean up root.f90 and make it easier to implement changes to the output. 
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
!! @par Revision History
!! Brought from the womb by Philipp Griewank, IMPRS (2010-10-11) \n
!! add some output values by Niels Fuchs, MPIMET (2017-03-01)
!!
MODULE mo_output

  USE mo_parameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output_begin,output_begin_bgc,output_raw,output_raw_snow,output_raw_lay,output,output_bgc,output_settings!,output_end

CONTAINS
  !>
  !! Settings output.
  !!
  !! Writes important values to latter identify run.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2011-02-12)
  !!
  SUBROUTINE output_settings(description,testcase,N_top,N_bottom,Nlayer,fl_q_bottom,T_bottom,S_bu_bottom,thick_0,time_out,    &
    time_total,dt,boundflux_flag,atmoflux_flag,albedo_flag,grav_flag,flush_flag,flood_flag,grav_heat_flag,flush_heat_flag,    &
    harmonic_flag,prescribe_flag,salt_flag,turb_flag,bottom_flag,tank_flag,precip_flag,bgc_flag,N_bgc,k_snow_flush)
    INTEGER,         INTENT(in) :: testcase,N_top,N_bottom,Nlayer, boundflux_flag,atmoflux_flag,albedo_flag,grav_flag,flush_flag, &
                                   flood_flag,grav_heat_flag , flush_heat_flag,harmonic_flag,prescribe_flag,salt_flag,turb_flag,  &
                                   bottom_flag,tank_flag,precip_flag,bgc_flag,N_bgc
    REAL(wp),        INTENT(in) :: fl_q_bottom,T_bottom,S_bu_bottom,thick_0,time_out,time_total,dt,k_snow_flush
    CHARACTER*12000, INTENT(in) :: description
    
    
    OPEN(1234,file='./output/dat_settings.dat' ,STATUS='replace')
    WRITE(1234,*)'################  Description  ###############'
    WRITE(1234,*) TRIM(description)
    WRITE(1234,*)'#################  Testcase  #################'
    WRITE(1234,'(A16,I9)')    'testcase        =',testcase
       
    !*************************************************************************************************************************
    !Layer settings and allocation
    !*************************************************************************************************************************
    WRITE(1234,*) '##############  Basic settings  ##############'
    WRITE(1234,'(A16,F15.3)') 'dt              =',dt
    WRITE(1234,'(A16,F15.3)') 'thick_0         =',thick_0
    WRITE(1234,'(A16,F15.3)') 'time_out        =',time_out
    WRITE(1234,'(A16,F15.3)') 'time_total      =',time_total
    
    
    WRITE(1234,'(A16,F15.3)') 'fl_q_bottom     =',fl_q_bottom
    WRITE(1234,'(A16,F15.3)') 'T_bottom        =',T_bottom
    WRITE(1234,'(A16,F15.3)') 'S_bu_bottom     =',S_bu_bottom
    
    
    WRITE(1234,'(A16,I9.0)')  'N_top           =',N_top
    WRITE(1234,'(A16,I9.0)')  'N_middle        =',Nlayer-N_top-N_bottom
    WRITE(1234,'(A16,I9.0)')  'N_bottom        =',N_bottom
    WRITE(1234,'(A16,I9.0)')  'Nlayer          =',Nlayer
    
    !*************************************************************************************************************************
    !Flags
    !*************************************************************************************************************************
    WRITE(1234,*) '##################  Flags  ###################'
    !________________________top heat flux____________
    WRITE(1234,'(A16,I9.0)')  'boundflux_flag  =',boundflux_flag
    WRITE(1234,'(A16,I9.0)')  'atmoflux_flag   =',atmoflux_flag
    WRITE(1234,'(A16,I9.0)')  'albedo_flag     =',albedo_flag
    !________________________brine_dynamics____________
    WRITE(1234,'(A16,I9.0)')  'grav_flag       =',grav_flag
    WRITE(1234,'(A16,I9.0)')  'flush_flag      =',flush_flag
    WRITE(1234,'(A16,I9.0)')  'flood_flag      =',flood_flag
    WRITE(1234,'(A16,I9.0)')  'grav_heat_flag  =',grav_heat_flag
    WRITE(1234,'(A16,I9.0)')  'flush_heat_flag =',flush_heat_flag
    WRITE(1234,'(A16,I9.0)')  'harmonic_flag   =',harmonic_flag
    WRITE(1234,'(A16,F15.3)')  'k_snow_flush    =',k_snow_flush !< Niels, 2017
    !________________________Salinity____________
    WRITE(1234,'(A16,I9.0)')  'prescribe_flag  =',prescribe_flag
    WRITE(1234,'(A16,I9.0)')  'salt_flag       =',salt_flag
    !________________________bottom setting______________________
    WRITE(1234,'(A16,I9.0)')  'turb_flag       =',turb_flag
    WRITE(1234,'(A16,I9.0)')  'bottom_flag     =',bottom_flag
    WRITE(1234,'(A16,I9.0)')  'tank_flag       =',tank_flag
    WRITE(1234,'(A16,I9.0)')  'precip_flag     =',precip_flag
    WRITE(1234,'(A16,I9.0)')  'bgc_flag        =',bgc_flag
    WRITE(1234,'(A16,I9.0)')  'N_bgc           =',N_bgc


    CLOSE(1234)
  END SUBROUTINE output_settings

  !>
  !! Standard output.
  !!
  !! For time=n*time_out data is exported.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-10-11)
  !!
  SUBROUTINE output(Nlayer,T,psi_s,psi_l,thick,S_bu,ray,format_T,format_psi, &
       format_thick,format_snow,freeboard,thick_snow,T_snow,psi_l_snow,psi_s_snow,           &
       energy_stored,freshwater,total_resist,thickness,bulk_salin,                &
       grav_drain,grav_salt,grav_temp,T2m,T_top,perm,format_perm,flush_v,flush_h,psi_g,melt_thick_output,format_melt)
    INTEGER,                       INTENT(in) :: Nlayer
    REAL(wp), DIMENSION(Nlayer),   INTENT(in) :: T,psi_s,psi_l,thick,S_bu,perm,flush_v,flush_h,psi_g
    REAL(wp), DIMENSION(Nlayer-1), INTENT(in) :: ray
    REAL(wp),                      INTENT(in) :: freeboard,thick_snow,T_snow,psi_l_snow,psi_s_snow,energy_stored,&
                                                 freshwater,thickness,bulk_salin, &
                                                 total_resist,grav_drain,grav_salt,grav_temp,T2m,T_top
    REAL(wp), DIMENSION(3),        INTENT(in) :: melt_thick_output  !< Niels, 2017: 1: accumulated melt_thick, 2: accumulated melt_thick_snow, 3: accumulated top ice thickness variations (recheck 3: in mo_layer_dynamics)
    CHARACTER*12000,               INTENT(in) :: format_T,format_psi,format_thick,format_snow,format_perm,format_melt

    WRITE(30,format_T)     T
    WRITE(31,format_psi)   psi_s
    WRITE(32,format_thick) thick
    WRITE(33,format_T)     S_bu
    WRITE(34,format_T)     ray
    WRITE(35,format_psi)   psi_l
    WRITE(40,'(F9.3)')     freeboard
    WRITE(41,format_snow)  thick_snow,T_snow,psi_l_snow,psi_s_snow
    WRITE(42,'(F15.1,2x,F10.5,2x,F10.5,2x,F10.5,2x,F10.5)')energy_stored,freshwater,total_resist,thickness,bulk_salin
    WRITE(43,'(F9.6,2X,F9.5,2X,F7.3)')grav_drain,grav_salt,grav_temp
    WRITE(45,*)            T2m,T_top
    WRITE(46,format_perm) perm  !< Niels, 2017 add: output permeability
    WRITE(47,format_perm) flush_v   !< Niels, 2017 add: output vertical flushing
    WRITE(48,format_perm) flush_h   !< Niels, 2017 add: output horizontal flushing
    WRITE(49,format_psi) psi_g  !< Niels, 2017 add: output gas fraction !OBS: not simulated physically in SAMSIM
    WRITE(50,format_melt) melt_thick_output !< Niels, 2017
  
  END SUBROUTINE output

  !>
  !! Standard bgc output.
  !!
  !! For time=n*time_out data is exported.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2014-02-06)
  !!
  SUBROUTINE output_bgc(Nlayer,N_active,bgc_bottom,N_bgc,bgc_abs,psi_l,thick,m,format_bgc)
    INTEGER,                             INTENT(in) :: Nlayer, N_bgc, N_active
    REAL(wp), DIMENSION(N_bgc),          INTENT(in) :: bgc_bottom
    REAL(wp), DIMENSION(Nlayer),         INTENT(in) :: psi_l,m,thick
    REAL(wp), DIMENSION(Nlayer,N_bgc),   INTENT(in) :: bgc_abs
    REAL(wp), DIMENSION(Nlayer)                     :: bgc_bu,bgc_br
    CHARACTER*12000,                     INTENT(in) :: format_bgc
    INTEGER :: k,kk
    
    DO k=1,N_bgc
      DO kk=1,Nlayer
        IF (kk<=N_active) THEN
           IF (m(kk).NE.0._wp) THEN
              bgc_bu(kk) = bgc_abs(kk,k)/m(kk)
              IF (psi_l(kk).NE.0._wp .AND. thick(kk).NE.0) THEN
                 bgc_br(kk) = bgc_abs(kk,k)/psi_l(kk)/thick(kk)/rho_l
              ELSE
                 bgc_br(kk) = 0000._wp
              END IF
           ELSE
              bgc_bu(kk) = 0._wp
              bgc_br(kk) = 0._wp
           END IF
        ELSE
              bgc_bu(kk) = bgc_bottom(k)
              bgc_br(kk) = bgc_bottom(k)
        END IF
      END DO
      WRITE(2*k+400,format_bgc)bgc_bu
      WRITE(2*k+401,format_bgc)bgc_br
    END DO
  
  END SUBROUTINE output_bgc

  !>
  !! Output for debugging purposes.
  !!
  !! Data for each layer is written out each time step to aid in finding errors or understanding model behavior.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-10-11)
  !!
  SUBROUTINE output_raw(Nlayer,N_active,time,T,thick,S_bu,psi_s,psi_l,psi_g)
    INTEGER,                     INTENT(in) :: Nlayer,N_active
    INTEGER                                 :: k
    REAL(wp),                    INTENT(in) :: time
    REAL(wp), DIMENSION(Nlayer), INTENT(in) :: T,thick,S_bu,psi_s,psi_l,psi_g

    DO k = 1,Nlayer
       IF ( k<=N_active ) THEN
          WRITE(k+200, '(F8.4,2x,f10.3,2x,f7.5,2x,f9.5,2x,f4.2,2x,f4.2,2x,f4.2)') &
               time/(86400._wp),T(k),thick(k),S_bu(k),psi_s(k),psi_l(k),psi_g(k)
       ELSE
          WRITE(k+200,'(F6.2,2x,f7.3,2x,f5.3,2x,f8.5,2x,f4.2,2x,f4.2,2x,f4.2)')time/(86400._wp),0.0,0.0,0.0,0.0,0.0,0.0
       END IF
    END DO

  END SUBROUTINE output_raw

  !>
  !! Output for debugging purposes.
  !!
  !! Data of snow layer is written out at each time step to aid in finding errors or understanding model behavior.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-10-11)
  !!
  SUBROUTINE output_raw_snow(time,T_snow,thick_snow,S_abs_snow,m_snow,psi_s_snow,psi_l_snow,psi_g_snow)
    REAL(wp), INTENT(in) :: time
    REAL(wp), INTENT(in) :: T_snow,thick_snow,S_abs_snow,m_snow,psi_s_snow,psi_l_snow,psi_g_snow

    IF (thick_snow>0.0_wp) THEN
       WRITE(66, '(F8.4,2x,f10.3,2x,f5.3,2x,f4.1,2x,f4.2,2x,f4.2,2x,f4.2)') &
            time/(30._wp*86400._wp),T_snow,thick_snow,S_abs_snow/MAX(m_snow,0.001_wp),psi_s_snow,psi_l_snow,psi_g_snow
    ELSE 
       WRITE(66, '(F8.4,2x,f10.3,2x,f5.3,2x,f4.1,2x,f4.2,2x,f4.2,2x,f4.2)') &
            time/(30._wp*86400._wp),0.0,0.0,0.0,0.0,0.0,0.0
    END IF

  END SUBROUTINE output_raw_snow

  !>
  !! Output for debugging layer dynamics..
  !!
  !! Is used when debug_flag = 2 to track when which layer dynamics occur (see mo_layer_dynamics).
  !!
  !!

  SUBROUTINE output_raw_lay(Nlayer,N_active,H_abs,m,S_abs,thick,string)
    INTEGER,                     INTENT(in) :: Nlayer,N_active
    INTEGER                                 :: k
    REAL(wp), DIMENSION(Nlayer), INTENT(in) :: H_abs,S_abs,thick,m 
    REAL(wp)                                :: mm 
    CHARACTER*6,                 INTENT(in) :: string

    DO k = 1,Nlayer
       IF ( k<=N_active ) THEN
          IF (m(k)==0.0_wp) THEN        
             mm = 99999999._wp
          ELSE
             mm = m(k)
          END IF
          WRITE(k+200, '(A6,2x,f11.0,2x,f6.3,2x,f8.5,2x,f9.3,2x,I2)') &
               string,h_abs(k),thick(k),s_abs(k)/mm,mm/MAX(thick(k),0.0000000000000000001_wp),N_active
       ELSE 
          WRITE(k+200, '(A6,2x,f11.0,2x,f6.3,2x,f8.1,2x,f9.3,2x,I2)')string,0.0,0.0,0.0,0.0,N_active
       END IF
    END DO
  END SUBROUTINE output_raw_lay

  !>
  !! Output files are opened and format strings are created.
  !!
  !! Format strings are defined according to the number of layers used which define the output format.
  !! Files are opened.
  !!
  !! @par Revision History
  !! created by Philipp Griewank, IMPRS (2010-10-11)
  !! moved by Philipp Griewank, IMPRS (2011-03-09)
  !!
  SUBROUTINE output_begin(Nlayer,debug_flag,format_T,format_psi,format_thick,format_snow,format_T2m_top,format_perm,&
                          &format_melt)
    INTEGER,         INTENT(in)  :: Nlayer,debug_flag
    CHARACTER*12000, INTENT(out) :: format_T,format_psi,format_thick,format_snow,format_T2m_top,format_perm,&
                                    &format_melt
    INTEGER                      :: k
    CHARACTER*12000              :: output_string


    IF (debug_flag==2) THEN
       DO k=1,Nlayer
          IF (k<10) THEN
             WRITE(output_string,'(A18,I1,A4)')'./output/thermo00',k,'.txt'
             OPEN(k+200,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
          ELSE IF (k<100) THEN
             WRITE(output_string,'(A17,I2,A4)')'./output/thermo0',k,'.txt'
             OPEN(k+200,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
          ELSE
             WRITE(output_string,'(A16,I3,A4)')'./output/thermo',k,'.txt'
             OPEN(k+200,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
          END IF
       END DO
    END IF

    format_T     = '(F9.3,2x'
    format_perm  = '(ES14.7,2x'
    format_psi   = '(F9.3,2x'
    format_thick = '(F9.5,2x'

    format_snow    = '(F9.3,2x,F9.3,2x,F9.3,2x,F9.3)'
    format_T2m_top = '(F9.3,2x,F9.3,2x,F9.3)'
    format_melt    = '(ES14.7,2x,ES14.7,2x,ES14.7)'

    DO k=1,Nlayer
       format_T     = TRIM(format_T)//',f9.3,2x'
       format_psi   = TRIM(format_psi)//',f9.3,2x'
       format_thick = TRIM(format_thick)//',f9.5,2x'
       format_perm  = TRIM(format_perm)//',ES14.7,2x'
    END DO

    format_T     = TRIM(format_T)//')'
    format_psi   = TRIM(format_psi)//')'
    format_thick = TRIM(format_thick)//')'
    format_perm  = TRIM(format_perm)//')'

    
    OPEN(30,file='./output/dat_T.dat'           ,STATUS='replace',Recl=12288)
    OPEN(31,file='./output/dat_psi_s.dat'       ,STATUS='replace',Recl=12288)
    OPEN(32,file='./output/dat_thick.dat'       ,STATUS='replace',Recl=12288)
    OPEN(33,file='./output/dat_S_bu.dat'        ,STATUS='replace',Recl=12288)
    OPEN(34,file='./output/dat_ray.dat'         ,STATUS='replace',Recl=12288)
    OPEN(35,file='./output/dat_psi_l.dat'       ,STATUS='replace',Recl=12288)
    OPEN(40,file='./output/dat_freeboard.dat'   ,STATUS='replace',Recl=12288)
    OPEN(41,file='./output/dat_snow.dat'        ,STATUS='replace',Recl=12288)
    OPEN(42,file='./output/dat_vital_signs.dat' ,STATUS='replace',Recl=12288)
    OPEN(43,file='./output/dat_grav_drain.dat'  ,STATUS='replace',Recl=12288)
    OPEN(45,file='./output/dat_T2m_T_top.dat'   ,STATUS='replace',Recl=12288)
    OPEN(46,file='./output/dat_perm.dat'        ,STATUS='replace',Recl=12288)
    OPEN(47,file='./output/dat_flush_v.dat'        ,STATUS='replace',Recl=12288)
    OPEN(48,file='./output/dat_flush_h.dat'        ,STATUS='replace',Recl=12288)
    OPEN(49,file='./output/dat_psi_g.dat'        ,STATUS='replace',Recl=12288)
    OPEN(50,file='./output/dat_melt.dat'        ,STATUS='replace',Recl=12288)

    OPEN(66,file='./output/snow.txt'            ,STATUS='replace')

  END SUBROUTINE output_begin
  
  
  !>
  !! Output files for bgc are opened and format strings are created.
  !!
  !! Same thing as out_begin but for bgc
  !! Each tracer is outputted in bulk and in brine concentration in a separate file.
  !! Added ADJUSTL to the output strings because they got wierd 
  !!
  !! @par Revision History
  !! created by Dr. Philipp Griewank, MPI (2014-02-07)
  !! fix by Dr. Philipp Griewank, UniK (2018-05-18)
  SUBROUTINE output_begin_bgc(Nlayer,N_bgc,format_bgc)
    INTEGER,         INTENT(in)  :: Nlayer,N_bgc
    CHARACTER*12000, INTENT(out) :: format_bgc
    INTEGER                      :: k
    CHARACTER*12000              :: output_string

     

    format_bgc = '(F16.8,2x'

    DO k=1,Nlayer
      format_bgc = TRIM(format_bgc)//',f16.8,2x'
    END DO
    
    format_bgc = TRIM(format_bgc)//')'

    DO k=1,N_bgc
       IF (k<10) THEN
          WRITE(output_string,'(A18,I1,A7)')'./output/dat_bgc0',k,'.bu.dat'
          OPEN(2*k+400,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
          WRITE(output_string,'(A18,I1,A7)')'./output/dat_bgc0',k,'.br.dat'
          OPEN(2*k+401,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
       ELSE 
          WRITE(output_string,'(A17,I2,A7)')'./output/dat_bgc',k,'.bu.dat'
          OPEN(2*k+400,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
          WRITE(output_string,'(A17,I2,A7)')'./output/dat_bgc',k,'.br.dat'
          OPEN(2*k+401,file=TRIM(ADJUSTL(output_string)),STATUS='replace',Recl=12288)
       END IF
    END DO

  END SUBROUTINE output_begin_bgc

END MODULE mo_output
