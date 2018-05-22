!>
!! Module contains changes specific testcases require during the main timeloop. 
!!
!! Most settings related to the testcases are defined in mo_init, but if changes to the code need to applied after the timestepping has begun they are located here. 
!! Changes were initially simply implemented in the main timeloop, but things got confusing. 
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
!!
!!
!! @par Revision History
!! Removed from mo_grotz by Philipp Griewank, IMPRS (2014-04-16)
MODULE mo_testcase_specifics

  USE mo_parameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC::sub_test1,sub_test2,sub_test3,sub_test4,sub_test6,sub_test9,sub_test34


CONTAINS

  !>
  !! Subroutine for changing T_top for testcase 1.
  !!
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, IMPRS (2014-04-16)
  !!
  SUBROUTINE sub_test1 (time,T_top)
    REAL(wp),                    INTENT(in) :: time
    REAL(wp),                   INTENT(inout) :: T_top

          IF (ABS(time-12.*3600.)<0.01) THEN
             PRINT*,'T_top switching commences'
             T_top = -10._wp
          ELSE IF (ABS(time-24.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-36.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-48.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-60.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-72.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-84.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-96.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-108.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-120.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-132.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-144.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-156.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-168.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-180.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-192.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-204.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-216.*3600.)<0.01) THEN
             T_top = -5._wp
          ELSE IF (ABS(time-228.*3600.)<0.01) THEN
             T_top = -10._wp
          ELSE IF (ABS(time-240.*3600.)<0.01) THEN
             T_top = -5._wp
          END IF

 END SUBROUTINE sub_test1

  !>
  !! Subroutine for changing T_top for testcase 2.
  !!
  !! T2m is adjusted over time.
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, IMPRS (2014-04-17)
  !!
  SUBROUTINE sub_test2 (time,T2m)
    REAL(wp),                    INTENT(in) :: time
    REAL(wp),                   INTENT(inout) :: T2m

   
 
   IF      (time>86400._wp*25._wp) THEN
      T2m = 15._wp
   Else IF      (time>86400._wp*15._wp) THEN
      T2m = 1._wp
   END IF
   
   END SUBROUTINE sub_test2
   
  !>
  !! Subroutine for changing T2m for testcase 9.
  !!
  !! T2m is adjusted over time.
  !!
  !! @par Revision History
  !! Formed by Niels Fuchs, MPI (2016-01-18)
  !!
  SUBROUTINE sub_test9 (time,T2m)
    REAL(wp),                    INTENT(in) :: time
    REAL(wp),                   INTENT(inout) :: T2m

   
 
   IF (time<(19.75_wp*3600._wp)) THEN
      T2m = 0._wp
   ELSE IF (time<(86400._wp*3._wp+2.25_wp*3600._wp)) THEN
      T2m = -15._wp
   ELSE
      T2m = 1._wp
   END IF
   
   
 END SUBROUTINE sub_test9
 
 !>
  !! Subroutine for changing T2m for testcase 34.
  !!
  !! T2m is adjusted over time.
  !!
  !! @par Revision History
  !! adjusted by Niels Fuchs, MPI (2016-01-18)
  !!
  SUBROUTINE sub_test34 (time,T2m)
    REAL(wp),                    INTENT(in) :: time
    REAL(wp),                   INTENT(inout) :: T2m

   
 
   IF (time<2._wp*3600._wp) THEN
      T2m = 0._wp
   ELSE IF (time<(86400._wp*5._wp)) THEN
      T2m = -15._wp
   ELSE IF (time<(86400._wp*7._wp)) THEN
      T2m = -5._wp
   ELSE
      T2m = 1._wp
   END IF
   
 END SUBROUTINE sub_test34

  !>
  !! Subroutine for setting snow for testcase 3.
  !!
  !! Precipitation rates are set
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, (2014-04-18)
  !!
  SUBROUTINE sub_test3 (time,liquid_precip,solid_precip)
    REAL(wp),  INTENT(in)    :: time
    REAL(wp)                 :: day
    REAL(wp),  INTENT(inout) :: liquid_precip, solid_precip 

     day    = time/86400._wp
    DO WHILE (day>360)
       day = day-360
    END DO
   
 liquid_precip  = 0._wp     
 solid_precip   = 0.15_wp/86400._wp/356._wp

   
   
 END SUBROUTINE sub_test3

  !>
  !! Subroutine for setting snow for testcase 4.
  !!
  !! 
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, (2014-04-18)
  !!
  SUBROUTINE sub_test4 (time,fl_q_bottom)
    REAL(wp),  INTENT(in)    :: time
    REAL(wp),  INTENT(inout) :: fl_q_bottom 
    fl_q_bottom = -7._wp*SIN(time*(2._wp*pi)/(86400._wp*365._wp))+7._wp 

  END SUBROUTINE sub_test4

  !>
  !! Subroutine for changing T_top for testcase 6 which seeks to reproduce lab measurements of Roni Glud.
  !!
  !!
  !! @par Revision History
  !! Formed by Philipp Griewank, IMPRS (2014-04-38)
  !!
  SUBROUTINE sub_test6 (time,T2m)
    REAL(wp),              INTENT(in)    :: time
    REAL(wp),              INTENT(inout) :: T2m
          IF      (time>1714._wp*60._wp) THEN
             T2m = -19._wp
          ELSE IF (time>1676._wp*60._wp) THEN
             T2m = -5._wp
          ELSE IF (time>1525._wp*60._wp) THEN
             T2m = -18._wp
          ELSE IF (time>1483._wp*60._wp) THEN
             T2m = -5._wp
          ELSE IF (time>1385._wp*60._wp) THEN
             T2m = -18._wp
          ELSE IF (time>1349._wp*60._wp) THEN
             T2m = -5._wp
          ELSE IF (time>1160._wp*60._wp) THEN
             T2m = -18._wp
          ELSE IF (time>1100._wp*60._wp) THEN
             T2m = -5._wp
          END IF

 END SUBROUTINE sub_test6

END MODULE mo_testcase_specifics

