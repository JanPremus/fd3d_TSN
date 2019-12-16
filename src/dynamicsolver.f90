! Main body of the program - , initialization, 
! reading input files and running FD subroutine. 
!-------------------------------------------------------
! Authors: Frantisek Gallovic (8/2019)
! Charles University in Prague, Faculty of Mathematics and Physics

! This code is published under the GNU General Public License. To any
! licensee is given permission to modify the work, as well as to copy
! and redistribute the work or any derivative version. Still we would
! like to kindly ask you to acknowledge the authors and don't remove
! their names from the code. This code is distributed in the hope
! that it will be useful, but WITHOUT ANY WARRANTY.
! ------------------------------------------------------   
    
    PROGRAM dynamicsolver
    USE source_com
    USE friction_com
    USE fd3dparam_com
    USE inversion_com
    USE pml_com
    USE SlipRates_com
#if defined GPUMPI
    USE openacc
#endif
    IMPLICIT NONE
    integer i,k,it,itot,ncpu
    real dum
    real,allocatable:: msfts(:,:)
    logical err
    
#if defined GPUMPI
    i=acc_get_num_devices(acc_device_nvidia)
#endif

    call fd3d_init()  !Reads FD parameters
    call inversion_init()  !Reads GFs, observed waveforms
    
    write(*,*)'Running forward modeling:'

#if defined SCEC
#if defined TPV5
    CALL forwardspecialTPV5()
#endif
#if defined TPV8
    CALL forwardspecialTPV8()
#endif
#if defined TPV9
    CALL forwardspecialTPV9()
#endif
#if defined TPV104
    CALL forwardspecialTPV104()
#endif
#if defined TPV103
	print *, 'Running TPV103 benchmark' 
    CALL forwardspecialTPV103()
#endif
#else
    print *,'reading forwardmodel.dat'
    CALL readinversionresult()
#endif

	ioutput=1
	
    write(*,*)'Running dynamic rupture simulation...'

    call fd3d()
    print *,'-----------------------------------------------------'
    print *,'Average speed:       ',output_param(1)
    print *,'Seismic moment:      ',output_param(2)
    print *,'Surface of rupture:  ',output_param(3)
    print *,'Average stress drop: ',output_param(4)
 !  print *,'Energy release rate: ',output_param(5)
 !  print *,'Available energy:    ',output_param(6)
 !  print *,'kappa:               ',output_param(6)/output_param(5)
    print *,'-----------------------------------------------------'

    END PROGRAM
