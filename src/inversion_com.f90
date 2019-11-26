! Routines required for reading undersampled solutions 
! from dynamic inversions using fd3d_pt.
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
      MODULE inversion_com
      INTEGER:: RUNI,NLI,NWI
      REAL,ALLOCATABLE,DIMENSION(:,:):: DcI,TsI,T0I     !Test variables (for which misfit is calculated)
      real,allocatable,dimension(:,:,:):: DcA,TsA,T0A   !Array of variables in MC chains:
      REAL,ALLOCATABLE,DIMENSION(:,:,:):: ruptimeA,riseA,slipA,schangeA
      REAL,ALLOCATABLE :: pgaA(:,:,:),mwA(:),ruptdistA(:,:),MomentRateA(:,:)
      real,allocatable,dimension(:):: VRA
      INTEGER randseed,StepType
      REAL StepSizeT0,StepSizeTs,StepSizeD
      END MODULE
