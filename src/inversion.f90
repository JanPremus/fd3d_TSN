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
    
    SUBROUTINE inversion_init()
    USE inversion_com
    USE fd3dparam_com
    USE source_com
    USE SlipRates_com
    IMPLICIT NONE
    
    open(10,FILE='inputinv.dat')
    
    read(10,*)RUNI
    read(10,*)NLI,NWI
	dtseis=1
    allocate(DcI(NLI,NWI),TsI(NLI,NWI),T0I(NLI,NWI))
	allocate(MSRX(NL*NW*ntfd),MSRZ(NL*NW*ntfd),MomentRate(ntfd))
    END

    SUBROUTINE inversion_modeltofd3d() ! inversion from model controll points to fd3d grid
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    USE pml_com
    IMPLICIT NONE
    real,dimension(:),allocatable:: xintrpl,yintrpl,xnew,ynew
    integer i,k,ii,kk
    real DW,DL,dum,xs,zs,t,u

! Bilinear interpolation
    DL=dh*(nxt-2*nabc)/real(NLI-1)
    DW=dh*(nzt-nfs-nabc)/real(NWI-1)
    do k=nabc+1,nzt-nfs
      ZS=dh*(k-1-nabc)+dh/2.
      kk=int(ZS/DW)+1
      u=(ZS-DW*(kk-1))/DW
      do i=nabc+1,nxt-nabc
        XS=dh*(i-1-nabc)+dh/2.
        ii=int(XS/DL)+1
        t=(XS-DL*(ii-1))/DL
#if defined DIPSLIP
        striniZ(i,k)=(1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)
        striniX(i,k)=0.
#else
        striniX(i,k)=(1.-t)*(1.-u)*T0I(ii,kk)+t*(1.-u)*T0I(ii+1,kk)+t*u*T0I(ii+1,kk+1)+(1.-t)*u*T0I(ii,kk+1)
        striniZ(i,k)=0.
#endif
        peak_xz(i,k)=((1.-t)*(1.-u)*TsI(ii,kk)+t*(1.-u)*TsI(ii+1,kk)+t*u*TsI(ii+1,kk+1)+(1.-t)*u*TsI(ii,kk+1))*normstress(k)
        Dc(i,k)     =(1.-t)*(1.-u)*DcI(ii,kk)+t*(1.-u)*DcI(ii+1,kk)+t*u*DcI(ii+1,kk+1)+(1.-t)*u*DcI(ii,kk+1)
      enddo
    enddo
    dyn_xz=0.
    coh=0.e6

    END
    
    SUBROUTINE readinversionresult()
    USE inversion_com
    USE fd3dparam_com
    USE friction_com
    IMPLICIT NONE
    real x,z,rr,DL,DW
    integer i,j
    real dum
    
    open(244,FILE='forwardmodel.dat')
    read(244,*)dum,dum,T0I(:,:),TsI(:,:),DcI(:,:)
    close(244)
    
    CALL inversion_modeltofd3d()
    
    END



    

