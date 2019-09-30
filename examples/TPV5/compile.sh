##pgfortran -DSCEC -DTPV5 -Mcuda=cc30,cc35,cc50,cc60 -ta=tesla:cc70 -fast -Mpreprocess -acc -Mbackslash -ofd3d_pt_GPU_TSN dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 inversion.f90
ifort -DSCEC -DTPV5 -fast -fpp -o fd3d_pt_CPU_TSN dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 inversion.f90

