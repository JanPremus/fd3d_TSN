# FD3D_TSN compile script
# tested compilers: pgfortran 19.10, gfortran 9.2.1, ifort 19.0.5
# tested environment: Intel CPU & Nvidia GPU compute capability 3+

#!/bin/bash

[ ! -d result ] && mkdir result

D_MACRO=' -DDIPSLIP -DSCEC -DTPV9 '
SRC_FILES='fd3d_init.f90 fd3d_deriv.f90 fd3d_theo.f90  dynamicsolver.f90'

PFC=pgfortran
PFC_OPTS='-fast -Mpreprocess -acc -ta=tesla:ccall -ofd3d_pgi_TSN'
# PGI target-accelerator/processor flags: -acc -ta=tesla:ccall -ta=multicore -ta=host -tp=haswell -tp=skylake -mavx512f
# PGI informational flags: -Minfo=accel -Mcuda=ptxinfo
# PGI debugging flags: 
# PGI environment variables: PGI_ACC_TIME, PGI_ACC_NOTIFY, ACC_NUM_CORES, ACC_BIND

GFC=gfortran-8
GFC_OPTS='-Ofast -march=native -cpp -o fd3d_gnu_TSN'
# gfortran accelerator/processor flags: -foffload=nvptx-none -fopenacc-dim=:4:32 -march=native
# gfortran informational flags: -fopt-info-optimized-omp
# gfortran debugging flags: -g -fbacktrace
# gfortran-9 environment variables: GOMP_DEBUG
# gfortran-9 on Ubuntu: apt install gfortran-9 gcc-9-offload-nvptx

IFC=ifort
IFC_OPTS='-Ofast -xHost -fpp -o fd3d_if_TSN'
# ifort processor flags: -xHost -xSKYLAKE-AVX512 -mtune=skylake-avx512

echo $PFC...
# -- pgfortran OpenACC (GPU)
  $PFC $D_MACRO $PFC_OPTS -acc $SRC_FILES
# $PFC $D_MACRO $PFC_OPTS -acc -ta=tesla:ccall $SRC_FILES
# -- pgfortran serial (CPU)
# $PFC $D_MACRO $PFC_OPTS $SRC_FILES
# -- pgfortran OpenMP (CPU)
# $PFC $D_MACRO $PFC_OPTS -ta=host -mp $SRC_FILES
# -- pgfortran OpenACC (CPU)
# $PFC $D_MACRO $PFC_OPTS -ta=multicore $SRC_FILES
# $PFC $D_MACRO $PFC_OPTS -ta=multicore -tp=skylake -mavx512f $SRC_FILES

echo $GFC...
# -- gfortran OpenACC (GPU) 
  $GFC $D_MACRO $GFC_OPTS -fopenacc $SRC_FILES
# -- gfortran serial (CPU)
# $GFC $D_MACRO $GFC_OPTS $SRC_FILES
# -- gfortran OpenMP (CPU)
# $GFC $D_MACRO $GFC_OPTS -fopenmp $SRC_FILES

echo $IFC...
# -- ifort serial (CPU)
  $IFC $D_MACRO $IFC_OPTS $SRC_FILES
# -- ifort OpenMP (CPU)
# $IFC $D_MACRO $IFC_OPTS -qopenmp $SRC_FILES
