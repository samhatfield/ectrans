####################################################################
# COMPILER
####################################################################

set( CMAKE_C_COMPILER nvc )
set( CMAKE_CXX_COMPILER nvc++ )
set( CMAKE_Fortran_COMPILER nvfortran )

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-mp -mp=bind,allcores,numa" )
set( OpenMP_CXX_FLAGS           "-mp -mp=bind,allcores,numa" )
set( OpenMP_Fortran_FLAGS       "-mp -mp=bind,allcores,numa" )

####################################################################
# OpenAcc FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc=gpu -gpu=cc90,lineinfo,fastmath,rdc" )

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 90)
endif()

####################################################################
# COMMON FLAGS
####################################################################

set( ECBUILD_Fortran_FLAGS "-fpic" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mframe" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee" )
set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz" )

set( ECBUILD_Fortran_FLAGS_BIT "-O2 -gopt" )

set( ECBUILD_C_FLAGS "-O2 -gopt -traceback" )

set( ECBUILD_CXX_FLAGS "-O2 -gopt -DEIGEN_DONT_VECTORIZE" )

# We use the GCC perl installation as something is broken in the NVHPC perl installation
# PATH=/usr/local/apps/eb/env/release/2025/software/Perl/5.40.2-GCCcore-15.1.0/bin:$PATH
# PERL5LIB=/usr/local/apps/eb/env/release/2025/software/Perl/5.40.2-GCCcore-15.1.0/lib/perl5/5.40.2

# Necessary when CUDA math libs are installed in a different location to cudart
set( CMAKE_EXE_LINKER_FLAGS "-L/$ENV{NVHPC_ROOT}/math_libs/13.0/lib64" )

# A workaround for a quirky interaction between FIAT and ecWAM
set( HAVE_MPI_F08 OFF )
