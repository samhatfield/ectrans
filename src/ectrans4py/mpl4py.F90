SUBROUTINE MPL_INIT4PY(KRANK, KSIZE)
! ** PURPOSE
!    Initialise FIAT's MPL (Message Passing Library) for distributed ecTrans.
!    MPL_INIT detects an already-initialised MPI (e.g. from mpi4py) and attaches
!    to MPI_COMM_WORLD, so this must be called *after* MPI is up.
!
! ** RETURNS
!    KRANK : this task's MPL rank (1-based)
!    KSIZE : total number of MPI tasks
USE ISO_FORTRAN_ENV, ONLY: INT64
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRANK, KSIZE
CALL MPL_INIT(LDINFO=.FALSE.)
KRANK = MPL_MYRANK()
KSIZE = MPL_NPROC()
END SUBROUTINE MPL_INIT4PY


SUBROUTINE MPL_END4PY()
! ** PURPOSE
!    Finalise FIAT's MPL. (Does not finalise MPI itself.)
USE MPL_MODULE, ONLY: MPL_END
IMPLICIT NONE
CALL MPL_END(LDMEMINFO=.FALSE.)
END SUBROUTINE MPL_END4PY


SUBROUTINE MPL_INIT_COMM4PY(KCOMM, KRANK, KSIZE)
! ** PURPOSE
!    Initialise FIAT's MPL bound to a *user* communicator instead of
!    MPI_COMM_WORLD, via FIAT's LMPLUSERCOMM seam (the same mechanism the
!    OOPS-coupled IFS uses). This makes every MPL collective -- and therefore all
!    of ecTrans -- local to that communicator: when the world is split into
!    compute + IO-server groups, the transforms live on the compute group and a
!    server rank never has to join them.
!
! ** ARGUMENTS
!    KCOMM : Fortran handle of the communicator (mpi4py: comm.py2f())
! ** RETURNS
!    KRANK : this task's MPL rank within KCOMM (1-based)
!    KSIZE : number of tasks in KCOMM
USE ISO_FORTRAN_ENV, ONLY: INT64
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC
USE MPL_DATA_MODULE, ONLY: LMPLUSERCOMM, MPLUSERCOMM
USE PARKIND1, ONLY: JPIM
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(IN) :: KCOMM
INTEGER(KIND=INT64), INTENT(OUT) :: KRANK, KSIZE
LMPLUSERCOMM = .TRUE.
MPLUSERCOMM = INT(KCOMM, KIND=JPIM)
CALL MPL_INIT(LDINFO=.FALSE.)
KRANK = MPL_MYRANK()
KSIZE = MPL_NPROC()
END SUBROUTINE MPL_INIT_COMM4PY
