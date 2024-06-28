! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_MOD
CONTAINS
SUBROUTINE FSC(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE
INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PSCALAR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PEWDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PUVDERS(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(PSCALAR(IR,JF)**2 + PSCALAR(II,JF)**2)
        ZPHASE = ATAN2(PSCALAR(II,JF),PSCALAR(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        
        PSCALAR(IR,JF) =  ZAMP*COS(ZPHASE)
        PSCALAR(II,JF) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS(IR,JF)**2 + PNSDERS(II,JF)**2)
          ZPHASE = ATAN2(PNSDERS(II,JF),PNSDERS(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
          PNSDERS(IR,JF) =  ZAMP*COS(ZPHASE)
          PNSDERS(II,JF) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(PUV(IR,JF)**2 + PUV(II,JF)**2)
        ZPHASE = ATAN2(PUV(II,JF),PUV(IR,JF)) + REAL(JM,JPRB)*ZSHIFT
        PUV(IR,JF) =  ZAMP*COS(ZPHASE)
        PUV(II,JF) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,2*KF_UV
      PUV(JLON,JF) = PUV(JLON,JF)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,KF_SCALARS
      PNSDERS(JLON,JF) = PNSDERS(JLON,JF)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,2*KF_UV
      PUVDERS(IR,JF) = -PUV(II,JF)*ZMUL
      PUVDERS(II,JF) =  PUV(IR,JF)*ZMUL
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,KF_SCALARS
      PEWDERS(IR,JF) = -PSCALAR(II,JF)*ZMUL
      PEWDERS(II,JF) =  PSCALAR(IR,JF)*ZMUL
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
