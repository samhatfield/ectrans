! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GET_CURRENT(KRESOL, LDLAM)

USE PARKIND1,                    ONLY: JPIM
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_BOOL, C_PTR, C_LOC, C_NULL_PTR

IMPLICIT NONE

INTERFACE
SUBROUTINE GET_CURRENT_CPP_BINDING(KRESOL, KLAM) BIND(C, NAME="get_current")
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
  TYPE(C_PTR), VALUE, INTENT(IN) :: KRESOL
  TYPE(C_PTR), VALUE, INTENT(IN) :: KLAM
END SUBROUTINE GET_CURRENT_CPP_BINDING
END INTERFACE

INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KRESOL
LOGICAL,            OPTIONAL, INTENT(OUT) :: LDLAM

! C-interoperable handles for the optional output arguments. A null C pointer indicates the argument
! was not present in the Fortran call; otherwise the pointer targets a local temporary which the
! C++ layer writes into and which is copied back after the call.
TYPE(C_PTR) :: KRESOL_PTR, KLAM_PTR
INTEGER(KIND=C_INT),  TARGET :: IRESOL
LOGICAL(KIND=C_BOOL), TARGET :: LLAM

IF (PRESENT(KRESOL)) THEN
  KRESOL_PTR = C_LOC(IRESOL)
ELSE
  KRESOL_PTR = C_NULL_PTR
ENDIF
IF (PRESENT(LDLAM)) THEN
  KLAM_PTR = C_LOC(LLAM)
ELSE
  KLAM_PTR = C_NULL_PTR
ENDIF

CALL GET_CURRENT_CPP_BINDING(KRESOL_PTR, KLAM_PTR)

IF (PRESENT(KRESOL)) KRESOL = IRESOL
IF (PRESENT(LDLAM))  LDLAM = LLAM

END SUBROUTINE GET_CURRENT
