! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TRANS_END(CDMODE)

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_PTR, C_LOC, C_NULL_PTR, C_NULL_CHAR

IMPLICIT NONE

INTERFACE
SUBROUTINE TRANS_END_CPP_BINDING(CDMODE) BIND(C, NAME="trans_end")
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
  TYPE(C_PTR), VALUE, INTENT(IN) :: CDMODE
END SUBROUTINE TRANS_END_CPP_BINDING
END INTERFACE

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDMODE

! Fixed-size, null-terminated C-string buffer for the optional mode string. A null C pointer
! indicates the argument was not present in the Fortran call.
INTEGER, PARAMETER :: JP_MAXLEN = 32
CHARACTER(KIND=C_CHAR, LEN=1), TARGET :: ZMODE(JP_MAXLEN+1)
TYPE(C_PTR) :: CDMODE_PTR
INTEGER :: ILEN, JCHAR

IF (PRESENT(CDMODE)) THEN
  ILEN = MIN(LEN(CDMODE), JP_MAXLEN)
  DO JCHAR = 1, ILEN
    ZMODE(JCHAR) = CDMODE(JCHAR:JCHAR)
  ENDDO
  ZMODE(ILEN+1) = C_NULL_CHAR
  CDMODE_PTR = C_LOC(ZMODE)
ELSE
  CDMODE_PTR = C_NULL_PTR
ENDIF

CALL TRANS_END_CPP_BINDING(CDMODE_PTR)

END SUBROUTINE TRANS_END
