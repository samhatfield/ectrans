! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TRANS_RELEASE(KRESOL)

USE PARKIND1,                    ONLY: JPIM
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT

IMPLICIT NONE

INTERFACE
SUBROUTINE TRANS_RELEASE_CPP_BINDING(KRESOL) BIND(C, NAME="trans_release")
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
  INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: KRESOL
END SUBROUTINE TRANS_RELEASE_CPP_BINDING
END INTERFACE

INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL

CALL TRANS_RELEASE_CPP_BINDING(KRESOL)

END SUBROUTINE TRANS_RELEASE
