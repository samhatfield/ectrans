program cpp_test_program

use parkind1, only: jpim

implicit none

#include "setup_trans0.h"
#include "setup_trans.h"
#include "dir_trans.h"

integer(kind=jpim), parameter :: nflevg = 137
integer(kind=jpim), parameter :: nsmax = 79
integer(kind=jpim), parameter :: ndgl = 2 * (nsmax + 1)

integer(kind=jpim) :: ivsetuv(nflevg)
integer(kind=jpim) :: iresol
integer(kind=jpim) :: iloen(ndgl)
integer(kind=jpim) :: i

ivsetuv = 1

call setup_trans0(kout=6, kprintlev=1)
call setup_trans(ksmax=nsmax, kdgl=ndgl, kresol=iresol)

write(6,*) "iresol = ", iresol

do i = 1, ndgl / 2
    iloen(i) = 20 + 4 * (i - 1)
    iloen(ndgl - i + 1) = iloen(i)
enddo

call setup_trans(ksmax=nsmax, kdgl=ndgl, kresol=iresol, kloen=iloen)

write(6,*) "iresol = ", iresol

call dir_trans(kvsetuv=ivsetuv)

end program cpp_test_program