
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine nuclei
use modmain
implicit none
! local variables
integer is,ir
! external functions
real(8), external :: radnucl
do is=1,nspecies
! approximate nuclear radius
  rnucl(is)=radnucl(spzn(is))
! nuclear volume
  volnucl(is)=(4.d0/3.d0)*pi*rnucl(is)**3
! number of radial mesh points to nuclear radius
  nrnucl(is)=1
  do ir=1,nrmt(is)
    if (rsp(ir,is) > rnucl(is)) then
      nrnucl(is)=ir
      exit
    end if
  end do
end do
end subroutine

