
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalsv(fext,ikp,vpl,evalsvp)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: evalsvp(nstsv)
! local variables
integer isym,ik,ist
integer recl,nstsv_
real(8) vkl_(3),t1
if (ikp > 0) then
  ik=ikp
else
! find the k-point number
  call findkpt(vpl,isym,ik)
end if
! find the record length
inquire(iolength=recl) vkl_,nstsv_,evalsvp
!$OMP CRITICAL(u204)
open(204,file='EVALSV'//trim(fext),form='UNFORMATTED',access='DIRECT',recl=recl)
read(204,rec=ik) vkl_,nstsv_,evalsvp
close(204)
!$OMP END CRITICAL(u204)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1 > epslat) then
  write(*,*)
  write(*,'("Error(getevalsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALSV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv /= nstsv_) then
  write(*,*)
  write(*,'("Error(getevalsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVALSV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! apply scissor correction if required
if (tscissor) then
  do ist=1,nstsv
    if (evalsvp(ist).gt.efermi) evalsvp(ist)=evalsvp(ist)+scissor
  end do
end if
end subroutine

