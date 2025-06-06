
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init3
use modmain
use modgw
use modtddft
use modvars
use modmpi
implicit none
! local variables
integer ig,iw
real(8) w1,w2,t1,t2

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! G-vectors for response functions
ngrf=1
do ig=2,ngvec
  if (gc(ig) > gmaxrf) then
    ngrf=ig-1
    exit
  end if
end do
ngrf=min(ngrf,ngvc)
! write the G-vectors to file
if (mp_mpi) call writegvecrf
! frequencies for reponse functions
nwrf=1
if (allocated(wrf)) deallocate(wrf)
if (any(task == [188,320,330,331])) then
  nwrf=nwplot
  allocate(wrf(nwrf))
  w1=wplot(1)
  w2=max(wplot(2),w1)
  t1=(w2-w1)/dble(nwplot)
  do iw=1,nwplot
    t2=w1+t1*dble(iw-1)
    wrf(iw)=cmplx(t2,swidth,8)
  end do
! set the first frequency to zero for the bootstrap functional
  if ((fxctype(1) == 210).or.(fxctype(1) == 211)) then
    wrf(1)=cmplx(0.d0,swidth,8)
  end if
else if (any(task == [600,610,620,640]).or.ksgwrho) then
! GW Matsubara frequencies
  call genwgw
else
  nwrf=1
  allocate(wrf(nwrf))
  wrf(1)=cmplx(0.d0,swidth,8)
end if
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('gmaxrf',rv=gmaxrf)
  call writevars('ngrf',iv=ngrf)
  call writevars('nwrf',iv=nwrf)
  call writevars('wrf',nv=nwrf,zva=wrf)
end if

end subroutine

