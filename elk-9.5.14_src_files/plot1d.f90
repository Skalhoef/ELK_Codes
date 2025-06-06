
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot1d
! !INTERFACE:
subroutine plot1d(fnum1,fnum2,nf,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum1 : plot file number (in,integer)
!   fnum2 : vertex location file number (in,integer)
!   nf    : number of functions (in,integer)
!   rfmt  : real muffin-tin function (in,real(npmtmax,natmtot,nf))
!   rfir  : real intersitial function (in,real(ngtot,nf))
! !DESCRIPTION:
!   Produces a 1D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} along the lines connecting the vertices in the global array
!   {\tt vvlp1d}. See routine {\tt rfplot}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum1,fnum2,nf
real(8), intent(in) :: rfmt(npmtmax,natmtot,nf),rfir(ngtot,nf)
! local variables
integer jf,ip,iv
real(8) fmin,fmax,t1
! allocatable arrays
real(8), allocatable :: fp(:,:)
if ((nf < 1).or.(nf > 4)) then
  write(*,*)
  write(*,'("Error(plot1d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
allocate(fp(npp1d,nf))
! connect the 1D plotting vertices
call plotpt1d(avec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
do jf=1,nf
! evaluate function at each point
  call rfplot(npp1d,vplp1d,rfmt(:,:,jf),rfir(:,jf),fp(:,jf))
end do
do ip=ip01d,npp1d
! write the point distances and function to file
  write(fnum1,'(5G18.10)') dpp1d(ip),(fp(ip,jf),jf=1,nf)
end do
! write the vertex location lines
fmin=minval(fp(:,:))
fmax=maxval(fp(:,:))
t1=0.5d0*(fmax-fmin)
fmin=fmin-t1
fmax=fmax+t1
do iv=1,nvp1d
  write(fnum2,'(2G18.10)') dvp1d(iv),fmin
  write(fnum2,'(2G18.10)') dvp1d(iv),fmax
  write(fnum2,*)
end do
deallocate(fp)
end subroutine
!EOC

