
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: brzint
! !INTERFACE:
subroutine brzint(nsm,ngridk,nsk,ivkik,nw,wint,n,ld,e,f,g)
! !USES:
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   nsm    : level of smoothing for output function (in,integer)
!   ngridk : k-point grid size (in,integer(3))
!   nsk    : k-point subdivision grid size (in,integer(3))
!   ivkik  : map from (i1,i2,i3) to k-point index
!            (in,integer(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
!   nw     : number of energy divisions (in,integer)
!   wint   : energy interval (in,real(2))
!   n      : number of functions to integrate (in,integer)
!   ld     : leading dimension (in,integer)
!   e      : array of energies as a function of k-points (in,real(ld,*))
!   f      : array of weights as a function of k-points (in,real(ld,*))
!   g      : output function (out,real(nw))
! !DESCRIPTION:
!   Given energy and weight functions, $e$ and $f$, on the Brillouin zone and a
!   set of equidistant energies $\omega_i$, this routine computes the integrals
!   $$ g(\omega_i)=\frac{\Omega}{(2\pi)^3}\int_{\rm BZ} f({\bf k})
!    \delta(\omega_i-e({\bf k}))d{\bf k}, $$
!   where $\Omega$ is the unit cell volume. This is done by first interpolating
!   $e$ and $f$ on a finer $k$-point grid using the trilinear method. Then for
!   each $e({\bf k})$ on the finer grid the nearest $\omega_i$ is found and
!   $f({\bf k})$ is accumulated in $g(\omega_i)$. If the output function is
!   noisy then either {\tt nsk} should be increased or {\tt nw} decreased.
!   Alternatively, the output function can be artificially smoothed up to a
!   level given by {\tt nsm}. See routine {\tt fsmooth}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!   Improved efficiency, May 2007 (Sebastian Lebegue)
!   Added parallelism, March 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nsm,ngridk(3),nsk(3)
integer, intent(in) :: ivkik(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
integer, intent(in) :: nw
real(8), intent(in) :: wint(2)
integer, intent(in) :: n,ld
real(8), intent(in) :: e(ld,*),f(ld,*)
real(8), intent(out) :: g(nw)
! local variables
integer nk,i1,i2,i3,j1,j2,j3,k1,k2,k3,i,iw,nthd
integer i000,i001,i010,i011,i100,i101,i110,i111
real(8) wd,dw,dwi,w1,t1,t2
! automatic arrays
real(8) f0(n),f1(n),e0(n),e1(n)
real(8) f00(n),f01(n),f10(n),f11(n)
real(8) e00(n),e01(n),e10(n),e11(n)
if ((ngridk(1) < 1).or.(ngridk(2) < 1).or.(ngridk(3) < 1)) then
  write(*,*)
  write(*,'("Error(brzint): ngridk < 1 : ",3I8)') ngridk
  write(*,*)
  stop
end if
if ((nsk(1) < 1).or.(nsk(2) < 1).or.(nsk(3) < 1)) then
  write(*,*)
  write(*,'("Error(brzint): nsk < 1 : ",3I8)') nsk
  write(*,*)
  stop
end if
! total number of k-points
nk=ngridk(1)*ngridk(2)*ngridk(3)
! length of interval
wd=wint(2)-wint(1)
! energy step size
dw=wd/dble(nw)
dwi=1.d0/dw
w1=wint(1)*dwi
g(:)=0.d0
call holdthd(nk,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(f0,f1,e0,e1) &
!$OMP PRIVATE(f00,f01,f10,f11) &
!$OMP PRIVATE(e00,e01,e10,e11) &
!$OMP PRIVATE(k1,k2,k3,i1,i2,i3) &
!$OMP PRIVATE(i000,i001,i010,i011) &
!$OMP PRIVATE(i100,i101,i110,i111) &
!$OMP PRIVATE(t1,t2,i,iw) &
!$OMP REDUCTION(+:g) &
!$OMP NUM_THREADS(nthd)
!$OMP DO COLLAPSE(3)
do j1=0,ngridk(1)-1
  do j2=0,ngridk(2)-1
    do j3=0,ngridk(3)-1
      k1=mod(j1+1,ngridk(1))
      k2=mod(j2+1,ngridk(2))
      k3=mod(j3+1,ngridk(3))
      i000=ivkik(j1,j2,j3); i001=ivkik(j1,j2,k3)
      i010=ivkik(j1,k2,j3); i011=ivkik(j1,k2,k3)
      i100=ivkik(k1,j2,j3); i101=ivkik(k1,j2,k3)
      i110=ivkik(k1,k2,j3); i111=ivkik(k1,k2,k3)
      do i1=0,nsk(1)-1
        t2=dble(i1)/dble(nsk(1))
        t1=1.d0-t2
        f00(:)=f(:,i000)*t1+f(:,i100)*t2
        f01(:)=f(:,i001)*t1+f(:,i101)*t2
        f10(:)=f(:,i010)*t1+f(:,i110)*t2
        f11(:)=f(:,i011)*t1+f(:,i111)*t2
        t1=t1*dwi
        t2=t2*dwi
        e00(:)=e(:,i000)*t1+e(:,i100)*t2-w1
        e01(:)=e(:,i001)*t1+e(:,i101)*t2-w1
        e10(:)=e(:,i010)*t1+e(:,i110)*t2-w1
        e11(:)=e(:,i011)*t1+e(:,i111)*t2-w1
        do i2=0,nsk(2)-1
          t2=dble(i2)/dble(nsk(2))
          t1=1.d0-t2
          f0(:)=f00(:)*t1+f10(:)*t2
          f1(:)=f01(:)*t1+f11(:)*t2
          e0(:)=e00(:)*t1+e10(:)*t2
          e1(:)=e01(:)*t1+e11(:)*t2
          do i3=0,nsk(3)-1
            t2=dble(i3)/dble(nsk(3))
            t1=1.d0-t2
            do i=1,n
              iw=nint(e0(i)*t1+e1(i)*t2)+1
              if ((iw >= 1).and.(iw <= nw)) g(iw)=g(iw)+f0(i)*t1+f1(i)*t2
            end do
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! normalise function
t1=dw*dble(nk)*dble(nsk(1)*nsk(2)*nsk(3))
t1=1.d0/t1
g(:)=t1*g(:)
! smooth output function if required
if (nsm > 0) call fsmooth(nsm,nw,g)
end subroutine
!EOC

