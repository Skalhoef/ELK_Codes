
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gradrfmt
! !INTERFACE:
subroutine gradrfmt(nr,nri,ri,wcr,rfmt,ld,grfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on inner part of muffin-tin (in,integer)
!   ri    : 1/r on the radial mesh (in,real(nr))
!   wcr   : weights for spline coefficients on radial mesh (in,real(12,nr))
!   rfmt  : real muffin-tin function (in,real(*))
!   ld    : leading dimension (in,integer)
!   grfmt : gradient of rfmt (out,real(ld,3))
! !DESCRIPTION:
!   Calculates the gradient of a real muffin-tin function. In other words, given
!   the real spherical harmonic expansion coefficients $f_{lm}(r)$ of a function
!   $f({\bf r})$, the routine returns ${\bf F}_{lm}$ where
!   $$ \sum_{lm}{\bf F}_{lm}(r)R_{lm}(\hat{\bf r})=\nabla f({\bf r}), $$
!   and $R_{lm}$ is a real spherical harmonic function. This is done by first
!   converting the function to a complex spherical harmonic expansion and then
!   using the routine {\tt gradzfmt}. See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created August 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: ri(nr),wcr(12,nr),rfmt(*)
integer, intent(in) :: ld
real(8), intent(out) :: grfmt(ld,3)
! local variables
integer i
! automatic arrays
complex(8) zfmt(ld),gzfmt(ld,3)
! convert real to complex spherical harmonic expansion
call rtozfmt(nr,nri,rfmt,zfmt)
! compute the gradient
call gradzfmt(nr,nri,ri,wcr,zfmt,ld,gzfmt)
! convert complex to real spherical harmonic expansion
do i=1,3
  call ztorfmt(nr,nri,gzfmt(:,i),grfmt(:,i))
! improve stability by smoothing the gradient
  call rfmtsm(msmgmt,nr,nri,grfmt(:,i))
end do
end subroutine
!EOC

