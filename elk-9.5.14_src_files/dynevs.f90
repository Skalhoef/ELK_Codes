
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynevs(ev,a,w)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(in) :: ev(nbph,nbph)
complex(8), intent(inout) :: a(nbph,nbph)
real(8), intent(out) :: w(nbph)
! local variables
integer i,j,k
real(8) t1,t2
complex(8) z1
! automatic arrays
real(8) wt(nbph)
! find the eigenvalues and eigenvectors of the matrix a
call eveqnzh(nbph,nbph,a,w)
! reorder eigenvalues so that the eigenvectors maximally overlap with ev
wt(:)=w(:)
do i=1,nbph
  j=1
  t1=0.d0
  do k=1,nbph
    z1=dot_product(ev(:,i),a(:,k))
    t2=dble(z1)**2+aimag(z1)**2
    if (t2 > t1) then
      j=k
      t1=t2
    end if
  end do
  w(i)=wt(j)
end do
end subroutine

