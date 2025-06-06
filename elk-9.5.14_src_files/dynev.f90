
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynev(dynp,w,ev)
use modmain
use modphonon
implicit none
! arguments
complex(8), intent(in) :: dynp(nbph,nbph)
real(8), intent(out) :: w(nbph)
complex(8), intent(out) :: ev(nbph,nbph)
! local variables
integer is,ia,js,ja
integer ip,jp,i,j
real(8) t1
ev(:,:)=0.d0
i=0
do is=1,nspecies
  do ia=1,natoms(is)
    do ip=1,3
      i=i+1
      j=0
      do js=1,nspecies
! mass factor
        if ((spmass(is) <= 0.d0).or.(spmass(js) <= 0.d0)) then
! infinite mass
          t1=0.d0
        else
          t1=1.d0/sqrt(spmass(is)*spmass(js))
        end if
        do ja=1,natoms(js)
          do jp=1,3
            j=j+1
            if (i <= j) then
! use Hermitian average of dynamical matrix
              ev(i,j)=0.5d0*t1*(dynp(i,j)+conjg(dynp(j,i)))
            end if
          end do
        end do
      end do
    end do
  end do
end do
! find the eigenvalues and eigenvectors of the dynamical matrix
call eveqnzh(nbph,nbph,ev,w)
do i=1,nbph
  if (w(i) >= 0.d0) then
    w(i)=sqrt(w(i))
  else
    w(i)=-sqrt(abs(w(i)))
  end if
end do
end subroutine

