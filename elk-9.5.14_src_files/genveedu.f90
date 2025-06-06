
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genveedu
! !INTERFACE:
subroutine genveedu(idu,u,j,vee)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   idu : DFT+U entry (in,integer)
!   u   : parameter U (out,real)
!   j   : parameter J (out,real)
!   vee : Coulomb matrix elements (out,real(-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,
!          -lmaxdm:lmaxdm,-lmaxdm:lmaxdm))
! !DESCRIPTION:
!   Calculates the Coulomb matrix elements used in DFT+U calculations. See
!   {\it Phys. Rev. B} {\bf 52}, 5467 (1995).
!
! !REVISION HISTORY:
!   Created November 2007 (FC,JKD,FB,LN)
!   Modified July 2009 (FC)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: idu
real(8), intent(out) :: u,j
real(8), intent(out) :: vee(-lmaxdm:lmaxdm,-lmaxdm:lmaxdm,-lmaxdm:lmaxdm, &
 -lmaxdm:lmaxdm)
! local variables
integer l,m1,m2,m3,m4,k,q
real(8) sm1,sm2,t1
! automatic arrays
real(8) :: f(0:2*lmaxdm)
! external functions
real(8), external :: gaunt
l=isldu(2,idu)
! calculate Slater integrals
call genfdu(idu,u,j,f)
do m1=-l,l
  do m2=-l,l
    do m3=-l,l
      do m4=-l,l
        sm1=0.d0
        do k=0,2*l,2
          sm2=0.d0
          do q=-k,k
            t1=gaunt(l,k,l,m1,q,m2)*gaunt(l,k,l,m3,-q,m4)
            if (mod(q,2) == 0) then
              sm2=sm2+t1
            else
              sm2=sm2-t1
            end if
          end do
          sm1=sm1+f(k)*sm2/dble(2*k+1)
        end do
        vee(m1,m3,m2,m4)=fourpi*sm1
      end do
    end do
  end do
end do
end subroutine
!EOC

