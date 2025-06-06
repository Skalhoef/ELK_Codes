
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
real(8) function wigner3j(j1,j2,j3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $3j$-symbol. There are many equivalent formulae for
!   the $3j$-symbols, the following provides high accuracy for $j\le 50$
!   \begin{align*}
!    &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}= \\
!    &(-1)^{j1+j2+m3}\sqrt{\frac{(j_1+m_1)!\,(j_2+m_2)!\,(j_3+m_3)!\,
!    (j_3-m_3)!\,(j_1-m_1)!\,(j_2-m_2)!}{(j_2-j_1+j_3)!\,(j_1-j_2+j_3)!\,
!    (j_1+j_2-j_3)!\,(1+j_1+j_2+j_3)!}}\,\sum_k(-1)^k \\
!    &\frac{(j_2-j_1+j_3)!\,(j_1-j_2+j_3)!\,(j_1+j_2-j_3)!}{(j_3-j_1-m_2+k)!\,
!    (j_3-j_2+m_1+k)!\,(j_1+j_2-j_3-k)!\,k!\,(j_1-m_1-k)!\,(j_2+m_2-k)!},
!   \end{align*}
!   where the sum is over all integers $k$ for which the factorials in the
!   summand are non-negative.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1,j2,j3
integer, intent(in) :: m1,m2,m3
! local variables
integer k,k1,k2,l1,l2,l3,n1,n2
real(8) sgn,sm,t1
! external functions
real(8), external :: factn,factr
! check input variables
if ((j1 < 0).or.(j2 < 0).or.(j3 < 0).or.(abs(m1) > j1).or.(abs(m2) > j2) &
 .or.(abs(m3) > j3)) then
  write(*,*)
  write(*,'("Error(wigner3j): invalid arguments :")')
  write(*,'("j1 = ",I8," j2 = ",I8," j3 = ",I8)') j1,j2,j3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((j1 == 0).and.(j2 == 0).and.(j3 == 0)) then
  wigner3j=1.d0
  return
end if
if ((j1 > 50).or.(j2 > 50).or.(j3 > 50)) then
  write(*,*)
  write(*,'("Error(wigner3j): angular momenta out of range : ",3I8)') j1,j2,j3
  write(*,*)
  stop
end if
l1=j2-j1+j3
l2=j1-j2+j3
l3=j1+j2-j3
if ((m1+m2+m3 /= 0).or.(l1 < 0).or.(l2 < 0).or.(l3 < 0)) then
  wigner3j=0.d0
  return
end if
n1=j1-m1
n2=j2+m2
k1=max(0,n1-l2,n2-l1)
k2=min(l3,n1,n2)
if (mod(k1-j1+j2+m3,2) /= 0) then
  sgn=-1.d0
else
  sgn=1.d0
end if
sm=0.d0
do k=k1,k2
  t1=sgn*factr(l1,l1-n2+k)*factr(l2,l2-n1+k)*factr(l3,l3-k)
  sm=sm+t1/(factn(k)*factn(n1-k)*factn(n2-k))
  sgn=-sgn
end do
t1=factr(j1+m1,l1)*factr(j2+m2,l2)*factr(j3+m3,l3)
t1=t1*factr(j3-m3,1+j1+j2+j3)*factn(j1-m1)*factn(j2-m2)
wigner3j=sm*sqrt(t1)
end function
!EOC

