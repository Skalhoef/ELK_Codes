
! Copyright (C) 2002-2012 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gengvec
! !INTERFACE:
subroutine gengvec
! !USES:
use modmain
! !DESCRIPTION:
!   Generates a set of ${\bf G}$-vectors used for the Fourier transform of the
!   charge density and potential and sorts them according to length. Integers
!   corresponding to the vectors in lattice coordinates are stored, as well as
!   the map from these integer coordinates to the ${\bf G}$-vector index. A map
!   from the ${\bf G}$-vector set to the standard FFT array structure is also
!   generated. Finally, the number of ${\bf G}$-vectors with magnitude less than
!   {\tt gmaxvr} is determined.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig,jg,ifg,n1
integer i1,i2,i3,j1,j2,j3
real(8) v1(3),v2(3),v3(3)
! allocatable arrays
integer, allocatable :: idx(:),ivg1(:,:)
real(8), allocatable :: vgc1(:,:),gc1(:)
if (gmaxvr < 0.d0) gmaxvr=abs(gmaxvr)*gkmax
! ensure |G| cut-off is at least twice |G+k| cut-off
gmaxvr=max(gmaxvr,2.d0*gkmax+epslat)
! find the G-vector grid sizes
call gridsize(avec,gmaxvr,ngridg,ngtot,intgv)
! number of complex FFT elements for real-complex transforms
n1=ngridg(1)/2+1
nfgrz=n1*ngridg(2)*ngridg(3)
! allocate global G-vector arrays
if (allocated(ivg)) deallocate(ivg)
allocate(ivg(3,ngtot))
if (allocated(ivgig)) deallocate(ivgig)
allocate(ivgig(intgv(1,1):intgv(2,1),intgv(1,2):intgv(2,2), &
 intgv(1,3):intgv(2,3)))
if (allocated(igfft)) deallocate(igfft)
allocate(igfft(ngtot))
if (allocated(igrzf)) deallocate(igrzf)
allocate(igrzf(nfgrz))
if (allocated(vgc)) deallocate(vgc)
allocate(vgc(3,ngtot))
if (allocated(gc)) deallocate(gc)
allocate(gc(ngtot))
! allocate local arrays
allocate(idx(ngtot),ivg1(3,ngtot))
allocate(vgc1(3,ngtot),gc1(ngtot))
ig=0
do i1=intgv(1,1),intgv(2,1)
  v1(:)=dble(i1)*bvec(:,1)
  do i2=intgv(1,2),intgv(2,2)
    v2(:)=v1(:)+dble(i2)*bvec(:,2)
    do i3=intgv(1,3),intgv(2,3)
      v3(:)=v2(:)+dble(i3)*bvec(:,3)
      ig=ig+1
! map from G-vector to (i1,i2,i3) index
      ivg1(1,ig)=i1
      ivg1(2,ig)=i2
      ivg1(3,ig)=i3
! G-vector in Cartesian coordinates
      vgc1(:,ig)=v3(:)
! length of each G-vector
      gc1(ig)=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
    end do
  end do
end do
! sort by vector length
call sortidx(ngtot,gc1,idx)
! reorder arrays
do ig=1,ngtot
  jg=idx(ig)
  ivg(1:3,ig)=ivg1(1:3,jg)
  gc(ig)=gc1(jg)
  vgc(1:3,ig)=vgc1(1:3,jg)
end do
! find the number of vectors with |G| < gmaxvr
ngvec=1
do ig=2,ngtot
  if (gc(ig) > gmaxvr) then
    ngvec=ig-1
    exit
  end if
end do
! generate index arrays
do ig=1,ngtot
  i1=ivg(1,ig)
  i2=ivg(2,ig)
  i3=ivg(3,ig)
! map from (i1,i2,i3) to G-vector index
  ivgig(i1,i2,i3)=ig
! Fourier transform index
  if (i1 >= 0) then
    j1=i1
  else
    j1=ngridg(1)+i1
  end if
  if (i2 >= 0) then
    j2=i2
  else
    j2=ngridg(2)+i2
  end if
  if (i3 >= 0) then
    j3=i3
  else
    j3=ngridg(3)+i3
  end if
  igfft(ig)=j3*ngridg(2)*ngridg(1)+j2*ngridg(1)+j1+1
! map from G-vector index to real-complex FFT index and vice versa
  if (i1 >= 0) then
    ifg=j3*ngridg(2)*n1+j2*n1+j1+1
    igrzf(ifg)=ig
  end if
end do
deallocate(idx,ivg1,vgc1,gc1)
end subroutine
!EOC

