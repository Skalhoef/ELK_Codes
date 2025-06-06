
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdisp
use modmain
use modphonon
implicit none
! local variables
integer iq,i,iv
real(8) wmin,wmax,t1
! allocatable arrays
complex(8), allocatable :: dynq(:,:,:),dynr(:,:,:)
complex(8), allocatable :: dynp(:,:),ev(:,:)
! initialise universal variables
call init0
call init2
! allocate local arrays
allocate(dynq(nbph,nbph,nqpt))
allocate(dynr(nbph,nbph,nqptnr))
allocate(dynp(nbph,nbph))
allocate(ev(nbph,nbph))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! generate a set of q-point vectors along a path in the Brillouin zone
call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
wmin=0.d0
wmax=0.d0
! compute the phonon frequencies along the path
do iq=ip01d,npp1d
! compute the dynamical matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),dynr,dynp)
! find the phonon frequencies and eigenvectors
  call dynev(dynp,wphq(:,iq),ev)
  wmin=min(wmin,wphq(1,iq))
  wmax=max(wmax,wphq(nbph,iq))
end do
t1=(wmax-wmin)*0.5d0
wmin=wmin-t1
wmax=wmax+t1
! output the vertex location lines
open(50,file='PHDLINES.OUT',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),wmin
  write(50,'(2G18.10)') dvp1d(iv),wmax
  write(50,*)
end do
close(50)
! output the phonon dispersion
open(50,file='PHDISP.OUT',form='FORMATTED')
do i=1,nbph
  do iq=ip01d,npp1d
    write(50,'(2G18.10)') dpp1d(iq),wphq(i,iq)
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(phdisp):")')
write(*,'(" phonon dispersion written to PHDISP.OUT")')
write(*,'(" vertex location lines written to PHDLINES.OUT")')
deallocate(dynq,dynr,dynp,ev)
end subroutine

