
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine alpha2f
use modmain
use modphonon
use modtest
implicit none
! local variables
integer ik,iq,i,j
integer i1,i2,i3,iw
real(8) wmin,wmax,wd,dw
real(8) wlog,wrms,lambda,tc
real(8) v(3),t1
! allocatable arrays
real(8), allocatable :: gq(:,:),wq(:),a2fp(:),w(:),a2f(:)
complex(8), allocatable :: dynq(:,:,:),dynr(:,:,:),dynp(:,:)
complex(8), allocatable :: a2fmq(:,:,:),a2fmr(:,:,:),a2fmp(:,:)
complex(8), allocatable :: ev(:,:),b(:,:)
! initialise universal variables
call init0
call init1
call init2
! allocate local arrays
allocate(gq(nbph,nqpt),wq(nbph),a2fp(nbph))
allocate(w(nwplot),a2f(nwplot))
allocate(dynq(nbph,nbph,nqpt))
allocate(dynr(nbph,nbph,nqptnr))
allocate(dynp(nbph,nbph))
allocate(a2fmq(nbph,nbph,nqpt))
allocate(a2fmr(nbph,nbph,nqptnr))
allocate(a2fmp(nbph,nbph))
allocate(ev(nbph,nbph),b(nbph,nbph))
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
end do
! compute the density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! read in the phonon linewidths for each q-point
call readgamma(gq)
! loop over phonon q-points
do iq=1,nqpt
! find the eigenvalues and eigenvectors of the dynamical matrix
  call dynev(dynq(:,:,iq),wphq(:,iq),ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues squared are the phonon linewidths divided by the frequency
  do i=1,nbph
    if (wphq(i,iq) > 1.d-8) then
      t1=sqrt(abs(gq(i,iq)/wphq(i,iq)))
    else
      t1=0.d0
    end if
    do j=1,nbph
      b(i,j)=t1*conjg(ev(j,i))
    end do
  end do
  call zgemm('N','N',nbph,nbph,nbph,zone,ev,nbph,b,nbph,zzero,a2fmq(:,:,iq), &
   nbph)
end do
! Fourier transform the matrices to real-space
call dynqtor(a2fmq,a2fmr)
! find the minimum and maximum frequencies
wmin=0.d0
wmax=0.d0
do iq=1,nqpt
  wmin=min(wmin,wphq(1,iq))
  wmax=max(wmax,wphq(nbph,iq))
end do
t1=(wmax-wmin)*0.1d0
wmin=wmin-t1
wmax=wmax+t1
wd=wmax-wmin
if (wd < 1.d-8) wd=1.d0
dw=wd/dble(nwplot)
! generate energy grid
do iw=1,nwplot
  w(iw)=dw*dble(iw-1)+wmin
end do
a2f(:)=0.d0
do i1=0,ngrkf-1
  v(1)=dble(i1)/dble(ngrkf)
  do i2=0,ngrkf-1
    v(2)=dble(i2)/dble(ngrkf)
    do i3=0,ngrkf-1
      v(3)=dble(i3)/dble(ngrkf)
! compute the dynamical matrix at this particular q-point
      call dynrtoq(v,dynr,dynp)
! find the phonon eigenvalues and eigenvectors
      call dynev(dynp,wq,ev)
! compute the α²F matrix at this particular q-point
      call dynrtoq(v,a2fmr,a2fmp)
! diagonalise the α²F matrix simultaneously with the dynamical matrix
! (thanks to Matthieu Verstraete and Ryotaro Arita for correcting this)
      call dynevs(ev,a2fmp,a2fp)
! square the eigenvalues to recover the linewidths divided by the frequency
      a2fp(:)=a2fp(:)**2
      do i=1,nbph
        t1=(wq(i)-wmin)/dw+1.d0
        iw=nint(t1)
        if ((iw >= 1).and.(iw <= nwplot)) then
          a2f(iw)=a2f(iw)+a2fp(i)
        end if
      end do
    end do
  end do
end do
t1=twopi*(fermidos/2.d0)*dw*dble(ngrkf)**3
if (t1 > 1.d-8) then
  t1=1.d0/t1
else
  t1=0.d0
end if
a2f(:)=t1*a2f(:)
! smooth Eliashberg function if required
if (nswplot > 0) call fsmooth(nswplot,nwplot,a2f)
! write Eliashberg function to file
open(50,file='ALPHA2F.OUT',form='FORMATTED')
do iw=1,nwplot
  write(50,'(2G18.10)') w(iw),a2f(iw)
end do
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Eliashberg function α²F written to ALPHA2F.OUT")')
! compute lambda, logarithmic average frequency, RMS average frequency and
! McMillan-Allen-Dynes superconducting critical temperature
call mcmillan(w,a2f,lambda,wlog,wrms,tc)
open(50,file='MCMILLAN.OUT',form='FORMATTED')
write(50,*)
write(50,'("Electron-phonon coupling constant, λ : ",G18.10)') lambda
write(50,*)
write(50,'("Logarithmic average frequency : ",G18.10)') wlog
write(50,*)
write(50,'("RMS average frequency : ",G18.10)') wrms
write(50,*)
write(50,'("Coulomb pseudopotential, μ* : ",G18.10)') mustar
write(50,*)
write(50,'("McMillan-Allen-Dynes superconducting critical temperature")')
write(50,'(" [Eq. 34, Phys. Rev. B 12, 905 (1975)] (kelvin) : ",G18.10)') tc
write(50,*)
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Electron-phonon coupling constant, λ,")')
write(*,'(" logarithmic and RMS average frequencies,")')
write(*,'(" and McMillan-Allen-Dynes superconducting critical temperature")')
write(*,'(" written to MCMILLAN.OUT")')
! write lambda to test file
call writetest(251,'electron-phonon coupling constant, lambda',tol=5.d-2, &
 rv=lambda)
deallocate(wq,gq,a2fp,w,a2f)
deallocate(dynq,dynr,dynp)
deallocate(a2fmq,a2fmr,a2fmp)
deallocate(ev,b)
end subroutine

