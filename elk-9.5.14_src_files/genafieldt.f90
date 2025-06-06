
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genafieldt
! !INTERFACE:
subroutine genafieldt
! !USES:
use modmain
use modtddft
! !DESCRIPTION:
!   Generates a time-dependent vector potential, ${\bf A}(t)$, representing a
!   laser pulse and stores it in {\tt AFIELDT.OUT}. The vector potential is
!   constructed from a sum of sinusoidal waves, each modulated with a Gaussian
!   envelope function:
!   $$ {\bf A}(t)={\bf A}_0
!    \frac{e^{-(t-t_0)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
!    \sin(\omega(t-t_0)+\phi). $$
!   Seven real numbers have to be specified for each pulse, namely the vector
!   amplitude ${\bf A}_0$, peak time $t_0$, full-width at half-maximum
!   $d=2\sqrt{2\ln 2}\sigma$, frequency $\omega$ and phase $\phi$.
!
! !REVISION HISTORY:
!   Created May 2012 (K. Krieger)
!   Modified, January 2014 (S. Sharma)
!   Modified, February 2014 (JKD)
!   Added spin-dependent A-fields, January 2023 (E. Harris-Lee)
!EOP
!BOC
implicit none
! local variables
integer its,i,j
real(8) av0(3),t0,d,w,phi,rc
real(8) ft,ppd,s,t
real(8) av(3),s0,sv(3)
real(8) t1,t2,t3,t4
! conversion factor of power density to W/cm^2
real(8), parameter :: cpd=ha_si/(t_si*(100.d0*br_si)**2)
! generate the time step grid
call gentimes
open(50,file='TD_INFO.OUT',form='FORMATTED')
write(50,*)
write(50,'("(All units are atomic unless otherwise specified)")')
write(50,*)
write(50,'("1 atomic unit of time is ",G18.10," attoseconds")') t_si*1.d18
write(50,*)
write(50,'("Total simulation time : ",G18.10)') tstime
write(50,'(" in attoseconds       : ",G18.10)') tstime*t_si*1.d18
write(50,*)
write(50,'("Time step length : ",G18.10)') dtimes
write(50,'(" in attoseconds  : ",G18.10)') dtimes*t_si*1.d18
write(50,*)
write(50,'("Number of time steps : ",I8)') ntimes
write(50,*)
write(50,'("Number of laser pulses : ",I6)') npulse
write(50,'("Number of ramps : ",I6)') nramp
write(50,'("Number of steps : ",I6)') nstep
! allocate and zero time-dependent A-field array
if (allocated(afieldt)) deallocate(afieldt)
allocate(afieldt(3,ntimes))
afieldt(:,:)=0.d0
! allocate and zero spin- and time-dependent A-field array
if (tafspt) then
  if (allocated(afspt)) deallocate(afspt)
  allocate(afspt(3,3,ntimes))
  afspt(:,:,:)=0.d0
end if
!----------------------!
!     laser pulses     !
!----------------------!
do i=1,npulse
! vector amplitude
  av0(1:3)=pulse(1:3,i)
! frequency
  w=pulse(4,i)
! phase
  phi=pulse(5,i)
! chirp rate
  rc=pulse(6,i)
! peak time
  t0=pulse(7,i)
! full-width at half-maximum
  d=pulse(8,i)
! Gaussian sigma
  s=d/(2.d0*sqrt(2.d0*log(2.d0)))
! write information to TD_INFO.OUT
  write(50,*)
  write(50,'("Pulse : ",I6)') i
  write(50,'(" vector amplitude : ",3G18.10)') av0(:)
  write(50,'(" laser frequency : ",G18.10)') w
  write(50,'("  in eV          : ",G18.10)') w*ha_ev
  write(50,'(" laser wavelength (Angstroms) : ",G18.10)') 1.d10/(w*ha_im)
  write(50,'(" phase (degrees) : ",G18.10)') phi
  write(50,'(" chirp rate : ",G18.10)') rc
  write(50,'(" peak time : ",G18.10)') t0
  write(50,'(" full-width at half-maximum : ",G18.10)') d
  write(50,'(" Gaussian σ = FWHM / 2√(2ln2) : ",G18.10)') s
  t1=av0(1)**2+av0(2)**2+av0(3)**2
  ppd=t1*(w**2)/(8.d0*pi*solsc)
  write(50,'(" peak laser power density : ",G18.10)') ppd
  write(50,'("  in W/cm²                : ",G18.10)') ppd*cpd
  if (tafspt) then
    s0=pulse(9,i)
    sv(1:3)=pulse(10:12,i)
    write(50,'(" spin components for σ_0, σ_x, σ_y, σ_z : ")')
    write(50,'(4G18.10)') s0,sv(:)
  end if
! loop over time steps
  do its=1,ntimes
    t=times(its)
    t1=t-t0
    t2=-0.5d0*(t1/s)**2
    t3=w*t1+phi*pi/180.d0+0.5d0*rc*t**2
    ft=exp(t2)*sin(t3)
    if (abs(ft) < 1.d-20) ft=0.d0
    av(:)=ft*av0(:)
! spin-polarised vector potential
    if (tafspt) then
      do j=1,3
        afspt(:,j,its)=afspt(:,j,its)+av(:)*sv(j)
      end do
      av(:)=s0*av(:)
    end if
    afieldt(:,its)=afieldt(:,its)+av(:)
  end do
end do
!---------------!
!     ramps     !
!---------------!
do i=1,nramp
! vector amplitude
  av0(1:3)=ramp(1:3,i)
! ramp start time
  t0=ramp(4,i)
! linear coefficient
  t1=ramp(5,i)
! quadratic coefficient
  t2=ramp(6,i)
! cubic coefficient
  t3=ramp(7,i)
! quartic coefficient
  t4=ramp(8,i)
! write information to TD_INFO.OUT
  write(50,*)
  write(50,'("Ramp : ",I6)') i
  write(50,'(" vector amplitude : ",3G18.10)') av0(:)
  write(50,'(" ramp start time : ",G18.10)') t0
  write(50,'(" coefficients : ",4G18.10)') t1,t2,t3,t4
  if (tafspt) then
    s0=ramp(9,i)
    sv(1:3)=ramp(10:12,i)
    write(50,'(" spin components for σ_0, σ_x, σ_y, σ_z : ")')
    write(50,'(4G18.10)') s0,sv(:)
  end if
! loop over time steps
  do its=1,ntimes
    t=times(its)-t0
    if (t > 0.d0) then
      ft=t*(t1+t*(t2+t*(t3+t*t4)))
      av(:)=ft*av0(:)
      if (tafspt) then
        do j=1,3
          afspt(:,j,its)=afspt(:,j,its)+av(:)*sv(j)
        end do
        av(:)=s0*av(:)
      end if
      afieldt(:,its)=afieldt(:,its)+av(:)
    end if
  end do
end do
!---------------!
!     steps     !
!---------------!
do i=1,nstep
! vector amplitude
  av(1:3)=step(1:3,i)
! step start time
  t0=step(4,i)-1.d-14
! step stop time
  t1=step(5,i)
! write information to TD_INFO.OUT
  write(50,*)
  write(50,'("Step : ",I6)') i
  write(50,'(" vector amplitude : ",3G18.10)') av(:)
  write(50,'(" step start and stop times : ",2G18.10)') t0,t1
  if (tafspt) then
    s0=step(6,i)
    sv(1:3)=step(7:9,i)
    write(50,'(" spin components for σ_0, σ_x, σ_y, σ_z : ")')
    write(50,'(4G18.10)') s0,sv(:)
  end if
! loop over time steps
  do its=1,ntimes
    t=times(its)
    if (t > t1) exit
    if (t >= t0) then
      afieldt(:,its)=afieldt(:,its)+av(:)
      if (tafspt) then
        do j=1,3
          afspt(:,j,its)=afspt(:,j,its)+av(:)*sv(j)
        end do
      end if
    end if
  end do
end do
close(50)
! write the vector potential to AFIELDT.OUT
open(50,file='AFIELDT.OUT',form='FORMATTED')
write(50,'(I8," : number of time steps")') ntimes
do its=1,ntimes
  write(50,'(I8,4G18.10)') its,times(its),afieldt(:,its)
end do
close(50)
! write the spin-polarised vector potential to AFSPT.OUT
if (tafspt) then
  open(50,file='AFSPT.OUT',form='FORMATTED')
  write(50,'(I8," : number of time steps")') ntimes
  do its=1,ntimes
    write(50,'(I8,10G18.10)') its,times(its),afspt(:,:,its)
  end do
  close(50)
end if
write(*,*)
write(*,'("Info(genafieldt):")')
write(*,'(" Time-dependent A-field written to AFIELDT.OUT")')
if (tafspt) then
  write(*,'(" Time- and spin-dependent A-field written to AFSPT.OUT")')
end if
write(*,'(" Laser pulse, ramp and step parameters written to TD_INFO.OUT")')
write(*,*)
write(*,'(" 1 atomic unit of time is ",G18.10," attoseconds")') t_si*1.d18
write(*,'(" Total simulation time : ",G18.10)') tstime
write(*,'("  in attoseconds       : ",G18.10)') tstime*t_si*1.d18
deallocate(times,afieldt)
end subroutine
!EOC

