
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvcl(vclcv,vclvv)
use modmain
use modmpi
use modomp
implicit none
! arguments
complex(8), intent(out) :: vclcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(out) :: vclvv(nstsv,nstsv,nkpt)
! local variables
integer ik,ncv,nvv
integer lp,nthd
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
!$OMP CRITICAL(oepvcl_)
  write(*,'("Info(oepvcl): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(oepvcl_)
  call oepvclk(ik,vclcv(:,:,:,ik),vclvv(:,:,ik))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! broadcast matrix elements to all other processes
ncv=ncrmax*natmtot*nstsv
nvv=nstsv*nstsv
do ik=1,nkpt
  lp=mod(ik-1,np_mpi)
  call mpi_bcast(vclcv(:,:,:,ik),ncv,mpi_double_complex,lp,mpicom,ierror)
  call mpi_bcast(vclvv(:,:,ik),nvv,mpi_double_complex,lp,mpicom,ierror)
end do
end subroutine

