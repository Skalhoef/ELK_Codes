
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwbandstr
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ip,iw
real(8) dw
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:)
real(8), allocatable :: bmt(:,:,:),bir(:,:)
real(8), allocatable :: wr(:),sf(:)
complex(8), allocatable :: se(:,:,:)
! store original parameters
vkloff0(:)=vkloff(:)
! initialise universal variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! generate k-points along a path for band structure plots
call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! compute the matrix elements of -V_xc and -B_xc
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
if (spinpol) then
  allocate(bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag))
end if
call gwlocal(vmt,vir,bmt,bir)
! real axis frequencies
allocate(wr(nwplot))
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  wr(iw)=dw*dble(iw-1)+wplot(1)
end do
allocate(sf(nwplot),se(nstsv,nstsv,0:nwfm))
if (mp_mpi) then
  open(85,file='GWBAND.OUT',form='FORMATTED',action='WRITE')
  write(85,'(2I6," : grid size")') nwplot,npp1d
end if
! loop over plot points along path
do ip=ip01d,npp1d
  if (mp_mpi) then
    write(*,'("Info(gwbandstr): ",I6," of ",I6," plot points")') ip,npp1d
  end if
! reset the OpenMP thread variables
  call omp_reset
! change the k-point offset
  vkloff(:)=vplp1d(:,ip)*ngridk(:)
! generate the new k-point set
  call init1
! determine the Kohn-Sham ground-state for this k-point offset
  call linengy
  call genapwlofr
  call gensocfr
  call genevfsv
  call occupy
! write the momentum matrix elements to file
  call genpmat
! generate the inverse dielectric function and write to file
  call epsinv
! determine the self-energy for the first k-point
  if (mp_mpi) then
    write(*,'("Info(gwbandstr): calculating self-energy for first k-point")')
  end if
  call gwsefmk(1,vmt,vir,bmt,bir,se)
! solve the Dyson equation on the real axis
  call dysonr(1,wr,se,sf)
  if (mp_mpi) then
    do iw=1,nwplot
      write(85,'(3G18.10)') dpp1d(ip),wr(iw),sf(iw)
    end do
    flush(85)
  end if
! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)
end do
deallocate(vmt,vir,wr,sf,se)
if (spinpol) deallocate(bmt,bir)
if (mp_mpi) then
  close(85)
  write(*,*)
  write(*,'("Info(gwbandstr):")')
  write(*,'(" GW spectral function band structure written to GWBAND.OUT")')
end if
! restore original input parameters
vkloff(:)=vkloff0(:)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

