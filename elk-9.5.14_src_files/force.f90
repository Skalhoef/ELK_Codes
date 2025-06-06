
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: force
! !INTERFACE:
subroutine force
! !USES:
use modmain
use modtddft
use modtest
use modmpi
use modomp
! !DESCRIPTION:
!   Computes the various contributions to the atomic forces. In principle, the
!   force acting on a nucleus is simply the gradient at that site of the
!   classical electrostatic potential from the other nuclei and the electronic
!   density. This is a result of the Hellmann-Feynman theorem. However because
!   the basis set is dependent on the nuclear coordinates and is not complete,
!   the Hellman-Feynman force is inaccurate and corrections to it are required.
!   The first is the core correction which arises because the core wavefunctions
!   were determined by neglecting the non-spherical parts of the Kohn-Sham
!   potential $v_s$. Explicitly this is given by
!   $$ {\bf F}_{\rm core}^{\alpha}=\int_{\rm MT_{\alpha}} v_s({\bf r})
!    \nabla\rho_{\rm core}^{\alpha}({\bf r})\,d{\bf r} $$
!   for atom $\alpha$. The second, which is the incomplete basis set (IBS)
!   correction, is due to the position dependence of the APW functions, and is
!   derived by considering the change in total energy if the eigenvector
!   coefficients were fixed and the APW functions themselves were changed. This
!   would result in changes to the first-variational Hamiltonian and overlap
!   matrices given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G-G'})
!    \left(H^{\alpha}_{\bf G+k,G'+k}-\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G-G'})\left(O^{\alpha}_{\bf G+k,G'+k}
!    -\tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)
!   \end{align*}
!   where both ${\bf G}$ and ${\bf G'}$ run over the APW indices;
!   $\tilde{\Theta}_{\alpha}$ is the form factor of the smooth step function for
!   muffin-tin $\alpha$; and $H^{\alpha}$ and $O^{\alpha}$ are the muffin-tin
!   Hamiltonian and overlap matrices, respectively. The APW-local-orbital part
!   is given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G+k})H^{\alpha}_{\bf G+k,G'+k}\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G+k})O^{\alpha}_{\bf G+k,G'+k}
!   \end{align*}
!   where ${\bf G}$ runs over the APW indices and ${\bf G'}$ runs over the
!   local-orbital indices. There is no contribution from the
!   local-orbital-local-orbital part of the matrices. We can now write the IBS
!   correction in terms of the basis of first-variational states as
!   \begin{align*}
!    {\bf F}_{ij}^{\alpha{\bf k}}=\sum_{\bf G,G'}
!    b^{i{\bf k}*}_{\bf G}b^{j{\bf k}}_{\bf G'}\left(
!    \delta H_{\bf G,G'}^{\alpha}-\epsilon_j\delta O_{\bf G,G'}^{\alpha}\right),
!   \end{align*}
!   where $b^{i{\bf k}}$ is the first-variational eigenvector.
!   Finally, the ${\bf F}_{ij}^{\alpha{\bf k}}$ matrix elements can be
!   multiplied by the second-variational coefficients, and contracted over all
!   indices to obtain the IBS force:
!   \begin{align*}
!    {\bf F}_{\rm IBS}^{\alpha}=\sum_{\bf k}w_{\bf k}\sum_{l\sigma}n_{l{\bf k}}
!    \sum_{ij}c_{\sigma i}^{l{\bf k}*}c_{\sigma j}^{l{\bf k}}
!    {\bf F}_{ij}^{\alpha{\bf k}}
!    +\int_{\rm MT_{\alpha}}v_s({\bf r})\nabla\left[\rho({\bf r})
!    -\rho^{\alpha}_{\rm core}({\bf r})\right]\,d{\bf r},
!   \end{align*}
!   where $c^{l{\bf k}}$ are the second-variational coefficients, $w_{\bf k}$
!   are the $k$-point weights, $n_{l{\bf k}}$ are the occupation numbers.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Fixed problem with second-variational forces, May 2008 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,idm,is,ias
integer nr,nri,i,j,nthd
real(8) fav(3),ca,t1,t2
real(8) ts0,ts1
! automatic arrays
real(8) grfmt(npmtmax,3)
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
! external functions
real(8), external :: rfmtinp,rfmtint
call timesec(ts0)
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
!---------------------------------!
!     Hellmann-Feynman forces     !
!---------------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),vclmt(:,ias),npmtmax,grfmt)
! force from Coulomb potential
  forcehf(:,ias)=-spzn(is)*grfmt(1,:)*y00
! force on nuclei from time-dependent E-field
  if (tafieldt) then
    do i=1,3
      t1=(chgsmt(ias,i)+spzn(is))*efieldt(i)
      forcehf(i,ias)=forcehf(i,ias)+t1
    end do
  end if
end do
! symmetrise Hellmann-Feynman forces
call symveca(forcehf)
!----------------------------------!
!     IBS correction to forces     !
!----------------------------------!
! set the IBS forces to zero
forceibs(:,:)=0.d0
! compute k-point dependent contribution to the IBS forces
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi) /= lp_mpi) cycle
  call forcek(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! add IBS forces from each process and redistribute
if (np_mpi > 1) then
  call mpi_allreduce(mpi_in_place,forceibs,3*natmtot,mpi_double_precision, &
   mpi_sum,mpicom,ierror)
end if
if (tafieldt) then
  t2=afieldt(1,itimes)**2+afieldt(2,itimes)**2+afieldt(3,itimes)**2
  t2=t2*ca**2
end if
! integral of Kohn-Sham potential with gradient of density
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),rhomt(:,ias),npmtmax,grfmt)
  do i=1,3
    t1=rfmtinp(nr,nri,wrmt(:,is),vsmt(:,ias),grfmt(:,i))
! remove contribution from gauge correction to the current density
    if (tafieldt) then
      t1=t1-t2*rfmtint(nr,nri,wrmt(:,is),grfmt(:,i))
    end if
    forceibs(i,ias)=forceibs(i,ias)+t1
  end do
end do
! integral of Kohn-Sham magnetic field with magnetisation gradient
if (spinpol) then
  allocate(rfmt(npmtmax,natmtot))
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfsht(nrcmt(is),nrcmti(is),bsmt(:,ias,idm),rfmt(:,ias))
    end do
    call rfmtctof(rfmt)
    do ias=1,natmtot
      is=idxis(ias)
      nr=nrmt(is)
      nri=nrmti(is)
      call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),magmt(:,ias,idm), &
       npmtmax,grfmt)
      do i=1,3
        t1=rfmtinp(nr,nri,wrmt(:,is),rfmt(:,ias),grfmt(:,i))
        forceibs(i,ias)=forceibs(i,ias)+t1
      end do
    end do
  end do
  deallocate(rfmt)
end if
! time-dependent vector potential times integral of the gradient of the current
! density over the muffin-tin
if (tafieldt.and.tjr) then
  do j=1,3
    t1=ca*afieldt(j,itimes)
    do ias=1,natmtot
      is=idxis(ias)
      nr=nrmt(is)
      nri=nrmti(is)
      call gradrfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),jrmt(:,ias,j),npmtmax, &
       grfmt)
      do i=1,3
        t2=rfmtint(nr,nri,wrmt(:,is),grfmt(:,i))
        forceibs(i,ias)=forceibs(i,ias)+t1*t2
      end do
    end do
  end do
end if
! symmetrise IBS forces
call symveca(forceibs)
! total force on each atom
do ias=1,natmtot
  forcetot(1:3,ias)=forcehf(1:3,ias)+forceibs(1:3,ias)
end do
! symmetrise total forces
call symveca(forcetot)
! remove the average force, if required, to prevent translation of atomic basis
if (tfav0) then
  fav(:)=0.d0
  do ias=1,natmtot
    fav(:)=fav(:)+forcetot(:,ias)
  end do
  fav(:)=fav(:)/dble(natmtot)
  do ias=1,natmtot
    forcetot(:,ias)=forcetot(:,ias)-fav(:)
  end do
end if
! zero force on atoms with negative mass
do ias=1,natmtot
  is=idxis(ias)
  if (spmass(is) <= 0.d0) forcetot(1:3,ias)=0.d0
end do
! compute maximum force magnitude over all atoms
forcemax=0.d0
do ias=1,natmtot
  t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
  if (t1 > forcemax) forcemax=t1
end do
! restrict maximum force magnitude if required
if (maxforce >= 0.d0) then
  if (forcemax > maxforce) then
    t1=maxforce/forcemax
    forcetot(1:3,1:natmtot)=t1*forcetot(1:3,1:natmtot)
  end if
end if
call timesec(ts1)
timefor=timefor+ts1-ts0
! write total forces to test file
call writetest(750,'total forces',nv=3*natmtot,tol=1.d-3,rva=forcetot)
end subroutine
!EOC

