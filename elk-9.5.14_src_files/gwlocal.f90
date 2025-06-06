
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwlocal(vmt,vir,bmt,bir)
use modmain
use modomp
implicit none
! arguments
real(8), intent(out) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(out) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! local variables
integer idm,is,ias,nthd
integer nrc,nrci,npc
! automatic arrays
real(8) rfmt(npcmtmax)
call holdthd(natmtot+1,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,ias,is) &
!$OMP PRIVATE(nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! convert exchange-correlation potential to a coarse radial mesh
  call rfmtftoc(nrc,nrci,vxcmt(:,ias),rfmt)
! negate because V_xc should be removed from the self-energy
  rfmt(1:npc)=-rfmt(1:npc)
! convert to spherical coordinates
  call rbsht(nrc,nrci,rfmt,vmt(:,ias))
! multiply by radial integration weights
  call rfcmtwr(nrc,nrci,wrcmt(:,is),vmt(:,ias))
end do
!$OMP END DO NOWAIT
! negate and multiply the interstitial V_xc by the characteristic function
!$OMP SINGLE
vir(:)=-vxcir(:)*cfunir(:)
!$OMP END SINGLE NOWAIT
! do the same for B_xc in the spin-polarised case
if (spinpol) then
  do idm=1,ndmag
!$OMP DO
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      call rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),rfmt)
      rfmt(1:npc)=-rfmt(1:npc)
      call rbsht(nrc,nrci,rfmt,bmt(:,ias,idm))
      call rfcmtwr(nrc,nrci,wrcmt(:,is),bmt(:,ias,idm))
    end do
!$OMP END DO NOWAIT
  end do
!$OMP SINGLE
  do idm=1,ndmag
    bir(:,idm)=-bxcir(:,idm)*cfunir(:)
  end do
!$OMP END SINGLE
end if
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

