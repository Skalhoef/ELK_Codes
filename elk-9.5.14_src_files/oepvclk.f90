
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvclk(ikp,vclcv,vclvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vclcv(ncrmax,natmtot,nstsv)
complex(8), intent(out) :: vclvv(nstsv,nstsv)
! local variables
integer ik,jk,nst,ist1,ist2,ist3
integer is,ia,ias,nrc,nrci,npc
integer iv(3),ig,iq,ic,jc,m1,m2
real(8) vc(3)
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(4), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(4), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(4), allocatable :: wfcr1(:,:),wfcr2(:,:)
complex(4), allocatable :: crhomt1(:,:,:),crhomt2(:,:),crhoir1(:,:)
complex(4), allocatable :: cvclmt(:,:),cvclir(:)
! external functions
complex(8), external :: zcfinp,zcfmtinp
! allocate local arrays
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(wfcr1(npcmtmax,2),wfcr2(npcmtmax,2))
allocate(crhomt1(npcmtmax,natmtot,nstsv),crhoir1(ngtc,nstsv))
allocate(crhomt2(npcmtmax,nstcr))
allocate(cvclmt(npcmtmax,natmtot),cvclir(ngtc))
! zero the Coulomb matrix elements
vclcv(:,:,:)=0.d0
vclvv(:,:)=0.d0
! get the eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(filext,ikp,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),vgkc(:,:,1,ikp),gkc(:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states of the input k-point
call genwfsv_sp(.false.,.false.,nstsv,[0],ngdgc,igfc,ngk(1,ikp),igkig(:,1,ikp),&
 apwalm,evecfv,evecsv,wfmt1,ngtc,wfir1)
! loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine the q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
! check if the q-point is in user-defined set
  iv(:)=iv(:)*ngridq(:)
  if (any(mod(iv(:),ngridk(:)) /= 0)) cycle
  iv(:)=iv(:)/ngridk(:)
  iq=ivqiq(iv(1),iv(2),iv(3))
  vc(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvc
! determine the G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length
    gqc(ig)=sqrt(vgqc(1,ig)**2+vgqc(2,ig)**2+vgqc(3,ig)**2)
! spherical harmonics for G+q-vectors
    call genylmv(lmaxo,vgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvc,vgqc,ngvc,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngvc,gqc,gclgq)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ngvc,gqc,ngvc,jlgqrmt)
! find the matching coefficients
  call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-points
  call getevecfv(filext,0,vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(filext,0,vkl(:,ik),evecsv)
! count and index occupied states
  nst=0
  do ist3=1,nstsv
    if (evalsv(ist3,jk) > efermi) cycle
    nst=nst+1
    idx(nst)=ist3
  end do
! calculate the wavefunctions for occupied states
  call genwfsv_sp(.false.,.false.,nst,idx,ngdgc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
  do ist3=1,nst
! compute the complex overlap densities for all valence-valence states
    do ist1=1,nstsv
      call gencrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist3),wfir2(:,:,ist3), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),crhomt1(:,:,ist1),crhoir1(:,ist1))
    end do
! compute the complex overlap densities for all valence-core states
    jc=0
    do is=1,nspecies
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ist1=1,nstsp(is)
          if (spcore(ist1,is)) then
            do m1=-ksp(ist1,is),ksp(ist1,is)-1
              jc=jc+1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
              call wavefcr(.false.,lradstp,is,ia,ist1,m1,npcmtmax,wfcr1)
              if (spinpol) then
                call crho2(npc,wfmt2(:,ias,1,ist3),wfmt2(:,ias,2,ist3),wfcr1, &
                 wfcr1(:,2),crhomt2(:,jc))
              else
                call crho1(npc,wfmt2(:,ias,1,ist3),wfcr1,crhomt2(:,jc))
              end if
! convert to spherical harmonics
              call cfshtip(nrc,nrci,crhomt2(:,jc))
            end do
          end if
        end do
      end do
    end do
    do ist2=1,nstsv
      if (evalsv(ist2,ikp) > efermi) then
! calculate the Coulomb potential
        call gencvclmt(nrcmt,nrcmti,nrcmtmax,rlcmt,wprcmt,npcmtmax, &
         crhomt1(:,:,ist2),cvclmt)
        call cpotcoul(nrcmt,nrcmti,npcmt,nrcmtmax,rlcmt,ngdgc,igfc,ngvc,gqc, &
         gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,crhoir1(:,ist2),npcmtmax,cvclmt,cvclir)
        cvclir(:)=cvclir(:)*cfrc(:)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
        do ist1=1,nstsv
          if (evalsv(ist1,ikp) < efermi) then
            z1=zcfinp(crhomt1(:,:,ist1),crhoir1(:,ist1),cvclmt,cvclir)
            vclvv(ist1,ist2)=vclvv(ist1,ist2)-wqptnr*z1
          end if
        end do
!-------------------------------------------!
!     core-valence-valence contribution     !
!-------------------------------------------!
        jc=0
        do is=1,nspecies
          nrc=nrcmt(is)
          nrci=nrcmti(is)
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            ic=0
            do ist1=1,nstsp(is)
              if (spcore(ist1,is)) then
                do m1=-ksp(ist1,is),ksp(ist1,is)-1
                  ic=ic+1
                  jc=jc+1
                  z1=zcfmtinp(nrc,nrci,wrcmt(:,is),crhomt2(:,jc),cvclmt(:,ias))
                  vclcv(ic,ias,ist2)=vclcv(ic,ias,ist2)-wqptnr*z1
                end do
! end loop over ist1
              end if
            end do
! end loops over atoms and species
          end do
        end do
! end loop over ist2
      end if
    end do
! end loop over ist3
  end do
! end loop over non-reduced k-point set
end do
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m1=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m1,npcmtmax,wfcr1)
! compute the complex overlap densities for the core-valence states
          do ist1=1,nstsv
            if (spinpol) then
              call crho2(npc,wfcr1,wfcr1(:,2),wfmt1(:,ias,1,ist1), &
               wfmt1(:,ias,2,ist1),crhomt1(:,ias,ist1))
            else
              call crho1(npc,wfcr1,wfmt1(:,ias,1,ist1),crhomt1(:,ias,ist1))
            end if
            call cfshtip(nrc,nrci,crhomt1(:,ias,ist1))
          end do
! compute the complex overlap densities for the core-core states
          ic=0
          do ist1=1,nstsp(is)
            if (spcore(ist1,is)) then
              do m2=-ksp(ist1,is),ksp(ist1,is)-1
                ic=ic+1
                call wavefcr(.false.,lradstp,is,ia,ist1,m2,npcmtmax,wfcr2)
                call crho2(npc,wfcr1,wfcr1(:,2),wfcr2,wfcr2(:,2),crhomt2(:,ic))
                call cfshtip(nrc,nrci,crhomt2(:,ic))
              end do
            end if
          end do
          do ist2=1,nstsv
            if (evalsv(ist2,ikp) > efermi) then
! calculate the Coulomb potential
              call cpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
               crhomt1(:,ias,ist2),cvclmt)
!-------------------------------------------!
!     valence-core-valence contribution     !
!-------------------------------------------!
              do ist1=1,nstsv
                if (evalsv(ist1,ikp) < efermi) then
                  z1=zcfmtinp(nrc,nrci,wrcmt(:,is),crhomt1(:,ias,ist1),cvclmt)
                  vclvv(ist1,ist2)=vclvv(ist1,ist2)-z1
                end if
              end do
!----------------------------------------!
!     core-core-valence contribution     !
!----------------------------------------!
              ic=0
              do ist1=1,nstsp(is)
                if (spcore(ist1,is)) then
                  do m2=-ksp(ist1,is),ksp(ist1,is)-1
                    ic=ic+1
                    z1=zcfmtinp(nrc,nrci,wrcmt(:,is),crhomt2(:,ic),cvclmt)
                    vclcv(ic,ias,ist2)=vclcv(ic,ias,ist2)-z1
                  end do
! end loop over ist1
                end if
              end do
! end loop over ist2
            end if
          end do
! end loops over ist3 and m1
        end do
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(vgqc,gqc,gclgq,jlgqrmt)
deallocate(apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfir1,wfmt2,wfir2,wfcr1,wfcr2)
deallocate(crhomt1,crhomt2,crhoir1)
deallocate(cvclmt,cvclir)
return

contains

pure subroutine crho1(n,wf1,wf2,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf1(n),wf2(n)
complex(4), intent(out) :: crho(n)
crho(:)=conjg(wf1(:))*wf2(:)
end subroutine

pure subroutine crho2(n,wf11,wf12,wf21,wf22,crho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(4), intent(out) :: crho(n)
crho(:)=conjg(wf11(:))*wf21(:)+conjg(wf12(:))*wf22(:)
end subroutine

end subroutine
!EOC

