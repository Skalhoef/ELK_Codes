
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdvs(iq,is,ia,ip,dvsmt,dvsir)
use modmain
implicit none
! arguments
integer, intent(in) :: iq,is,ia,ip
complex(8), intent(out) :: dvsmt(npmtmax,natmtot),dvsir(ngtot)
! local variables
integer js,jas,ios
integer version_(3),nspecies_
integer lmmaxo_,natoms_
integer nrmt_,ngridg_(3)
character(256) fext,fname
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfmt(lmmaxo,nrmtmax,natmtot))
call dynfext(iq,is,ia,ip,fext)
fname='DVS'//trim(fext)
open(150,file=trim(fname),form='UNFORMATTED',action='READ',status='OLD', &
 iostat=ios)
if (ios /= 0) then
  write(*,*)
  write(*,'("Error(readdvs): error opening ",A)') trim(fname)
  write(*,*)
  stop
end if
read(150) version_
if ((version(1) /= version_(1)).or.(version(2) /= version_(2)) &
 .or.(version(3) /= version_(3))) then
  write(*,*)
  write(*,'("Warning(readdvs): different versions")')
  write(*,'(" current : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" file    : ",I3.3,".",I3.3,".",I3.3)') version_
  write(*,'(" in file ",A)') trim(fname)
end if
read(150) nspecies_
if (nspecies /= nspecies_) then
  write(*,*)
  write(*,'("Error(readdvs): differing nspecies")')
  write(*,'(" current : ",I4)') nspecies
  write(*,'(" file    : ",I4)') nspecies_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(150) lmmaxo_
if (lmmaxo /= lmmaxo_) then
  write(*,*)
  write(*,'("Error(readdvs): differing lmmaxo")')
  write(*,'(" current : ",I4)') lmmaxo
  write(*,'(" file    : ",I4)') lmmaxo_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
do js=1,nspecies
  read(150) natoms_
  if (natoms(js) /= natoms_) then
    write(*,*)
    write(*,'("Error(readdvs): differing natoms for species ",I4)') js
    write(*,'(" current : ",I4)') natoms(js)
    write(*,'(" file    : ",I4)') natoms_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
  read(150) nrmt_
  if (nrmt(js) /= nrmt_) then
    write(*,*)
    write(*,'("Error(readdvs): differing nrmt for species ",I4)') js
    write(*,'(" current : ",I6)') nrmt(js)
    write(*,'(" file    : ",I6)') nrmt_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
end do
read(150) ngridg_
if ((ngridg(1) /= ngridg_(1)).or.(ngridg(2) /= ngridg_(2)).or. &
 (ngridg(3) /= ngridg_(3))) then
  write(*,*)
  write(*,'("Error(readdvs): differing ngridg")')
  write(*,'(" current : ",3I6)') ngridg
  write(*,'(" file    : ",3I6)') ngridg_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(150) zfmt,dvsir
do jas=1,natmtot
  js=idxis(jas)
  call zfmtpack(.true.,nrmt(js),nrmti(js),zfmt(:,:,jas),dvsmt(:,jas))
end do
close(150)
deallocate(zfmt)
end subroutine

