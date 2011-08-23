
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genshtmat
! !INTERFACE:
subroutine genshtmat
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the forward and backward spherical harmonic transformation (SHT)
!   matrices using the spherical covering set produced by the routine
!   {\tt sphcover}. These matrices are used to transform a function between its
!   $(l,m)$-expansion coefficients and its values at the $(\theta,\phi)$ points
!   on the sphere.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer itp,info
! automatic arrays
integer ipiv(lmmaxvr)
real(8) tp(2,lmmaxvr),rlm(lmmaxvr),work(lmmaxvr)
complex(8) ylm(lmmaxvr),zwork(lmmaxvr)
! allocate real SHT matrices for lmaxvr
if (allocated(rbshtvr)) deallocate(rbshtvr)
allocate(rbshtvr(lmmaxvr,lmmaxvr))
if (allocated(rfshtvr)) deallocate(rfshtvr)
allocate(rfshtvr(lmmaxvr,lmmaxvr))
! allocate complex SHT matrices for lmaxvr
if (allocated(zbshtvr)) deallocate(zbshtvr)
allocate(zbshtvr(lmmaxvr,lmmaxvr))
if (allocated(zfshtvr)) deallocate(zfshtvr)
allocate(zfshtvr(lmmaxvr,lmmaxvr))
! generate spherical covering set for lmaxvr
call sphcover(lmmaxvr,tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxvr
  call genrlm(lmaxvr,tp(:,itp),rlm)
  rbshtvr(itp,1:lmmaxvr)=rlm(1:lmmaxvr)
  call genylm(lmaxvr,tp(:,itp),ylm)
  zbshtvr(itp,1:lmmaxvr)=ylm(1:lmmaxvr)
end do
! find the forward SHT arrays
! real
rfshtvr(:,:)=rbshtvr(:,:)
call dgetrf(lmmaxvr,lmmaxvr,rfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) goto 10
call dgetri(lmmaxvr,rfshtvr,lmmaxvr,ipiv,work,lmmaxvr,info)
if (info.ne.0) goto 10
! complex
zfshtvr(:,:)=zbshtvr(:,:)
call zgetrf(lmmaxvr,lmmaxvr,zfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) goto 10
call zgetri(lmmaxvr,zfshtvr,lmmaxvr,ipiv,zwork,lmmaxvr,info)
if (info.ne.0) goto 10
return
10 continue
write(*,*)
write(*,'("Error(genshtmat): unable to find inverse spherical harmonic &
 &transform")')
write(*,'(" => improper spherical covering")')
write(*,*)
stop
end subroutine
!EOC

