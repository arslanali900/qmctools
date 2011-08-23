
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradwf2(ik,evecfv,evecsv,gwf2mt,gwf2ir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
real(8), intent(inout) :: gwf2mt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: gwf2ir(ngrtot)
! local variables
integer ispn,jspn,ist
integer is,ia,ias,i,j,n
integer nr,ir,itp,igk,ifg
real(8) wo,t1
complex(8) zq(2),zt1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: gwfmt(:,:,:)
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: zfft1(:,:)
complex(8), allocatable :: zfft2(:)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(wfmt1(lmmaxvr,nrmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrmtmax,nspinor))
allocate(gwfmt(lmmaxvr,nrmtmax,3))
allocate(zfmt(lmmaxvr,nrmtmax))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
do is=1,nspecies
  nr=nrmt(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! de-phasing factor for spin-spirals
    if (spinsprl.and.ssdph) then
      t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1)=cmplx(cos(t1),sin(t1),8)
      zq(2)=conjg(zq(1))
    end if
    done(:,:)=.false.
    do j=1,nstsv
      wo=wkpt(ik)*occsv(j,ik)
      if (abs(wo).gt.epsocc) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt2(:,:,:)=0.d0
          i=0
          do ispn=1,nspinor
            if (spinsprl) then
              jspn=ispn
            else
              jspn=1
            end if
            do ist=1,nstfv
              i=i+1
              zt1=evecsv(i,j)
              if (spinsprl.and.ssdph) zt1=zt1*zq(ispn)
              if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                if (.not.done(ist,jspn)) then
                  call wavefmt(1,lmaxvr,is,ia,ngk(jspn,ik), &
                   apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr, &
                   wfmt1(:,:,ist,jspn))
                  done(ist,jspn)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt1(:,:,ist,jspn),1,wfmt2(:,:,ispn),1)
              end if
            end do
          end do
        else
! spin-unpolarised wavefunction
          call wavefmt(1,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,j,1),lmmaxvr, &
           wfmt2)
        end if
! compute the gradient of the wavefunction
        do ispn=1,nspinor
          call gradzfmt(lmaxvr,nr,spr(:,is),lmmaxvr,nrmtmax,wfmt2(:,:,ispn), &
           gwfmt)
! convert gradient from spherical harmonics to spherical coordinates
          do i=1,3
            call zgemm('N','N',lmmaxvr,nr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
             gwfmt(:,:,i),lmmaxvr,zzero,zfmt,lmmaxvr)
            do ir=1,nr
              do itp=1,lmmaxvr
                t1=wo*(dble(zfmt(itp,ir))**2+aimag(zfmt(itp,ir))**2)
                gwf2mt(itp,ir,ias)=gwf2mt(itp,ir,ias)+t1
              end do
            end do
          end do
        end do
      end if
    end do
  end do
end do
deallocate(apwalm,wfmt1,wfmt2,gwfmt,zfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfft1(ngrtot,nspinor))
allocate(zfft2(ngrtot))
do j=1,nstsv
  wo=wkpt(ik)*occsv(j,ik)
  if (abs(wo).gt.epsocc) then
    t1=wo/omega
    zfft1(:,:)=0.d0
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        if (spinsprl) then
          jspn=ispn
        else
          jspn=1
        end if
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            do igk=1,ngk(jspn,ik)
              ifg=igfft(igkig(igk,jspn,ik))
              zfft1(ifg,ispn)=zfft1(ifg,ispn)+zt1*evecfv(igk,ist,jspn)
            end do
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      do igk=1,ngk(1,ik)
        ifg=igfft(igkig(igk,1,ik))
        zfft1(ifg,1)=evecfv(igk,j,1)
      end do
    end if
! compute gradient of wavefunction
    do ispn=1,nspinor
      if (spinsprl) then
        jspn=ispn
      else
        jspn=1
      end if
      do i=1,3
        zfft2(:)=0.d0
        do igk=1,ngk(jspn,ik)
          ifg=igfft(igkig(igk,jspn,ik))
          zfft2(ifg)=zi*vgkc(i,igk,jspn,ik)*zfft1(ifg,ispn)
        end do
! Fourier transform gradient to real-space
        call zfftifc(3,ngrid,1,zfft2)
        do ir=1,ngrtot
          gwf2ir(ir)=gwf2ir(ir)+t1*(dble(zfft2(ir))**2+aimag(zfft2(ir))**2)
        end do
      end do
    end do
  end if
end do
deallocate(zfft1,zfft2)
return
end subroutine
