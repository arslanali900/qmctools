!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   subroutine core_charge_ftr( tpre )
!=----------------------------------------------------------------------------=!
     !
     !  Compute the fourier trasform of the core charge, from the radial
     !  mesh to the reciprocal space
     !
     use kinds,              ONLY : DP
     use ions_base,          ONLY : nsp
     use atom,               ONLY : rgrid
     use uspp_param,         ONLY : upf
     use gvecb,              ONLY : ngb, gb
     use small_box,          ONLY : omegab, tpibab
     use pseudo_base,        ONLY : compute_rhocg
     use cp_interfaces,      ONLY : build_cctab, chkpstab
     use pseudopotential,    ONLY : tpstab, rhoc1_sp, rhocp_sp
     use cell_base,          ONLY : omega, tpiba2, tpiba
     USE splines,            ONLY : spline
     use reciprocal_vectors, ONLY : ngm, g, gstart
     USE core,               ONLY : rhocb, rhocg, drhocg, nlcc_any
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: tpre
     !
     INTEGER :: is, ig
     REAL(DP) :: xg, cost1
     !
     !
     IF( .NOT. nlcc_any ) RETURN
     !
     IF( .NOT. ALLOCATED( rgrid ) ) &
        CALL errore( ' core_charge_ftr ', ' rgrid not allocated ', 1 )
     IF( .NOT. ALLOCATED( upf ) ) &
        CALL errore( ' core_charge_ftr ', ' upf not allocated ', 1 )
     !
     IF( tpstab ) THEN
        !
        CALL build_cctab( )
        !
     END IF
     !
     do is = 1, nsp
        !
        if( upf(is)%nlcc ) then
           !
              CALL compute_rhocg( rhocb(:,is), rhocb(:,is), rgrid(is)%r, &
                  rgrid(is)%rab, upf(is)%rho_atc(:), gb, omegab, tpibab**2, &
                  rgrid(is)%mesh, ngb, 0 )
              !
           IF( tpre ) THEN
              !
              IF( tpstab ) THEN
                 !
                 cost1 = 1.0d0/omega
                 !
                 IF( gstart == 2 ) THEN
                    rhocg (1,is) = rhoc1_sp(is)%y( 1 ) * cost1
                    drhocg(1,is) = 0.0d0
                 END IF
                 DO ig = gstart, SIZE( rhocg, 1 )
                    xg = SQRT( g(ig) ) * tpiba
                    rhocg (ig,is) = spline( rhoc1_sp(is), xg ) * cost1
                    drhocg(ig,is) = spline( rhocp_sp(is), xg ) * cost1
                 END DO
                 !
              ELSE

                 CALL compute_rhocg( rhocg(:,is), drhocg(:,is), rgrid(is)%r, &
                                     rgrid(is)%rab, upf(is)%rho_atc(:), g, &
                                     omega, tpiba2, rgrid(is)%mesh, ngm, 1 )

              END IF
              !
           END IF
           !
        endif
        !
     end do

     return
   end subroutine core_charge_ftr



!-----------------------------------------------------------------------
      subroutine add_cc( rhoc, rhog, rhor )
!-----------------------------------------------------------------------
!
! add core correction to the charge density for exch-corr calculation
!
      USE kinds,              ONLY: DP
      use electrons_base,     only: nspin
      use control_flags,      only: iprsta
      use io_global,          only: stdout
      use mp_global,          only: intra_image_comm
      use cell_base,          only: omega
      use recvecs_indexes,    only: np
      USE mp,                 ONLY: mp_sum

      ! this isn't really needed, but if I remove it, ifc 7.1
      ! gives an "internal compiler error"
      use reciprocal_vectors, only: gstart
      use gvecp,              only: ngm
      use grid_dimensions,    only: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnrx
      USE cp_interfaces,      ONLY: fwfft
      USE fft_base,           ONLY: dfftp
!
      implicit none
      !
      REAL(DP),    INTENT(IN)   :: rhoc( nnrx )
      REAL(DP),    INTENT(INOUT):: rhor( nnrx, nspin )
      COMPLEX(DP), INTENT(INOUT):: rhog( ngm,  nspin )
      !
      COMPLEX(DP), ALLOCATABLE :: wrk1( : )
!
      integer :: ig, ir, iss, isup, isdw
      REAL(DP) :: rsum
      !
      IF( iprsta > 2 ) THEN
         rsum = SUM( rhoc ) * omega / DBLE(nr1*nr2*nr3)
         CALL mp_sum( rsum, intra_image_comm )
         WRITE( stdout, 10 ) rsum 
10       FORMAT( 3X, 'Core Charge = ', D14.6 )
      END IF
      !
      ! In r-space:
      !
      if ( nspin .eq. 1 ) then
         iss=1
         call daxpy(nnrx,1.d0,rhoc,1,rhor(1,iss),1)
      else
         isup=1
         isdw=2
         call daxpy(nnrx,0.5d0,rhoc,1,rhor(1,isup),1)
         call daxpy(nnrx,0.5d0,rhoc,1,rhor(1,isdw),1)
      end if 
      !
      ! rhoc(r) -> rhoc(g)  (wrk1 is used as work space)
      !
      allocate( wrk1( nnrx ) )

      wrk1(:) = rhoc(:)

      call fwfft('Dense',wrk1, dfftp )
      !
      ! In g-space:
      !
      if (nspin.eq.1) then
         do ig=1,ngm
            rhog(ig,iss)=rhog(ig,iss)+wrk1(np(ig))
         end do
      else
         do ig=1,ngm
            rhog(ig,isup)=rhog(ig,isup)+0.5d0*wrk1(np(ig))
            rhog(ig,isdw)=rhog(ig,isdw)+0.5d0*wrk1(np(ig))
         end do
      end if

      deallocate( wrk1 )
!
      return
      end subroutine add_cc


!
!-----------------------------------------------------------------------
      subroutine force_cc(irb,eigrb,vxc,fion1)
!-----------------------------------------------------------------------
!
!     core correction force: f = \int V_xc(r) (d rhoc(r)/d R_i) dr
!     same logic as in newd - uses box grid. For parallel execution:
!     the sum over node contributions is done in the calling routine
!
      USE kinds,           ONLY: DP
      use electrons_base,  only: nspin
      use gvecb,           only: gxb, ngb, npb, nmb
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx
      use cell_base,       only: omega
      use ions_base,       only: nsp, na, nat
      use small_box,       only: tpibab
      use uspp_param,      only: upf
      use core,            only: rhocb
      use cp_interfaces,   only: invfft
      use fft_base,        only: dfftb
      use reciprocal_vectors, only: gstart
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none

! input
      integer, intent(in)        :: irb(3,nat)
      complex(8), intent(in):: eigrb(ngb,nat)
      real(8), intent(in)   :: vxc(nnr,nspin)
! output
      real(8), intent(inout):: fion1(3,nat)
! local
      integer iss, ix, ig, is, ia, nfft, isa
      real(8) fcc(3,nat), fac, boxdotgrid
      complex(8) ci, facg
      complex(8), allocatable :: qv(:)
      external  boxdotgrid
!
      call start_clock( 'forcecc' )
      ci = (0.d0,1.d0)
      fac = omega/DBLE(nr1*nr2*nr3*nspin)
      fcc = 0.d0

      allocate( qv( nnrb ) )

      isa = 0

      do is=1,nsp
         if( .not. upf(is)%nlcc ) go to 10
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            if ( dfftb%np3( ia + isa ) <= 0 ) go to 15
#else
         do ia=1,na(is),2
!
! two fft's on two atoms at the same time (when possible)
!
            nfft=2
            if(ia.eq.na(is)) nfft=1
#endif
            do ix=1,3
               qv(:) = (0.d0, 0.d0)
               if (nfft.eq.2) then
                  do ig=1,ngb
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa  )*facg                 &
     &                      + ci * eigrb(ig,ia+isa+1)*facg
                     qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa  )*facg)          &
     &                      + ci * CONJG(eigrb(ig,ia+isa+1)*facg)
                  end do
               else
                  do ig=1,ngb
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                     qv(npb(ig)) = eigrb(ig,ia+isa)*facg
                     qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa)*facg)
                  end do
               end if
!
               call invfft('Box',qv,dfftb,ia+isa)
               !
               ! note that a factor 1/2 is hidden in fac if nspin=2
               !
               do iss=1,nspin
                  fcc(ix,ia+isa) = fcc(ix,ia+isa) + fac *               &
     &                 boxdotgrid(irb(1,ia  +isa),1,qv,vxc(1,iss))
                  if (nfft.eq.2)                                         &
     &               fcc(ix,ia+1+isa) = fcc(ix,ia+1+isa) + fac *           &
     &                    boxdotgrid(irb(1,ia+1+isa),2,qv,vxc(1,iss))
               end do
            end do
15          continue
         end do
10       continue
         isa = isa + na(is)
      end do
!
      do ia = 1, nat
        fion1(:,ia) = fion1(:,ia) + fcc(:,ia)
      end do

      deallocate( qv )
!
      call stop_clock( 'forcecc' )
      return
      end subroutine force_cc


!
!-----------------------------------------------------------------------
      subroutine set_cc( irb, eigrb, rhoc )
!-----------------------------------------------------------------------
!
!     Calculate core charge contribution in real space, rhoc(r)
!     Same logic as for rhov: use box grid for core charges
!
      use ions_base,       only: nsp, na, nat
      use uspp_param,      only: upf
      use grid_dimensions, only: nr3, nnr => nnrx
      use gvecb,           only: ngb, npb, nmb
      use control_flags,   only: iprint
      use core,            only: rhocb
      use cp_interfaces,   only: invfft
      use fft_base,        only: dfftb
      use smallbox_grid_dimensions, only: nr1b, nr2b, nr3b, &
            nr1bx, nr2bx, nr3bx, nnrb => nnrbx

      implicit none
! input
      integer, intent(in)        :: irb(3,nat)
      complex(8), intent(in):: eigrb(ngb,nat)
! output
      real(8), intent(out)  :: rhoc(nnr)
! local
      integer nfft, ig, is, ia, isa
      complex(8) ci
      complex(8), allocatable :: wrk1(:)
      complex(8), allocatable :: qv(:)
!
      call start_clock( 'set_cc' )
      ci=(0.d0,1.d0)
!
      allocate( qv ( nnrb ) )
      allocate( wrk1 ( nnr ) )
      wrk1 (:) = (0.d0, 0.d0)
!
      isa = 0
      do is=1,nsp
         if (.not.upf(is)%nlcc) go to 10
#ifdef __PARA
         do ia=1,na(is)
            nfft=1
            if ( dfftb%np3( ia + isa ) <= 0 ) go to 15
#else
         do ia=1,na(is),2
            nfft=2
            if( ia.eq.na(is) ) nfft=1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
#endif
            qv(:) = (0.d0, 0.d0)
            if(nfft.eq.2)then
               do ig=1,ngb
                  qv(npb(ig))= eigrb(ig,ia  +isa)*rhocb(ig,is)          &
     &                    + ci*eigrb(ig,ia+1+isa)*rhocb(ig,is)
                  qv(nmb(ig))= CONJG(eigrb(ig,ia  +isa)*rhocb(ig,is))   &
     &                    + ci*CONJG(eigrb(ig,ia+1+isa)*rhocb(ig,is))
               end do
            else
               do ig=1,ngb
                  qv(npb(ig)) = eigrb(ig,ia+isa)*rhocb(ig,is)
                  qv(nmb(ig)) = CONJG(eigrb(ig,ia+isa)*rhocb(ig,is))
               end do
            endif
!
            call invfft('Box',qv,dfftb,isa+ia)
!
            call box2grid(irb(1,ia+isa),1,qv,wrk1)
            if (nfft.eq.2) call box2grid(irb(1,ia+1+isa),2,qv,wrk1)
!
15          continue
         end do
10       continue
         isa = isa + na(is)
      end do
!
      call dcopy(nnr,wrk1,2,rhoc,1)

      deallocate( qv  )
      deallocate( wrk1 )
!
      call stop_clock( 'set_cc' )
!
      return
   end subroutine set_cc
