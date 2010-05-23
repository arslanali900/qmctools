!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_base
!=----------------------------------------------------------------------=!

        USE kinds, ONLY: DP
        USE parallel_include

        USE fft_types, ONLY: fft_dlay_descriptor

        IMPLICIT NONE

        ! ... data structure containing all information
        ! ... about fft data distribution for a given
        ! ... potential grid, and its wave functions sub-grid.

        TYPE ( fft_dlay_descriptor ) :: dfftp ! descriptor for dense grid
        TYPE ( fft_dlay_descriptor ) :: dffts ! descriptor for smooth grid
        TYPE ( fft_dlay_descriptor ) :: dfftb ! descriptor for box grids

        SAVE

        PRIVATE

        PUBLIC :: fft_scatter, grid_gather, grid_scatter
        PUBLIC :: dfftp, dffts, dfftb, fft_dlay_descriptor
        PUBLIC :: cgather_sym, cgather_smooth, cscatter_sym, cscatter_smooth



!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!
#if defined __NONBLOCKING_FFT
!
!   NON BLOCKING SCATTER, should be better on switched network
!   like infiniband, ethernet, myrinet
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign, use_tg )
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (sign > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nrx3 elements for a fft of order nr3
  !    nrx3 can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !
  !  b) From planes to columns (sign < 0)
  !
  !  Quite the same in the opposite direction
  !
  !  The output is overwritten on f_in ; f_aux is used as work space
  !
  !  If optional argument "use_tg" is true the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
#ifdef __PARA
  USE parallel_include
#endif
  USE mp_global,   ONLY : nproc_pool, me_pool, intra_pool_comm, nproc, &
                          my_image_id, nogrp, pgrp_comm, nplist, me_pgrp, npgrp
  USE kinds,       ONLY : DP

  IMPLICIT NONE

  INTEGER, INTENT(in)           :: nrx3, nxx_, sign, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  LOGICAL, OPTIONAL, INTENT(in) :: use_tg

#ifdef __PARA

  INTEGER :: dest, from, k, ip, proc, ierr, me, ipoffset, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: sendcount(nproc_pool), sdispls(nproc_pool), recvcount(nproc_pool), rdispls(nproc_pool)
  INTEGER :: offset(nproc_pool)
  INTEGER :: sh(nproc_pool), rh(nproc_pool)
  !
  LOGICAL :: use_tg_ , lrcv, lsnd
  LOGICAL :: tsts(nproc_pool), tstr(nproc_pool)
  INTEGER :: istat( MPI_STATUS_SIZE )

#if defined __HPM
     !       CALL f_hpmstart( 10, 'scatter' )
#endif

  !
  !  Task Groups

  use_tg_ = .false.

  IF( present( use_tg ) ) use_tg_ = use_tg

  me     = me_pool + 1
  !
  IF( use_tg_ ) THEN
    !  This is the number of procs. in the plane-wave group
     nprocp   = npgrp
     ipoffset = me_pgrp
     gcomm    = pgrp_comm
  ELSE
     nprocp   = nproc_pool
     ipoffset = me_pool
     gcomm    = intra_pool_comm
  ENDIF
  !
  IF ( nprocp == 1 ) RETURN
  !
  CALL start_clock ('fft_scatter')
  !
  ! sendcount(proc): amount of data processor "me" must send to processor
  ! recvcount(proc): amount of data processor "me" must receive from
  !
  ! offset is used to locate the slices to be sent to proc
  ! sdispls+1 is the beginning of data that must be sent to proc
  ! rdispls+1 is the beginning of data that must be received from pr
  !
  IF( use_tg_ ) THEN
     DO proc = 1, nprocp
        gproc = nplist( proc ) + 1
        sendcount (proc) = npp_ ( gproc ) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ ( gproc )
     ENDDO
     offset(1) = 0
     DO proc = 2, nprocp
        gproc = nplist( proc - 1 ) + 1
        offset(proc) = offset(proc - 1) + npp_ ( gproc )
     ENDDO
  ELSE
     DO proc = 1, nprocp
        sendcount (proc) = npp_ (proc) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ (proc)
     ENDDO
     offset(1) = 0
     DO proc = 2, nprocp
        offset(proc) = offset(proc - 1) + npp_ (proc - 1)
     ENDDO
  ENDIF
  !
  sdispls (1) = 0
  rdispls (1) = 0
  DO proc = 2, nprocp
     sdispls (proc) = sdispls (proc - 1) + sendcount (proc - 1)
     rdispls (proc) = rdispls (proc - 1) + recvcount (proc - 1)
  ENDDO
  !
  ierr = 0
  !
  IF ( sign > 0 ) THEN
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices and send
     !
     DO ip = 1, nprocp

        ! the following two lines make the loop iterations different on each
        ! proc in order to avoid that all procs send a msg at the same proc
        ! at the same time.
        !
        proc = ipoffset + 1 + ip
        IF( proc > nprocp ) proc = proc - nprocp

        gproc  = proc
        IF( use_tg_ ) gproc  = nplist( proc ) + 1
        !
        from = 1 + offset( proc )
        dest = 1 + sdispls( proc )
        !
        !  optimize for large parallel execution, where npp_ ( gproc ) ~ 1
        !
        SELECT CASE ( npp_ ( gproc ) )
        CASE ( 1 )
           DO k = 1, ncp_ (me)
              f_aux (dest + (k - 1) ) =  f_in (from + (k - 1) * nrx3 )
           ENDDO
        CASE ( 2 )
           DO k = 1, ncp_ (me)
              f_aux ( dest + (k - 1) * 2 - 1 + 1 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 1 )
              f_aux ( dest + (k - 1) * 2 - 1 + 2 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 2 )
           ENDDO
        CASE ( 3 )
!$omp parallel do
           DO k = 1, ncp_ (me)
              f_aux ( dest + (k - 1) * 3 - 1 + 1 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 1 )
              f_aux ( dest + (k - 1) * 3 - 1 + 2 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 2 )
              f_aux ( dest + (k - 1) * 3 - 1 + 3 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 3 )
           ENDDO
        CASE ( 4 )
!$omp parallel do
           DO k = 1, ncp_ (me)
              f_aux ( dest + (k - 1) * 4 - 1 + 1 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 1 )
              f_aux ( dest + (k - 1) * 4 - 1 + 2 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 2 )
              f_aux ( dest + (k - 1) * 4 - 1 + 3 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 3 )
              f_aux ( dest + (k - 1) * 4 - 1 + 4 ) =  f_in ( from + (k - 1) * nrx3 - 1 + 4 )
           ENDDO
        CASE DEFAULT
!$omp parallel do default(shared), private(i, kdest, kfrom)
           DO k = 1, ncp_ (me)
              kdest = dest + (k - 1) * npp_ ( gproc ) - 1
              kfrom = from + (k - 1) * nrx3 - 1
              DO i = 1, npp_ ( gproc )
                 f_aux ( kdest + i ) =  f_in ( kfrom + i )
              ENDDO
           ENDDO
        END SELECT
        !
        ! post the non-blocking send, f_aux can't be overwritten until operation has completed
        !
        CALL mpi_isend( f_aux( sdispls( proc ) + 1 ), sendcount( proc ), MPI_DOUBLE_COMPLEX, &
             proc-1, me, gcomm, sh( proc ), ierr )
        !
        IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', ' forward send info<>0', abs(ierr) )
        !
        !
     ENDDO
     !
     ! step two: receive
     !
     DO ip = 1, nprocp
        !
        proc = ipoffset + 1 - ip
        IF( proc < 1 ) proc = proc + nprocp
        !
        ! now post the receive
        !
        CALL mpi_irecv( f_in( rdispls( proc ) + 1 ), recvcount( proc ), MPI_DOUBLE_COMPLEX, &
             proc-1, MPI_ANY_TAG, gcomm, rh( proc ), ierr )
        !
        IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', ' forward receive info<>0', abs(ierr) )
        !
        tstr( proc )  = .false.
        tsts( proc )  = .false.
        !
     ENDDO
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in( rdispls( nprocp ) + recvcount( nprocp ) + 1 : size( f_in )  ) = 0.0_DP
     !
     lrcv = .false.
     lsnd = .false.
     !
     ! exit only when all test are true: message operation have completed
     !
     DO WHILE ( .not. lrcv .or. .not. lsnd )
        lrcv = .true.
        lsnd = .true.
        DO proc = 1, nprocp
           !
           IF( .not. tstr( proc ) ) THEN
              CALL mpi_test( rh( proc ), tstr( proc ), istat, ierr )
           ENDIF
           !
           IF( .not. tsts( proc ) ) THEN
              CALL mpi_test( sh( proc ), tsts( proc ), istat, ierr )
           ENDIF
           !
           lrcv = lrcv .and. tstr( proc )
           lsnd = lsnd .and. tsts( proc )
           !
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     DO ip = 1, nprocp

        !  post the non blocking send

        proc = ipoffset + 1 + ip
        IF( proc > nprocp ) proc = proc - nprocp

        CALL mpi_isend( f_in( rdispls( proc ) + 1 ), recvcount( proc ), MPI_DOUBLE_COMPLEX, &
             proc-1, me, gcomm, sh( proc ), ierr )
        IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', ' backward send info<>0', abs(ierr) )

        !  post the non blocking receive

        proc = ipoffset + 1 - ip
        IF( proc < 1 ) proc = proc + nprocp

        CALL mpi_irecv( f_aux( sdispls( proc ) + 1 ), sendcount( proc ), MPI_DOUBLE_COMPLEX, &
             proc-1, MPI_ANY_TAG, gcomm, rh(proc), ierr )
        IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', ' backward receive info<>0', abs(ierr) )

        tstr( ip )  = .false.
        tsts( ip )  = .false.

     ENDDO
     !
     lrcv = .false.
     lsnd = .false.
     !
     ! exit only when all test are true: message hsve been sent and received
     !
     DO WHILE ( .not. lsnd )
        !
        lsnd = .true.
        !
        DO proc = 1, nprocp
           !
           IF( .not. tsts( proc ) ) THEN
              CALL mpi_test( sh( proc ), tsts( proc ), istat, ierr )
           ENDIF

           lsnd = lsnd .and. tsts( proc )

        ENDDO

     ENDDO
     !
     lrcv = .false.
     !
     DO WHILE ( .not. lrcv )
        !
        lrcv = .true.
        !
        DO proc = 1, nprocp

           gproc = proc
           IF( use_tg_ ) gproc = nplist(proc)+1

           IF( .not. tstr( proc ) ) THEN

              CALL mpi_test( rh( proc ), tstr( proc ), istat, ierr )

              IF( tstr( proc ) ) THEN

                 from = 1 + sdispls( proc )
                 dest = 1 + offset( proc )
                 !
                 !  optimize for large parallel execution, where npp_ ( gproc ) ~ 1
                 !
                 SELECT CASE ( npp_ ( gproc ) )
                 CASE ( 1 )
                    DO k = 1, ncp_ (me)
                       f_in ( dest + (k - 1) * nrx3 ) = f_aux ( from + k - 1 )
                    ENDDO
                 CASE ( 2 )
                    DO k = 1, ncp_ ( me )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 1 ) = f_aux( from + (k - 1) * 2 - 1 + 1 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 2 ) = f_aux( from + (k - 1) * 2 - 1 + 2 )
                    ENDDO
                 CASE ( 3 )
!$omp parallel do
                    DO k = 1, ncp_ ( me )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 1 ) = f_aux( from + (k - 1) * 3 - 1 + 1 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 2 ) = f_aux( from + (k - 1) * 3 - 1 + 2 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 3 ) = f_aux( from + (k - 1) * 3 - 1 + 3 )
                    ENDDO
                 CASE ( 4 )
!$omp parallel do
                    DO k = 1, ncp_ ( me )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 1 ) = f_aux( from + (k - 1) * 4 - 1 + 1 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 2 ) = f_aux( from + (k - 1) * 4 - 1 + 2 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 3 ) = f_aux( from + (k - 1) * 4 - 1 + 3 )
                       f_in ( dest + (k - 1) * nrx3 - 1 + 4 ) = f_aux( from + (k - 1) * 4 - 1 + 4 )
                    ENDDO
                 CASE DEFAULT
!$omp parallel do default(shared), private(i, kdest, kfrom)
                    DO k = 1, ncp_ ( me )
                       kdest = dest + (k - 1) * nrx3 - 1
                       kfrom = from + (k - 1) * npp_ ( gproc ) - 1
                       DO i = 1, npp_ ( gproc )
                          f_in ( kdest + i ) = f_aux( kfrom + i )
                       ENDDO
                    ENDDO
                 END SELECT

              ENDIF

           ENDIF

           lrcv = lrcv .and. tstr( proc )

        ENDDO

     ENDDO

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

#if defined __HPM
     !       CALL f_hpmstop( 10 )
#endif

  RETURN

END SUBROUTINE fft_scatter
!
!
!
#else
!
!   ALLTOALL based SCATTER, should be better on network
!   with a defined topology, like on bluegene and cray machine
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign, use_tg )
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (sign > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nrx3 elements for a fft of order nr3
  !    nrx3 can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !
  !  b) From planes to columns (sign < 0)
  !
  !  Quite the same in the opposite direction
  !
  !  The output is overwritten on f_in ; f_aux is used as work space
  !
  !  If optional argument "use_tg" is true the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
#ifdef __PARA
  USE parallel_include
#endif
  USE mp_global,   ONLY : nproc_pool, me_pool, intra_pool_comm, nproc, &
                          my_image_id, nogrp, pgrp_comm, nplist
  USE kinds,       ONLY : DP

  IMPLICIT NONE

  INTEGER, INTENT(in)           :: nrx3, nxx_, sign, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  LOGICAL, OPTIONAL, INTENT(in) :: use_tg

#ifdef __PARA

  INTEGER :: dest, from, k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: sendcount (nproc_pool), sdispls (nproc_pool), recvcount (nproc_pool), rdispls (nproc_pool)
  !
  LOGICAL :: use_tg_

#if defined __HPM
     !       CALL f_hpmstart( 10, 'scatter' )
#endif

  !
  !  Task Groups

  use_tg_ = .false.

  IF( present( use_tg ) ) use_tg_ = use_tg

  me     = me_pool + 1
  !
  IF( use_tg_ ) THEN
    !  This is the number of procs. in the plane-wave group
     nprocp = nproc_pool / nogrp
  ELSE
     nprocp = nproc_pool
  ENDIF
  !
  IF (nprocp.eq.1) RETURN
  !
  CALL start_clock ('fft_scatter')
  !
  ! sendcount(proc): amount of data processor "me" must send to processor
  ! recvcount(proc): amount of data processor "me" must receive from
  ! offset1(proc) is used to locate the slices to be sent to proc
  ! sdispls(proc)+1 is the beginning of data that must be sent to proc
  ! rdispls(proc)+1 is the beginning of data that must be received from pr
  !
  !
  IF( use_tg_ ) THEN
     DO proc = 1, nprocp
        gproc = nplist( proc ) + 1
        sendcount (proc) = npp_ ( gproc ) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ ( gproc )
     ENDDO
  ELSE
     DO proc = 1, nprocp
        sendcount (proc) = npp_ (proc) * ncp_ (me)
        recvcount (proc) = npp_ (me) * ncp_ (proc)
     ENDDO
  ENDIF
  !
  sdispls (1) = 0
  rdispls (1) = 0
  DO proc = 2, nprocp
     sdispls (proc) = sdispls (proc - 1) + sendcount (proc - 1)
     rdispls (proc) = rdispls (proc - 1) + recvcount (proc - 1)
  ENDDO
  !

  ierr = 0
  IF (sign.gt.0) THEN
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 1

     DO proc = 1, nprocp
        from = offset
        dest = 1 + sdispls (proc)
        IF( use_tg_ ) THEN
           gproc = nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
        !
        !  optimize for large parallel execution, where npp_ ( gproc ) ~ 1
        !
        IF( npp_ ( gproc ) > 128 ) THEN
           DO k = 1, ncp_ (me)
              CALL DCOPY (2 * npp_ ( gproc ), f_in (from + (k - 1) * nrx3), &
                   1, f_aux (dest + (k - 1) * npp_ ( gproc ) ), 1)
           ENDDO
        ELSEIF( npp_ ( gproc ) == 1 ) THEN
           DO k = 1, ncp_ (me)
              f_aux (dest + (k - 1) ) =  f_in (from + (k - 1) * nrx3 )
           ENDDO
        ELSE
           DO k = 1, ncp_ (me)
              kdest = dest + (k - 1) * npp_ ( gproc ) - 1
              kfrom = from + (k - 1) * nrx3 - 1
              DO i = 1, npp_ ( gproc )
                 f_aux ( kdest + i ) =  f_in ( kfrom + i )
              ENDDO
           ENDDO
        ENDIF
        offset = offset + npp_ ( gproc )
     ENDDO
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in = 0.0_DP
     !
     ! step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = pgrp_comm
     ELSE
        gcomm = intra_pool_comm
     ENDIF

     CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

     CALL mpi_alltoallv (f_aux(1), sendcount, sdispls, MPI_DOUBLE_COMPLEX, f_in(1), &
          recvcount, rdispls, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', 'info<>0', abs(ierr) )
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     !  step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = pgrp_comm
     ELSE
        gcomm = intra_pool_comm
     ENDIF

     CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

     CALL mpi_alltoallv (f_in(1), recvcount, rdispls, MPI_DOUBLE_COMPLEX, f_aux(1), &
          sendcount, sdispls, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     f_in = 0.0_DP
     !
     offset = 1
     !
     DO proc = 1, nprocp
        from = 1 + sdispls (proc)
        dest = offset
        IF( use_tg_ ) THEN
           gproc = nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
        !
        !  optimize for large parallel execution, where npp_ ( gproc ) ~ 1
        !
        IF( npp_ ( gproc ) > 128 ) THEN
           DO k = 1, ncp_ (me)
              CALL DCOPY ( 2 * npp_ ( gproc ), f_aux (from + (k - 1) * npp_ ( gproc ) ), 1, &
                                            f_in  (dest + (k - 1) * nrx3 ), 1 )
           ENDDO
        ELSEIF ( npp_ ( gproc ) == 1 ) THEN
           DO k = 1, ncp_ (me)
              f_in ( dest + (k - 1) * nrx3 ) = f_aux ( from + (k - 1) )
           ENDDO
        ELSE
           DO k = 1, ncp_ (me)
              kdest = dest + (k - 1) * nrx3 - 1
              kfrom = from + (k - 1) * npp_ ( gproc ) - 1
              DO i = 1, npp_ ( gproc )
                 f_in ( kdest + i ) = f_aux( kfrom + i )
              ENDDO
           ENDDO
        ENDIF
        !
        offset = offset + npp_ ( gproc )
        !
     ENDDO

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

#if defined __HPM
     !       CALL f_hpmstop( 10 )
#endif

  RETURN

END SUBROUTINE fft_scatter

#endif

!----------------------------------------------------------------------------
SUBROUTINE grid_gather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers nproc_pool distributed data on the first processor of every pool
  !
  ! ... REAL*8  f_in  = distributed variable (nxx)
  ! ... REAL*8  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
  USE kinds,     ONLY : DP
  USE parallel_include
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in( : ), f_out( : )
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  IF( size( f_in ) < dfftp%nnr ) &
     CALL errore( ' grid_gather ', ' f_in too small ', dfftp%nnr - size( f_in ) )
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  info = size( f_out ) - displs( nproc_pool - 1 ) - recvcount( nproc_pool - 1 )
  !
  IF( info < 0 ) &
     CALL errore( ' grid_gather ', ' f_out too small ', -info )
  !
  info = 0
  !
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, root_pool,    &
                    intra_pool_comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE grid_gather


!----------------------------------------------------------------------------
SUBROUTINE grid_scatter( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... REAL*8  f_in  = gathered variable (nrx1*nrx2*nrx3)
  ! ... REAL*8  f_out = distributed variable (nxx)
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in( : ), f_out( : )
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  IF( size( f_out ) < dfftp%nnr ) &
     CALL errore( ' grid_scatter ', ' f_out too small ', dfftp%nnr - size( f_in ) )
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  info = size( f_in ) - displs( nproc_pool - 1 ) - sendcount( nproc_pool - 1 )
  !
  IF( info < 0 ) &
     CALL errore( ' grid_scatter ', ' f_in too small ', -info )
  !
  info = 0
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /= dfftp%nnr ) f_out(sendcount(me_pool)+1:dfftp%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE grid_scatter
!
! ... "gather"-like subroutines
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (in phonon code)
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
  USE mp_global, ONLY : intra_pool_comm, intra_image_comm, &
                        nproc_pool, me_pool
  USE mp,        ONLY : mp_barrier
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in( : ), f_out(:)
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_ALLGATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       intra_pool_comm, info )
  !
  CALL errore( 'cgather_sym', 'info<>0', info )
  !
!  CALL mp_barrier( intra_image_comm )
  !
  CALL stop_clock( 'cgather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cgather_smooth ( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the smooth AND complex fft grid
  !
  ! ... gathers nproc_pool distributed data on the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = distributed variable ( dffts%nnr )
  ! ... COMPLEX*16  f_out = gathered variable (nrx1s*nrx2s*nrx3s)
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, root_pool,    &
                    intra_pool_comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_smooth
!
! ... "scatter"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nrx1*nrx2*nrx3)
  ! ... COMPLEX*16  f_out = distributed variable (nxx)
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'cscatter_sym' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'cscatter_sym', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /=  dfftp%nnr  ) f_out(sendcount(me_pool)+1: dfftp%nnr ) = 0.D0
  !
  CALL stop_clock( 'cscatter_sym' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_smooth( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the smooth AND complex fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nrx1s*nrx2s*nrx3s)
  ! ... COMPLEX*16  f_out = distributed variable ( dffts%nnr)
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__PARA)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /=  dffts%nnr  ) f_out(sendcount(me_pool)+1: dffts%nnr ) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_smooth

!=----------------------------------------------------------------------=!
   END MODULE fft_base
!=----------------------------------------------------------------------=!
