!
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM pw2qmcpack
  !----------------------------------------------------------------------- 

  ! This subroutine writes the file "prefix".pwscf.xml and "prefix".pwscf.h5
  ! containing the  plane wave coefficients and other stuff needed by QMCPACK. 

#include "f_defs.h"

  USE io_files,  ONLY : nd_nmbr, prefix, outdir, tmp_dir, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, npool, nimage, nogrp, npgrp
  USE environment,ONLY : environment_start
  !
  IMPLICIT NONE
  INTEGER :: ios
  LOGICAL :: write_psir

  NAMELIST / inputpp / prefix, outdir, write_psir
#ifdef __PARA
  CALL mp_startup ( )
#endif
CALL environment_start ( 'pw2qmcpack' )
  IF ( npool > 1 .or. nimage > 1) THEN
     CALL errore('pw2qmcpack', 'pool or image parallelization not (yet) implemented',1)
  ENDIF
!   CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  write_psir = .true.
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  ios = 0
  IF ( ionode )  THEN 
     !
     !READ (5, inputpp, err=200, iostat=ios)
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  CALL mp_bcast( ios, ionode_id ) 
  IF ( ios/=0 ) CALL errore('pw2qmcpack', 'reading inputpp namelist', ABS(ios))
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast(prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  CALL mp_bcast(write_psir, ionode_id ) 
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL compute_qmcpack(write_psir)
  !
  CALL stop_pp
  STOP

END PROGRAM pw2qmcpack


SUBROUTINE compute_qmcpack(write_psir)

  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba2, at, bg
  USE printout_base, ONLY: title
  USE constants, ONLY: tpi
  USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, g, gg, ecutwfc, gcutm, nl, igtongl
  USE klist , ONLY: nks, nelec, nelup, neldw, xk
  USE lsda_mod, ONLY: lsda, nspin
  !!V3.0
  !!USE scf, ONLY: rho, rho_core, rhog, rhog_core
  !!USE vlocal, ONLY: vloc, vnew, strf
  !!USE wvfct, ONLY: npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  USE scf, ONLY: rho, rho_core, rhog_core, vnew
  USE vlocal, ONLY: vloc, strf
  USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
!   USE becmod,   ONLY: becp, calbec
  USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_global, ONLY: stdout
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc, iun => iunsat, tmp_dir, prefix
  USE wavefunctions_module, ONLY : evc, psic
  use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  use iotk_module
  use iotk_xtox_interf
  USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY: mp_sum
  USE dfunct, ONLY : newd
  

  IMPLICIT NONE
  LOGICAL :: write_psir
  INTEGER :: ig, ibnd, ik, io, na, j, ispin, nbndup, nbnddown, &
       nk, ngtot, ig7, ikk, nt, ijkb0, ikb, ih, jh, jkb, at_num, &
       nelec_tot, nelec_up, nelec_down, ii, igx, igy, igz, n_rgrid(3)
  INTEGER, ALLOCATABLE :: INDEX(:), igtog(:), igtomin(:)
  LOGICAL :: exst, found
  REAL(DP) :: ek, eloc, enl, charge, etotefield
  REAL(DP) :: bg_qmc(3,3), g_red(3), lattice_real(3,3)
  COMPLEX(DP), ALLOCATABLE :: phase(:),aux(:), hpsi(:,:), eigpacked(:)
  COMPLEX(DP), ALLOCATABLE :: psitr(:)
  REAL(DP), ALLOCATABLE ::  tau_r(:,:), g_cart(:,:),psireal(:),eigval(:)
  INTEGER :: ios, ierr, h5len,oldh5,ig_c,save_complex, nup,ndown
  INTEGER, EXTERNAL :: atomic_number, is_complex
  REAL(DP), ALLOCATABLE :: g_qmc(:,:)
  INTEGER, ALLOCATABLE :: gint_den(:,:), gint_qmc(:,:)
  REAL (DP), EXTERNAL :: ewald

  CHARACTER(256)          :: tmp,h5name,eigname
  CHARACTER(iotk_attlenx) :: attr

  CALL init_us_1
  CALL newd
  !####io = 77

  !####WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')

  !####CALL seqopn( 77, 'pwfn.data', 'formatted',exst)  

  !ALLOCATE (hpsi(npwx, nbnd))
  ALLOCATE (aux(nrxx))
  CALL allocate_bec_type ( nkb, nbnd, becp )
!   ALLOCATE (becp (nkb,nbnd))
  ! four times npwx should be enough
  ALLOCATE (INDEX (4*npwx) )
  ALLOCATE (igtog (4*npwx) )
  ALLOCATE (igtomin(4*npwx) )


  !hpsi (:,:) = (0.d0, 0.d0)
  INDEX(:) = 0
  igtog(:) = 0
  igtomin(:) = 0

  IF( lsda ) THEN
     nbndup = nbnd
     nbnddown = nbnd
     nk = nks/2
     !     nspin = 2
  ELSE
     nbndup = nbnd
     nbnddown = 0
     nk = nks
     !     nspin = 1
  ENDIF

  !  if(nks > 1) rewind(iunigk)
  !  do ik=1,nks
  !     if(nks > 1) read(iunigk) npw, igk
  !     
  !  if(nks > 1) rewind(iunigk)
  ek  = 0.d0
  eloc= 0.d0
  enl = 0.d0
  !
  DO ispin = 1, nspin 
     !
     !     calculate the local contribution to the total energy
     !
     !      bring rho to G-space
     !
     !aux(:) = CMPLX ( rho(:,ispin), 0.d0)
     aux(:) = CMPLX ( rho%of_r(:,ispin), 0.d0)
     CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     DO nt=1,ntyp
        DO ig = 1, ngm
           eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                * CONJG(aux(nl(ig)))
        ENDDO
     ENDDO

     DO ik = 1, nk
        ikk = ik + nk*(ispin-1)
        CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, ikk, - 1)
        CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
        !CALL ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
        CALL calbec ( npw, vkb, evc, becp )

        DO ig =1, npw
           IF( igk(ig) > 4*npwx ) & 
                CALL errore ('pw2qmcpack','increase allocation of index', ig)
           INDEX( igk(ig) ) = 1
        ENDDO
     ENDDO
  ENDDO

!#ifdef __PARA
!  CALL reduce(1,eloc)
!  CALL reduce(1,ek)
!  CALL poolreduce(1,ek)
!  CALL poolreduce(1,enl)
!#endif
#ifdef __PARA
call mp_sum( eloc,  intra_pool_comm )
call mp_sum( ek,    intra_pool_comm )
call mp_sum( ek,    inter_pool_comm )
call mp_sum( enl,   inter_pool_comm )
!call mp_sum( demet, inter_pool_comm )
#endif

  eloc = eloc * omega 
  ek = ek * tpiba2

  ngtot = 0
! igtomin maps indices from the full set of G-vectors to the
! minimal set which includes the G-spheres of all k-points
  DO ig = 1, 4*npwx
     IF( INDEX(ig) == 1 ) THEN
        ngtot = ngtot + 1
        igtog(ngtot) = ig
        igtomin(ig) = ngtot
     ENDIF
  ENDDO

  ALLOCATE (g_qmc(3,ngtot))
  ALLOCATE (gint_qmc(3,ngtot))
  ALLOCATE (gint_den(3,ngm))
  ALLOCATE (g_cart(3,ngtot))
  ALLOCATE (tau_r(3,nat))

  ! get the number of electrons
  nelec_tot= NINT(nelec)
  nup=NINT(nelup)
  ndown=NINT(neldw)

  if(nup .eq. 0) then
    ndown=nelec_tot/2
    nup=nelec_tot-ndown
  endif

  h5name = TRIM( prefix ) // '.pwscf.h5'
  eigname = "eigenstates_"//trim(iotk_itoa(nr1s))//'_'//trim(iotk_itoa(nr2s))//'_'//trim(iotk_itoa(nr3s))

  tmp = TRIM( tmp_dir )//TRIM( prefix ) //'.pwscf.h5'
  h5len = LEN_TRIM(tmp)

  ! writing to xml and hdf5
  ! open hdf5 file 
  oldh5=0
  CALL esh5_open_file(tmp,h5len,oldh5)

  bg_qmc(:,:)=bg(:,:)/alat

  !! create a file for particle set
  tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.ptcl.xml'
  CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )

  CALL iotk_write_attr (attr,"name","global",first=.true.)
  CALL iotk_write_begin(iun, "simulationcell",ATTR=attr)
  CALL iotk_write_attr (attr,"name","lattice",first=.true.)
  CALL iotk_write_attr (attr,"units","bohr")
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
   
  lattice_real=alat*at
  WRITE(iun,100) lattice_real(1,1), lattice_real(2,1), lattice_real(3,1)
  WRITE(iun,100) lattice_real(1,2), lattice_real(2,2), lattice_real(3,2)
  WRITE(iun,100) lattice_real(1,3), lattice_real(2,3), lattice_real(3,3)

  CALL esh5_write_supercell(lattice_real)

  !WRITE(iun,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
  !WRITE(iun,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
  !WRITE(iun,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_attr (attr,"name","reciprocal",first=.true.)
  CALL iotk_write_attr (attr,"units","2pi/bohr")
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,100) bg_qmc(1,1), bg_qmc(2,1), bg_qmc(3,1)
  WRITE(iun,100) bg_qmc(1,2), bg_qmc(2,2), bg_qmc(3,2)
  WRITE(iun,100) bg_qmc(1,3), bg_qmc(2,3), bg_qmc(3,3)
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_attr (attr,"name","bconds",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,'(a)') 'p p p'
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_attr (attr,"name","LR_dim_cutoff",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  WRITE(iun,'(a)') '15'
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "simulationcell")

  ! <particleset name="ions">
  CALL iotk_write_attr (attr,"name","ion0",first=.true.)
  CALL iotk_write_attr (attr,"size",nat)
  CALL iotk_write_begin(iun, "particleset",ATTR=attr)

  CALL esh5_open_atoms(nat,ntyp)

  ! ionic species --> group
  DO na=1,ntyp

  tmp=TRIM(atm(na))
  h5len=LEN_TRIM(tmp)
  CALL esh5_write_species(na,tmp,h5len,atomic_number(tmp),zv(na))

  CALL iotk_write_attr (attr,"name",TRIM(atm(na)),first=.true.)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) zv(na)
  CALL iotk_write_end(iun, "parameter")

  CALL iotk_write_end(iun, "group")
  ENDDO


  ! <attrib name="ionid"/>
  CALL iotk_write_attr (attr,"name","ionid",first=.true.)
  CALL iotk_write_attr (attr,"datatype","stringArray")
  CALL iotk_write_begin(iun, "attrib",ATTR=attr)
  write(iun,'(a)') (TRIM(atm(ityp(na))),na=1,nat)
  CALL iotk_write_end(iun, "attrib")

  ! <attrib name="position"/>
  CALL iotk_write_attr (attr,"name","position",first=.true.)
  CALL iotk_write_attr (attr,"datatype","posArray")
  CALL iotk_write_attr (attr,"condition","0")
  CALL iotk_write_begin(iun, "attrib",ATTR=attr)
  ! write in cartesian coordinates in bohr
  ! problem with xyz ordering inrelation to real-space grid
  DO na = 1, nat
  !tau_r(1)=alat*(tau(1,na)*bg_qmc(1,1)+tau(2,na)*bg_qmc(1,2)+tau(3,na)*bg_qmc(1,3))
  !tau_r(2)=alat*(tau(1,na)*bg_qmc(2,1)+tau(2,na)*bg_qmc(2,2)+tau(3,na)*bg_qmc(2,3))
  !tau_r(3)=alat*(tau(1,na)*bg_qmc(3,1)+tau(2,na)*bg_qmc(3,2)+tau(3,na)*bg_qmc(3,3))
  tau_r(1,na)=alat*tau(1,na)
  tau_r(2,na)=alat*tau(2,na)
  tau_r(3,na)=alat*tau(3,na)
  WRITE(iun,100) (tau_r(j,na),j=1,3)
  ENDDO
  !write(iun,100) tau
  CALL iotk_write_end(iun, "attrib")
  CALL iotk_write_end(iun, "particleset")

  !cartesian positions
  CALL esh5_write_positions(tau_r)
  CALL esh5_write_species_ids(ityp)

  CALL esh5_close_atoms()
  ! </particleset>

  ! <particleset name="e">
  CALL iotk_write_attr (attr,"name","e",first=.true.)
  CALL iotk_write_attr (attr,"random","yes")
  CALL iotk_write_begin(iun, "particleset",ATTR=attr)

  ! <group name="u" size="" >
  CALL iotk_write_attr (attr,"name","u",first=.true.)
  CALL iotk_write_attr (attr,"size",nup)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) -1
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "group")

  ! <group name="d" size="" >
  CALL iotk_write_attr (attr,"name","d",first=.true.)
  CALL iotk_write_attr (attr,"size",ndown)
  CALL iotk_write_begin(iun, "group",ATTR=attr)
  CALL iotk_write_attr (attr,"name","charge",first=.true.)
  CALL iotk_write_begin(iun, "parameter",ATTR=attr)
  write(iun,*) -1
  CALL iotk_write_end(iun, "parameter")
  CALL iotk_write_end(iun, "group")
  CALL iotk_write_end(iun, "particleset")
  CALL iotk_close_write(iun)

  !! close the file
  DO ik = 0, nk-1
  ! create a xml input file for each k-point
     IF(nk .gt. 1) THEN
       tmp = TRIM( tmp_dir ) // TRIM( prefix ) //TRIM(iotk_index(ik))// '.wfs.xml'
     ELSE
       tmp = TRIM( tmp_dir ) // TRIM( prefix )// '.wfs.xml'
     ENDIF
     CALL iotk_open_write(iun, FILE=TRIM(tmp), ROOT="qmcsystem", IERR=ierr )
     ! <wavefunction name="psi0">
     CALL iotk_write_attr (attr,"name","psi0",first=.true.)
     CALL iotk_write_attr (attr,"target","e")
     CALL iotk_write_begin(iun, "wavefunction",ATTR=attr)
       write(iun,'(a)') '<!-- Uncomment this out to use plane-wave basis functions'
       CALL iotk_write_attr (attr,"type","PW",first=.true.)
       CALL iotk_write_attr (attr,"href",TRIM(h5name))
       CALL iotk_write_attr (attr,"version","1.10")
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
       write(iun,'(a)') '--> '
       CALL iotk_write_attr (attr,"type","bspline",first=.true.)
       CALL iotk_write_attr (attr,"href",TRIM(h5name))
       CALL iotk_write_attr (attr,"sort","1")
       CALL iotk_write_attr (attr,"tilematrix","1 0 0 0 1 0 0 0 1")
       CALL iotk_write_attr (attr,"source","ion0")
       CALL iotk_write_attr (attr,"version","0.10")
       CALL iotk_write_begin(iun, "determinantset",ATTR=attr)
          CALL iotk_write_attr (attr,"ecut",ecutwfc/2,first=.true.)
          ! basisset to overwrite cutoff to a smaller value
          CALL iotk_write_begin(iun, "basisset",ATTR=attr)
             ! add grid to use spline on FFT grid
             CALL iotk_write_attr (attr,"dir","0",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr1s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
             CALL iotk_write_attr (attr,"dir","1",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr2s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
             CALL iotk_write_attr (attr,"dir","2",first=.true.)
             CALL iotk_write_attr (attr,"npts",nr3s)
             CALL iotk_write_attr (attr,"closed","no")
             CALL iotk_write_empty(iun, "grid",ATTR=attr)
          CALL iotk_write_end(iun, "basisset")
          
          !CALL iotk_write_attr (attr,"href",TRIM(h5name),first=.true.)
          !CALL iotk_write_empty(iun, "coefficients",ATTR=attr)
  
          ! write the index of the twist angle
          CALL iotk_write_attr (attr,"name","twistIndex",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          write(iun,*) ik
          CALL iotk_write_end(iun, "h5tag")

          CALL iotk_write_attr (attr,"name","twistAngle",first=.true.)
          CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          g_red(1)=at(1,1)*xk(1,ik+1)+at(2,1)*xk(2,ik+1)+at(3,1)*xk(3,ik+1)
          g_red(2)=at(1,2)*xk(1,ik+1)+at(2,2)*xk(2,ik+1)+at(3,2)*xk(3,ik+1)
          g_red(3)=at(1,3)*xk(1,ik+1)+at(2,3)*xk(2,ik+1)+at(3,3)*xk(3,ik+1)
          !write(iun,100) xk(1,ik+1),xk(2,ik+1),xk(3,ik+1)
          write(iun,100) g_red(1),g_red(2),g_red(3)
          CALL iotk_write_end(iun, "h5tag")
          !write(iun,'(a)') '<!-- Uncomment this out for bspline wavefunctions '
          !!CALL iotk_write_attr (attr,"name","eigenstates",first=.true.)
          !!CALL iotk_write_begin(iun, "h5tag",ATTR=attr)
          !!write(iun,'(a)') TRIM(eigname)
          !!CALL iotk_write_end(iun, "h5tag")
          !write(iun,'(a)') '--> '

  
          CALL iotk_write_begin(iun, "slaterdeterminant")
             ! build determinant for up electrons
             CALL iotk_write_attr (attr,"id","updet",first=.true.)
             CALL iotk_write_attr (attr,"size",nup)
             CALL iotk_write_begin(iun, "determinant",ATTR=attr)
                CALL iotk_write_attr (attr,"mode","ground",first=.true.)
                CALL iotk_write_attr (attr,"spindataset",0)
                CALL iotk_write_begin(iun, "occupation",ATTR=attr)
                CALL iotk_write_end(iun, "occupation")
             CALL iotk_write_end(iun, "determinant")
  
             ! build determinant for down electrons
             CALL iotk_write_attr (attr,"id","downdet",first=.true.)
             CALL iotk_write_attr (attr,"size",ndown)
             IF( lsda ) CALL iotk_write_attr (attr,"ref","updet")
             CALL iotk_write_begin(iun, "determinant",ATTR=attr)
               CALL iotk_write_attr (attr,"mode","ground",first=.true.)
               IF( lsda ) THEN
                 CALL iotk_write_attr (attr,"spindataset",1)
               ELSE
                 CALL iotk_write_attr (attr,"spindataset",0)
               ENDIF
               CALL iotk_write_begin(iun, "occupation",ATTR=attr)
               CALL iotk_write_end(iun, "occupation")
             CALL iotk_write_end(iun, "determinant")
          CALL iotk_write_end(iun, "slaterdeterminant")
  
       CALL iotk_write_end(iun, "determinantset")
     CALL iotk_write_end(iun, "wavefunction")
  
     CALL iotk_close_write(iun)
  ENDDO

  DO ig=1, ngtot
    ig_c =igtog(ig)
    g_cart(1,ig)=tpi/alat*g(1,ig_c)
    g_cart(2,ig)=tpi/alat*g(2,ig_c)
    g_cart(3,ig)=tpi/alat*g(3,ig_c)
    g_qmc(1,ig)=at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c)
    g_qmc(2,ig)=at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c)
    g_qmc(3,ig)=at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c)
    gint_qmc(1,ig)=NINT(at(1,1)*g(1,ig_c)+at(2,1)*g(2,ig_c)+at(3,1)*g(3,ig_c))
    gint_qmc(2,ig)=NINT(at(1,2)*g(1,ig_c)+at(2,2)*g(2,ig_c)+at(3,2)*g(3,ig_c))
    gint_qmc(3,ig)=NINT(at(1,3)*g(1,ig_c)+at(2,3)*g(2,ig_c)+at(3,3)*g(3,ig_c))
  ENDDO
  DO ig=1,ngm
     gint_den(1,ig)=NINT(at(1,1)*g(1,ig)+at(2,1)*g(2,ig)+at(3,1)*g(3,ig))
     gint_den(2,ig)=NINT(at(1,2)*g(1,ig)+at(2,2)*g(2,ig)+at(3,2)*g(3,ig))
     gint_den(3,ig)=NINT(at(1,3)*g(1,ig)+at(2,3)*g(2,ig)+at(3,3)*g(3,ig))
  ENDDO


  n_rgrid(1)=nr1s
  n_rgrid(2)=nr2s
  n_rgrid(3)=nr3s

  save_complex=0
  DO ik = 1, nk
     !! evaluate the phase
     !phase(:) = (0.d0,0.d0)
     !if ( ig_(ik,ib)>0) phase( nls(ig_(ik,ib)) ) = (1.d0,0.d0)
     !call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
     g_red(1)=at(1,1)*xk(1,ik)+at(2,1)*xk(2,ik)+at(3,1)*xk(3,ik)
     g_red(2)=at(1,2)*xk(1,ik)+at(2,2)*xk(2,ik)+at(3,2)*xk(3,ik)
     g_red(3)=at(1,3)*xk(1,ik)+at(2,3)*xk(2,ik)+at(3,3)*xk(3,ik)

     IF(g_red(1)*g_red(1)+g_red(2)*g_red(2)+g_red(3)*g_red(3)>1e-12) THEN
        save_complex=1
     END IF
  END DO

  CALL esh5_open_electrons(nup, ndown,nspin,nk,nbnd,n_rgrid, save_complex)

  !!NOT YET DECIDED
  !!CALL esh5_write_basis(g_qmc,g_cart,ngtot)
  !!CALL esh5_write_parameters(nelec_tot,nspin,nbnd,nk,ecutwfc/2,alat,at)
  !
100 FORMAT (3(1x,f20.15))

  ALLOCATE (eigpacked(ngtot))
  ALLOCATE (eigval(nbnd))

  ! start a main section to save eigen values and vector
  DO ik = 1, nk
     g_red(1)=at(1,1)*xk(1,ik)+at(2,1)*xk(2,ik)+at(3,1)*xk(3,ik)
     g_red(2)=at(1,2)*xk(1,ik)+at(2,2)*xk(2,ik)+at(3,2)*xk(3,ik)
     g_red(3)=at(1,3)*xk(1,ik)+at(2,3)*xk(2,ik)+at(3,3)*xk(3,ik)

     !! open kpoint 
     CALL esh5_open_kpoint(ik)
     CALL esh5_write_kpoint_data(g_red(1), 1.0d0, ngtot, gint_qmc)

     DO ispin = 1, nspin 

        CALL esh5_open_spin(ispin)

        ikk = ik + nk*(ispin-1)
        !!IF( nks > 1 ) THEN
        !!   CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        !!   CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        !!ENDIF
        DO ibnd = 1, nbnd
!!           DO ig=1, ngtot
!!              ! now for all G vectors find the PW coefficient for this k-point
!!              found = .FALSE.
!!              DO ig7 = 1, npw
!!                 IF( igk(ig7) == igtog(ig) )THEN
!!                    eigpacked(ig)=evc(ig7,ibnd)
!!                    found = .TRUE.
!!                    GOTO 17
!!                 ENDIF
!!              ENDDO
!!              ! if can't find the coefficient this is zero
!!17            IF( .NOT. found ) eigpacked(ig)=(0.d0,0.d0) 
!!           ENDDO
!!           CALL esh5_write_psi_g(ibnd,eigpacked,ngtot)
           eigval(ibnd)=0.5*et(ibnd,ikk)
        ENDDO
        CALL esh5_write_eigvalues(eigval)
        CALL esh5_close_spin()
     ENDDO
     CALL esh5_close_kpoint()
  ENDDO

  !!CALL esh5_close_eigg

  !ALLOCATE (phase(nrxxs))

  ALLOCATE(psireal(nrx1s*nrx2s*nrx3s))
  ALLOCATE(psitr(nrx1s*nrx2s*nrx3s))

  ! open real-space wavefunction on FFT grid
  !!CALL esh5_open_eigr(nr1s,nr2s,nr3s)
  DO ik = 1, nk
     CALL esh5_open_kpoint(ik)
     DO ispin = 1, nspin 
        ikk = ik + nk*(ispin-1)
        IF( nks > 1 ) THEN
           CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
           CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
        ENDIF
        CALL esh5_open_spin(ispin)
        DO ibnd = 1, nbnd !!transform G to R
           eigpacked(:)=(0.d0,0.d0)
           eigpacked(igtomin(igk(1:npw)))=evc(1:npw,ibnd)
           CALL esh5_write_psi_g(ibnd,eigpacked,ngtot)
           IF (write_psir) THEN
              psic(:)=(0.d0,0.d0)
              psic(nls(igk(1:npw)))=evc(1:npw,ibnd)
              call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
              

              IF(save_complex .eq. 1) THEN
                 !psic(:)=psic(:)/omega
                 ii=1
                 DO igx=1,nr1s
                    DO igy=0,nr2s-1
                       DO igz=0,nr3s-1
                          psitr(ii)=psic(igx+nr1s*(igy+igz*nr2s))/omega
                          ii=ii+1
                       ENDDO
                    ENDDO
                 ENDDO
                 CALL esh5_write_psi_r(ibnd,psitr,save_complex)
              ELSE
                 DO ig=1,nr1s*nr2s*nr3s
                    psireal(ig)=real(psic(ig))
                 ENDDO
                 psireal(1:nr1s*nr2s*nr3s)=real(psic(1:nr1s*nr2s*nr3s))/omega
                 ii=1
                 DO igx=1,nr1s
                    DO igy=0,nr2s-1
                       DO igz=0,nr3s-1
                          psireal(ii)=real(psic(igx+nr1s*(igy+igz*nr2s)))/omega
                          ii=ii+1
                       ENDDO
                    ENDDO
                 ENDDO
                 CALL esh5_write_psi_r(ibnd,psireal,save_complex)
              ENDIF
           ENDIF
           !! conversion and output complete for each band
        ENDDO
        CALL esh5_close_spin()
     ENDDO
     CALL esh5_close_kpoint()
  ENDDO

  ! write charge density
  ! ignore spin for the time being
  !CALL esh5_write_rho(rho,rhog(1,1),ngm)

  CALL esh5_open_density(gint_den,ngm,nr1s,nr2,nr3s)
  ! CALL esh5_open_density(g_qmc,igtog,ngtot,nr1s,nr2,nr3s)
    DO ispin = 1, nspin
       !! ii=1
       !! DO igx=1,nr1s
       !! DO igy=0,nr2s-1
       !! DO igz=0,nr3s-1
       !! !psireal(ii)=rho(igx+nr1s*(igy+igz*nr2s),1)
       !! psireal(ii)=rho%of_r(igx+nr1s*(igy+igz*nr2s),ispin)
       !! ii=ii+1
       !! ENDDO
       !! ENDDO
       !! ENDDO
       !! CALL esh5_write_rho(psireal,rhog(1,1),ngm)
       !! what are the G????
       CALL esh5_write_density_g(ispin,rho%of_g(1,ispin))
       !! CALL esh5_write_density_r(ispin,psireal)
    ENDDO
  CALL esh5_close_density()
  CALL esh5_close_electrons()
  CALL esh5_close_file()

  DEALLOCATE (igtog)
  DEALLOCATE (index)
!   DEALLOCATE (becp)
  CALL deallocate_bec_type (becp)
  DEALLOCATE (aux)
  DEALLOCATE (eigpacked)
  DEALLOCATE (g_qmc)
  DEALLOCATE (g_cart)
  !DEALLOCATE (psireal)
  !DEALLOCATE (psitr)
  !DEALLOCATE (phase)

END SUBROUTINE compute_qmcpack


