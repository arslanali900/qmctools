
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

function group(name, num)
  character(len=200) group
  character(*) name
  integer      num
  character(len=50) numbuff
  
  write(numbuff,'(I)'), num
  write(group,'(A,"_",A)'),name,trim(adjustl(numbuff))
  return 
end function


subroutine writehdf5
use modmain
use HDF5
use h5lt
implicit none
! local variables
integer ik,ist
integer is,ia,ias,nr,ir,l,ilo,ispn,ord,ig,lm,m,iatom,i,ng,icore
integer, dimension(1) :: n
CHARACTER(LEN=7), PARAMETER :: fname = "lapw.h5"
CHARACTER(LEN=200) group_name, group
!CHARACTER(LEN=20) group_num
real(8) x,t1
! HDF5 file_id, error
integer(HID_T) h5id, radial_group, eig_group, ion_group, twist_group, &
               band_group, species_group, atom_group, l_group, local_group,&
               spin_group, param_group, core_group, state_group
integer(HSIZE_T), dimension(1) :: dim1 
integer(HSIZE_T), dimension(2) :: dim2
integer(HSIZE_T), dimension(3) :: dim3
integer(HSIZE_T), dimension(4) :: dim4
integer(HSIZE_T), dimension(5) :: dim5
integer h5err
real(8),          dimension(1) :: val1
! allocatable arrays
complex(8), allocatable :: evecfv(:,:,:)
real(8),    allocatable :: evalfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
integer,    allocatable :: atom_types(:)
integer,    allocatable :: species_local(:), atom_local(:), &
     l_local(:), m_local(:), ilo_local(:)
! Indices are: (real/imag, g-point, apw order, m, spin)
real(8), allocatable :: apw_match(:,:,:,:,:)
real(8), allocatable :: apw_match_real(:,:,:,:)
real(8), allocatable :: gvecs(:,:)
real(8), allocatable :: atom_pos(:,:)
! Indices are real/complex, state
real(8), allocatable :: phi(:,:)
! external functions
real(8) sdelta
external sdelta
! initialise universal variables
call init0
call init1
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evalfv(nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(phi(2,nmatmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(apw_match(2,ngkmax,apwordmax,lmmaxapw,nspnfv))
allocate(apw_match_real(2,ngkmax,apwordmax,2*lmaxapw+1))
allocate(gvecs(3,ngkmax))
allocate(atom_types(natmtot))
allocate(atom_pos(3,natmtot))
allocate(species_local(nlotot))
allocate(atom_local(nlotot))
allocate(l_local(nlotot))
allocate(m_local(nlotot))
allocate(ilo_local(nlotot))

! read the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! set the occupancies

do ispn=1,nspnfv
  do ik=1,nkpt
    call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
         sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
  end do
end do


write(*,'("Writing HFD5file for QMC import...")')

! Initialized HDF5
call h5open_f(h5err)
! Create the file
call h5fcreate_f(fname, H5F_ACC_TRUNC_F, h5id, h5err, H5P_DEFAULT_F, &
                 H5P_DEFAULT_F)

! Write the file type
call h5ltmake_dataset_string_f (h5id, 'type', 'LAPW', h5err)

! Create the ions group
call h5gcreate_f (h5id, 'ions', ion_group, h5err)
! Write the geometry


! Close the group
call h5gclose_f (ion_group, h5err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write the APW  and local functions !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create the apw group
call h5gcreate_f (h5id, 'radialfuncs', radial_group, h5err)
! Write out lmaxapw
n(1) = lmaxapw
dim1(1) = 1
call h5ltmake_dataset_int_f (radial_group, 'lmax', 1, dim1, n, h5err)
! Loop over species
do is=1,nspecies
!  write(group_name, '("species_",I1)'),is-1
  group_name = group("species",is-1);
  call h5gcreate_f (radial_group, group_name, species_group, h5err)
  nr=nrmt(is)
  dim1(1) = nr
  n(1) = natoms(is)
  call h5ltmake_dataset_double_f (species_group, 'r', 1, dim1, spr(:,is), h5err)
  dim1(1) = 1
  ! Write the number of atoms in the species
  call h5ltmake_dataset_int_f    (species_group, 'num_atoms', 1, dim1,&
                                  n, h5err)
  ! Write the muffin tin radius
  val1(1) = rmt(is)
  call h5ltmake_dataset_double_f (species_group, 'radius', 1, dim1,&
       val1, h5err)
  do ia=1,natoms(is)
     ! create atom group
!     write(group_name, '("atom_",I1)'),ia-1
     group_name = group("atom", ia-1)
     call h5gcreate_f (species_group, group_name, atom_group, h5err)
     ! Write the atom position
     dim1(1)=3
     call h5ltmake_dataset_double_f (atom_group, 'pos_reduced', 1, dim1, &
          atposl(:,ia,is), h5err)
     call h5ltmake_dataset_double_f (atom_group, 'pos_cart', 1, dim1, &
          atposc(:,ia,is), h5err)
     ias=idxas(ia,is)
     ! Loop over local orbitals, v(r)
     do ilo=1,nlorb(is)
        ! write(*,'("lorbord=",I1)'),lorbord(ilo,is)
!        write(group_name, '("local_",I1)'),ilo-1
        group_name = group ("local",ilo-1)
        call h5gcreate_f (atom_group, group_name, local_group, h5err)
        n(1) = lorbl(ilo,is)
        dim1(1) = 1
        call h5ltmake_dataset_int_f (local_group, 'l', 1, dim1, n, h5err)
        dim1(1) = nr
        call h5ltmake_dataset_double_f (local_group, 'v_of_r', 1, dim1,&
             lofr(:,1,ilo,ias), h5err)
        call h5gclose_f (local_group, h5err)
     end do
     
     ! Loop over the l channels for APW functions, u(r)
     do l=0,lmaxapw
        dim2(1) = nr
        dim2(2)=apword(l,is)
!        write(group_name, '("APW_l_",I1)'),l
        group_name = group("APW_l",l)
        call h5gcreate_f (atom_group, group_name, l_group, h5err) 
        call h5ltmake_dataset_double_f (l_group, 'u_of_r', 2, dim2,&
             apwfr(:,1:1,:,l,ias),h5err)
        dim1(1)=apword(l,is)
        call h5ltmake_dataset_double_f (l_group, 'du_dr_final', 1, dim1,&
             apwdfr(:,l,ias), h5err)
        ! loop over spins
        do ispn=1,nspnfv
!           write(group_name, '("spin_",I1)'),ispn-1
           group_name = group("spin",ispn-1)
           call h5gcreate_f (l_group, group_name, spin_group, h5err)
           ! loop over k-points
           do ik=1,nkpt
!              write(group_name, '("twist_",I1)'),ik-1
              group_name = group("twist",ik-1)
              call h5gcreate_f (spin_group, group_name, twist_group, h5err)
              ! Compute matching condition
              call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
                   sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
              ! Collect matching functions into an array
              do m=-l,l 
                 lm = idxlm(l,m)
                 do ord=1,apword(l,is)
                    do ig=1,ngk(ik,ispn)
                       apw_match(1,ig,ord,m+l+1,ispn) = &
                            real(apwalm(ig,ord,lm,ias,ispn))
                       apw_match(2,ig,ord,m+l+1,ispn) = &
                            dimag(apwalm(ig,ord,lm,ias,ispn))
                    end do ! ig=1,nkg(ik,ispn)
                 end do ! ord=1,apword(l,is)
              end do ! m=-l,l
              dim3(1)=ngk(ik,ispn)
              dim3(2)=apword(l,is)
              dim3(3)=2*l+1
              call h5ltmake_dataset_double_f(twist_group, 'match_real', 3, dim3, &
                   apw_match(1,1:ngk(ik,ispn),1:apword(l,is),1:2*l+1,ispn), h5err)
              call h5ltmake_dataset_double_f(twist_group, 'match_imag', 3, dim3, &
                   apw_match(2,1:ngk(ik,ispn),1:apword(l,is),1:2*l+1,ispn), h5err)
              call h5gclose_f (twist_group, h5err)
           end do ! ik=1,nkpt
           call h5gclose_f (spin_group, h5err)
        end do ! ispn=1,nspnfv
        call h5gclose_f (l_group, h5err)        
     end do ! l=0,lmaxapw
     call h5gclose_f (atom_group, h5err)
  end do ! ia=1,natoms(is)
  call h5gclose_f (species_group, h5err)
end do ! is=1,nspecies
! Close the APW group
call h5gclose_f (radial_group, h5err)

! Create the eignestates group
call h5gcreate_f (h5id, 'eigenstates', eig_group, h5err)

! Loop over species
!do is=1,nspecies
!  nr=nrmt(is)
!  do ir=1,nr
!    write(*,'(G18.10)'),spr(ir,is)
!  end do
!  do ia=1,natoms(is)
!    ias=idxas(ia,is)
!  end do
!end do


! Write the mapping between local orbital index and species,
! atom, l, m
dim1=nlotot
do is=1,nspecies
   do ia=1,natoms(is)
      ias = idxas(ia,is)
      do ilo=1,nlorb(is)
         l=lorbl(ilo,is)
         do m=-l,l
            lm=idxlm(l,m)
            i=idxlo(lm,ilo,ias)
            species_local(i) = is
            atom_local(i) = ias-1
            l_local(i) = l
            m_local(i) = m
            ilo_local(i) = ilo-1
         end do ! m=-l,l
      end do ! ilo=1,nlorb(is)
   end do ! ia=1,natoms(is)
end do !ia=1,nspecies
call h5ltmake_dataset_int_f(eig_group,'local_species', 1, dim1,&
     species_local, h5err)
call h5ltmake_dataset_int_f(eig_group,'local_atom'   , 1, dim1,&
     atom_local, h5err)
call h5ltmake_dataset_int_f(eig_group,'local_l'      , 1, dim1,&
     l_local, h5err)
call h5ltmake_dataset_int_f(eig_group,'local_m'      , 1, dim1,&
     m_local, h5err)
call h5ltmake_dataset_int_f(eig_group,'local_ilo'    , 1, dim1,&
     ilo_local, h5err)

! Loop over k-points
do ispn=1,nspnfv
!   write(group_name,'("spin_",I1)'),ispn-1
   group_name = group("spin",ispn-1)
   call h5gcreate_f (eig_group, group_name, spin_group, h5err)
   do ik=1,nkpt
      ng = ngk(ik,ispn)
      call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
      call getevalfv(vkl(1,ik),evalfv)
      ! call getevecsv(vkl(1,ik),evecsv)
      ! Create twist group
!      write(group_name, '("twist_",I1)'),ik-1
      group_name = group("twist",ik-1)
      ! Write the twist_angle
      call h5gcreate_f (spin_group, group_name, twist_group, h5err)
      dim1(1) = 3
      call h5ltmake_dataset_double_f(twist_group, 'twist_angle', 1, dim1,&
           vkl(:,ik), h5err)
      ! Write out the G-vectors in lattice coordinates
      dim2(1) = 3
      dim2(2) = ng
      do ig=1,ng
         gvecs(:,ig) = vgkl(:,ig,ik,ispn)-vkl(:,ik)
      end do
      call h5ltmake_dataset_double_f(twist_group, 'gvecs_reduced', 2, dim2,&
           gvecs, h5err)
      do ig=1,ngk(ik,ispn)
         gvecs(:,ig) = vgkc(:,ig,ik,ispn)-vkc(:,ik)
      end do
      call h5ltmake_dataset_double_f(twist_group, 'gvecs_cart', 2, dim2,&
           gvecs, h5err)

      ! Loop over bands
      do ist=1,nstfv
!         write(group_num,'(I)'),ist-1
!         write(group_name,'("band_",A)'),trim(adjustl(group_num))
         group_name = group("band",ist-1)
         call h5gcreate_f(twist_group, group_name, band_group, h5err)
         dim1(1) = 1
         val1(1) = evalfv(ist,ispn)
         call h5ltmake_dataset_double_f(band_group, 'eigenvalue', 1, dim1,&
              val1, h5err)
         dim2(1) = 2
         dim2(2) = ngk(ik,ispn)
         do ig=1,ng
            phi(1,ig) = dreal(evecfv(ig,ist,ispn))
            phi(2,ig) = dimag(evecfv(ig,ist,ispn))
         end do
         ! Loop over species
         do is=1,nspecies
            do ia=1,natoms(is)
               ias = idxas(ia,is)
               do ilo=1,nlorb(is)
                  l=lorbl(ilo,is)
                  do m=-l,l
                     lm=idxlm(l,m)
                     i=ng+idxlo(lm,ilo,ias)
                     phi(1,i) = dreal(evecfv(i,ist,ispn))
                     phi(2,i) = dimag(evecfv(i,ist,ispn))
                  end do ! m=-l,l
               end do ! ilo=1,nlorb(is)
            end do ! ia=1,natoms(is)
         end do ! is=1,nspecies
         call h5ltmake_dataset_double_f(band_group, 'phi_apw', 2, dim2,&
              phi(:,1:ng), h5err)
         dim2(2) = nmat(ik,ispn)-ng
         call h5ltmake_dataset_double_f(band_group, 'phi_local', 2, dim2,&
              phi(:,ng+1:nmat(ik,ispn)), h5err)
         
         call h5gclose_f(band_group, h5err)
      end do
      
      ! Close twist group
      call h5gclose_f (twist_group, h5err)
   end do ! ik=1,nkpt
   call h5gclose_f (spin_group, h5err)
end do ! isp=1,nspnf
! Close the group
call h5gclose_f (eig_group, h5err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write the core states !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Create the core_states group
call h5gcreate_f (h5id, 'corestates', core_group, h5err)
! Loop over species
do is=1,nspecies
!   write (group_name, '("species_",I1)'),is-1
   group_name = group("species",is-1)
   call h5gcreate_f (core_group, group_name, species_group, h5err)
   dim1(1) = spnr(is) - 4
   call h5ltmake_dataset_double_f (species_group, 'r', 1, dim1,&
        spr(5:spnr(is),is), h5err)
   do ia=1,natoms(is)
      ias = idxas(ia,is)
!      write (group_name, '("atom_",I1)'),ia-1
      group_name = group("atom", ia-1)
      call h5gcreate_f (species_group, group_name, atom_group, h5err)
      icore=0
      do ist=1,spnst(is)
         if (spcore(ist,is)) then
!            write (group_name, '("state_",I1)'),icore
            group_name = group("state",icore)
            call h5gcreate_f (atom_group, group_name, state_group, h5err)
            dim1(1) = spnr(is) - 4
            call h5ltmake_dataset_double_f(state_group, 'g0', 1, dim1,&
                 rwfcr(5:spnr(is),1,ist,ias), h5err)
            call h5ltmake_dataset_double_f(state_group, 'f0', 1, dim1,&
                 rwfcr(5:spnr(is),2,ist,ias), h5err)
            dim1(1) = 1
            val1(1) = evalcr(ist,ias)
            call h5ltmake_dataset_double_f(state_group, 'eigenvalue',1,dim1,&
                 val1, h5err)
            n(1) = spn(ist,is)
            call h5ltmake_dataset_int_f(state_group, 'n', 1, dim1, n, h5err)
            n(1) = spl(ist,is)
            call h5ltmake_dataset_int_f(state_group, 'l', 1, dim1, n, h5err)
            n(1) = spk(ist,is)
            call h5ltmake_dataset_int_f(state_group, 'k', 1, dim1, n, h5err)
            call h5gclose_f  (state_group, h5err)
            icore = icore + 1
         end if
      end do ! ist=spnst(is)
      call h5gclose_f (atom_group, h5err)
   end do ! ia=1,natoms(is)
   call h5gclose_f (species_group, h5err)
end do ! is=1,nspecies
call h5gclose_f (core_group, h5err)


! Create the parameters group
call h5gcreate_f (h5id, 'parameters', param_group, h5err)
dim2(1)=3
dim2(2)=3
call h5ltmake_dataset_double_f(param_group, 'lattice', 2, dim2, avec, h5err)
call h5ltmake_dataset_double_f(param_group, 'reciprocal_lattice', &
     2, dim2, bvec, h5err)
do is=1,nspecies
   do ia=1,natoms(is)
      atom_types(idxas(ia,is)) = int(-spzn(is))
      atom_pos(:,idxas(ia,is)) = atposc(:,ia,is)
   end do
end do
dim1(1)=natmtot
call h5ltmake_dataset_int_f(param_group, 'atom_types', 1, dim1, &
     atom_types, h5err)
dim2(1)=3
dim2(2)=natmtot
call h5ltmake_dataset_double_f(param_group, 'atom_pos', 2, dim2, &
     atom_pos, h5err)


call h5gclose_f (param_group, h5err)
! Close the file
call h5fclose_f(h5id, h5err)


deallocate(apwalm,apw_match,apw_match_real,gvecs,atom_types,atom_pos,phi)
deallocate(species_local,atom_local,l_local,m_local,ilo_local)
deallocate(evecfv,evalfv,evecsv)
write(*,'("Done.")')

return
end subroutine

