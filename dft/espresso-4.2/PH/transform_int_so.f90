!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int1_so(int1,na,iflag)
!----------------------------------------------------------------------------
!
! This routine multiply int1 by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in int1_nc.  
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE phus,                 ONLY : int1_nc
!
IMPLICIT NONE

INTEGER :: na, iflag
COMPLEX(DP) :: int1(nhm,nhm,3,nat,nspin_mag)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ipol, np, is1, is2, ijs
COMPLEX(DP) :: fact(4)
LOGICAL :: same_lj

np=ityp(na)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO ipol=1,3
                     ijs=0
                     DO is1=1,npol
                        DO is2=1,npol
                           ijs=ijs+1
                           IF (iflag==0) THEN
                              fact(1)=int1(kh,lh,ipol,na,1)
                           ELSE
                              fact(1)=CONJG(int1(kh,lh,ipol,na,1))
                           ENDIF
                           int1_nc(ih,jh,ipol,na,ijs)=                       &
                               int1_nc(ih,jh,ipol,na,ijs) +                  &
                               fact(1)*                       &
                             (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)  + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)   ) 
                           IF (domag) THEN
                              IF (iflag==0) THEN
                                 fact(2)=int1 (kh,lh,ipol,na,2)
                                 fact(3)=int1 (kh,lh,ipol,na,3)
                                 fact(4)=int1 (kh,lh,ipol,na,4)
                              ELSE
                                 fact(2)=CONJG(int1 (kh,lh,ipol,na,2))
                                 fact(3)=CONJG(int1 (kh,lh,ipol,na,3))
                                 fact(4)=CONJG(int1 (kh,lh,ipol,na,4))
                              ENDIF
                              int1_nc(ih,jh,ipol,na,ijs)=                     &
                                 int1_nc(ih,jh,ipol,na,ijs) +                 &
                                 fact(2)*                       &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)+ &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 (0.D0,-1.D0) * fact(3)*        &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 fact(4)*                      &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE transform_int1_so
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int2_so(int2,nb,iflag)
!----------------------------------------------------------------------------
!
! This routine rotates int2 as appropriate for the spin-orbit case
! and saves it in int2_so.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm 
USE noncollin_module,     ONLY : npol
USE spin_orb,             ONLY : fcoef
USE phus,                 ONLY : int2_so
!
IMPLICIT NONE
INTEGER :: nb, iflag
COMPLEX(DP) :: int2(nhm,nhm,3,nat,nat)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijs, np, is1, is2, na, ipol
COMPLEX(DP) :: fact
LOGICAL :: same_lj

np=ityp(nb)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO na=1,nat 
                     DO ipol=1,3
                        IF (iflag==0) THEN
                           fact=int2(kh,lh,ipol,na,nb) 
                        ELSE
                           fact=CONJG(int2(kh,lh,ipol,na,nb))
                        ENDIF 
                        ijs=0
                        DO is1=1,npol
                           DO is2=1,npol
                              ijs=ijs+1
                              int2_so(ih,jh,ipol,na,nb,ijs)= &
                              int2_so(ih,jh,ipol,na,nb,ijs)+ &
                              fact* &
                            (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np) + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)  )
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE transform_int2_so
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int3_so(int3,na,npert)
!----------------------------------------------------------------------------
!
! This routine multiply int3 by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in int3_nc.  
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE phus,                 ONLY : int3_nc
!
IMPLICIT NONE

COMPLEX(DP) :: int3(nhm,nhm,npert,nat,nspin_mag)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ipol, np, npert, is1, is2, ijs
LOGICAL :: same_lj

np=ityp(na)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO ipol=1,npert
                     ijs=0
                     DO is1=1,npol
                        DO is2=1,npol
                           ijs=ijs+1
                           int3_nc(ih,jh,ipol,na,ijs)=                       &
                               int3_nc(ih,jh,ipol,na,ijs) +                  &
                               int3 (kh,lh,ipol,na,1)*                       &
                             (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)  + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)   ) 
                           IF (domag) THEN
                              int3_nc(ih,jh,ipol,na,ijs)=                     &
                                 int3_nc(ih,jh,ipol,na,ijs) +                 &
                                 int3(kh,lh,ipol,na,2)*                       &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)+ &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 (0.D0,-1.D0) * int3(kh,lh,ipol,na,3)*        &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 int3 (kh,lh,ipol,na,4)*                      &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE transform_int3_so
!

!----------------------------------------------------------------------------
SUBROUTINE transform_int4_so(int4,na)
!----------------------------------------------------------------------------
!
! This routine multiply int4 by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in int4_nc.  
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE phus,                 ONLY : int4_nc
!
IMPLICIT NONE

INTEGER :: na
COMPLEX(DP) :: int4(nhm*(nhm+1)/2,3,3,nat,nspin_mag)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ipol, jpol, np, is1, is2, ijs
INTEGER, ALLOCATABLE :: ijh_save(:,:)
INTEGER :: find_ijh, ijh_l
LOGICAL :: same_lj

ALLOCATE(ijh_save(nhm,nhm))
np=ityp(na)
DO ih=1,nh(np)
   DO jh=1,nh(np)
      ijh_save(ih,jh)=find_ijh(ih,jh,nh(np))
   END DO
END DO

DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  ijh_l=ijh_save(kh,lh)
                  DO ipol=1,3
                     DO jpol=1,3
                        ijs=0
                        DO is1=1,npol
                           DO is2=1,npol
                              ijs=ijs+1
                              int4_nc(ih,jh,ipol,jpol,na,ijs)=               &
                                   int4_nc(ih,jh,ipol,jpol,na,ijs) +         &
                                   int4(ijh_l,ipol,jpol,na,1) *              &
                                  (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)+&
                                   fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                              IF (domag) THEN
                                 int4_nc(ih,jh,ipol,jpol,na,ijs)=             &
                                    int4_nc(ih,jh,ipol,jpol,na,ijs) +         &
                                    int4(ijh_l,ipol,jpol,na,2)*               &
                                  (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)+&
                                  fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 (0.D0,-1.D0) * int4(ijh_l,ipol,jpol,na,3) *   &
                                  (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)-&
                                  fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                   int4(ijh_l,ipol,jpol,na,4)*                &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
DEALLOCATE(ijh_save)
       !
RETURN
END SUBROUTINE transform_int4_so

!----------------------------------------------------------------------------
SUBROUTINE transform_int5_so(int5,nb)
!----------------------------------------------------------------------------
!
! This routine rotates int5 as appropriate for the spin-orbit case
! and saves it in int5_so.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm 
USE noncollin_module,     ONLY : npol
USE spin_orb,             ONLY : fcoef
USE phus,                 ONLY : int5_so
!
IMPLICIT NONE
INTEGER :: nb
COMPLEX(DP) :: int5(nhm*(nhm+1)/2,3,3,nat,nat)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijs, np, is1, is2, na, ipol, jpol

INTEGER, ALLOCATABLE :: ijh_save(:,:)
INTEGER :: find_ijh, ijh_l
LOGICAL :: same_lj

ALLOCATE(ijh_save(nhm,nhm))
np=ityp(nb)
DO ih=1,nh(np)
   DO jh=1,nh(np)
      ijh_save(ih,jh)=find_ijh(ih,jh,nh(np))
   END DO
END DO

DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  ijh_l=ijh_save(kh,lh)
                  DO na=1,nat 
                     DO ipol=1,3
                        DO jpol=1,3
                           ijs=0
                           DO is1=1,npol
                              DO is2=1,npol
                                 ijs=ijs+1
                                 int5_so(ih,jh,ipol,jpol,na,nb,ijs)= &
                                 int5_so(ih,jh,ipol,jpol,na,nb,ijs)+ &
                                 int5(ijh_l,ipol,jpol,na,nb)* &
                               (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np) + &
                                fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)  )
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
DEALLOCATE(ijh_save)
       !
RETURN
END SUBROUTINE transform_int5_so