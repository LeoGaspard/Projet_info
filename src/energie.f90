!! THIS MODULE CONTAINS THE SUBROUTINES THAT COMPUTE
!! THE UFF ENERGIES, ACCORDING TO THE UFF IMPLEMENTATION
!!
!! Rappé A. K. et al, J. Am. Chem. Soc. 1992, 114 10024-10035
MODULE energie
        IMPLICIT NONE
        CONTAINS
        
        ! INPUT     :
        !               - a, b : list of reals, UFF properties for atom a and b
        !               - n    : integer, the bond order between atom a and b
        ! OPERATION :
        !               Compute the req according to UFF
        ! RETURN    :
        !               Real, the equilibrium distance
        REAL FUNCTION req(a,b,n)
                REAL, DIMENSION(11), INTENT(IN) :: a,b
                INTEGER, INTENT(IN)             :: n

                req = a(1)+b(1)-0.1332*(a(1)+b(1))*LOG(REAL(n))+a(1)*b(1)*(SQRT(a(9))-SQRT(b(9)))**2/(a(9)*a(1)+b(9)*b(1))
        END FUNCTION

        ! INPUT     :
        !               - a, b : list of reals, UFF properties for atom a and b
        !               - n    : integer, the bond order between atom a and b
        ! OPERATION :
        !               Compute the stretching force constant according to UFF
        ! RETURN    :
        !               Real, the force constant
        REAL FUNCTION kbond(a,b,n)
                REAL, DIMENSION(11),INTENT(IN) :: a,b
                INTEGER,INTENT(IN)             :: n

                kbond = 664.12*(a(6)*b(6))/(req(a,b,n)**3)
        END FUNCTION

        ! INPUT     :
        !               - a, b : list of real, UFF properties for atom a and b
        !               - n    : integer, the bond order between atom a and b
        !               - d    : real, the distance between atom a and b
        ! OPERATION :
        !               Compute the UFF energy for a bond strech, with and
        !               harmonic potential
        ! RETURN    :
        !               Real, the energy
        REAL FUNCTION ebond(a,b,n,d)

                REAL, DIMENSION(11), INTENT(IN) :: a,b
                REAL, INTENT(IN)                :: d
                INTEGER, INTENT(IN)             :: n

                ebond = 0.5*kbond(a,b,n)*(d-req(a,b,n))**2
        END FUNCTION

        ! INPUT     :
        !               - positions : matrix of real, contain the atomic positions
        !               - B         : matrix of integer, bond order matrix
        !               - names     : list of character, atomic UFF type
        !               - types     : list of character, all UFF types in the molecule
        !               - prop      : matrix of real, the UFF properties for the atom types
        ! OPERATION :
        !               Loop over all atom couples to compute bond energy if bonded
        ! RETURN    :
        !               Real, sum of all the bond energy within the molecule     
        REAL FUNCTION bond_energy(positions,B,names,types,prop)
                USE math

                CHARACTER(len=5), DIMENSION(:), INTENT(IN) :: names,types
                INTEGER, DIMENSION(:,:), INTENT(IN)        :: B
                REAL, DIMENSION(:,:), INTENT(IN)           :: positions, prop
                REAL, DIMENSION(11)                        :: a,c 
                REAL                                       :: d
                INTEGER                                    :: i,j,k

                DO i=1,SIZE(names)
                        DO j=i,SIZE(names)
                                IF(B(i,j)/=0) THEN
                                        d = distance(positions(i,:),positions(j,:))
                                        DO k=1,SIZE(types)
                                                IF(types(k)==names(i)) THEN
                                                        a=prop(k,:)
                                                        EXIT
                                                END IF
                                        END DO                                        
                                        DO k=1,SIZE(types)
                                                IF(types(k)==names(j)) THEN
                                                        c=prop(k,:)
                                                        EXIT
                                                END IF
                                        END DO                                        
                                        bond_energy = bond_energy + ebond(a,c,B(i,j),d)
                                END IF
                        END DO
                END DO
        END FUNCTION

        ! INPUT     :
        !               - a, b, c : list of reals, UFF properties for atom a and b
        !               - nab,nbc : integer, the bond order between atoms ab and bc
        ! OPERATION :
        !               Compute the bending force constant according to UFF
        ! RETURN    :
        !               Real, the force constant
        REAL FUNCTION kangle(a,b,c,nab,nbc)
                REAL, DIMENSION(11), INTENT(IN) :: a,b,c
                INTEGER, INTENT(IN)             :: nab,nbc

                kangle = (664.12/(req(a,b,nab)*req(b,c,nbc)))*((a(6)*b(6))/(req(a,c,1)**5))&
                        &*req(a,b,nab)*req(b,c,nbc)*(req(a,b,nab)*req(b,c,nbc)*(1-COS(b(2)*3.14159265359/180)**2)-req(a,c,1)**2*&
                        &COS(b(2)*3.14159265359/180))
        END FUNCTION

        ! INPUT     :
        !               - a, b, c  : list of real, UFF properties for atom a and b
        !               - nab, nbc : integer, the bond order between atom ab and bc
        !               - theta    : real, the angle abc
        ! OPERATION :
        !               Compute the UFF energy for a bend abc
        ! RETURN    :
        !               Real, the energy
        REAL FUNCTION eangle(a,b,c,nab,nbc,theta)
                REAL, DIMENSION(11), INTENT(IN)  :: a,b,c
                INTEGER, INTENT(IN)              :: nab,nbc
                REAL, INTENT(IN)                 :: theta
                REAL                             :: c0, c1, c2

                c2 = 1/(4*SIN(b(2)*3.14159265359/180)**2)
                c1 = -4*c2*COS(b(2)*3.14159265359/180)
                c0 = c2*(2*COS(b(2)*3.14159265359/180)**2+1)

                eangle = kangle(a,b,c,nab,nbc)*(c0+c1*COS(theta)+c2*COS(theta)**2)
        END FUNCTION

        ! INPUT     :
        !               - positions : matrix of real, contain the atomic positions
        !               - B         : matrix of integer, bond order matrix
        !               - names     : list of character, atomic UFF type
        !               - types     : list of character, all UFF types in the molecule
        !               - prop      : matrix of real, the UFF properties for the atom types
        ! OPERATION :
        !               Loop over all atom trouples to compute angle bend energy
        ! RETURN    :
        !               Real, sum of all the angle bend energy   
        REAL FUNCTION angle_energy(positions,B,names,types,prop)
                USE math

                CHARACTER(len=5), DIMENSION(:), INTENT(IN) :: names,types
                REAL, DIMENSION(:,:), INTENT(IN)           :: positions, prop
                INTEGER, DIMENSION(:,:), INTENT(IN)        :: B
                REAL, DIMENSION(11)                        :: a,d,c 
                REAL                                       :: theta
                INTEGER                                    :: i,j,k,l


                DO i=1,SIZE(names)
                        DO j=i,SIZE(names)
                                IF(B(i,j)/=0) THEN
                                        DO k=1,SIZE(names)
                                                IF(k/=i .AND. B(j,k)/=0) THEN
                                                        DO l=1,SIZE(types)
                                                                IF(types(l)==names(i)) THEN
                                                                        a=prop(l,:)
                                                                        EXIT
                                                                END IF
                                                        END DO                                        
                                                        DO l=1,SIZE(types)
                                                                IF(types(l)==names(j)) THEN
                                                                        d=prop(l,:)
                                                                        EXIT
                                                                END IF
                                                        END DO                                        
                                                        DO l=1,SIZE(types)
                                                                IF(types(l)==names(k)) THEN
                                                                        c=prop(l,:)
                                                                        EXIT
                                                                END IF
                                                        END DO                                        
                                                        angle_energy = angle_energy + eangle(a,d,c,B(i,j),B(j,k)&
                                                                &,angle(positions(i,:),positions(j,:),positions(k,:)))
                                                END IF
                                        END DO                                        
                                END IF
                        END DO
                END DO

        END FUNCTION
END MODULE
