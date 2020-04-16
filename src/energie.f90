!! THIS MODULE CONTAINS THE SUBROUTINES THAT COMPUTE
!! THE UFF ENERGIES, ACCORDING TO THE UFF IMPLEMENTATION
!!
!! Rappé A. K. et al, J. Am. Chem. Soc. 1992, 114 10024-10035
MODULE energie
        IMPLICIT NONE
        CONTAINS
        
        ! INPUT     :
        !               - a, b : list of reals, UFF properties for atom a and b
        !               - n    : real, the bond order between atom a and b
        ! OPERATION :
        !               Compute the req according to UFF
        ! RETURN    :
        !               Real, the equilibrium distance
        REAL FUNCTION req(a,b,n)
                REAL, DIMENSION(11), INTENT(IN) :: a,b
                REAL,INTENT(IN)                 :: n
                REAL                            :: ri,rj,rBO,rEN,xi,xj

                ri = a(1)
                rj = b(1)
                xi=a(9)
                xj=b(9)
                rBO = -0.1332*(ri+rj)*LOG(REAL(n))
                rEN = ri*rj*((SQRT(xi)-SQRT(xj))**2)/(xi*ri+xj*rj)

                req = ri + rj + rBO - rEN
        END FUNCTION

        ! INPUT     :
        !               - a, b : list of reals, UFF properties for atom a and b
        !               - n    : real, the bond order between atom a and b
        ! OPERATION :
        !               Compute the stretching force constant according to UFF
        ! RETURN    :
        !               Real, the force constant
        REAL FUNCTION kbond(a,b,n)
                REAL, DIMENSION(11),INTENT(IN) :: a,b
                REAL,INTENT(IN)                :: n

                kbond = 664.12*(a(6)*b(6))/(req(a,b,n)**3)
        END FUNCTION

        ! INPUT     :
        !               - a, b : list of real, UFF properties for atom a and b
        !               - n    : real, the bond order between atom a and b
        !               - d    : real, the distance between atom a and b
        ! OPERATION :
        !               Compute the UFF energy for a bond strech, with and
        !               harmonic potential
        ! RETURN    :
        !               Real, the energy
        REAL FUNCTION ebond(a,b,n,d)

                REAL, DIMENSION(11), INTENT(IN) :: a,b
                REAL, INTENT(IN)                :: d,n

                ebond = 0.5*kbond(a,b,n)*(d-req(a,b,n))**2
        END FUNCTION

        ! INPUT     :
        !               - positions : matrix of real, contain the atomic positions
        !               - B         : matrix of reals, bond order matrix
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
                REAL, DIMENSION(:,:), INTENT(IN)           :: B
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
        !               - nab,nbc : real, the bond order between atoms ab and bc
        ! OPERATION :
        !               Compute the bending force constant according to UFF
        ! RETURN    :
        !               Real, the force constant
        REAL FUNCTION kangle(a,b,c,nab,nbc)
                REAL, DIMENSION(11),INTENT(IN) :: a,b,c
                REAL,INTENT(IN)                :: nab,nbc
                REAL                           :: beta,prefactor,endfactor,rterm,rab,rbc,rac,theta0

                theta0 = b(2)*3.14159265359/180
                rab = req(a,b,nab)
                rbc = req(b,c,nbc)
                rac = SQRT(rab**2+rbc**2-2*rab*rbc*COS(theta0))
                rterm=rab*rbc


                beta = 664.12/(rab*rbc)
                prefactor = beta*(a(6)*c(6))/(rac**5)
                endfactor = 3*rterm*(1-COS(theta0)**2)-rac*rac*COS(theta0)

                kangle =prefactor*rterm*endfactor
        END FUNCTION

        ! INPUT     :
        !               - a, b, c  : list of real, UFF properties for atom a and b
        !               - nab, nbc : real, the bond order between atom ab and bc
        !               - theta    : real, the angle abc
        ! OPERATION :
        !               Compute the UFF energy for a bend abc
        ! RETURN    :
        !               Real, the energy
        REAL FUNCTION eangle(a,b,c,nab,nbc,theta)
                REAL, DIMENSION(11), INTENT(IN) :: a,b,c
                REAL, INTENT(IN)                :: nab,nbc
                REAL, INTENT(IN)                :: theta
                REAL                            :: theta0, c0, c1, c2

                theta0 = b(2)*3.14159265359/180

                c2 = 1/(4*SIN(theta0)**2)
                c1 = -4*c2*COS(theta0)
                c0 = c2*(2*COS(theta0)**2+1)

                eangle = kangle(a,b,c,nab,nbc)*(c0+c1*COS(theta)+c2*COS(2*theta))
        END FUNCTION

        ! INPUT     :
        !               - positions : matrix of real, contain the atomic positions
        !               - B         : matrix of real, bond order matrix
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
                REAL, DIMENSION(:,:), INTENT(IN)           :: B
                REAL, DIMENSION(11)                        :: a,d,c 
                REAL                                       :: theta
                INTEGER                                    :: i,j,k,l


                DO i=1,SIZE(names)
                        DO j=1,SIZE(names)
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
                angle_energy = angle_energy/2
        END FUNCTION

        ! INPUT :
        !		- a, b : list of reals, UFF properties for atom a and b
        ! OPERATION :
        !		Compute the well depth according to UFF
        ! RETURN :
        !		Real, the well depth in kcal/mol
        REAL FUNCTION wd(a,b)
        	REAL, DIMENSION(11), INTENT(IN) :: a,b
        
                wd = SQRT(a(4)*b(4))
        END FUNCTION
        
        ! INPUT :
        !		- a, b : list of reals, UFF properties for atom a and b 
        ! OPERATION :
        !		Compute the Van Der Waals bond length
        ! RETURN :
        !		Real, the Van Der Waals bond lenght in Angström
        REAL FUNCTION vdwl(a,b)
        	REAL, DIMENSION(11), INTENT(IN) :: a,b
                
                vdwl = SQRT(a(3)*b(3))
        END FUNCTION
        
        ! INPUT : 
        !		- a, b : list of reals, UFF properties for atom a and b
        !		- d    : real, distance between a and b
        ! OPERATION :
        !		Compute the Van Der Waals energy between two atoms
        ! RETURN :
        !		Real, the Van Der Waals energy 
        REAL FUNCTION evdw(a,b,d)
        
        	REAL, DIMENSION(11), INTENT(IN) :: a,b
        	REAL, INTENT(IN)                :: d
                REAL :: Dab,xab,r6,r12

                Dab = wd(a,b)
                xab = vdwl(a,b)
                r6 = (xab/d)**6
                r12 = r6*r6 


                evdw = Dab*(r12-2*r6)
        END FUNCTION
        
        ! INPUT :
        !		- positions	: matrix of real, contains the atomic positions
        !		- names		: list of charater, atomic UFF type
        !		- types		: list of character, all UFF types in the molecule
        !		- prop		: matrix of real, the UFF properties for the atom types
        ! OPERATION :
        !		Loop over all atom couples to compute the Van Der Waals energy 
        ! RETURN :
        !		Real, sum of all atoms couples to compute Van Der Waals energy
        REAL FUNCTION vdw_energy(positions,names,types,prop,adj)
                USE math
        
                CHARACTER(len=5), DIMENSION(:), INTENT(IN) :: names, types
        	REAL, DIMENSION(:,:), INTENT(IN)           :: positions, prop
                INTEGER, DIMENSION(:,:), INTENT(IN)        :: adj
        	REAL, DIMENSION(11)                        :: a, b
        	REAL                                       :: d
                INTEGER                                    :: i, j, k, l

                vdw_energy = 0
                DO i=1,SIZE(names)-1
                        DO j=i+1,SIZE(names)
                                d = distance(positions(i,:),positions(j,:))
                                IF (d.LT.20 .AND. adj(i,j)==0) THEN
                                        DO l=1,SIZE(names)
                                                IF(adj(i,l)==1 .AND. adj(j,l)==1) THEN
                                                        EXIT
                                                ELSE IF(l==SIZE(names)) THEN
                                                        DO k=1, SIZE(types)
                                                                IF(types(k)==names(i)) THEN
                                                                        a=prop(k,:)
                                                                        EXIT
                                                                END IF
                                                        END DO
                                                        DO k=1, SIZE(types)
                                                                IF(types(k)==names(j)) THEN
                                                                        b=prop(k,:)
                                                                        EXIT
                                                                END IF
                                                        END DO
                                                        vdw_energy = vdw_energy + evdw(a,b,d)
                                                END IF
                                        END DO
                                END IF
                        END DO
               END DO
        END FUNCTION

! INPUT		:
!			- j     : the name (len=5) of the atom j
!			- k     : the name (len=5) of the atom k
!			- bjk   : the j-k bond order
!                       - n     : the periodicity
!                       - phi0  : the equilibrium dihedral angle
!                       - v     : the barrier
!                       - nj    : the number of neighbors to j
!                       - nk    : the number of neighbors to k
!                       - i     : the name (len=5) of the atom i
!                       - l     : the name (len=5) of the atom l
!                       - vj,vk : the V parameter for the atom
!                       - uj,uk : the U parameter for the atom
! OPERATION	:
!			Compute the barrier to rotation
! RETURN	:
!			- v, real, the value of barrier to rotation in kcal/mol
!			- n, integer, the periodicity
!			- phi0, real, the equilibrium angle
SUBROUTINE find_barrier(j,k,bjk,n,phi0,v,nj,nk,i,l,vj,vk,uj,uk)
        CHARACTER(len=5), INTENT(IN)                    :: j,k,i,l
        INTEGER, INTENT(IN)                             :: nj,nk
        INTEGER, INTENT(INOUT)                          :: n
        REAL, INTENT(IN)                                :: vj,vk,uj,uk,bjk
        REAL                                            :: nvj,nvk
        REAL, INTENT(INOUT)                             :: phi0,v
        REAL                                            :: pi
        CHARACTER(len=2)                                :: tj,tk
        CHARACTER(len=1)                                :: hj,hk,hi,hl !Hybridation   3 = sp3   2 = sp2   R = resonant

        pi = 3.14159265359
        READ(i,'(2a,1a)') tk,hi
        READ(l,'(2a,1a)') tk,hl
        READ(j,'(2a,1a)') tj,hj
        READ(k,'(2a,1a)') tk,hk

        nvj = vj
        nvk = vk

        !Check for sp3-sp3 generic case
        IF(hj=='3' .AND. hk=='3') THEN
                n = 3
                phi0 = pi/3
                IF((tj=='O_' .OR. tj=='Te' .OR. tj=='Se' .OR. tj=='Po' .OR. tj=='S_') .AND.&
                 & (tk=='O_' .OR. tk=='S_' .OR. tk=='Se' .OR. tk=='Te'.OR. tk=='Po')) THEN
                        nvj = 6.8
                        n = 2
                        phi0 = pi/2
                        IF(tj=='O_') THEN
                                nvj = 2
                        END IF
                        IF(tk=='O_') THEN
                                nvk = 2
                        END IF
                END IF
                v = SQRT(vj*vk)
        !Check for sp2-sp2 generic case
        ELSE IF((hj=='2' .OR. hj=='R').AND. (hk=='2'.OR.hk=='R')) THEN
                v = 5*SQRT(uj*uk)*(1+4.18*LOG(bjk))
                n = 2
                phi0 = pi
        !Check for the sp2-sp3 case
        ELSE IF(((hj=='2'.OR.hj=='R') .AND. hk=='3') .OR. (hj=='3' .AND. (hk=='2'.OR.hj=='R'))) THEN
                n = 6
                v = 1
                phi0 = 0
                IF((tj=='O_' .OR. tj=='Te' .OR. tj=='Se' .OR. tj=='Po' .OR. tj=='S_').AND.hj=='3') THEN
                        IF(tk/='O_' .OR. tk/='Te' .OR. tk/='Se' .OR. tk/='Po' .OR. tk/='S_') THEN
                                n = 2
                                phi0 = pi/2
                                V = 5*SQRT(uj*uk)*(1+4.18*LOG(bjk))
                        END IF
                END IF
                IF((tk=='O_' .OR. tk=='Te' .OR. tk=='Se' .OR. tk=='Po' .OR. tk=='S_').AND.hk=='3') THEN
                        IF(tj/='O_' .OR. tj/='Te' .OR. tj/='Se' .OR. tj/='Po' .OR. tj/='S_') THEN
                                n = 2
                                phi0 = pi/2
                                V = 5*SQRT(uj*uk)*(1+4.18*LOG(bjk))
                        END IF
                END IF
                IF(hj=='2') THEN
                        IF(hi=="2") THEN
                                v = 2
                                n = 3
                                phi0 = pi
                        END IF
                END IF
                IF(hk=='2') THEN
                        IF(hl=="2") THEN
                                v = 2
                                n = 3
                                phi0 = pi
                        END IF
                END IF
        ELSE
                v = 1
                n = 6
                phi0 = 0
        END IF
END SUBROUTINE

! INPUT		:
!			- phi   : the dihedral angle
!			- j,k   : the name(len=5) of the atom
!			- bjk   : the j-k bond order
!			- nj,nk : the number of neighbors to j,k
!                       - vj,vk : the V parameter for the atom
!                       - uj,uk : the U parameter for the atom
! OPERATION 	: 
!			Compute the torsion energy
! RETURN	:
!			Real, the torsion energy

REAL FUNCTION t_energy(phi,j,k,bjk,nj,nk,i,l,vj,vk,uj,uk)
        REAL, INTENT(IN)             :: vj,vk,uj,uk
        INTEGER,INTENT(IN)           :: nj,nk
        REAL, INTENT(IN)             :: phi,bjk
	REAL                         :: v, phi0
        INTEGER                      :: n
        CHARACTER(len=5), INTENT(IN) :: j,k,i,l
        CALL find_barrier(j,k,bjk,n,phi0,v,nj,nk,i,l,vj,vk,uj,uk) 

        t_energy = v*(1-COS(n*phi)*COS(n*phi0))/2
        t_energy = t_energy/((nj-1)*(nk-1))
END FUNCTION

! INPUT		:
!                       - names		: list of character, the names of the atoms with the format 'xxV  ' xx represent the atom and V the numver of neighbours
!                       - B		: matrix of real, bond order matrix
!			- positions 	: matrix of real, contain the atomic positions
!                       - D		: matrix of integer, the adjency matrix
!                       - prop          : the list of the atomic properties
!                       - types         : the list of the atom types
! OPERATION	:
!			Loop over all atom quadrouple to compute the torsion energy
! RETURN	:
!			Real,sum of all the torsion energy

REAL FUNCTION torsion_energy(names,B,positions, D,prop,types)
        USE math

	REAL, DIMENSION(:,:), INTENT(IN)                :: positions,prop,B
        INTEGER, DIMENSION(:,:), INTENT(IN)             :: D
        CHARACTER(len=5), DIMENSION(:), INTENT(IN)      :: names,types
	REAL                                            :: v, phi0, phi
        INTEGER                                         :: n,i,j,k,l,p
        REAL, DIMENSION(11)                             :: propj,propk

        torsion_energy = 0

        DO i=1,SIZE(names)
                DO j=1,SIZE(names)
                        IF (B(i,j)/=0) THEN
                                DO k=1,SIZE(names)
                                        IF(B(k,j)/=0 .AND. k/=i) THEN
                                                DO l=1, SIZE(names)
                                                        IF (B(k,l)/=0 .AND. l/=j .AND. l/=i) THEN
                                                                DO p=1, SIZE(types)
                                                                        IF(types(p)==names(j)) THEN
                                                                                propj=prop(p,:)
                                                                                EXIT
                                                                        END IF
                                                                END DO
                                                                DO p=1, SIZE(types)
                                                                        IF(types(p)==names(k)) THEN
                                                                                propk=prop(p,:)
                                                                                EXIT
                                                                        END IF
                                                                END DO
                                                                phi = dihedral(positions(i,:),positions(j,:)&
                                                                        &,positions(k,:),positions(l,:))
                                                                torsion_energy = torsion_energy + &
                                                                &t_energy(phi,names(j),names(k),B(j,k)&
                                                                &,SUM(D(:,j)),SUM(D(:,k)),names(i),names(l),&
                                                                &propj(7),propj(7),propk(8),propk(8))
                                                        END IF
                                                END DO
                                        END IF
                                END DO
                        END IF
                END DO
        END DO
        torsion_energy = torsion_energy/2
END FUNCTION
        !INPUT :
        !		-Kf : real, force constant
        !		-C0,C1,C2 : integers, parameters
        !		-omega : real, variable, improper angle.
        !OPERATION :
        !		computes the improper torsion energy
        !RESULTS :
        !		real, the energy
	REAL FUNCTION e_improper(i,j,k,l,omega)

                CHARACTER(LEN=5)                ::i,j,k,l
		REAL, INTENT(IN)                ::omega
                REAL                            ::C0,C1,C2
		REAL                            ::Kf,eimproper
                CALL improper_parameters(i,j,k,l,c0,c1,c2,Kf)

                eimproper=Kf*(C0+C1*COS(omega)+C2*COS(2*omega))
        END FUNCTION
        
        !INPUT 	   :
        !	   -i,j,k,l : name (len=5) of the atoms 1,j,k,l, i being the central atom bonded to the 3 others
        !OPERATION :
        !	   gives the values of the parameters Kf,C0,C1 and C2 in function of the nature of the atoms
        !RESUTS    :
        !	   the parameters		
        SUBROUTINE improper_parameters(i,j,k,l,c0,c1,c2,Kf)
                
                CHARACTER(len=5),INTENT(IN)              :: i,j,k,l
                CHARACTER(len=2)                         :: ti,tj,tk,tl
                CHARACTER(len=1)                         :: hi
                REAL                                     :: Kf
		REAL                                     :: pi,C0,C1,C2
                pi=3.14159265359

                READ(i,'(2a,1a)') ti,hi               !ti: type for the atom i 
                READ(j,'(2a)') tj
                READ(k,'(2a)') tk
                READ(l,'(2a)') tl


                IF(ti=='C_' .AND. (hi=='R' .OR. hi=='2')) THEN
                        IF(tj=='O_' .OR. tk=='O_' .OR. tl=='O_' .OR. tj=='N_' .OR. tk=='N_' .OR. tl=='N_') THEN
                                Kf=50.0
                                C0=1
                                C1=-1
                                C2=0
                        ELSE
                                Kf=6.0
                                C0=1
                                C1=-1
                                C2=0
                        END IF

                ELSE IF (ti=='P_') THEN
                        C2=1
                        C1=-4*COS(84.4339*pi/180)
                        C0=-(C1*COS(84.4339*pi/180)+C2*COS(2*84.4339*pi/180))
                        Kf=22.0/(C0+C1+C2)

                ELSE IF (ti=='As') THEN
                        C2=1
                        C1=-4*COS(86.9735*pi/180)
                        C0=-(C1*COS(86.9735*pi/180)+C2*COS(2*86.9735*pi/180))
                        Kf=22.0/(C0+C1+C2)

                ELSE IF (ti=='Sb') THEN
                        C2=1
                        C1=-4*COS(87.7047*pi/180)
                        C0=-(C1*COS(87.7047*pi/180)+C2*COS(2*87.7047*pi/180))
                        Kf=22/(C0+C1+C2)

                ELSE IF (ti=='Bi') THEN
                        C2=1
                        C1=-4*COS(pi/2)
                        C0=-(C1*COS(pi/2)+C2*COS(pi))
                        Kf=22.0/(C0+C1+C2)
                ELSE
                        Kf=0.0
                        C0=0
                        C1=0
                        C2=0
                END IF

        END SUBROUTINE

        ! INPUT     :
        !               - positions : matrix of real, contain the atomic positions
        !               - B         : matrix of integer, bond order matrix
        !               - names     : list of character, atomic UFF type
        ! OPERATION :
        !               Loop over all atom quadruplets to compute improper torsion energy energy
        ! RETURN    :
        !               Real, sum of all the angle bend energy   
        REAL FUNCTION improper_energy(positions,B,names)
                USE math

                CHARACTER(len=5), DIMENSION(:), INTENT(IN) :: names
                REAL, DIMENSION(:,:), INTENT(IN)           :: B,positions
		REAL                                       :: omega
                INTEGER                                    :: i,j,k,l
                
                improper_energy=0
                DO i=1,SIZE(names)
                        DO j=1,SIZE(names)
                                IF(B(i,j)/=0) THEN
                                        DO k=1,SIZE(names)
                                                IF(B(i,k)/=0 .AND. j/=k) THEN
                                                        DO l=1,SIZE(names)
                                                                IF(B(i,l)/=0 .AND. j/=l .AND. l/=k) THEN
                                                                        omega=improper(positions(i,:)&
                                                                        &,positions(j,:),positions(k,:)&
                                                                        &,positions(l,:))
                                                                        improper_energy=improper_energy &
                                                                        & +e_improper(names(i),names(j),&
                                                                        &names(k),names(l),omega)
                                                                END IF
                                                        END DO
                                                END IF
                                        END DO
                                END IF
                        END DO
                END DO
        improper_energy=improper_energy/6
        END FUNCTION


SUBROUTINE compute_energy(positions,B,names,types,prop,D,outfile)
        USE ecriture
        CHARACTER(len=5), DIMENSION(:), INTENT(IN)    :: names,types
        REAL, DIMENSION(:,:), INTENT(IN)              :: B,positions,prop
        INTEGER, DIMENSION(:,:), INTENT(IN)           :: D
        INTEGER, INTENT(IN)                           :: outfile
        REAL                                          :: benergy,aenergy,tenergy,ienergy,vdwenergy


        benergy = bond_energy(positions,B,names,types,prop)
        aenergy = angle_energy(positions,B,names,types,prop)
        tenergy = torsion_energy(names,B,positions,D,prop,types)
        ienergy = improper_energy(positions,B,names)
        vdwenergy = vdw_energy(positions,names,types,prop,D)

        CALL write_energies(outfile,benergy,aenergy,tenergy,ienergy,vdwenergy)
END SUBROUTINE
END MODULE

