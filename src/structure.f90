! This module contains all the subroutine that define the structure 
! of the molecule (bond and bond order) from the atomic positions


MODULE structure
        IMPLICIT NONE
        CONTAINS

        ! INPUT     :
        !               - D         : allocatable matrix of integer, void at the begining
        !               - positions : matrix of real number, the atomic coordinates
        !               - names     : matrix of character, the atomic names
        !               - nAtom     : integer, the number of atoms
        ! OPERATION :
        !               For each atom in the list, looks at the distance to every other atom,
        !               if it's inferior to a certain value (sum of the covalent radii + 0.2)
        !               the two atoms are considered bonded
        ! OUTPUT    :
        !               - D        : adjacency matrix
        !                       D(i,j) = 0 if atom i and j are not bonded
        !                       D(i,j) = 0 if atom i and j are bonded
        !                       i and j from 1 to nAtom 
        SUBROUTINE connectivity(D,positions,names,nAtom)
                USE math
                USE lecture
                INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: D
                REAL, DIMENSION(:,:), INTENT(IN)                    :: positions
                CHARACTER(len=5), DIMENSION(:), INTENT(IN)          :: names
                INTEGER, INTENT(IN)                                 :: nAtom
                REAL                                                :: ra,rb
                INTEGER                                             :: i,j,codat


                ! Initializing the adjacency matrix with zeros
                ALLOCATE(D(nAtom,nAtom))
                D = 0
                
                ! Opening the file with the covalent radius informations
                codat=13
                OPEN(UNIT=codat,FILE="data/covalent.dat")

                ! Loop over all atoms to determine connectivity
                DO i=1,nAtom-1
                        DO j=i+1,nAtom
                                ra = get_covalent_radius(names(i),codat)
                                rb = get_covalent_radius(names(j),codat)
                                IF(distance(positions(i,:),positions(j,:)) < 0.8) THEN
                                        PRINT *,"Two atoms on top of each other"&
                                                &,distance(positions(i,:),positions(j,:)),names(i),names(j)
                                ELSE IF(distance(positions(i,:),positions(j,:)) < ra+rb+0.2) THEN
                                        D(i,j) = 1
                                        D(j,i) = 1
                                END IF 
                        END DO
                END DO
                CLOSE(codat)
        END SUBROUTINE

        ! INPUT     :
        !               - D     : matrix of integer, the adjacency matrix
        !               - names : list of character, the names of the atom
        !               - nAtom : integer, the number of atoms
        ! OPERATION :
        !               For every atom, write the number of neighbours in its
        !               name. Will be used to identify the UFF type
        ! OUTPUT    :
        !               - names : list of character, the names of the atom
        !                         and the number of neighbours
        !               
        !                 names(i) = 'C 3  ' for a carbon with 3 neighbours
        SUBROUTINE atom_neighbour(D,names,nAtom)
                INTEGER, DIMENSION(:,:), INTENT(IN)          :: D
                CHARACTER(len=5), DIMENSION(:),INTENT(INOUT) :: names
                INTEGER, INTENT(IN)                          :: nAtom
                INTEGER                                      :: i,j,nv

                DO i=1,nAtom
                        nv = SUM(D(:,i))
                        WRITE(names(i),'(a2,I1)') names(i),nv
                END DO
        END SUBROUTINE


        ! INPUT     :
        !               - D          : matrix of integer, the adjacency matrix
        !               - nAtom      : integer, the number of atoms
        !               - discovered : a list containing the index of the already 
        !                              seen atoms
        !               - start      : integer, the index of the starting atom for
        !                              the recursive DFS part
        !               - parent     : integer, the index of the starting atom for
        !                              the whole DFS
        !               - step       : integer, the number of the DFS current step
        !               - ncycle     : integer, the size of the cycle we are 
        !                              looking for
        !               - nelec      : integer, the number of electrons in the path
        !               - narom      : integer, the number aromatic cycle containing
        !                              parent
        !               - names      : list of character, the names of the atoms
        ! OPERATION :
        !               Use Depth-First Search into the molecular graph (D) to determine
        !               wether the "parent" atom is part of a number "ncycle" of cycles,
        !               a number "narom" of them being aromatic.
        ! OUTPUT    :
        !               - discovered : variable size array of integer, containing the 
        !                              indexes of previously seen atoms
        !                       discovered = [4,2,12] if the path was 4 -> 2 -> 12
        !               - ncycle     : integer, the number of cycle in which the
        !                              "parent" takes part
        !                       ncycle = 2 for a central atom in naphtalene
        !               - nelec      : the number of pi electrons in the current path
        !                       nelec = 5 for the path C-N-C-C
        !               - narom      : the number of aromatic cycles in which the
        !                              "parent" takes part
        RECURSIVE SUBROUTINE find_cycle(D,nAtom,discovered,start,parent,step,ncycle,csize,nelec,narom,names)
                USE ArrayManip

                INTEGER, DIMENSION(:,:),INTENT(IN)                :: D
                INTEGER, INTENT(IN)                               :: nAtom,start,parent,step,csize
                CHARACTER(len=5), DIMENSION(:), INTENT(IN)        :: names
                INTEGER, INTENT(INOUT)                            :: ncycle,nelec,narom
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: discovered
                INTEGER                                           :: i,j

                ! Add the starting atom to the discovered list
                CALL add_element(discovered,start)

                ! Add the number of electrons corresponding to the 
                ! atomic type
                SELECT CASE(names(start))
                        CASE('C 3')
                                nelec = nelec+1
                        CASE('C_2')
                                nelec = nelec+1
                        CASE('C_1')
                                nelec = nelec+1
                        CASE('C_R')
                                nelec = nelec+1
                        CASE('C 2')
                                nelec = nelec+1
                        CASE('N 3')
                                nelec = nelec+2
                        CASE('N_R')
                                nelec = nelec+2
                        CASE('O 2')
                                nelec = nelec+2
                        CASE('O_2')
                                nelec = nelec+2
                        CASE('S 2')
                                nelec = nelec+2
                        CASE('S_R')
                                nelec = nelec+2
                        CASE DEFAULT
                                nelec = nelec
                END SELECT

                DO j=1,nAtom
                        ! If the "start" is connected to the "parent" we have found
                        ! a cycle
                        IF(D(start,j)==1 .AND. j==parent .AND. step==csize) THEN
                                ncycle = ncycle+1
                                IF(MOD(nelec-2,4)==0 .AND. nelec > 4) THEN
                                        narom = narom + 1
                                END IF
                        ! If the "start" isn't connected to any atom, 
                        ! we remove the added electron and the "start"
                        ! from "discovered"
                        ELSE IF(j == nAtom) THEN
                                SELECT CASE(names(start))
                                        CASE('C 3')
                                                nelec = nelec-1
                                        CASE('C_R')
                                                nelec = nelec-1
                                        CASE('C_2')
                                                nelec = nelec-1
                                        CASE('C_1')
                                                nelec = nelec-1
                                        CASE('C 2')
                                                nelec = nelec-1
                                        CASE('N 3')
                                                nelec = nelec-2
                                        CASE('N_R')
                                                nelec = nelec-2
                                        CASE('O 2')
                                                nelec = nelec-2
                                        CASE('O_2')
                                                nelec = nelec-2
                                        CASE('S 2')
                                                nelec = nelec-2
                                        CASE('S_R')
                                                nelec = nelec-2
                                        CASE DEFAULT
                                                nelec = nelec
                                END SELECT
                                discovered = discovered(:SIZE(discovered)-1)
                        ! If the "start" is connected to an atom which is not
                        ! the "parent", we start DFS again from this atom
                        ELSE IF(D(start,j)==1 .AND. step+1 <= csize .AND. .NOT. ANY(j==discovered)) THEN
                                CALL find_cycle(D,nAtom,discovered,j,parent,step+1,ncycle,csize,nelec,narom,names)
                        END IF
                END DO
        END SUBROUTINE
        
        ! INPUT     : 
        !               - D, matrix of integer, the adjacency matrix
        !               - names, list of charcater, the names of the atoms with the
        !                 format 'xxV  ' xx represent the atom and V the number of
        !                 neighbours
        !               - nAtom : integer, number of atoms
        ! OPERATION :
        !               For each atom, change the name to assigne an UFF atom type
        ! OUTPUT    :
        !               - names, the list of atom UFF types
        SUBROUTINE atom_type_assign(D,names,nAtom)
                INTEGER, DIMENSION(:,:),INTENT(IN)           :: D
                INTEGER, INTENT(IN)                          :: nAtom
                CHARACTER(len=5), DIMENSION(:),INTENT(INOUT) :: names
                INTEGER                                      :: i,j,nv,ncycle,narom,nelec,dox
                CHARACTER(len=2)                             :: na
                INTEGER, DIMENSION(:), ALLOCATABLE           :: discovered
                LOGICAL                                      :: arom
                
                !Put the number of neighbours in the name
                CALL atom_neighbour(D,names,nAtom)


                ! Loop over all the atoms and assigne the UFF type
                DO i=1,nAtom
                        READ(names(i),'(a2,i1)') na,nv
                        SELECT CASE(na)
                                CASE('H')
                                        IF(nv==1) THEN
                                                names(i) = 'H_   '
                                        ELSE
                                                names(i) = 'Hb   '
                                        END IF
                                CASE('He')
                                        names(i) = 'He4+4'
                                CASE('Li')
                                        names(i) = 'Li   '
                                CASE('Be')
                                        names(i) = 'Be3+2'
                                CASE('B')
                                        IF(nv==3) THEN
                                                names(i) = 'B_2  '
                                        ELSE IF(nv==4) THEN
                                                names(i) = 'B_3  '
                                        ELSE
                                                PRINT *,"Unknown boron (B) type : ",i
                                        END IF
                                CASE('C')
                                        arom = .FALSE.
                                        IF(nv==4) THEN
                                                names(i) = 'C_3  '
                                        ELSE
                                                DO j=3,7
                                                narom=0
                                                ncycle=0
                                                nelec=0
                                                CALL find_cycle(D,nAtom,discovered,i,i,1,ncycle,j,nelec,narom,names)
                                                DEALLOCATE(discovered)
                                                IF(narom>0) THEN
                                                        arom = .TRUE.
                                                END IF
                                                END DO
                                                IF(arom) THEN
                                                        names(i)='C_R  '
                                                ELSE IF(nv==3) THEN
                                                        names(i)='C_2  '
                                                ELSE IF(nv==2) THEN
                                                        names(i)='C_1  '
                                                END IF
                                        END IF
                                CASE('N')
                                        arom = .FALSE.
                                        DO j=3,7
                                        narom = 0
                                        ncycle = 0
                                        nelec = 0
                                        CALL find_cycle(D,nAtom,discovered,i,i,1,ncycle,j,nelec,narom,names)
                                        DEALLOCATE(discovered)
                                        IF(narom>0) THEN
                                                arom = .TRUE.
                                        END IF
                                        END DO
                                        IF(arom) THEN
                                                names(i)='N_R  '
                                        ELSE IF(nv==1) THEN
                                                names(i)='N_1  '
                                        ELSE IF(nv==2) THEN
                                                names(i)='N_2  '
                                        ELSE IF(nv==3) THEN
                                                names(i)='N_3  '
                                        END IF
                                CASE('O')
                                        arom = .FALSE.
                                        DO j=3,7
                                        narom = 0
                                        ncycle = 0
                                        nelec = 0
                                        CALL find_cycle(D,nAtom,discovered,i,i,1,ncycle,j,nelec,narom,names)
                                        DEALLOCATE(discovered)
                                        IF(narom>0) THEN
                                                arom = .TRUE.
                                        END IF
                                        END DO
                                        IF(arom) THEN
                                                names(i)='O_R  '
                                        ELSE IF(nv==1) THEN
                                                names(i)='0_1  '
                                        ELSE IF(nv==2) THEN
                                                names(i)='0_3  '
                                        END IF
                                CASE('F')
                                        names(i) = 'F_   '
                                CASE('Ne')
                                        names(i) = 'Ne4+4'
                                CASE('Na')
                                        names(i) = 'Na   '
                                CASE('Mg')
                                        names(i) = 'Mg3+2'
                                CASE('Al')
                                        names(i) = 'Al3  '
                                CASE('Si')
                                        names(i) = 'Si3  '
                                CASE('P')
                                        dox=0
                                        DO j=1,nAtom
                                                IF(D(i,j)==1) THEN
                                                        IF(names(j) == 'H' .OR. names(j) == 'H_   ' .OR. names(j) == 'Hb   ') THEN
                                                                dox = dox - 1
                                                        ELSE IF(names(j) == 'O' .OR. names(j) == 'O_1  ' .OR. names(j) == 'O_2  '&
                                                        &.OR. names(j) == 'O_3  ' .OR. names(j) == 'O_R  ') THEN
                                                        dox = dox + 2
                                                       END IF
                                                END IF
                                        END DO
                                        IF(dox == 3) THEN
                                                names(i) = 'P_3+3'
                                        ELSE
                                                names(i) = 'P_3+5'
                                        END IF
                                CASE('S')
                                        dox = 0
                                        arom = .FALSE.
                                        DO j=3,7
                                        narom = 0
                                        ncycle = 0
                                        nelec = 0
                                        CALL find_cycle(D,nAtom,discovered,i,i,1,ncycle,j,nelec,narom,names)
                                        DEALLOCATE(discovered)
                                        IF(narom>0) THEN
                                                arom = .TRUE.
                                        END IF
                                        END DO
                                        DO j=1,nAtom
                                                IF(D(i,j)==1) THEN
                                                        IF(names(j) == 'H' .OR. names(j) == 'H_   ' .OR. names(j) == 'Hb   ') THEN
                                                                dox = dox - 1
                                                        ELSE IF(names(j) == 'O' .OR. names(j) == 'O_1  ' .OR. names(j) == 'O_2  '&
                                                        &.OR. names(j) == 'O_3  ' .OR. names(j) == 'O_R  ') THEN
                                                        dox = dox + 2
                                                       END IF
                                                END IF
                                        END DO
                                        IF(arom) THEN
                                                names(i) = 'S_R  '
                                        ELSE IF(dox == 2) THEN
                                                names(i) = 'S_3+2'
                                        ELSE IF(dox == 6) THEN
                                                names(i) = 'S_3+6'
                                        ELSE
                                                names(i) = 'S_3+4'
                                        END IF
                                CASE('Cl')
                                        names(i) = 'Cl   '
                                CASE('Ar')
                                        names(i) = 'Ar4+4'
                                CASE('K')
                                        names(i) = 'K_   '
                                CASE('Ca')
                                        names(i) = 'Ca6+2'
                                CASE('Sc')
                                        names(i) = 'Sc3+3'
                                CASE('Ti')
                                        IF(nv==6) THEN
                                                names(i) = 'Ti6+4'
                                        ELSE
                                                names(i) = 'Ti3+4'
                                        END IF
                                CASE('V')
                                        names(i) = 'V_3+5'
                                CASE('Cr')
                                        names(i) = 'Cr6+3'
                                CASE('Mn')
                                        names(i) = 'Mn6+2'
                                CASE('Fe')
                                        IF(nv==6) THEN
                                                names(i) = 'Fe6+2'
                                        ELSE
                                                names(i) = 'Fe3+2'
                                        END IF
                                CASE('Co')
                                        names(i) = 'Co6+3'
                                CASE('Ni')
                                        names(i) = 'Ni4+2'
                                CASE('Cu')
                                        names(i) = 'Cu3+1'
                                CASE('Zn')
                                        names(i) = 'Zn3+2'
                                CASE('Ga')
                                        names(i) = 'Ga3+3'
                                CASE('Ge')
                                        names(i) = 'Ge3  '
                                CASE('As')
                                        names(i) = 'As3+3'
                                CASE('Se')
                                        names(i) = 'Se3+2'
                                CASE('Br')
                                        names(i) = 'Br   '
                                CASE('Kr')
                                        names(i) = 'Kr4+4'
                                CASE('Rb')
                                        names(i) = 'Rb   '
                                CASE('Sr')
                                        names(i) = 'Sr6+2'
                                CASE('Y')
                                        names(i) = 'Y_3+3'
                                CASE('Zr')
                                        names(i) = 'Zr3+4'
                                CASE('Nb')
                                        names(i) = 'Nb3+5'
                                CASE('Mo')
                                        IF(nv==6) THEN
                                                names(i) = 'Mo6+6'
                                        ELSE
                                                names(i) = 'Mo3+6'
                                        END IF
                                CASE('Tc')
                                        names(i) = 'Tc6+5'
                                CASE('Ru')
                                        names(i) = 'Ru6+2'
                                CASE('Rh')
                                        names(i) = 'Rh6+3'
                                CASE('Pd')
                                        names(i) = 'Pd4+2'
                                CASE('Ag')
                                        names(i) = 'Ag1+1'
                                CASE('Cd')
                                        names(i) = 'Cd3+2'
                                CASE('In')
                                        names(i) = 'In3+3'
                                CASE('Sn')
                                        names(i) = 'Sn3  '
                                CASE('Sb')
                                        names(i) = 'Sb3+3'
                                CASE('Te')
                                        names(i) = 'Te3+2'
                                CASE('I')
                                        names(i) = 'I_   '
                                CASE('Xe')
                                        names(i) = 'Xe4+4'
                                CASE('Cs')
                                        names(i) = 'Cs   '
                                CASE('Ba')
                                        names(i) = 'Ba6+2'
                                CASE('La')
                                        names(i) = 'La3+3'
                                CASE('Ce')
                                        names(i) = 'Ce6+3'
                                CASE('Pr')
                                        names(i) = 'Pr6+3'
                                CASE('Nd')
                                        names(i) = 'Nd6+3'
                                CASE('Pm')
                                        names(i) = 'Pm6+3'
                                CASE('Sm')
                                        names(i) = 'Sm6+3'
                                CASE('Eu')
                                        names(i) = 'Eu6+3'
                                CASE('Gd')
                                        names(i) = 'Gd6+3'
                                CASE('Tb')
                                        names(i) = 'Tb6+3'
                                CASE('Dy')
                                        names(i) = 'Dy6+3'
                                CASE('Ho')
                                        names(i) = 'Ho6+3'
                                CASE('Er')
                                        names(i) = 'Er6+3'
                                CASE('Tm')
                                        names(i) = 'Tm6+3'
                                CASE('Yb')
                                        names(i) = 'Yb6+3'
                                CASE('Lu')
                                        names(i) = 'Lu6+3'
                                CASE('Hf')
                                        names(i) = 'Hf3+4'
                                CASE('Ta')
                                        names(i) = 'Ta3+5'
                                CASE('W')
                                        dox = 0
                                        IF(nv==6) THEN
                                                names(i) = 'W_6+6'
                                        ELSE
                                                DO j=1,nAtom
                                                        IF(D(i,j)==1) THEN
                                                                IF(names(j) == 'H' .OR. names(j) == 'H_   ' .OR. names(j) == 'Hb   &
                                                                        &') THEN
                                                                        dox = dox - 1
                                                                ELSE IF(names(j) == 'O' .OR. names(j) == 'O_1  ' .OR. names(j) == &
                                                                        &'O_2  ' .OR. names(j) == 'O_3  ' .OR. names(j) == &
                                                                        &'O_R  ') THEN
                                                                dox = dox + 2
                                                               END IF
                                                        END IF
                                                END DO
                                                IF(dox==4) THEN
                                                        names(i)='W_3+4'
                                                ELSE
                                                        names(i)='W_3+6'
                                                END IF
                                        END IF
                                CASE('Re')
                                        IF(nv==6) THEN
                                                names(i) = 'Re6+5'
                                        ELSE
                                                names(i) = 'Re3+7'
                                        END IF
                                CASE('Os')
                                        names(i) = 'Os6+6'
                                CASE('Ir')
                                        names(i) = 'Ir6+3'
                                CASE('Pt')
                                        names(i) = 'Pt4+2'
                                CASE('Au')
                                        names(i) = 'Au4+3'
                                CASE('Hg')
                                        names(i) = 'Hg1+2'
                                CASE('Tl')
                                        names(i) = 'Tl3+3'
                                CASE('Pb')
                                        names(i) = 'Pb3  '
                                CASE('Bi')
                                        names(i) = 'Bi3+3'
                                CASE('Po')
                                        names(i) = 'Po3+2'
                                CASE('At')
                                        names(i) = 'At   '
                                CASE('Rn')
                                        names(i) = 'Rn4+4'
                                CASE('Fr')
                                        names(i) = 'Fr   '
                                CASE('Ra')
                                        names(i) = 'Ra6+2'
                                CASE('Ac')
                                        names(i) = 'Ac6+3'
                                CASE('Th')
                                        names(i) = 'Th6+4'
                                CASE('Pa')
                                        names(i) = 'Pa6+4'
                                CASE('U')
                                        names(i) = 'U_6+4'
                                CASE('Np')
                                        names(i) = 'Np6+4'
                                CASE('Pu')
                                        names(i) = 'Pu6+4'
                                CASE('Am')
                                        names(i) = 'Am6+4'
                                CASE('Cm')
                                        names(i) = 'Cm6+3'
                                CASE('Bk')
                                        names(i) = 'Bk6+3'
                                CASE('Cf')
                                        names(i) = 'Cf6+3'
                                CASE('Es')
                                        names(i) = 'Es6+3'
                                CASE('Fm')
                                        names(i) = 'Fm6+3'
                                CASE('Md')
                                        names(i) = 'Md6+3'
                                CASE('No')
                                        names(i) = 'No6+3'
                                CASE('Lw')
                                        names(i) = 'Lw6+3'
                                CASE DEFAULT
                                        PRINT *,"Don't know this atom : ",na
                        END SELECT
                END DO
        END SUBROUTINE

        ! INPUT     :
        !               - D, matrix of integer, the adjacency matrix
        !               - names, list of characters, the names of the atoms
        !               - B, matrix of real, bond order matrix
        ! OPERATION :
        !               Use the "standard" atomic valence to compute the
        !               bond order matrix, inspired by :
        !       Kim Y., Kim W. Y.; Bull. Koeran Chem. Soc. 2015, 36, 1769-1777 
        ! OUTPUT    :
        !               - B, matrix of integer, bond order matrix
        !       B(i,j) = 0 if i and j are not connected
        !       B(i,j) = bond order if i and j are connected
        !       i,j from 1 to the number of atoms
        SUBROUTINE bond_order(D,names,B)
                USE ArrayManip

                CHARACTER(len=5), DIMENSION(:), INTENT(IN)    :: names
                INTEGER, DIMENSION(:,:), INTENT(INOUT)        :: D
                INTEGER, DIMENSION(:), ALLOCATABLE            :: UA,DU
                REAL, DIMENSION(:,:), ALLOCATABLE             :: B
                INTEGER                                       :: i, j,nv

                ALLOCATE(B(SIZE(D),SIZE(D)))
                ALLOCATE(DU(SIZE(D)))

                DU = 0
                B = D

                DO i=1,SIZE(names)
                        nv = SUM(D(:,i))
                        SELECT CASE(names(i))
                                CASE('H')
                                        DU(i) = 1 - nv
                                CASE('B')
                                        DU(i) = 3 - nv
                                CASE('C')
                                        DU(i) = 4 - nv
                                CASE('O')
                                        DU(i) = 2 - nv
                                CASE('N')
                                        DU(i) = 3 - nv
                                CASE('P')
                                        DU(i) = 3 - nv
                                CASE('S')
                                        DU(i) = 2 - nv
                                CASE DEFAULT
                                        DU(i) = 0
                        END SELECT
                END DO
                
                i=1
                DO WHILE(SUM(DU)/=0)
                       IF(DU(i)==0) THEN
                               i = i+1
                       ELSE 
                               DO j=1,SIZE(DU)
                                       IF(DU(i)>0 .AND. DU(j)>0 .AND. D(i,j) > 0) THEN
                                               B(i,j) = B(i,j) + 1
                                               B(j,i) = B(j,i) + 1
                                               DU(i) = DU(i)-1
                                               DU(j) = DU(j)-1
                                       END IF
                               END DO
                       END IF
                       IF(i>SIZE(DU)) THEN
                               EXIT
                       END IF
                END DO
        END SUBROUTINE

        ! INPUT : 
        !               -names : character(len=5), names of the atoms
        !               -B     : matrix of real, bond order matrix
        ! OPERATION :
        !               Changes the bond order between resonant atoms
        !               to 1.5
        ! OUTPUT    :
        !               -B : bond order matrix        

        SUBROUTINE bond_order_arom(names,B)
                CHARACTER(len=5), DIMENSION(:), INTENT(IN)        :: names
                REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)  :: B
                CHARACTER(len=1)                                  :: hi,hj
                CHARACTER(len=2)                                  :: dummy
                INTEGER                                           :: i,j

                DO i=1,SIZE(names)
                        DO j=1,SIZE(names)
                                READ(names(i),'(2a,1a)') dummy,hi
                                READ(names(j),'(2a,1a)') dummy,hj
                                IF(hi=="R".AND.hj=="R".AND.B(i,j)/=0.0) THEN
                                        B(i,j)=1.5
                                END IF
                        END DO
                END DO
        END SUBROUTINE
END MODULE


