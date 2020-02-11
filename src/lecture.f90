MODULE lecture
        IMPLICIT NONE
        CONTAINS
              
!! Variables :
!!              - fName     = name of the .xyz file
!!              - positions = dynamic 2D array which size is the number of atoms
!!              - names     = dynamic 1D array which size is the number of atoms
!!              - nAtom     = the number of atoms (read from the .xyz file)
!!
!! Operation :
!!              Opens the corresponding .xyz file  and read it to store the values
!!              of the positions and names of the different atoms
!!
!! Output         :
!!              - positions(i,j)
!!                      i c [1,nAtom]   ith atom
!!                      j c [1,3]       j is the x,y or z coordinate (respectively 1,2,3)
!!              - names(i)
!!                      i c [1,nAtom]   ith atom name (N,C,H,O...)

        SUBROUTINE read_xyz(fName,positions,names, nAtom)
        CHARACTER(len=32)                             :: fName
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: positions
        CHARACTER(len=2), DIMENSION(:), ALLOCATABLE   :: names
        INTEGER                                       :: i,nAtom

        OPEN(UNIT=9,FILE=fName)

        READ(9,*) nAtom
        ALLOCATE(positions(nAtom,3))
        ALLOCATE(names(nAtom))
        READ(9,*)
        DO i=1,nAtom
                READ(9,*) names(i),positions(i,1),positions(i,2),positions(i,3)
        END DO
        CLOSE(9)
        END SUBROUTINE

!! Variables      :
!!              -atomList = 1D array of 10 components which are the possible atom names
!!              -minMax   = 2D array of 2x10 components which are the minimum and maximum 
!!                          bond length
!!              
!! Operation :
!!              Opens the bond.dat file in the data directory, extracts the values for the 
!!              bond length of the 10 possible atoms
!!
!! Output         :
!!              - atomList(i)
!!                      i c [1,10] the ith atom name
!!              - minMax(i,j)
!!                      i c [1,10] the ith atom
!!                      j c [1,2]  the minimum (j=1) or minimum (j=2) value for the bond length

        SUBROUTINE read_atom_sphere(atomList,minMax)
        REAL, DIMENSION(10,2)           :: minMax
        CHARACTER(len=2), DIMENSION(10) :: atomList
        INTEGER                         :: n,i

        n = 10
        OPEN(UNIT=9,FILE="data/bond.dat") 

        DO i=1,n
                READ(9,*) atomList(i),minMax(i,1),minMax(i,2)                                                
        END DO
        CLOSE(9)
        END SUBROUTINE

!! Variables :
!!              - atomName : string of length 2, the atom name within gaff formalism
!!              - mass     : real number, will contain the mass of the atom
!!              - charge   : real number, will contain the charge of the atom
!!              - gaffdat  : integer, the number of the unit corresponding to the gaff2.dat file
!!
!! Operation :
!!              Parse the gaff2.dat file in order to get the mass and charge of the atom corresponding to the input
!!
!! Output :
!!              - mass   : real number, mass of the atom, 0 if the atom was not found
!!              - charge : real number, charge of the atom, 0 if the atom was not found

        SUBROUTINE get_atom_properties(atomName,mass,charge,gaffdat)
                CHARACTER(len=2) :: atomName, atomTest
                REAL             :: mass, charge
                INTEGER          :: i,gaffdat
                LOGICAL          :: success

                REWIND(gaffdat)
                success = .FALSE.

                READ(gaffdat,*)
                DO i=1,97
                        READ(gaffdat,*) atomTest, mass, charge
                        IF(atomTest==atomNAme) THEN
                                success = .TRUE.
                                EXIT
                        END IF
                END DO 
                IF(success.EQV..FALSE.) THEN
                        mass = 0
                        charge = 0
                END IF
        END SUBROUTINE

!! Variable :
!!              - atoma, atomb : string of length 2, name of the two atoms within the gaff formalism
!!              - k            : real, will contain the force constant
!!              - r            : real, will contain the equilibrium distance
!!              -gaffdat       : integer, the number of the unit corresponding to the gaff2.dat file
!!
!! Operation :
!!              Parses the gaff2.dat file in order to find the value of the force constant and equilibrium 
!!              values for the atoma-atomb bond
!!
!! Output :
!!              - k : real, the value of the force constant, 0 if the bond was not found
!!              - r : real, the value of the equilibrium distance, 0 if the bond was not found

        SUBROUTINE get_spring_properties(atoma, atomb, k, r,gaffdat)
                CHARACTER(len=2) :: atoma, atomb
                CHARACTER(len=5) :: spring,test
                REAL             :: k,r
                INTEGER          :: i,gaffdat
                LOGICAL          :: success 

                REWIND(gaffdat)
                success = .FALSE.
                spring = atoma//"-"//atomb

                READ(gaffdat,*) test
                DO i=0,97
                        READ(gaffdat,*)
                END DO
                READ(gaffdat,*) test
                DO i=1,935
                        READ(gaffdat,'(a5,f7.1,f10.4)') test,k,r
                        IF(test==spring) THEN
                                success = .TRUE.
                                EXIT
                        END IF
                END DO
                IF(success.EQV..FALSE.) THEN
                        REWIND(gaffdat)
                        DO i=0,97
                                READ(gaffdat,*)
                        END DO
                        READ(gaffdat,*)

                        spring = atomb//"-"//atoma
                        DO i=1,935
                                READ(gaffdat,'(a5,f7.1,f10.4)') test,k,r
                                IF(test==spring) THEN
                                        success = .TRUE.
                                        EXIT
                                END IF
                        END DO
                END IF

                IF(success.EQV..FALSE.) THEN
                        k = 0.0
                        r = 0.0
                END IF
        END SUBROUTINE

!! Variables :
!!              - atoma,atomb,atomc : string of length 2, the names of the atoms
!!              - k                 : real, will contain the value of the force constant
!!              - theta             : real, will contain the value of the equilibrium angle
!!              - gaffdat           : integer, the number of the unit corresponding to the gaff2.dat file
!!
!! Operation :
!!              Parses the gaff2.dat file in order to find the value of the force constant and equilibrium angle
!!              corresponding to the atoma-atomb-atomc angle
!!
!! Output :
!!              - k     : real, force constant, 0 if the angle was not found
!!              - theta : real, equilibrium angle, 0 if the angle was not found  
        SUBROUTINE get_angle_properties(atoma,atomb,atomc,k,theta,gaffdat)
                CHARACTER(len=2) :: atoma,atomb,atomc
                CHARACTER(len=8) :: angle,test
                REAL             :: k,theta
                INTEGER          :: i,gaffdat
                LOGICAL          :: success
                
                REWIND(gaffdat)
                success = .FALSE.
                angle = atoma//"-"//atomb//"-"//atomc

                DO i=0,1034
                        READ(gaffdat,*)
                END DO

                DO i=0,5312
                        READ(gaffdat,'(a8,f8.1,f12.2)') test,k,theta
                        IF(test==angle) THEN
                                success = .TRUE.
                                EXIT
                        END IF
                END DO
                IF(success.EQV..FALSE.) THEN

                        REWIND(gaffdat)
                        angle = atomc//"-"//atomb//"-"//atoma

                        DO i=0,1034
                                READ(gaffdat,*)
                        END DO

                        DO i=0,5312
                                READ(gaffdat,'(a8,f8.1,f12.2)') test,k,theta
                                IF(test==angle) THEN
                                        EXIT
                                END IF
                        END DO
                END IF
                IF(success.EQV..FALSE.) THEN
                        k = 0
                        theta = 0
                END IF
        END SUBROUTINE


END MODULE
