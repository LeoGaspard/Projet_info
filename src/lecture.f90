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
        CHARACTER(len=32), INTENT(IN)                                :: fName
        REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)             :: positions
        CHARACTER(len=5), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: names
        INTEGER, INTENT(INOUT)                                       :: nAtom
        INTEGER                                                      :: i

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

!! Variables :
!!              - Atom : CHARACTER, the atom name before the atomtype assignment (C,B,N,H...)
!!              - codat : The unit number for the file containing the covalent radius data
!! Operation :
!!              Read the codat file to find the covalent radius of a given atome
!! Output :
!!              - REAL the value of the covalent radius


        REAL FUNCTION get_covalent_radius(atom,codat)
                CHARACTER(len=2), INTENT(IN) :: atom
                CHARACTER(len=2)             :: test
                INTEGER, INTENT(IN)          :: codat
                INTEGER                      :: i

                REWIND(codat)

                DO i=1,96
                        READ(codat,*) test,get_covalent_radius
                        IF(test==atom) THEN
                                EXIT
                        END IF
                END DO
        END FUNCTION

        ! INPUT     :
        !               - types : list of character, the different atomic types
        !               - prop  : matrix of real, to store the properties
        !               - uff   : integer, the unit number of the uff file
        ! OPERATION :
        !               Read in the UFF data file the lines corresponding 
        !               to the different 'types'
        ! OUTPUT    :
        !               -prop, a nAtom x 11 matrix containing the properties
        SUBROUTINE get_uff_properties(types,prop,uff)
                INTEGER, INTENT(IN)                        :: uff
                CHARACTER(len=5), DIMENSION(:), INTENT(IN) :: types
                REAL, DIMENSION(:,:), INTENT(INOUT)        :: prop
                INTEGER                                    :: i,j
                CHARACTER(len=5)                           :: test
                LOGICAL                                    :: success

                success = .FALSE.


                DO i=1,SIZE(types)
                        REWIND(uff)
                        READ(uff,*)
                        DO j=1,127
                                READ(uff,*) test,prop(i,1),prop(i,2)&
                                        &,prop(i,3),prop(i,4),prop(i,5),prop(i,6),prop(i,7),prop(i,8),prop(i,9),&
                                &prop(i,10),prop(i,11)
                                IF(test==types(i)) THEN
                                        success = .TRUE.
                                        EXIT
                                END IF
                        END DO
                        IF(success.EQV..FALSE.) THEN
                                prop(i,:) = 0
                                PRINT *,"Don't know this atom type : ",types(i)
                        END IF
                END DO
        END SUBROUTINE
END MODULE
