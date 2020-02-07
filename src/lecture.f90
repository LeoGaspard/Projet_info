MODULE lecture
        IMPLICIT NONE
        CONTAINS
        SUBROUTINE read_xyz(fName,positions,names)
        CHARACTER(len=32) :: fName
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: positions
        CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: names
        INTEGER :: i,nAtom

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
END MODULE
