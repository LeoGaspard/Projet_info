MODULE ArrayManip
        IMPLICIT NONE
        ! INPUT     :
        !               - element : str or int, the element to be added
        !               - list    : allocatable, str or int, the list
        ! OPERATION :
        !               Add the "element" to the "list"
        ! OUTPUT    :
        !               - list, the list containing the old + new elements
        INTERFACE add_element
                MODULE PROCEDURE add_str,add_int
        END INTERFACE
        CONTAINS
                SUBROUTINE add_str(list,element)
                        INTEGER                                                    :: i,lsize
                        CHARACTER(len=5), INTENT(IN)                               :: element
                        CHARACTER(len=5), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: list
                        CHARACTER(len=5), DIMENSION(:), ALLOCATABLE                :: tmpList

                        IF(ALLOCATED(list)) THEN
                                lsize = SIZE(list)
                                ALLOCATE(tmpList(lsize+1))
                                DO i=1,lsize
                                        tmpList(i) = list(i)
                                END DO
                                tmpList(lsize+1) = element
                                DEALLOCATE(list)
                                CALL MOVE_ALLOC(tmpList,list)
                        ELSE
                                ALLOCATE(list(1))
                                list(1) = element
                        END IF
                END SUBROUTINE
                SUBROUTINE add_int(list,element)
                        INTEGER                                           :: i,lsize
                        INTEGER, INTENT(IN)                               :: element
                        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: list
                        INTEGER, DIMENSION(:), ALLOCATABLE                :: tmpList

                        IF(ALLOCATED(list)) THEN
                                lsize = SIZE(list)
                                ALLOCATE(tmpList(lsize+1))
                                DO i=1,lsize
                                        tmpList(i) = list(i)
                                END DO
                                tmpList(lsize+1) = element
                                DEALLOCATE(list)
                                CALL MOVE_ALLOC(tmpList,list)
                        ELSE
                                ALLOCATE(list(1))
                                list(1) = element
                        END IF
                END SUBROUTINE
END MODULE
