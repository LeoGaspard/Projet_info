!! THIS MODULE CONTAINS ALL THE SUBROUTINE THAT WILL WRITE
!! INTO THE OUTPUT FILE
MODULE ecriture
        IMPLICIT NONE
        CONTAINS
                SUBROUTINE write_positions(outfile,positions,names)
                        INTEGER :: outfile,i
                        REAL, DIMENSION(:,:) :: positions
                        CHARACTER(len=5), DIMENSION(:) :: names

                        WRITE (outfile,*) "Printout of the atomic coordinates (Angstrom) : " 
                        WRITE (outfile,*) 
                        DO i=1,SIZE(names)
                                WRITE (outfile,'(a2,f13.4,f13.4,f13.4)') names(i),positions(i,1),positions(i,2),positions(i,3)
                        END DO
                        WRITE (outfile,*) 
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                END SUBROUTINE
                SUBROUTINE write_connectivity(outfile,B)
                        REAL, DIMENSION(:,:) :: B
                        INTEGER :: i,j,outfile

                        WRITE (outfile,*) "Printout of the calculated connectivity :"
                        WRITE (outfile,*)
                        WRITE (outfile,*) "               Atom                 Atom             Bond order"
                        DO i=1,SIZE(B(1,:))
                                DO j=i,SIZE(B(1,:))
                                        IF(B(i,j)/=0) THEN
                                                WRITE (outfile,'(i20,i20,f20.1)') i,j,B(i,j)
                                        END IF
                                END DO
                        END DO
                        WRITE (outfile,*)
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                END SUBROUTINE
                SUBROUTINE write_uff_properties(outfile,types,prop)
                        REAL, DIMENSION(:,:), ALLOCATABLE :: prop
                        CHARACTER(len=5), DIMENSION(:), ALLOCATABLE :: types
                        INTEGER :: outfile,i

                        WRITE (outfile,*) "Printout of the used Universal Force Field parameters"
                        WRITE (outfile,*)
                        WRITE (outfile,*) "Atom          r1        theta0          x1           D1          zeta   &
                              &       Z1          Vi            Uj            Xi       Hardness     Radius"

                        DO i=1,SIZE(types)
                                WRITE (outfile,'(a5,5X,11(f8.3,5X))') types(i),prop(i,:)
                        END DO
                        WRITE (outfile,*)
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                END SUBROUTINE

END MODULE
