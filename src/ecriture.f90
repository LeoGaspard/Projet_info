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
                SUBROUTINE write_header(outfile)
                        INTEGER :: outfile,i
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,'(52a)') 'MOLECULAR MECHANICS SIMULATIONS FOR A SCHOOL PROJECT'
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'Implementation made by :'
                        WRITE(outfile,*) '                        - Yann Damour : yann.damour@univ-tlse3.fr'
                        WRITE(outfile,*) '                        - Pierre Racine : pierre.racine@univ-tlse3.fr'
                        WRITE(outfile,*) '                        - Léo Gaspard : leo.gaspard@univ-tlse3.fr'                    
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)

                        WRITE(outfile,*) 'UFF implementation following :'
                        WRITE(outfile,*) 'A. K. Rappé, C. J. Casewit, K. S. Colwell,'
                        WRITE(outfile,*) 'W. A. Goddard III, W. M. Skiff;'
                        WRITE(outfile,*) 'J. Am. Chem. Soc. 1992, 114, 10024-10035'
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                END SUBROUTINE
                SUBROUTINE write_energies(outfile,benergy,aenergy,tenergy,ienergy,vdwenergy)
                        REAL    :: benergy,aenergy,tenergy,ienergy,vdwenergy
                        INTEGER :: outfile,i
                 

                        WRITE(outfile,*) 'Bond stretch contribution to the energy          :'
                        WRITE(outfile,'(f15.5)') benergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'Angle bend contribution to the energy            :'
                        WRITE(outfile,'(f15.5)') aenergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'Torsion contribution to the energy               :'
                        WRITE(outfile,'(f15.5)') tenergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'Inversion contribution to the energy             :'
                        WRITE(outfile,'(f15.5)') ienergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'Van der Waals contribution to the energy         :'
                        WRITE(outfile,'(f15.5)') vdwenergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                        WRITE(outfile,*) 'TOTAL ENERGY                                     :'
                        WRITE(outfile,'(f15.5)') benergy+aenergy+tenergy+ienergy+vdwenergy
                        DO i=1,150
                                WRITE (outfile,'(a1)',advance='no') "*"
                        END DO
                        WRITE (outfile,*)
                END SUBROUTINE

END MODULE
