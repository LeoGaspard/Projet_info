MODULE math
        IMPLICIT NONE
        CONTAINS

        ! INPUT     : 
        !               -a, b  : Real vectors, positions of objects a and b
        ! OPERATION :
        !               Compute the distance between a and b
        ! RETURN    :
        !               Real, the distance
        REAL FUNCTION distance(a,b)
                REAL, DIMENSION(3), INTENT(IN) :: a,b

                distance =  SQRT((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
                RETURN 
        END FUNCTION

        ! INPUT     : 
        !               -posA, posB,posC  : Real vectors, positions of 
        !                objects a,b and c
        ! OPERATION :
        !               Compute the angle a-b-c
        ! RETURN    :
        !               Real, the angle
        REAL FUNCTION angle(posA,posB,posC)
                REAL, DIMENSION (3), INTENT(IN) :: posA,posB,posC
                REAL                            :: ab,ac,bc,x

                ab= distance(posA,posB)
                ac= distance(posA,posC)
                bc= distance(posB,posC)

                x=((-ac**2+ab**2+bc**2)/(2*ab*bc))
                angle= ACOS(x)
       END FUNCTION
        
       ! INPUT     :
       !                - vect_a, vect_b : Real vectors, positions of 
       !                  objects a and b
       ! OPERATION :
       !                Computes the cross product between vectors a
       !                and b
       ! RETURN    :
       !                Real, the cross product vector
       FUNCTION cross_product(vect_a,vect_b) RESULT(cp)
                REAL, DIMENSION(3), INTENT(IN) :: vect_a,vect_b
                REAL, DIMENSION(3)             :: cp

                cp=(/ vect_a(2)*vect_b(3)-vect_a(3)*vect_b(2),vect_a(3)*vect_b(1)-vect_a(1)*vect_b(3), &
                        &vect_a(1)*vect_b(2)-vect_a(2)*vect_b(1) /)
       END FUNCTION

       ! INPUT     :
       !                -a,b,c,d : Real vectors, positions of 
       !                 objects a,c,b,d
       ! OPERATION :
       !                Compute the dihedral angle a-b-c-d
       ! RETURN    :
       !                Real, the dihedral angle
       REAL FUNCTION dihedral(a,b,c,d)
                REAL, DIMENSION (3),INTENT(IN)  :: a,b,c,d
                REAL, DIMENSION (3)             :: vect_ab,vect_bc,vect_cd
                REAL, DIMENSION (3)             :: vect_n1,vect_n2
                REAL                            :: doprod

                vect_ab=(/ b(1)-a(1),b(2)-a(2),b(3)-a(3) /)
                vect_bc=(/ c(1)-b(1),c(2)-b(2),c(3)-b(3) /)
                vect_cd=(/ d(1)-c(1),d(2)-c(2),d(3)-c(3) /)

                vect_n1=cross_product(vect_ab,vect_bc)
                vect_n2=cross_product(vect_bc,vect_cd)
                doprod = REAL(NINT(DOT_PRODUCT(vect_n1,vect_n2)/(NORM2(vect_n1)*NORM2(vect_n2))*10000))/10000
                IF (DOT_PRODUCT(vect_bc,cross_product(vect_n1,vect_n2)) .LT. 0) THEN
                        dihedral=ACOS(doprod)
                ELSE
                        dihedral=-ACOS(doprod)
                END IF
       END FUNCTION
END MODULE
