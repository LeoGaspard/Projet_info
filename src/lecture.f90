MODULE lecture
        IMPLICIT NONE
        CONTAINS
              
!! Variables :
!!              - fName     = nom du fichier .xyz
!!              - positions = tableau dynamique 2D dont la taille correspondra avec le nombre d'atomes
!!              - names     = tableau dynamique 1D dont la taille correspondra avec le nombre d'atomes
!!              - nAtom     = contiendra le nombre d'atomes
!!
!! Fonctionnement :
!!              Ouvre le fichier .xyz correspondant et le lis pour stocker les valeurs correspondantes
!!              aux positions ainsi qu'aux noms des atomes
!!
!! Output         :
!!              - positions(i,j)
!!                      i c [1,nAtom]   correspond à l'atome en question
!!                      j c [1,3]       correspond aux coordonnées x,y,z respectivement 1,2,3
!!              - names(i)
!!                      i c [1,nAtom]   correspond au nom de l'atom en question (N,C,H,O...)

        SUBROUTINE read_xyz(fName,positions,names, nAtom)
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

!! Variables      :
!!              -atomList = tableau 1D de dimension 10 qui contiendra les labels des atomes possibles
!!              -minMax   = tableau 2D de dimension 2,10 qui contiendra les valeurs minimales et 
!!                          maximales possibles pour les liaisons
!!              
!! Fonctionnement :
!!              Ouvre le fichier bonds.dat situé dans le répertoire data, en extrait les valeurs
!!              minimales et maximales de liaisons pour les 10 types d'atomes
!!
!! Output         :
!!              - atomList(i)
!!                      i c [1,10] correspond au nom de l'atome
!!              - minMax(i,j)
!!                      i c [1,10] correspond à l'atome en question
!!                      j c [1,2]  correspond aux valeurs minimales et maximales (respectivement i et j)
!!                      de liaisons faites

        SUBROUTINE read_atom_sphere(atomList,minMax)
        REAL, DIMENSION(10,2) :: minMax
        CHARACTER(len=2), DIMENSION(10) :: atomList
        INTEGER :: n,i

        n = 10
        OPEN(UNIT=9,FILE="data/bond.dat") 

        DO i=1,n
                READ(9,*) atomList(i),minMax(i,1),minMax(i,2)                                                
        END DO
        CLOSE(9)
        END SUBROUTINE
END MODULE
