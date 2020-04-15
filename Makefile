CC = gfortran
EXEC = FF_U_NOOB
LIBS = 
FLAGS = 


all: main.o
	$(CC) *.o -o $(EXEC) $(LIBS)

main.o : structure.o lecture.o ecriture.o math.o energie.o
	$(CC) src/main.f90 -c $(FLAGS)

structure.o : math.o lecture.o arraymanip.o
	$(CC) src/structure.f90 -c $(FLAGS)

lecture.o : 
	$(CC) src/lecture.f90 -c $(FLAGS)

math.o :
	$(CC) src/math.f90 -c $(FLAGS)

arraymanip.o : 
	$(CC) src/ArrayManip.f90 -c $(FLAGS)

ecriture.o :
	$(CC) src/ecriture.f90 -c $(FLAGS)

energie.o : math.o
	$(CC) src/energie.f90 -c $(FLAGS)


clear :
	rm -f *.o *.mod

mr_proper :
	rm -f *.o $(EXEC) *.mod
