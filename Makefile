all : simulate

simulate : simulate.o rigidbodyeoms.o savepng.o
	gcc -Wall -O3 -funroll-loops -lGL -lGLU -lglut -lgsl -lpng -lcblas -latlas -lm -o simulate simulate.o rigidbodyeoms.o savepng.o

rigidbodyeoms.o : rigidbodyeoms.c
	gcc -Wall -O3 -funroll-loops -c rigidbodyeoms.c

simulate.o : simulate.c
	gcc -Wall -O3 -funroll-loops -c simulate.c

savepng.o : savepng.c
	gcc -Wall -O3 -funroll-loops -c savepng.c

clean :
	rm -f simulate *.o *.in *.dir
