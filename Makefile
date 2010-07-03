all : simulate

simulate : simulate.o rigidbodyeoms.o
	gcc -Wall -O3 -funroll-loops -lGL -lGLU -lglut -lgsl -lcblas -latlas -lm -o simulate simulate.o rigidbodyeoms.o

rigidbodyeoms.o : rigidbodyeoms.c
	gcc -Wall -O3 -funroll-loops -c rigidbodyeoms.c

simulate.o : simulate.c
	gcc -Wall -O3 -funroll-loops -c simulate.c

clean :
	rm -f simulate *.o *.in *.dir
