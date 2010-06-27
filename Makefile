rigidbodyeoms.o: rigidbodyeoms.c
	gcc -Wall -O3 -funroll-loops -c rigidbodyeoms.c

clean :
	rm *.o
