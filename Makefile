
default: main

main: partition-sfc.o partition-sfc.h main.o
	mpicc  -Wall -Wextra -o main main.c partition-sfc.c

partition-sfc.o: partition-sfc.c partition-sfc.h
	mpicc  -Wall -Wextra -c partition-sfc.c

main.o: main.c partition-sfc.h
	mpicc  -Wall -Wextra -c main.c
clean:
	rm -fv *.o main
