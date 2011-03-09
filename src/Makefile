CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

all: cg

cg: bspedupack.o bspcg.o
	$(CC) $(CFLAGS) -o cg bspedupack.o bspcg.o $(LFLAGS)

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o cg
