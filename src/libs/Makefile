CC= mpicc
CFLAGS= -DDEBUG -g -Wall # -flibrary-level 2 -bspfifo 10000  -fcombine-puts -fcombine-puts-buffer 256K,128M,4K
LFLAGS= -lm -lbsponmpi

all: bspinprod.o bspmv.o vecio.o matsort.o paullib.o

matsort.o: matsort.h matsort.c
	$(CC) $(CFLAGS) -c matsort.c

paullib.o: paullib.h paullib.c
	$(CC) $(CFLAGS) -c paullib.c

vecio.o: vecio.c vecio.h
	$(CC) $(CFLAGS) -c vecio.c

bspinprod.o: bspinprod.c bspedupack.h bspfuncs.h
	$(CC) $(CFLAGS) -c bspinprod.c

bspmv.o: bspmv.c bspedupack.h bspfuncs.h
	$(CC) $(CFLAGS) -c bspmv.c

clean:
	rm -f *.o
