CC= mpicc
CFLAGS= -O3 # -flibrary-level 2 -bspfifo 10000  -fcombine-puts -fcombine-puts-buffer 256K,128M,4K
LFLAGS= -lm -lbsponmpi

all: bspinprod.o bspmv.o

bspinprod.o: bspinprod.c bspedupack.h
	$(CC) $(CFLAGS) -c bspinprod.c

bspmv.o: bspmv.c bspedupack.h
	$(CC) $(CFLAGS) -c bspmv.c

clean:
	rm -f *.o
