CC= mpicc
CFLAGS= -Wall # -DDEBUG -g # -flibrary-level 2 -bspfifo 10000  -fcombine-puts -fcombine-puts-buffer 256K,128M,4K
LFLAGS= -lm -lbsponmpi

all: bspinprod.o bspmv.o vecio.o matsort.o paullib.o vecalloc-seq.o bspedupack.o

matsort.o: matsort.h matsort.c
	$(CC) $(CFLAGS) -c matsort.c

paullib.o: paullib.h paullib.c
	$(CC) $(CFLAGS) -c paullib.c

vecio.o: vecio.c vecio.h
	$(CC) $(CFLAGS) -c vecio.c

bspinprod.o: bspinprod.c bspedupack.h bspfuncs.h
	$(CC) $(CFLAGS) -c bspinprod.c

bspmv.o: bspmv.c bspedupack.o bspfuncs.h
	$(CC) $(CFLAGS) -c bspmv.c

vecalloc-seq.o: vecalloc-seq.h vecalloc-seq.c
	gcc -c vecalloc-seq.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o
