CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspedupack.o bspcg.o
OBJS_SEQ=seq.o
OBJS_GEN=genmat.o vecalloc-seq.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o libs/vecio.o libs/matsort.o libs/paullib.o
BINS=cg genmat seq

all: lib $(BINS)

lib: force-look
	cd libs; make all

seq: $(OBJS_SEQ) lib
	gcc -o seq $(OBJS_SEQ) $(LIB_OBJS) -lm

genmat: lib $(OBJS_GEN)
	gcc -o genmat $(OBJS_GEN) $(LIB_OBJS) -lm

cg: $(OBJS) lib
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LIBOBJS) $(LFLAGS)

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: lib
	$(CC) $(CFLAGS) -c libs/bspedupack.c

vecalloc-seq.o:
	gcc -c libs/vecalloc-seq.c

clean:
	rm -f *.o $(BINS)
	cd libs; make clean

force-look:
	@true
