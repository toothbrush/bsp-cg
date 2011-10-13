CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspcg.o
OBJS_SEQ=seq.o
OBJS_GEN=genmat.o libs/vecalloc-seq.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o libs/vecio.o libs/matsort.o libs/paullib.o libs/bspedupack.o
BINS=cg genmat seq

all: lib $(BINS)

lib: force-look
	cd libs; make all

seq: $(OBJS_SEQ) $(LIBOBJS)
	gcc -o seq $(OBJS_SEQ) $(LIB_OBJS) -lm

genmat: $(OBJS_GEN) $(LIBOBJS)
	gcc -o genmat $(OBJS_GEN) $(LIB_OBJS) -lm

cg: $(OBJS) $(LIBOBJS)
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LIBOBJS) $(LFLAGS)

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

clean:
	rm -f *.o $(BINS)
	cd libs; make clean

force-look:
	@true
