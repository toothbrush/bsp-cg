CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspedupack.o bspcg.o
OBJS2=bspedupack.o cg_test.o
OBJS_SEQ=seq.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o libs/vecio.o libs/matsort.o libs/paullib.o

all: lib cg cg2 seq

lib: force-look
	cd libs; make all

seq: $(OBJS_SEQ) lib
	gcc -o seq $(OBJS_SEQ) $(LIB_OBJS) -lm

cg2: $(OBJS2) lib
	$(CC) $(CFLAGS) -o cg2 $(OBJS2) $(LIBOBJS) $(LFLAGS)

cg: $(OBJS) lib
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LIBOBJS) $(LFLAGS)

cg_test.o: cg_test.c
	$(CC) $(CFLAGS) -c cg_test.c

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: lib
	$(CC) $(CFLAGS) -c libs/bspedupack.c

clean:
	rm -f *.o cg
	cd libs; make clean

force-look:
	@true
