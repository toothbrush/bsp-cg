CC= mpicc
CFLAGS= -pg -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspedupack.o bspcg.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o libs/vecio.o libs/matsort.o libs/paullib.o

all: lib cg

lib: force-look
	cd libs; make all

cg: $(OBJS) lib
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LIBOBJS) $(LFLAGS)

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: lib
	$(CC) $(CFLAGS) -c libs/bspedupack.c

clean:
	rm -f *.o cg
	cd libs; make clean

force-look:
	@true
