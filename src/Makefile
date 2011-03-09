CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspedupack.o bspcg.o paullib.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o

all: lib cg

lib: force-look
	cd libs; make all

cg: $(OBJS)
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LIBOBJS) $(LFLAGS)

paullib.o: paullib.c
	$(CC) $(CFLAGS) -c paullib.c

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: lib
	$(CC) $(CFLAGS) -c libs/bspedupack.c

clean:
	rm -f *.o cg
	cd libs; make clean

force-look:
	@true
