CC= mpicc
CFLAGS= -g -Wall
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspedupack.o bspcg.o paullib.o 

all: cg

cg: $(OBJS)
	$(CC) $(CFLAGS) -o cg $(OBJS) $(LFLAGS)

paullib.o: paullib.c
	$(CC) $(CFLAGS) -c paullib.c

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o cg
