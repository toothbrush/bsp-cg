include cc.mk
LFLAGS= -lm -lbsponmpi #-Wl,-rpath -Wl,LIBDIR

# the objects required to build the final executable CG
OBJS=bspcg.o
OBJS_SEQ=seq.o
OBJS_GEN=genmat.o libs/vecalloc-seq.o libs/paullib.o
LIBOBJS=libs/bspmv.o libs/bspinprod.o libs/vecio.o libs/matsort.o libs/paullib.o libs/bspedupack.o
BINDIR=../bin
BINS=cg genmat seq

all: lib $(BINS)

lib: force-look
	cd libs; make all

seq: $(OBJS_SEQ) $(LIBOBJS) $(BINDIR)
	gcc -o $(BINDIR)/seq $(OBJS_SEQ) $(LIB_OBJS) -lm

genmat: $(OBJS_GEN) $(LIBOBJS) $(BINDIR)
	gcc $(CFLAGS) -o $(BINDIR)/genmat $(OBJS_GEN) $(LIB_OBJS) -lm

cg: $(OBJS) $(LIBOBJS) $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/cg $(OBJS) $(LIBOBJS) $(LFLAGS)

$(BINDIR):
	mkdir $(BINDIR)

genmat.o: genmat.c genmat.h $(LIBOBJS)
	gcc $(CFLAGS) -c -o genmat.o genmat.c

bspcg.o: bspcg.c
	$(CC) $(CFLAGS) -c bspcg.c

clean:
	for i in $(BINS) ; do rm -fv $(BINDIR)/$$i; done
	rm -fv *.o
	cd libs; make clean

force-look:
	@true
