
.SILENT:

#CFLAGS=-ansi -Wall -pedantic -O0 -g -m64 -fprefetch-loop-arrays -funroll-loops
#CMSG="Debug flags"
CFLAGS=-ansi -Wall -pedantic -O3 -m64 -fprefetch-loop-arrays -funroll-loops -ffast-math -march=native -mtune=native -DNDEBUG
CMSG="Performance flags"
SCFLAGS=${CFLAGS} -static
DCFLAGS=${CFLAGS} -fPIC
LIBNAME=sparselib
CURRENTVERSION=1.4.0
SHORTCURRENTVERSION=${CURRENTVERSION:%.=%}

FILE_SRC=FileToVT.hpp FileToVT.cpp mmio.cpp mmio.h
UTIL_SRC=util.cpp util.hpp
HILBERT_SRC=Matrix2HilbertCoordinates.cpp Matrix2HilbertCoordinates.hpp
BLOCKORDER_SRC=BlockOrderer.hpp SBDTree.cpp SBDTree.hpp MinCCS.hpp

FILE_C=$(filter %.cpp,${FILE_SRC})
UTIL_C=$(filter %.cpp,${UTIL_SRC})
HILBERT_C=$(filter %.cpp,${HILBERT_SRC})
BLOCKORDER_C=$(filter %.cpp,${BLOCKORDER_SRC})

FILE_DEPS=${FILE_SRC} ${FILE_C:%.cpp=%.o}
UTIL_DEPS=${UTIL_SRC} ${UTIL_C:%.cpp=%.o}
HILBERT_DEPS=${HILBERT_SRC} ${HILBERT_C:%.cpp=%.o}
BLOCKORDER_DEPS=${BLOCKORDER_SRC} ${BLOCKORDER_C:%.cpp=%.o}
DEPS=${FILE_DEPS} ${UTIL_DEPS} ${HILBERT_DEPS} ${BLOCKORDER_DEPS}

CLEAN=$(filter %.o,${DEPS})
CLEAN:=${CLEAN:%.o=%.static.o} ${CLEAN}

CC=g++

all: msg mm2htrp driver sbd2trp sbdmv mtx2crstrp

msg:
	echo "***Compiling using ${CMSG}***"

help:
	echo ""
	echo "DEFAULT TARGETS (make all)"
	echo "=========================="
	echo "The default target makes, next to the library (see next section), the following utilities:"
	echo "   1. 'mm2htrp'    Converts a .mtx file into a Hilbert-ordered .trp file"
	echo "   2. 'driver'     Calculates Ax (a SpMV); with A loaded in from an .mtx file"
	echo "   3. 'sbdmv'      Loads a .mtx file, containing a matrix in 2D SBD, in a hierarchical format and times SpMVs"
	echo "   4. 'sbd2trp'    Converts a .mtx file, containing a matrix in 2D SBD format, to an optimised .trp file"
	echo "   5. 'mtx2crstrp' Loads a matrix file, and reorders the triplets in CRS ordering"
	echo ""
	echo "Execute each utility without arguments for details."
	echo ""
	echo "LIBRARY COMPILATION (make all/libs)"
	echo "==============================="
	echo "'make libs' attempts to compile both a static and shared library."
	echo "('make all' implies 'make libs')"
	echo "The libraries (.a / .so) are put into the ./lib/ directory."
	echo "The header files required for the library are put into ./include/"
	echo ""
	echo "OTHER MAKE TARGETS:"
	echo "==================="
	echo "'make driver-static' - Compiles the static library, as well as the driver utility using the static library."
	echo "'make driver-shared' - Compiles the shared library, as well as the driver utility using the shared library."
	echo "'make mtx2crstrp'    - Compiles only the mtx2crstrp utility."
	echo "'make driver'        - Compiles the static and shared libraries, as well as the driver utilities."
	echo "'make mm2htrp'       - Compiles only the mm2htrp utility."
	echo "'make sbd2trp'       - Compiles only the sbd2trp utility."
	echo "'make sbdmv'         - Compiles only the sbdmv utility."
	echo "'make shared'        - Only compiles the shared library (.so)"
	echo "'make static'        - Only compiles the static library (.a)"
	echo "'make intel'         - Does 'make all libs' using the intel compiler (icc) targetting intel architectures."
	echo "'make core2'         - Does 'make all libs' using the intel compiler (icc) targetting core2 architectures."
	echo "'make docs'          - Compiles javadoc documentation; will be generated in the ./doxy_out/ directory."
	echo ""
	echo "CLEANUP TARGETS:"
	echo "================"
	echo "'make clean'     - Cleans all .o files, and any compiled test executables."
	echo "'make veryclean' - Implies 'make clean', also removes all compiled executables, the library directories './lib' and './include', and the doxygen directory."
	echo ""

intel:
	make CC=icc CFLAGS="-O2 -DNDEBUG -fprefetch-loop-arrays -funroll-loops -ansi -fPIC -m64 -msse3 -lrt" all libs

core2:
	make CC=icc CFLAGS="-O2 -DNDEBUG -fprefetch-loop-arrays -funroll-loops -ansi -fPIC -m64 -march=core2 -mtune=core2 -lrt" all libs

libs: clean
	make CFLAGS="${CFLAGS} -DNDEBUG" static shared

driver-static: static driver.static.o
	${CC} ${CFLAGS} -o driver-static driver.static.o lib/lib${LIBNAME}${CURRENTVERSION}.a -lrt

driver-shared: shared driver.o
	${CC} ${CFLAGS} -o driver-shared driver.o -Llib/ lib/lib${LIBNAME}.so.${CURRENTVERSION} -lrt

driver:
	make driver-static
	make driver-shared

%.static.o: %.c %.h
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.static.o: %.cpp %.h
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.static.o: %.cpp %.hpp
	${CC} ${SCFLAGS} -c -o $@ ${^:%.hpp=}

%.static.o: %.cpp
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.c %.h
	${CC} ${DCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.cpp %.h
	${CC} ${DCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.cpp %.hpp
	${CC} ${DCFLAGS} -c -o $@ ${^:%.hpp=}

%.o: %.cpp
	${CC} ${DCFLAGS} -c -o $@ $^

static: ${DEPS:%.o=%.static.o}
	mkdir include || true
	cp *.hpp include/
	cp *.h include/
	mkdir lib || true;
	ar rcs lib${LIBNAME}${CURRENTVERSION}.a $(filter %.o,$^)
	mv lib${LIBNAME}${CURRENTVERSION}.a ./lib/

shared: ${DEPS}
	mkdir include || true
	cp *.hpp include/
	cp *.h include/
	mkdir lib || true;
	gcc 	${CFLAGS} \
		-shared -Wl,-soname,lib${LIBNAME}.so.${SHORTCURRENTVERSION} \
		-o ./lib/lib${LIBNAME}.so.${CURRENTVERSION} $(filter %.o,$^) -lc

docs: doxy.conf
	doxygen doxy.conf

test: test_crs test_icrs test_ts test_hts test_read mm2htrp test_zz_crs test_zz_icrs test_bicrs test_ddmatrix m2h_test
	echo "*** Running on randomly generated matrix, same seed"
	echo "CRS:"
	/usr/bin/time -f "%S %U" ./test_crs
	echo "ICRS:"
	/usr/bin/time -f "%S %U" ./test_icrs
	echo "TS:"
	/usr/bin/time -f "%S %U" ./test_ts
	echo "Hilbert TS:"
	/usr/bin/time -f "%S %U" ./test_hts
	echo "Zig-zag ICRS:"
	/usr/bin/time -f "%S %U" ./test_zz_icrs
	echo "TS read:"
	/usr/bin/time -f "%S %U" ./test_read
	echo
	echo "*** Running on fixed, hardcoded matrix"
	echo "Bi-ICRS:"
	/usr/bin/time -f "%S %U" ./test_bicrs
	echo "Dense Diagonal:"
	/usr/bin/time -f "%S %U" ./test_ddmatrix
	echo
	echo "*** Running on test.mtx"
	echo "Matrix Market to Hilbert Triplet test:"
	/usr/bin/time -f "%S %U" ./mm2htrp test.mtx
	mv test.mtx.hilbert.trp test.trp
	echo "TS read:"
	/usr/bin/time -f "%S %U" ./test_read
	ls -lh temp.trp
	rm -f  temp.trp
	make clean
	echo "Execute ./m2h_test to manually check 2D integer to 1D Hilbert conversion."

test_crs: tests/test_crs.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_crs.cpp -o test_crs

test_icrs: tests/test_icrs.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_icrs.cpp -o test_icrs

test_ts: tests/test_ts.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_ts.cpp -o test_ts

test_read: tests/test_read.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_read.cpp -o test_read

test_bicrs: tests/test_bicrs.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_bicrs.cpp -o test_bicrs

test_ddmatrix: tests/test_ddmatrix.cpp
	${CC} -Wall -O2 -ansi -g -I. tests/test_ddmatrix.cpp -o test_ddmatrix

test_hts: tests/test_hts.cpp ${HILBERT_DEPS}
	${CC} -Wall -O2 -ansi -g -I. -o test_hts tests/test_hts.cpp $(filter %.o,$^)

test_zz_crs: tests/test_zz_crs.cpp ${FILE_DEPS}
	g++ -Wall -O2 -ansi -I. -o test_zz_crs tests/test_zz_crs.cpp $(filter %.0,$^)

test_zz_icrs: tests/test_zz_icrs.cpp ${FILE_DEPS}
	g++ -Wall -O3 -ansi -g -I. -o test_zz_icrs tests/test_zz_icrs.cpp $(filter %.o,$^)

m2h_test: tests/m2h_test.cpp
	${CC} -Wall -O2 -ansi -I. -o m2h_test tests/m2h_test.cpp Matrix2HilbertCoordinates.cpp

mm2htrp: mm2htrp.cpp
	${CC} -Wall -O3 -ansi -o mm2htrp mm2htrp.cpp Matrix2HilbertCoordinates.cpp FileToVT.cpp mmio.cpp

sbd2trp: sbd2trp-main.cpp ${FILE_DEPS} Triplet.hpp ${BLOCKORDER_DEPS} ${UTIL_DEPS}
	${CC} ${CFLAGS} -o sbd2trp sbd2trp-main.cpp $(filter %.o,$^)

sbdmv: sbdmv.cpp ${FILE_DEPS} Triplet.hpp ${HILBERT_DEPS} ${UTIL_DEPS} ${BLOCKORDER_DEPS}
	${CC} ${CFLAGS} -o sbdmv sbdmv.cpp $(filter %.o,$^) -lrt

mtx2crstrp: mtx2crstrp.cpp ${FILE_DEPS}
	${CC} ${CFLAGS} -o $@ mtx2crstrp.cpp $(filter %.o,$^)

clean:
	rm -rf test_crs test_icrs test_ts test_hts test_read test.trp test_zz_icrs test_bicrs test_zz_crs test_ddmatrix
	rm -rf driver*.o
	rm -rf $(filter %.o,${CLEAN})

veryclean: clean
	rm -f mm2htrp
	rm -f sbd2trp
	rm -f m2h_test
	rm -f mtx2crstrp
	rm -rf sbdmv
	rm -rf include
	rm -rf lib
	rm -rf doxy_out
	rm -f driver-static
	rm -f driver-shared

