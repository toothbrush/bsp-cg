

/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#define FIXED_SEED 13925

/**Define OUTPUT_Z to have the output vector printed to stdout */
//#define OUTPUT_Z

#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>

#include "SparseMatrix.hpp"
#include "TS.hpp"
#include "CRS.hpp"
#include "ICRS.hpp"
#include "ZZ_CRS.hpp"
#include "ZZ_ICRS.hpp"
#include "SVM.hpp"
#include "HTS.hpp"
#include "BICRS.hpp"
#include "CCSWrapper.hpp"
#include "Hilbert.hpp"
#include "BlockHilbert.hpp"
#include "BisectionHilbert.hpp"

double checksum( double* z, unsigned long int m ) {
	if( m==0 ) return 0.0;
	double sum = z[ 0 ];
	for( unsigned long int i=1; i<m; i++ )
		sum += z[ i ];
	return sum / static_cast< double >( m );
}

Matrix< double >* selectMatrix( const int scheme, const int ccs, const std::string file ) {
	double zero = 0.0;
	size_t pos = file.find_last_of( '.' );
	std::string ext = file.substr( pos + 1, file.length() );
	std::cout << "Matrix file extension: \"" << ext << "\"." << std::endl;
	if( ext.compare( "trp" ) == 0 ) {
		std::cout << "Expecting binary triplet format, reading in..." << std::endl;
		unsigned long int m, n;
		std::vector< Triplet< double > > input = Triplet< double >::load( file, m, n );
		if( ccs ) {
			switch( scheme ) {
				case 0: { return new CCSWrapper< double, TS< double > >     ( input, m, n, zero ); }
				case 1: { return new CCSWrapper< double, CRS< double > >    ( input, m, n, zero ); }
				case 2: { return new CCSWrapper< double, ICRS< double > >   ( input, m, n, zero ); }
				case 3: { return new CCSWrapper< double, ZZ_CRS< double > > ( input, m, n, zero ); }
				case 4: { return new CCSWrapper< double, ZZ_ICRS< double > >( input, m, n, zero ); }
				case 5: { return new CCSWrapper< double, SVM< double > >    ( input, m, n, zero ); }
				case 6: { return new CCSWrapper< double, HTS< double > >    ( input, m, n, zero ); }
				case 7: { return new CCSWrapper< double, BICRS< double > >  ( input, m, n, zero ); }
				case 8: { return new CCSWrapper< double, Hilbert< double > >( input, m, n, zero ); }
				case 9: { return new CCSWrapper< double, BlockHilbert< double > >( input, m, n, zero ); }
				case 10:{ return new CCSWrapper< double, BisectionHilbert< double > >( input, m, n, zero ); }
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		} else {
			switch( scheme ) {
				case 0: { return new TS< double >     ( input, m, n, zero ); }
				case 1: { return new CRS< double >    ( input, m, n, zero ); }
				case 2: { return new ICRS< double >   ( input, m, n, zero ); }
				case 3: { return new ZZ_CRS< double > ( input, m, n, zero ); }
				case 4: { return new ZZ_ICRS< double >( input, m, n, zero ); }
				case 5: { return new SVM< double >    ( input, m, n, zero ); }
				case 6: { return new HTS< double >    ( input, m, n, zero ); }
				case 7: { return new BICRS< double >  ( input, m, n, zero ); }
				case 8: { return new Hilbert< double >( input, m, n, zero ); }
				case 9: { return new BlockHilbert< double >( input, m, n, zero ); }
				case 10:{ return new BisectionHilbert< double >( input, m, n, zero ); }
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		}
	} else /*if( ext.compare( "mtx" ) == 0 )*/ {
		std::cout << "Matrix-market format expected, reading in..." << std::endl;
		if( ccs ) {
			switch( scheme ) {
				case 0: { return new CCSWrapper< double, TS< double > >     ( file, zero ); }
				case 1: { return new CCSWrapper< double, CRS< double > >    ( file, zero ); }
				case 2: { return new CCSWrapper< double, ICRS< double > >   ( file, zero ); }
				case 3: { return new CCSWrapper< double, ZZ_CRS< double > > ( file, zero ); }
				case 4: { return new CCSWrapper< double, ZZ_ICRS< double > >( file, zero ); }
				case 5: { return new CCSWrapper< double, SVM< double > >    ( file, zero ); }
				case 6: { return new CCSWrapper< double, HTS< double > >    ( file, zero ); }
				case 7: { return new CCSWrapper< double, BICRS< double > >  ( file, zero ); }
				case 8: { return new CCSWrapper< double, Hilbert< double > >( file, zero ); }
				case 9: { return new CCSWrapper< double, BlockHilbert< double > >( file, zero ); }
				case 10:{ return new CCSWrapper< double, BisectionHilbert< double > >( file, zero ); }
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		} else {
			switch( scheme ) {
				case 0: { return new TS< double >     ( file, zero ); }
				case 1: { return new CRS< double >    ( file, zero ); }
				case 2: { return new ICRS< double >   ( file, zero ); }
				case 3: { return new ZZ_CRS< double > ( file, zero ); }
				case 4: { return new ZZ_ICRS< double >( file, zero ); }
				case 5: { return new SVM< double >    ( file, zero ); }
				case 6: { return new HTS< double >    ( file, zero ); }
				case 7: { return new BICRS< double >  ( file, zero ); }
				case 8: { return new Hilbert< double >( file, zero ); }
				case 9: { return new BlockHilbert< double >( file, zero ); }
				case 10:{ return new BisectionHilbert< double >( file, zero ); }
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		}
	}
}

int main( int argc, char** argv ) {
	struct timespec start, stop; 
	double time;

#ifndef NDEBUG
	std::cout << "-->WARNING: COMPILED *WITH* ASSERTIONS!<--" << std::endl;
#endif
	
	if( argc<=3 ) {
		std::cout << "Usage: " << argv[0] << " <mtx> <scheme> <x> <REP>" << std::endl << std::endl;
		std::cout << "calculates Ax=y and reports average time taken as well as the mean of y." << std::endl;
		std::cout << "with\t\t <mtx> filename of the matrix A in matrix-market or binary triplet format." << std::endl;
		std::cout << "    \t\t <scheme> number of a sparse scheme to use, see below." << std::endl;
		std::cout << "    \t\t <x> 0 for taking x to be the 1-vector, 1 for taking x to be random (fixed seed)." << std::endl;
		std::cout << "    \t\t <REP> (optional, default is 1) number of repititions of the in-place multiplication." << std::endl;
		std::cout << std::endl << "Possible schemes:" << std::endl;
		std::cout << " 0: TS (triplet scheme)" << std::endl;
		std::cout << " 1: CRS (also known as CSR)" << std::endl;
		std::cout << " 2: ICRS (Incremental CRS)" << std::endl;
		std::cout << " 3: ZZ-CRS (Zig-zag CRS)" << std::endl;
		std::cout << " 4: ZZ-ICRS (Zig-zag ICRS)" << std::endl;
		std::cout << " 5: SVM (Sparse vector matrix)" << std::endl;
		std::cout << " 6: HTS (Hilbert-ordered triplet scheme)" << std::endl;
		std::cout << " 7: BICRS (Bi-directional Incremental CRS)" << std::endl;
		std::cout << " 8: Hilbert (Hilbert-ordered triplets backed by BICRS)" << std::endl;
		std::cout << " 9: Block Hilbert (Sparse matrix blocking, backed by Hilbert and HBICRS)" << std::endl;
		std::cout << "10: Bisection Hilbert (Sparse matrix blocking by bisection, backed by Hilbert and HBICRS)" << std::endl;
		std::cout << std::endl << "The in-place Ax=y calculation is preceded by a quasi pre-fetch." << std::endl;
		std::cout << "Add a minus sign before the scheme number to enable use of the CCS wrapper (making each CRS-based structure CCS-based instead)" << std::endl;
		std::cout << "Note: binary triplet format is machine dependent. ";
		std::cout << "Take care when using the same binary files on different machine architectures." << std::endl;
		return EXIT_FAILURE;
	}

	std::string file = std::string( argv[1] );
	int scheme = atoi( argv[2] );
	int ccs    = scheme < 0 ? 1 : 0;
	if( ccs ) scheme = -scheme;
	int x_mode = atoi( argv[3] );
	int rep = 1;
	if( argc >= 5 )
		rep = atoi( argv[4] );

	std::cout << argv[0] << " called with matrix input file " << file << ", scheme number ";
	std::cout << scheme << " and x being " << (x_mode?"random":"the 1-vector") << "." << std::endl;
	std::cout << "Number of repititions of in-place zax is " << rep << std::endl;

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start);
	Matrix< double >* matrix = selectMatrix( scheme, ccs, file );
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop);
	time  = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( matrix == NULL ) {
		std::cerr << "Error during sparse scheme loading, exiting." << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Matrix dimensions: " << matrix->m() << " times " << matrix->n() << "." << std::endl;
	std::cout << "Datastructure loading time: " << time << " ms." << std::endl << std::endl;

	srand( FIXED_SEED );
	double* x = new double[ matrix->n() ];
	for( unsigned long int j=0; j<matrix->n(); j++ )
		x[ j ] = x_mode?rand():1.0;

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );
	double* z = matrix->mv( x );
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop);
	time = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
#ifdef OUTPUT_Z
	for( unsigned long int j=0; j<matrix->m(); j++ )
		std::cout << z[ j ] << std::endl;
#endif
	std::cout << "out-of-place z=Ax: mean=" << checksum( z, matrix->m() ) << ", ";
	std::cout << "time= " << time << " ms." << std::endl;

	double *times = new double[ rep ];

	//Run rep*rep instances
	for( int run = 0; run < rep; run++ ) {
		//"prefetch"
		matrix->zax( x, z );
		clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start);
		for( int j=0; j<rep; j++ )
			matrix->zax( x, z );
		clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop);
		time  = (stop.tv_sec-start.tv_sec)*1000;
		time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
		time /= static_cast<double>( rep );
		times[ run ] = time;
	}

	//calculate statistics
	double meantime, mintime, vartime;
	meantime = vartime = 0.0;
	mintime = times[ 0 ];
	for( int run = 0; run < rep; run++ ) {
		if( times[ run ] < mintime ) mintime = times[ run ];
		meantime += times[ run ] / static_cast< double >( rep );
	}
	for( int run = 0; run < rep; run++ ) {
		vartime += ( times[ run ] - meantime ) * ( times[ run ] - meantime ) / static_cast< double >( rep - 1 );
	}

	std::cout << "In-place:" << std::endl;
	std::cout << "Mean =" << checksum( z, matrix->m() ) << std::endl;
	std::cout << "Time = " << meantime << " (average), \t" <<  mintime << " (fastest), \t" << vartime << " (variance) ms. " << std::endl;

	delete [] times;
	delete [] x;
	delete [] z;
	delete matrix;

	return EXIT_SUCCESS;
}

