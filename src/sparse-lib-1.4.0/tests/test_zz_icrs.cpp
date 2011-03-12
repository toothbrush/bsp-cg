#include "ZZ_ICRS.hpp"
#include "FileToVT.hpp"

#include<vector>
#include<iostream>

unsigned long int randuli( const unsigned long int m ) {
	return static_cast< unsigned long int >(  (static_cast< double >( rand() ) / static_cast< double >( RAND_MAX ) ) * static_cast< double >( m ) );
}

double randd() {
	return static_cast< double >( rand() ) / static_cast< double >( RAND_MAX );
}

int main( int argc, char ** argv ) {
	unsigned long int m, n, nnz;
	double zero = 0.0;
	std::vector< Triplet< double > > naive;

	if ( argc <= 1 ) {
		m = 5000;
		n = 5000;
		nnz = 15000;

		srand( 543213 );
		for( unsigned long int i=0; i<nnz; i++ ) {
			const unsigned long int _i = randuli( m );
			const unsigned long int _j = randuli( n );
			//naive.push_back( Triplet< double >( i, i, 1.0 ) );
			//naive.push_back( Triplet< double >( _i == m ? m-1 : _i , _j == n ? n-1 : _j, 1.0 ) );
			naive.push_back( Triplet< double >( _i == m ? m-1 : _i , _j == n ? n-1 : _j, randd() ) );
		}
	} else
		naive = FileToVT::parse( std::string( argv[1] ), m, n, nnz );

	Triplet< double >::save( std::string( "temp.trp" ), naive, m, n );

	ZZ_ICRS< double > test( naive, m, n, zero );
	
	double* x = new double[ n ];
	for( unsigned long int i=0; i<n; i++ )
		x[i]=1.0;
	std::cout << "Init done, starting mv..." << std::endl;
	double* z = test.mv( x );
	//std::cout << "[";
	double mean = 0;
	for( unsigned long int i=0; i<( m-1 ); i++ ) {
		mean += z[ i ] / m;
		//std::cout << z[ i ] << ",";
	}
	mean += z[ m-1 ] / m;
	//std::cout << z[ m-1 ] << "]" << std::endl;
	std::cout << "Mean: " << mean << std::endl;
	delete [] z;
	delete [] x;
}

