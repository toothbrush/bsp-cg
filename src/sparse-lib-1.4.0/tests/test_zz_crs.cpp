#include "ZZ_CRS.hpp"

#include<vector>
#include<iostream>
#include<cstdlib>

unsigned long int randuli( const unsigned long int m ) {
	return static_cast< unsigned long int >(  (static_cast< double >( rand() ) / static_cast< double >( RAND_MAX ) ) * static_cast< double >( m ) );
}

double randd() {
	return static_cast< double >( rand() ) / static_cast< double >( RAND_MAX );
}

int main() {
	unsigned long int m = 5000;
	unsigned long int n = 5000;
	unsigned long int nnz = 15000;
	double zero = 0.0;

	std::vector< Triplet< double > > naive;

	srand( 543213 );
	for( unsigned long int i=0; i<nnz; i++ ) {
		const unsigned long int _i = randuli( m );
		const unsigned long int _j = randuli( n );
		//naive.push_back( Triplet< double >( i, i, 1.0 ) );
		//naive.push_back( Triplet< double >( _i == m ? m-1 : _i , _j == n ? n-1 : _j, 1.0 ) );
		naive.push_back( Triplet< double >( _i == m ? m-1 : _i , _j == n ? n-1 : _j, randd() ) );
	}

	ZZ_CRS< double > test( naive, m, n, zero );

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

