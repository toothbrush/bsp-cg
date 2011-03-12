#include "TS.hpp"

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
	unsigned long int m = 0;
	unsigned long int n = 0;
	double zero = 0.0;

	std::vector< Triplet< double > > naive = Triplet< double >::load( "test.trp", m, n );

	TS< double > test( naive, m, n, zero );
	
	double* x = new double[ n ];
	for( unsigned long int i=0; i<n; i++ )
		x[i]=1.0;
	std::cout << "Init done, starting MV..." << std::endl;
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

