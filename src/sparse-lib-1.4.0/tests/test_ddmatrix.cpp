#include "DD_MATRIX.hpp"

#include<vector>
#include<iostream>
#include<cstdlib>

//must be defined global, so that it is ensured constant at compiletime (necessary when template is used)
int offsets[] = { -2, 0, 2 };

int main() {
	unsigned long int n = 5;
	unsigned long int m = 5;
	double zero = 0.0;
	double n1[] = { -1, -1, -1 };
	double n2[] = { 1, 1, 1, 1, 1 };
	double n3[] = { -1, -1, -1 };
	double *nonzeroes[] = { n1, n2, n3 };
	DD_MATRIX< double, 3, offsets > test( (double **)nonzeroes, static_cast< unsigned long int >( 5 ) , static_cast< unsigned long int >( 5 ), zero );
	double* x = new double[ n ];
	for( unsigned long int i=0; i<n; i++ )
		x[i]=1.0;
	std::cout << "Init done, starting MV..." << std::endl;
	double* z = test.mv( x );
	std::cout << "y = [";
	double mean = 0;
	for( unsigned long int i=0; i<( m-1 ); i++ ) {
		mean += z[ i ] / m;
		std::cout << z[ i ] << ",";
	}
	mean += z[ m-1 ] / m;
	std::cout << z[ m-1 ] << "]" << std::endl;
	std::cout << "Mean: " << mean << std::endl;
	delete [] z;
	delete [] x;
}

