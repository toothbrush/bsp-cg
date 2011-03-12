
#include "BICRS.hpp"

int main( int argc, char **argv ) {
	
	int row[] = {0,1,2,3,1,1,2,0};
	int col[] = {0,1,2,3,3,2,3,1};
	double vals[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5};
	int m = 4;
	int n = 4;
	int nz = 8;
	BICRS< double > bicrs( row, col, vals, m, n, nz, 0.0 );
	double x[n];
	for( int i=0; i<n; i++ ) x[ i ] = 1.0;
	double *y;
	y = bicrs.mv( x );

	std::cout << "y-vector: " << std::endl;
	std::cout << "( " << y[ 0 ];
	for( int i=1; i<m; i++ ) std::cout << " , " << y[ i ];
	std::cout << " )" << std::endl;

	delete [] y;
}

