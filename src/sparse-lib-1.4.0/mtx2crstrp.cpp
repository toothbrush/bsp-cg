
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>

#include "SparseMatrix.hpp"
#include "TS.hpp"

bool compare( const Triplet< double > one, const Triplet< double > two ) {
	if ( one.i() < two.i() )
		return true;
	if ( one.i() > two.i() )
		return false;
	if ( one.j() < two.j() )
		return true;
	return false;
}

int main( int argc, char** argv ) {
	
	if( argc != 3 ) {
		std::cout << "Usage: " << argv[0] << " <input-matrix> <output-trp>" << std::endl << std::endl;
		std::cout << "sorts a matrix market or .trp file into CRS order and writes to output binary tiplets." << std::endl;
		std::cout << "Note: binary triplet format is machine dependent;";
		std::cout << "      take care when using the same binary files on different machine architectures." << std::endl;
		return EXIT_FAILURE;
	}

	std::string file = std::string( argv[1] );
	std::string out  = std::string( argv[2] );

	unsigned long int m, n;
	std::vector< Triplet< double > > matrix = FileToVT::parse( file, m, n );

	std::cout << "Matrix dimensions: " << m << " times " << n << "." << std::endl;
	std::cout << "Sorting..." << std::endl;

	std::sort( matrix.begin(), matrix.end(), compare );

	std::cout << "Saving..." << std::endl;
	Triplet< double >::save( out, &(matrix[0]), m, n, matrix.size() );

	std::cout << "Done!" << std::endl;
	return EXIT_SUCCESS;
}

