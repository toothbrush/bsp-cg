#include "HTS.hpp"
#include "FileToVT.hpp"

#include<vector>
#include<sstream>
#include<iostream>
#include<algorithm>

unsigned long int randuli( const unsigned long int m ) {
	return static_cast< unsigned long int >(  (static_cast< double >( rand() ) / static_cast< double >( RAND_MAX ) ) * static_cast< double >( m ) );
}

double randd() {
	return static_cast< double >( rand() ) / static_cast< double >( RAND_MAX );
}

inline
unsigned long int string_to_ulint( std::string s ) {
        std::istringstream iss(s);
        unsigned long int n;
        iss >> n;
	if ( iss.fail() ) {
		std::cerr << "Warning: String to int conversion failure!" << std::endl;
	}
        return n;
}

int main( int argc, char** argv ) {

	if( argc != 3 && argc != 4  && argc != 2 ) {
		std::cout << "Usage: " << argv[0] << " <m> <n> <nnz>, OR" << std::endl;
		std::cout << "       " << argv[0] << " <filename> (which reads a matrix from file (matrix market format)" << std::endl << std::endl;
		std::cout << "Writes the (generated) matrix in Triplet format to \"<filename>.hilbert.trp\". The order of triplets is that of the Hilbert ordering scheme." << std::endl;
		exit( !EXIT_SUCCESS );
	}

	std::vector< Triplet< double > > naive;
	unsigned long int m;
	unsigned long int n;
	double zero = 0.0;
	bool justload=false;

	if( argc == 3 )
		justload = true;

	if( argc == 4 ) {
		m = string_to_ulint( std::string( argv[1] ) );
		n = string_to_ulint( std::string( argv[2] ) );
		unsigned long int nnz = string_to_ulint( std::string( argv[3] ) );;
		zero = 0.0;

		srand( 543213 );
		for( unsigned long int i=0; i<nnz; i++ ) { //note that double indices can occur, but only one will fall through in actual zax calculation.
			const unsigned long int _i = randuli( m );
			const unsigned long int _j = randuli( n );
			const double temp = randd();
			if( temp != zero )
				naive.push_back( Triplet< double >( _i == m ? m-1 : _i , _j == n ? n-1 : _j, randd() ) );
		}
	} else if( argc== 2 || argc==3 ) {
		naive = FileToVT::parse( std::string( argv[1] ), m, n );
		zero = static_cast< double >( 0 );
	}

#ifdef _DEBUG
	for( unsigned long int i = 0; i<naive.size(); i++ ) {
		const Triplet< double > t = naive[ i ];
		std::cout << t.i() << " " << t.j() << " " << t.value << std::endl;
	}
#endif
		
	if( justload ) return EXIT_SUCCESS;

	HTS< double > hts( naive, m, n, zero );

	hts.saveBinary( std::string( argv[1] ) + ".hilbert.trp" );

	return EXIT_SUCCESS;

}

