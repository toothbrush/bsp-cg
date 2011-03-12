
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#include "util.hpp"

int parseULIVector( std::string filename, std::vector< std::vector< unsigned long int >* > &vector ) {
	std::ifstream in( filename.c_str(), std::ios::in );
	unsigned long int temp;
	unsigned long int swch = 0;

	if( !in ) {
		std::cerr << "Error reading '" << filename << "'!" << std::endl;
		return EXIT_FAILURE;
	}

	while( true ) {
		in >> temp;
		if( !in ) break;
		vector[swch]->push_back( temp );
		swch = (swch+1)%vector.size();
	}
	
	return EXIT_SUCCESS;
}

