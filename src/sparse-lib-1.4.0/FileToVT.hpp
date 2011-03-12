
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_FTVT
#define _H_FTVT

#include <cstdlib>
#include <iostream>
#include <vector>

#include "mmio.h"
#include "Triplet.hpp"

/** 
 *  Flag indicating if we support cache-simulated triplets, as defined
 *  in the CACHE_SIM library.
 */
#ifdef _SUPPORT_CS
	#include "CS_Triplet.hpp"
#endif

/** Class responsible for reading in matrix market files and converting them to vector< Triplet > format. */
class FileToVT {

   public:

	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename );
	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename, unsigned long int &m, unsigned long int &n );
	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename, unsigned long int &m, unsigned long int &n, unsigned long int &nnz );
#ifdef _SUPPORT_CS
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename );
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename, unsigned long int &m, unsigned long int &n );
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename, unsigned long int &m, unsigned long int &n, unsigned long int &nnz );
#endif

};

#endif //_H_FTVT

