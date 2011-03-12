
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include "BlockHilbert.hpp"

#ifndef _H_BISECTIONHILBERT
#define _H_BISECTIONHILBERT

/** The Bisection Hilbert triplet scheme. In effect similar to (HierarchicalBICRS),
    but uses Hilbert coordinates to determine the order of the blocks,
    and a bisection algorithm to construct the individual blocks.
    Wraps around the BisectionHilbert class which already implements this scheme. */
template< typename T >
class BisectionHilbert: public BlockHilbert< T > {

   private:
	
   protected:

   public:

	/** Base deconstructor. */
	virtual ~BisectionHilbert() {}

	/** Base constructor. */
	BisectionHilbert() {
		this->bisection = 1;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @param bisect Whether bisection-based blocking should be used.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	BisectionHilbert( std::string file, T zero = 0 ) {
		this->bisection = 1;
		loadFromFile( file, zero );
	}

	/** 
	 *  Base constructor.
	 *  Warning: the zero parameter is currently NOT USED!
	 *  @param input Raw input of normal triplets.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero What elements is considered to-be zero.
	 */
	BisectionHilbert( std::vector< Triplet< T > >& input, unsigned long int m, unsigned long int n, T zero ) {
		this->bisection = 1;
		load( input, m, n, zero );
	}

};

#endif

