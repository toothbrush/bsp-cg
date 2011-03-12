
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2010.
 * Released under GPL, see the file LICENSE.txt
 */

#include <vector>
#include "SparseMatrix.hpp"
#include "HilbertTriplet.hpp"
#include "Triplet.hpp"
#include "BICRS.hpp"

#ifdef _DEBUG
	#include <iostream>
#endif

/** The Hilbert triplet scheme. In effect similar to the Hilbert Triplet scheme (HTS),
    but uses BICRS to store the nonzeroes. */
template< typename T >
class Hilbert: public SparseMatrix< T, unsigned long int > {

   private:
	
	/** Convenience typedef */
	typedef unsigned long int ULI;

   protected:

	/** Minimum number of expansions */
	ULI minexp;

	/** Vector storing the non-zeros and their locations. */
	std::vector< HilbertTriplet< T > > ds;

	/** Actual data structure. */
	BICRS< T > *ads;

   public:

	/** Base deconstructor. */
	virtual ~Hilbert() {
		if( ads != NULL )
			delete ads;
	}

	/** Base constructor. */
	Hilbert() {
		ads = NULL;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	Hilbert( std::string file, T zero = 0 ) {
		ads = NULL;
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
	Hilbert( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
		ads = NULL;
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = 0;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();
		for( ULI i=0; i<this->nnz; i++ ) {
			ds.push_back( HilbertTriplet< T >( input[ i ].i(), input[ i ].j(), input[ i ].value ) );
			ds[ i ].calculateHilbertCoordinate();
		}
		HilbertTripletCompare< T > compare;
		std::sort( ds.begin(), ds.end(), compare );
#ifdef _DEBUG
	        for( ULI i=0; i<this->nnz; i++ )
			std::cout << ds[i].getMostSignificantHilbertBits() << " " << ds[i].getLeastSignificantHilbertBits() << std::endl;
#endif
		//create temporary structure
		std::vector< Triplet< T > > tds;
		typename std::vector< HilbertTriplet< T > >::iterator it = ds.begin();
		for( ; it!=ds.end(); ++it )
			tds.push_back( Triplet< T >( it->i(), it->j(), it->value ) );
		//create actual structure
		ads = new BICRS< T >( tds, m, n, zero );
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = ds[ 0 ].i();
		col = ds[ 0 ].j();
	}

	/**
	 * Calculates z=xA.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zxa( T* x, T* z ) {
		ads->zxa( x, z );
	}

	/**
	 * Calculates z=Ax.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zax( T* x, T* z ) {
		ads->zax( x, z );
	}

	/** Saves the current Hilbert structure in binary triplet form.
	 *  @param fn Filename to save to.
	 *  @see HilbertTriplet< T >::save
`	 */
	void saveBinary( const std::string fn ) {
		HilbertTriplet< T >::save( fn, ds, this->nor, this->noc );
	}

	/** Loads from binary triplets, assumes Hilbert ordering already done */
	void loadBinary( const std::string fn ) {
		std::cerr << "Warning: assuming binary file was saved by a HTS scheme, i.e., that it is pre-ordered using Hilbert coordinates." << std::endl;
		std::vector< Triplet< T > > tds = Triplet< T >::load( fn, this->nor, this->noc );
		ads = new BICRS< T >( tds, this->nor, this->noc, 0 );
	}

};

