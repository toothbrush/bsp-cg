
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2010.
 * Released under GPL, see the file LICENSE.txt
 */

#include <vector>
#include <algorithm>
#include "SparseMatrix.hpp"
#include "HilbertTriplet.hpp"
#include "HilbertTripletCompare.hpp"
#include "Triplet.hpp"

#ifdef _DEBUG
	#include <iostream>
#endif

#ifndef _H_HTS
#define _H_HTS

/** The Hilbert triplet scheme. In effect similar to the triplet scheme (TS),
    but uses Hilbert coordinates to determine the order of storage. */
template< typename T >
class HTS: public SparseMatrix< T, unsigned long int > {

   private:
	
	/** Convenience typedef */
	typedef unsigned long int ULI;

   protected:

	/** Minimum number of expansions */
	ULI minexp;

	/** Vector storing the non-zeros and their locations. */
	std::vector< HilbertTriplet< T > > ds;

	/** HilbertCoordinate comparison function. */
	bool cmp( HilbertTriplet< T > &left, HilbertTriplet< T > &right ) {

		if( left.i() == right.i() && left.j() == right.j() ) return true;

		const std::vector< unsigned long int > h_one = left.getHilbert();
		const std::vector< unsigned long int > h_two = right.getHilbert();

		unsigned long int max = h_one.size();
		bool max_is_one = true;
		if ( h_two.size() < max ) { max = h_two.size(); max_is_one = false; }
		for( unsigned long int i=0; i<max; i++ )
			if( h_one[ i ] != h_two[ i ] )
				return h_one[ i ] < h_two[ i ];
#ifdef _DEBUG		
		std::cout << "expand, ";
#endif
		if ( max_is_one )
			left.morePrecision( this->nor, this->noc );
		else
			right.morePrecision( this->nor, this->noc );

		return cmp( left, right );
	}

	/**
	 *  Binary search for finding a given triplet in a given range.
	 *  @param x triplet to-be found.
	 *  @param left Left bound of the range to search in.
	 *  @param right Right bound of the range to search in.
	 *  @return Index of the triplet searched.
	 */	
	unsigned long int find( HilbertTriplet< T > &x, ULI &left, ULI &right ) {
#ifdef _DEBUG
		std::cout << "left: " << left << ", right: " << right << std::endl;
#endif
		if( ds[left].getHilbert().size() < ds[minexp].getHilbert().size() ) minexp = left;
		if( ds[right].getHilbert().size() < ds[minexp].getHilbert().size() ) minexp = right;

		if ( left == right ) return left;
		if ( left+1 == right ) return right;
		if ( cmp( x, ds[ left ] ) ) return left;
		if ( cmp( ds[ right ], x ) ) return right+1;

		ULI mid = static_cast< unsigned long int >( ( left + right ) / 2 );
#ifdef _DEBUG
		std::cout << "mid: " << mid << std::endl;
#endif
		if ( cmp( x, ds[ mid ] ) )
			return find( x, left, mid );
		else
			return find( x, mid, right );
	}

   public:

	/** Base constructor. */
	HTS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	HTS( std::string file, T zero = 0 ) {
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
	HTS( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
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
		for( ULI i=0; i<this->nnz; i++ )
			z[ ds[i].j() ] += ds[i].value * x[ ds[i].i() ];
	}

	/**
	 * Calculates z=Ax.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zax( T* x, T* z ) {
		for( ULI i=0; i<this->nnz; i++ )
			z[ ds[i].i() ] += ds[i].value * x[ ds[i].j() ];
	}

	/** Saves the current HTS.
	 *  @param fn Filename to save to.
	 *  @see HilbertTriplet< T >::save
`	 */
	void saveBinary( const std::string fn ) {
		HilbertTriplet< T >::save( fn, ds, this->nor, this->noc );
	}

	/** Loads from binary triplets, assumes Hilbert ordering already done */
	void loadBinary( const std::string fn ) {
		std::cerr << "Warning: assuming binary file was saved by HTS scheme, i.e., that it is pre-ordered using Hilbert coordinates." << std::endl;
		ds = Triplet< T >::load( fn, this->nor, this->noc );
	}

};

#endif

