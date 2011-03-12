
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include "SparseMatrix.hpp"
#include "Triplet.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>
#include<iostream>

//#define _DEBUG

#ifndef _H_ZZ_CRS
#define _H_ZZ_CRS

/**
 *  The zig-zag compressed row storage sparse matrix data structure.
 */
template< typename T >
class ZZ_CRS: public SparseMatrix< T, long int > {

  private:

  protected:

	/** Array containing the actual this->nnz non-zeros. */
	T* ds;

	/** Array containing the column jumps. */
	long int* col_ind;

	/** Array containing the row jumps. */
	long int* row_start;

	/** Comparison function used for sorting input data. */
	static bool compareTripletsLTR( const Triplet< T >* one, const Triplet< T >* two ) {
		if( one->i() < two->i() )
			return true;
		if ( one->i() > two->i() )
			return false;
		return one->j() < two->j();
	}

	/** Comparison function used for sorting input data. */
	static bool compareTripletsRTL( const Triplet< T >* one, const Triplet< T >* two ) {
		if( one->i() < two->i() )
			return true;
		if ( one->i() > two->i() )
			return false;
		return one->j() > two->j();
	}

  public:

	/** Base constructor. */
	ZZ_CRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	ZZ_CRS( std::string file, T zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid CRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns to be stored.
	 *  @param zero Which element is considered to be the zero element.
	 */
	ZZ_CRS( const long int number_of_nonzeros, const long int number_of_rows, const long int number_of_cols, T zero ):
		SparseMatrix< T, unsigned long int >( number_of_nonzeros, number_of_rows, number_of_cols, zero ) {
		ds = new T[ this->nnz ];
		col_ind = new long int[ this->nnz ];
		row_start = new long int[ this->nor + 1 ];
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	ZZ_CRS( ZZ_CRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->noc = toCopy.noc;
		this->nor = toCopy.nor;
		ds = new T[ this->nnz ];
		col_ind = new long int[ this->nnz ];
		row_start = new long int[ this->nor ];
		for( long int i=0; i<this->nnz; i = i + 1 ) {
			ds[ i ] = toCopy.ds[ i ];
			col_ind[ i ]= toCopy.col_ind[ i ];
		}
		for( long int i=0; i<=this->nor; i++ )
			row_start[ i ] = toCopy.row_start[ i ];
	}

	/**
	 *  Constructor which transforms a collection of input triplets to CRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero Which element is considered zero.
	 */
	ZZ_CRS( std::vector< Triplet< T > >& input, const long int m, const long int n, const T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const long int m, const long int n, const T zero ) {
		this->zero_element = zero;
		this->nnz = input.size();
		this->nor = m;
		this->noc = n;

		std::vector< std::vector< Triplet< T >* > > ds( this->nor );
	
		//move input there
		typename std::vector< Triplet< T > >::iterator in_it = input.begin();
		for( ; in_it != input.end(); in_it++ ) {
			Triplet< T >* cur = &(*in_it);
			const long int currow = cur->i();
			const T value = cur->value;
			if( value == this->zero_element ) { this->nnz--; continue; }
			ds.at( currow ).push_back( cur );
		}

		//allocate arrays
		this->ds = new T[ this->nnz ];
		col_ind = new long int[ this->nnz ];
		row_start = new long int[ this->nor + 1 ];
		
		//make ZZ-CRS
		long int index = 0;
		bool LTR       = true;
		for( long int currow = 0; currow < this->nor; currow++ ) {
			row_start[ currow ] = index;
			if( ds.at( currow ).size() == 0 ) continue;
			if( LTR )
				sort( ds.at( currow ).begin(), ds.at( currow ).end(), compareTripletsLTR );
			else
				sort( ds.at( currow ).begin(), ds.at( currow ).end(), compareTripletsRTL );
			LTR = !LTR;
			typename std::vector< Triplet< T >* >::iterator row_it = ds.at( currow ).begin();
			for( ; row_it!=ds.at( currow ).end(); row_it++ ) {
				const Triplet< T > cur = *(*row_it);
				this->ds[ index ] = cur.value;
				col_ind[ index ] = cur.j();
				index++;
			}
		}
		row_start[ this->nor ] = this->nor;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = static_cast< unsigned long int >( this->row_start[ 0 ] );
		col = static_cast< unsigned long int >( this->col_ind[ 0 ] );
	}

	/**
	 *  In-place z=xA multiplication algorithm.
	 *  TODO optimise!
	 *  @param x The vector x supplied for left-multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zxa( T* x, T* z ) {
		long int index=0, row=0;
		for( ; row < this->nor; row++ )
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				z[ col_ind[ index ] ] += ds[ index ] * x[ row ];
	}

	/**
	 *  In-place z=Ax multiplication algorithm.
	 *  TODO optimise!
	 *  @param x The vector x supplied for multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zax( T* x, T* z ) {
		long int index=0,row=0;
		for( ; row < this->nor; row++ )
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				z[ row ] += ds[ index ] * x[ col_ind[ index ] ];
	}

	/** Base deconstructor. */
	~ZZ_CRS() {
		delete [] ds;
		delete [] col_ind;
		delete [] row_start;
	}

};

#endif

