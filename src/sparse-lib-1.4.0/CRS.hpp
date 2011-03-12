/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */


#include "Triplet.hpp"
#include "SparseMatrix.hpp"
#include <vector>
#include <assert.h>

//#define _DEBUG

#ifndef _H_CRS
#define _H_CRS

#ifdef _DEBUG
#include<iostream>
#endif

/**
 *  The compressed row storage sparse matrix data structure.
 */
template< typename T >
class CRS: public SparseMatrix< T, unsigned long int > {

  private:

	/** Convenience typedef */
	typedef unsigned long int ULI;
	
  protected:

	/** Array keeping track of individual row starting indices. */
	ULI* row_start;

	/** Array containing the actual nnz non-zeros. */
	T* ds;

	/** Array containing the column indeces corresponding to the elements in ds. */
	ULI* col_ind;

 	/** Sorts 1D columnwise */
        static int compareTriplets( const void * left, const void * right ) {
                const Triplet< T > one = **( (Triplet< T > **)left );
                const Triplet< T > two = **( (Triplet< T > **)right );
                if ( one.j() < two.j() )
                        return -1;
                if ( one.j() > two.j() )
                        return 1;
                return 0;
        }     

	/** 
         *  Helper function which finds a value with a given column index on a given subrange of indices.
         *  @param col_index The given column index.
	 *  @param search_start The start index of the subrange (inclusive).
	 *  @param search_end The end index of the subrange (exlusive).
	 *  @param ret Reference to the variable where the return *index* is stored.
	 *  @return Whether or not a non-zero value should be returned.
	 */
	bool find( const ULI col_index, const ULI search_start, const ULI search_end, ULI &ret ) {
		for( ULI i=search_start; i<search_end; i++ )
			if( col_ind[ i ] == col_index ) {
				ret = i;
				return true;
			}
		return false;
	}

  public:

	/** Base constructor. */
	CRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	CRS( std::string file, T zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid CRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns of the matrix.
	 *  @param zero The element considered to be zero.
	 */
	CRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ) {
		this->nnz = number_of_nonzeros;
		this->nor = number_of_rows;
		this->noc = number_of_cols;
		this->zero_element = zero;
		row_start = new ULI[ this->nor + 1 ];
		ds = new T[ this->nnz ];
		col_ind = new ULI[ this->nnz ];
	}

	/** Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	CRS( CRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->nor = toCopy.nor;
		row_start = new ULI[ this->nor + 1 ];
		ds = new T[ this->nnz ];
		col_ind = new ULI[ this->nnz ];
		for( ULI i=0; i<this->nnz; i++ ) {
			ds[ i ] = toCopy.ds[ i ];
			col_ind[ i ] = toCopy.col_ind[ i ];
		}
		for( ULI i=0; i<this->nor; i++ )
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
	 *  @param zero The element considered to be zero.
	 */
	CRS( std::vector< Triplet< T > > input, unsigned long int m, unsigned long int n, T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, unsigned long int m, unsigned long int n, T zero ) {
		std::cout << "\tLoading in a vector of " << input.size() << " triplets into CRS..." << std::endl;

		this->zero_element = zero;
		//find nnz
		this->nnz = input.size();
	
		this->nor = m;
		this->noc = n;

		//build better datastructure
		std::vector< std::vector< Triplet< T >* > > ds( this->nor );
		
		//move input there
		typename std::vector< Triplet< T > >::iterator in_it = input.begin();
		for( ; in_it != input.end(); ++in_it ) {
			//Triplet< T >* cur = &(*in_it);
			const ULI currow = in_it->i();
			const T value = in_it->value;
			if( value == this->zero_element ) { 
				this->nnz--;
				continue;
			}
			ds.at( currow ).push_back( &(*in_it) );
		}

		//allocate arrays
		row_start = new ULI[ this->nor + 1 ];
		this->ds = new T[ this->nnz ];
		col_ind = new ULI[ this->nnz ];

                //make CRS
                ULI index = 0;
                for( ULI currow = 0; currow < this->nor; currow++ ) {
                        row_start[ currow ] = index;
                        if( ds.at( currow ).size() == 0 ) continue;
                        qsort( &( ds.at( currow )[ 0 ] ), ds.at( currow ).size(), sizeof( Triplet< T >* ), &compareTriplets );
                        typename std::vector< Triplet< T >* >::iterator row_it = ds.at( currow ).begin();
                        for( ; row_it!=ds.at( currow ).end(); row_it++ ) {
                                const Triplet< T > cur = *(*row_it);
                                this->ds[ index ] = cur.value;
                                col_ind[ index ] = cur.j();
                                index++;
                        }
                }
		row_start[ this->nor ] = this->nnz;
		std::cout << "\t" << index << " nonzeroes loaded into CRS structure." << std::endl;
		assert( index == this->nnz );
	}

	/**
	 *  Method which provides random matrix access to the stored sparse matrix.
	 *  @param i Row index.
	 *  @param j Column index.
	 *  @return Matrix valuei at (i,j).
	 */
	T& random_access( ULI i, ULI j ) {
		ULI found_index;
		if ( find( j, row_start[ i ], row_start[ i+1 ], found_index ) ) {
#ifdef _DEBUG
			std::cout << "Searched col_ind between " << row_start[ i ] << " and " << row_start[ i + 1 ] << ", found: " << std::endl;
			std::cout << "Element (" << i << "," << j << ") found on index " << found_index << ", returning " << ds[ found_index ] << std::endl;
#endif
			return ds[ found_index ];
		} else
			return this->zero_element;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = this->row_start[ 0 ];
		col = this->col_ind[ 0 ];
	}

	/** 
	 *  In-place z=xA function.
	 */
        virtual void zxa( T*__restrict__ x, T*__restrict__ z ) {
		T*__restrict__ ds_p = ds;
		ULI*__restrict__ col_ind_p = col_ind;
		ULI row, index;
		for( row = 0; row < this->nor; row++, x++ ) {
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				z[ *col_ind_p++ ] += *ds_p++ * *x;
		}
        }

	/** 
	 *  In-place z=Ax function.
	 *  
	 *  @param x The x vector to multiply current matrix with.
	 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zax( T*__restrict__ x, T*__restrict__ z ) {
		T*__restrict__ ds_p = ds;
		ULI*__restrict__ col_ind_p = col_ind;
		ULI row, index;
		for( row = 0; row < this->nor; row++, z++ ) {
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				*z += *ds_p++ * x[ *col_ind_p++ ];
		}
        }

	/** Returns pointer to the row_start vector. */
	ULI* rowJump() { return row_start; }

	/** Returns pointer to the column index vector. */
	ULI* columnIndices() { return col_ind; }

	/** Returns pointer to the matrix nonzeros vector. */
	T* values() { return ds; } 

	/** Base deconstructor. */
	virtual ~CRS() {
		delete [] row_start;
		delete [] ds;
		delete [] col_ind;
	}

};

#undef _DEBUG
#endif

