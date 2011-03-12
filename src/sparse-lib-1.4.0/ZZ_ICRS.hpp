
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include "SparseMatrix.hpp"
#include "Triplet.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>

//#define _DEBUG

#ifndef _H_ZZ_ICRS
#define _H_ZZ_ICRS

//#define _UNSIGNED

//#ifdef _DEBUG
#include<iostream>
//#endif

/**
 *  The zig-zag incremental compressed row storage sparse matrix data structure.
 */
template< typename T >
#ifdef _UNSIGNED
class ZZ_ICRS: public SparseMatrix< T, unsigned long int > {
#else
class ZZ_ICRS: public SparseMatrix< T, long int > {
#endif

  private:

#ifdef _UNSIGNED
	/** Conveniece typedef */
	typedef unsigned long int ULI;
#else
	/**Convenience typedef */
	typedef long int ULI;
#endif

  protected:

	/** Array containing the actual this->nnz non-zeros. */
	T* ds;

	/** Array containing the column jumps. */
	ULI* c_ind;

	/** Array containing the row jumps. */
	ULI* r_ind;

	/** Comparison function used for sorting input data. */
	static bool compareTriplets( const Triplet< T >& one, const Triplet< T >& two ) {
		if( one.i() < two.i() )
			return true;
		if ( one.i() > two.i() )
			return false;
		return one.j() < two.j();
	}

  public:

	/** Base constructor. */
	ZZ_ICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	ZZ_ICRS( std::string file, T zero = 0 ) {
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
	ZZ_ICRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ):
		SparseMatrix< T, ULI >( number_of_nonzeros, number_of_rows, number_of_cols, zero ) {
		ds = new T[ this->nnz ];
		c_ind = new ULI[ this->nnz ];
		r_ind = new ULI[ this->nnz ];
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	ZZ_ICRS( ZZ_ICRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->noc = toCopy.noc;
		this->nor = toCopy.nor;
		ds = new T[ this->nnz ];
		c_ind = new ULI[ this->nnz ];
		r_ind = new ULI[ this->nnz ];
		for( ULI i=0; i<this->nnz; i = i + 1 ) {
			ds[ i ] = toCopy.ds[ i ];
			c_ind[ i ]= toCopy.c_ind[ i ];
			r_ind[ i ] = toCopy.r_ind[ i ];
		}
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
	ZZ_ICRS( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = zero;

		this->nor = m;
		this->noc = n;
	
		typename std::vector< Triplet< T > >::iterator in_it;

		//WARNING: this->noc*this->nor typically overflows on 32 bits!
		//         immediately start recording differences
		//         instead of trying to directly calculate
		//         the index as i*this->noc+j.
	
		sort( input.begin(), input.end(), compareTriplets );

		//filter out zeros and count the number of row jumps
		std::vector< ULI > r_ind_temp;
		typename std::vector< Triplet< T > >::iterator it = input.begin();
		unsigned long int prev_row = (*it).i();
		r_ind_temp.push_back( prev_row );
		it++;
		for( ; it!=input.end(); it++ ) {
			if( (*it).i() > prev_row ) {
				r_ind_temp.push_back( (*it).i() - prev_row );	
				prev_row = (*it).i();
			}
		}
		this->nnz = input.size();

		//allocate arrays
		ds = new T[ this->nnz ];
		c_ind = new ULI[ this->nnz + 1 ];
		r_ind = new ULI[ r_ind_temp.size() + 1 ];

		//set last entry
		c_ind[ this->nnz ] = this->noc;
		//r_ind does not have to be set; altough the last element is read, it is actually never used.

		//copy row-jump vector
		for( ULI i=0; i<static_cast< ULI >( r_ind_temp.size() ); i++ )
			r_ind[ i ] = r_ind_temp[ i ];

		//make ICRS
		ULI prev_col = 0;
		prev_row = 0;

		for( ULI i=0; i<this->nnz; ++i ) {
			const Triplet< T > cur = input[ i ];
			const unsigned long int currow = cur.i();
			const unsigned long int curcol = cur.j();
			if( currow == prev_row ) {
				c_ind[ i ] = curcol - prev_col;
				ds[ i ] = cur.value;
				prev_col = curcol;
			} else {
				prev_row = currow;
				ULI stop = i;
				ULI j = i;
				
				while( j < this->nnz && input[ j ].i() == static_cast< unsigned long int >( currow ) )
					j++; //go to next row
				//ULI next = j; //where the next row starts, OR next = j = this->nnz // unused
				j--; // j now points at last element at new row

				//process the first new element
				ds[ i ] = input[ j ].value;
				c_ind[ i ] = this->noc + ( input[ j ].j() - prev_col );
				prev_col = input[ j ].j();

				i++; //go to next indices
				j--;

				for( ; j>=stop; j-- ) { //reversively add the rest to data structure
					ds[ i ] = input[ j ].value;
					c_ind[ i ] = prev_col - input[ j ].j();
					prev_col = input[ j ].j();
					i++; //increment target index
				}

				if( i >= this->nnz ) break;

				//at this point, the reverse row is completely added
				//and i should be at the next row index.
#ifdef _UNSIGNED
				c_ind[ i ] = this->noc + input[ i ].j();
#else
				c_ind[ i ] = this->noc + ( prev_col - input[ i ].j() );
#endif

				ds[ i ] = input[ i ].value;
				prev_col = input[ i ].j();
				prev_row = input[ i ].i();
			}
		}
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = static_cast< unsigned long int >( this->r_ind[ 0 ] );
		col = static_cast< unsigned long int >( this->c_ind[ 0 ] );
	}

	/**
	 *  In-place z=xA multiplication algorithm.
	 *
	 *  @param x The vector x supplied for left-multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zxa( T*__restrict__ pDataX, T*__restrict__ pDataZ ) {

		//get some pointers
		T   *__restrict__ pDataA	= ds;
		ULI *__restrict__ pIncRow   	= r_ind;
		ULI *__restrict__ pIncCol   	= c_ind;
		const T * const pDataZbeg 	= pDataZ;
		const T * const pDataZend 	= pDataZ + this->noc;
		const T *__restrict__ const pDataAend = ds + this->nnz;
		//const T * const pDataXbeg = pDataX;
		//const T * const pDataXend = pDataX + this->nor;

		//should be optimised out by compiler when NDEBUG is defined 
		ULI rowshadow  = 0;

		//go to first column
		pDataZ += *pIncCol;
		pIncCol++;

		//walk through all matrix entries
		while( pDataA < pDataAend ) {

			//jump to correct row
			pDataX += *pIncRow;
			rowshadow += *pIncRow;

			//walk left-to-right
			while( pDataZ < pDataZend ) {

				*pDataZ  += *pDataA  * *pDataX;
				 pDataA++;
				 pDataZ  += *pIncCol          ;
				 pIncCol++;

			}

			//check if there is no more data
			if( pDataA >= pDataAend ) break;

			//jump to correct row and column
			pDataZ  -=  this->noc  ;
			pIncRow++;
			pDataX  += *pIncRow    ;
			rowshadow+=*pIncRow;
			assert( rowshadow >= 0  );
			assert( rowshadow < this->nor );

			//walk right-to-left
			while( pDataZbeg <= pDataZ ) {
				
				*pDataZ += *pDataA  * *pDataX;
				 pDataA++;
				 pDataZ -= *pIncCol          ;
 				 pIncCol++;

			}

			//jump to crrect column
			pDataZ += this->noc;

			//set row increment pointer to crrect location
			pIncRow++;
			
		}

	}


	/**
	 *  In-place z=Ax multiplication algorithm.
	 *
	 *  @param x The vector x supplied for multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zax( T*__restrict__ pDataX, T*__restrict__ pDataZ ) {

		//get some pointers
		T   *__restrict__ pDataA    = ds;
		ULI *__restrict__ pIncRow   = r_ind;
		ULI *__restrict__ pIncCol   = c_ind;
		const T * const pDataXbeg = pDataX;
		const T * const pDataXend = pDataX + this->noc;
		const T *__restrict__ const pDataAend = ds + this->nnz;

		//should be optimised out by compiler when NDEBUG is defined 
		ULI rowshadow  = 0;

		//go to first column
		pDataX += *pIncCol;
		pIncCol++;

		//walk through all matrix entries
		while( pDataA < pDataAend ) {

			//jump to correct row
			pDataZ += *pIncRow;
			rowshadow += *pIncRow;

			//walk left-to-right
			while( pDataX < pDataXend ) {

				*pDataZ  += *pDataA  * *pDataX;
				 pDataA++;
				 pDataX  += *pIncCol          ;
				 pIncCol++;

			}

			//check if there is no more data
			if( pDataA >= pDataAend ) break;

			//jump to correct row and column
			pDataX   -=  this->noc  ;
			pIncRow++;
			pDataZ   += *pIncRow    ;
			rowshadow+=*pIncRow;
			assert( rowshadow >= 0  );
			assert( rowshadow < this->nor );

			//walk right-to-left
			while( pDataXbeg <= pDataX ) {
				
				*pDataZ += *pDataA  * *pDataX;
				 pDataA++;
				 pDataX -= *pIncCol          ;
 				 pIncCol++;

			}

			//jump to crrect column
			pDataX += this->noc;

			//set row increment pointer to crrect location
			pIncRow++;
			
		}

	}

	/** Base deconstructor. */
	~ZZ_ICRS() {
		delete [] ds;
		delete [] c_ind;
		delete [] r_ind;
	}

};

#endif

