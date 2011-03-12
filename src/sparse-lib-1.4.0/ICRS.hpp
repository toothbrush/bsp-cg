
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include "Triplet.hpp"
#include "SparseMatrix.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>

//#define _DEBUG

#ifndef _H_CS_CRS
#define _H_CS_CRS

#ifdef _DEBUG
#include<iostream>
#endif

/**
 *  The *incremental* compressed row storage sparse matrix data structure.
 */
template< typename T >
class ICRS: public SparseMatrix< T, unsigned long int > {

  private:

	/** Convenience typedef */
	typedef unsigned long int ULI;

  protected:

	/** Array containing the actual nnz non-zeros. */
	T* ds;

	/** Array containing the column jumps. */
	ULI* c_ind;

	/** Array containing the row jumps. */
	ULI* r_ind;

	/** Comparison function used for sorting input data. */
	static int compareTriplets( const void * left, const void * right ) {
                const Triplet< T > one = *( (Triplet< T > *)left );
                const Triplet< T > two = *( (Triplet< T > *)right );
                if( one.i() < two.i() )
                        return -1;
                if ( one.i() > two.i() )
                        return 1;
                if ( one.j() < two.j() )
                        return -1;
                if ( one.j() > two.j() )
                        return 1;
                return 0;
        }

  public:

	/** Base constructor. */
	ICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	ICRS( std::string file, T zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid ICRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows of the matrix to be stored.
	 *  @param number_of_cols The number of columns of the matrix to be stored.
	 *  @param zero Which element is considered to be the zero element.
	 */
	ICRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ):
		SparseMatrix< T, ULI >( number_of_nonzeros, number_of_cols, number_of_rows, zero ) {
		ds = new T[ this->nnz ];
		c_ind = new ULI[ this->nnz ];
		r_ind = new ULI[ this->nnz ];
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy Reference to the CRS datastructure to copy.
	 */
	ICRS( ICRS< T >& toCopy ) {
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
	 *  pairs of indices. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *
	 *  @param input The input collection of triplets (i,j,val).
	 *  @param m The number of rows of the input matrix.
	 *  @param n The number of columns of the input matrix.
	 *  @param zero Which element is considered zero.
	 */
	ICRS( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = zero;

		this->nor = m;
		this->noc = n;
	
		typename std::vector< Triplet< T > >::iterator in_it;

		//WARNING: noc*nor typically overflows on 32 bits!
		//         immediately start recording differences
		//         instead of trying to directly calculate
		//         the index as i*noc+j.

		//Complexity compiler package-dependent. Probably O(nnz log(nnz)) average, O(nnz^2) worst-case.
		//for standard C++ sort:
		//sort( input.begin(), input.end(), compareTriplets );
		//for old school C quicksort:
		qsort( &input[ 0 ], input.size(), sizeof( Triplet< T > ), &compareTriplets );

		//filtering out zeros is skipped for now.
		
		//Count the number of row jumps
		std::vector< ULI > r_ind_temp;
		typename std::vector< Triplet< T > >::iterator it = input.begin();
		ULI prev_row = (*it).i();
		r_ind_temp.push_back( prev_row );
		it++;

		//O(nnz log(nnz)); push_back on a vector uses amortised log(n) array growing algorithm on inserts
		for( ; it!=input.end(); it++ ) {
			assert( (*it).i() >= prev_row );
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
		//O(m) worst case
		for( ULI i=0; i<r_ind_temp.size(); i++ ) {
			r_ind[ i ] = r_ind_temp[ i ];
		}

		//make ICRS
		ULI prev_col = 0;
		prev_row = r_ind[ 0 ];

		//O(nnz)
		unsigned long int check_jumps = 0;
		unsigned long int check_row   = r_ind_temp[0];
		for( ULI i=0; i<this->nnz; ++i ) {
			const Triplet< T > cur = input[ i ];
			const ULI currow = cur.i();
			const ULI curcol = cur.j();
			if( currow == prev_row ) {
				c_ind[ i ] = curcol - prev_col;
				assert( currow == check_row );
			} else {
				assert( currow > prev_row );
				c_ind[ i ] = this->noc + ( curcol - prev_col );
				prev_row = currow;
				check_jumps++;
				check_row += r_ind_temp[ check_jumps ];
				assert( currow == check_row );
			}
			ds[ i ] = cur.value;
			prev_col = curcol;

#ifdef _DEBUG
			std::cout << currow << "," << curcol << "(" << cur.value << ") maps to " << c_ind[ i ] << std::endl;
#endif
		}

		//assert row jumps is equal to r_ind_temp's size
		assert( check_jumps == r_ind_temp.size()-1 );
		
		//clear temporary r_ind vector
		r_ind_temp.clear();
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = this->r_ind[ 0 ];
		col = this->c_ind[ 0 ];
	}

	/**
	 *  In-place z=xA function. Adapted from the master thesis of Joris Koster, Utrecht University.
	 *
	 *  @param pDataX Pointer to array x to multiply by the current matrix (Ax).
	 *  @param pDataZ Pointer to result array. Must be pre-allocated and its elements set to zero for correct results.
	 */
        virtual void zxa( T*__restrict__ pDataX, T*__restrict__ pDataZ ) {
		         T *__restrict__ pDataA    = ds;
		const    T *__restrict__ pDataAend = ds + this->nnz;
		//const    T *__restrict__ const pDataXend = pDataX + this->nor;
		const T * const pDataZend = pDataZ + this->noc;
		//const    T *pDataZend = z + nor; //unused

		//register ULI *__restrict__ pIncRow   = r_ind; //use of register keyword is discouraged
		//register ULI *__restrict__ pIncCol   = c_ind; //use of register keyword is discouraged
		ULI *__restrict__ pIncRow   = r_ind;
		ULI *__restrict__ pIncCol   = c_ind;

		//go to first column
		pDataZ += *pIncCol++;
		while( pDataA < pDataAend ) {
			//jump to correct row
			pDataX += *pIncRow;
			while( pDataZ < pDataZend ) {
				*pDataZ += *pDataA * *pDataX;
				pDataA++;
				pDataZ  += *pIncCol;
				pIncCol++;
			}
			pDataZ -= this->noc;
			pIncRow++;
		}
        }

	/**
	 *  In-place z=Ax function. Adapted from the master thesis of Joris Koster, Utrecht University.
	 *
	 *  @param pDataX Pointer to array x to multiply by the current matrix (Ax).
	 *  @param pDataZ Pointer to result array. Must be pre-allocated and its elements set to zero for correct results.
	 */
        virtual void zax( T*__restrict__ pDataX, T*__restrict__ pDataZ ) {
		         T *__restrict__ pDataA    = ds;
		const    T *__restrict__ pDataAend = ds + this->nnz;
		const T * const pDataXend = pDataX + this->noc;
		//const T * const pDataZend = pDataZ + nor; //unused

		//register ULI *__restrict__ pIncRow   = r_ind; //use of register keyword is discouraged
		//register ULI *__restrict__ pIncCol   = c_ind; //use of register keyword is discouraged
		ULI *__restrict__ pIncRow   = r_ind;
		ULI *__restrict__ pIncCol   = c_ind;

		//T *__restrict__ const x = pDataX;
		//T *__restrict__ const z = pDataZ;

		//go to first column
		pDataX += *pIncCol++;
		while( pDataA < pDataAend ) {
			//jump to correct row
			pDataZ += *pIncRow;
			while( pDataX < pDataXend ) {
				*pDataZ += *pDataA * *pDataX;
				pDataA++;
				pDataX  += *pIncCol;
				pIncCol++;
			}
			pDataX -= this->noc;
			pIncRow++;
		}
        }

	/** Base deconstructor. */
	~ICRS() {
		delete [] ds;
		delete [] c_ind;
		delete [] r_ind;
	}

};

#endif

