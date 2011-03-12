
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include <vector>
#include <assert.h>
#include "Triplet.hpp"
#include "SparseMatrix.hpp"

#ifndef _H_TS
#define _H_TS

/** The triplet scheme; a storage scheme for sparse matrices using triplets. */
template< typename T >
class TS: public SparseMatrix< T, unsigned long int > {

   private:
	
	/** Convenience typedef */
	typedef unsigned long int ULI;

   protected:

	/** The row indices of the nonzeros. */
	ULI* i;
	
	/** The column indices of the nonzeros. */
	ULI* j;

	/** The values of the nonzeros. */
	T* ds;

   public:

	/** Base constructor. */
	TS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	TS( std::string file, T zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/** Base constructor.
	 *  @param input std::vector of triplets to be stored in this scheme.
	 *  @param m total number of rows.
	 *  @param n total number of columns.
	 *  @param zero what is to be considered the zero element.
	 */
	TS( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		ULI offset = 0;

		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();
		ds = new T[ this->nnz ];
		i = new ULI[ this->nnz ];
		j = new ULI[ this->nnz ];
		for( ULI r=0; r<this->nnz; r++ )
			if( input[ r ].value != this->zero_element ) {
				//ds[ r ] = input[ r ];
				ds[ r - offset ] = input[ r ].value;
				i[ r - offset ] = input[ r ].i();
				j[ r - offset ] = input[ r ].j();
			} else {
				offset++;
			}
		this->nnz -= offset;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = this->i[ 0 ];
		col = this->j[ 0 ];
	}

	/**
	 *  In-place z=xA calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of xA.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
        virtual void zxa( T* x, T* z ) {
                for( ULI r=0; r<this->nnz; r++ ) {
                        assert( i[ r ] >= 0  );
                        assert( i[ r ] < this->nor );
                        assert( j[ r ] >= 0  );
                        assert( j[ r ] < this->noc );
                        z[ j[ r ] ] += ds[ r ] * x[ i[ r ] ];
                }
        }

	/**
	 *  In-place z=Ax calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of Ax.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
        virtual void zax( T* x, T* z ) {
                for( ULI r=0; r<this->nnz; r++ ) {
                        assert( i[ r ] >= 0  );
                        assert( i[ r ] < this->nor );
                        assert( j[ r ] >= 0  );
                        assert( j[ r ] < this->noc );
                        z[ i[ r ] ] += ds[ r ] * x[ j[ r ] ];
                }
        }

	/** Base destructor. */
	~TS() {
		delete [] i;
		delete [] j;
		delete [] ds;
	}

};

#endif

