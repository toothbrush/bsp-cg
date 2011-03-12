
#include "SparseMatrix.hpp"
#include <assert.h>
#include <iostream>

/** Set when in debug mode */
//#define _DEBUG

#ifndef _H_BICRS
#define _H_BICRS

/**
 *  Bi-directional Incremental Compressed Row Storage scheme.
 *  Supports jumping back and forward within columns.
 *  Supports jumping back and forward between rows.
 *  Main storage direction in column-wise.
 *  Storage requirements are 2nz plus the number of row jumps required.
 *  Many row jumps are disadvantageous to storage as well as speed.
 *  @param _t_value The type of the nonzeros in the matrix.
 *
 *  Warning: this class uses assertions! For optimal performance,
 *           define the NDEBUG flag (e.g., pass -DNDEBUG as a compiler
 *	     flag).
 */
template< typename _t_value >
class BICRS: public SparseMatrix< _t_value, unsigned long int > {

     protected:
	/** Stores the row jumps; size is _at maximum_ the number of nonzeros. */
	int* r_inc;

	/** Stores the column jumps; size is exactly the number of nonzeros. */
	int* c_inc;

	/** Stores the values of the individual nonzeros. Size is exactly the number of nonzeros. */
	_t_value* vals;

	/** Caches n times two. */
	int ntt;
	
     public:

	/** Base deconstructor. */
	virtual ~BICRS() {
		delete [] r_inc;
		delete [] c_inc;
		delete [] vals;
	}

	/** Base constructor. */
	BICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	BICRS( std::string file, _t_value zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor. Stores triplets in exactly the same order as passed
	 *  to this constructor.
	 *  @param row The row numbers of the individual nonzeros.
	 *  @param col The column numbers of the individual nonzeros.
	 *  @param val The values of the nonzeros.
	 *  @param m Number of matrix rows.
	 *  @param n Number of matrix columns.
	 *  @param nz Number of nonzeros.
	 *  @param zero Which value is to be regarded zero here.
	 */
	BICRS( int* row, int* col, _t_value* val, int m, int n, int nz, _t_value zero ) {
		load( row, col, val, m, n, nz, zero );
	}

	/** Base constructor.
	 *  @see SparseMatrix::SparseMatrix( input, m, n, zero )
	 */
	BICRS( std::vector< Triplet< _t_value > >& input, unsigned long int m, unsigned long int n, _t_value zero ) {
		load( input, m, n, zero );
	}

	/**
	 *  This function will rewrite the std::vector< Triplet > structure to one suitable
	 *  for the other load function.
	 *  @see load( row, col, val, m, n, nz ) 
	 *  @see SparseMatrix::load 
	 */
	virtual void load( std::vector< Triplet< _t_value > >& input, unsigned long int m, unsigned long int n, _t_value zero ) {
		int nz = input.size();
		int* row = new int[ nz ];
		int* col = new int[ nz ];
		_t_value* val = new _t_value[ nz ];
		unsigned long int c = 0;
		typename std::vector< Triplet< _t_value > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, c++ ) {
			row[ c ] = (*it).i();
			col[ c ] = (*it).j();
			val[ c ] = (*it).value;
		}
		load( row, col, val, m, n, nz, zero );
	}

	/** @see BICRS( row, col, val, m, n, nz ) */
	void load( int* row, int* col, _t_value* val, int m, int n, int nz, _t_value zero ) {
#ifdef _DEBUG
		std::cerr << "Warning: _DEBUG flag set." << std::endl;
#endif
		this->zero_element = zero;
		this->nnz = nz;
		this->nor = m;
		this->noc = n;
		this->ntt=2*n;
		int jumps = 0;
		int prevrow = row[ 0 ];
		for( unsigned long int i=1; i<this->nnz; i++ ) {
			if( row[ i ] != prevrow )
				jumps++;
			prevrow = row[ i ];
		}
#ifdef _DEBUG
		std::cout << jumps << " row jumps found." << std::endl;
#endif
		r_inc = new int[ jumps + 2 ];
		c_inc = new int[ this->nnz + 1 ];
		vals  = new double[ this->nnz ];
		for( unsigned long int i=0; i<this->nnz; ++i ) vals[i] = val[i];
	
		r_inc[ 0 ] = row[ 0 ];
		prevrow = row[ 0 ];	
		int prevcol = c_inc[ 0 ] = col[ 0 ];
#ifdef _DEBUG
		std::cout << "c_inc: " << prevcol << std::endl;
		std::cout << "r_inc: " << prevrow << std::endl;
#endif
		int c = 1;
		for( unsigned long int i=1; i<this->nnz; i++ ) {
			this->c_inc[ i ] = col[ i ] - prevcol;
			if( row[ i ] != prevrow ) {
				this->c_inc[ i ] += ntt;
				this->r_inc[ c++ ] = row[ i ] - prevrow;
#ifdef _DEBUG
				std::cout << "c_inc: " << ntt << std::endl;
				std::cout << "r_inc: " << row[ i ] - prevrow << std::endl;
#endif
				prevrow = row[ i ];
			}
#ifdef _DEBUG
			else
				std::cout << "c_inc: " << col[ i ] - prevcol << std::endl;
#endif
			prevcol = col[ i ];
		}
		//overflow so to signal end of matrix
		c_inc[ this->nnz ] = ntt;
		//initialise last two row jumps to zero (prevent undefined jump)
		r_inc[ c ] = 0;

#ifdef _DEBUG
		std::cout << "Construction done." << std::endl;
#endif
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = this->r_inc[ 0 ];
		col = this->c_inc[ 0 ];
	}

	/**
	 *  Calculates y=xA, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zxa( _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * y		= y_p;
		int *__restrict__ c_inc_p	= c_inc;
		int *__restrict__ r_inc_p	= r_inc;
		_t_value *__restrict__ v_p	= vals;

#ifndef NDEBUG
		const _t_value * x				= x_p;
		const _t_value * const x_end			= x+this->nor;
		const int *__restrict__      const c_inc_end	= c_inc+this->nnz+1;
#endif
		const _t_value * const y_end			= y+this->noc;
		const _t_value *__restrict__ const v_end	= vals+this->nnz;

		y_p += *c_inc_p++;
		x_p += *r_inc_p++;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( c_inc_p <  c_inc_end );
			assert( r_inc_p >= r_inc );
			while( y_p < y_end ) {
#ifdef _DEBUG
				std::cout << (x_p-x) << "," << (y_p-y) << " next increment: " << (*(c_inc_p+1))<< std::endl;
#endif
				*y_p += *v_p++ * *x_p;
				y_p += *c_inc_p++;
			}
			y_p -= ntt;
			x_p += *r_inc_p++;
		}
	}

	/**
	 *  Calculates y=Ax, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zax( _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * x		= x_p;
		int *__restrict__ c_inc_p	= c_inc;
		int *__restrict__ r_inc_p	= r_inc;
		_t_value *__restrict__ v_p	= vals;

#ifndef NDEBUG
		const _t_value * y				= y_p;
		const _t_value * const y_end			= y+this->nor;
		const int *__restrict__ const c_inc_end		= c_inc+this->nnz+1;
#endif
		const _t_value * const x_end			= x+this->noc;
		const _t_value *__restrict__ const v_end	= vals+this->nnz;

		x_p += *c_inc_p++;
		y_p += *r_inc_p++;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( c_inc_p <  c_inc_end );
			assert( r_inc_p >= r_inc );
			while( x_p < x_end ) {
#ifdef _DEBUG
				std::cout << (y_p-y) << "," << (x_p-x) << " next increment: " << (*(c_inc_p+1))<< std::endl;
#endif
				*y_p += *v_p++ * *x_p;
				 x_p += *c_inc_p++;
			}
			x_p -= ntt;
			y_p += *r_inc_p++;
		}
	}

};

#endif

