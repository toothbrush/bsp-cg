/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_SM
#define _H_SM

#include "Matrix.hpp"
#include "FileToVT.hpp"

/**
 *  Interface common to all sparse matrix storage schemes.
 */
template< typename T, typename IND >
class SparseMatrix: public Matrix< T > {
	private:
	
	protected:
		/** Number of rows. */
		IND nor;

		/** Number of columns */
		IND noc;

		/** Number of non-zeros. */
		IND nnz;

	public:

		/** The element considered to be zero. */
		T zero_element;

		/** Base constructor. */
		SparseMatrix() {}

		/** Base constructor. */
		SparseMatrix( const IND nzs, const IND nr, const IND nc, const T zero ):
			nnz( nzs ), nor( nr ), noc( nc ), zero_element( zero ) {}

		/** Base deconstructor. */
		virtual ~SparseMatrix() {}

		/** 
		 *  Function reading in from std::vector< Triplet< T > > format.
		 *  @param input The input matrix in triplet format.
		 *  @param m The number of rows of the input matrix.
		 *  @param n The number of columns of the input matrix.
		 *  @param zero Which element is to be considered zero.
		 */
		virtual void load( std::vector< Triplet< T > >& input, const IND m, const IND n, const T zero ) = 0;

		/**
		 *  Function which loads a matrix from a matrix market file.
		 *  @param file Filename (including path) to the matrix market file.
		 *  @param zero Which element is to be considered zero.
		 */
		void loadFromFile( const std::string file, const T zero=0 ) {
			unsigned long int m, n;
			std::vector< Triplet< T > > VT = FileToVT::parse( file, m, n );
			this->load( VT, m, n, zero );
		}

		/** @return The number of rows. */
		virtual unsigned long int m() {
			return static_cast< unsigned long int >( nor );
		}

		/** @return The number of columns. */
		virtual unsigned long int n() {
			return static_cast< unsigned long int >( noc );
		}

		/** @return The number of nonzeros. */
		virtual unsigned long int nzs() {
			return static_cast< unsigned long int >( nnz );
		}

		/** Returns the first nonzero index, per reference. */
		virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) = 0;

		/**
		 *  Calculates and returns z=Ax. The vectors x should have appropiately set length; this is not
		 *  enforced. Memory leaks, segmentation faults, the works; one or more of these will occur if dimensions
		 *  do not make sense.
		 *  The return array is allocated to length equal to the number of rows in this function.
		 *  @param x The input (dense) x-vector.
		 *  @return The matrix-vector multiplication Ax, where A corresponds to the currently stored matrix.
		 *  @see mv
		 */
		virtual T* mv( T* x ) {
			T* ret = new T[ nor ];
			for( IND i=0; i<nor; i++ )
				ret[ i ] = zero_element;
			zax( x, ret );
			return ret;
		}

		/**
		 *  In-place z=Ax function.
		 *
		 *  @param x The x vector to multiply current matrix with.
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zax( T*__restrict__ x, T*__restrict__ z ) = 0;

		/**
		 *  In-place z=xA function.
		 *
		 *  @param x The x vector to apply in left-multiplication to A
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zxa( T*__restrict__x, T*__restrict__ z ) = 0;

};

#endif //_H_SM

