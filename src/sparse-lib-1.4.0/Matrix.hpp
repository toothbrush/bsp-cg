/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_M
#define _H_M

/**
 *  Operations common to all matrices.
 */
template< typename T >
class Matrix {
	private:
	
	protected:

	public:

		/** Base constructor. */
		Matrix() {}

		/** Base deconstructor. */
		virtual ~Matrix() {}

		/** @return The number of rows. */
		virtual unsigned long int m() = 0;

		/** @return The number of columns. */
		virtual unsigned long int n() = 0;

		/** @return The number of nonzeros. */
		virtual unsigned long int nzs() {
			return m() * n();
		}

		/**
		 *  Calculates z=Ax (where A is this matrix).
		 *  @param x The input vector.
		 *  @return The output vector z.
		 */
		virtual T* mv( T* x ) = 0;

		/**
		 *  In-place z=xA function.
		 *
		 *  @param x The input vector to left-multiply to the current matrix.
		 *  @param z The output vector. Must be pre-allocated and initialised to zero for correct results.
		 */
		virtual void zxa( T*__restrict__ x, T*__restrict__ z ) = 0;

		/**
		 *  In-place z=Ax function.
		 *
		 *  @param x The x vector to multiply current matrix with.
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zax( T*__restrict__ x, T*__restrict__ z ) = 0;

};

#endif //_H_M

