/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#include <vector>

#include "Triplet.hpp"
#include "SparseMatrix.hpp"

#ifndef _H_CCSWRAPPER
#define _H_CCSWRAPPER

template< typename T, typename SparseMatrixType >
class CCSWrapper : public SparseMatrix< T, unsigned long int > {

	private:
		SparseMatrixType *ds;

	protected:
		void transposeVector( std::vector< Triplet< T > > &input ) {
			typename std::vector< Triplet< T > >::iterator it = input.begin();
			for( ; it != input.end(); ++it )
				it->transpose();
		}

	public:
		CCSWrapper() { ds = NULL; }

		CCSWrapper( std::string file, T zero = 0 ) { loadFromFile( file, zero ); }

		CCSWrapper( const unsigned long int nnz, const unsigned long int nor, const unsigned long int noc, T zero ) {
			ds = new SparseMatrixType( nnz, noc, nor, zero );
		}

		CCSWrapper( std::vector< Triplet< T > > &input, unsigned long int m, unsigned long int n, T zero ) {
			load( input, m, n, zero );
		}

		virtual ~CCSWrapper() {
			if ( ds != NULL ) delete ds;
		}

		virtual void load( std::vector< Triplet< T > > &input, unsigned long int m, unsigned long int n, T zero ) {
			transposeVector( input );
			ds = new SparseMatrixType( input, n, m, zero );
		}

		virtual void loadFromFile( const std::string file, const T zero = 0 ) {
			unsigned long int m, n;
			std::vector< Triplet< T > > VT = FileToVT::parse( file, m, n );
			load( VT, m, n, zero );
		}

		virtual unsigned long int m() { return ds->n(); }

		virtual unsigned long int n() { return ds->m(); }

		virtual unsigned long int nzs() { return ds->nzs(); }

		virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) { ds->getFirstIndexPair( col, row ); }
		
		virtual T* mv( T* x ) {
			T * const ret = new T[ m() ];
			for( unsigned long int i = 0; i < m(); i ++ ) ret[ i ] = ds->zero_element;
			ds->zxa( x, ret );
			return ret;
		}

		virtual void zax( T*__restrict__ x, T*__restrict__ z ) {
			ds->zxa( x, z );
		}

		virtual void zxa( T*__restrict__ x, T*__restrict__ z ) {
			ds->zax( x, z );
		}

};

#endif

