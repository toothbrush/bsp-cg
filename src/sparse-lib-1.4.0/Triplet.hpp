
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include<cstdlib>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>

//#define _DEBUG

#ifndef _H_TRIPLET
#define _H_TRIPLET

/** The triplet scheme. Stores sparse matrices by storing in an array format
 *  triplets of the form (i,j,v), with (i,j) the position of the nonzero and
 *  v the nonzero value. Triplets can be stored in any order.
 */
template< typename T >
class Triplet {

   protected:
	
	/** The row coordinate of this triplet. */
	unsigned long int row;
	
	/** The column coordinate of this triplet. */
	unsigned long int column;

   public:

#ifdef TRIPLET_META
	/** Stores meta data */
	TRIPLET_META_TYPE meta;
#endif

	/** @return Row index of this triplet. */	
	const unsigned long int i() const { return row; }

	/** @return Column index of this triplet. */
	const unsigned long int j() const { return column; }

	/** Value stored at this triplet. */
	T value;

	/** Base constructor.
	 *  @param i row index.
	 *  @param j column index.
	 *  @param val nonzero value.
	 */
	Triplet( unsigned long int ii, unsigned long int ij, T val ): row( ii ), column( ij ), value( val ) {}

	/** Base constructor. Sets all values to zero. */
	Triplet(): row( 0 ), column( 0 ), value( 0 ) {}

	/** Loads an array of triplets from a binary file.
	 *  Warning: there is a difference between 32 and 64 bits files!
	 *  @param fn Filename of the file to load from.
	 *  @param m Reference to where the total number of rows is to-be stored.
	 *  @param n Reference to where the total number of columns is to-be stored.
	 *  @return A std::vector containing the triplets in the order they were loaded in.
	 */
	static std::vector< Triplet< T > > load( const std::string fn, unsigned long int &m, unsigned long int &n ) {
		std::fstream file( fn.c_str(), std::ios::in | std::ios::binary );
		unsigned long int i; unsigned long int j; double v;
		std::vector< Triplet< T > > ret;
		if( !file.is_open() ) {
			std::cerr << "Error while opening file" << std::endl;
			exit( 1 );
		}
		file.read( (char*) &i, sizeof( unsigned long int ) );
		file.read( (char*) &j, sizeof( unsigned long int ) );
		m = i;
		n = j;
#ifdef _DEBUG
		std::cout << "m: " << m << ", n: " << n << std::endl;
#endif
		while( true ) {
			file.read( (char*) &i, sizeof( unsigned long int ) );
			if( !file ) break;
			file.read( (char*) &j, sizeof( unsigned long int ) );
			file.read( (char*) &v, sizeof( T ) );
#ifdef _DEBUG
			std::cout << "Pushed back: ( " << i << " , " << j << " , " << v << " )" << std::endl;
#endif
			ret.push_back( Triplet< T >( i, j, v ) );
		}
		file.close();
		return ret;
	}

	/**
	 *  Saves an array of triplets to a file, in binary format.
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Array of triplets to write.
	 *  @param m Total number of rows in the matrix.+
	 *  @param n Total number of columns in the matrix.
	 *  @param s Size of the array toWrite.
	 */
        static void save( std::string fn, Triplet< T >* toWrite, const unsigned long int m, const unsigned long int n, const unsigned long int s ) {
                std::fstream myFile ( fn.c_str(), std::ios::out | std::ios::binary);
                myFile.write( (char*) &m, sizeof( unsigned long int ) );
                myFile.write( (char*) &n, sizeof( unsigned long int ) );
                for( unsigned long int i = 0; i<s; i++ ) {
                        const unsigned long int wi = toWrite[ i ].i();
                        const unsigned long int wj = toWrite[ i ].j();
                        const double wv = toWrite[ i ].value;
#ifdef _DEBUG
                        std::cout << "Wrote: ( " << wi << " , " << wj << " , " << wv << " ) " << std::endl;
#endif
                        myFile.write( (char*) &( wi ), sizeof( unsigned long int ) );
                        myFile.write( (char*) &( wj ), sizeof( unsigned long int ) );
                        myFile.write( (char*) &( wv ), sizeof( T ) );
                }
                myFile.close();
        }

	/** Transposes this triplet, i.e., swapping the row and column value. */
	void transpose() { const unsigned long int t = row; row = column; column = t; }

	/**
	 *  Saves a std::vector of triplets to a file, in binary format.
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Vector of triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 */
        static void save( std::string fn, std::vector< Triplet< T > > &toWrite, const unsigned long int m, const unsigned long int n ) {
                save( fn, &( toWrite[ 0 ] ), m, n, toWrite.size() );
        }

};

#endif

