
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2008.
 * Released under GPL, see the file LICENSE.txt
 */

#include "Matrix2HilbertCoordinates.hpp"
#include<vector>
#include<string>
#include<fstream>

#ifndef _H_HILBERT_TRIPLET
#define _H_HILBERT_TRIPLET

/**
 *  Class which upon reading in a matrix, first determines the Hilbert-coordinate
 *  of the nonzero at (i,j). Triplets (i,j,v), where v is the nonzero value at
 *  (i,j), are then stored in order of their respective Hilbert coordinate.
 */
template< typename T >
class HilbertTriplet {

   protected:

	/** The row coordinate of this triplet. */
	unsigned long int row;

	/** The column coordinate of this triplet. */
	unsigned long int column;

	/** Most significant part of a 128-bits Hilbert coordinate, for one-shot, non-iterative, calculation. */
	unsigned long int hilbert1;

	/** Least significant part of a 128-bits Hilbert coordinate, for one-shot, non-iterative, calculation. */
	unsigned long int hilbert2;

   public:

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
	HilbertTriplet( unsigned long int i, unsigned long int j, T val ): row( i ), column( j ), hilbert1( 0 ), hilbert2( 0 ), value( val ) {}

	/** Base constructor. Sets all values to zero. */
	HilbertTriplet(): row( 0 ), column( 0 ), hilbert1( 0 ), hilbert2( 0 ), value( 0 ) {}

	/** Calculates the full Hilbert coordinate */
	void calculateHilbertCoordinate() {
		Matrix2HilbertCoordinates::IntegerToHilbert( row, column, hilbert1, hilbert2 );
	}

	/**
	 *	Gets the Hilbert coordinates.
	 *	Does not check if calculateHilbertCoordinate() was called first, otherwise (0,0) will be returned.
	 *	Note that the Hilbert coordinate is a 1D 128-bits unsigned integer.
	 *
	 *   @param h1 The first (most significant) 64-bits of the Hilbert coordinate
	 *   @param h2 the remainder of the Hilbert coordinate.
	 */
	void getHilbertCoordinate( unsigned long int &h1, unsigned long int &h2 ) {
		h1 = hilbert1;
		h2 = hilbert2;
	}

	/** @return h1 of HilbertTriplet<T>::getHilbertCoordinate() */
	const unsigned long int getMostSignificantHilbertBits() {
		return hilbert1;
	}

	/** @return h2 of HilbertTriplet<T>::getHilbertCoordinate() */
	const unsigned long int getLeastSignificantHilbertBits() {
		return hilbert2;
	}

	/**
	 *  Saves an array of Hilbert triplets to a file, in binary format.
	 *  Does NOT save the Hilbert coordinate!
	 *  (For reading-in the written file, use the regular Triplet scheme)
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Array of Hilbert triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 *  @param s Size of the array toWrite.
	 */
	static void save( std::string fn, HilbertTriplet< T >* toWrite, const unsigned long int m, const unsigned long int n, const unsigned long int s ) {
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

	/**
	 *  Saves a std::vector of Hilbert triplets to a file, in binary format.
	 *  Does NOT save the Hilbert coordinate!
	 *  (For reading-in the written file, use the regular Triplet scheme)
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Vector of Hilbert triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 */
	static void save( std::string fn, std::vector< HilbertTriplet< T > > &toWrite, const unsigned long int m, const unsigned long int n ) {
		save( fn, &( toWrite[ 0 ] ), m, n, toWrite.size() );
	}

};

#endif

