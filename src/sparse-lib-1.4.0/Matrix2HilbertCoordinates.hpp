/*
 *  Copyright (C) 2010-2007  A.N. Yzelman
 *
 *  Last modified at 13th of October, 2010, by A.N. Yzelman
 */


#ifndef _H_MATRIX2HILBERTCOORDINATES
#define _H_MATRIX2HILBERTCOORDINATES

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>

/**
 *  Class which maps coordinates to 1D Hilbert Coordinates.
 */
class Matrix2HilbertCoordinates {

   private:

	/** Base constructor. Does nothing. */
	Matrix2HilbertCoordinates() {}
   
   protected:

	static const unsigned char BITWIDTH = 64; //Warning: this should not exceed the bit length
						  //         of an unsigned long int on the system
						  //         used. Understand the algorithm before
						  //         changing this value.

   public:

	/** 
	 *  New method, October 2010. Maps any 2D coordinate (i,j),
	 *  with i and j 64-bits unsigned integers,
	 *  to a 1D 128-bits unsigned integer.
	 *
	 *  @param i  A 64-bits unsigned integer value in one dimension
	 *  @param j  A 64-bits unsigned integer value in the other dimension
	 *  @param h1 First part of the 128-bit Hilbert coordinate, unsigned integer format (most significant, first 64 bits)
	 *  @param h2 Second part of the 128-bit Hilbert coordinate (least significant, last 64 bits)
	 */
	static void IntegerToHilbert( const unsigned long int i, const unsigned long int j,
					unsigned long int &h1, unsigned long int &h2 );

};

#endif

