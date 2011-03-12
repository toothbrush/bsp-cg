
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2010.
 * Released under GPL, see the file LICENSE.txt
 */

#include "HilbertTriplet.hpp"

#ifndef _H_HILBERTTRIPLETCOMPARE
#define _H_HILBERTTRIPLETCOMPARE

template< typename T >
class HilbertTripletCompare {

	public:

	bool operator() ( HilbertTriplet< T > i, HilbertTriplet< T > j ) {
		const unsigned long int ihilbert1 = i.getMostSignificantHilbertBits();
		const unsigned long int jhilbert1 = j.getMostSignificantHilbertBits();
		if( ihilbert1 < jhilbert1 ) return true;
		if( ihilbert1 > jhilbert1 ) return false;
		const unsigned long int ihilbert2 = i.getLeastSignificantHilbertBits();
		const unsigned long int jhilbert2 = j.getLeastSignificantHilbertBits();
		if( ihilbert2 < jhilbert2 ) return true;
		return false;
	}
};

#endif

