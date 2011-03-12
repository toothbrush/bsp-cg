
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_SBDTREE
#define _H_SBDTREE

#include <cstdlib>
#include <vector>
#include <iostream>
#include <limits.h>
#include <assert.h>

class SBDTree {

	protected:
		
		unsigned long int *parent;
		unsigned long int *left_child;
		unsigned long int *right_child;
		unsigned long int *r_lo;
		unsigned long int *r_hi;
		unsigned long int *c_lo;
		unsigned long int *c_hi;
		unsigned long int root;
		unsigned long int nodes;
		char root_is_set;

		static const unsigned long int NO_SUCH_NODE = ULONG_MAX;

		void build( std::vector< unsigned long int > &hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );
	public:

		/** Base constructor */
		SBDTree( std::vector< unsigned long int > &r_hierarchy, std::vector< unsigned long int > &c_hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );

		/** Base constructor. Warning: avoids some assertions! */
		SBDTree( std::vector< unsigned long int > &hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );

		/** Base deconstructor. */
		~SBDTree();

		/** Gets, from a separator node, the bounding box of the nonzeroes contained in the separator.
		    Note that this is *not* the r_lo/hi,c_lo/hi of the node itself; those are the bounds of the
		    row-wise and column-wise separators.
		
		    This is a logarithmic operation.
		  */
		void getSeparatorBB( const unsigned long int index,
					unsigned long int &r_lo, unsigned long int &r_hi,
					unsigned long int &c_lo, unsigned long int &c_hi );

		/** Returns the parent of a given node. */
		unsigned long int up( const unsigned long int index );

		/** Returns the left child of a given node. */
		unsigned long int left( const unsigned long int index );

		/** Returns the right child of a given node. */
		unsigned long int right( const unsigned long int index );

		/** Returns the row bounds corresponding to a given node. */
		void rowBounds( const unsigned long int index,
				unsigned long int &r_lo,
				unsigned long int &r_hi );

		/** Returns the column bounds corresponding to a given node. */
		void columnBounds( const unsigned long int index,
					unsigned long int &c_lo,
					unsigned long int &c_hi );

		/** Whether the given node is a leaf node */
		char isLeaf( const unsigned long int index );

		/** Gets the number of nodes */
		unsigned long int size();

		/** Gets the root index */
		unsigned long int getRoot();
};

#endif

