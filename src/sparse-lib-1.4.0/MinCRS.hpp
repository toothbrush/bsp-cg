
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_MINCRS
#define _H_MINCRS

#include "BlockOrderer.hpp"

/** Codes the Minimal CRS block order */
template< typename T >
class MinCRS: public BlockOrderer< T > {
	protected:
		virtual void pre_readout ( const unsigned long int index ) {
			//switch between internal nodes and external nodes
			if( this->tree->isLeaf( index ) ) {
				//immediately add everything in range
				this->output->push_back( this->items[ index ] );
				if( this->datatype != NULL )
					this->datatype->push_back( NORMAL_DS );
			} else {
				//initialise
				typename std::vector< Triplet< T > >::iterator it = this->items[ index ].begin();
				std::vector< Triplet< T > > cur;
				std::vector< Triplet< T > > rep;
				//loop over this node's triplets
				for( ; it != this->items[ index ].end(); ++it ) {
					//left horizontal separator first
					if( upper_vertical( index, *it ) )
						cur.push_back( *it );
					else
						rep.push_back( *it );
				}
				if( cur.size() > 0 )
					this->items[ index ] = rep; //replace with smaller vector
				this->output->push_back( cur );
				if( this->datatype != NULL )
					this->datatype->push_back( VERTICAL_DS );
			}
		}

		virtual void in_readout  ( const unsigned long int index ) {
			//skip leaf node
			if( this->tree->isLeaf( index ) ) return;
			//initialise
			typename std::vector< Triplet< T > >::iterator it = this->items[ index ].begin();
			std::vector< Triplet< T > > cur1;
			std::vector< Triplet< T > > cur2;
			std::vector< Triplet< T > > cur3;
			std::vector< Triplet< T > > rep;
			//loop over this node's triplets
			for( ; it != this->items[ index ].end(); ++it ) {
				//now the 3 vertical separators (incl middle) in row order
				if( left_horizontal( index, *it ) )
					cur1.push_back( *it );
				else if( middle( index, *it ) )
					cur2.push_back( *it );
				else if( right_horizontal( index, *it ) )
					cur3.push_back( *it );
				else
					rep.push_back( *it );
			}
			if( cur1.size() + cur2.size() + cur3.size() > 0 )
				this->items[ index ] = rep; //replace with smaller vector
			this->output->push_back( cur1 );
			this->output->push_back( cur2 );
			this->output->push_back( cur3 );
			if( this->datatype != NULL ) {
				this->datatype->push_back( HORIZONTAL_DS );
				this->datatype->push_back( NORMAL_DS );
				this->datatype->push_back( HORIZONTAL_DS );
			}
		}

		virtual void post_readout( const unsigned long int index ) {
			//skip leaf node
			if( this->tree->isLeaf( index ) ) return;
			//initialise
			typename std::vector< Triplet< T > >::iterator it = this->items[ index ].begin();
			std::vector< Triplet< T > > cur;
			std::vector< Triplet< T > > rep;
			//loop over this node's triplets
			for( ; it != this->items[ index ].end(); ++it ) {
				//right horizontal separator last
				if( lower_vertical( index, *it ) )
					cur.push_back( *it );
				else
					rep.push_back( *it );
			}
			if( cur.size() > 0 )
				this->items[ index ] = rep; //replace with smaller vector
			assert( rep.size() == 0 );
			this->output->push_back( cur );
			if( this->datatype != NULL )
				this->datatype->push_back( VERTICAL_DS );
		}

	public:
		//(de)constructor of superclass gets called automagically
};

#endif

