
/* Copyright by A.N. Yzelman, Dept. of Mathematics, Utrecht University, 2007-2009.
 * Released under GPL, see the file LICENSE.txt
 */

#ifndef _H_BLOCKORDERER
#define _H_BLOCKORDERER

#include "SBDTree.hpp"
#include "Triplet.hpp"
#include <vector>
#include <map>

#define NORMAL_DS 2 //ICRS
#define HORIZONTAL_DS -2 //ICCS
#define VERTICAL_DS 2 //ICRS

/** Induces a block order by fully traversing an SBDTree */
template< typename T >
class BlockOrderer {
	
	protected:

		SBDTree *tree;
		void prefix( const unsigned long int index ) {
			switch( traverse_mode ) {
			case TRAVERSE_HEIGHT:
				++cur_height;
				break;
			case READOUT:
				pre_readout( index );
				break;
			default:
				std::cerr << "BlockOrderer:: Unknown or unset traverse mode!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		void infix( const unsigned long int index ) {
			switch( traverse_mode ) {
			case TRAVERSE_HEIGHT:
				height[ index ] = cur_height;
				break;
			case READOUT:
				in_readout( index );
				break;
			default:
				std::cerr << "BlockOrderer:: Unknown or unset traverse mode!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		void postfix( const unsigned long int index ) {
			switch( traverse_mode ) {
			case TRAVERSE_HEIGHT:
				--cur_height;
				break;
			case READOUT:
				post_readout( index );
				break;
			default:
				std::cerr << "BlockOrderer:: Unknown or unset traverse mode!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		void traverse() {
			char *leftpass  = new char[ tree->size() ];
			char *rightpass = new char[ tree->size() ];
			for( unsigned long int i=0; i<tree->size(); i++ ) leftpass[i] = rightpass[i] = 0;

			const unsigned long int root = tree->getRoot();
			unsigned long int cur = root;
		
			while( true ) { //traverse! :)
				//if this is a leaf node, call all the functions and move up
				if( tree->isLeaf( cur ) ) {
					prefix ( cur );
					infix  ( cur );
					postfix( cur );
					if( cur == root ) {
						delete [] leftpass;
						delete [] rightpass;
						return; //done
					} else
						cur = tree->up( cur ); //move up
				} else {
					//this is an internal node
					//if not passed left yet
					if( !leftpass[ cur ] ) {
						prefix( cur ); //execute prefix function
						leftpass[ cur ] = 1; //set left pass
						cur = tree->left( cur ); //go to left subtree
					} else {
						//if not passed right yet
						if( !rightpass[ cur ] ) {
							infix( cur ); //execute infix function
							rightpass[ cur ] = 1; //set right pass
							cur = tree->right( cur ); //go to right subtree
						} else {
							//passed children
							postfix( cur ); //execute post function
							if( cur == root ) {
								delete [] leftpass;
								delete [] rightpass;
								return; //done
							} else
								cur = tree->up( cur ); //move up
						}
					}
				}
			}
			std::cerr << "BlockOrderer::traverse: Reached what should be unreachable code!" << std::endl;
			assert( false );
			exit( EXIT_FAILURE );
		}

		char traverse_mode;

		/** Helpers for determining node heights */
		static const char TRAVERSE_HEIGHT = 0;
		unsigned long int *height;
		unsigned long int cur_height;
		void pre_height ( const unsigned long int index );
		void in_height  ( const unsigned long int index );
		void post_height( const unsigned long int index );

		/** Helpers for final tree read-out */
		static const char READOUT = 1;
		/** Following depends on exact block order */
		virtual void pre_readout ( const unsigned long int index ) = 0;
		virtual void in_readout  ( const unsigned long int index ) = 0;
		virtual void post_readout( const unsigned long int index ) = 0;
		std::vector< Triplet< T > > *items;
		std::vector< std::vector< Triplet< T > > > *output;
		std::vector< signed char > *datatype;

		/** Helper functions for determining place of a nonzero within a separator cross */
		inline bool left_horizontal( const unsigned long int index, const Triplet< T > triplet ) {
			unsigned long int bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi, r_lo, r_hi, c_lo, c_hi;
			this->tree->getSeparatorBB( index, bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi );
			this->tree->rowBounds( index, r_lo, r_hi );
			this->tree->columnBounds( index, c_lo, c_hi );
			return ( triplet.i() >= r_lo && triplet.i() < r_hi && triplet.j() >= bb_c_lo && triplet.j() < c_lo );
		}

		inline bool right_horizontal( const unsigned long int index, const Triplet< T > triplet ) {
			unsigned long int bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi, r_lo, r_hi, c_lo, c_hi;
			this->tree->getSeparatorBB( index, bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi );
			this->tree->rowBounds( index, r_lo, r_hi );
			this->tree->columnBounds( index, c_lo, c_hi );
			return ( triplet.i() >= r_lo && triplet.i() < r_hi && triplet.j() >= c_hi && triplet.j() < bb_c_hi );
		}

		inline bool upper_vertical( const unsigned long int index, const Triplet< T > triplet ) {
			unsigned long int bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi, r_lo, r_hi, c_lo, c_hi;
			this->tree->getSeparatorBB( index, bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi );
			this->tree->rowBounds( index, r_lo, r_hi );
			this->tree->columnBounds( index, c_lo, c_hi );
			return ( triplet.i() >= bb_r_lo && triplet.i() < r_lo && triplet.j() >= c_lo && triplet.j() < c_hi );
		}

		inline bool lower_vertical( const unsigned long int index, const Triplet< T > triplet ) {
			unsigned long int bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi, r_lo, r_hi, c_lo, c_hi;
			this->tree->getSeparatorBB( index, bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi );
			this->tree->rowBounds( index, r_lo, r_hi );
			this->tree->columnBounds( index, c_lo, c_hi );
			return ( triplet.i() >= r_hi && triplet.i() < bb_r_hi && triplet.j() >= c_lo && triplet.j() < c_hi );
		}

		inline bool middle( const unsigned long int index, const Triplet< T > triplet ) {
			unsigned long int bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi, r_lo, r_hi, c_lo, c_hi;
			this->tree->getSeparatorBB( index, bb_r_lo, bb_r_hi, bb_c_lo, bb_c_hi );
			this->tree->rowBounds( index, r_lo, r_hi );
			this->tree->columnBounds( index, c_lo, c_hi );
			return ( triplet.i() >= r_lo && triplet.i() < r_hi && triplet.j() >= c_lo && triplet.j() < c_hi );
		}
		//Note: I'm perhaps expecting too much from the compiler; the inline keyword hopefully
		//      hints enough that these functions should really be integrated into the pre/in/post
		//	readout functions, outside the for loop along the nonzeroes...

	public:

		/** Base constructor */
		BlockOrderer(): tree( NULL ), height( NULL ), items( NULL ), output( NULL ), datatype( NULL ) {}

		/** Base destructor */
		~BlockOrderer() {
			if( tree != NULL ) delete tree;
			if( height != NULL ) delete [] height;
			if( items != NULL ) delete [] items;
			if( output != NULL ) {
				std::cerr << "Warning: BlockOrderer::output was somehow not reset to NULL!" << std::endl;
			}
		}

		/** Induces this ordering on a set of nonzeroes */
		std::vector< std::vector< Triplet< T > > > induce( const std::vector< Triplet< T > > input,
			std::vector< unsigned long int > &r_hierarchy, std::vector< unsigned long int > &c_hierarchy,
			std::vector< unsigned long int > &r_bounds, std::vector< unsigned long int > &c_bounds,
			std::vector< signed char > *datatype ) {
	
			//Note:	
			//this function only passes around Triplet< T >'s, and does not have to do anything specific
			//depending on T. Therefore, this function is templated, and not the entire class.
			std::vector< std::vector< Triplet< T > > > ret;

#ifdef _INFO
			std::cout << "BlockOrderer, phase I (build tree)" << std::endl;
#endif
			//phase running time: O( p log p )?

			tree = new SBDTree( r_hierarchy, c_hierarchy, r_bounds, c_bounds );

#ifdef _INFO
			std::cout << "BlockOrderer, phase II (build height map and dictionary)" << std::endl;
#endif
			//phase running time: O( 2p-1 )

			//get height map
			height = new unsigned long int[ tree->size() ];
			for( unsigned long int i=0; i<tree->size(); ++i ) height[ i ] = ULONG_MAX;
			std::map< unsigned long int, unsigned long int > row_translate, col_translate;
			traverse_mode = TRAVERSE_HEIGHT;
			cur_height = 0;
			traverse();
			//check height map
			for( unsigned long int i=0; i<tree->size(); ++i ) assert( height[ i ] != ULONG_MAX );
			
			//build dictionaries
			unsigned long int lo, hi;
			for( unsigned long int i=0; i<tree->size(); ++i ) {
				tree->rowBounds( i, lo, hi );
				//prevent translating indices to empty blocks
				if( lo != hi )
					row_translate[ hi ] = i;
				tree->columnBounds( i, lo, hi );
				if( lo != hi )
					col_translate[ hi ] = i;
			}

#ifdef _INFO
			std::cout << "BlockOrderer, phase III (assign triplets)" << std::endl;
#endif
			//phase running time: O( nnz log p )

			items = new std::vector< Triplet< T > >[ tree->size() ];
			std::vector< Triplet< double > >::const_iterator it = input.begin();
			for( ; it != input.end(); ++it ) {
				const unsigned long int i = row_translate.upper_bound( it->i() )->second;
				const unsigned long int j = col_translate.upper_bound( it->j() )->second;
				if( i == j ) items[ i ].push_back( *it );
				else {
					//not on a leaf node; on a separator cross (but which?)
					//well, the one with the lowest height
					if( height[ i ] == height[ j ] ) assert( false );
					if( height[ i ] < height [ j ] )
						items[ i ].push_back( *it );
					else
						items[ j ].push_back( *it );
				}
			}
#ifdef _INFO
			std::cout << "BlockOrderer, phase IV (read out tree & return)" << std::endl;
#endif
			//phase running time: O( 2p-1 )

			output = &ret;
			this->datatype = datatype;
			traverse_mode = READOUT;
			traverse();
			output = NULL;
			this->datatype = NULL;

			return ret;
		}
};

#endif

