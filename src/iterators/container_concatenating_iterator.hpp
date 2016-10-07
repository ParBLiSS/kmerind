/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    container_concatenating_iterator.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   an iterator template for concatenating multiple nested iterators.  Meant to be used to create specialized iterators.
 * @details specialization is via functors that check for end of a top level iterator instance, and for initializing a lower level iterator
 *
 */
#ifndef CONTAINER_CONCATENATING_ITERATOR_HPP_
#define CONTAINER_CONCATENATING_ITERATOR_HPP_

#include <vector>
#include <iterator>
#include <utility>
#include "utils/logging.h"
#include <type_traits>  // for testing const iterator

namespace bliss
{
  namespace iterator
  {

  	  /// generic container adaptor for ContainerConcatenatingIterator.  See KmerParser for other exmples of ContainerAdaptor.
  	  template <typename DUMMY>
  	  class ConcatenatingIteratorContainerAdapter {
  	  public:
  		  template <typename Container>
  		  using iterator_type = typename Container::const_iterator;

  		  template <typename Container>
  		  iterator_type<Container> begin(Container const & c) const {
  			  return c.begin();
  		  }
  		  template <typename Container>
  		  iterator_type<Container> end(Container const & c) const {
  			  return c.end();
  		  }

  	  };


    /**
     * @class    bliss::iterator::ConcatenatingIterator
     * @brief    this class presents a single/sequential view of a series of underlying iterator ranges
     * @details  forward iterator only.  does not support going backward therefore not bidir or random access.
     *			TODO: InnerContainerAdaptor is currently mostly KmerParserTypes, which has begin and end operations that take a SeqType object.
     *					ContainerConcatenatingIterator is therefore not yet compatible with stl containers as inner containers directly.
     *					an adaptor would be needed.
     */
    template<typename OuterIter, typename InnerContainerAdapter >
    class ContainerConcatenatingIterator :
    		public ::std::iterator<
			::std::forward_iterator_tag,
			typename ::std::iterator_traits<
				typename InnerContainerAdapter::template iterator_type<
					typename ::std::iterator_traits<OuterIter>::value_type
				>
    		>::value_type,
			std::ptrdiff_t
        >
    {
      protected:
    	/// outer iterator's data type
    	using outer_value_type = typename ::std::iterator_traits<OuterIter>::value_type;
    	using InnerIter = typename InnerContainerAdapter::template iterator_type<outer_value_type>;
//    	using iterator_category = std::forward_iterator_tag;
//    	using value_type = typename ::std::iterator_traits<InnerIter>::value_type;
//    	using pointer = typename ::std::iterator_traits<InnerIter>::pointer;
//    	using reference = typename ::std::iterator_traits<InnerIter>::reference;
//    	using difference_type = std::ptrdiff_t;

    	using base_traits = std::iterator_traits<InnerIter>;

    	InnerContainerAdapter adapter;

    	/// current outer iterator
    	OuterIter o_curr;
    	/// maximum outer iterator
    	OuterIter o_end;

    	/// current inner iterator (associated with o_curr)
    	InnerIter i_curr;

    	/// end inner iterator (associated with o_curr).  mostly here so we don't call adapter.end() every call to operator++
    	InnerIter i_curr_end;

    	/// end inner iterator (associated with o_end)
    	InnerIter i_end;

    	/// indicate if i_curr is set, implies i_curr_end is set.  not set if o_curr is undereferenceable.
    	bool o_curr_dereferenceable;
    	/// indicate if i_end is set.  not set if o_end is undereferenceable.
    	bool o_end_dereferenceable;

    	inline bool are_same(
    			OuterIter const & ox, OuterIter const & oy,
				InnerIter const & ix, InnerIter const & iy,
    			bool bx, bool by) const {

    		// if outer iters are the same, then check if both inneriters are set, if so compare inneriters,
    		//       else return true if both are unset, false if one of them is unset.
    		// else return false;
    		//return (ox == oy) && ((bx && by) ? (ix == iy) : !(bx || by));


            // problem:  2 ways of specifying end (dereferenceable, or not).   if the other iterator is not specified the same way,
            // then they are not comparable and this could cause infinite loop.
            // the comparison needs to be symmetric....
            // probably have to use dereferenceability to indicate that inner iterators should not be checked.
    		// below fixes the problem
    		// outer have to match.  then either inner iterator match, or at least one is not dereferenceble.
    		return (ox == oy) && (!(bx && by) || (ix == iy));
    	}


    	inline bool is_contained(OuterIter oiter, InnerIter iiter) const {
    		return (std::distance(adapter.begin(*oiter), iiter) >= 0) &&
    			   (std::distance(iiter, adapter.end(*oiter)) >= 0);
    	}

        /// enforce that iterator is at a dereferenceable position
    	/// this means that the
        void ensure_dereferenceable() {
        	bool at_end = are_same(o_curr, o_end, i_curr, i_end, o_curr_dereferenceable, o_end_dereferenceable);

        	/*
        	 * check to see if position is dereferenceable.
        	 * dereferenceability table.  requires that o_curr_deref is properly updated.
        	 * 	| deref	| o_curr_deref	|  i_curr != i_curr_end
        	 * 	|   N	| 	n
        	 * 	|   Y 	|   y			|  		!=
        	 * 	|	N	|   y			|		==
        	 */
        	bool dereferenceable = (o_curr_dereferenceable && (i_curr != i_curr_end));

        	while (!dereferenceable && !at_end) {
        		// if current is not dereferenceable, and we are not at end, then the i_curr must be at i_curr_end.

        		// move to next outer iter.
        		++o_curr;

        		// if o_curr at o_end, use its deref state, else assumed dereferenceable.
        		// so in other words, o_curr not dereferenceable if o_curr is at o_end and o_end is not dereferenceable.
        		//   o_curr NOT deref if ( (o_curr == o_end) && !o_end_deref )
        		o_curr_dereferenceable = (o_curr != o_end) || o_end_dereferenceable;

        		// if dereferenceable, get the new i_curr and i_curr_end
        		if (o_curr_dereferenceable) {
        			i_curr = adapter.begin(*o_curr);
        			i_curr_end = adapter.end(*o_curr);
        		}

        		// dereferenceable if o_curr dereferecenable and inner iter not at end.
        		dereferenceable = o_curr_dereferenceable && (i_curr != i_curr_end);

        		// o_curr has to be at o_end, and in that case we are at end if o_curr not deref, or i_curr == i_end
        		at_end = (o_curr == o_end) && (!o_curr_dereferenceable || (i_curr == i_end));
        	}
        }


      public:

        /// default constructor for ForwardIterator.
    	ContainerConcatenatingIterator() :
    		o_curr(), o_end(o_curr),
    		o_curr_dereferenceable(false), o_end_dereferenceable(false) {};

    	/// Default destructor
    	~ContainerConcatenatingIterator() {};

        /// constructor for start concatenating iterator using copy semantic
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter,
        		OuterIter obegin, InnerIter ibegin, OuterIter oend, InnerIter iend)
        : adapter(_adapter), o_curr(obegin), o_end(oend),
		  i_curr(ibegin), i_curr_end(adapter.end(*o_curr)), i_end(iend),
		  o_curr_dereferenceable(true), o_end_dereferenceable(true) {

        	if (!is_contained(o_curr, i_curr)) throw std::invalid_argument("ibegin parameter is not within obegin container");
        	if (!is_contained(o_end, i_end)) throw std::invalid_argument("iend parameter is not within oend container");

        	ensure_dereferenceable();
        };


        /// constructor for start concatenating iterator using copy semantic.  outer end iter is assumed to be not dereferenceable.
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter,
        		OuterIter obegin, InnerIter ibegin, OuterIter oend)
        : adapter(_adapter), o_curr(obegin), o_end(oend),
		  i_curr(ibegin), i_curr_end(adapter.end(*o_curr)),
		  o_curr_dereferenceable(true), o_end_dereferenceable(false) {
        	if (!is_contained(o_curr, i_curr)) throw std::invalid_argument("ibegin parameter is not within obegin container");

        	ensure_dereferenceable();
        };

        /// constructor for start concatenating iterator using copy semantic.  outer begin iterator is assumed to be dereferenceable.
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter,
        		OuterIter obegin, OuterIter oend, InnerIter iend)
        : adapter(_adapter), o_curr(obegin), o_end(oend),
		  i_curr(adapter.begin(*o_curr)), i_curr_end(adapter.end(*o_curr)), i_end(iend),
		  o_curr_dereferenceable(true), o_end_dereferenceable(true) {
        	if (!is_contained(o_end, i_end)) throw std::invalid_argument("iend parameter is not within oend container");

        	ensure_dereferenceable();
        };


        /// constructor for start concatenating iterator using copy semantic.  assume outer end iter is not deferenceable.
        /// assume outer begin iter is dereferenceable unless it's same as outer end iter.
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter,
        		OuterIter obegin, OuterIter oend)
        : adapter(_adapter), o_curr(obegin), o_end(oend),
		  o_curr_dereferenceable(false), o_end_dereferenceable(false) {
        	if (o_curr != o_end) {
        		// not at end, can get inner begin and end
        		i_curr = adapter.begin(*o_curr);
        		i_curr_end = adapter.end(*o_curr);
        		o_curr_dereferenceable = true;

            	ensure_dereferenceable();
        	}
        };

        /// constructor for end concatenating iterator using copy semantic
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter, OuterIter oend, InnerIter iend)
        : adapter(_adapter), o_curr(oend), o_end(oend),
		  i_curr(iend), i_curr_end(iend), i_end(iend),
		  o_curr_dereferenceable(true), o_end_dereferenceable(true) {
        	if (!is_contained(o_end, i_end)) throw std::invalid_argument("iend parameter is not within oend container");

        };

        /// constructor for end concatenating iterator using copy semantic
        ContainerConcatenatingIterator(InnerContainerAdapter const & _adapter, OuterIter oend)
        : adapter(_adapter), o_curr(oend), o_end(oend), o_curr_dereferenceable(false), o_end_dereferenceable(false) {};


        /// constructor for end iterator, copy semantic
        ContainerConcatenatingIterator(ContainerConcatenatingIterator const & other) :
        	adapter(other.adapter), o_curr(other.o_curr), o_end(other.o_end),
			i_curr(other.i_curr), i_curr_end(other.i_curr_end), i_end(other.i_end),
			o_curr_dereferenceable(other.o_curr_dereferenceable), o_end_dereferenceable(other.o_end_dereferenceable) {};


        // note that explicit keyword cannot be on copy and move constructors else the constructors are not defined/found.

        /// copy constructor.  respects multi-pass requirement of forward iterator.
        ContainerConcatenatingIterator& operator=(ContainerConcatenatingIterator const & other)
        {
        	adapter 	  = other.adapter;
        	o_curr		  = other.o_curr;
            o_end		  = other.o_end;
            i_curr		  = other.i_curr;
            i_curr_end    = other.i_curr_end;
            i_end		  = other.i_end;
            o_curr_dereferenceable = other.o_curr_dereferenceable;
            o_end_dereferenceable  = other.o_end_dereferenceable;

          return *this;
        }

        /**
         * @brief increment:  move to the next position in the concatenating iterator, which may cross range boundaries.
         * @return reference to self.
         */
        ContainerConcatenatingIterator& operator++()
		{
          //=== first set of conditions are when we are at the end.
          if (are_same(o_curr, o_end, i_curr, i_end, o_curr_dereferenceable, o_end_dereferenceable)) return *this;

          // if not at end, then o_curr must be dereferenceable from constructor and previous calls to ++
          // this also means that i_curr is not at i_curr_end.

          // so increment i_curr, then ensure dereferenceable.
          ++i_curr;
          ensure_dereferenceable();

          return *this;
		}


        /**
         * post increment.  make a copy then increment that.
         */
        ContainerConcatenatingIterator operator++(int)
		{
          ContainerConcatenatingIterator output(*this);
          this->operator++();
          return output;
		}

        //=== input iterator specific

        /// comparison operator
        inline bool operator==(const ContainerConcatenatingIterator& rhs) const
		{
        	return are_same(o_curr, rhs.o_curr, i_curr, rhs.i_curr, o_curr_dereferenceable, rhs.o_curr_dereferenceable);
//        	if (o_curr == rhs.o_curr) {
//        		if (o_curr_dereferenceable && rhs.o_curr_dereferenceable)
//        			return i_curr == rhs.i_curr;
//        		else if (! (o_curr_dereferenceable || rhs.o_curr_dereferenceable))
//        			return true;
//        	}
//        	return false;

		}

        /// comparison operator
        inline bool operator!=(const ContainerConcatenatingIterator& rhs) const
		{
          return !(this->operator==(rhs));
		}

        /// get pointer.  should be dereferenceable, unless we're at end.
        inline typename base_traits::pointer operator->() const
        {
          return &(*i_curr);
        }

        /// dereference operator.  normal (not an output iterator) - returns a rvalue
        //  should be dereferenceable unless we're at end.
        inline typename base_traits::value_type operator*() const
        {
          return *i_curr;
        }



    };

  } /* namespace iterator */
} /* namespace bliss */

#endif /* CONTAINER_CONCATENATING_ITERATOR_HPP_ */
