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
 * @file    compound_iterators.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   an iterator that returns only elements that match a predicate, and then applies a transform to the output.
 * @details mostly a convenience class.
 */

#ifndef COMPOUND_ITERATOR_HPP_
#define COMPOUND_ITERATOR_HPP_

#include <iterator>
#include "iterators/filter_iterator.hpp"
#include "iterators/transform_iterator.hpp"

namespace bliss
{

  namespace iterator
  {

  	  /**
  	   * @brief		filter_transform_iterator filters an iterator and transforms the output.
  	   * @details	e.g. get kmers from kmer_count pair whose count is above threshold.
  	   */
  	  template<typename Iterator,
			 typename Filter = ::fsc::TruePredicate,
			 typename Transform = ::bliss::kmer::transform::Identity<typename ::std::iterator_traits<Iterator>::value_type>
  		>
  	  using filter_transform_iterator = ::bliss::iterator::transform_iterator<
  			  ::bliss::iterator::filter_iterator<Filter, Iterator>,
			   Transform>;

  	  // make start iterator
  	  template<typename Iterator,
			 typename Filter = ::fsc::TruePredicate,
			 typename Transform = ::bliss::kmer::transform::Identity<typename ::std::iterator_traits<Iterator>::value_type>
  		>
  	  filter_transform_iterator<Iterator, Filter, Transform>
  	  make_filter_transform_iterator(Iterator curr, Iterator end,
  			  Filter const & f = Filter(), Transform const & t = Transform()) {
  		  return ::bliss::iterator::transform_iterator<Iterator, Transform>(
  				  ::bliss::iterator::filter_iterator<Iterator, Filter>(f, curr, end), t);
  	  }
  	  // make end iterator
  	  template<typename Iterator,
			 typename Filter = ::fsc::TruePredicate,
			 typename Transform = ::bliss::kmer::transform::Identity<typename ::std::iterator_traits<Iterator>::value_type>
  		>
  	  filter_transform_iterator<Iterator, Filter, Transform>
  	  make_filter_transform_iterator(Iterator end,
  			  Filter const & f = Filter(), Transform const & t = Transform()) {
  		  return ::bliss::iterator::transform_iterator<Iterator, Transform>(
  				  ::bliss::iterator::filter_iterator<Iterator, Filter>(f, end), t);
  	  }


  	  /**
  	   * @brief		elemental_filter_iterator - filter out elements in the iterator based on the content of another iterator.
  	   * @details	e.g. get every kmer which has both edges.
  	   * 			TODO: do this later.  It might be better to do it using filter_transfrom_iterator anyways.
  	   *
  	   * @tparam Predicate   operates on PredicateInputIterator and return true or false.
  	   */
  	  template <typename Iterator, typename PredicateInputIterator,
	  	  	  	typename Predicate = ::fsc::TruePredicate,
				typename Transform = ::bliss::kmer::transform::Identity<
					std::pair<typename ::std::iterator_traits<Iterator>::value_type,
					          typename ::std::iterator_traits<PredicateInputIterator>::value_type
							 >
  	  	  	  	  	>
  	  	  >
  	  using elemental_filter_transform_iterator =
  			  ::bliss::iterator::filter_transform_iterator<
			    ::bliss::iterator::zip_iterator<Iterator, PredicateInputIterator>,
				tuple_element_predicate<Predicate, 1, typename std::iterator_traits<Iterator>::value_type,
									   typename std::iterator_traits<PredicateInputIterator>::value_type
									  >,
				Transform
			>;

  	  /// make start elemental_filter_transform_iterator
  	  template<typename Iterator, typename PredicateInputIterator,
			 typename Predicate = ::fsc::TruePredicate,
				typename Transform = ::bliss::kmer::transform::Identity<
					std::pair<typename ::std::iterator_traits<Iterator>::value_type,
					          typename ::std::iterator_traits<PredicateInputIterator>::value_type
							 >
	  	  	  	  	>
  		>
  	  elemental_filter_transform_iterator<Iterator, PredicateInputIterator, Predicate, Transform>
  	  make_elemental_filter_transform_iterator(Iterator start, Iterator end,
  			  PredicateInputIterator pstart, PredicateInputIterator pend,
  			  Predicate const & f = Predicate(), Transform const & t = Transform()) {
		  using ZipIter = ::bliss::iterator::ZipIterator<Iterator, PredicateInputIterator>;

  		  return make_filter_transform_iterator(
  				  ZipIter(start, pstart),ZipIter(end, pend),
				  [&f](typename ::std::iterator_traits<ZipIter>::value_type const & x){
  			  	  	  return f(x.second);
  		  	  	  },
				  t);
  	  }


  	  /// make end elemental_filter_transform_iterator
  	  template<typename Iterator, typename PredicateInputIterator,
			 typename Predicate = ::fsc::TruePredicate,
				typename Transform = ::bliss::kmer::transform::Identity<
					std::pair<typename ::std::iterator_traits<Iterator>::value_type,
					          typename ::std::iterator_traits<PredicateInputIterator>::value_type
							 >
	  	  	  	  	>
  		>
  	  elemental_filter_transform_iterator<Iterator, PredicateInputIterator, Predicate, Transform>
  	  make_elemental_filter_transform_iterator(Iterator end,
  			  PredicateInputIterator pend,
  			  Predicate const & f = Predicate(), Transform const & t = Transform()) {
		  using ZipIter = ::bliss::iterator::ZipIterator<Iterator, PredicateInputIterator>;

  		  return make_filter_transform_iterator(
  				  ZipIter(end, pend),
				  [&f](typename ::std::iterator_traits<ZipIter>::value_type const & x){
  			  	  	  return f(x.second);
  		  	  	  },
				  t);
  	  }

  	  template <int I>
  	  struct tuple_element_extractor {
  		  template <typename tuple>
  		  typename std::tuple_element<I, tuple>::type operator()(tuple const & x) {
  			  return std::get<I>(x);
  		  }
  	  };


  	  /// make start elemental_filter_transform_iterator
  	  template<typename Iterator, typename PredicateInputIterator,
			 typename Predicate = ::fsc::TruePredicate
  		>
  	  elemental_filter_transform_iterator<Iterator, PredicateInputIterator, Predicate,
	  	  ::bliss::iterator::tuple_element_extractor<0> >
  	  make_elemental_filter_iterator(Iterator start, Iterator end,
  			  PredicateInputIterator pstart, PredicateInputIterator pend,
  			  Predicate const & f = Predicate()) {
		  using ZipIter = ::bliss::iterator::ZipIterator<Iterator, PredicateInputIterator>;

  		  return make_filter_transform_iterator(
  				  ZipIter(start, pstart),ZipIter(end, pend),
				  [&f](typename ::std::iterator_traits<ZipIter>::value_type const & x){
  			  	  	  return f(x.second);
  		  	  	  },
				  tuple_element_extractor<0>() );
  	  }


  	  /// make end elemental_filter_iterator
  	  template<typename Iterator, typename PredicateInputIterator,
			 typename Predicate = ::fsc::TruePredicate
  		>
  	  elemental_filter_transform_iterator<Iterator, PredicateInputIterator, Predicate,
	  	  ::bliss::iterator::tuple_element_extractor<0> >
  	  make_elemental_filter_iterator(Iterator end,
  			  PredicateInputIterator pend,
  			  Predicate const & f = Predicate()) {
		  using ZipIter = ::bliss::iterator::ZipIterator<Iterator, PredicateInputIterator>;

  		  return make_filter_transform_iterator(
  				  ZipIter(end, pend),
				  [&f](typename ::std::iterator_traits<ZipIter>::value_type const & x){
  			  	  	  return f(x.second);
  		  	  	  },
				  tuple_element_extractor<0>());
  	  }


  } // iterator
} // bliss
#endif /* COMPOUND_ITERATORS_HPP_ */
