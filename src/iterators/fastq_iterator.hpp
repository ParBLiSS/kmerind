/**
 * fastq_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 *
 * special iterator to get fastq reads.
 *
 *
 *  TODO:
 *  commenting.
 */

#ifndef FASTQ_ITERATOR_HPP_
#define FASTQ_ITERATOR_HPP_

#include <iterator>
#include "iterators/function_traits.hpp"

namespace bliss
{

namespace iterator
{

  struct fastq_sequence {
      // assume seq_name is the start, end is the beginning of next.  each line has a EOL at the end.
      char* seq_name;
      char* seq;
      char* qual_name;
      char* qual;
      char* end;
  };


  /**
   * has the following characteristics:
   * transform  (changes data type)
   * delimiter  (trivial search for the end of the "record" / beginning of next record.)   != random access iterator.
   * finite  (needs to know base iterator's end.)
   *
   * implement forward only for now.  searching backwards is not easy.
   */
  template<typename Parser,
           typename Iterator
           >
  class fastq_iterator
    : public std::iterator<typename std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                                                  std::random_access_iterator_tag>::value ||
                                                    std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                                          std::bidirectional_iterator_tag>::value,
                                                          std::forward_iterator_tag,
                                                    typename std::iterator_traits<Iterator>::iterator_category>::type,
                           typename std::remove_reference<typename bliss::functional::function_traits<Parser, typename std::iterator_traits<Iterator>::value_type>::return_type>::type,
                           typename std::iterator_traits<Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename bliss::functional::function_traits<Parser, typename std::iterator_traits<Iterator>::value_type>::return_type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename bliss::functional::function_traits<Parser, typename std::iterator_traits<Iterator>::value_type>::return_type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before fastq_iterator::Parser, etc are defined.
        typedef std::iterator_traits<Iterator>                                                         base_traits;
        typedef bliss::functional::function_traits<Parser, typename std::iterator_traits<Iterator>::value_type>                                    functor_traits;

        Iterator                                                                                          _base;
        Parser                                                                                             _f;

      public:
        typedef fastq_iterator< Parser, Iterator >                                                          type;
        typedef std::iterator_traits<type>                                                                  traits;

//TODO
        typedef typename base_traits::iterator_category                                                     iterator_category;
        typedef typename std::remove_reference<typename functor_traits::return_type>::type                  value_type;
        typedef typename base_traits::difference_type                                                       difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type                                        reference_type;
        typedef typename std::add_pointer<value_type>::type                                                 pointer_type;

        typedef typename base_traits::value_type                                                            base_value_type;


        // class specific constructor
        explicit
        fastq_iterator(const Iterator& base_iter, const Parser & f)
          : _base(base_iter), _f(f) {};


        // for all classes
        // no EXPLICIT so can copy assign.
        fastq_iterator(const type& Other) : _base(Other._base), _f(Other._f) {};

        type& operator=(const type& Other) {
          _f = Other._f;
          _base = Other._base;
          return *this;
        }

        type& operator++() {

          // walk the base and construct the output.

          ++_base;
          return *this;
        }

        type operator++(int) {
          type output(*this);
          return ++output;
        }


        // accessor functions for internal state.
        Parser& getParser() {
          return _f;
        }
        const Parser& getParser() const {
          return _f;
        }
        Iterator& getBaseIterator() {
          return _base;
        }
        const Iterator& getBaseIterator() const {
          return _base;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        { return _base == rhs._base; }

        inline bool operator!=(const type& rhs)
        { return _base != rhs._base; }

        inline value_type operator*() {
          // return the built or build here?


          return _f(*_base);
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator


        // forward iterator
        template<typename = std::enable_if<std::is_same<iterator_category, std::forward_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::random_access_iterator_tag>::value > >
        explicit
        fastq_iterator() : _f(), _base() {};


    };


} // iterator
} // bliss
#endif /* FASTQ_ITERATOR_HPP_ */
