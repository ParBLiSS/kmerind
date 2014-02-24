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

#include <cassert>


namespace bliss
{

namespace iterator
{

  template<typename Iterator>
  struct fastq_sequence {
      // assume seq_name is the start, end is the beginning of next.  each line has a EOL at the end.
      union id_type {
          uint64_t id;
          struct {
            uint16_t pos;
            char file;
            char seq[5];
          } ids;
      };
      Iterator name;
      Iterator name_end;
      Iterator seq;
      Iterator seq_end;
      Iterator qual;
      Iterator qual_end;
  };

  /**
   * Iterator: base iterator from which to get the data for transform.
   */
  template<typename Iterator>
  struct fastq_parser {
      // internal state

      typedef fastq_sequence<Iterator> seq_type;

      seq_type output;
      Iterator _start;
      size_t   _global_offset;

      fastq_parser(Iterator const & start, size_t const & global_offset)
        : _start(start), _global_offset()
      {}

      /**
       * parses with an iterator, so as to have complete control over the increment.
       */

      size_t operator()(Iterator &iter, Iterator const& end) {
        // first initialize  (on windowed version, will need to have a separate way of initializing., perhaps with an overloaded operator.
        output = seq_type();

        size_t offset = _global_offset + (iter - _start);

        size_t dist = 0;

        if (iter == end)  // at end, terminate.
        {
          printf("not walked. %ld\n", dist);
          return dist;
        }
        // do some computation.

        // to simplify initial, get rid of leading \n
        while ((*iter == '\n') && (iter != end)) {
          ++dist;
          ++iter;
        }
        if (iter == end)  // if the range consists of \n only. at end, terminate
        {
          printf("walked (trimmed) %ld\n", dist);

          return dist;
        }

        // store the "pointers"
        Iterator starts[4];
        Iterator ends[4];

        // init current state to 1 (first char of line
        bool parsing = true;
        int line_num = 0;  // increment after the "end" of a line, so as to allow immediate termination.
        bool isEOL = false;
        bool wasEOL = true;   // beginning of first line, so prev char must be newline.


        //TODO:  optimize further.  even grep is faster (below takes about 50ms.  grep is at 30ms to search @
        // walk through the data range
        while (parsing && (iter != end))  // kind of slow
        {
          isEOL = (*iter == '\n');  // slow

          // early termination of iteration logic.  either 2 \n, or 2 non-\n
          if (isEOL != wasEOL)  // kind of slow
          {

            // otherwise we have \nX or X\n
            if (isEOL)  // a new eol.  X\n case
            {
              ends[line_num] = iter;

              ++line_num;
              if (line_num >= 4)
              {
  //              printf("got 4 lines.  %ld\n", dist);
                parsing = false;
              }
            }
            else  // first char.  \nX case
            {
              starts[line_num] = iter;
            }
            // iter == \n and newline Char.  keep going.  or iter != \n and !newline char.  in the middle.  keep going.

          }

          wasEOL = isEOL;
          ++dist;
          ++iter;   // kind of slow
        }
        //printf("walked total %ld\n", dist);

        // check to make sure we finished okay  - if not, don't update teh fastq_sequence object.
        assert(*(starts[0]) == '@');
        assert(*(starts[2]) == '+');


        if ((iter == end) && (line_num < 4))
          return dist;


        // now populate the output
        output.id = offset;
        output.name = starts[0];
        output.name_end = ends[0];
        output.seq = starts[1];
        output.seq_end = ends[1];
        output.qual = starts[3];
        output.qual_end = ends[3];

        return dist;
      }

      seq_type& operator()()
      {
        return output;
      }
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
                                                     typename std::iterator_traits<Iterator>::iterator_category
                                                    >::type,
                           typename std::remove_reference<typename bliss::functional::function_traits<Parser>::return_type>::type,
                           typename std::iterator_traits<Iterator>::difference_type,
                           typename std::add_pointer<typename std::remove_reference<typename bliss::functional::function_traits<Parser>::return_type>::type>::type,
                           typename std::add_rvalue_reference<typename std::remove_reference<typename bliss::functional::function_traits<Parser>::return_type>::type>::type
                           >
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before fastq_iterator::Parser, etc are defined.
        typedef std::iterator_traits<Iterator>                                                           base_traits;
        typedef bliss::functional::function_traits<Parser> functor_traits;

        Iterator                                                                                          _curr;
        Iterator                                                                                          _next;
        Iterator                                                                                          _end;
        Parser                                                                                             _f;
        typename std::remove_reference<typename functor_traits::return_type>::type                         empty_output;

      public:
        typedef fastq_iterator< Parser, Iterator >                                                          type;
        typedef std::iterator_traits<type>                                                                  traits;

//TODO
        typedef typename std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                                      std::random_access_iterator_tag>::value ||
                                          std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                                      std::bidirectional_iterator_tag>::value,
                                          std::forward_iterator_tag,
                                          typename std::iterator_traits<Iterator>::iterator_category
                                          >::type                                                     iterator_category;
        typedef typename std::remove_reference<typename functor_traits::return_type>::type                  value_type;
        typedef typename base_traits::difference_type                                                       difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type                                        reference_type;
        typedef typename std::add_pointer<value_type>::type                                                 pointer_type;

        typedef typename base_traits::value_type                                                            base_value_type;



        // class specific constructor
        explicit
        fastq_iterator(const Parser & f, const Iterator& curr, const Iterator& end)
          : _curr(curr), _next(curr), _end(end), _f(f), empty_output() {};
        explicit
        fastq_iterator(const Parser & f, const Iterator& end)
          : _curr(end), _next(end), _end(end), _f(f), empty_output() {};


        // for all classes
        // no EXPLICIT so can copy assign.
        fastq_iterator(const type& Other) : _curr(Other._curr), _next(Other._next), _end(Other._end), _f(Other._f), empty_output() {};

        type& operator=(const type& Other) {
          _f = Other._f;
          _curr = Other._curr;
          _next = Other._next;
          _end = Other._end;
          return *this;
        }

        type& operator++() {

          // special case, at end of iter.
          if (_curr == _end)
          {
            // no movement
            return *this;
          }

          // property:  to set _next to the appropriate position, need to parse.

          // ++ parses so that _curr moves.  _next is pushed forward.  if _next has already been moved (via *), then just move curr there.

          if (_curr == _next)  // at beginning.  so need to parse, and right away
          {
            // walk the base iterator until function is done with construction or at end.
            // doing the construction here because we have to walk anyways.
            _f(_next, _end);    // _next is moved.
          }
          // else next has already been moved (by * operator)

          _curr = _next;

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
          return _curr;
        }
        const Iterator& getBaseIterator() const {
          return _curr;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        { return _curr == rhs._curr; }

        inline bool operator!=(const type& rhs)
        { return _curr != rhs._curr; }

        inline value_type operator*() {
          // boundary case: at end of iterators
          if (_curr == _end) {
            return empty_output;
          }

          // special case, have not yet parsed the content
          if (_curr == _next)
          {
            // after parse, _next has been moved but curr stays where it is.
            _f(_next, _end);
          }
          // else result was already computed and stored in the functor.

          return _f();
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator


        // forward iterator
        template<typename = std::enable_if<std::is_same<iterator_category, std::forward_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::bidirectional_iterator_tag>::value ||
                                           std::is_same<iterator_category, std::random_access_iterator_tag>::value > >
        explicit
        fastq_iterator() : _f(), _curr(), _next(), empty_output() {};


    };


} // iterator
} // bliss
#endif /* FASTQ_ITERATOR_HPP_ */
