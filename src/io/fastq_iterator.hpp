
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

// C includes
#include <cassert>
#include <cstdint>

// C++ STL includes
#include <iterator>
#include <type_traits>

// own includes
#include <utils/logging.h>
#include <iterators/function_traits.hpp>

namespace bliss
{

  namespace io
  {
    /**
     * @class     bliss::io::fastq_sequence_id
     * @brief     data structure to represent a fastq sequence id.
     * @detail    this is set up as a union to allow easy serialization
     *            and parsing of the content.
     *            this keeps a 40 bit sequence ID, a file id, and a position within the file.
     *
     *            A separate FASTA version will have a different partitioning.
     */
    union fastq_sequence_id
    {
        uint64_t composite;
        struct
        {
            uint32_t seq_id;
            uint8_t seq_msb;   // seq_id and seq_msb: 40 bits allow for 1 trillion entries.
            uint8_t file_id;   // file id
            uint16_t pos;      // position value
        } components;
    };



    /**
     * @class     bliss::io::fastq_sequence
     * @brief     represents a fastq read.
     * @details   Iterator allows walking through the fastq data.
     *            Alphabet allows interpretation of the sequence
     *
     *
     *            a separate one would be used for FASTA sequence,
     *            where there is no quality score.
     */
    template<typename Iterator, typename Alphabet>
    struct fastq_sequence
    {
        typedef typename std::remove_pointer<
            typename std::remove_reference<Iterator>::type>::type ValueType;
        typedef Alphabet AlphabetType;
        typedef Iterator IteratorType;

        Iterator name;
        Iterator name_end;
        Iterator seq;
        Iterator seq_end;

        fastq_sequence_id id;

        static void allocCopy(const fastq_sequence<Iterator, Alphabet>& src, fastq_sequence<Iterator, Alphabet>& dest) {
          dest.id.composite = src.id.composite;

          size_t length = src.name_end - src.name;
          dest.name = new ValueType[length];
          memcpy(dest.name, src.name, length);
          dest.name_end = dest.name + length;

          length = src.seq_end - src.seq;
          dest.seq = new ValueType[length];
          memcpy(dest.seq, src.seq, length);
          dest.seq_end = dest.seq + length;

        }
        static void deleteCopy(fastq_sequence<Iterator, Alphabet>& dest) {
          delete [] dest.name;
          delete [] dest.seq;
        }
    };

    /**
     * @class     bliss::io::fastq_sequence
     * @brief     represents a fastq read without quality score.
     * @details   Iterator allows walking through the fastq data.
     *            Alphabet allows interpretation of the sequence
     *            Quality allows interpretation of the quality score.
     *
     *            has to be derived class because can't conditionally define member variables.
     *
     *            a separate one would be used for FASTA sequence,
     *            where there is no quality score.
     */
    template<typename Iterator, typename Alphabet, typename Quality>
    struct fastq_sequence_quality : public fastq_sequence<Iterator, Alphabet>
    {
        typedef fastq_sequence<Iterator, Alphabet> base_class_t;
        typedef typename base_class_t::ValueType ValueType;
        typedef Alphabet AlphabetType;
        typedef Iterator IteratorType;
        typedef Quality ScoreType;

        Iterator qual;
        Iterator qual_end;

        static void allocCopy(const fastq_sequence_quality<Iterator, Alphabet, Quality>& src, fastq_sequence_quality<Iterator, Alphabet, Quality>& dest) {
          base_class_t::allocCopy(src, dest);

          size_t length = src.qual_end - src.qual;
          dest.qual = new ValueType[length];
          memcpy(dest.qual, src.qual, length);
          dest.qual_end = dest.qual + length;

        }
        static void deleteCopy(fastq_sequence_quality<Iterator, Alphabet, Quality>& dest) {
          base_class_t::deleteCopy(dest);
          delete [] dest.qual;
        }
    };



    /**
     * Automatically
     */
    template<typename Iterator, typename Alphabet, typename Quality = void>
    struct fastq_parser
    {
        // internal state

        typedef typename std::conditional<std::is_same<Quality, void>::value,
              fastq_sequence<Iterator, Alphabet>,
              fastq_sequence_quality<Iterator, Alphabet, Quality> >::type   SeqType;
        typedef Alphabet                                                    AlphabetType;
        typedef std::enable_if<!std::is_same<Quality, void>::value, Quality>
                                                                            QualityType;

        SeqType output;

        fastq_parser() : output()
        {
        }

        template <typename Q = Quality>
        typename std::enable_if<!std::is_same<Q, void>::value>::type
        populateQuality(const Iterator & start, const Iterator & end) {
          output.qual = start;
          output.qual_end = end;
        }
        template <typename Q = Quality>
        typename std::enable_if<std::is_same<Q, void>::value>::type
        populateQuality(const Iterator & start, const Iterator & end) {
        }

        /**
         * parses with an iterator, so as to have complete control over the increment.
         */
        Iterator operator()(const Iterator &it, const Iterator &end, const RangeType &coordinates)
        {
          // first initialize  (on windowed version, will need to have a separate
          // way of initializing., perhaps with an overloaded operator.
          Iterator iter = it;
          output = SeqType();

          // trim leading \n
          while ((*iter == '\n') && (iter != end))
          {
            ++iter;
          }

          if (iter == end) // if the range consists of \n only. at end, terminate
          {
            printf("ERROR: nothing was parsed. %lu to %lu.  iter %p, end %p\n", coordinates.start, coordinates.end, iter, end);
            return iter;
          }

          // store the "pointers"
          Iterator starts[4];
          Iterator ends[4];

          int line_num = 0;  // first line has num 0.
          bool isEOL;

          // increment after the "end" of a line, so as to allow immediate termination

          // first char is already known (and not \n).  also, not at end.
          starts[line_num] = iter;
          ++iter;
          bool parsing = (iter != end);
          bool wasEOL = false;


          //TODO:  optimize further.  even grep is faster (below takes about 50ms.
          // grep is at 30ms to search @
          // walk through the data range
          while (parsing)
          {
            isEOL = (*iter == '\n');  // slow

            // early termination of iteration logic.  either 2 \n, or 2 non-\n
            if (isEOL != wasEOL)  // kind of slow
            {

              // we have \nX or X\n
              if (isEOL)  // a new eol.  X\n case
              {
                ends[line_num] = iter;

                ++line_num;
                if (line_num >= 4)
                {
                  parsing = false;
                }
              }
              else  // first char.  \nX case
              {
                starts[line_num] = iter;
              }
              // iter == \n and newline Char.  keep going.
              // or iter != \n and !newline char.  in the middle.  keep going.

              wasEOL = isEOL;  // only toggle if isEOL != wasEOL.
            }

            ++iter;   // kind of slow
            if (iter == end) {
              parsing = false;
              ends[line_num] = iter;
              ++line_num;
            }
          }

          // check to make sure we finished okay  - if not, don't update the
          // fastq_sequence object.
          if ((iter == end) && (line_num < 4)) {
            printf("ERROR: parsing failed. lines %d, %lu to %lu.  start %p, curr %p, end %p\n", line_num, coordinates.start, coordinates.end, it, iter, end);
            unsigned char* offending = new unsigned char[(end - it) + 1];
            memcpy(offending, it, (end - it));
            memset(offending + (end - it), 0, sizeof(unsigned char));
            printf("  offending string is %s\n", offending);
            delete [] offending;
            return iter;
          }

          assert(*(starts[0]) == '@');
          assert(*(starts[2]) == '+');

          // now populate the output
          output.id.composite = coordinates.start + (starts[0] - it);
          output.name = starts[0];
          output.name_end = ends[0];
          output.seq = starts[1];
          output.seq_end = ends[1];

          populateQuality(starts[3], ends[3]);
          return iter;
        }

        SeqType& operator()()
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
     *
     *
     */
    template<typename Parser, typename Iterator>
    class fastq_iterator :
        public std::iterator<
            typename std::conditional<
                std::is_same<
                    typename std::iterator_traits<Iterator>::iterator_category,
                    std::random_access_iterator_tag>::value
                || std::is_same<
                    typename std::iterator_traits<Iterator>::iterator_category,
                    std::bidirectional_iterator_tag>::value,
                std::forward_iterator_tag,
                typename std::iterator_traits<Iterator>::iterator_category>::type,
            typename std::remove_reference<
                typename bliss::functional::function_traits<Parser>::return_type>::type,
            typename std::iterator_traits<Iterator>::difference_type>
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before fastq_iterator::Parser, etc are defined.
        typedef std::iterator_traits<Iterator> base_traits;
        typedef bliss::functional::function_traits<Parser> functor_traits;

        Iterator _curr;
        Iterator _next;
        Iterator _end;
        Iterator _temp;
        Parser _f;
        typename std::remove_reference<typename functor_traits::return_type>::type empty_output;
        RangeType range;

      public:
        typedef fastq_iterator<Parser, Iterator> type;
        typedef std::iterator_traits<type> traits;

//TODO
        typedef typename std::conditional<
            std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::random_access_iterator_tag>::value
            || std::is_same<
                typename std::iterator_traits<Iterator>::iterator_category,
                std::bidirectional_iterator_tag>::value,
            std::forward_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::type iterator_category;
        typedef typename std::remove_reference<
            typename functor_traits::return_type>::type ValueType;
        typedef typename base_traits::difference_type difference_type;
        typedef typename std::add_rvalue_reference<ValueType>::type reference_type;
        typedef typename std::add_pointer<ValueType>::type pointer_type;

        typedef typename base_traits::value_type BaseValueType;

        // class specific constructor
        fastq_iterator(const Parser & f, const Iterator& curr,
                       const Iterator& end, const RangeType &_range)
            : _curr(curr), _next(curr), _end(end), _f(f), empty_output(), range(_range)
        {
        }

        fastq_iterator(const Parser & f, const Iterator& end, const RangeType &_range)
            : _curr(end), _next(end), _end(end), _f(f), empty_output(), range(_range.end, _range.end)
        {
        }

        // for all classes
        fastq_iterator(const type& Other)
            : _curr(Other._curr), _next(Other._next), _end(Other._end),
              _f(Other._f), empty_output(), range(Other.range)
        {
        }

        type& operator=(const type& Other)
        {
          _f = Other._f;
          _curr = Other._curr;
          _next = Other._next;
          _end = Other._end;
          range = Other.range;
          return *this;
        }

        type& operator++()
        {

          // special case, at end of iter.
          if (_curr == _end)
          {
            // no movement
            return *this;
          }

          // property:  to set _next to the appropriate position, need to parse.

          // ++ parses so that _curr moves.  _next is pushed forward.  if _next
          //has already been moved (via *), then just move curr there.

          if (_curr == _next) // at beginning.  so need to parse, and right away
          {
            // walk the base iterator until function is done with construction or
            // at end.
            // doing the construction here because we have to walk anyways.
            _temp = _next;
            _next = _f(_next, _end, range);    // _next is moved.
            range.start += (_next - _temp) / sizeof(BaseValueType);
          }
          // else next has already been moved (by * operator)

          _curr = _next;

          return *this;
        }

        type operator++(int)
        {
          type output(*this);
          return ++output;
        }

        // accessor functions for internal state.
        Parser& getParser()
        {
          return _f;
        }
        const Parser& getParser() const
        {
          return _f;
        }
        Iterator& getBaseIterator()
        {
          return _curr;
        }
        const Iterator& getBaseIterator() const
        {
          return _curr;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        {
          return _curr == rhs._curr;
        }

        inline bool operator!=(const type& rhs)
        {
          return _curr != rhs._curr;
        }

        inline ValueType operator*()
        {
          // boundary case: at end of iterators
          if (_curr == _end)
          {
            return empty_output;
          }

          // special case, have not yet parsed the content
          if (_curr == _next)
          {
            // after parse, _next has been moved but curr stays where it is.
            _temp = _next;
            _next = _f(_next, _end, range);
            range.start += (_next - _temp) / sizeof(BaseValueType);
          }
          // else result was already computed and stored in the functor.

          return _f();
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator

        // forward iterator
        template<
            typename = std::enable_if<
                std::is_same<iterator_category, std::forward_iterator_tag>::value || std::is_same<
                    iterator_category, std::bidirectional_iterator_tag>::value
                || std::is_same<iterator_category,
                    std::random_access_iterator_tag>::value> >
        explicit fastq_iterator()
            : _f(), _curr(), _next(), empty_output()
        {
        }


    };

  } // iterator
} // bliss
#endif /* FASTQ_ITERATOR_HPP_ */
