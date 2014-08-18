/**
 * @file    fastq_iterator.hpp
 * @ingroup bliss::io
 * @author  Tony Pan <tpan@gatech.edu>
 * @date    Feb 5, 2014
 * @brief   Contains an iterator for traversing a FASTQ formatted DataBlock (see bliss::io::DataBlock) sequence record by record, and support classes.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add Licence
 *
 */


#ifndef FASTQ_ITERATOR_HPP_
#define FASTQ_ITERATOR_HPP_

// C includes
#include <cassert>
#include <cstdint>

// C++ STL includes
#include <iterator>
#include <algorithm>
#include <sstream>
#include <type_traits>

// own includes
#include <utils/logging.h>
#include <iterators/function_traits.hpp>
#include <partition/range.hpp>

namespace bliss
{

  namespace io
  {
    /**
     * @class     bliss::io::FASTQSequenceId
     * @brief     represents a fastq sequence's id, also used for id of the FASTQ file, and for position inside a fastq sequence..
     * @detail    this is set up as a union to allow easy serialization
     *            and parsing of the content.
     *            this keeps a 40 bit sequence ID, broken up into a 32 bit id and 8 bit significant bit id
     *                       an 8 bit file id
     *                       a 16 bit position within the sequence.
     *
     *            A separate FASTA version will have a different fields but keeps the same 64 bit total length.
     *
     *            file size is at the moment limited to 1TB (40 bits) in number of bytes.
     */
    union FASTQSequenceId
    {
        /// the concatenation of the id components as a single unsigned 64 bit field
        uint64_t composite;

        /// the id field components
        struct
        {
            /// sequence's id, lower 32 of 40 bits (potentially as offset in the containing file)
            uint32_t seq_id;
            /// sequence's id, upper 8 of 40 bits (potentially as offset in the containing file)
            uint8_t seq_msb;
            /// id of fastq file
            uint8_t file_id;
            /// offset within the read.  Default 0 refers to the whole sequence
            uint16_t pos;      // position value
        } components;
    };


    /// TODO: change this to use a single iterator, with ranges for the other fields?


    /**
     * @class     bliss::io::Sequence
     * @brief     represents a biological sequence, and provides iterators for traversing the sequence.
     *
     * @tparam Iterator   allows walking through the fastq data.
     * @tparam Alphabet   allows interpretation of the sequence
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.
     */
    template<typename Iterator, typename Alphabet, typename IdType=FASTQSequenceId>
    struct Sequence
    {
        /// Type of the sequence elements
        using ValueType = typename std::iterator_traits<Iterator>::value_type;
        /// The alphabet this sequence uses
        typedef Alphabet AlphabetType;
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;

        /// begin iterator for the sequence
        Iterator seqBegin;
        /// end iterator for the sequence.
        Iterator seqEnd;
        /// id of the sequence.
        IdType id;
    };

    /**
     * @class     bliss::io::SequenceWithQuality
     * @brief     extension of bliss::io::Sequence to include quality scores for each position.
     *
     * @tparam Iterator   allows walking through the fastq data.
     * @tparam Alphabet   allows interpretation of the sequence
     * @tparam Quality    data type for quality scores for each base
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.
     */
    template<typename Iterator, typename Alphabet, typename Quality, typename IdType=FASTQSequenceId>
    struct SequenceWithQuality : public Sequence<Iterator, Alphabet, IdType>
    {
        /// Type of the sequence elements
        using ValueType = typename std::iterator_traits<Iterator>::value_type;
        /// The alphabet this sequence uses
        typedef Alphabet AlphabetType;
        /// Iterator type for traversing the sequence.
        typedef Iterator IteratorType;
        /// the desired quality score's type, e.g. double
        typedef Quality ScoreType;

        /// begin iterator for the quality scores
        Iterator qualBegin;
        /// end iterator for the quality scores
        Iterator qualEnd;
    };



    /**
     * @class bliss::io::FASTQParser
     * @brief Functoid encapsulating functionality for increment and dereference in an iterator.
     * @details   The purpose of this class is to abstract the increment and dereference operations in
     *    an iterator that transforms a block of data into a collection of FASTQ records.  this allows us to reuse
     *    a transform iterator.
     *
     *    This class assumes that the iterator at start is pointing to the beginning of a FASTQ record already.
     *    Increment is done by walking through the data and recording the positions of the start and end of each
     *    of the 4 lines in a FASTQ record.  because we are always searching to the end of a record, successive
     *    calls to increment always are aligned correctly with record boundaries
     *
     *
     * @tparam Iterator   The underlying iterator to be traversed to generate a Sequence
     * @tparam Alphabet   allows interpretation of the sequence
     * @tparam Quality    data type for quality scores for each base.  defaults to void to mean No Quality Score.
     */
    template<typename Iterator, typename Alphabet, typename Quality = void>
    class FASTQParser
    {
        // internal state
      protected:

        /// Sequence type, conditionally set based on Quality template param to either normal Sequence or SequenceWithQuality.
        typedef typename std::conditional<std::is_void<Quality>::value,
              Sequence<Iterator, Alphabet>,
              SequenceWithQuality<Iterator, Alphabet, Quality> >::type      SeqType;
        /// The alphabet this sequence uses
        typedef Alphabet                                                    AlphabetType;
        /// the desired quality score's type, e.g. double.  conditionally defined only when needed
        typedef typename std::enable_if<!std::is_void<Quality>::value, Quality>::type
                                                                            QualityType;
        /// Range Type, set to use unsigned 64 bit for internal range representation
        typedef bliss::partition::range<size_t>                              RangeType;

        /// state of the functoid , which is the current sequence the iterator is pointing to.
        SeqType output;

      public:
        /// default constructor.
        FASTQParser() : output() {};

        /**
         * @brief function to populate the quality score. defined only when Quality type is not void
         * @param start   beginning of the quality score character sequence.
         * @param end     end of the quality score character sequence.
         */
        template <typename Q = Quality>
        inline typename std::enable_if<!std::is_void<Q>::value>::type
        populateQuality(const Iterator & start, const Iterator & end) {
          output.qual = start;
          output.qual_end = end;
        }

        /**
         * @brief function to populate the quality score. defined only when Quality type is void
         * @param start   beginning of the quality score character sequence.
         * @param end     end of the quality score character sequence.
         */
        template <typename Q = Quality>
        inline typename std::enable_if<std::is_void<Q>::value>::type
        populateQuality(const Iterator & start, const Iterator & end) {}

        /**
         * TODO:  change to use range coordinates for comparison instead of iterators.  measure performance...
         *
         * parses with an iterator, so as to have complete control over the increment.
         */
        /**
         * @brief Called by the containing iterator during iterator increment (++)
         * @details   Internally, this function parses the data by line, tracking the begin
         *            and end of the lines, and produces a single sequence record, stored in
         *            a member variable.
         *
         *            The source iterator is not modified.
         *
         * @param it            source iterator, pointing to data to traverse
         * @param end           end of the source iterator - not to go past.
         * @param coordinates   coordinates in units of source character types (e.g. char) from the beginning of the file.
         * @return              source iterator, advanced forward by the amount parsed.
         */
        Iterator increment(const Iterator &it, const Iterator &end, RangeType &coordinates)
        {
          // first make a copy of it so we can change it.
          Iterator iter = it;

          size_t count = 0;

          // trim leading \n
          while ((*iter == '\n') && (iter != end))
          {
            ++iter;  ++count;
          }

          // if the range consists of \n only, then we'd be at end/ terminate
          if (iter == end)
          {
            ERRORF("ERROR: nothing was parsed. %lu to %lu.\n", coordinates.start, coordinates.end);
            output = SeqType();
            return iter;
          }  // else we have some real data.

          //== generate the sequence instance.

          // store the "pointers"
          Iterator starts[4];
          Iterator ends[4];

          int line_num = 0;  // first line has num 0.

          //== do the first char
          // first char is already known (and not \n).  also, not at end.
          starts[line_num] = iter;
          output.id.composite = coordinates.start + count;
          bool isEOL = false;  // already trimmed all leading \n
          bool found = false;  // not yet found
          ++iter;  ++count;
          // now wasEOL is set since we've moved position.
          bool wasEOL = isEOL;

          //TODO:  optimize further.  even grep is faster (below takes about 50ms.  grep is at 30ms to search @


          //== walk through the data range until the end.
          while ((iter != end))
          {
            // check if current char is EOL
            isEOL = (*iter == '\n');  // slow

            // we have \nX or X\n.   skip this logic if 2 \n or 2 non-\n.
            if (isEOL != wasEOL)  // kind of slow
            {
              // a new EOL.  X\n case
              if (isEOL)
              {
                // found the end of a line
                ends[line_num] = iter;
                ++line_num;

                // once reached 4 lines, can stop
                if (line_num >= 4)
                {
                  found = true;

                  // increment past the \n just before break.
                  ++iter; ++count;
                  break;
                }
              }
              // first char of a new line.  \nX case
              else
              {
                starts[line_num] = iter;
              }

              wasEOL = isEOL;  // only toggle if isEOL != wasEOL.
            } // else we have \n\n or [^\n][^\n], continue.

            ++iter; ++count;  // kind of slow
          }


          //== if reached the end by this time
          if (iter == end) {

            // at end of iterator, but still missing a \n, then marked the last char as end of a line.
            if (!found) {
              ends[line_num] = iter;
              ++line_num;
            }

            // check to make sure we finished okay (found 4 lines) - if not, don't update the Sequence object.
            if (line_num < 4) {
              ERRORF("ERROR: parsing failed. lines %d, %lu to %lu. \n", line_num, coordinates.start, coordinates.end);
              std::stringstream ss;
              std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);

              std::copy(it, end, oit);
              ERRORF("  offending string is %s\n", ss.str().c_str());

              output = SeqType();
              return iter;
            }
          }

          assert(*(starts[0]) == '@');
          assert(*(starts[2]) == '+');

          //== now populate the output (instance variable)
          output.seq = starts[1];
          output.seq_end = ends[1];
          // this part is dependent on whether Quality template param is specified.
          populateQuality(starts[3], ends[3]);

          //== update the coordinates.
          coordinates.start += count;

          //== return the updated iter.
          return iter;
        }

        /// for dereferencing the parent iterator and return the result, which is stored in the local state variable: output.
        SeqType& dereference()
        {
          return output;
        }
    };

TODO: here.

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
    template<typename Parser, typename Iterator, typename base_traits = std::iterator_traits<Iterator> >
    class fastq_iterator :
        public std::iterator<
            typename std::conditional<
                  std::is_same<typename base_traits::iterator_category, std::random_access_iterator_tag>::value ||
                  std::is_same<typename base_traits::iterator_category, std::bidirectional_iterator_tag>::value,
                  std::forward_iterator_tag,
                  typename base_traits::iterator_category
                >::type,
            typename std::remove_reference<typename bliss::functional::function_traits<Parser>::return_type>::type,
            typename base_traits::difference_type
          >
    {
      public:


      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before fastq_iterator::Parser, etc are defined.
        typedef fastq_iterator<Parser, Iterator> type;

        typedef std::iterator_traits<type> traits;
        typedef typename Parser::RangeType RangeType;

        Iterator _curr;
        Iterator _next;
        Iterator _end;
        Iterator _temp;
        Parser _f;
        typename traits::value_type empty_output;
        RangeType range;     // the offsets in the base iterator.  if base iterator is pointer, then units are in bytes.

      public:

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
            _next = _f.increment(_temp, _end, range);    // _next is moved.
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

        inline typename traits::value_type operator*()
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
            _next = _f.increment(_temp, _end, range);
          }
          // else result was already computed and stored in the functor.

          return _f();
        }

        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator

        // forward iterator
        template<
            typename = typename std::enable_if<
                std::is_same<typename traits::iterator_category, std::forward_iterator_tag>::value ||
                std::is_same<typename traits::iterator_category, std::bidirectional_iterator_tag>::value ||
                std::is_same<typename traits::iterator_category, std::random_access_iterator_tag>::value>::type
               >
        explicit fastq_iterator()
            : _f(), _curr(), _next(), empty_output()
        {
        }


    };

  } // iterator
} // bliss
#endif /* FASTQ_ITERATOR_HPP_ */
