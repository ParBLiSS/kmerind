/**
 * @file    SequencesIterator.hpp
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


#ifndef SequencesIterator_HPP_
#define SequencesIterator_HPP_

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
     * @note  NOT THREAD SAFE.
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

      public:
        /// default constructor.
        FASTQParser() {};

        /**
         * @brief function to populate the quality score. defined only when Quality type is not void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Q = Quality>
        typename std::enable_if<!std::is_void<Q>::value>::type
        populateQuality(const Iterator & start, const Iterator & end, SeqType& output) {
          output.qualStart = start;
          output.qualEnd = end;
        }

        /**
         * @brief function to populate the quality score. defined only when Quality type is void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Q = Quality>
        typename std::enable_if<std::is_void<Q>::value>::type
        populateQuality(const Iterator & start, const Iterator & end, SeqType& output) {}


        /**
         * @brief function to populate the sequence iterators.
         * @param[in] start   beginning of the sequence.
         * @param[in] end     end of the sequence.
         * @param[out] output  The output sequence record type to modify
         */
        void populateSequence(const Iterator & start, const Iterator & end, SeqType& output) {
          output.seqStart = start;
          output.seqEnd = end;
        }

        /**
         * @brief Called by the containing iterator during iterator increment (++)
         * @details   Internally, this function parses the data by line, tracking the begin
         *            and end of the lines, and produces a single sequence record, stored in
         *            a member variable.
         *
         *            The source iterator is not modified.  A copy of the iterator is modified.
         *
         *            Iterator can be a forward iterator or better.
         *
         * @param[in/out] iter          source iterator, pointing to data to traverse
         * @param[in]     end           end of the source iterator - not to go past.
         * @param[in/out] coordinates   coordinates in units of source character types (e.g. char) from the beginning of the file.
         * @param[in/out] output        updated sequence type, values updated.
         */
        void increment(Iterator &iter, const Iterator &end, size_t &offset, SeqType& output)
        {
          //== first make a copy of iter so we can change it.
          Iterator orig_iter = iter;
          size_t orig_offset = offset;

          //== trim leading \n
          while ((*iter == '\n') && (iter != end))
          {
            ++iter;  ++offset;
          }

          // if at this point we are at end, then this is an incomplete record.
          if (iter == end)
          {
            ERRORF("ERROR: nothing was parsed. %lu to %lu.\n", orig_offset, offset);
            output = SeqType();
            return;
          }
          // else we have some real data.

          //== generate the sequence instance.

          // store the "pointers"
          Iterator starts = end;
          int line_num = 0;  // first line has num 0.

          //== do the first char
          // first char is already known (and not \n), so isEOL is false.  also, not at end.
//          starts[line_num] = iter;
          output.id.composite = offset;
          bool isEOL = false;  // already trimmed all leading \n

          // move to next char, and set wasEOL.
          ++iter;  ++offset;
          bool wasEOL = isEOL;

          //TODO:  optimize further.  even grep is faster (below takes about 50ms.  grep is at 30ms to search @


          //== walk through the data range until the end.
          while ((iter != end))
          {
            //== check if current char is EOL
            isEOL = (*iter == '\n');  // slow

            //== check if this is start of a line or end of a line
            //  where we have \nX or X\n.   skip this logic if 2 \n or 2 non-\n.
            if (isEOL != wasEOL)  // kind of slow
            {
              // a new EOL.  X\n case
              if (isEOL)
              {
                // found the end of a line
                if (line_num == 1) populateSequence(starts, iter, output);
                else if (line_num == 3) populateQuality(starts, iter, output);
                ++line_num;
                // reset start
                starts = end;

                // once reached 4 lines, can stop the loop
                if (line_num >= 4)
                {
                  // increment past the \n just before break.
                  ++iter; ++offset;
                  break;
                }
              }
              // first char of a new line.  \nX case
              else
              {
                starts = iter;  // regardless of the line number, save it.
              }

              wasEOL = isEOL;  // only toggle if isEOL != wasEOL.
            } // else we have \n\n or [^\n][^\n], continue.

            ++iter; ++offset;  // kind of slow
          }

          //== at this point, got 4 lines (iter == end or not), or don't have 4 lines (iter has to be  at end)
          // if iter reached end with less than 4 lines, then mark the last char as the end of a line
          if (line_num < 4) {
            if (line_num == 1) populateSequence(starts, iter, output);
            else if (line_num == 3) populateQuality(starts, iter, output);
            ++line_num;
          }

          TODO:

          //== after completing the last line, check to see if we have the needed number of lines.
          //== if no quality, need 2 lines, and if quality needed, need 4 lines
          if (std::is_void<Quality>::value) { // not requiring quality scores, so need at least 2 lines.
            if (line_num < 2) {
              // not enough lines, so raise error and change output

              output = SeqType();
            }
            // else have enough lines, everything should be updated already.


          } else {                            // requiring quality scores, so need 4 lines

            if (line_num < 4) {
              // not enough lines, so raise error and change output
              output = SeqType();
            }
            // else have enough lines, everything should be updated already.

          }



          //== if reached the end by this time
          if (iter == end) {

            // at end of iterator, but still missing a \n, then marked the last char as end of a line.
            if (!found) {
            }

            // check to make sure we finished okay (found 4 lines) - if not, don't update the Sequence object.
            if (line_num < 4) {
              ERRORF("ERROR: parsing failed. lines %d, %lu to %lu. \n", line_num, orig_offset, offset);

              std::stringstream ss;
              std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
              std::copy(orig_iter, end, oit);
              ERRORF("  offending string is %s\n", ss.str().c_str());

              output = SeqType();
              return;
            }
          }

          // else, we have 4 lines.

          //== now populate the output (instance variable)
          output.seqStart = starts[1];
          output.seqEnd = ends[1];
          // this part is dependent on whether Quality template param is specified.
          populateQuality(starts[3], ends[3], output);

          // coordinates already updated.

        }
    };

    /**
     * @class bliss::io::SequencesIterator
     * @brief Iterator for parsing and traversing a block of data to access individual FASTQ sequence records.
     * @details The iterator is initialized with start and end iterators pointing to the source data.
     *    The source data may be of primitive type such as unsigned char.  The input iterators should be aligned
     *    to record boundaries.
     *
     *    The iterator uses FASTQParser to walk through the data one sequence record at a time
     *
     *    The iterator has the following characteristics.
     *      transform  (changes data type)
     *      delimiter  (trivial search for the end of the "record" / beginning of next record.)
     *      finite  (needs to know base iterator's end.)
     *
     *    Design considerations:
     *      1. repeated calls to operator*() returns the same output, WITHOUT reparsing the underlying data
     *      2. repeated calls to operator++() moves the current pointer forward each time
     *      3. comparison of 2 SequencesIterators compares the underlying input iterators at the start positions of the current record
     *      4. iterators have states:  before, during, after.
     *          before:  operator* points to the first output, no need to call ++ first.  so this is same case as "during".
     *          during:  operator* gets the most current output that's a result of calling ++, ++ moves pointer forward
     *          after:   operator* returns empty result, ++ is same as no-op.
     *     Based on these, we need to
     *      1. cache the output
     *      2. keep a current iterator (for comparison between SequencesIterators)
     *      3. keep a next iterator (for parsing)
     *      3. have constructor parse the first entry, so that dereference functions properly
     *      4. check for end condition (next == end) and return if at end.
     *
     *      The functoid then only has 1 simple task:  to increment, as it does not need to store locally.
     *
     * @note this iterator is not a forward or bidirectional iterator, as traversing in reverse is not implemented.
     *       NOT SPECIFIC TO FASTQ FORMAT
     *
     * @tparam Parser     Functoid type to supply specific logic for increment and dereference functionalities
     * @tparam Iterator   Type of the input iterator to parse/traverse
     * @tparam base_traits  default to std::iterator_traits of Iterator, defined for convenience.
     */
    template<typename Parser, typename Iterator >
    class SequencesIterator :
        public std::iterator<
            typename std::conditional<
                  std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::random_access_iterator_tag>::value ||
                  std::is_same<typename std::iterator_traits<Iterator>::iterator_category, std::bidirectional_iterator_tag>::value,
                  std::forward_iterator_tag,
                  typename std::iterator_traits<Iterator>::iterator_category
                >::type,
            typename std::remove_reference<typename bliss::functional::function_traits<Parser>::return_type>::type,
            typename std::iterator_traits<Iterator>::difference_type
          >
    {
      public:

        /// type of data elements in the input iterator
        typedef typename std::iterator_traits<Iterator>::value_type BaseValueType;

        /// type of data elements in the input iterator
        typedef typename std::iterator_traits<SequencesIterator<Parser, Iterator> >::value_type ValueType;

      protected:
        /// type of SequencesIterator class
        typedef SequencesIterator<Parser, Iterator> type;

        /// type traits of SequencesIterator class.
        typedef std::iterator_traits<type> traits;

        /// cache of the current result.
        ValueType seq;

        /// iterator at the current starting point in the input data from where the current FASTQ record was parsed.
        Iterator _curr;

        /// iterator at the next starting point in the input data from where the next FASTQ record should be parsed.
        Iterator _next;

        /// iterator pointing to the end of the input data, not to go beyond.
        Iterator _end;

        /// instance of the parser functoid, used during record parsing.
        Parser parser;

        /**
         * @brief  the offset from the start of the file.  This is used as sequenceId.
         * @details  for each file with file id (fid) between 0 and 255, the seq id starts at (fid << 40).
         *          if base iterator is pointer, then units are in bytes.
         */
        size_t seqId;


      public:
        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         * @details   at construction, the internal _curr, _next iterators are set to the input start iterator, and _end is set to the input end.
         *            then immediately the first record is parsed, with output populated with the first entry and next pointer advanced.
         *
         *
         * @param f       parser functoid for parsing the data.
         * @param start   beginning of the data to be parsed
         * @param end     end of the data to be parsed.
         * @param _range  the Range associated with the start and end of the source data.  coordinates relative to the file being processed.
         */
        SequencesIterator(const Parser & f, const Iterator& start,
                       const Iterator& end, const size_t &_offset)
            : seq(), _curr(start), _next(start), _end(end), parser(f), seqId(_offset)
        {
          // parse the first entry, if start != end.
          // effect: _next incremented to next parse start point.
          //         seqId updated to next start point,
          //         output updated with either empty (if _next at end)  or the new output
          if (_next != _end) parser.increment(_next, _end, seqId, seq);;
        }

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing SequencesIterator.
         * @param end     end of the data to be parsed.
         */
        SequencesIterator(const Iterator& end)
            : seq(), _curr(end), _next(end), _end(end), parser(), seqId(std::numeric_limits<size_t>::max())
        {
        }

        /**
         * @brief default copy constructor
         * @param Other   The SequencesIterator to copy from
         */
        SequencesIterator(const type& Other)
            : seq(Other.seq), _curr(Other._curr), _next(Other._next), _end(Other._end),
              parser(Other.parser),  seqId(Other.seqId)
        {}

        /**
         * @brief  default copy assignment operator
         * @param Other   The SequencesIterator to copy from
         * @return  modified copy of self.
         */
        type& operator =(const type& Other)
        {
          seq = Other.seq;
          _curr = Other._curr;
          _next = Other._next;
          _end = Other._end;
          parser = Other.parser;
          seqId = Other.seqId;
          return *this;
        }

        /// default constructor deleted.
        SequencesIterator() = delete;


        /**
         * @brief   iterator's increment operation
         * @details to increment, the class uses the parser functoid to parse and traverse the input iterator
         *          At the end of the traversal, a FASTQ sequence record has been made and is kept in the parser functoid
         *          as state variable, and the input iterator has advanced to where parsing would start for the next record.
         *          The corresponding starting of the current record is saved as _curr iterator.
         *
         *          conversely, to get the starting position of the next record, we need to parse the current record.
         *
         * @note    for end iterator, because _curr and _end are both pointing to the input's end, "increment" does not do anything.
         * @return  self, updated with internal position.
         */
        type &operator ++()
        {
          //== at end of the data block, can't parse any further, no movement.  (also true for end iterator)
          if (_next == _end)
          {
            //== check if reaching _end for the first time
            if (_curr != _end) {
              // advance to end
              _curr = _next;
              // and set output to empty.
              seq = ValueType();
            }
          } else {
            //== else we can parse, so do it.

            //== move the current forward
            _curr = _next;

            //== then try parsing.  this advances _next and seqId, ready for next ++ call, and updates output cache.
            parser.increment(_next, _end, seqId, seq);
          }

          //== states updated.  return self.
          return *this;
        }

        type operator ++(int)
        {
          type output(*this);
          return ++output;
        }

        // accessor functions for internal state.
        Parser& getParser()
        {
          return parser;
        }
        const Parser& getParser() const
        {
          return parser;
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
        bool operator ==(const type& rhs) const
        {
          return _curr == rhs._curr;
        }

        bool operator !=(const type& rhs) const
        {
          return _curr != rhs._curr;
        }

        ValueType &operator *() const
        {
          return seq;
        }

        ValueType *operator ->() const {
          return &seq;
        }


    };

  } // iterator
} // bliss
#endif /* SequencesIterator_HPP_ */
