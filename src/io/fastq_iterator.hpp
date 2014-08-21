/**
 * @file    SequencesIterator.hpp
 * @ingroup bliss::io
 * @author  Tony Pan <tpan@gatech.edu>
 * @date    Feb 5, 2014
 * @brief   Contains an iterator for traversing a sequence data, record by record, a FASTQ format specific parser for use by the iterator,
 *          Sequence Id datatype definition, and 2 sequence types (with and without quality score iterators)
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
#include <io/io_exception.hpp>
#include <utils/logging.h>
#include <iterators/function_traits.hpp>
#include <partition/range.hpp>

namespace bliss
{

  namespace io
  {
    /**
     * @class     bliss::io::FASTQSequenceId
     * @brief     represents a fastq sequence's id, also used for id of the FASTQ file, and for position inside a FASTQ sequence..
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


    /**
     * @class     bliss::io::Sequence
     * @brief     represents a biological sequence, and provides iterators for traversing the sequence.
     * @details   internally we store 2 iterators, start and end, instead of start + length.
     *            reason is that using length is best for random access iterators and we don't have guaranty of that.
     * @tparam Iterator   allows walking through the sequence data.
     * @tparam Alphabet   allows interpretation of the sequence
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQSequenceId
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
     * @tparam Iterator   allows walking through the sequence data.
     * @tparam Alphabet   allows interpretation of the sequence
     * @tparam Quality    data type for quality scores for each base
     * @tparam IdType     the format and components of the id of a sequence. Differs from FASTA to FASTQ and based on length of reads.  defaults to FASTQSequenceId
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
     * @brief Functoid encapsulating increment functionality traversing data in an iterator record by record, by parsing and searching for record boundaries.
     * @details   The purpose of this class is to abstract the increment operations in
     *    an iterator that transforms a block of data into a collection of FASTQ records. Using this class allows us to reuse
     *    an sequence parsing iterator, but potentially for different sequence record formats.
     *
     *    This class assumes that the iterator at start is pointing to the beginning of a FASTQ record already.
     *    Increment is done by walking through the data and recording the positions of the start and end of each
     *    of the 4 lines in a FASTQ record.  Because we are always searching to the end of a record (by reading 4 lines always),
     *    successive calls to increment always are aligned correctly with record boundaries
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

      protected:
        /// constant representing the EOL character
        static const typename std::iterator_traits<Iterator>::value_type eol = '\n';

        /// alias for this type, for use inside this class.
        using type = FASTQParser<Iterator, Alphabet, Quality>;

      public:

        /// Sequence type, conditionally set based on Quality template param to either normal Sequence or SequenceWithQuality.
        typedef typename std::conditional<std::is_void<Quality>::value,
              Sequence<Iterator, Alphabet>,
              SequenceWithQuality<Iterator, Alphabet, Quality> >::type      SequenceType;


        /// default constructor.
        FASTQParser() {};

        /// default destructor
        virtual ~FASTQParser() {};

      protected:
        /**
         * @brief  search for first non-EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
         * @details       iter can point to the previous EOL, or a nonEOL character.
         * @param[in/out] iter    iterator to advance
         * @param[in]     end     position to stop the traversal
         * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
         * @return        iterator at the new position, where the Non EOL char is found, or end.
         */
        inline Iterator& findNonEOL(Iterator& iter, const Iterator& end, size_t &offset) {
          while ((*iter == type::eol)  &&  (iter != end)) {
            ++iter;
            ++offset;
          }
          return iter;
        }

        /**
         * @brief  search for first EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
         * @details       iter can point to a nonEOL char, or an EOL char.
         * @param[in/out] iter    iterator to advance
         * @param[in]     end     position to stop the traversal
         * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
         * @return        iterator at the new position, where the EOL char is found, or end.
         */
        inline Iterator& findEOL(Iterator& iter, const Iterator& end, size_t &offset) {
          while ((*iter != type::eol)  &&  (iter != end)) {
            ++iter;
            ++offset;
          }
          return iter;
        }


        /**
         * @brief function to populate the quality score. defined only when Quality type is not void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Q = Quality>
        inline typename std::enable_if<!std::is_void<Q>::value>::type
        setQualityIterators(const Iterator & start, const Iterator & end, SequenceType& output) {
          output.qualBegin = start;
          output.qualEnd = end;
        }

        /**
         * @brief function to populate the quality score. defined only when Quality type is void
         * @param[in] start   beginning of the quality score character sequence.
         * @param[in] end     end of the quality score character sequence.
         * @param[out] output  The output sequence record type to modify
         */
        template <typename Q = Quality>
        inline typename std::enable_if<std::is_void<Q>::value>::type
        setQualityIterators(const Iterator & start, const Iterator & end, SequenceType& output) {}


        /**
         * @brief function to populate the sequence iterators.
         * @param[in] start   beginning of the sequence.
         * @param[in] end     end of the sequence.
         * @param[out] output  The output sequence record type to modify
         */
        inline void setSequenceIterators(const Iterator & start, const Iterator & end, SequenceType& output) {
          output.seqBegin = start;
          output.seqEnd = end;
        }

        /**
         * @brief constructs an IOException object with the relevant debug data.
         * @param errType     string indicating the source of the error.  user choice
         * @param start       iterator pointing to beginning of the data in question
         * @param end         iterator pointing to end of the data in question
         * @param startOffset offset for the beginning of the data in question
         * @param endOffset   offset for the end of the data in question
         */
        void handleError(const std::string& errType, const Iterator &start, const Iterator &end, const size_t& startOffset, const size_t& endOffset) throw (bliss::io::IOException) {
          std::stringstream ss;
          ss << "ERROR: did not find "<< errType << " in " << startOffset << " to " << endOffset << std::endl;
          ss << "  offending string is \"" << std::endl;
          std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type> oit(ss);
          std::copy(start, end, oit);
          ss << "\".";
          throw bliss::io::IOException(ss.str());

        }

      public:
        /**
         * @brief increments the iterator to the beginning of the next record, while saving the current record in memory and update the record id.
         * @details   Called by the containing iterator during operator++ or operator+=
         *            Internally, this function finds the begin and end of the text lines 4 times,
         *            and populate the output sequence type with start and end positions of lines 2 and 4.
         *
         *            Each call to findNonEOL and findEOL may not move the iterator, if iter is already at end, or if iter already points to
         *            a nonEOL or EOL character, respectively.  This means that the output may have Sequence or Quality Score start and end iterators
         *            point to the block's end.  We check these conditions and throw an exception so the incomplete Sequence record does not get used
         *            inadvertently.
         *
         *            Iterator can be a forward iterator or better.
         *
         * @param[in/out] iter          source iterator, pointing to data to traverse
         * @param[in]     end           end of the source iterator - not to go past.
         * @param[in/out] offset        starting position in units of source character types (e.g. char) from the beginning of the file.
         * @param[in/out] output        updated sequence type, values updated.
         * @throw IOException           if parse failed, throw exception
         */
        void increment(Iterator &iter, const Iterator &end, size_t &offset, SequenceType& output) throw (bliss::io::IOException)
        {
          //== first make a copy of iter so we can later use the original in exception handling.
          Iterator orig_iter = iter;
          size_t orig_offset = offset;

          // local variables, temp storage.
          Iterator lstart = end, lend = end;

          //==== parse 4 lines (always, because FASTQ record has 4 lines) for start and end.
          //==== okay to call findNonEOL or findEOL 4 times - if iter at end, won't advance. and okay for output to have end.
          // each call to findNonEOL will find first char, and trim the leading \n

          //== find 1st line, also the starting point.
          // find start of line, \nX or end.
          findNonEOL(iter, end, offset);
          // offset pointing to first start, save it as the id.  - either \n or end
          output.id.composite = offset;
          // then find the end of that line  - either \n or end.
          findEOL(iter, end, offset);

          // == find 2nd line, and save it.
          lstart = findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = findEOL(iter, end, offset);
          setSequenceIterators(lstart, lend, output);

          // == find 3rd line, and discard it.
          findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          findEOL(iter, end, offset);

          // == find 4th line, and save it.  this depends on conditional definition of setQualityIterators.
          lstart = findNonEOL(iter, end, offset);
          // then find the end of that line  - either \n or end.
          lend = findEOL(iter, end, offset);
          // don't even bother to call if Quality type is void - avoids copy operations.
          if (!std::is_void<Quality>::value) setQualityIterators(lstart, lend, output);

          // lend at this point is pointing to the last \n or at end.
          if (iter != end) {  // if at \n, advance it by 1 position
            ++iter;  ++offset;
          }

          //== now check for error conditions.
          // if either sequence or if required quality score were not found, then failed parsing
          if (output.seqBegin == output.seqEnd) {
            handleError("sequence", orig_iter, iter, orig_offset, offset);
          }
          if (!std::is_void<Quality>::value && (output.qualBegin == output.qualEnd)) {
            output.seqBegin = output.seqEnd;  // force seq data not to be used.
            handleError("required quality score", orig_iter, iter, orig_offset, offset);
          }
          // leave iter and offset advanced.


        }


    };

    /// template class' static variable definition (declared and initialized in class)
    template<typename Iterator, typename Alphabet, typename Quality>
    const typename std::iterator_traits<Iterator>::value_type bliss::io::FASTQParser<Iterator, Alphabet, Quality>::eol;

    /**
     * @class bliss::io::SequencesIterator
     * @brief Iterator for parsing and traversing a block of data to access individual sequence records (of some file format)
     * @details The iterator is initialized with start and end iterators pointing to the source data.
     *    The source data should be of primitive type such as unsigned char.  The input iterators should be aligned
     *    to record boundaries.
     *
     *    The iterator uses a Parser Functoid to walk through the data one sequence record at a time.  Different Parser
     *    types may be used, e.g. FASTA, FASTQ, to parse different file formats.
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
     *      4. iterator states:  before first use, during use, after reaching the end.
     *          before:  operator* points to the first output, no need to call ++ first.  so this is same case as "during".
     *          during:  operator* gets the most current output that's a result of calling ++, ++ moves pointer forward
     *          after:   operator* returns empty result, ++ is same as no-op.
     *     Based on these, we need to
     *      1. cache the output
     *      2. keep a current iterator (for comparison between SequencesIterators)
     *      3. keep a next iterator (for parsing)
     *      4. have constructor parse the first entry, so that dereference functions properly
     *      5. check for end condition (next == end) and return if at end.
     *
     *      The Parser functoid then only has 1 simple task:  to increment, as it does not need to store locally.
     *
     * @note this iterator is not a forward or bidirectional iterator, as traversing in reverse is not implemented.
     *
     * @tparam Parser     Functoid type to supply specific logic for increment and dereference functionalities.
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
            typename Parser::SequenceType,
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
        size_t offset;


      public:
        /**
         * @brief constructor, initializes with a start and end.  this represents the start of the output iterator
         * @details   at construction, the internal _curr, _next iterators are set to the input start iterator, and _end is set to the input end.
         *            then immediately the first record is parsed, with output populated with the first entry and next pointer advanced.
         *
         * @param f       parser functoid for parsing the data.
         * @param start   beginning of the data to be parsed
         * @param end     end of the data to be parsed.
         * @param _range  the Range associated with the start and end of the source data.  coordinates relative to the file being processed.
         */
        SequencesIterator(const Parser & f, const Iterator& start,
                       const Iterator& end, const size_t &_offset)
            : seq(), _curr(start), _next(start), _end(end), parser(f), offset(_offset)
        {
          // parse the first entry, if start != end.
          // effect: _next incremented to next parse start point.
          //         offset updated to next start point,
          //         output updated with either empty (if _next at end)  or the new output
          if (_next != _end) parser.increment(_next, _end, offset, seq);;
        }

        /**
         * @brief constructor, initializes with only the end.  this represents the end of the output iterator.
         * @details   at construction, the internal _curr, _next, and _end iterators are set to the input end iterator.
         *          this instance therefore does not traverse and only serves as marker for comparing to the traversing SequencesIterator.
         * @param end     end of the data to be parsed.
         */
        SequencesIterator(const Iterator& end)
            : seq(), _curr(end), _next(end), _end(end), parser(), offset(std::numeric_limits<size_t>::max())
        {
        }

        /**
         * @brief default copy constructor
         * @param Other   The SequencesIterator to copy from
         */
        SequencesIterator(const type& Other)
            : seq(Other.seq), _curr(Other._curr), _next(Other._next), _end(Other._end),
              parser(Other.parser),  offset(Other.offset)
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
          offset = Other.offset;
          return *this;
        }

        /// default constructor deleted.
        SequencesIterator() = delete;


        /**
         * @brief   iterator's pre increment operation
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

            //== then try parsing.  this advances _next and offset, ready for next ++ call, and updates output cache.
            parser.increment(_next, _end, offset, seq);
          }

          //== states updated.  return self.
          return *this;
        }

        /**
         * @brief   iterator's post increment operation
         * @details to increment, the class uses the parser functoid to parse and traverse the input iterator
         *          At the end of the traversal, a FASTQ sequence record has been made and is kept in the parser functoid
         *          as state variable, and the input iterator has advanced to where parsing would start for the next record.
         *          The corresponding starting of the current record is saved as _curr iterator.
         *
         *          conversely, to get the starting position of the next record, we need to parse the current record.
         *
         * @note    for end iterator, because _curr and _end are both pointing to the input's end, "increment" does not do anything.
         * @param   unnamed
         * @return  copy of self incremented.
         */
        type operator ++(int)
        {
          type output(*this);
          return ++output;
        }

        /**
         * @brief   accessor for the internal base iterator, pointing to the current position
         * @return  base iterator
         */
        Iterator& getBaseIterator()
        {
          return _curr;
        }
        /**
         * @brief   const accessor for the internal base iterator, pointing to the current position
         * @return  const base iterator
         */
        const Iterator& getBaseIterator() const
        {
          return _curr;
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators match.
         */
        bool operator ==(const type& rhs) const
        {
          return _curr == rhs._curr;
        }

        /**
         * @brief     comparison operator.
         * @details   compares the underlying base iterator positions.
         * @param rhs the SequenceIterator to compare to.
         * @return    true if this and rhs SequenceIterators do not match.
         */
        bool operator !=(const type& rhs) const
        {
          return _curr != rhs._curr;
        }

        /**
         * @brief dereference operator
         * @return    a const reference to the cached sequence object
         */
        const ValueType &operator *() const
        {
          return seq;
        }

        /**
         * @brief pointer dereference operator
         * @return    a const reference to the cached sequence object
         */
        const ValueType *operator ->() const {
          return &seq;
        }


    };

  } // iterator
} // bliss
#endif /* SequencesIterator_HPP_ */
