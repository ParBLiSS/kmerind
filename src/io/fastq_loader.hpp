/**
 * fastq_loader.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FASTQLOADER_HPP_
#define FASTQLOADER_HPP_

#include <cstdint>


#include <vector>


#include "io/file_loader.hpp"
#include "io/fastq_iterator.hpp"

namespace bliss
{
  namespace io
  {

//    struct parseSequence {
//        fastq_sequence operator()()
//    };

    /**
     *
     */
    template <typename Alphabet, typename Quality = void>
    class fastq_loader : public file_loader
    {
      public:
        typedef bliss::io::fastq_parser<char*, Alphabet, Quality> ParserType;
        typedef typename ParserType::SeqType                      SequenceType;
        typedef bliss::io::fastq_iterator<char*, ParserType>      IteratorType;


        virtual ~fastq_loader();
        fastq_loader(const std::string &_filename,
                     const file_loader::range_type &range, const size_t &total,
                     bool preload = false);

        IteratorType begin();
        IteratorType end();
      protected:
        file_loader::range_type align_to_sequence(
            const std::string &_filename, const file_loader::range_type &input,
            const size_t &total) throw (io_exception);

        size_t find_sequence_start(char const* _data,
                                   const file_loader::range_type &range)
                                       throw (io_exception);

        // for parallel scan way of assigning ids.
//        std::vector< SequenceType > seqPositions;
        uint64_t seqIdStart;


    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQLOADER_HPP_ */
