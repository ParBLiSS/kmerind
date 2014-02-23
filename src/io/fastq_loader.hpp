/**
 * fastq_loader.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FASTQLOADER_HPP_
#define FASTQLOADER_HPP_

#include "io/file_loader.hpp"
#include <vector>

#include "iterators/fastq_iterator.hpp"

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
    class fastq_loader : public file_loader
    {
      public:
        virtual ~fastq_loader();
        fastq_loader(std::string const & _filename, file_loader::range_type const & range,
                     size_t const &total);

      protected:
        std::vector<bliss::iterator::fastq_sequence<char*>> seqPositions;
        uint64_t seqIdStart;

        file_loader::range_type align_to_sequence(std::string const & _filename,
                                                            file_loader::range_type const & input,
                                                            size_t const & total) throw(io_exception);

        size_t find_sequence_start(char const* _data, file_loader::range_type const& range) throw(io_exception);


      public:
        void get_sequence_positions() throw(io_exception);
        void test();

    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQLOADER_HPP_ */
