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
    class fastq_loader : public file_loader
    {
      public:
        virtual ~fastq_loader();
        fastq_loader(const std::string &_filename,
                     const file_loader::range_type &range, const size_t &total,
                     bool preload = false);

      protected:

        file_loader::range_type align_to_sequence(
            const std::string &_filename, const file_loader::range_type &input,
            const size_t &total) throw (io_exception);

        size_t find_sequence_start(char const* _data,
                                   const file_loader::range_type &range)
                                       throw (io_exception);

        // for parallel scan way of assigning ids.
        std::vector<bliss::iterator::fastq_sequence<char*>> seqPositions;
        uint64_t seqIdStart;

      public:

        // testing.
//        void test();
//        void test2();
//        void get_sequence_positions() throw (io_exception);
//        void assign_sequence_ids() throw (io_exception);

        typedef bliss::iterator::fastq_iterator<
            bliss::iterator::fastq_parser<char*>, char*> iterator;

        iterator begin();
        iterator end();

    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQLOADER_HPP_ */
