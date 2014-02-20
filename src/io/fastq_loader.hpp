/**
 * fastq_loader.hpp
 *
 *  Created on: Feb 18, 2014
 *      Author: tpan
 */

#ifndef FASTQLOADER_HPP_
#define FASTQLOADER_HPP_

#include "io/file_loader.hpp"

namespace bliss
{
  namespace io
  {

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
        file_loader::range_type adjust_to_record_boundaries(std::string const & _filename,
                                                            file_loader::range_type const & input,
                                                            size_t const & total) throw(io_exception);

        size_t search_for_record_boundary(char const* _data, file_loader::range_type const& range) throw(io_exception);
    };

  } /* namespace io */
} /* namespace bliss */
#endif /* FASTQLOADER_HPP_ */
