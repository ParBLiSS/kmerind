/**
 * @file		io_exception.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef IO_EXCEPTION_HPP_
#define IO_EXCEPTION_HPP_

namespace bliss
{
  namespace io
  {

    /**
     * @class			bliss::io::IOException
     * @brief
     * @details
     *            Note: have move constructor and assignment operators.
     */
    class IOException : public std::exception
    {
      protected:
        std::string message;

      public:
        IOException(const char* _msg)
        {
          message = _msg;
        }

        IOException(const std::string &_msg)
        {
          message = _msg;
        }

        IOException(IOException&& other) {
          message.swap(other.message);
        }

        IOException& operator=(IOException&& other) {
          if (this != &other) {
            std::string().swap(message);
            message.swap(other.message);
          }
          return *this;
        }

        virtual ~IOException() {};

        virtual const char* what() const throw ()
        {
          return message.c_str();
        }
    };
  } /* namespace io */
} /* namespace bliss */

#endif /* IO_EXCEPTION_HPP_ */
