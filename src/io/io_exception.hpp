/**
 * @file		io_exception.hpp
 * @ingroup io
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   derived from std exception to represent an IO exception
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
     * @class			IOException
     * @brief     derived from std exception to represent an IO exception
     * @details   this will be retired in the near future.
     */
    class IOException : public std::exception
    {
      protected:
        /// internal error message string
        std::string message;

      public:
        /**
         * @brief constructor
         * @param _msg    const char* error message for the exception
         */
        IOException(const char* _msg)
        {
          message = _msg;
        }

        /**
         * @brief constructor
         * @param _msg    const string error message for the exception
         */
        IOException(const std::string &_msg)
        {
          message = _msg;
        }

        /**
         * @brief move constructor. takes ownership of the error message content.
         * @param other   the exception to move from.  empty after the move
         */
        IOException(IOException&& other) {
          message.swap(other.message);
        }

        /**
         * @brief move assignment operator..  swaps the internal error messages.
         * @param other   the exception to move from.   empty after the move
         * @return        updated io exception object with internal data of the other exception
         */
        IOException& operator=(IOException&& other) {
          if (this != &other) {
            std::string().swap(message);
            message.swap(other.message);
          }
          return *this;
        }

        /**
         * @brief destructor
         */
        virtual ~IOException() {};

        /**
         * @brief   overridden method to retrieve the content of the error message.
         * @return  actual error message.
         */
        virtual const char* what() const throw ()
        {
          return message.c_str();
        }
    };
  } /* namespace io */
} /* namespace bliss */

#endif /* IO_EXCEPTION_HPP_ */
