/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file		io_exception.hpp
 * @ingroup io
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief   derived from std exception to represent an IO exception
 * @details
 *
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
