/**
 * @file    CopyableAtomic.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef COPYABLEATOMIC_HPP_
#define COPYABLEATOMIC_HPP_

namespace bliss
{
  namespace concurrent
  {

    /**
     * @class    bliss::concurrent::CopyableAtomic
     * @brief
     * @details
     *
     */
    template <typename T>
    struct CopyableAtomic
    {
      std::atomic<T> val;

      CopyableAtomic() : val(T()) {}


      explicit CopyableAtomic(const T& v) : val(v) {}
      explicit CopyableAtomic(const std::atomic<T> & v) : val(v.load()) {}
      explicit CopyableAtomic(const CopyableAtomic<T> & other) :
          val(other.val.load()) {}
      explicit CopyableAtomic(CopyableAtomic<T>&& other) : val(other.val.load())
      {
        other.val.store(T());
      }
      CopyableAtomic& operator=(const CopyableAtomic<T> & other) {
        val = other.val.load();
        return *this;
      }

      CopyableAtomic& operator=(CopyableAtomic<T> && other) {
        val = other.val.load();
        other.val.store(T());
        return *this;
      }
    };

  } /* namespace concurrent */
} /* namespace bliss */

#endif /* COPYABLEATOMIC_HPP_ */
