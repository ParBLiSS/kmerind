/**
 * @file    distributed_map_base.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_
#define SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_

#include <containers/container_utils.hpp>
#include <functional>
#include <iterator>

namespace dsc
{

  /**
   * @brief base class for a IteratorFilter.  uses CRTP
   * @details   goal of the class is to provide subclasses that adhere to the same interface without introducing dynamic polymorphism.
   *        4 types of subclasses are possible, the first is provided here.
   *          1. static so that every evaluation returns the same "true" or "false"
   *          2. depend on element only.  operator(iter, iter) always return 0, indicating dependence on element.
   *          3. depend on range only.  constructor evaluates, and the value is return for both operators.
   *          4. combination of both.  developer should implement the appropriate logic.
   */
  template<typename Iter, typename Derived>
  struct IteratorFilter {


      /// default constructor
      IteratorFilter() {};

      /// constructor with start and end iterator.  This is to allow for evaluation based on range alone, or range + element.
      IteratorFilter(Iter b, Iter e) {};

      /// return result for an element.  value may be based on the range.  subclass define eval function.
      bool operator()(typename ::std::iterator_traits<Iter>::value_type const & x) {
        return static_cast<Derived>(*this).eval(x);
      }

      /// return result for a range.  subclass define eval function
      int operator()() {
        return static_cast<Derived>(*this).eval();
      }
  };

  struct Identity {
      template <typename T>
      bool operator()(T const & x) const { return true; }
      template <typename Iter>
      bool operator()(Iter b, Iter e) const { return true; }
  };


  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less = ::std::less<Key>,
      class Equal = ::std::equal_to<Key>,
      class Alloc = ::std::allocator< ::std::pair<Key, T> >
  >  class map_base {

    protected:

      static KeyTransform<Key> trans;

      template <typename Comparator>
      struct TransformedComp {
          Comparator comp;
          inline bool operator()(Key const & x, Key const & y) const {
            return comp(trans(x), trans(y));
          }
          template<typename V>
          inline bool operator()(::std::pair<Key, V> const & x, Key const & y) const {
            return this->operator()(x.first, y);
          }
          template<typename V>
          inline bool operator()(::std::pair<const Key, V> const & x, Key const & y) const {
            return this->operator()(x.first, y);
          }
          template<typename V>
          inline bool operator()(Key const & x, ::std::pair<Key, V> const & y) const {
            return this->operator()(x, y.first);
          }
          template<typename V>
          inline bool operator()(Key const & x, ::std::pair<const Key, V> const & y) const {
            return this->operator()(x, y.first);
          }
          template<typename V>
          inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) const {
            return this->operator()(x.first, y.first);
          }
          template<typename V>
          inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) const {
            return this->operator()(x.first, y.first);
          }

      };


//      struct Greater {
//          Less lt;
//          inline bool operator()(Key const &x, Key const &y) const {
//            return lt(y, x);
//          }
//      };

      using TransformedEqual = TransformedComp<Equal>;
      using TransformedLess = TransformedComp<Less>;
//      using TransformedGreater = TransformedComp<Greater>;
      static TransformedEqual equal;
      static TransformedLess less;
//      static TransformedGreater greater;

      mutable size_t key_multiplicity;

      // communication stuff...
      MPI_Comm comm;
      int comm_size;
      int comm_rank;


      // defined Communicator as a friend
      friend Comm;



      template <class Iter>
      static void sort_ascending(Iter b, Iter e) {
          if (!::std::is_sorted(b, e, less)) {
            ::std::sort(b, e, less);
          }
      }

      virtual void local_clear() noexcept = 0;
      //virtual void local_reserve(size_t n) = 0;
      virtual bool local_empty() const noexcept = 0;
      virtual size_t local_size() const noexcept = 0;


      map_base(MPI_Comm _comm, int _comm_size) :
          key_multiplicity(1), comm(_comm), comm_size(_comm_size) {
        MPI_Comm_rank(comm, &comm_rank);
      }

    public:
      virtual ~map_base() {};

      virtual size_t update_multiplicity() = 0;
      virtual std::vector<std::pair<Key, T> > to_vector() const = 0;
      virtual void to_vector(std::vector<std::pair<Key, T> > & result) const  = 0;
      virtual std::vector<Key> keys() const = 0;
      virtual void keys(std::vector<Key> & result) const = 0;

      /// check if empty.
      virtual bool empty() const noexcept {
        if (comm_size == 1)
          return this->local_empty();
        else // all reduce
          return mxx::test_all(this->local_empty(), comm);
      }

      /// get size of distributed container
      virtual size_t size() const noexcept {
        size_t s = this->local_size();
        if (comm_size == 1)
          return s;
        else
          return mxx::allreduce(s, comm);
      }


      /// clears the unordered_map
      void clear() noexcept {
        // clear + barrier.
        this->local_clear();
        if (comm_size > 1) MPI_Barrier(comm);
      }

  };



  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  KeyTransform<Key> map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::trans;

  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedLess
      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::less;

  template<typename Key, typename T,
      class Comm,
      template <typename> class KeyTransform,
      class Less,
      class Equal,
      class Alloc>
  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedEqual
      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::equal;

//  template<typename Key, typename T,
//      class Comm,
//      template <typename> class KeyTransform,
//      class Less,
//      class Equal,
//      class Alloc>
//  typename map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::TransformedGreater
//      map_base<Key, T, Comm, KeyTransform, Less, Equal, Alloc>::greater;

}



#endif /* SRC_CONTAINERS_DISTRIBUTED_MAP_BASE_HPP_ */
