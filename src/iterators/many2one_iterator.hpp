/**
 * @file    many2one_iterator.hpp
 * @ingroup iterators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the many2one iterator.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */


#ifndef BLISS_ITERATORS_MANY2ONE_ITERATOR_HPP
#define BLISS_ITERATORS_MANY2ONE_ITERATOR_HPP

#include <iterator>

#include "common/bit_ops.hpp"
#include "iterators/iterator_utils.hpp"
#include "utils/function_traits.hpp"
#include "iterators/transform_iterator.hpp"

namespace bliss
{
namespace iterator
{

template<typename Iterator, typename Functor>
class _shared_many2one_iterator
  : public _shared_transforming_iterator<Iterator, Functor, Iterator, Iterator>
{

protected:

  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<_shared_many2one_iterator> traits;

  /// The difference type for the step size
  typedef typename traits::difference_type difference_type;

  /// The type of the derived class
  typedef _shared_transforming_iterator<Iterator, Functor, Iterator, Iterator> derived_type;


public:
  /*******************************
   *  Member accessor functions  *
   *******************************/

  /**
   * @brief  Returns the `m` in `m -> 1`, i.e. the step size of this iterator.
   *
   * @return The step size m.
   */
  const difference_type& getStepSize() const
  {
    return this->_m;
  }

protected:
  /**********************
   *  Member variables  *
   **********************/

  /// The end of the input sequence
  Iterator _end;

  /// the number of base elements per compressed element (i.e the `m` in m->1)
  difference_type _m;

  /**************************************************
   *  Private constructor (no direct construction)  *
   **************************************************/
  _shared_many2one_iterator() {}

  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The many2one functor.
   * @param m           The number of elements combined into one.
   */
  _shared_many2one_iterator(const Iterator& base_iter, const Iterator& end_iter,
                           const Functor& f, const difference_type m)
      : derived_type(base_iter, f), _end(end_iter), _m(m)
  {
  }
};


template<typename BaseIterator, typename Compressor, typename DerivedIterator>
class _many2one_iterator_dir : public _shared_many2one_iterator<BaseIterator, Compressor>
{
protected:
  // the base iterator traits
  typedef std::iterator_traits<BaseIterator> base_traits;

  // the type of the derived iterator, so that operator functions do not have
  // to be overloaded via polymorphism
  typedef DerivedIterator type;

  // assert that the base iterator is NOT a random access iterator
  static_assert(!std::is_same<typename base_traits::iterator_category,
                             std::random_access_iterator_tag>::value,
                "The (Bi)DirectedIterator implementation must be"
                "based on a (Bi)DirectedIterator (i.e. Input, FWD or Bidir)");


protected:
  /**********************
   *  Member variables  *
   **********************/

  // The next iterator (read ahead), only used for bidirectional, forward and
  // input iterators
  BaseIterator _next;

  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<type> traits;

  // base class type
  typedef _shared_many2one_iterator<BaseIterator, Compressor> base_class_type;

public:

  /**********************************
   *  Typedefs for iterator traits  *
   **********************************/

  /// The iterator category of this iterator
  typedef typename traits::iterator_category iterator_category;
  /// The value type of this iterator
  typedef typename traits::value_type value_type;
  /// The difference type of this iterator
  typedef typename traits::difference_type difference_type;
  /// The reference type of a value of this iterator
  typedef typename traits::reference reference_type;
  /// The pointer type of a value of this iterator
  typedef typename traits::pointer pointer_type;

public:
  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The many2one functor.
   * @param m           The number of elements combined into one.
   */
  _many2one_iterator_dir(const BaseIterator& base_iter, const BaseIterator& end_iter,
                            const Compressor& f, const difference_type m)
      : base_class_type(base_iter, end_iter, f, m), _next(base_iter)
  {
  }

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This returns the compressed value (m -> 1 combination) from the
   * next `m` elements of the base iterator.
   *
   * @return    The value at the current position.
   */
  value_type operator*()
  {
    if (this->_base != this->_next)
    {
      this->_next = this->_base;
    }
    // the functor iteraters _next forward by at most _m steps
    return this->_f(this->_next, this->_end);
  }

  /**
   * @brief     Pre-increment operator.
   *
   * Advances this iterator by one position.
   *
   * @return    A reference to this.
   */
  type& operator++()
  {
    //  advance the value by m
    if (this->_next == this->_base)
    {
      // advance the min of m and the end of the input sequence
      iter_tools<BaseIterator>::advance_at_most(this->_next, this->_end, this->_m);
    }
    this->_base = this->_next;
    return *dynamic_cast<type*>(this);
  }

  /****************************
   *  Bidirectional Iterator  *
   ****************************/

  static constexpr bool is_bidir = std::is_same<iterator_category,
                                std::bidirectional_iterator_tag>::value;
  // bidirectional iterator
  /**
   * @brief     Pre-decrement operator.
   *
   * Reduces this iterator by one position.
   *
   * @return    A reference to this.
   */
  template<typename T = type>
  typename std::enable_if<is_bidir, T>::type&
  operator--()
  {
    // go back `m`, cast to signed int before negating
    // TODO: what happens if we are at an uneven end?
    std::advance(this->_base, - static_cast<int>(this->_m));
    return *dynamic_cast<type*>(this);
  }

  /**
   * @brief     Post-decrement operator.
   *
   * Reduces this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-decremented iterator.
   */
  template<typename T = type>
  typename std::enable_if<is_bidir, type>::type
  operator--(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type*>(this));
    this->operator--();
    return tmp;
  }

  /**
   * @brief     Post-increment operator.
   *
   * Advances this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-incremented iterator.
   */
  type operator++(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type*>(this));
    this->operator++();
    return tmp;
  }
};

template<typename BaseIterator, typename Compressor, typename DerivedIterator>
class _many2one_iterator_ra : public _shared_many2one_iterator<BaseIterator, Compressor>
{
protected:
  // the base iterator traits
  typedef typename std::iterator_traits<BaseIterator> base_traits;

  // the base class
  typedef _shared_many2one_iterator<BaseIterator, Compressor> base_class_type;

  // the type of the derived iterator, so that operator functions do not have
  // to be overloaded via polymorphism
  typedef DerivedIterator type;

  // assert that the base iterator is actually a random access iterator
  static_assert(std::is_same<typename base_traits::iterator_category,
                             std::random_access_iterator_tag>::value,
                "The RandomAccessIterator implementation must be"
                "based on a RandomAccessIterator");

  /// The std::iterator_traits of this iterator
  typedef typename std::iterator_traits<base_class_type> traits;

public:

  /**********************************
   *  Typedefs for iterator traits  *
   **********************************/

  /// The iterator category of this iterator
  typedef typename traits::iterator_category iterator_category;
  /// The value type of this iterator
  typedef typename traits::value_type value_type;
  /// The difference type of this iterator
  typedef typename traits::difference_type difference_type;
  /// The reference type of a value of this iterator
  typedef typename traits::reference reference_type;
  /// The pointer type of a value of this iterator
  typedef typename traits::pointer pointer_type;

public:
  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The many2one functor.
   * @param m           The number of elements combined into one.
   */
  _many2one_iterator_ra(const BaseIterator& base_iter, const BaseIterator& end_iter,
                           const Compressor& f, const difference_type m)
      : base_class_type(base_iter, end_iter, f, m)
  {
  }

  /************************************************
   *  RandomAccess Iterator supported functions:  *
   ************************************************/

  /************************
   *  specific functions  *
   ************************/

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This returns the compressed value (m -> 1 combination) from the
   * next `m` elements of the base iterator.
   *
   * @return    The value at the current position.
   */
  // same function, but for RandomAccessIterators
  value_type operator*()
  {
    BaseIterator _next = this->_base;
    return this->_f(_next, this->_end);
  }

  /**
   * @brief     Pre-increment operator.
   *
   * Advances this iterator by one position.
   *
   * @return    A reference to this.
   */
  type& operator++()
  {
    iter_tools<BaseIterator>::advance_at_most(this->_base, this->_end, this->_m);
    return *dynamic_cast<type*>(this);
  }

  /**
   * @brief     Pre-decrement operator.
   *
   * Reduces this iterator by one position.
   *
   * @return    A reference to this.
   */
  type& operator--()
  {
    // go back `m`, cast to signed int before negating
    // TODO: what happens if we are at an uneven end?
    std::advance(this->_base, - static_cast<int>(this->_m));
    return *dynamic_cast<type*>(this);
  }

  /* advancing iterator:  operator+ */

  /**
   * @brief     Advances this iterator by `n` positions.
   *
   * @param n   The number of positions to advance.
   * @return    A reference to this after advancing.
   */
  type& operator+=(difference_type n)
  {
    // advance in steps of _m
    if (n < 0)
      *this -= (-n);
    else
      iter_tools<BaseIterator>::advance_at_most(this->_base, this->_end, n * this->_m);
    return *dynamic_cast<type*>(this);
  }

  /* decrementing iterator:  operator- */

  /**
   * @brief     Decrements the iterator by `n` positions.
   *
   * @param n   The number of positions to decrement.
   * @return    A reference to this after decrementing.
   */
  type& operator-=(difference_type n)
  {
    // go back in steps of size `m`
    if (n < 0)
      *this += (-n);
    else
      std::advance(this->_base, -n * static_cast<int>(this->_m));
    return *dynamic_cast<type*>(this);
  }

  /**
   * @brief     Returns the difference between this and the given iterator.
   *
   * The difference is the number of elements between the two iterator
   * positions.
   *
   * @param other   The iterator to subtract from `this`.
   *
   * @return    The difference between this and the given iterator.
   */
  difference_type operator-(const type& other)
  {
    // divide base difference by `m`
    return intCeil<difference_type>(this->_base - other._base, this->_m);
  }

  /***********************
   *  derived functions  *
   ***********************/
  // i.e functions that are purely based on previosly declared functions
  //     where no specialization to the exact iterator is necessary

  /**
   * @brief     Post-increment operator.
   *
   * Advances this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-incremented iterator.
   */
  type operator++(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type*>(this));
    this->operator++();
    return tmp;
  }

  /**
   * @brief     Post-decrement operator.
   *
   * Reduces this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-decremented iterator.
   */
  type operator--(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type*>(this));
    this->operator--();
    return tmp;
  }

  /**
   * @brief     Returns the n'th element as seen from the current iterator
   *            position.
   *
   * @param n   The offset.
   *
   * @return    The element at offset `n` from the current position.
   */
  value_type operator[](difference_type n)
  {
    // reduce to: *(*this + n)
    // create tmp copy
    type tmp(*dynamic_cast<type*>(this));
    // reduce to already implemented operators
    tmp += n;
    return *tmp;
  }

  /**
   * @brief     Advances a copy of this iterator by `n` positions.
   *
   * @param n   The number of positions to advance.
   * @return    The advanced iterator.
   */
  type operator+(difference_type n)
  {
    type output(*dynamic_cast<type*>(this));
    // reduced to += operator
    output += n;
    return output;
  }

  /**
   * @brief     Advances a copy of the `right` iterator by `n` positions.
   *
   * @param n   The number of positions to advance.
   * @return    The advanced iterator.
   */
  friend type operator+(difference_type n, const type& right)
  {
    // reduced to + operator
    return right + n;
  }

  /**
   * @brief     Decrements a copy of this iterator by `n` positions.
   *
   * @param n   The number of positions to decrement by.
   * @return    The decremented iterator.
   */
  type operator-(difference_type n)
  {
    type output(*dynamic_cast<type*>(this));
    // reduced to -= operator
    output -= n;
    return output;
  }

  /* ordered comparison operators */

  /**
   * @brief  Compares this iterator to the given iterator for the relation <.
   */
  bool operator<(const type& rhs)
  {
    return this->_base < rhs._base;
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation >.
   */
  bool operator>(const type& rhs)
  {
    return this->_base > rhs._base;
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation <=.
   */
  bool operator<=(const type& rhs)
  {
    return this->_base <= rhs._base;
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation >=.
   */
  bool operator>=(const type& rhs)
  {
    return this->_base >= rhs._base;
  }
};



/**
 * @brief many2one Iterator class (m -> 1 mapping)
 *
 * This iterator wraps around a base iterator and combines a fixed number of
 * elements (e.g. `m`) from the base iterator into a single element of a
 * different type.
 *
 * Advancing this iterator by a single position advances the base iterator
 * by `m` positions.
 *
 * The given compression functor takes two parameters as input: a `begin`
 * and an `end` iterator of the iterator base class. The functor must `pack`
 * at most `m` or `end - begin` elements (which ever is smaller) from the
 * underlying iterator sequence into a single value and return this value.
 * The value_type of this iterator is determined from the given functors
 * return type.
 *
 * This iterator "inherits" the iterator properties from the base iterator.
 * I.e. if the base iterator is a forward_iterator, this iterator will also be
 * a forward iterator. Similarly for the iterator tag types:
 *  - input_iterator
 *  - forward_iterator
 *  - bidirectional_iterator
 *  - random_access_iterator
 *
 * Any output iterator concepts are not defined, as the functor can not be
 * reversed.
 *
 * Note: this is loosely patterned after boost's transform iterator.
 *
 * @tparam Compressor   The type of the many2one functor.
 * @tparam Iterator     The type of the base iterator.
 */
template<typename Iterator, typename Compressor>
class many2one_iterator : public std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                std::random_access_iterator_tag>::value,
                                _many2one_iterator_ra<Iterator, Compressor, many2one_iterator<Iterator, Compressor> >,
                                _many2one_iterator_dir<Iterator, Compressor, many2one_iterator<Iterator, Compressor> > >::type
{
protected:

  /***********************
   *  Private type defs  *
   ***********************/

  /// The iterator traits of the wrapped base iterator.
  typedef std::iterator_traits<Iterator> base_traits;



  // the traits (e.g. return type) of the compression function
  typedef bliss::functional::function_traits<Compressor,
      Iterator, Iterator> functor_traits;
  /// The type of this class
  typedef many2one_iterator type;

  /// Whether the iterator is a RandomAccessIterator
  static constexpr bool _is_ra = std::is_same<typename base_traits::iterator_category,
                                 std::random_access_iterator_tag>::value;

  /// Iterator base class type in case the wrapped iterator is RandomAccess
  typedef _many2one_iterator_ra<Iterator, Compressor, many2one_iterator> _ra_base_t;
  /// Iterator base class type in case the wrapped iterator is a directed iterator
  typedef _many2one_iterator_dir<Iterator, Compressor, many2one_iterator> _dir_base_t;
  /// The type of the base class
  typedef typename std::conditional<_is_ra, _ra_base_t, _dir_base_t>::type base_class_type;


public:

  /// The difference type of the iterator
  typedef typename base_traits::difference_type difference_type;

  /******************
   *  Constructors  *
   ******************/

  // class specific constructor
  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The many2one functor.
   * @param m           The number of elements combined into one.
   */
  many2one_iterator(const Iterator& base_iter, const Iterator& end_iter,
                       const Compressor& f, const difference_type m)
    : base_class_type(base_iter, end_iter, f, m)
  {
  }

  /// default constructor
  many2one_iterator()
    : base_class_type(Iterator(), Iterator(), Compressor(), 0) {}

  /// copy constructor
  many2one_iterator(const many2one_iterator& other)
    : base_class_type(other._base, other._end, other._f, other._m)
  {
  }

  /// copy assign operator
  many2one_iterator& operator=(const many2one_iterator& other)
  {
    // check for self assignment
    if (this != &other)
    {
      this->_base = other._base;
      this->_end = other._end;
      this->_f = other._f;
      this->_m = other._m;
    }
    return *this;
  }

  /*
  // specific to at least fwd iterators: default contructable
  template< typename = typename std::enable_if<is_min_fwd && !is_ra>::type >
  many2one_iterator()
      :  _base(), _next(), _end(), _f(), _m(0)
  {
  }
  */

  /*
  // specific to RandomAccessIterators
  template< typename = typename std::enable_if<is_ra>::type >
  many2one_iterator()
      :  _base(), _end(), _f(), _m(0)
  {
  }
  */

};

} // iterator
} // bliss
#endif /* BLISS_ITERATORS_MANY2ONE_ITERATOR_HPP */
