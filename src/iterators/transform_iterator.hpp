/**
 * @file    transform_iterator.hpp
 * @ingroup iterators
 * @author  Tony Pan <tpan7@gatech.edu>
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the transforming iterator.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */

#ifndef BLISS_ITERATORS_TRANSFORM_ITERATOR_HPP
#define BLISS_ITERATORS_TRANSFORM_ITERATOR_HPP

#include <iterator>
#include "utils/function_traits.hpp"

namespace bliss
{

namespace iterator
{

template<typename Iterator, typename Functor, typename ... FunctorArgs>
class _shared_transforming_iterator : public std::iterator<
    /* inherit interator traits from std::iterator: */
    // 1) iterator type tag (same as base)
    typename std::iterator_traits<Iterator>::iterator_category,
    // 2) value type = return type of compression function
    typename std::remove_reference<
      typename bliss::functional::function_traits<Functor, FunctorArgs ... >::return_type>::type,
    // 3) difference_type = same as base iterator
    typename std::iterator_traits<Iterator>::difference_type>
{
protected:
  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<_shared_transforming_iterator> traits;

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
  /*******************************
   *  Member accessor functions  *
   *******************************/

  /**
   * @brief     Returns the functor.
   *
   * @return    The functor.
   */
  Functor& getFunctor()
  {
    return _f;
  }

  /**
   * @brief     Returns the functor.
   *
   * @return    The functor.
   */
  const Functor& getFunctor() const
  {
    return _f;
  }

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  Iterator& getBaseIterator()
  {
    return _base;
  }

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  const Iterator& getBaseIterator() const
  {
    return _base;
  }

  /**
   * @brief     Returns whether this and the given iterator point to the same
   *            positon.
   *
   * @return    Whether this and the given iterator point to the same position.
   */
  inline bool operator==(const _shared_transforming_iterator& rhs) const
  {
    return _base == rhs._base;
  }

  /**
   * @brief     Returns whether this and the given iterator are different.
   *
   * @return    Whether this and the given iterator are different.
   */
  inline bool operator!=(const _shared_transforming_iterator& rhs) const
  {
    return _base != rhs._base;
  }

protected:
  /**********************
   *  Member variables  *
   **********************/

  /// The current base iterator position
  Iterator _base;

  /// The functor
  Functor _f;

  /**************************************************
   *  Private constructor (no direct construction)  *
   **************************************************/
  _shared_transforming_iterator() {}

  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The functor.
   */
  _shared_transforming_iterator(const Iterator& base_iter, const Functor& f)
      : _base(base_iter), _f(f)
  {
  }
  /// destructor
  virtual ~_shared_transforming_iterator() {}
};


template<typename BaseIterator, typename Functor, typename DerivedIterator>
class _transforming_iterator_dir : public _shared_transforming_iterator<BaseIterator, Functor, typename std::iterator_traits<BaseIterator>::value_type>
{
protected:
  // the base iterator traits
  typedef std::iterator_traits<BaseIterator> base_traits;

  // the type of the derived iterator, so that operator functions do not have
  // to be overloaded via polymorphism
  typedef DerivedIterator type;


protected:
  /**********************
   *  Member variables  *
   **********************/

  // base class type
  typedef _shared_transforming_iterator<BaseIterator, Functor, typename std::iterator_traits<BaseIterator>::value_type> base_class_type;

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

protected:

  // iterator tags
  static constexpr bool _is_ra = std::is_same<typename base_traits::iterator_category,
                                 std::random_access_iterator_tag>::value;
  static constexpr bool _is_bidir = std::is_same<iterator_category,
                                std::bidirectional_iterator_tag>::value;
  static constexpr bool _is_min_bidir = _is_ra || _is_bidir;

public:
  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The functor.
   */
   _transforming_iterator_dir(const BaseIterator& base_iter, const Functor& f)
      : base_class_type(base_iter, f)
  {
  }

  virtual ~_transforming_iterator_dir() {}

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This is the value at the current base iterator position, transformed with
   * the transformation functor.
   *
   * @return    The value at the current position.
   */
  value_type operator*()
  {
    return this->_f(*this->_base);
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
    ++this->_base;
    return *dynamic_cast<type*>(this);
  }

  /****************************
   *  Bidirectional Iterator  *
   ****************************/
  // bidirectional iterator
  /**
   * @brief     Pre-decrement operator.
   *
   * Reduces this iterator by one position.
   *
   * @return    A reference to this.
   */
  template<typename T = type>
  typename std::enable_if<_is_min_bidir, T>::type&
  operator--()
  {
    --this->_base;
    return *dynamic_cast<type*>(this);
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

  /**
   * @brief     Post-decrement operator.
   *
   * Reduces this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-decremented iterator.
   */
  template<typename T = type>
  typename std::enable_if<_is_min_bidir, T>::type
  operator--(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type>(this));
    this->operator--();
    return tmp;
  }
};

template<typename BaseIterator, typename Functor, typename DerivedIterator>
class _transforming_iterator_ra : public _transforming_iterator_dir<BaseIterator, Functor, DerivedIterator>
{
protected:
  // the base iterator traits
  typedef typename std::iterator_traits<BaseIterator> base_traits;

  // the base class
  typedef _transforming_iterator_dir<BaseIterator, Functor, DerivedIterator> base_class_type;

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
   * @brief Constructor, taking the base iterator and the functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The functor.
   */
  _transforming_iterator_ra(const BaseIterator& base_iter, const Functor& f)
      : base_class_type(base_iter, f)
  {
  }

  virtual ~_transforming_iterator_ra() {}

  /************************************************
   *  RandomAccess Iterator supported functions:  *
   ************************************************/

  /************************
   *  specific functions  *
   ************************/

  /* advancing iterator:  operator+ */

  /**
   * @brief     Advances this iterator by `n` positions.
   *
   * @param n   The number of positions to advance.
   * @return    A reference to this after advancing.
   */
  type& operator+=(difference_type n)
  {
    std::advance(this->_base, n);
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
    std::advance(this->_base, -n);
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
    return this->_base - other._base;
  }

  /***********************
   *  derived functions  *
   ***********************/
  // i.e functions that are purely based on previosly declared functions
  //     where no specialization to the exact iterator is necessary

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
 * @brief Transform Iterator class
 *
 * This iterator wraps around a base iterator and transforms each value of the
 * base iterator with the given transformation function.
 * Dereferencing this iterator is analogous to: `f(*base_iterator)`
 * where `f` is the transformation function
 *
 * The functor has to take the base iterator's value_type as argument
 * and return any type. The value_type of this iterator is set to the return
 * type of the transformation functor.
 *
 * This iterator "inherits" the iterator properties from the base iterator.
 * I.e. if the base iterator is a forward_iterator, this iterator will also be
 * a forward iterator. The same goes for the iterator tag types:
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
 * @tparam Iterator     The type of the base iterator.
 * @tparam Transformer  The type of the tranformer/transformation functor.
 */
template<typename Iterator, typename Transformer>
class transform_iterator : public std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                std::random_access_iterator_tag>::value,
                                _transforming_iterator_ra<Iterator, Transformer, transform_iterator<Iterator, Transformer> >,
                                _transforming_iterator_dir<Iterator, Transformer, transform_iterator<Iterator, Transformer> > >::type
{
protected:

  /***********************
   *  Private type defs  *
   ***********************/

  // the base iterator traits
  typedef std::iterator_traits<Iterator> base_traits;
  // the traits (e.g. return type) of the transformer function
  typedef bliss::functional::function_traits<Transformer,
      typename base_traits::value_type> functor_traits;
  /// The type of this class
  typedef transform_iterator<Iterator, Transformer> type;
  /// Whether the iterator is a RandomAccessIterator
  static constexpr bool _is_ra = std::is_same<typename base_traits::iterator_category,
                                 std::random_access_iterator_tag>::value;

  /// Iterator base class type in case the wrapped iterator is RandomAccess
  typedef _transforming_iterator_ra<Iterator, Transformer, transform_iterator> _ra_base_t;
  /// Iterator base class type in case the wrapped iterator is a directed iterator
  typedef _transforming_iterator_dir<Iterator, Transformer, transform_iterator> _dir_base_t;
  /// The type of the base class
  typedef typename std::conditional<_is_ra, _ra_base_t, _dir_base_t>::type base_class_type;


public:


  /******************
   *  Constructors  *
   ******************/

  // class specific constructor
  /**
   * @brief Constructor, taking the base iterator and the transformation
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The transforming functor.
   */
  transform_iterator(const Iterator& base_iter, const Transformer & f)
      : base_class_type(base_iter, f)
  {
  }

  /**
   * @brief     Default contructor.
   */
  // specific to at least fwd iterators: default contructable
  transform_iterator()
      : base_class_type()
  {
  }

  /// default destructor
  virtual ~transform_iterator() {}

  /// copy assignment operator
  transform_iterator& operator=(const transform_iterator& other)
  {
    if (this != &other)
    {
      this->_base = other._base;
      this->_f = other._f;
    }
    return *this;
  }


};

} // iterator
} // bliss
#endif /* BLISS_ITERATORS_TRANSFORM_ITERATOR_HPP */
