/**
 * @file    one2many_iterator.hpp
 * @ingroup iterators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the one2many iterator.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */


#ifndef BLISS_ITERATORS_ONE2MANY_ITERATOR_HPP
#define BLISS_ITERATORS_ONE2MANY_ITERATOR_HPP

#include <iterator>

#include <common/bit_ops.hpp>
#include <utils/function_traits.hpp>
#include <iterators/transform_iterator.hpp>

namespace bliss
{
namespace iterator
{

template<typename Iterator, typename Functor>
class _shared_one2many_iterator
  : public _shared_transforming_iterator<Iterator, Functor, Iterator, Iterator>
{

protected:

  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<_shared_one2many_iterator> traits;

  /// The difference type for the step size
  typedef typename traits::difference_type difference_type;

  /// The type of the derived class
  typedef _shared_transforming_iterator<Iterator, Functor, Iterator, Iterator> derived_type;


public:
  /*******************************
   *  Member accessor functions  *
   *******************************/

  /**
   * @brief  Returns the `m` in `1 -> m`, i.e. the step size of this iterator.
   *
   * @return The step size m.
   */
  const difference_type& getStepSize() const
  {
    return this->_m;
  }


  /**
   * @brief Returns the current offset within the current element.
   *
   * @return The current offset.
   */
  const difference_type& getOffset() const
  {
    return this->_offset;
  }

  /**
   * @brief Returns the current offset within the current element.
   *
   * @return The current offset.
   */
  difference_type& getOffset()
  {
    return this->_offset;
  }

  /**
   * @brief     Returns whether this and the given iterator point to the same
   *            positon.
   *
   * @return    Whether this and the given iterator point to the same position.
   */
  inline bool operator==(const _shared_one2many_iterator& rhs) const
  {
    return this->_base == rhs._base && this->_offset == rhs._offset;
  }

  /**
   * @brief     Returns whether this and the given iterator are different.
   *
   * @return    Whether this and the given iterator are different.
   */
  inline bool operator!=(const _shared_one2many_iterator& rhs) const
  {
    return this->_base != rhs._base || this->_offset != rhs._offset;
  }

protected:
  /**********************
   *  Member variables  *
   **********************/

  /// the number of base elements per input element (i.e the `m` in 1->m)
  difference_type _m;

  /// The current offset inside the current word
  difference_type _offset;

  /**************************************************
   *  Private constructor (no direct construction)  *
   **************************************************/
  _shared_one2many_iterator() {}

  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements extracted from one.
   */
  _shared_one2many_iterator(const Iterator& base_iter,
                           const Functor& f, const difference_type m)
      : derived_type(base_iter, f), _m(m), _offset(0)
  {
  }

  /**
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements extracted from one.
   * @param offset      The initial offset in the current base iterator. The
   *                    state of the iterator will be equivalent to initalizing
   *                    to `base_iter` and iterating `offset` times.
   */
  _shared_one2many_iterator(const Iterator& base_iter,
                           const Functor& f, const difference_type m,
                           const difference_type offset)
      : derived_type(base_iter, f), _m(m), _offset(offset)
  {
    if (this->_offset >= this->_m)
    {
      difference_type steps = this->_offset / this->_m;
      difference_type rem = this->_offset % this->_m;
      std::advance(this->_base, steps);
      this->_offset = rem;
    }
  }
};


template<typename BaseIterator, typename Functor, typename DerivedIterator>
class _one2many_iterator_dir : public _shared_one2many_iterator<BaseIterator, Functor>
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
  typedef _shared_one2many_iterator<BaseIterator, Functor> base_class_type;

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
  static constexpr bool _is_bidir = std::is_same<typename base_traits::iterator_category,
                                std::bidirectional_iterator_tag>::value;
  static constexpr bool _is_min_bidir = _is_ra || _is_bidir;

public:
  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the one2many functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements combined into one.
   */
  _one2many_iterator_dir(const BaseIterator& base_iter,
                         const Functor& f, const difference_type m)
      : base_class_type(base_iter, f, m)
  {
  }

  /**
   * @brief Constructor, taking the base iterator and the one2many functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements combined into one.
   * @param offset      The initial offset in the current base iterator. The
   *                    state of the iterator will be equivalent to initalizing
   *                    to `base_iter` and iterating `offset` times.
   */
  _one2many_iterator_dir(const BaseIterator& base_iter,
                         const Functor& f, const difference_type m,
                         const difference_type offset)
      : base_class_type(base_iter, f, m, offset)
  {
  }

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This returns the extracted value (1->m) at offset `_offset` from the
   * current base iterator position.
   *
   * @return    The value at the current position.
   */
  value_type operator*()
  {
    // the functor extracts the relevant part and returns the type
    return this->_f(*this->_base, this->_offset);
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
    //  advance the offset by one
    ++this->_offset;
    if (this->_offset == this->_m)
    {
      ++this->_base;
      ++this->_offset = 0;
    }
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
    // at beginning of current element go back in base iterator
    if (this->_offset == 0)
    {
      --this->_base;
      this->_offset = this->m - 1;
    }
    // otherwise just decrease the offset
    else
    {
      --this->_offset;
    }
    return *dynamic_cast<type*>(this);
  }


  /********************************
   *  derived operator functions  *
   ********************************/

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
  typename std::enable_if<is_bidir, T>::type
  operator--(int)
  {
    // create a copy
    type tmp(*dynamic_cast<type*>(this));
    this->operator--();
    return tmp;
  }
};

template<typename BaseIterator, typename Functor, typename DerivedIterator>
class _one2many_iterator_ra : public _one2many_iterator_dir<BaseIterator, Functor, DerivedIterator>
{
protected:
  // the base iterator traits
  typedef typename std::iterator_traits<BaseIterator> base_traits;

  // the base class
  typedef _one2many_iterator_dir<BaseIterator, Functor, DerivedIterator> base_class_type;

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
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements per base element.
   */
  _one2many_iterator_ra(const BaseIterator& base_iter,
                        const Functor& f, const difference_type m)
      : base_class_type(base_iter, f, m)
  {
  }

  /**
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements per base element.
   * @param offset      The initial offset in the current base iterator. The
   *                    state of the iterator will be equivalent to initalizing
   *                    to `base_iter` and iterating `offset` times.
   */
  _one2many_iterator_ra(const BaseIterator& base_iter,
                        const Functor& f, const difference_type m,
                        const difference_type offset)
      : base_class_type(base_iter, f, m, offset)
  {
  }

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
    if (n < 0)
      *this -= (-n);
    else
    {
      // get base iterator steps and new offset
      difference_type steps = (n + this->_offset) / this->_m;
      difference_type offset = (n + this->_offset) % this->m;
      // advance base and set new offset
      std::advance(this->_base, steps);
      this->_offset = offset;
    }

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
    if (n < 0)
      *this += (-n);
    else
    {
      // get base iterator steps and new offset
      if (n <= this->_offset)
      {
        this->_offset -= n;
      }
      else
      {
        // remove offset
        n -= this->_offset;
        this->_offset = 0;
        // then get number of further steps to go back
        difference_type steps = (n - 1) / this->_m + 1;
        difference_type offset = steps*this->_m - n;
        // advance base and set new offset
        std::advance(this->_base, -steps);
        this->_offset = offset;
      }
    }
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
    assert(this->_m == other._m);
    // get base difference
    difference_type base_diff = this->_base - other._base;
    return base_diff * this->_m + (this->_offset - other._offset);
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
    return this->_base < rhs._base
      || (this->_base == rhs._base && this->_offset < rhs._offset);
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation >.
   */
  bool operator>(const type& rhs)
  {
    return this->_base > rhs._base
      || (this->_base == rhs._base && this->_offset > rhs._offset);
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation <=.
   */
  bool operator<=(const type& rhs)
  {
    return this->_base < rhs._base
      || (this->_base == rhs._base && this->_offset <= rhs._offset);
  }

  /**
   * @brief  Compares this iterator to the given iterator for the relation >=.
   */
  bool operator>=(const type& rhs)
  {
    return this->_base > rhs._base
      || (this->_base == rhs._base && this->_offset >= rhs._offset);
  }
};



/**
 * @brief one2many Iterator class (1 -> m extraction)
 *
 * This iterator wraps around a base iterator and returns `m` elements per base
 * element.
 *
 * Advancing this iterator by a single position advances the base iterator
 * only every `m` iterations.
 *
 * The given functor takes two parameters as input: The current base iterator
 * position and the current offset for retrieving the `offset`th element within
 * the one base element.
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
 * @tparam Functor   The type of the one2many functor.
 * @tparam Iterator     The type of the base iterator.
 */
template<typename Iterator, typename Functor>
class one2many_iterator : public std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                std::random_access_iterator_tag>::value,
                                _one2many_iterator_ra<Iterator, Functor, one2many_iterator<Iterator, Functor> >,
                                _one2many_iterator_dir<Iterator, Functor, one2many_iterator<Iterator, Functor> > >::type
{
protected:

  /***********************
   *  Private type defs  *
   ***********************/

  /// The iterator traits of the wrapped base iterator.
  typedef std::iterator_traits<Iterator> base_traits;



  // the traits (e.g. return type) of the compression function
  typedef bliss::functional::function_traits<Functor,
      Iterator, Iterator> functor_traits;
  /// The type of this class
  typedef one2many_iterator type;

  /// Whether the iterator is a RandomAccessIterator
  static constexpr bool _is_ra = std::is_same<typename base_traits::iterator_category,
                                 std::random_access_iterator_tag>::value;

  /// Iterator base class type in case the wrapped iterator is RandomAccess
  typedef _one2many_iterator_ra<Iterator, Functor, one2many_iterator> _ra_base_t;
  /// Iterator base class type in case the wrapped iterator is a directed iterator
  typedef _one2many_iterator_dir<Iterator, Functor, one2many_iterator> _dir_base_t;
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
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements extracted from one.
   */
  one2many_iterator(const Iterator& base_iter,
                       const Functor& f, const difference_type m)
    : base_class_type(base_iter, f, m)
  {
  }

  /**
   * @brief Constructor, taking the base iterator and the one2many
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param f           The one2many functor.
   * @param m           The number of elements extracted from one.
   */
  one2many_iterator(const Iterator& base_iter,
                       const Functor& f, const difference_type m,
                       const difference_type offset)
    : base_class_type(base_iter, f, m, offset)
  {
  }

  /// default constructor
  one2many_iterator()
    : base_class_type(Iterator(), Functor(), 0, 0) {}

  /// copy constructor
  one2many_iterator(const one2many_iterator& other)
    : base_class_type(other._base, other._f, other._m, other._offset)
  {
  }

  /// copy assignment iterator
  one2many_iterator& operator=(const one2many_iterator& other)
  {
    // check for self assignment
    if (this != &other)
    {
      this->_base = other._base;
      this->_f = other._f;
      this->_m = other._m;
      this->_offset = other._offset;
    }
    return *this;
  }

  /*
  // specific to at least fwd iterators: default contructable
  template< typename = typename std::enable_if<is_min_fwd && !is_ra>::type >
  one2many_iterator()
      :  _base(), _next(), _end(), _f(), _m(0)
  {
  }
  */

  /*
  // specific to RandomAccessIterators
  template< typename = typename std::enable_if<is_ra>::type >
  one2many_iterator()
      :  _base(), _end(), _f(), _m(0)
  {
  }
  */

};

} // iterator
} // bliss
#endif /* BLISS_ITERATORS_ONE2MANY_ITERATOR_HPP */
