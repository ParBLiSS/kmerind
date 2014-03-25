/**
 * @file    compression_iterator.hpp
 * @ingroup interators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the compressing iterator.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */


#ifndef BLISS_ITERATORS_COMPRESSING_ITERATOR_HPP
#define BLISS_ITERATORS_COMPRESSING_ITERATOR_HPP

#include <iterator>

#include <common/bit_ops.hpp>
#include <iterators/iterator_tools.hpp>
#include <iterators/function_traits.hpp>

namespace bliss
{
namespace iterator
{


template<typename Functor, typename Iterator>
class _shared_compressing_iterator : public std::iterator<
    /* inherit interator traits from std::iterator: */
    // 1) iterator type tag (same as base)
    typename std::iterator_traits<Iterator>::iterator_category,
    // 2) value type = return type of compression function
    typename std::remove_reference<
      typename bliss::functional::function_traits<Functor,
        Iterator,  Iterator>::return_type>::type,
    // 3) difference_type = same as base iterator
    typename std::iterator_traits<Iterator>::difference_type>
{
protected:
  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<_shared_compressing_iterator> traits;

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
   * @brief     Returns the compression functor.
   *
   * @return    The compression functor.
   */
  Functor& getCompressor()
  {
    return _f;
  }

  /**
   * @brief     Returns the compression functor.
   *
   * @return    The compression functor.
   */
  const Functor& getCompressor() const
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
  inline bool operator==(const _shared_compressing_iterator& rhs) const
  {
    return _base == rhs._base;
  }

  /**
   * @brief     Returns whether this and the given iterator are different.
   *
   * @return    Whether this and the given iterator are different.
   */
  inline bool operator!=(const _shared_compressing_iterator& rhs) const
  {
    return _base != rhs._base;
  }

protected:
  /**********************
   *  Member variables  *
   **********************/

  /// The current base iterator position
  Iterator _base;

  // The end of the input sequence
  Iterator _end;

  /// The compressor function
  const Functor _f;

  /// the number of base elements per compressed element (i.e the `m` in m->1)
  const difference_type _m;

  /**************************************************
   *  Private constructor (no direct construction)  *
   **************************************************/
  _shared_compressing_iterator() {}

  /******************
   *  Constructors  *
   ******************/

  /**
   * @brief Constructor, taking the base iterator and the compressing
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The compressing functor.
   * @param m           The number of elements combined into one.
   */
  _shared_compressing_iterator(const Iterator& base_iter, const Iterator& end_iter,
                           const Functor& f, const difference_type m)
      : _base(base_iter), _end(end_iter), _f(f), _m(m)
  {
  }
};



template<typename Compressor, typename BaseIterator, typename DerivedIterator>
// TODO: template specialization for RA vs FWD
class _compressing_iterator_ra : _shared_compressing_iterator<Compressor, BaseIterator>
{
protected:
  // the base iterator traits
  typedef typename std::iterator_traits<BaseIterator> base_traits;

  // the base class
  typedef _shared_compressing_iterator<Compressor, BaseIterator> base_class_type;

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
   * @brief Constructor, taking the base iterator and the compressing
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The compressing functor.
   * @param m           The number of elements combined into one.
   */
  _compressing_iterator_ra(const BaseIterator& base_iter, const BaseIterator& end_iter,
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
    return *this;
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
      iter_tools<BaseIterator>::advance_at_most(this->_base, n * this->_m, this->_end);
    return *this;
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
    return *this;
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
    type tmp(*this);
    this->operator++();
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
    type tmp(*this);
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
    type output(*this);
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
    type output(*this);
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

template<typename Compressor, typename BaseIterator, typename DerivedIterator>
// TODO: template specialization for RA vs FWD
class _compressing_iterator_dir : _shared_compressing_iterator<Compressor, BaseIterator>
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
  typedef _shared_compressing_iterator<Compressor, BaseIterator> base_class_type;

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
   * @brief Constructor, taking the base iterator and the compressing
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The compressing functor.
   * @param m           The number of elements combined into one.
   */
  _compressing_iterator_dir(const BaseIterator& base_iter, const BaseIterator& end_iter,
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
    return *this;
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
  typename std::enable_if<is_bidir, type>::type&
  operator--()
  {
    // go back `m`, cast to signed int before negating
    // TODO: what happens if we are at an uneven end?
    std::advance(this->_base, - static_cast<int>(this->_m));
    return *this;
  }

  /**
   * @brief     Post-decrement operator.
   *
   * Reduces this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-decremented iterator.
   */
  typename std::enable_if<is_bidir, type>::type
  operator--(int)
  {
    // create a copy
    type tmp(*this);
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
    type tmp(*this);
    this->operator++();
    return tmp;
  }
};


/**
 * @brief Compressing Iterator class (m -> 1 mapping)
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
 * @tparam Compressor   The type of the compressing functor.
 * @tparam Iterator     The type of the base iterator.
 */
template<typename Compressor, typename Iterator>
class compressing_iterator : std::conditional<std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                                std::random_access_iterator_tag>::value,
                                _compressing_iterator_ra<Compressor, Iterator, compressing_iterator<Compressor, Iterator> >,
                                _compressing_iterator_dir<Compressor, Iterator, compressing_iterator<Compressor, Iterator> > >::type
{
protected:

  /***********************
   *  Private type defs  *
   ***********************/

  // the base iterator traits
  typedef std::iterator_traits<Iterator> base_traits;

  typedef typename base_traits::difference_type difference_type;

  // the traits (e.g. return type) of the compression function
  typedef bliss::functional::function_traits<Compressor,
      Iterator, Iterator> functor_traits;
  /// The type of this class
  typedef compressing_iterator type;
  /// The type of the base class
  typedef typename std::conditional<std::is_same<typename base_traits::iterator_category,
            std::random_access_iterator_tag>::value,
            _compressing_iterator_ra<Compressor, Iterator, compressing_iterator>,
            _compressing_iterator_dir<Compressor, Iterator, compressing_iterator> >::type base_class_type;



protected:
  /**************************************************
   *  Helper expressions for the iterator category  *
   **************************************************/

  /*
  /// Whether the iterator is an InputIterator
  static constexpr bool is_input = std::is_same<iterator_category,
                                std::input_iterator_tag>::value;
  /// Whether the iterator is a ForwardIterator
  static constexpr bool is_fwd = std::is_same<iterator_category,
                              std::forward_iterator_tag>::value;
  /// Whether the iterator is a BidirectionalIterator
  static constexpr bool is_bidir = std::is_same<iterator_category,
                                std::bidirectional_iterator_tag>::value;
  /// Whether the iterator is a RandomAccessIterator
  static constexpr bool is_ra = std::is_same<iterator_category,
                                std::random_access_iterator_tag>::value;

  /// Whether the iterator is at least a ForwardIterator
  static constexpr bool is_min_fwd = is_fwd || is_bidir || is_ra;
  /// Whether the iterator is at least a BidirectionalIterator
  static constexpr bool is_min_bidir = is_bidir || is_ra;
  */


public:

  /******************
   *  Constructors  *
   ******************/

  // class specific constructor
  /**
   * @brief Constructor, taking the base iterator and the compressing
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param end_iter    The end of the input sequence. This is needed in case
   *                    the length of the input sequence is not perfectly
   *                    dividable by m.
   * @param f           The compressing functor.
   * @param m           The number of elements combined into one.
   */
  compressing_iterator(const Iterator& base_iter, const Iterator& end_iter,
                       const Compressor& f, const difference_type m)
    : base_class_type(base_iter, end_iter, f, m)
  {
  }

  /*
  // specific to at least fwd iterators: default contructable
  template< typename = std::enable_if<is_min_fwd && !is_ra> >
  compressing_iterator()
      :  _base(), _next(), _end(), _f(), _m(0)
  {
  }
  */

  /*
  // specific to RandomAccessIterators
  template< typename = std::enable_if<is_ra> >
  compressing_iterator()
      :  _base(), _end(), _f(), _m(0)
  {
  }
  */

};

} // iterator
} // bliss
#endif /* BLISS_ITERATORS_COMPRESSING_ITERATOR_HPP */
