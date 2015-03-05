/**
 * @file    sliding_window_iterator.hpp
 * @ingroup iterators
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the sliding window iterators.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */


#ifndef BLISS_ITERATORS_SLIDING_WINDOW_ITERATOR_HPP
#define BLISS_ITERATORS_SLIDING_WINDOW_ITERATOR_HPP

#include <iterator>

#include "iterators/iterator_utils.hpp"

namespace bliss
{
namespace iterator
{

template<typename BaseIterator, typename Window>
class sliding_window_iterator : public std::iterator<
    // 1) iterator type tag: at most forward, but if the base iterator is less,
    //    then inherit from it
    typename iterator_tags<BaseIterator>::inherit_forward,
    // 2) value type = return type of the getValue function in the Window class
    typename std::remove_reference<decltype(std::declval<Window>().getValue())>::type,
    // 3) difference_type = same as base iterator
    typename std::iterator_traits<BaseIterator>::difference_type>
{
protected:
  // the base iterator traits
  typedef std::iterator_traits<BaseIterator> base_traits;

  /// This type
  typedef sliding_window_iterator type;

protected:
  /**********************
   *  Member variables  *
   **********************/

  // input iterators, _base is used for comparison and represents the current
  // position, while _next reads ahead in order to fill the sliding window.
  BaseIterator _base;
  BaseIterator _next;

  // the window
  Window _window;

  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<type> traits;

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


  /**********************
   *  Member accessors  *
   **********************/

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  BaseIterator& getBaseIterator()
  {
    return _base;
  }

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  const BaseIterator& getBaseIterator() const
  {
    return _base;
  }

  /**
   * @brief Returns a reference to the current window structure.
   *
   * @return The current sliding window structure.
   */
  Window& getWindow()
  {
    return _window;
  }

  /**
   * @brief Returns a reference to the current window structure.
   *
   * @return The current sliding window structure.
   */
  const Window& getWindow() const
  {
    return _window;
  }


  /*****************
   *  Comparators  *
   *****************/

  /**
   * @brief     Returns whether this and the given iterator point to the same
   *            positon.
   *
   * @return    Whether this and the given iterator point to the same position.
   */
  inline bool operator==(const sliding_window_iterator& rhs) const
  {
    return _base == rhs._base;
  }

  /**
   * @brief     Returns whether this and the given iterator are different.
   *
   * @return    Whether this and the given iterator are different.
   */
  inline bool operator!=(const sliding_window_iterator& rhs) const
  {
    return _base != rhs._base;
  }

public:
  /******************
   *  Constructors  *
   ******************/

  /// Default contructor
  sliding_window_iterator() = default;

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param window      The sliding window structure
   */
  sliding_window_iterator(const BaseIterator& base_iter, const Window& window, bool init_first=true)
      : _base(base_iter), _next(base_iter), _window(window)
  {
    if (init_first)
    {
      do_init();
    }
  }

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   */
  sliding_window_iterator(const BaseIterator& base_iter, bool init_first=true)
      : _base(base_iter), _next(base_iter), _window()
  {
    if (init_first)
    {
      do_init();
    }
  }

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This returns the buffered value of the current sliding window position.
   *
   * @return    The value at the current position.
   */
  inline value_type operator*()
  {
    // in case the current element has not been read yet: do so now
    if (this->_base == this->_next)
    {
      // add the current element to the sliding window, this iterates the _next
      // pointer by one
      this->_window.next(this->_next);
    }
    // the result was already computed and stored in the functor
    return this->_window.getValue();
  }

  /**
   * @brief     Pre-increment operator.
   *
   * Advances this iterator by one position.
   *
   * @return    A reference to this.
   */
  inline type& operator++()
  {
    // if we have not read the current element (thus iterating one postion
    // ahead with _next), then do it now (this happens of operator*() is not
    // called)
    if (this->_base == this->_next)
    {
      // slide the window, inserting the current element and iterating ahead
      // with _next
      this->_window.next(this->_next);
    }

    // set both iterators to the same position
    this->_base = this->_next;

    // return a reference to this
    return *this;
  }

  /**
   * @brief     Post-increment operator.
   *
   * Advances this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-incremented iterator.
   */
  inline type operator++(int)
  {
    // create a copy
    type tmp(*this);
    this->operator++();
    return tmp;
  }
protected:
  /**
   * @brief Initializes the first window.
   */
  inline void do_init()
  {
    this->_window.init(this->_base);
    this->_next = this->_base;
    ++this->_next;
  }
};



template<typename BaseIterator, typename Window>
class one2many_sliding_window_iterator : public std::iterator<
    // 1) iterator type tag: at most forward, but if the base iterator is less,
    //    then inherit from it
    typename iterator_tags<BaseIterator>::inherit_forward,
    // 2) value type = return type of the getValue function in the Window class
    typename std::remove_reference<decltype(std::declval<Window>().getValue())>::type,
    // 3) difference_type = same as base iterator
    typename std::iterator_traits<BaseIterator>::difference_type>
{
protected:
  // the base iterator traits
  typedef std::iterator_traits<BaseIterator> base_traits;

  /// This type
  typedef one2many_sliding_window_iterator type;

  /// The std::iterator_traits of this iterator
  typedef std::iterator_traits<type> traits;

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
  /**********************
   *  Member variables  *
   **********************/

  // input iterators, _base is used for comparison and represents the current
  // position, while _next reads ahead in order to fill the sliding window.
  BaseIterator _base;
  BaseIterator _next;

  // the offsets in each word (many2one)
  difference_type _base_offset;
  difference_type _next_offset;

  // the window
  Window _window;

public:

  /**********************
   *  Member accessors  *
   **********************/

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  BaseIterator& getBaseIterator()
  {
    return _base;
  }

  /**
   * @brief     Returns the current base iterator.
   *
   * @return    The current base iterator.
   */
  const BaseIterator& getBaseIterator() const
  {
    return _base;
  }

  /**
   * @brief Returns the current offset within the current iterator position.
   *
   * @return the current offset.
   */
  difference_type getOffset() const
  {
    return _base_offset;
  }

  /**
   * @brief Returns a reference to the current window structure.
   *
   * @return The current sliding window structure.
   */
  Window& getWindow()
  {
    return _window;
  }

  /**
   * @brief Returns a reference to the current window structure.
   *
   * @return The current sliding window structure.
   */
  const Window& getWindow() const
  {
    return _window;
  }


  /*****************
   *  Comparators  *
   *****************/

  /**
   * @brief     Returns whether this and the given iterator point to the same
   *            positon.
   *
   * @return    Whether this and the given iterator point to the same position.
   */
  inline bool operator==(const one2many_sliding_window_iterator& rhs) const
  {
    return _base == rhs._base && _base_offset == rhs._base_offset;
  }

  /**
   * @brief     Returns whether this and the given iterator are different.
   *
   * @return    Whether this and the given iterator are different.
   */
  inline bool operator!=(const one2many_sliding_window_iterator& rhs) const
  {
    return !this->operator==(rhs);
  }

public:
  /******************
   *  Constructors  *
   ******************/

  /// default constructor
  one2many_sliding_window_iterator()
    :_base(), _next(), _base_offset(0), _next_offset(0), _window()
  {
  }

  /// constructor, initialize with a base iterator.  optionally initialize the first sliding window result
  one2many_sliding_window_iterator(const BaseIterator& base_iter, bool init_first = true)
    :_base(base_iter), _next(base_iter), _base_offset(0), _next_offset(0), _window()
  {
    if (init_first)
    {
      do_init();
    }
  }

  /// constructor, initialize the sliding window to offset.
  one2many_sliding_window_iterator(const BaseIterator& base_iter, difference_type offset)
    :_base(base_iter), _next(base_iter), _base_offset(0), _next_offset(0), _window()
  {
    this->advance(offset);
  }

  /**
   * @brief Constructor, taking the base iterator and the many2one
   *        functor.
   *
   * @param base_iter   The base iterator that is wrapped via this iterator.
   * @param window      The sliding window structure
   */
  one2many_sliding_window_iterator(const BaseIterator& base_iter, const Window& window, bool init_first = true)
      : _base(base_iter), _next(base_iter),
        _base_offset(0), _next_offset(0), _window(window)
  {
    if (init_first)
    {
      do_init();
    }
  }

  one2many_sliding_window_iterator(const BaseIterator& base_iter, const Window& window, difference_type offset)
      : _base(base_iter), _next(base_iter),
        _base_offset(0), _next_offset(0), _window(window)
  {
    this->advance(offset);
  }

  /**
   * @brief     Returns the value at the current iterator position.
   *
   * This returns the buffered value of the current sliding window position.
   *
   * @return    The value at the current position.
   */
  inline value_type operator*()
  {
    // in case the current element has not been read yet: do so now
    if (this->_base == this->_next && this->_base_offset == this->_next_offset)
    {
      // add the current element to the sliding window, this sets the next
      // and offset pointers to the next readable element/offset
      this->_window.next(this->_next, this->_next_offset);
    }
    // the result was already computed and stored in the functor
    return this->_window.getValue();
  }

  /**
   * @brief     Pre-increment operator.
   *
   * Advances this iterator by one position.
   *
   * @return    A reference to this.
   */
  inline type& operator++()
  {
    // if we have not read the current element (thus iterating one postion
    // ahead with _next), then do it now (this happens of operator*() is not
    // called)
    if (this->_base == this->_next && this->_base_offset == this->_next_offset)
    {
      // add the current element to the sliding window, this sets the next
      // and offset pointers to the next readable element/offset
      this->_window.next(this->_next, this->_next_offset);
    }

    // set both iterators to the same position
    this->_base = this->_next;
    this->_base_offset = this->_next_offset;

    // return a reference to this
    return *this;
  }

  /**
   * @brief     Post-increment operator.
   *
   * Advances this iterator by one position, but returns the old iterator state.
   *
   * @return    A copy to the non-incremented iterator.
   */
  inline type operator++(int)
  {
    // create a copy
    type tmp(*this);
    this->operator++();
    return tmp;
  }

protected:
  /**
   * @brief Initializes the first window.
   */
  inline void do_init()
  {
    this->_window.init(this->_base, this->_base_offset);
    this->_next = this->_base;
    this->_next_offset = this->_base_offset;
    // make a copy of the current window to get the next iterator and offset
    Window wincpy = this->_window;
    // iterate by one, simply to get the correct `_next` pointers
    wincpy.next(this->_next, this->_next_offset);
  }

  /**
   * @brief Advances the iterator by the given number of times.
   *
   * Used for initialization, this function skips over all content and
   * invalidates the iterator for dereferencing or further advancing.
   * This is purely used to initialize an `end` style iterator that is
   * used only for comparison.
   */
  inline void advance(difference_type by)
  {
    this->_window.skip(this->_base, this->_base_offset, by);
    this->_next = this->_base;
    this->_next_offset = this->_base_offset;
    // TODO: set this iterator as non-readable?
  }
};



} // iterator
} // bliss
#endif /* BLISS_ITERATORS_SLIDING_WINDOW_ITERATOR_HPP */
