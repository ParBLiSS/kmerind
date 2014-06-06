/**
 * buffered_transform_iterator.hpp
 *
 *  Created on: Feb 5, 2014
 *      Author: tpan
 *
 *  COMPATIBLE WITH A BUFFERING FUNCTOR with
 *    operator()(Iterator curr) and
 *    TO operator()()
 *
 *  patterned after boost's transform iterator.  no multiple inheritance.  uses functor class.  allows variadic parameters
 *  support functor (possibly multiple operator() ), function pointer.
 *
 *  Does NOT support member function pointer.  for those, recommend wrap in a functor.
 *  NOTE: function pointer performs far worse then functor.  probably due to compiler optimization
 *
 *
 *  Relies on compiler to cast the input parameter to the specific types needed by the functor, including constness and references.
 *
 *     - no std::function - general feeling is that it's slower.
 *
 *  TODO:
 *  commenting.
 */

#ifndef BUFFERED_TRANSFORM_ITERATOR_HPP_
#define BUFFERED_TRANSFORM_ITERATOR_HPP_

#include <iterator>
#include "iterators/function_traits.hpp"

namespace bliss
{

  namespace iterator
  {

    /**
     * transform iterator class
     *
     * transforms each element in the list. otherwise it's the same kind of iterator as the base.
     *
     * inheriting from std::iterator ONLY to get iterator_traits support.
     *
     * careful with the use of enable_if.  false causes the "type" typedef inside not to be defined, creating an error.
     *   this may be a SFINAE usage  or implementation issue.
     */
    template<typename Transformer, typename Iterator>
    class buffered_transform_iterator :
        public std::iterator<
            /* 1.) iterator type (tag) */
            typename std::conditional<
                // if the base iterator a RandomAccess OR Bidirectional?
                std::is_same<
                    typename std::iterator_traits<Iterator>::iterator_category,
                    std::random_access_iterator_tag>::value
                || std::is_same<
                    typename std::iterator_traits<Iterator>::iterator_category,
                    std::bidirectional_iterator_tag>::value,
                // then: become a forward iterator
                std::forward_iterator_tag,
                // else: `steal` the underlying iterator type
                typename std::iterator_traits<Iterator>::iterator_category>::type,
            // the value type ( = functor return type)
            typename std::remove_reference<
                typename bliss::functional::function_traits<Transformer>::return_type>::type,
            // difference type = underlying diff type
            typename std::iterator_traits<Iterator>::difference_type>
    {
      protected:
        // define first, to avoid -Wreorder error (where the variables are initialized before buffered_transform_iterator::Transformer, etc are defined.
        typedef std::iterator_traits<Iterator> base_traits;
        typedef bliss::functional::function_traits<Transformer> functor_traits;

        Iterator _curr;
        Iterator _next;
        Transformer _f;

      public:
        typedef buffered_transform_iterator<Transformer, Iterator> type;
        typedef std::iterator_traits<type> traits;

        typedef typename base_traits::iterator_category iterator_category;
        typedef typename std::remove_reference<
        typename functor_traits::return_type>::type value_type;
        typedef typename base_traits::difference_type difference_type;
        typedef typename std::add_rvalue_reference<value_type>::type reference_type;
        typedef typename std::add_pointer<value_type>::type pointer_type;

        typedef typename base_traits::value_type base_value_type;

        // class specific constructor
        buffered_transform_iterator(const Iterator& base_iter,
                                    const Transformer & f)
            : _curr(base_iter), _next(base_iter), _f(f)
        {
        }
        ;

        // for all classes
        // no EXPLICIT so can copy assign.
        buffered_transform_iterator(const type& Other)
            : _curr(Other._curr), _next(Other._next), _f(Other._f)
        {
        }

        type& operator=(const type& Other)
        {
          _f = Other._f;
          _next = Other._next;
          _curr = Other._curr;
          return *this;
        }

        inline value_type operator*()
        {

          // special case, have not yet parsed the content
          if (_curr == _next)
          {
            // after parse, _next has been moved but curr stays where it is.
            _f(_next);
          }
          // else result was already computed and stored in the functor.

          return _f();
        }

        type& operator++()
        {

          // property:  to set _next to the appropriate position, need to parse.

          // ++ parses so that _curr moves.  _next is pushed forward.  if _next has already been moved (via *), then just move curr there.

          if (_curr == _next) // at beginning.  so need to parse, and right away
          {
            // walk the base iterator until function is done with construction or at end.
            // doing the construction here because we have to walk anyways.
            _f(_next);    // _next is moved.
          }
          // else next has already been moved (by * operator)

          _curr = _next;

          return *this;
        }

        type operator++(int)
        {
          type output(*this);
          return ++output;
        }

        // accessor functions for internal state.
        Transformer& getTransformer()
        {
          return _f;
        }
        const Transformer& getTransformer() const
        {
          return _f;
        }
        Iterator& getBaseIterator()
        {
          return _curr;
        }
        const Iterator& getBaseIterator() const
        {
          return _curr;
        }

        // input iterator specific
        inline bool operator==(const type& rhs)
        {
          return _curr == rhs._curr;
        }

        inline bool operator!=(const type& rhs)
        {
          return _curr != rhs._curr;
        }


        // no -> operator.  -> returns a pointer to the value held in iterator.
        // however, in transform iterator, we are not holding data, so pointer does not make sense.

        // NO support for output iterator

        // forward iterator
        template<
            typename = typename std::enable_if<
                std::is_same<iterator_category, std::forward_iterator_tag>::value || std::is_same<
                    iterator_category, std::bidirectional_iterator_tag>::value
                || std::is_same<iterator_category,
                    std::random_access_iterator_tag>::value>::type >
        explicit buffered_transform_iterator()
            : _f(), _next()
        {
        }
        ;

    };

  } // iterator
} // bliss
#endif /* BUFFERED_BUFFERED_TRANSFORM_ITERATOR_HPP_ */
