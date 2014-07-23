/**
 * @file		container_traits.hpp
 * @ingroup bliss::utils
 * @author	tpan
 * @brief   helper functions to check if a container has the required methods.
 * @details  adopted from  http://stackoverflow.com/questions/9530928/checking-a-member-exists-possibly-in-a-base-class-c11-version
 *            updated to use macro to generate the simple cases.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONTAINER_TRAITS_HPP_
#define CONTAINER_TRAITS_HPP_

#include <utility>
#include <type_traits>

    // preprocessor directive to save some typing.
#define BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(NAME) \
        template<typename CC = C> \
        static auto test_##NAME##(CC*) -> decltype(std::declval<CC>().##NAME##(), std::true_type()); \
        template<typename> \
        static std::false_type test_##NAME##(...);


namespace bliss {
  namespace utils {

    /**
     * container traits class, templated to container's data type.
     */
    template<typename C>
    class container_traits {

      private:
        BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(begin);


      public:
        /**
         * @brief   determine if the container has an begin method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last.
         *
         *          method is defined only if begin exists.
         *
         * @tparam  T   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        static constexpr bool hasBeginMethod = std::is_same<decltype(has_begin<C>(nullptr)), std::true_type>::value;

        /**
         * @brief   determine if a type has an end method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last.
         *
         *          method is defined only if end exists.
         *
         * @tparam  T   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if end method is defined, for all else return false.
         */
        //BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(end);
        static constexpr bool hasEndMethod = !std::is_same<decltype(std::declval<C>().end()), void>::value;

        /**
         * @brief   determine if a type has an cbegin method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last.
         *
         *          method is defined only if cbegin exists.
         *
         * @tparam  T   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        //BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(cbegin);
        static constexpr bool hasCBeginMethod = !std::is_same<decltype(std::declval<C>().cbegin()), void>::value;

        /**
         * @brief   determine if a type has an cend method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last.
         *
         *          method is defined only if cend exists.
         *
         * @tparam  T   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        //BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(cend);
        static constexpr bool hasCEndMethod = !std::is_same<decltype(std::declval<C>().cend()), void>::value;
        /**
         * @brief     check to see if container supports constant iteration.
         * @tparam T  container type to check
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return    true if container has cbegin and cend.
         */
        template<typename CC = C>
        static constexpr auto is_const_iterable(CC*) -> decltype(std::declval<CC>().cend(), std::declval<C>().cbegin(), bool())
        {
          return true;
        }
        template<typename CC = C>
        static constexpr auto is_const_iterable(CC*) -> typename std::enable_if<std::is_same<decltype(std::declval<CC>().cend()), void>::value, bool>::type
        {
          return false;
        }

        /**
         * @brief     check to see if container supports iteration.
         * @tparam T  container type to check
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return    true if container has begin and end.
         */
        static constexpr bool is_iterable() {
          return hasBeginMethod && hasEndMethod;
        }


        /**
         * @brief   determine if a type has an assign method that takes start and end iterators
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last.
         *
         *          method is defined only if assign exists.
         *
         *          adopted from  http://stackoverflow.com/questions/9530928/checking-a-member-exists-possibly-in-a-base-class-c11-version
         *
         * @tparam  T              type to check, should be a container.
         * @tparam  InputIterator  type of iterator to use with "assign"
         * @param       unamed container pointer parameter for automatic type deduction for the templates
         * @param       unamed input iterator pointer parameter for automatic type deduction for the templates
         * @return  true (if method is defined, which is when T has "assign"
         */
        template <typename InputIterator>
        static constexpr bool is_assignable()
        {
          return !std::is_same<decltype(std::declval<C>().assign(std::declval<InputIterator>(), std::declval<InputIterator>())), void>::value;
        }




    };

    template<typename C>
    constexpr bool container_traits<C>::hasEndMethod;
    template<typename C>
    constexpr bool container_traits<C>::hasCBeginMethod;
    template<typename C>
    constexpr bool container_traits<C>::hasCEndMethod;

  } //namespace utils
} //namespace bliss


#endif /* CONTAINER_TRAITS_HPP_ */
