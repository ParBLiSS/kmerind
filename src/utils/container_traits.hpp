/**
 * @file	container_traits.hpp
 * @ingroup bliss::utils
 * @author	Tony Pan
 * @brief   helper functions to check if a container has the required methods.
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef CONTAINER_TRAITS_HPP_
#define CONTAINER_TRAITS_HPP_

#include <utility>
#include <type_traits>

    // preprocessor macro to save some typing
#define BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(NAME) \
        template<typename CC = C> \
        static constexpr auto has_##NAME(CC*) -> decltype(std::declval<CC>().NAME(), std::true_type()); \
        template<typename> \
        static constexpr std::false_type has_##NAME(...);


namespace bliss {
  namespace utils {

    /**
     * @class container_traits
     * @brief This class is used to determine if certain methods are present.
     * @tparam C    container type, whose traits this class extracts.
     */
    template<typename C>
    class container_traits {

      private:

        /**
         * @brief   determine if the container has an begin method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last (true/false).
         *
         *          method is defined only if begin exists.
         *          adopted from  http://dev.krzaq.cc/checking-whether-a-class-has-a-member-function-with-a-given-signature/
         *
         * @tparam  CC   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(begin);

        /**
         * @brief   determine if a type has an end method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last (true/false).
         *
         *          method is defined only if end exists.
         *          adopted from  http://dev.krzaq.cc/checking-whether-a-class-has-a-member-function-with-a-given-signature/
         *
         * @tparam  CC   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if end method is defined, for all else return false.
         */
        BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(end);

        /**
         * @brief   determine if a type has an cbegin method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last (true/false).
         *
         *          method is defined only if cbegin exists.
         *          adopted from  http://dev.krzaq.cc/checking-whether-a-class-has-a-member-function-with-a-given-signature/
         *
         * @tparam  CC   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(cbegin);

        /**
         * @brief   determine if a type has an cend method
         * @detail  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last (true/false).
         *
         *          method is defined only if cend exists.
         *          adopted from  http://dev.krzaq.cc/checking-whether-a-class-has-a-member-function-with-a-given-signature/
         *
         * @tparam  CC   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        BLISS_UTILS_CONTAINER_HAS_METHOD_WITH_NAME(cend);


        /**
         * @brief   determine if the container has an assign method
         * @details  uses declval to get lvalues, and use decltype and SFINAE to do "enable_if".
         *          decltype evaluates all expressions, but returns type of last (true/false).
         *
         *          method is defined only if assign exists.
         *          test with the container's own iterator, so also depend on begin() and end() existing.
         *
         *          adopted from  http://dev.krzaq.cc/checking-whether-a-class-has-a-member-function-with-a-given-signature/
         *
         * @tparam  CC   type to check, should be a container.
         * @param       unamed pointer parameter for automatic type deduction for the templates
         * @return      true if method is defined, for all else return false.
         */
        template<typename CC = C>
        static constexpr auto has_assign(CC*) ->
            decltype(std::declval<CC>().assign(std::declval<decltype(std::declval<CC>().begin())>(),
                                               std::declval<decltype(std::declval<CC>().end())>()
                                               ), std::true_type());
        template<typename>
        static constexpr std::false_type has_assign(...);

      public:

        /**
         * @var       isConstIterable
         * @brief     static bool indicating if container supports constant iteration (has cbegin and cend)
         */
        static constexpr bool isConstIterable = std::is_same<decltype(has_cbegin<C>(nullptr)), std::true_type>::value
                  && std::is_same<decltype(has_cend<C>(nullptr)), std::true_type>::value;


        /**
         * @var       isIterable
         * @brief     static bool indicating if container supports iteration (has begin and end)
         */
        static constexpr bool isIterable = std::is_same<decltype(has_begin<C>(nullptr)), std::true_type>::value
            && std::is_same<decltype(has_end<C>(nullptr)), std::true_type>::value;


        /**
         * @var       isAssignable
         * @brief     static bool indicating if container supports assignment.
         */
        static constexpr bool isAssignable = std::is_same<decltype(has_assign<C>(nullptr)), std::true_type>::value;

    };


  } //namespace utils
} //namespace bliss


#endif /* CONTAINER_TRAITS_HPP_ */
