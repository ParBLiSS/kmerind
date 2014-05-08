/**
 * @file		template_helper.hpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef TEMPLATE_HELPER_HPP_
#define TEMPLATE_HELPER_HPP_

// to support the "curiously recursive template pattern "static polymorphism"
// usage:  Declare_Base(X) { basic() };  Declare_Base_CTRP(X) { interface() };  Derive_Base_CTRP(Y, X) { impl()  };
//    see http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
#define Declare_Base(BASE)      class BASE
#define Declare_Base_CTRP(BASE) class BASE_CTRP : public BASE
#define Derive_Base_CTRP(Type, BASE)  class Type : public BASE_CTRP


#endif /* TEMPLATE_HELPER_HPP_ */
