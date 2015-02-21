//
// metatest.cpp: This program produces a nice set of
//				 errors to flex the muscles of STLFilt's
//				 template metaprogramming-specific features.
//
// Note: You must have a recent version of the Boost libs installed.

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/reverse_fold.hpp>
#include <boost/mpl/plus.hpp>

namespace mpl = boost::mpl;

mpl::reverse_fold<
        mpl::range_c<long, 0, 50>
      , void
      , mpl::plus<>
>::type x;
