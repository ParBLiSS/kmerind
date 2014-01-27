#include <common/dummy.hpp>
#include <gtest/gtest.h>

TEST(BlissCommonSuite, Dummy1) {
    ASSERT_EQ(getSomeNumber(12, 4), 48) << "something went wrong!";
    ASSERT_EQ(getSomeNumber(12, 3), 36) << "whoops, broken test case :)";
}


TEST(BlissCommonSuite, Dummy2) {
    EXPECT_EQ(getSomeNumber(-3, -3), 9) << "something went wrong!";
    EXPECT_EQ(getSomeNumber(0, 1), 0) << "whoops, broken test case :)";
}
