#include <common/dummy.hpp>
#include <gtest/gtest.h>

TEST(MyTestSuitName, MyTestCaseName) {
    ASSERT_EQ(getSomeNumber(12, 4), 48) << "something went wrong!";
    ASSERT_EQ(getSomeNumber(12, 3), 36) << "whoops, broken test case :)";
}