#include <common/dummy.hpp>
#include <gtest/gtest.h>

TEST(BlissCommonSuite2, DummyBlah) {
    ASSERT_EQ(getSomeNumber(12, 4), 48) << "something went wrong!";
    ASSERT_EQ(getSomeNumber(12, 3), 36) << "whoops, broken test case :)";
}


TEST(BlissCommonSuite2, DummyFoo) {
    EXPECT_EQ(getSomeNumber(1, -3), -3) << "something went wrong!";
    EXPECT_EQ(getSomeNumber(0, 1), 0) << "whoops, broken test case :)";
}
