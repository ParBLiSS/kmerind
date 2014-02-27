// include google test
#include <gtest/gtest.h>

// include C stdlib
#include <cstdint>

// include files to test
#include <common/padding.hpp>



TEST(PaddingTest, RemovePaddingSameType)
{
  uint32_t paddedseq[4] = {0x0000dead, 0x0000beef,0x00001234,0x0000abba};
  uint32_t unpadded[2];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 4*16, 16);

  EXPECT_EQ(unpadded[0], 0xbeefdead) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0xabba1234) << std::hex << "wrong result is: 0x"  << unpadded[1];
}

TEST(PaddingTest, RemovePaddingSameType2)
{
  uint32_t paddedseq[6] = {0x0012dead, 0x0034beef,0xffff1234,0x0056abba, 0x00c0ffee, 0x00decaff};
  uint32_t unpadded[(6*6)/8 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 6*6*4, 8);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xef12dead) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0x123434be) << std::hex << "wrong result is: 0x"  << unpadded[1];
  EXPECT_EQ(unpadded[2], 0x56abbaff) << std::hex << "wrong result is: 0x"  << unpadded[2];
  EXPECT_EQ(unpadded[3], 0xffc0ffee) << std::hex << "wrong result is: 0x"  << unpadded[3];
  EXPECT_EQ(unpadded[4], 0x0000deca) << std::hex << "wrong result is: 0x"  << unpadded[4];
}



TEST(PaddingTest, RemovePaddingSameType3)
{
  uint32_t paddedseq[4] = {0x0000dead, 0x000beef,0xffff1234,0x0056abba};
  uint32_t unpadded[(13*4)/32 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 4*13, 32-13);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xd3ddfead) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0x0005dd48) << std::hex << "wrong result is: 0x"  << unpadded[1];
}

TEST(PaddingTest, RemovePadding32to16by8)
{
  uint32_t paddedseq[6] = {0x0012dead, 0x0034beef,0xffff1234,0x0056abba, 0x00c0ffee, 0x00decaff};
  uint16_t unpadded[(6*6*4)/16 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 6*6*4, 8);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xdead) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0xef12) << std::hex << "wrong result is: 0x"  << unpadded[1];
  EXPECT_EQ(unpadded[2], 0x34be) << std::hex << "wrong result is: 0x"  << unpadded[2];
  EXPECT_EQ(unpadded[3], 0x1234) << std::hex << "wrong result is: 0x"  << unpadded[3];
  EXPECT_EQ(unpadded[4], 0xbaff) << std::hex << "wrong result is: 0x"  << unpadded[4];
  EXPECT_EQ(unpadded[5], 0x56ab) << std::hex << "wrong result is: 0x"  << unpadded[5];
  EXPECT_EQ(unpadded[6], 0xffee) << std::hex << "wrong result is: 0x"  << unpadded[6];
  EXPECT_EQ(unpadded[7], 0xffc0) << std::hex << "wrong result is: 0x"  << unpadded[7];
  EXPECT_EQ(unpadded[8], 0xdeca) << std::hex << "wrong result is: 0x"  << unpadded[8];
}

TEST(PaddingTest, RemovePadding32to8by28)
{
  uint32_t paddedseq[6] = {0x0012dead, 0x0034beef,0xffff1234,0x0056abba, 0x00c0ffee, 0x00decaff};
  uint8_t unpadded[(4*6)/8 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 4*6, 28);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xfd) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0xa4) << std::hex << "wrong result is: 0x"  << unpadded[1];
  EXPECT_EQ(unpadded[2], 0xfe) << std::hex << "wrong result is: 0x"  << unpadded[2];
}

TEST(PaddingTest, RemovePadding32to16by13)
{
  uint32_t paddedseq[4] = {0x0000dead, 0x000beef,0xffff1234,0x0056abba};
  uint16_t unpadded[(13*4)/16 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 4*13, 32-13);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xfead) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0xd3dd) << std::hex << "wrong result is: 0x"  << unpadded[1];
  EXPECT_EQ(unpadded[2], 0xdd48) << std::hex << "wrong result is: 0x"  << unpadded[2];
  EXPECT_EQ(unpadded[3], 0x0005) << std::hex << "wrong result is: 0x"  << unpadded[3];
}

TEST(PaddingTest, RemovePadding32to8by13)
{
  uint32_t paddedseq[4] = {0x0000dead, 0x000beef,0xffff1234,0x0056abba};
  uint8_t unpadded[(13*4)/8 + 1];

  uint32_t* begin = paddedseq;
  bliss::removePadding(begin, unpadded, 4*13, 32-13);

  // check for correct depadding (entered manually)
  EXPECT_EQ(unpadded[0], 0xad) << std::hex << "wrong result is: 0x"  << unpadded[0];
  EXPECT_EQ(unpadded[1], 0xfe) << std::hex << "wrong result is: 0x"  << unpadded[1];
  EXPECT_EQ(unpadded[2], 0xdd) << std::hex << "wrong result is: 0x"  << unpadded[2];
  EXPECT_EQ(unpadded[3], 0xd3) << std::hex << "wrong result is: 0x"  << unpadded[3];
  EXPECT_EQ(unpadded[4], 0x48) << std::hex << "wrong result is: 0x"  << unpadded[4];
  EXPECT_EQ(unpadded[5], 0xdd) << std::hex << "wrong result is: 0x"  << unpadded[5];
  EXPECT_EQ(unpadded[6], 0x05) << std::hex << "wrong result is: 0x"  << unpadded[6];
}
