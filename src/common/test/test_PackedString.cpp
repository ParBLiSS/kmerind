// include google test
#include <gtest/gtest.h>

// include files to test
#include <common/PackedString.hpp>
#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>

TEST(BlissCommonSuite, TestPackedString1) {
  std::string dna = "ACTGGGCCATAATCTCTCATGGATGCTACGAGCTGATCGTAGCTGACTAGTCGA";
  // first step: translate (in place)
  AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());
  // check that every letter's value is smaller than the total alphabet size
  for (char c : dna)
  {
    EXPECT_LT(static_cast<unsigned int>(c), AlphabetTraits<DNA>::getSize());
  }

  // then pack it
  PackedString<DNA> packedStr(dna);

  ASSERT_EQ(packedStr.size(), dna.size()) << "The packed sequence should be of same length as the original string.";

  // unpack it
  std::vector<DNA> unpacked_dna(packedStr.size());
  packedStr.unpackSequence(unpacked_dna.begin());

  // compare the vector and the given string
  for (unsigned int i = 0; i < dna.size(); ++i)
  {
    EXPECT_EQ(dna[i], unpacked_dna[i]) << "The unpacked char was not the same as the original char.";
  }
}