// include google test
#include <gtest/gtest.h>

// include files to test
#include <common/PackedString.hpp>
#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>
#include <common/PackingIterator.hpp>

class PackingTest : public ::testing::Test
{
public:
  std::vector<std::string> dna_seqs = {"", "A", "T", "C", "AC", "AAAAAAAA",
    "TTTTTT", "GG", "CACACACA", "ACTGGGCCATAATCTCTCATGGATGCTACGAGCTGATCGTAGCTGACTAGTCGA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCTGCTCGTGATCGATCGATCGATCGATCGATCGACTAGCTAGCTACGTACGTACGTACGATCGATCGTACGATCGATCGATCGATCGATCGTACGTACGTACGTCAGCTAGCTCGATCGATATGCTAGATACACACACTACATACTCACATACACATCAACTCATACTTCATCAGACGATCGATCGATCTCGCTGATCGACTGTCGATCGTCGATCGTAGCTCGTAGCTCGTAGCTCG"
  };
};




TEST_F(PackingTest, TestPackedString1) {
  for (std::string dna : dna_seqs)
  {
    // first step: translate (in place)
    bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());
    // check that every letter's value is smaller than the total alphabet size
    for (char c : dna)
    {
      EXPECT_LT(static_cast<unsigned int>(c), bliss::AlphabetTraits<DNA>::getSize()) << "The value of the translated chars must be smaller than the size of the alphabet";
    }

    // then pack it
    bliss::PackedString<DNA> packedStr(dna);

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
}

TEST_F(PackingTest, TestPackedString2) {
  for (std::string dna : dna_seqs)
  {
    // translate into different container
    std::vector<char> dna_vec(dna.size());
    bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna_vec.begin());
    // check that every letter's value is smaller than the total alphabet size
    for (char c : dna_vec)
    {
      EXPECT_LT(static_cast<unsigned int>(c), bliss::AlphabetTraits<DNA>::getSize()) << "The value of the translated chars must be smaller than the size of the alphabet";
    }

    // then pack it
    bliss::PackedString<DNA> packedStr(dna_vec.begin(), dna_vec.end());

    ASSERT_EQ(packedStr.size(), dna.size()) << "The packed sequence should be of same length as the original string.";
    ASSERT_EQ(packedStr.size(), dna_vec.size()) << "The packed sequence should be of same length as the original string.";

    // unpack it into a basic_string<DNA> container
    std::basic_string<DNA> unpacked_dna(packedStr.size(),'X');
    packedStr.unpackSequence(unpacked_dna.begin());

    // compare the values using the [] operator
    for (unsigned int i = 0; i < packedStr.size(); ++i)
    {
      EXPECT_EQ(packedStr[i], unpacked_dna[i]) << "packed and unpacked strings have to compare equal";
    }

    // translate back to ASCII (in place)
    bliss::AlphabetTraits<DNA>::translateToAscii(unpacked_dna.begin(), unpacked_dna.end(), unpacked_dna.begin());

    // compare the vector and the given string
    for (unsigned int i = 0; i < dna.size(); ++i)
    {
      EXPECT_EQ(dna[i], unpacked_dna[i]) << "The unpacked char was not the same as the original char.";
    }
  }
}


TEST_F(PackingTest, TestPackingIterator) {
  for (std::string dna : dna_seqs)
  {
    // first step: translate (in place)
    bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());
    // check that every letter's value is smaller than the total alphabet size
    for (char c : dna)
    {
      EXPECT_LT(static_cast<unsigned int>(c), bliss::AlphabetTraits<DNA>::getSize()) << "The value of the translated chars must be smaller than the size of the alphabet";
    }

    bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packIt(dna.begin(), dna.end());
    bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packItEnd(dna.end(), dna.end());


    std::vector<WordType> packedStr2(packIt, packItEnd);
    // TODO continue testing the iterators

    // then pack it
    bliss::PackedString<DNA> packedStr(dna);

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
}