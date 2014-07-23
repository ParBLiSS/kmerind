// include google test
#include <gtest/gtest.h>

// include files to test
#include <common/PackedString.hpp>
#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>
#include <common/Kmer.hpp>
#include <common/packing_iterators.hpp>
#include <common/kmer_iterators.hpp>

class PackingTest : public ::testing::Test
{
public:
  std::vector<std::string> dna_seqs = {"", "A", "T", "C", "AC", "AAAAAAAA",
    "TTTTTT", "GG", "CACACACA", "ACTGGGCCATAATCTCTCATGGATGCTACGAGCTGATCGTAGCTGACTAGTCGA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCTGCTCGTGATCGATCGATCGATCGATCGATCGACTAGCTAGCTACGTACGTACGTACGATCGATCGTACGATCGATCGATCGATCGATCGTACGTACGTACGTCAGCTAGCTCGATCGATATGCTAGATACACACACTACATACTCACATACACATCAACTCATACTTCATCAGACGATCGATCGATCTCGCTGATCGACTGTCGATCGTCGATCGTAGCTCGTAGCTCGTAGCTCG"
  };

  std::vector<std::string> huge_dna_seqs;

  virtual void SetUp()
  {
    // fill huge array with random sequences
    /*
    std::stringstream ss;
    unsigned int dna_huge_length = 10000000;
    for (unsigned int i = 0; i < dna_huge_length; ++i)
    {
      // TODO randomize?
      ss << 'C';
    }
    dna_seqs.push_back(ss.str());
    */
  }
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

    // define a packing iterator to wrap around the string's iterators
    typedef bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packit_t;
    packit_t packIt(dna.begin(), dna.end());
    packit_t packItEnd(dna.end());
    // define an unpacking iterator to wrap around the packing iterator
    typedef bliss::UnpackingIterator<packit_t, bliss::AlphabetTraits<DNA>::getBitsPerChar()> unpackit_t;
    unpackit_t unpackIt(packIt);
    // pass the total unpacked size, so that the end is properly defined
    unpackit_t unpackEnd(packIt, dna.size());

    // simply pass the iterators to the constructor of a vector:
    //std::vector<char> unpacked_dna(unpackIt, unpackEnd);
    // pre allocate to total size
    std::vector<char> unpacked_dna(dna.size());

    std::copy(unpackIt, unpackEnd, unpacked_dna.begin());
    //auto outit = unpacked_dna.begin();
    //while(unpackIt != unpackEnd)
    //{
    //  *(outit++) = *(unpackIt++);
    //}

    // compare the vector and the given string
    for (unsigned int i = 0; i < dna.size(); ++i)
    {
      EXPECT_EQ(dna[i], unpacked_dna[i]) << "The unpacked char was not the same as the original char.";
    }
  }
}

// TODO: put this somewhere else
template<typename T>
std::string getTypeName()
{
  // NOTE: this is compiler dependent (works with gcc and clang though)
  std::string function_type(__PRETTY_FUNCTION__);
  std::size_t eq = function_type.find('=');
  std::size_t sem = function_type.find(';', eq);
  return function_type.substr(eq + 2, (sem - eq - 2));
}

TEST_F(PackingTest, TestKmerGeneration1) {
  for (std::string dna : dna_seqs)
  {
    // first step: translate (in place)
    bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());
    // check that every letter's value is smaller than the total alphabet size
    for (char c : dna)
    {
      EXPECT_LT(static_cast<unsigned int>(c), bliss::AlphabetTraits<DNA>::getSize()) << "The value of the translated chars must be smaller than the size of the alphabet";
    }

    // define a packing iterator to wrap around the string's iterators
    typedef bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packit_t;
    packit_t packIt(dna.begin(), dna.end());
    packit_t packItEnd(dna.end());

    // TODO generate Kmers
    typedef bliss::Kmer<21, 2, uint8_t> Kmer;
    typedef bliss::PackedKmerGenerationIterator< packit_t, Kmer > kmer_gen_it_t;
    if (dna.size() >= 21)
    {
      kmer_gen_it_t kmerGenIt(packIt);
      kmer_gen_it_t kmerGenEnd(packIt, dna.length());


      std::vector<Kmer> kmers(kmerGenIt, kmerGenEnd);

      EXPECT_EQ(dna.size() - 21 + 1, kmers.size()) << "the number of generated kmers is incorrect";
    }



  }
}
