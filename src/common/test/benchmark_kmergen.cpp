#include <gtest/gtest.h>

#include <chrono>

#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>
#include <common/PackingIterator.hpp>
#include <utils/generator.hpp>




TEST(Benchmark_KmerGeneration, BenchmarkKmer1)
{
  // generate a random piece of DNA
  std::string dna = bliss::utils::random_dna(100000000);

  auto start = std::chrono::high_resolution_clock::now();
  // first step: translate (in place)
  // TODO: do this as transformation iterator
  bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());
  auto stop = std::chrono::high_resolution_clock::now();
  std::cerr << "Duration of translation: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << "ms" << std::endl;

  // define a packing iterator to wrap around the string's iterators
  typedef bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packit_t;
  packit_t packIt(dna.begin(), dna.end());
  packit_t packItEnd(dna.end());

  /* Benchmark packing */

  start = std::chrono::high_resolution_clock::now();
  std::array<packit_t::value_type, 1000> packarr;
  packit_t begin = packIt;
  unsigned int j = 0;
  while (begin != packItEnd)
  {
    packarr[j % 1000] = *begin;
    ++begin;
    ++j;
  }
  stop = std::chrono::high_resolution_clock::now();
  std::cerr << "Duration of packing (j = " << j << ", packarr=" << packarr[0] << packarr[999] << "): " << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << "ms" << std::endl;


  /*****************************
   *  Packed k-mer generation  *
   *****************************/

  // generate Kmers
  typedef bliss::Kmer<35, 2, uint32_t> Kmer;
  std::cout << "size of kmer: " << sizeof(Kmer) << std::endl;
  typedef bliss::PackedKmerGenerationIterator< packit_t, Kmer > kmer_gen_it_t;

  kmer_gen_it_t kmerGenIt(packIt);
  kmer_gen_it_t kmerGenEnd(packIt, dna.length());
  std::size_t nKmers = dna.length() - 35 + 1;

  /* benchmark packing + kmer generation */
  std::array<Kmer, 1000> kmer_arr;
  unsigned int i = 0;

  start = std::chrono::high_resolution_clock::now();
  //std::copy(kmerGenIt, kmerGenEnd, kmers.begin());
  for (;kmerGenIt != kmerGenEnd; ++kmerGenIt)
  {
    kmer_arr[i % 1000] = *kmerGenIt;
    ++i;
  }
  stop = std::chrono::high_resolution_clock::now();

  EXPECT_EQ(i, nKmers);
  std::cerr << "kmer[0] = " << kmer_arr[0].toString() << std::endl;
  std::cerr << "Duration of packing + kmer generation (i = " << i << "): " << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << "ms" << std::endl;


  /**************************************************
   *  Direct k-mer construction (no prior packing)  *
   **************************************************/

  typedef bliss::KmerGenerationIterator<std::string::iterator, Kmer> kmer_char_gen_it_t;

  kmer_char_gen_it_t kmerGenIt2(dna.begin(), true);
  kmer_char_gen_it_t kmerGenEnd2(dna.end(), false);
  nKmers = dna.length() - 35 + 1;


  /* benchmark kmer generation */
  std::array<Kmer, 1000> kmer_arr2;
  i = 0;

  start = std::chrono::high_resolution_clock::now();
  for (;kmerGenIt2 != kmerGenEnd2; ++kmerGenIt2)
  {
    kmer_arr2[i % 1000] = *kmerGenIt2;
    ++i;
  }
  stop = std::chrono::high_resolution_clock::now();

  EXPECT_EQ(i, nKmers);
  std::cerr << "kmer[0] = " << kmer_arr2[0].toString() << std::endl;
  std::cerr << "Duration of direct kmer generation (i = " << i << "): " << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << "ms" << std::endl;


  // check that for both methods, the final 1000 kmers are identical
  for (unsigned int j = 0; j < 1000; ++j)
  {
    EXPECT_EQ(kmer_arr[j], kmer_arr2[j]);
  }
}
