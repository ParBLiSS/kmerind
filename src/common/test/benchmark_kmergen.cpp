#include <gtest/gtest.h>

#include <chrono>

#include <common/alphabets.hpp>
#include <common/AlphabetTraits.hpp>
#include <common/PackingIterator.hpp>
#include <utils/generator.hpp>




TEST(Benchmark_KmerGeneration, BenchmarkKmer1)
{
  // generate a random piece of DNA
  std::string dna = bliss::utils::random_dna(1000000);

  // first step: translate (in place)
  // TODO: do this as transformation iterator
  bliss::AlphabetTraits<DNA>::translateFromAscii(dna.begin(), dna.end(), dna.begin());

  // define a packing iterator to wrap around the string's iterators
  typedef bliss::PackingIterator<std::string::iterator, bliss::AlphabetTraits<DNA>::getBitsPerChar()> packit_t;
  packit_t packIt(dna.begin(), dna.end());
  packit_t packItEnd(dna.end());

  // generate Kmers
  typedef bliss::Kmer<21, 2, uint8_t> Kmer;
  typedef bliss::KmerGenerationIterator< packit_t, Kmer > kmer_gen_it_t;

  kmer_gen_it_t kmerGenIt(packIt);
  kmer_gen_it_t kmerGenEnd(packIt, dna.length());
  std::size_t nKmers = dna.length() - 21 + 1;

  // allocate vector
  std::vector<Kmer> kmers(nKmers);

  // benchmark!
  auto start = std::chrono::high_resolution_clock::now();
  std::copy(kmerGenIt, kmerGenEnd, kmers.begin());
  auto stop = std::chrono::high_resolution_clock::now();

  std::cerr << "Duration of kmer generation: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl;

}