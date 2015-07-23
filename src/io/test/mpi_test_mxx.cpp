
// include google test
#include <gtest/gtest.h>
#include <io/mxx_support.hpp>
#include <cstdlib>  // rand

#include <vector>
#include <random>
#include <iterator>  // for ostream iterator.
#include <algorithm>  // for transform.
#include <unordered_set>
#include <unordered_map>

#include "common/alphabets.hpp"

#include "utils/timer.hpp"

// include files to test
#include "utils/logging.h"


/*
 * test class holding some information.  Also, needed for the typed tests
 */
template <typename T>
class Mxx2MPITest : public ::testing::Test
{
  protected:
    virtual void SetUp() {};
};


// indicate this is a typed test
TYPED_TEST_CASE_P(Mxx2MPITest);

TYPED_TEST_P(Mxx2MPITest, all2all)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  std::vector<TypeParam> send_counts(p);
  std::vector<int> src(p * (p-1) / 2, rank);

  for (int i = 0; i < p; ++i) {
    send_counts[(i + rank) %p] =  i;
  }

  // create send count first
  auto recv_counts = mxx2::all2all(src, send_counts, comm);
  int k = 0;
  for (int i = 0; i < p; ++i) {
    EXPECT_EQ(recv_counts[i], (p + rank - i) % p);
    for (int j = 0; j < recv_counts[i]; ++j) {
      src[k++] = i;
    }
  }

}




TYPED_TEST_P(Mxx2MPITest, reduce)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  // create send count first
  int result = mxx2::reduce(rank, [](int const& x, int const& y) { return x + y; }, comm, root);

  if (rank == root)
    EXPECT_EQ(result, p * (p-1) / 2);
}


TYPED_TEST_P(Mxx2MPITest, reducen)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  std::vector<int> src(p);
  for (int i = 0; i < p; ++i) {
    src[i] = i;
  }

  // create send count first
  auto result = mxx2::reduce(src, [](int const& x, int const& y) { return x + y; }, comm, root);

  if (rank == root) {
    for (int i = 0; i < p; ++i) {
      EXPECT_EQ(result[i], p * i);
    }
  }
}


TYPED_TEST_P(Mxx2MPITest, reduce_loc)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  std::vector<int> src(p);
  for (int i = 0; i < p; ++i) {
    src[(i + rank) % p] = i;
  }

  // create send count first
  auto data = mxx2::reduce_loc(src, MPI_MINLOC, comm, root);
  if (rank == root) {

     for (int i = 0; i < p; ++i) {
       EXPECT_EQ(data.first[i], 0);
       EXPECT_EQ(data.second[i], i);
     }
  }

}





TYPED_TEST_P(Mxx2MPITest, gathern)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;
  int n = 5;
  std::vector<int> src(n, rank);

  // create send count first
  std::vector<int> data = mxx2::gathern(src, comm, root);
  if (rank == root) {

     for (int i = 0; i < p; ++i) {
       for (int j = 0; j < n; ++j)
         EXPECT_EQ(data[i * n + j], i);
     }
  }

}



TYPED_TEST_P(Mxx2MPITest, scattern)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;
  int n = 5;

  // create send count first
  std::vector<int> src;
  if (rank == root) {
    // source vector
    src.resize(p * n);
    auto begin = src.begin();
    for (int i = 0; i < p; ++i) {
       std::fill(begin, begin + n, i);
       begin += n;
    }
  }
  auto data =  mxx2::scattern(src, comm, root);

  for (int i = 0; i < n; ++i)
    EXPECT_EQ(data[i], rank);

}


TYPED_TEST_P(Mxx2MPITest, gather)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  // create send count first
  std::vector<int> data = mxx2::gather(rank, comm, root);
  if (rank == root) {

     for (int i = 0; i < p; ++i) {
       EXPECT_EQ(data[i], i);
     }
  }

}



TYPED_TEST_P(Mxx2MPITest, scatter)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  // create send count first
  std::vector<int> data;
  if (rank == root) {
    // source vector
    data.resize(p);
    for (int i = 0; i < p; ++i) {
       data[i] = i;
    }
  }
  int result = mxx2::scatter(data, comm, root);

  EXPECT_EQ(result, rank);

}



TYPED_TEST_P(Mxx2MPITest, scatterv)
{

  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  // create send count first
  std::vector<int> data;
  std::vector<TypeParam> send_counts;
  if (rank == root) {
    send_counts.resize(p, 0);
    std::srand(rank);
    for (int i = 0; i < p; ++i) {
      send_counts[i] = std::rand() % 1000;
    }
    send_counts[(std::rand() %p)] = 0;

    // total:
    size_t total = std::accumulate(send_counts.begin(), send_counts.end(), 0UL, std::plus<size_t>());

    // source vector
    data.resize(total);
    auto begin = data.begin();
    for (int i = 0; i < p; ++i) {
       std::fill(begin, begin + send_counts[i], i);
       begin += send_counts[i];
    }
  }
  std::vector<int> results = mxx2::scatterv(data.begin(), data.end(), send_counts, comm, root);

  mxx::datatype<TypeParam> count_dt;
  TypeParam count;
  MPI_Scatter(&(send_counts[0]), 1, count_dt.type(), &count, 1, count_dt.type(), root, comm);

  // count should be the same
  EXPECT_EQ(count, results.size());

  // compare to expected values
  TypeParam incorrect = results.size();
  for (int i = 0; i < results.size(); ++i) {
    EXPECT_EQ(results[i], rank);
    if (results[i] == rank) --incorrect;
  }

//  printf("size = %lu, incorrect = %lu\n", results.size(), incorrect);

  EXPECT_EQ(incorrect, 0);

  TypeParam all_incorrect = 0;
  MPI_Reduce(&incorrect, &all_incorrect, 1, count_dt.type(), MPI_SUM, root, comm);

  if (rank == root)
    EXPECT_EQ(all_incorrect, 0);
}





// now register the test cases
REGISTER_TYPED_TEST_CASE_P(Mxx2MPITest, all2all, reduce, reducen, reduce_loc, gathern, scattern, gather, scatter, scatterv);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<int, size_t> Mxx2MPITestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, Mxx2MPITest, Mxx2MPITestTypes);

//
///* below is not yet working.
//*/
//#include "common/kmer.hpp"
//#include "containers/distributed_map_base.hpp"
//#include "wip/kmer_hash.hpp"
//
//class Mxx2KmerMPITest : public ::testing::TestWithParam<::std::pair<int, int> >
//{
//  protected:
//    typedef bliss::common::Kmer<21, bliss::common::DNA, uint64_t> Key;
//    typedef bliss::kmer::transform::lex_less<Key> KeyTransform;
//    typedef bliss::kmer::hash::farm<Key, false> Hash;
//
//    struct TransformedHash {
//        Hash h;
//        KeyTransform trans;
//
//        inline uint64_t operator()(Key const& k) const {
//          return h(trans(k));
//        }
//        template<typename V>
//        inline uint64_t operator()(::std::pair<Key, V> const& x) const {
//          return this->operator()(x.first);
//        }
//        template<typename V>
//        inline uint64_t operator()(::std::pair<const Key, V> const& x) const {
//          return this->operator()(x.first);
//        }
//    };
//    template <typename Comparator>
//    struct TransformedComp {
//        Comparator comp;
//        KeyTransform trans;
//
//        inline bool operator()(Key const & x, Key const & y) const {
//          return comp(trans(x), trans(y));
//        }
//        template<typename V>
//        inline bool operator()(::std::pair<Key, V> const & x, Key const & y) const {
//          return this->operator()(x.first, y);
//        }
//        template<typename V>
//        inline bool operator()(::std::pair<const Key, V> const & x, Key const & y) const {
//          return this->operator()(x.first, y);
//        }
//        template<typename V>
//        inline bool operator()(Key const & x, ::std::pair<Key, V> const & y) const {
//          return this->operator()(x, y.first);
//        }
//        template<typename V>
//        inline bool operator()(Key const & x, ::std::pair<const Key, V> const & y) const {
//          return this->operator()(x, y.first);
//        }
//        template<typename V>
//        inline bool operator()(::std::pair<Key, V> const & x, ::std::pair<Key, V> const & y) const {
//          return this->operator()(x.first, y.first);
//        }
//        template<typename V>
//        inline bool operator()(::std::pair<const Key, V> const & x, ::std::pair<const Key, V> const & y) const {
//          return this->operator()(x.first, y.first);
//        }
//
//    };
//
//    typedef TransformedComp<::std::less<Key> > TransformedLess;
//    typedef TransformedComp<::std::equal_to<Key> > TransformedEqual;
//
//    typedef std::unordered_set<Key, TransformedHash, TransformedEqual> Set;
//    typedef std::unordered_map<Key, int, TransformedHash, TransformedEqual> Map;
//
//
//    ::std::vector<Key> input;
//    ::std::vector<Key> gold;
//    ::std::vector<Key> test;
//
//    typedef ::std::pair<Key, int> P;
//
//    ::std::vector<P > reduc_input;
//    ::std::vector<P > reduc_gold;
//    ::std::vector<P > reduc_test;
//
//    int p;
//    TransformedHash hash;
//
//    virtual void SetUp()
//    { // generate some inputs
//      p = GetParam().second;
//
//
//      std::default_random_engine generator;
//      std::uniform_int_distribution<uint64_t> distribution(0, p-1);
//      uint64_t v;
//      uint64_t *vp = &v;
//      for (int i = 0; i < GetParam().first; ++i) {
//        v = distribution(generator);
//        input.emplace_back(Key(vp));
//        reduc_input.emplace_back(Key(vp), 1);
//      }
//    }
//
//};
//
//
//TEST_P(Mxx2KmerMPITest, bucketing)
//{
//
//  TIMER_INIT(bucket);
//
//  auto pp = this->p;
//  auto h = this->hash;
//
//  typedef typename Mxx2KmerMPITest::Key TypeParam;
//
//
//  TIMER_START(bucket);
//  std::vector<int> gold_counts = ::mxx2::bucketing_copy<int>(this->input, this->gold, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(bucket, "gold", this->input.size());
//
//  this->test = this->input;
//
//  TIMER_START(bucket);
//  std::vector<int> test_counts = ::mxx2::bucketing<int>(this->test, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(bucket, "test", this->input.size());
//
//  printf("bucket: count %d, buckets %d\n", this->GetParam().first, this->GetParam().second);
//  TIMER_REPORT(bucket, 0);
//
//  bool same_counts = ::std::equal(gold_counts.begin(), gold_counts.end(), test_counts.begin());
//  if (!same_counts)
//  {
//    ::std::ostream_iterator<int> oit(::std::cout, ", ");
//    ::std::cout << "gold counts: " << ::std::endl;
//    ::std::copy(gold_counts.begin(), gold_counts.end(), oit);
//    ::std::cout << std::endl;
//    ::std::cout << "test counts: " << ::std::endl;
//    ::std::copy(test_counts.begin(), test_counts.end(), oit);
//    ::std::cout << std::endl;
//  }
//  EXPECT_TRUE(same_counts);
//
//  bool same_contents = true;
//
//
//  auto start = this->test.begin();
//  auto start_gold = this->gold.begin();
//  {
//  for (size_t i = 0; i < test_counts.size(); ++i) {
//    ::std::sort(start, start + test_counts[i], Mxx2KmerMPITest::TransformedLess());
//    ::std::sort(start_gold, start_gold + gold_counts[i], Mxx2KmerMPITest::TransformedLess());
//    same_contents = ::std::equal(start_gold, start_gold + gold_counts[i], start, Mxx2KmerMPITest::TransformedEqual());
//
//    if (!same_contents)
//    {
//      ::std::cout << i << " test: ";
//      for (int j = 0; j < test_counts[i]; ++start, ++j) {
//        ::std::cout << start->getData()[0] << "(" << (h(*start) %pp) << "),";
//      }
//      ::std::cout << ::std::endl;
//      ::std::cout << i << " gold: ";
//      for (int j = 0; j < gold_counts[i]; ++start_gold, ++j) {
//        ::std::cout << start_gold->getData()[0] << "(" << (h(*start_gold) %pp) << "),";
//      }
//      ::std::cout << ::std::endl;
//    } else {
//      start += test_counts[i];
//      start_gold += gold_counts[i];
//    }
//
//    EXPECT_TRUE(same_contents);
//
//  }
//  }
//}
//
//
//TEST_P(Mxx2KmerMPITest, unique)
//{
//
//  TIMER_INIT(unique);
//
//  auto pp = this->p;
//  auto h = this->hash;
//  typedef typename Mxx2KmerMPITest::Key TypeParam;
//
//
//  this->gold = this->input;
//
//  TIMER_START(unique);
//  ::std::sort(this->gold.begin(), this->gold.end(), Mxx2KmerMPITest::TransformedLess());
//  TIMER_END(unique, "gold_sort", this->gold.size());
//  TIMER_START(unique);
//  auto end = ::std::unique(this->gold.begin(), this->gold.end(), Mxx2KmerMPITest::TransformedEqual());
//  TIMER_END(unique, "gold_unique", this->gold.size());
//  TIMER_START(unique);
//  this->gold.erase(end, this->gold.end());
//  TIMER_END(unique, "gold_erase", this->gold.size());
//  TIMER_START(unique);
//  ::std::vector<int> gold_counts = ::mxx2::bucketing<int>(this->gold, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(unique, "gold_bucket", this->gold.size());
//
//
//  {
//  auto global_reduc = this->input;
//
//  TIMER_START(unique);
//  typename Mxx2KmerMPITest::Key key;
//  typename Mxx2KmerMPITest::Set set(global_reduc.begin(), global_reduc.end(), global_reduc.size());
//  auto newend = ::std::copy(set.begin(), set.end(), global_reduc.begin());
//  global_reduc.erase(newend, global_reduc.end());
//  TIMER_END(unique, "global_unique", global_reduc.size());
//
//  TIMER_START(unique);
//  ::std::vector<int> global_counts = ::mxx2::bucketing<int>(global_reduc, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(unique, "global_bucket", global_reduc.size());
//  }
//
//
//
//  this->test = this->input;
//
//  TIMER_START(unique);
//  std::vector<int> test_counts = ::mxx2::bucketing<int>(this->test, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(unique, "test_bucket", this->test.size());
//  TIMER_START(unique);
//  ::mxx2::retain_unique<typename Mxx2KmerMPITest::Set, typename Mxx2KmerMPITest::TransformedEqual >(this->test, test_counts);
//  TIMER_END(unique, "test_unique", this->test.size());
//
//  printf("unique: count %d, buckets %d\n", this->GetParam().first, this->GetParam().second);
//
//  TIMER_REPORT(unique, 0);
//
//  bool same_counts = ::std::equal(gold_counts.begin(), gold_counts.end(), test_counts.begin());
//  if (!same_counts)
//  {
//    ::std::ostream_iterator<int> oit(::std::cout, ", ");
//    ::std::cout << "gold counts: " << ::std::endl;
//    ::std::copy(gold_counts.begin(), gold_counts.end(), oit);
//    ::std::cout << ::std::endl;
//    ::std::cout << "test counts: " << ::std::endl;
//    ::std::copy(test_counts.begin(), test_counts.end(), oit);
//    ::std::cout << ::std::endl;
//  }
//  EXPECT_TRUE(same_counts);
//
//  bool same_contents = true;
//
//  auto start = this->test.begin();
//  auto start_gold = this->gold.begin();
//  for (size_t i = 0; i < test_counts.size(); ++i) {
//    ::std::sort(start, start + test_counts[i], Mxx2KmerMPITest::TransformedLess());
//    ::std::sort(start_gold, start_gold + gold_counts[i], Mxx2KmerMPITest::TransformedLess());
//    same_contents = ::std::equal(start_gold, start_gold + gold_counts[i], start, Mxx2KmerMPITest::TransformedEqual());
//
//    if (!same_contents)
//    {
//      ::std::cout << i << " test: ";
//      for (int j = 0; j < test_counts[i]; ++start, ++j) {
//        ::std::cout << start->getData()[0] << "(" << (h(*start) %pp) << "),";
//      }
//      ::std::cout << ::std::endl;
//      ::std::cout << i << " gold: ";
//      for (int j = 0; j < gold_counts[i]; ++start_gold, ++j) {
//        ::std::cout << start_gold->getData()[0] << "(" << (h(*start_gold) %pp) << "),";
//      }
//      ::std::cout << ::std::endl;
//    } else {
//      start += test_counts[i];
//      start_gold += gold_counts[i];
//    }
//    EXPECT_TRUE(same_contents);
//
//  }
//
//}
//
//
//TEST_P(Mxx2KmerMPITest, reduce)
//{
//
//  TIMER_INIT(reduc);
//
//  auto pp = this->p;
//  auto h = this->hash;
//
//  typedef typename Mxx2KmerMPITest::P TypeParam;
//
//
//  this->reduc_gold = this->reduc_input;
//
//  TIMER_START(reduc);
//  ::std::sort(this->reduc_gold.begin(), this->reduc_gold.end(), Mxx2KmerMPITest::TransformedLess());
//  TIMER_END(reduc, "gold_sort", this->reduc_gold.size());
//
//  if (this->reduc_gold.size() > 0) {
//    auto curr = this->reduc_gold.begin();
//    auto it = this->reduc_gold.begin(); ++it;
//    TIMER_START(reduc);
//    for (auto max = this->reduc_gold.end(); it != max; ++it) {
//      if (it->first == curr->first) {
//        curr->second += it->second;
//      } else {
//        ++curr;
//        *curr = *it;
//      }
//    }
//    ++curr;  // very end, move to next entry.
//    TIMER_END(reduc, "gold_reduce", this->reduc_gold.size());
//
//    TIMER_START(reduc);
//    this->reduc_gold.erase(curr, this->reduc_gold.end());
//    TIMER_END(reduc, "gold_erase", this->reduc_gold.size());
//  }
//
//  TIMER_START(reduc);
//  ::std::vector<int> gold_counts = ::mxx2::bucketing<int>(this->reduc_gold, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(reduc, "gold_bucket", this->reduc_gold.size());
//
//
//  {
//  auto global_reduc = this->reduc_input;
//
//  TIMER_START(reduc);
//  typename Mxx2KmerMPITest::Key key;
//  int val;
//  typename Mxx2KmerMPITest::Map map(global_reduc.size());
//  for (auto it = global_reduc.begin(); it != global_reduc.end(); ++it) {
//    key = it->first;
//    val = it->second;
//    if (map.count(key) == 0) map[key] = val;  // don't rely on initialization to set T to 0.
//    else map[key] += val;
//  }
//  auto newend = ::std::copy(map.begin(), map.end(), global_reduc.begin());
//  global_reduc.erase(newend, global_reduc.end());
//  TIMER_END(reduc, "global_reduce", global_reduc.size());
//
//  TIMER_START(reduc);
//  ::std::vector<int> global_counts = ::mxx2::bucketing<int>(global_reduc, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(reduc, "global_bucket", global_reduc.size());
//  }
//
//
//  this->reduc_test = this->reduc_input;
//
//  TIMER_START(reduc);
//  std::vector<int> test_counts = ::mxx2::bucketing<int>(this->reduc_test, [&](TypeParam const & x) { return h(x) % pp; }, pp);
//  TIMER_END(reduc, "test_bucket", this->reduc_test.size());
//  TIMER_START(reduc);
//  ::mxx2::bucket_reduce<typename Mxx2KmerMPITest::Map, ::std::plus<int> >(this->reduc_test, test_counts);
//  TIMER_END(reduc, "test_reduce", this->reduc_test.size());
//
//  printf("unique: count %d, buckets %d\n", this->GetParam().first, this->GetParam().second);
//
//  TIMER_REPORT(reduc, 0);
//
//  bool same_counts = ::std::equal(gold_counts.begin(), gold_counts.end(), test_counts.begin());
//  if (!same_counts)
//  {
//    ::std::ostream_iterator<int> oit(::std::cout, ", ");
//    ::std::cout << "gold counts: " << ::std::endl;
//    ::std::copy(gold_counts.begin(), gold_counts.end(), oit);
//    ::std::cout << ::std::endl;
//    ::std::cout << "test counts: " << ::std::endl;
//    ::std::copy(test_counts.begin(), test_counts.end(), oit);
//    ::std::cout << ::std::endl;
//  }
//  EXPECT_TRUE(same_counts);
//
//  bool same_contents = true;
//
//  auto start = this->reduc_test.begin();
//  auto start_gold = this->reduc_gold.begin();
//  for (size_t i = 0; i < test_counts.size(); ++i) {
//    ::std::sort(start, start + test_counts[i], Mxx2KmerMPITest::TransformedLess());
//    ::std::sort(start_gold, start_gold + gold_counts[i], Mxx2KmerMPITest::TransformedLess());
//    same_contents = ::std::equal(start_gold, start_gold + gold_counts[i], start, Mxx2KmerMPITest::TransformedEqual());
//
//    if (!same_contents) {
//      ::std::cout << "test: ";
//      for (int j = 0; j < test_counts[i]; ++start, ++j) {
//        ::std::cout << "[" << start->first << ":" << start->second << "],";
//      }
//      ::std::cout << ::std::endl;
//      ::std::cout << "gold: ";
//      for (int j = 0; j < test_counts[i]; ++start, ++j) {
//        ::std::cout << "[" << start_gold->first << ":" << start_gold->second << "],";
//      }
//      ::std::cout << ::std::endl;
//    }
//    start += test_counts[i];
//    start_gold += gold_counts[i];
//    EXPECT_TRUE(same_contents);
//
//  }
//
//}
//
//
//
//
//// now register the test cases
////REGISTER_TEST_CASE_P(Mxx2KmerMPITest, bucketing, unique);
//
//
////////////////////// RUN the tests with different types.
//
////INSTANTIATE_TEST_CASE_P(Bliss, Mxx2KmerMPITest, ::testing::Values(::std::make_pair(100, 16),
////                                                           ::std::make_pair(100, 32)
////                                                           ));
//INSTANTIATE_TEST_CASE_P(Bliss, Mxx2KmerMPITest, ::testing::Values(::std::make_pair(100000, 16),
//                                                               ::std::make_pair(100000, 32),
//                                                               ::std::make_pair(100000, 64),
//                                                               ::std::make_pair(100000, 128),
//                                                               ::std::make_pair(100000, 256),
//                                                               ::std::make_pair(100000, 512),
//                                                               ::std::make_pair(100000, 1024),
//                                                               ::std::make_pair(100000, 2048),
//                                                               ::std::make_pair(200000, 512),
//                                                               ::std::make_pair(400000, 512),
//                                                               ::std::make_pair(800000, 512),
//                                                               ::std::make_pair(1600000, 512),
//                                                               ::std::make_pair(3200000, 512),
//                                                               ::std::make_pair(6400000, 512)
//                                                           ));
//
//
//
int main(int argc, char* argv[]) {
    int result = 0;

    ::testing::InitGoogleTest(&argc, argv);


#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif
    result = RUN_ALL_TESTS();
#if defined(USE_MPI)
    MPI_Finalize();
#endif
    return result;
}
