
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
#include "bliss-config.hpp"

#ifdef USE_MPI

/*
 * test class holding some information.  Also, needed for the typed tests
 */
template <typename T>
class Mxx2MPITest : public ::testing::Test
{
  public:
    ~Mxx2MPITest() {};

  protected:
    virtual void SetUp() {};

    template <typename TT>
    void printvec(::std::vector<TT> x, std::string const & name) {

      std::cout << name.c_str() << ": ";
      for (auto e : x) {
        std::cout << e << ", ";
      }
      std::cout << std::endl;
    }
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


  MPI_Barrier(MPI_COMM_WORLD);

}




TYPED_TEST_P(Mxx2MPITest, reduce)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  int root = p-1;

  // create send count first
  int result = mxx2::reduce::reduce(rank, [](int const& x, int const& y) { return x + y; }, comm, root);

  if (rank == root)
    EXPECT_EQ(result, p * (p-1) / 2);

  MPI_Barrier(MPI_COMM_WORLD);

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
  auto result = mxx2::reduce::reduce(src, [](int const& x, int const& y) { return x + y; }, comm, root);

  if (rank == root) {
    for (int i = 0; i < p; ++i) {
      EXPECT_EQ(result[i], p * i);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

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
  auto data = mxx2::reduce::reduce_loc(src, MPI_MINLOC, comm, root);
  if (rank == root) {

     for (int i = 0; i < p; ++i) {
       EXPECT_EQ(data.first[i], 0);
       EXPECT_EQ(data.second[i], i);
     }
  }
  MPI_Barrier(MPI_COMM_WORLD);

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
  MPI_Barrier(MPI_COMM_WORLD);

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

  MPI_Barrier(MPI_COMM_WORLD);


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
  MPI_Barrier(MPI_COMM_WORLD);


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

  MPI_Barrier(MPI_COMM_WORLD);

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

  MPI_Barrier(MPI_COMM_WORLD);

}

TYPED_TEST_P(Mxx2MPITest, unique)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> values(p);

  for (int i = 0; i < rank; ++i) {
    values[i] = (i / 3) % 3;
  }
  for (int i = rank; i < p; ++i) {
    values[i] = (i / 3 + 1) % 3;
  }
  std::vector<TypeParam> allvals = mxx2::gathern(values, MPI_COMM_WORLD, 0);

  auto end = mxx2::unique_contiguous(values, MPI_COMM_WORLD, [](TypeParam & x, TypeParam & y){ return x == y; });
  values.erase(end, values.end());

  std::vector<TypeParam> result = mxx::gather_vectors(values, MPI_COMM_WORLD);

  if (rank == 0) {
    auto gold = allvals;
    auto end2 = std::unique(gold.begin(), gold.end(), [](TypeParam & x, TypeParam & y){ return x == y; });
    gold.erase(end2, gold.end());


    bool same = std::equal(result.begin(), result.end(), gold.begin());

    if (!same) {

      this->template printvec(allvals, "input");
      this->template printvec(result, "unique");
      this->template printvec(gold, "unique gold");
    }
    EXPECT_TRUE(same);
  }
}




TYPED_TEST_P(Mxx2MPITest, scan)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  TypeParam input = 1;



  // forward
  TypeParam scanned = mxx2::scan::scan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  //printf("R %d scan, %lu, ", rank, scanned);
  EXPECT_EQ(scanned, rank + 1);

  // reverse
  scanned = mxx2::scan::rscan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  //printf("R %d scan_rev, %lu, ", rank, scanned);
  EXPECT_EQ(scanned, p - rank);

  // exclusive
  scanned = mxx2::scan::exscan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank > 0) {
    //printf("R %d exscan, %lu, ", rank, scanned);
    EXPECT_EQ(scanned, rank);
  }


  // exclusive reverse
  scanned = mxx2::scan::rexscan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank < (p-1)) {
    //printf("R %d exscan_rev, %lu, ", rank, scanned);
    EXPECT_EQ(scanned, p - rank - 1);
  }

}

TYPED_TEST_P(Mxx2MPITest, scan_n)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < p; ++i) {
    inputs[i] = i;
  }



  // forward
  std::vector<TypeParam> scanned = mxx2::scan::scan(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], i * (rank + 1));
  }
  // reverse
  scanned = mxx2::scan::rscan(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan_rev, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], i * (p - rank));
  }
  // exclusive
  scanned = mxx2::scan::exscan(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank > 0)
    for (int i = 0; i < p; ++i) {
      //printf("R %d exscan, %d = %lu \n", rank, i, scanned[i]);
      EXPECT_EQ(scanned[i], i * rank);
    }

  // exclusive reverse
  scanned = mxx2::scan::rexscan(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank < (p-1))
    for (int i = 0; i < p; ++i) {
      //printf("R %d exscan_rev, %d = %lu \n", rank, i, scanned[i]);
      EXPECT_EQ(scanned[i], i * (p - rank - 1));
    }
}

TYPED_TEST_P(Mxx2MPITest, scan_vec)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < p; ++i) {
    inputs[i] = 1;
  }


  // forward
  std::vector<TypeParam> scanned = mxx2::scan::scan_contiguous(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * rank + (i + 1));
  }
  // reverse
  scanned = mxx2::scan::rscan_contiguous(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan_rev, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * (p - rank - 1) + (p - i));
  }
  // exclusive
  scanned = mxx2::scan::exscan_contiguous(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = (rank > 0) ? 0 : 1; i < p; ++i) {
    //printf("R %d exscan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * rank + i);
  }

  // exclusive reverse
  scanned = mxx2::scan::rexscan_contiguous(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < (rank < (p-1) ? p : p-1); ++i) {
    //printf("R %d exscan_rev, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * (p - rank - 1) + (p - i - 1));
  }
}


TYPED_TEST_P(Mxx2MPITest, seg_convert)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  TypeParam seg = (rank / 2) % 2 ;  // non unique

  // convert to start
  uint8_t start = mxx2::segment::is_start(seg, MPI_COMM_WORLD);
  //printf("R %d seg start, %d\n ", rank, start);
  EXPECT_EQ((rank % 2 == 0 ? 1 : 0), start);

  // convert to end
  uint8_t end = mxx2::segment::is_end(seg, MPI_COMM_WORLD);
  //printf("R %d seg end, %d\n ", rank, end);
  EXPECT_EQ(rank % 2 == 1 ? 1 : (rank == p-1) ? 1 : 0, end);


  // convert to unique from start
  TypeParam useg = mxx2::segment::to_unique_segment_id_from_start<TypeParam>(start, 0, MPI_COMM_WORLD);
  //printf("R %d uniq_seg, %d \n ", rank, useg);
  EXPECT_EQ((rank / 2 + 1), useg);


  // convert to unique from end
  useg = mxx2::segment::to_unique_segment_id_from_end<TypeParam>(end, 0, MPI_COMM_WORLD);
  //printf("R %d uniq_seg_rev, %d \n ", rank, useg);
  EXPECT_EQ((p + 1)/2 - rank/2, useg);


  // convert from non-unique to unique
  useg = mxx2::segment::to_unique_segment_id(seg, MPI_COMM_WORLD);
  //printf("R %d uniq_seg 2, %d \n ", rank, useg);
  EXPECT_EQ((rank / 2 + 1), useg);


}

TYPED_TEST_P(Mxx2MPITest, seg_convert_n)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < p; ++i) {
    inputs[i] = ((rank + i) / 2) % 2;
  }
//  for (int i = 0; i < p; ++i) {
//    printf("R %d input, %d = %lu \n", rank, i, inputs[i]);
//  }

  // convert to start
  auto starts = mxx2::segment::is_start(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
        //printf("R %d seg start, %d = %lu \n", rank, i, starts[i]);
        EXPECT_EQ(((rank + i) % 2 == 0 ? 1 : (rank == 0) ? 1 : 0), starts[i]);
  }
  // convert to end
  auto ends = mxx2::segment::is_end(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
          //printf("R %d seg end, %d = %lu \n", rank, i, ends[i]);
          EXPECT_EQ((rank+i) % 2 == 1 ? 1 : (rank == p-1) ? 1 : 0, ends[i]);
  }

  // convert to unique from start
  auto useg = mxx2::segment::to_unique_segment_id_from_start<TypeParam>(starts, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
        //printf("R %d useg, %d = %lu \n", rank, i, useg[i]);
        EXPECT_EQ(((rank + i % 2) / 2 + 1), useg[i]);
  }

  // convert to unique from end
  useg = mxx2::segment::to_unique_segment_id_from_end<TypeParam>(ends, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
      //  printf("R %d useg rev, %d = %lu \n", rank, i, useg[i]);
        EXPECT_EQ((p + 1 + i % 2) / 2 - (rank + i % 2)/2, useg[i]);
  }

  // convert from non-unique to unique
  useg = mxx2::segment::to_unique_segment_id(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d uniq_seg 2, %d \n ", rank, useg[i]);
    EXPECT_EQ(((rank + i % 2) / 2 + 1), useg[i]);
  }

}


TYPED_TEST_P(Mxx2MPITest, seg_convert_vec)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < rank; ++i) {
    inputs[i] = rank % 2;
  }
  for (int i = rank; i < p; ++i) {
    inputs[i] = (rank + 1) % 2;
  }
//  for (int i = 0; i < p; ++i) {
//    printf("R %d input, %d = %lu \n", rank, i, inputs[i]);
//  }


  // convert to start
  auto starts = mxx2::segment::is_start_contiguous(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d seg start, %d = %lu \n", rank, i, starts[i]);
    EXPECT_EQ(starts[i], (i == rank) ? 1 : 0);
  }

  // convert to end
  auto ends = mxx2::segment::is_end_contiguous(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//      printf("R %d seg end, %d = %lu \n", rank, i, ends[i]);
      EXPECT_EQ(ends[i], (i == (rank-1)) ? 1 : ((rank == p - 1) && (i == p-1 )) ? 1 : 0);
  }


  // convert to unique from start
  auto useg = mxx2::segment::to_unique_segment_id_from_start_contiguous<TypeParam>(starts, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d useg, %d = %lu \n", rank, i, useg[i]);
    EXPECT_EQ(useg[i], (i < rank) ? rank : rank + 1);
  }


  // convert to unique from end
  useg = mxx2::segment::to_unique_segment_id_from_end_contiguous<TypeParam>(ends, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d useg rev, %d = %lu \n", rank, i, useg[i]);
    EXPECT_EQ(useg[i], (i < rank) ? p - rank + 1 : p - rank);
  }


  // convert from non-unique to unique
  useg = mxx2::segment::to_unique_segment_id_contiguous(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d useg 2, %d = %lu \n", rank, i, useg[i]);
    EXPECT_EQ(useg[i], (i < rank) ? rank : rank + 1);
  }


}



TYPED_TEST_P(Mxx2MPITest, reduce2)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  TypeParam input = rank;



  // forward
  TypeParam scanned = mxx2::reduce::reduce(input, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0) EXPECT_EQ(scanned, p-1);

  // reverse
  scanned = mxx2::reduce::allreduce(input, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  EXPECT_EQ(scanned, p-1);

}

TYPED_TEST_P(Mxx2MPITest, reduce_n)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < p; ++i) {
    inputs[i] = (rank + i) % p;
  }



  // forward
  std::vector<TypeParam> scanned = mxx2::reduce::reduce(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0)
    for (int i = 0; i < p; ++i)
      EXPECT_EQ(scanned[i], p-1);

  // reverse
  scanned = mxx2::reduce::allreduce(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i)
    EXPECT_EQ(scanned[i], p-1);

}

TYPED_TEST_P(Mxx2MPITest, reduce_vec)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  std::vector<TypeParam> inputs(p);
  for (int i = 0; i < p; ++i) {
    inputs[i] = rank * p + i;
  }



  // forward
  TypeParam scanned = mxx2::reduce::reduce_contiguous(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0)
    EXPECT_EQ(scanned, p * p - 1);

  // reverse
  scanned = mxx2::reduce::allreduce_contiguous(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  EXPECT_EQ(scanned, p * p - 1);

}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(Mxx2MPITest, all2all, reduce, reducen, reduce_loc, gathern, scattern, gather, scatter, scatterv, scan, scan_n, scan_vec, reduce2, reduce_n, reduce_vec, seg_convert, seg_convert_n, seg_convert_vec, unique);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<int, size_t> Mxx2MPITestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, Mxx2MPITest, Mxx2MPITestTypes);



/*
 * test class holding some information.  Also, needed for the typed tests
 */
class Mxx2SegmentedMPITest : public ::testing::Test
{
  public:
    ~Mxx2SegmentedMPITest() {};

  protected:

    virtual void SetUp() {};

    template <typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_reduce(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;

      // forward.
      result[0] = x[0];
      for (int i = 1; i < x.size(); ++i) {
        if (seg[i-1] == seg[i]) {
          result[i] = f(result[i-1], x[i]);
        }
        else result[i] = x[i];
      }

      // then reverse
      for (int i = result.size() - 2; i >= 0; --i) {
        if (seg[i+1] == seg[i]) {
          result[i] = result[i+1];
        }
      }

      return result;
    }

    template <typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_scan(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;

      result[0] = x[0];
      for (int i = 1; i < x.size(); ++i) {
        if (seg[i-1] == seg[i]) {
          result[i] = f(result[i-1], x[i]);
        }
        else result[i] = x[i];
      }
      return result;
    }


    template <typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_exscan(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 2) return result;

      result[0] = TT();
      TT temp = x[0];
      for (int i = 1; i < x.size(); ++i) {
        if (seg[i-1] == seg[i]) {
          result[i] = temp;
          temp = f(temp, x[i]);
        }
        else {
          result[i] = TT();
          temp = x[i];
        }
      }
      return result;
    }

    template <typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_rscan(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;


      result.back() = x.back();
      for (int i = x.size() - 2; i >= 0; --i) {
        if (seg[i+1] == seg[i]) {
          result[i] = f(result[i+1], x[i]);
        }
        else result[i] = x[i];
      }
      return result;
    }

    template <typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_rexscan(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 2) return result;

      result.back() = TT();
      TT temp = x.back();
      for (int i = x.size()-2; i >= 0; --i) {
        if (seg[i+1] == seg[i]) {
          result[i] = temp;
          temp = f(temp, x[i]);
        }
        else {
          result[i] = TT();
          temp = x[i];
        }
      }
      return result;
    }

    template <bool reverse, bool exclusive, typename TT, typename ST>
    bool equal(std::vector<TT> test, std::vector<TT> gold, std::vector<ST> seg) {
      if (!exclusive) {
        return std::equal(test.begin(), test.end(), gold.begin());
      }
      // exclusive.

      auto starts = seg;

      // not in start/end form.
      if (reverse) {
        starts.back() = 1;
        std::transform(seg.rbegin(), seg.rend() - 1, seg.rbegin() + 1, starts.rbegin() + 1, [](ST &x, ST &y){ return (x == y) ? 0 : 1; });
      } else {
        starts.front() = 1;
        std::transform(seg.begin(), seg.end() - 1, seg.begin() + 1, starts.begin() + 1, [](ST &x, ST &y){ return (x == y) ? 0 : 1; });
      }

//        printvec(starts, "starts");

      bool same = true;

      for (int i = 0; i < test.size() && i < gold.size(); ++i) {
        if (starts[i] == 0) same &= (test[i] == gold[i]);
      }
      return same;


    }

    template <typename TT>
    void printvec(::std::vector<TT> x, std::string const & name) {

      printf("%s: ", name.c_str());
      for (auto e : x) {
        printf("%d, ", e);
      }
      printf("\n");

    }

};


TEST_F(Mxx2SegmentedMPITest, seg_scan)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  uint seg = rank / 3;

  auto allinputs = mxx2::gather(rank, comm, 0);

  auto allseg = mxx2::gather(seg, comm, 0);

  std::vector<int> gold;
  bool same;

  int result = mxx2::segmented_scan::scan(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  auto print = mxx2::gather(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_scan(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<false, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "scan");
      this->template printvec(gold, "scan gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::segmented_scan::exscan(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  MPI_Barrier(MPI_COMM_WORLD);
  print = mxx2::gather(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_exscan(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<false, true>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "exscan");
      this->template printvec(gold, "exscan gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::segmented_scan::rscan(rank, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  MPI_Barrier(MPI_COMM_WORLD);
  print = mxx2::gather(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_rscan(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });

    same = this->template equal<true, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rscan");
      this->template printvec(gold, "rscan gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::segmented_scan::rexscan(rank, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  MPI_Barrier(MPI_COMM_WORLD);
  print = mxx2::gather(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_rexscan(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });

    same = this->template equal<true, true>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rexscan");
      this->template printvec(gold, "rexscan gold");
    }
    EXPECT_TRUE(same);
  }
}

TEST_F(Mxx2SegmentedMPITest, seg_scan_n)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<uint> seg(3);
    seg[0] = rank / 3;
    seg[1] = (rank + 1) / 3;
    seg[2] = (rank + 2) / 3;

  std::vector<int> inputs(3);
  inputs[0] = rank;
  inputs[1] = rank;
  inputs[2] = rank;

  std::vector<int> gold;
  bool same;

  auto result = mxx2::segmented_scan::scan(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_scan(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

      same = this->template equal<false, false>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "scann");
        this->template printvec(golds, "scann gold");
      }


      EXPECT_TRUE(same);

    }

  }
  MPI_Barrier(MPI_COMM_WORLD);

  result = mxx2::segmented_scan::exscan(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_exscan(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

      same = this->template equal<false, true>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "exscann");
        this->template printvec(golds, "exscann gold");
      }


      EXPECT_TRUE(same);

    }

  }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::segmented_scan::rscan(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_rscan(vals, allseg, [](int const &x, int const &y){ return std::max(x, y); } );

      same = this->template equal<true, false>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "rscann");
        this->template printvec(golds, "rscann gold");
      }


      EXPECT_TRUE(same);

    }

  }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::segmented_scan::rexscan(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_rexscan(vals, allseg, [](int const &x, int const &y){ return std::max(x, y); } );

      same = this->template equal<true, true>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "rexscann");
        this->template printvec(golds, "rexscann gold");
      }


      EXPECT_TRUE(same);

    }

  }

  MPI_Barrier(MPI_COMM_WORLD);
}


TEST_F(Mxx2SegmentedMPITest, seg_scan_v)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector< uint> seg(p);
    for (int i = 0; i < rank; ++i) {
      seg[i] = rank - 1;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank;
    }
  std::vector<int> inputs(p);
  for (int i = 0; i < p; ++i)
    inputs[i] = i;

//  for (int i = 0; i < p; ++i) {
//    printf("R %d,%d input %lu, seg %d \n", rank, i, inputs[i], seg[i]);
//  }

  auto allinputs = mxx2::gathern(inputs, comm, 0);

  auto allseg = mxx2::gathern(seg, comm, 0);

  std::vector<int> gold;
  bool same;

  auto result = mxx2::segmented_scan::scan_contiguous(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);
  auto print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_scan(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<false, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "scanv");
      this->template printvec(gold, "scanv gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::segmented_scan::exscan_contiguous(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_exscan(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });


    same = this->template equal<false, true>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "exscanv");
      this->template printvec(gold, "exscanv gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::segmented_scan::rscan_contiguous(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_rscan(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });

    same = this->template equal<true, false>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rscanv");
      this->template printvec(gold, "rscanv gold");
    }

    EXPECT_TRUE(same);
  }

  result = mxx2::segmented_scan::rexscan_contiguous(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {

    gold = this->template seq_seg_rexscan(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });


    same = this->template equal<true, true>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rexscanv");
      this->template printvec(gold, "rexscanv gold");
    }

    EXPECT_TRUE(same);

  }
}

TEST_F(Mxx2SegmentedMPITest, seg_reduce)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  uint seg = rank / 3;

  auto allseg = mxx2::gather(seg, comm, 0);

  auto vals = mxx2::gather(rank, comm, 0);


  int result = mxx2::segmented_reduce::allreduce(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  auto outs = mxx2::gather(result, comm, 0);

  if (rank == 0) {

    auto golds = this->template seq_seg_reduce(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

    bool same = this->template equal<false, false>(outs, golds, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(vals, "input");
      this->template printvec(outs, "reduce");
      this->template printvec(golds, "reduce gold");
    }

    EXPECT_TRUE(same);

  }


  result = mxx2::segmented_reduce::allreduce(rank, seg, std::plus<int>(), comm);

  outs = mxx2::gather(result, comm, 0);

  if (rank == 0) {

    auto golds = this->template seq_seg_reduce(vals, allseg, std::plus<int>() );

    bool same = this->template equal<false, false>(outs, golds, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(vals, "input");
      this->template printvec(outs, "reduce+");
      this->template printvec(golds, "reduce+ gold");
    }


    EXPECT_TRUE(same);

  }
  MPI_Barrier(MPI_COMM_WORLD);

}


TEST_F(Mxx2SegmentedMPITest, seg_reduce_n)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<uint> seg(3);

    seg[0] = rank / 3;
    seg[1] = (rank + 1) / 3;
    seg[2] = (rank + 2) / 3;
  std::vector<int> inputs(3);
  inputs[0] = rank;
  inputs[1] = rank;
  inputs[2] = rank;


  bool same;
  auto result = mxx2::segmented_reduce::allreduce(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_reduce(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

      same = this->template equal<false, false>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "reducen");
        this->template printvec(golds, "reducen gold");
      }


      EXPECT_TRUE(same);

    }

  }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::segmented_reduce::allreduce(inputs, seg, std::plus<int>(), comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {
      auto golds = this->template seq_seg_reduce(vals, allseg, std::plus<int>() );

      same = this->template equal<false, false>(outs, golds, allseg);

      if (!same) {
        this->template printvec(allseg, "seg");
        this->template printvec(vals, "input");
        this->template printvec(outs, "reducen+");
        this->template printvec(golds, "reducen+ gold");
      }

      EXPECT_TRUE(same);

    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


TEST_F(Mxx2SegmentedMPITest, seg_reduce_v)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<uint> seg(p);
  auto allseg = seg;

    for (int i = 0; i < rank; ++i) {
      seg[i] = rank - 1;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank;
    }

  std::vector<int> inputs(p);
  for (int i = 0; i < p; ++i)
    inputs[i] = i;

//  for (int i = 0; i < p; ++i) {
//    printf("R %d,%d input %lu, seg %d \n", rank, i, inputs[i], seg[i]);
//  }

  auto allinputs = mxx2::gathern(inputs, comm, 0);

  allseg = mxx2::gathern(seg, comm, 0);

  std::vector<int> gold;
  bool same;

  auto result = mxx2::segmented_reduce::allreduce_contiguous(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);
  auto print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_reduce(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<false, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "reducev");
      this->template printvec(gold, "reducev gold");
    }
    EXPECT_TRUE(same);
  }

  result = mxx2::segmented_reduce::allreduce_contiguous(inputs, seg, std::plus<int>(), comm);

  MPI_Barrier(MPI_COMM_WORLD);
  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_reduce(allinputs, allseg, std::plus<int>());

    same = this->template equal<false, false>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "reducev+");
      this->template printvec(gold, "reducev+ gold");
    }

    EXPECT_TRUE(same);
  }

}

int main(int argc, char* argv[]) {
    int result = 0;

    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);

    result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#else

TEST(Mxx2MPITest, fail)
{
  ASSERT_TRUE(false) << "tests not compiled with USE_MPI";
}

int main(int argc, char* argv[]) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  result = RUN_ALL_TESTS();

  return result;
}
#endif
