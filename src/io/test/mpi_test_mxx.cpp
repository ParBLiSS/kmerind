
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
  int result = mxx2::reduce_op::reduce(rank, [](int const& x, int const& y) { return x + y; }, comm, root);

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
  auto result = mxx2::reduce_op::reduce_elementwise(src, [](int const& x, int const& y) { return x + y; }, comm, root);

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
  auto data = mxx2::reduce_op::reduce_loc_elementwise(src, MPI_MINLOC, comm, root);
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

TYPED_TEST_P(Mxx2MPITest, scan)
{
  int rank;
  int p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  TypeParam input = 1;



  // forward
  TypeParam scanned = mxx2::scan_op::scan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  //printf("R %d scan, %lu, ", rank, scanned);
  EXPECT_EQ(scanned, rank + 1);

  // reverse
  scanned = mxx2::scan_op::scan_reverse(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  //printf("R %d scan_rev, %lu, ", rank, scanned);
  EXPECT_EQ(scanned, p - rank);

  // exclusive
  scanned = mxx2::scan_op::exscan(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank > 0) {
    //printf("R %d exscan, %lu, ", rank, scanned);
    EXPECT_EQ(scanned, rank);
  }


  // exclusive reverse
  scanned = mxx2::scan_op::exscan_reverse(input, std::plus<TypeParam>(), MPI_COMM_WORLD);
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
  std::vector<TypeParam> scanned = mxx2::scan_op::scan_elementwise(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], i * (rank + 1));
  }
  // reverse
  scanned = mxx2::scan_op::scan_reverse_elementwise(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan_rev, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], i * (p - rank));
  }
  // exclusive
  scanned = mxx2::scan_op::exscan_elementwise(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  if (rank > 0)
    for (int i = 0; i < p; ++i) {
      //printf("R %d exscan, %d = %lu \n", rank, i, scanned[i]);
      EXPECT_EQ(scanned[i], i * rank);
    }

  // exclusive reverse
  scanned = mxx2::scan_op::exscan_reverse_elementwise(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
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
  std::vector<TypeParam> scanned = mxx2::scan_op::scanv(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * rank + (i + 1));
  }
  // reverse
  scanned = mxx2::scan_op::scanv_reverse(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
    //printf("R %d scan_rev, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * (p - rank - 1) + (p - i));
  }
  // exclusive
  scanned = mxx2::scan_op::exscanv(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
  for (int i = (rank > 0) ? 0 : 1; i < p; ++i) {
    //printf("R %d exscan, %d = %lu \n", rank, i, scanned[i]);
    EXPECT_EQ(scanned[i], p * rank + i);
  }

  // exclusive reverse
  scanned = mxx2::scan_op::exscanv_reverse(inputs, std::plus<TypeParam>(), MPI_COMM_WORLD);
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
  uint8_t start = mxx2::segment<TypeParam>::is_start(seg, MPI_COMM_WORLD);
  //printf("R %d seg start, %d\n ", rank, start);
  EXPECT_EQ((rank % 2 == 0 ? 1 : 0), start);

  // convert to end
  uint8_t end = mxx2::segment<TypeParam>::is_end(seg, MPI_COMM_WORLD);
  //printf("R %d seg end, %d\n ", rank, end);
  EXPECT_EQ(rank % 2 == 1 ? 1 : (rank == p-1) ? 1 : 0, end);


  // convert to unique from start
  TypeParam useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_start(start, 0, MPI_COMM_WORLD);
  //printf("R %d uniq_seg, %d \n ", rank, useg);
  EXPECT_EQ((rank / 2 + 1), useg);


  // convert to unique from end
  useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_end(end, 0, MPI_COMM_WORLD);
  //printf("R %d uniq_seg_rev, %d \n ", rank, useg);
  EXPECT_EQ((p + 1)/2 - rank/2, useg);


  // convert from non-unique to unique
  useg = mxx2::segment<TypeParam>::to_unique_segment_id(seg, MPI_COMM_WORLD);
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
  auto starts = mxx2::segment<TypeParam>::is_start_elementwise(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
        //printf("R %d seg start, %d = %lu \n", rank, i, starts[i]);
        EXPECT_EQ(((rank + i) % 2 == 0 ? 1 : (rank == 0) ? 1 : 0), starts[i]);
  }
  // convert to end
  auto ends = mxx2::segment<TypeParam>::is_end_elementwise(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
          //printf("R %d seg end, %d = %lu \n", rank, i, ends[i]);
          EXPECT_EQ((rank+i) % 2 == 1 ? 1 : (rank == p-1) ? 1 : 0, ends[i]);
  }

  // convert to unique from start
  auto useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_start_elementwise(starts, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
        //printf("R %d useg, %d = %lu \n", rank, i, useg[i]);
        EXPECT_EQ(((rank + i % 2) / 2 + 1), useg[i]);
  }

  // convert to unique from end
  useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_end_elementwise(ends, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
      //  printf("R %d useg rev, %d = %lu \n", rank, i, useg[i]);
        EXPECT_EQ((p + 1 + i % 2) / 2 - (rank + i % 2)/2, useg[i]);
  }

  // convert from non-unique to unique
  useg = mxx2::segment<TypeParam>::to_unique_segment_id_elementwise(inputs, MPI_COMM_WORLD);
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
  auto starts = mxx2::segment<TypeParam>::is_start_v(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d seg start, %d = %lu \n", rank, i, starts[i]);
    EXPECT_EQ(starts[i], (i == rank) ? 1 : 0);
  }

  // convert to end
  auto ends = mxx2::segment<TypeParam>::is_end_v(inputs, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//      printf("R %d seg end, %d = %lu \n", rank, i, ends[i]);
      EXPECT_EQ(ends[i], (i == (rank-1)) ? 1 : ((rank == p - 1) && (i == p-1 )) ? 1 : 0);
  }


  // convert to unique from start
  auto useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_start_v(starts, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d useg, %d = %lu \n", rank, i, useg[i]);
    EXPECT_EQ(useg[i], (i < rank) ? rank : rank + 1);
  }


  // convert to unique from end
  useg = mxx2::segment<TypeParam>::to_unique_segment_id_from_end_v(ends, 0, MPI_COMM_WORLD);
  for (int i = 0; i < p; ++i) {
//    printf("R %d useg rev, %d = %lu \n", rank, i, useg[i]);
    EXPECT_EQ(useg[i], (i < rank) ? p - rank + 1 : p - rank);
  }


  // convert from non-unique to unique
  useg = mxx2::segment<TypeParam>::to_unique_segment_id_v(inputs, MPI_COMM_WORLD);
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
  TypeParam scanned = mxx2::reduce_op::reduce(input, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0) EXPECT_EQ(scanned, p-1);

  // reverse
  scanned = mxx2::reduce_op::allreduce(input, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
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
  std::vector<TypeParam> scanned = mxx2::reduce_op::reduce_elementwise(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0)
    for (int i = 0; i < p; ++i)
      EXPECT_EQ(scanned[i], p-1);

  // reverse
  scanned = mxx2::reduce_op::allreduce_elementwise(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
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
  TypeParam scanned = mxx2::reduce_op::reducev(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  if (rank == 0)
    EXPECT_EQ(scanned, p * p - 1);

  // reverse
  scanned = mxx2::reduce_op::allreducev(inputs, [](TypeParam & x, TypeParam & y) { return x > y ? x : y; }, MPI_COMM_WORLD);
  EXPECT_EQ(scanned, p * p - 1);

}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(Mxx2MPITest, all2all, reduce, reducen, reduce_loc, gathern, scattern, gather, scatter, scatterv, scan, scan_n, scan_vec, reduce2, reduce_n, reduce_vec, seg_convert, seg_convert_n, seg_convert_vec);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<int, size_t> Mxx2MPITestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, Mxx2MPITest, Mxx2MPITestTypes);



/*
 * test class holding some information.  Also, needed for the typed tests
 */
template <typename T>
class Mxx2SegmentedMPITest : public ::testing::Test
{
  protected:
    virtual void SetUp() {};

    template <bool start, typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_reduce_v(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;

      // forward.
      result[0] = x[0];
      for (int i = 1; i < x.size(); ++i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i-1] == seg[i]))) {
          result[i] = f(result[i-1], x[i]);
        }
        else result[i] = x[i];
      }

      // then reverse
      for (int i = result.size() - 2; i >= 0; --i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i+1] == seg[i]))) {
          result[i] = result[i+1];
        }
      }

      return result;
    }

    template <bool start, typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_scan_v(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;

      result[0] = x[0];
      for (int i = 1; i < x.size(); ++i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i-1] == seg[i]))) {
          result[i] = f(result[i-1], x[i]);
        }
        else result[i] = x[i];
      }
      return result;
    }


    template <bool start, typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_exscan_v(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 2) return result;

      result[0] = TT();
      result[1] = x[0];
      for (int i = 2; i < x.size(); ++i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i-1] == seg[i]))) {
          result[i] = f(result[i-1], x[i-1]);
        }
        else {
          result[i] = TT();
          ++i;
          if (i < x.size()) result[i] = x[i-1];
        }
      }
      return result;
    }

    template <bool start, typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_scan_v_reverse(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 1) return result;


      result.back() = x.back();
      for (int i = x.size() - 2; i >= 0; --i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i+1] == seg[i]))) {
          result[i] = f(result[i+1], x[i]);
        }
        else result[i] = x[i];
      }
      return result;
    }

    template <bool start, typename TT, typename ST, typename Func>
    std::vector<TT> seq_seg_exscan_v_reverse(std::vector<TT> x, std::vector<ST> seg, Func f) {
      std::vector<TT> result(x.size());

      if (x.size() < 2) return result;

      int i = x.size()-1;
      result[i] = TT(); --i;
      result[i] = x[i+1];
      for (; i >= 0; --i) {
        if ((start && (seg[i] == 0)) || (!start && (seg[i+1] == seg[i]))) {
          result[i] = f(result[i+1], x[i+1]);
        }
        else {
          result[i] = TT();
          --i;
          if (i >= 0) result[i] = x[i+1];
        }
      }
      return result;
    }

    template <bool start, bool reverse, bool exclusive, typename TT, typename ST>
    bool equal(std::vector<TT> test, std::vector<TT> gold, std::vector<ST> seg) {
      if (!exclusive) {
        return std::equal(test.begin(), test.end(), gold.begin());
      }
      // exclusive.

      auto starts = seg;

      if (!start) {
        // not in start/end form.
        if (reverse) {
          starts.back() = 1;
          std::transform(seg.rbegin(), seg.rend() - 1, seg.rbegin() + 1, starts.rbegin() + 1, [](ST &x, ST &y){ return (x == y) ? 0 : 1; });
        } else {
          starts.front() = 1;
          std::transform(seg.begin(), seg.end() - 1, seg.begin() + 1, starts.begin() + 1, [](ST &x, ST &y){ return (x == y) ? 0 : 1; });
        }

//        printvec(starts, "starts");
      }  // already in start/end form


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


// indicate this is a typed test
TYPED_TEST_CASE_P(Mxx2SegmentedMPITest);

TYPED_TEST_P(Mxx2SegmentedMPITest, seg_scan)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  uint seg;
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) seg = rank / 3;
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) seg = (rank / 3) % 3;
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) seg = (rank % 3 == 0);

  int result = mxx2::seg_scan::scan<TypeParam>(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  if (result != rank - (rank % 3)) {
    printf("R %d seg %d scan result %d\n", rank, seg, result);
  }
  EXPECT_EQ(result, rank - (rank % 3));

  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::seg_scan::exscan<TypeParam>(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  int gold;
  if (rank % 3 > 0) {
    gold = std::max(0, rank - (rank % 3));
    if (result != gold) {
      printf("R %d seg %d scan result %d\n", rank, seg, result);
    }
    EXPECT_EQ(result, gold);
  }

  MPI_Barrier(MPI_COMM_WORLD);

// for reverse
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) seg = (rank % 3 == 2) || (rank == p-1);


  result = mxx2::seg_scan::scan_reverse<TypeParam>(rank, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  gold = std::min(p-1, rank + 2 - (rank % 3));
  if (rank < p-1) {
    if (result != gold) {
      printf("R %d seg %d scan result %d\n", rank, seg, result);
    }
    EXPECT_EQ(result, gold);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::seg_scan::exscan_reverse<TypeParam>(rank, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  if ((rank % 3 < 2) && (rank < p-1)) {
    gold = std::min(p-1, rank + 2 - (rank % 3));
    if (result != gold) {
      printf("R %d seg %d scan result %d\n", rank, seg, result);
    }
    EXPECT_EQ(result, gold);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

TYPED_TEST_P(Mxx2SegmentedMPITest, seg_scan_n)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<uint> seg(3);
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) {
    seg[0] = rank / 3;
    seg[1] = (rank + 1) / 3;
    seg[2] = (rank + 2) / 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) {
    seg[0] = (rank / 3) % 3;
    seg[1] = ((rank + 1) / 3) % 3;
    seg[2] = ((rank + 2) / 3) % 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    seg[0] = (rank % 3 == 0) ? 1 : 0;
    seg[1] = ((rank + 1) % 3 == 0) ? 1 : 0;
    seg[2] = ((rank + 2) % 3 == 0) ? 1 : 0;
  }
  std::vector<int> inputs(3);
  inputs[0] = rank;
  inputs[1] = rank;
  inputs[2] = rank;


  int gold;
  auto result = mxx2::seg_scan::scan_elementwise<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    gold = std::max(0, rank - ((rank + i) % 3));
    if (result[i] != gold) {
      printf("R %d seg %d scan result %d\n", rank, seg[i], result[i]);
    }
    EXPECT_EQ(gold, result[i]);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::seg_scan::exscan_elementwise<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    gold = std::max(0,rank - ((rank + i) % 3));
    if ((rank + i) % 3 > 0) {
      if (result[i] != gold) {
        printf("R %d seg %d scan result %d\n", rank, seg[i], result[i]);
      }
      EXPECT_EQ(result[i], gold);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

// for reverse
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    seg[0] = (rank % 3 == 2) || (rank == p-1) ? 1 : 0;
    seg[1] = ((rank + 1) % 3 == 2) || (rank == p-1) ? 1 : 0;
    seg[2] = ((rank + 2) % 3 == 2) || (rank == p-1) ? 1 : 0;
  }


  result = mxx2::seg_scan::scan_reverse_elementwise<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    gold = std::max(0, std::min(p-1, rank + 2 - ((rank + i) % 3)));
    if (rank < p-1) {
      if (result[i] != gold) {
        printf("R %d seg %d scan result %d\n", rank, seg[i], result[i]);
      }
      EXPECT_EQ(result[i], gold);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);


  result = mxx2::seg_scan::exscan_reverse_elementwise<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    gold = std::max(0, std::min(p-1, rank + 2 - ((rank + i) % 3)));
    if (((rank + i) % 3 < 2) && (rank < p-1)) {
      if (result[i] != gold) {
        printf("R %d seg %d scan result %d\n", rank, seg[i], result[i]);
      }
      EXPECT_EQ(result[i], gold);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}


TYPED_TEST_P(Mxx2SegmentedMPITest, seg_scan_v)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<typename std::conditional<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, uint8_t, uint>::type> seg(p);
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) {
    for (int i = 0; i < rank; ++i) {
      seg[i] = rank - 1;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank;
    }
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) {
    for (int i = 0; i < rank; ++i) {
      seg[i] = (rank - 1) % 2;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank % 2;
    }
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    for (int i = 0; i < p; ++i) {
      seg[i] = 0;
    }
    seg[rank] = 1;
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

  auto result = mxx2::seg_scan::scanv<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);
  auto print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_scan_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "scanv");
      this->template printvec(gold, "scanv gold");
    }
    EXPECT_TRUE(same);
  }


  result = mxx2::seg_scan::exscanv<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_exscan_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });


    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, true>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "exscanv");
      this->template printvec(gold, "exscanv gold");
    }
    EXPECT_TRUE(same);
  }

// for reverse
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    for (int i = 0; i < p - 1; ++i) {
      seg[i] = 0;
    }
    if (rank > 0) seg[rank-1] = 1;

    allseg = mxx2::gathern(seg, comm, 0);

  }



  result = mxx2::seg_scan::scanv_reverse<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_scan_v_reverse<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });

    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, true, false>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rscanv");
      this->template printvec(gold, "rscanv gold");
    }

    EXPECT_TRUE(same);
  }

  result = mxx2::seg_scan::exscanv_reverse<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::max(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {

    gold = this->template seq_seg_exscan_v_reverse<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, [](int const &x, int const &y){ return std::max(x, y); });


    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, true, true>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "rexscanv");
      this->template printvec(gold, "rexscanv gold");
    }

    EXPECT_TRUE(same);

  }
}

TYPED_TEST_P(Mxx2SegmentedMPITest, seg_reduce)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  typename std::conditional<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, uint8_t, uint>::type seg;
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) {
    seg = rank / 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) {
    seg = (rank / 3) % 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    seg = (rank % 3 == 0) ? 1 : 0;
  }

  auto allseg = mxx2::gather(seg, comm, 0);

  auto vals = mxx2::gather(rank, comm, 0);


  int result = mxx2::seg_reduce::reduce<TypeParam>(rank, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  auto outs = mxx2::gather(result, comm, 0);

  if (rank == 0) {

    auto golds = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

    bool same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(outs, golds, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(vals, "input");
      this->template printvec(outs, "reduce");
      this->template printvec(golds, "reduce gold");
    }

    EXPECT_TRUE(same);

  }


  result = mxx2::seg_reduce::reduce<TypeParam>(rank, seg, std::plus<int>(), comm);

  outs = mxx2::gather(result, comm, 0);

  if (rank == 0) {

    auto golds = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(vals, allseg, std::plus<int>() );

    bool same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(outs, golds, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(vals, "input");
      this->template printvec(outs, "reduce");
      this->template printvec(golds, "reduce gold");
    }


    EXPECT_TRUE(same);

  }
  MPI_Barrier(MPI_COMM_WORLD);

}


TYPED_TEST_P(Mxx2SegmentedMPITest, seg_reduce_n)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<typename std::conditional<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, uint8_t, uint>::type> seg(3);
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) {
    seg[0] = rank / 3;
    seg[1] = (rank + 1) / 3;
    seg[2] = (rank + 2) / 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) {
    seg[0] = (rank / 3) % 3;
    seg[1] = ((rank + 1) / 3) % 3;
    seg[2] = ((rank + 2) / 3) % 3;
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    seg[0] = (rank % 3 == 0) ? 1 : 0;
    seg[1] = ((rank + 1) % 3 == 0) ? 1 : 0;
    seg[2] = ((rank + 2) % 3 == 0) ? 1 : 0;
  }
  std::vector<int> inputs(3);
  inputs[0] = rank;
  inputs[1] = rank;
  inputs[2] = rank;


  int gold, len;
  bool same;
  auto result = mxx2::seg_reduce::reducen<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {

      auto golds = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(vals, allseg, [](int const &x, int const &y){ return std::min(x, y); } );

      same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(outs, golds, allseg);

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


  result = mxx2::seg_reduce::reducen<TypeParam>(inputs, seg, std::plus<int>(), comm);
  for ( int i = 0; i < 3; ++i) {
    auto allseg = mxx2::gather(seg[i], comm, 0);
    auto vals = mxx2::gather(inputs[i], comm, 0);
    auto outs = mxx2::gather(result[i], comm, 0);

    if (rank == 0) {
      auto golds = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(vals, allseg, std::plus<int>() );

      same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(outs, golds, allseg);

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
}


TYPED_TEST_P(Mxx2SegmentedMPITest, seg_reduce_v)
{
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // create segments
  std::vector<typename std::conditional<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, uint8_t, uint>::type> seg(p);
  auto allseg = seg;
  if (std::is_same<TypeParam, mxx2::seg_scan::SegmentUniquelyMarked>::value ) {
    for (int i = 0; i < rank; ++i) {
      seg[i] = rank - 1;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank;
    }
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentFullyMarked>::value ) {
    for (int i = 0; i < rank; ++i) {
      seg[i] = (rank - 1) % 2;
    }
    for (int i = rank; i < p; ++i) {
      seg[i] = rank % 2;
    }
  }
  else if (std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value ) {
    for (int i = 0; i < p; ++i) {
      seg[i] = 0;
    }
    seg[rank] = 1;
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

  auto result = mxx2::seg_reduce::reducev<TypeParam>(inputs, seg, [](int const &x, int const &y){ return std::min(x, y); }, comm);

  MPI_Barrier(MPI_COMM_WORLD);
  auto print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, [](int const &x, int const &y){ return std::min(x, y); });

    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(print, gold, allseg);
    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "reducev");
      this->template printvec(gold, "reducev gold");
    }
    EXPECT_TRUE(same);
  }

  result = mxx2::seg_reduce::reducev<TypeParam>(inputs, seg, std::plus<int>(), comm);

  MPI_Barrier(MPI_COMM_WORLD);
  print = mxx2::gathern(result, comm, 0);
  if (rank == 0) {
    gold = this->template seq_seg_reduce_v<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value>(allinputs, allseg, std::plus<int>());

    same = this->template equal<std::is_same<TypeParam, mxx2::seg_scan::SegmentStartMarked>::value, false, false>(print, gold, allseg);

    if (!same) {
      this->template printvec(allseg, "seg");
      this->template printvec(allinputs, "input");
      this->template printvec(print, "reducev");
      this->template printvec(gold, "reducev gold");
    }

    EXPECT_TRUE(same);
  }

}


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(Mxx2SegmentedMPITest, seg_scan, seg_scan_n, seg_scan_v, seg_reduce, seg_reduce_n, seg_reduce_v);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<mxx2::seg_scan::SegmentFullyMarked, mxx2::seg_scan::SegmentStartMarked, mxx2::seg_scan::SegmentUniquelyMarked > Mxx2SegmentedMPITestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, Mxx2SegmentedMPITest, Mxx2SegmentedMPITestTypes);


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
