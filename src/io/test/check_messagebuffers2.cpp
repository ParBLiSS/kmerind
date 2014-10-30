/**
 * @file		check_messagebuffers.cpp
 * @ingroup
 * @author	tpan
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

// TODO: replace content here to use messagebuffers.


#include <unistd.h>  // for usleep

#include "io/message_buffers.hpp"
#include "concurrent/threadsafe_queue.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand
#include <algorithm>

int repeats = 11;
int bufferSize = 2047;



template<typename BuffersType>
void testPool(BuffersType && buffers, const std::string &name, int nthreads) {

  printf("*** TESTING %s: ntargets = %lu, pool threads %d\n", name.c_str(), buffers.getSize(), nthreads);


  printf("TEST append until full: ");
  typedef typename BuffersType::BufferPtrType BufferPtrType;
  std::pair<bool, BufferPtrType > result(false, std::move(BufferPtrType()));
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType> fullBuffers;

  //printf("test string is \"%s\", length %lu\n", data.c_str(), sizeof(int));
  int id = 0;

  int i;
  int count = 0;
  int count2 = 0;
  int count3 = 0;
  int count4 = 0;
  int count5 = 0;


#pragma omp parallel for num_threads(nthreads) default(none) private(i, result) firstprivate(id) shared(buffers, fullBuffers, repeats, bufferSize) reduction(+ : count, count2, count3, count4, count5)
  for (i = 0; i < repeats; ++i) {
    //printf("insert %lu chars into %d\n", sizeof(int), id);
    int data = i;
    result = buffers.append(&data, sizeof(int), id);

    if (result.first) {
      ++count; // success
    } else {
      ++count2; // failure


		if (result.second) {
		  ++count3;     // full buffer
		} else {
			++count4;   // failed insert and no full buffer
		}
    }

    if (result.second)
      fullBuffers.waitAndPush(std::move(result.second));  // full buffer
  }
//  if ((count + count2) != repeats) printf("\nFAIL: number of successful inserts should be %d.  actual %d", repeats, count);
//  else if (count2 != 0) printf("\nFAIL: number of failed insert overall should be 0. actual %d", count2);
//  else
//    if (!(count5 <= count && count <= (count5 + count3))) printf("\nFAIL: number of successful inserts should be close to successful inserts without full buffers.");

  // compute 2 expected since append returns full buffer only on failed insert and if there are n inserts that brings it to just before full, then it depends on timing
  // as to when the buffer becomes full.  (other threads will fail on append until swap happens, but not return a full buffer.)
  int expectedFull = count / (bufferSize/sizeof(int)) - (count % (bufferSize/sizeof(int)) == 0 ? 1 : 0);
  int expectedFull2 = count / (bufferSize/sizeof(int));

  if (fullBuffers.getSize() != count3) printf("\nFAIL: number of full Buffers do not match: fullbuffer size %ld  full count %d", fullBuffers.getSize(), count3);
  // buffer at 23 entries (86 bytes each, 2048 bytes per buffer) will not show as full until the next iterator.
  else if (count3 != expectedFull && count3 != expectedFull2) printf("\nFAIL: number of full Buffers is not right: %d should be %d or %d", count3, expectedFull, expectedFull2);
  //else if (count4 != 0) printf("\nFAIL: number of failed insert due to no buffer should be 0. actual %d", count4);
  else printf("PASS");
  printf("\n");
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu.  numFullBuffers = %d.  num failed append due to no buffer = %d, successful insert and no buffer %d\n", count2, count, fullBuffers.getSize(), count3, count4, count5);



  printf("TEST release: ");
  int count1 = 0;
  count2 = 0;
  count4 = 0;
  count5 = 0;
//  printf("buffer ids 0 to %lu initially in use.\n", buffers.getSize() - 1);
//  printf("releasing: ");


  int iterations = count3 + 10;
#pragma omp parallel for num_threads(nthreads) default(none) private(id, result) shared(buffers, fullBuffers, repeats, bufferSize, iterations) reduction(+ : count,count1, count2, count3, count4, count5)
  for (id = 0; id < iterations; ++id) {
    try {
      result = fullBuffers.tryPop();
      if (result.first) {
//        printf("%d ", result.second);
        if (result.second) {

          buffers.releaseBuffer(std::move(result.second));
          ++count1;    // successful pop
        } else {
          ++count2;   // successful pop but no actual buffer to release
        }
      } else {
        ++count5;     // failed pop.
      }
    } catch(const std::invalid_argument & e)
    {
      printf("\nFAIL with %s", e.what());
      ++count4;       // error during pop
    }
  }

  expectedFull = count / (bufferSize/sizeof(int)) - (count % (bufferSize/sizeof(int)) == 0 ? 1 : 0);
  expectedFull2 = count / (bufferSize/sizeof(int));

  if (count4 != 0) printf("\nFAIL: invalid argument exception during pop.  count = %d", count4);
  else if (count5 != 10) printf("\nFAIL: failed on pop %d times, expected 10", count5);
  else if (count2 != 0) printf("\nFAIL: succeeded in pop but not full buffer. %d", count2);
  else if (count1 != expectedFull && count1 != expectedFull2) printf("FAIL: expected %d or %d full buffers, but received %d", expectedFull, expectedFull2, count1);
  else if (count1 != count3) printf("\nFAIL: successful pops. expected %d.  actual %d", count3, count1);
  else
    printf("PASS");
  printf("\n");




  buffers.reset();

  printf("TEST all operations together: ");
  count3 = 0;
  count1 = 0;
  count2 = 0;
  count4 = 0;
  count = 0;
  count5 = 0;
  int count6 = 0;
  int count7 = 0;
  //printf("full buffer: ");

  std::vector< std::vector<int> > stored(nthreads);
  std::vector< std::vector<int> > appended(nthreads);


  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, result) firstprivate(id) shared(stored, appended, buffers, repeats, bufferSize) reduction(+ : count, count1, count2, count3,count4, count6, count7)
  for (i = 0; i < repeats; ++i) {
    int data = i;
        result = buffers.append(&data, sizeof(int), id);

    if (result.first) {
      ++count1;
      appended[omp_get_thread_num()].push_back(data);


      if (result.second) {
        throw std::logic_error("ERROR: append result is true-true.  should not get here.\n");

      }
    } else {
      ++count2;

      if (result.second) {
        ++count4;
       // count7 = count1;  // save the number of successful inserts so far.

        result.second->waitForAllUpdates();

        count += result.second->getFinalSize();

        stored[omp_get_thread_num()].insert(stored[omp_get_thread_num()].end(), result.second->operator int*(), result.second->operator int*() + result.second->getFinalSize());
        buffers.releaseBuffer(std::move(result.second));
      }
    }
  }

  // concatenate the vectors
  std::vector<int> allstored;
  std::vector<int> allappended;
  for (int k = 0; k < nthreads; ++k) {
    allstored.insert(allstored.end(), stored[k].begin(), stored[k].end());
    allappended.insert(allappended.end(), appended[k].begin(), appended[k].end());
  }


  BufferPtrType final = buffers.flushBufferForRank(id);  // flush blocks buffer and waits for all updates., but need to set final size.
  count5 = final->getFinalSize();

  allstored.insert(allstored.end(), final->operator int*(),  final->operator int*() + final->getFinalSize());
  buffers.releaseBuffer(std::move(final));

  //if (count7 != count/sizeof(int)) printf("\nFAIL: append count = %d, actual data inserted is %ld", count7, count/sizeof(int) );

  final = buffers.flushBufferForRank(id);
  count6 = final->getFinalSize();
  if (count6 > 0) printf("\nFAIL: received %d data after flush.", count6);
  buffers.releaseBuffer(std::move(final));


  // now sort it an check to see if we are missing anything
  std::sort(allstored.begin(), allstored.end());
  std::sort(allappended.begin(), allappended.end());
  std::vector<int> appendedButNotStored; std::vector<int> storedButNotAppended;
  std::set_difference(allappended.begin(), allappended.end(), allstored.begin(), allstored.end(), std::inserter(appendedButNotStored, appendedButNotStored.begin()));
  std::set_difference(allstored.begin(), allstored.end(), allappended.begin(), allappended.end(), std::inserter(storedButNotAppended, storedButNotAppended.begin()));



  if (appendedButNotStored.size() > 0) {
    printf("\nFAIL!: %ld elements appended but not stored: ", appendedButNotStored.size());
    for (i = 0; i < appendedButNotStored.size(); ++i) {
      printf("%d, ", appendedButNotStored[i]);
    }
  }
  if (storedButNotAppended.size() > 0) {
    printf("\nFAIL?: %ld elements stored but not appended: ", storedButNotAppended.size());
    for (i = 0; i < storedButNotAppended.size(); ++i) {
      printf("%d, ", storedButNotAppended[i]);
    }
  }


  if ((count + count5) != count1 * sizeof(int)) {
    printf("\nFAIL: total bytes %d (%d + %d) for %ld entries.  expected %d entries.\n", (count + count5), count , count5, (count + count5)/sizeof(int), count1);

    throw std::logic_error("missing entries.");
  }




  expectedFull = (count1 - 1 + (bufferSize/sizeof(int))) / (bufferSize/sizeof(int)) - (count1 % (bufferSize/sizeof(int)) == 0 ? 1 : 0);
  expectedFull2 = count1 / (bufferSize/sizeof(int));

//  if (count1 != repeats) printf("\nFAIL: number of successful inserts should be %d.  actual %d", repeats, count1);
//  else if (count2 != 0) printf("\nFAIL: number of failed insert overall should be 0. actual %d", count2);
//  else
  if (count3 != 0) printf("\nFAIL: number of full Buffers from successful insert is not right: %d should be 0", count3);
  else if (count4 != expectedFull && count4 != expectedFull2) {
    printf("\nFAIL: number of full Buffers from failed insert is not right: %d should be %d or %d", count4, expectedFull, expectedFull2);
    fflush(stdout);
  }
  else printf("PASS");
  printf("\n");

  //printf("\n");
  printf("Number of failed attempt to append to buffer is %d, success %d. full buffers size: %lu, released successful appendss %d, released failed appends %d\n", count2, count1, fullBuffers.getSize(), count3, count4);

};


int main(int argc, char** argv) {

  // construct, acquire, access, release
  if (argc > 1) {
    repeats = atoi(argv[1]);
  }

  /// thread unsafe.  test in single thread way.
//while(true) {

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(1, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(2, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(3, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(4, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(5, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(6, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(7, 2047)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(8, 2047)), "thread unsafe buffers", 1);



  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2047)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2047)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2047)), "thread safe buffers", 4);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2047)), "thread safe buffers", 8);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2047)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2047)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2047)), "thread safe buffers", 4);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2047)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2047)), "thread safe buffers", 3);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2047)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(5, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(5, 2047)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(6, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(6, 2047)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(7, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(7, 2047)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(8, 2047)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(8, 2047)), "thread safe buffers", 2);


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
