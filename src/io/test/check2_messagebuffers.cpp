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

#include "io/locking_message_buffers.hpp"
#include "concurrent/threadsafe_queue.hpp"
#include "omp.h"
#include <cassert>
#include <chrono>
#include <cstdlib>   // for rand

int repeats = 11;
int bufferSize = 2048;
std::string data("this is a test.  this a test of the emergency broadcast system.  this is only a test. ");


template<typename BuffersType>
void testPool(BuffersType && buffers, bliss::concurrent::LockType poollt, bliss::concurrent::LockType bufferlt, int nthreads) {

  printf("*** TESTING pool lock %d buffer lock %d: ntargets = %lu, pool threads %d\n", poollt, bufferlt, buffers.getSize(), nthreads);


  printf("TEST append until full: ");
  typedef typename BuffersType::BufferPtrType BufferPtrType;
  std::pair<bool, BufferPtrType > result(false, std::move(BufferPtrType()));
  bliss::concurrent::ThreadSafeQueue<typename BuffersType::BufferPtrType> fullBuffers;

  //printf("test string is \"%s\", length %lu\n", data.c_str(), data.length());
  int id = 0;

  int i;
  int count = 0;
  int count2 = 0;
  int count3 = 0;
  int count4 = 0;
  int count5 = 0;


#pragma omp parallel for num_threads(nthreads) default(none) private(i, result) firstprivate(id) shared(buffers, data, fullBuffers, repeats, bufferSize) reduction(+ : count, count2, count3, count4, count5)
  for (i = 0; i < repeats; ++i) {
    //printf("insert %lu chars into %d\n", data.length(), id);

    result = buffers.append(data.c_str(), data.length(), id);

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
  int expectedFull = count / (bufferSize/data.length()) - (count % (bufferSize/data.length()) == 0 ? 1 : 0);
  int expectedFull2 = count / (bufferSize/data.length());

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
#pragma omp parallel for num_threads(nthreads) default(none) private(id, result) shared(buffers, data, fullBuffers, repeats, bufferSize, iterations) reduction(+ : count,count1, count2, count3, count4, count5)
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

  expectedFull = count / (bufferSize/data.length()) - (count % (bufferSize/data.length()) == 0 ? 1 : 0);
  expectedFull2 = count / (bufferSize/data.length());

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
//  int count7 = 0;
  //printf("full buffer: ");
  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, result) firstprivate(id) shared(buffers, data, repeats, bufferSize) reduction(+ : count, count1, count2, count3,count4, count6)
  for (i = 0; i < repeats; ++i) {
    result = buffers.append(data.c_str(), data.length(), id);

    if (result.first) {
      ++count1;

      if (result.second) {
        ++count3;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = result.second->is_writing();
        if (updating) printf("  FULLBUFFER1: size %ld updating? %s, blocked? %s\n", result.second->getSize(), (updating ? "Y" : "N"), (result.second->is_read_only() ? "Y" : "N"));

        count += result.second->getSize();
        count6 += strlen(result.second->operator char*());
        buffers.releaseBuffer(std::move(result.second));
      }
    } else {
      ++count2;

      if (result.second) {
        ++count4;
//        count7 = count1;  // save the number of successful inserts so far.
        bool updating = result.second->is_writing();
        if (updating) printf("  FULLBUFFER: size %ld blocked? %s\n", result.second->getSize(), (result.second->is_read_only() ? "Y" : "N"));

        count += result.second->getSize();
        count6 += strlen(result.second->operator char*());

        buffers.releaseBuffer(std::move(result.second));
      }
    }
  }

  buffers.at(id)->block_and_flush();
  bool updating = buffers.at(id)->is_writing();
  if (updating) printf("  PreFLUSH: size %ld updating? %s, blocked? %s\n", buffers.at(id)->getSize(), (updating ? "Y" : "N"), (buffers.at(id)->is_read_only() ? "Y" : "N"));

  BufferPtrType final = buffers.flushBufferForRank(id);
  count5 = final->getSize();
  if (count != count6) printf("\nFAIL: number of bytes written %d and number of bytes in FinalSize %d are not the same", count6, count);
  if ((count + count5) != count1 * data.length()) {
    printf("\nFAIL: total bytes %d (%d + %d) for %ld entries.  expected %d entries.\n", (count + count5), count , count5, (count + count5)/data.length(), count1);

    printf("    content length %ld\n", strlen(final->operator char*())); //, final->operator char*());
  }
  if ((count + count5)/data.length() + 1 < count1) {
    printf("\nFAIL: missing %ld entries, expected 1", count1 -(count + count5)/data.length());
    throw std::logic_error("missing more than 1 entry.");
  }
  buffers.releaseBuffer(std::move(final));

  //if (count7 != count/data.length()) printf("\nFAIL: append count = %d, actual data inserted is %ld", count7, count/data.length() );


  final = buffers.flushBufferForRank(id);
  count5 = final->getSize();
  if (count5 > 0) printf("\nFAIL: received %d data after flush.", count5);
  buffers.releaseBuffer(std::move(final));




  expectedFull = (count1 - 1 + (bufferSize/data.length())) / (bufferSize/data.length()) - (count1 % (bufferSize/data.length()) == 0 ? 1 : 0);
  expectedFull2 = count1 / (bufferSize/data.length());

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

#ifdef BLISS_MUTEX
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::MUTEX;
#endif
#ifdef BLISS_SPINLOCK
  constexpr bliss::concurrent::LockType lt = bliss::concurrent::LockType::SPINLOCK;
#endif


  /// thread unsafe.  test in single thread way.
//while(true) {

  for (int i = 1; i <= 8; ++i) {
    testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 2047>(i)), bliss::concurrent::LockType::NONE, bliss::concurrent::LockType::NONE, 1);

    for (int j = 1; j <= 8; ++j) {
      testPool(std::move(bliss::io::SendMessageBuffers<lt, bliss::concurrent::LockType::LOCKFREE, 2047>(i)), lt, bliss::concurrent::LockType::LOCKFREE, j);
    }
  }


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
