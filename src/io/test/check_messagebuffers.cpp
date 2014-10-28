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

int repeats = 11;
int bufferSize = 2048;
std::string data("this is a test.  this a test of the emergency broadcast system.  this is only a test. ");


template<typename BuffersType>
void testPool(BuffersType && buffers, const std::string &name, int nthreads) {

  printf("*** TESTING %s: ntargets = %lu, pool threads %d\n", name.c_str(), buffers.getSize(), nthreads);


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
#pragma omp parallel for num_threads(nthreads) default(none) private(id, result) shared(buffers, data, fullBuffers, repeats, bufferSize) reduction(+ : count,count1, count2, count3, count4, count5)
  for (id = 0; id < 350; ++id) {
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
  else if (count5 != 350 - expectedFull && count5 != 350 - expectedFull2) printf("\nFAIL: failed on pop %d times", count5);
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
  //printf("full buffer: ");
  id = 0;
#pragma omp parallel for num_threads(nthreads) default(none) private(i, result) firstprivate(id) shared(buffers, data, repeats, bufferSize) reduction(+ : count, count1, count2, count3,count4)
  for (i = 0; i < repeats; ++i) {
    result = buffers.append(data.c_str(), data.length(), id);

    if (result.first) {
      ++count1;

      if (result.second) {
        ++count3;
        count += result.second->getFinalSize();
        buffers.releaseBuffer(std::move(result.second));
      }
    } else {
      ++count2;

      if (result.second) {
        ++count4;
        bool updating = result.second->isUpdating();
        if (updating) printf("  FULLBUFFER: size %d updating ? %s, blocked? %s\n", result.second->getFinalSize(), (updating ? "Y" : "N"), (result.second->isBlocked() ? "Y" : "N"));

        count += result.second->getFinalSize();

        buffers.releaseBuffer(std::move(result.second));
      }
    }
  }
  BufferPtrType final = buffers.flushBufferForRank(id);
  count5 = final->getApproximateSize();
  if ((count + count5) != count1 * data.length()) {
    printf("\nFAIL: number of entries in full %d and final buffers %d is not the same as number successfully inserted %ld bytes for %d entries.\n", count , count5 , count1 * data.length(), count1);

    printf("    content length %ld: %s\n", strlen(reinterpret_cast<char*>(final->getData())), reinterpret_cast<char*>(final->getData()));
  }
  buffers.releaseBuffer(std::move(final));

  expectedFull = (count1 - 1 + (bufferSize/data.length())) / (bufferSize/data.length()) - (count1 % (bufferSize/data.length()) == 0 ? 1 : 0);
  expectedFull2 = count1 / (bufferSize/data.length());

//  if (count1 != repeats) printf("\nFAIL: number of successful inserts should be %d.  actual %d", repeats, count1);
//  else if (count2 != 0) printf("\nFAIL: number of failed insert overall should be 0. actual %d", count2);
//  else
  if (count3 != 0) printf("\nFAIL: number of full Buffers from successful insert is not right: %d should be 0", count3);
  else if (count4 != expectedFull && count4 != expectedFull2) {
    printf("\nFAIL: number of full Buffers from failed insert is not right: %d should be %d or %d", count4, expectedFull, expectedFull2);
    throw std::logic_error("GDB stop here!");
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

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(1, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(2, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(3, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(4, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(5, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(6, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(7, 2048)), "thread unsafe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_UNSAFE>(8, 2048)), "thread unsafe buffers", 1);



  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 4);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(1, 2048)), "thread safe buffers", 8);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 3);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(2, 2048)), "thread safe buffers", 4);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 2);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(3, 2048)), "thread safe buffers", 3);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(4, 2048)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(5, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(5, 2048)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(6, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(6, 2048)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(7, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(7, 2048)), "thread safe buffers", 2);

  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(8, 2048)), "thread safe buffers", 1);
  testPool(std::move(bliss::io::SendMessageBuffers<bliss::concurrent::THREAD_SAFE>(8, 2048)), "thread safe buffers", 2);


  // this one is not defined because it's not logical.  not compilable.
  // bliss::io::BufferPool<bliss::concurrent::THREAD_UNSAFE, bliss::concurrent::THREAD_SAFE> tsusPool(8192, 8);


//}



}
