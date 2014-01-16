/**
 * indexing.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: tpan
 *
 *  Description:  testing code to try different ways of generating the index key from a collection of strings
 *
 */

#include <vector>
#include <random>
#include <cmath>

#include "utils/logging.h"

/**
 * parameters
 */
const uint8_t K = 23;
const uint8_t BITS = 2;
const uint8_t C = 4;
const uint32_t R = 100000;
const double L_MEAN = 100.0f;
const double L_STDEV = 20.0f;

typedef std::vector<uint32_t> SequenceT;

/**
 * allocate input array
 */
std::vector<SequenceT> generateStrings() {
  // create the array of arrays
  INFO("arrays created");
  std::vector<SequenceT> reads;

  std::vector<uint32_t> lengths;

  std::default_random_engine generator1;
  std::normal_distribution<double> distribution1(L_MEAN, L_STDEV);

  std::default_random_engine generator2;
  std::uniform_int_distribution<double> distribution2;

  // use a random number generator to generate the "packed string"
  for (int i = 0; i < R; ++i)
  {
    // generate the random lengths of the reads
    lengths[i] = static_cast<uint32_t>(round(distribution1(generator1)));
    SequenceT read;
    int idx = 0;
    for (int j = 0; j < lengths[i]; j +=16)
    {
      // pack in 32 bit.
      read[idx] = distribution2(generator2);
      ++idx;
    }

    reads[i] = read;
  }

}


/**
 * output arrays
 */
std::vector<uint64_t> initializeKeyArray() {
  // initialize to 0
}



template<int CORE, int SIMD>
void computeKeys()
{
  FATAL("BASE TEMPLATE CALLED. NOT IMPLEMENTED.");
}

#define SERIAL 0
#define OPENMP 1

#define SCALAR 0
#define SSE2 1

/**
 * serial version
 */
template<>
void computeKeys<SERIAL, SCALAR>() {

}
template<>
void computeKeys<SERIAL, SSE2>() {

}

/**
 * openMP version
 */
template<>
void computeKeys<OPENMP, SCALAR>() {

}

template<>
void computeKeys<OPENMP, SSE2>() {

}



/**
 * main function.
 */
int main(int argc, char* argv[]) {




	return 0;
}
