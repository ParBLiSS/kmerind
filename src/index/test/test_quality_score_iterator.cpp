/**
 * @file    test_quality_score.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

// include classes to test
#include <index/quality_score_iterator.hpp>
#include <vector>
#include <algorithm>
#include <cassert>

//// Usable AlmostEqual function
//// from http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
//// set maxUlps to be 1M or less.
//bool AlmostEqual2sComplement(float A, float B, int maxUlps)
//{
//    // Make sure maxUlps is non-negative and small enough that the
//    // default NAN won't compare as equal to anything.
//
//    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//
//    int aInt = *(int*)&A;
//
//    // Make aInt lexicographically ordered as a twos-complement int
//    if (aInt < 0)
//        aInt = 0x80000000 - aInt;
//
//    // Make bInt lexicographically ordered as a twos-complement int
//    int bInt = *(int*)&B;
//
//    if (bInt < 0)
//        bInt = 0x80000000 - bInt;
//
//    int intDiff = abs(aInt - bInt);
//
//    if (intDiff <= maxUlps)
//        return true;
//
//    return false;
//}
//
//bool AlmostEqual2sComplement(double A, double B, int maxUlps)
//{
//    // Make sure maxUlps is non-negative and small enough that the
//    // default NAN won't compare as equal to anything.
//
//    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//
//    long int aInt = *(long int*)&A;
//
//    // Make aInt lexicographically ordered as a twos-complement int
//    if (aInt < 0)
//        aInt = 0x8000000000000000 - aInt;
//
//    // Make bInt lexicographically ordered as a twos-complement int
//    long int bInt = *(long int*)&B;
//
//    if (bInt < 0)
//        bInt = 0x8000000000000000 - bInt;
//
//    long int intDiff = abs(aInt - bInt);
//
//    if (intDiff <= maxUlps)
//        return true;
//
//    return false;
//}


template <typename T>
typename std::enable_if<std::is_same<T, double>::value, bool>::type compare_vectors(const std::vector<T> & first, const std::vector<T> & second) {
  if (first == second) {
    return true;   // if same object
  }
  if (first.size() != second.size()) {
    std::cout << "size not the same" << std::endl;
    return false;  // if size differ
  }

  if (std::equal(first.cbegin(), first.cend(), second.cbegin(), [](T x, T y) {
    return fabs(x - y) < 1.0e-14; }
  )) return true;
  else {
    std::cout << "Vectors not the same" << std::endl;
    for (int i = 0; i < first.size(); ++i) {
      std::cout << (first[i] - second[i]) << std::endl;
    }

    return false;
  }
}

template <typename T>
typename std::enable_if<std::is_same<T, float>::value, bool>::type compare_vectors(const std::vector<T> & first, const std::vector<T> & second) {
  if (first == second) {
    return true;   // if same object
  }
  if (first.size() != second.size()) {
    std::cout << "size not the same" << std::endl;
    return false;  // if size differ
  }

  if (std::equal(first.cbegin(), first.cend(), second.cbegin(), [](T x, T y) {
    return fabs(x - y) < 1.0e-6; }
  )) return true;
  else {
    std::cout << "Vectors not the same" << std::endl;
    for (int i = 0; i < first.size(); ++i) {
      std::cout << (first[i] - second[i]) << std::endl;
    }

    return false;
  }
}

template <typename T>
typename std::enable_if<!std::is_floating_point<T>::value, bool>::type compare_vectors(const std::vector<T> & first, const std::vector<T> & second) {
  if (first == second) return true;   // if same object
  if (first.size() != second.size()) return false;  // if size differ

  if (std::equal(first.cbegin(), first.cend(), second.cbegin())) return true;
  else {

    std::cout << "Vectors not the same" << std::endl;
    for (int i = 0; i < first.size(); ++i) {
      std::cout << first[i] - second[i] << std::endl;
    }

    return false;
  }
}


// templated test function
template<typename CODEC, unsigned int K>
void iter_decode(const std::vector<unsigned char>& data, // quality score value
                           std::vector<typename CODEC::value_type>& output) {

  using QualIter = bliss::index::QualityScoreGenerationIterator<std::vector<unsigned char>::const_iterator,
        K, CODEC>;

  QualIter it(data.begin(), true);
  QualIter end(data.end(), false);

  output.clear();

  for (; it != end; ++it) {
    output.push_back(*it);
  }

}


// templated test function
template<typename CODEC, unsigned int K>
void codec_decode(const std::vector<unsigned char> & data, // quality score value
                           std::vector<typename CODEC::value_type>& output) {

  CODEC decoder;

//  std::cout << "decoder LUT SIZE: " << (int)CODEC::size << std::endl;
//  for (int i = 0; i < CODEC::size; ++i) {
//    std::cout << CODEC::lut[i] << ", ";
//  }
//  std::cout << std::endl;

  output.clear();

  using OT = typename CODEC::value_type;

  int invalid = 0;
  OT sum = 0.0;
  for (int i = 0; i < K; ++i) {
    auto tmp = decoder.decode(data[i]);

    if (tmp >= CODEC::min_log_prob)
      sum += tmp;
    else
      ++invalid;
  }
  if (invalid > 0) output.push_back(0);
  else output.push_back(std::exp2(sum));

  for (int i = K; i < data.size(); ++i) {
    auto tmp = decoder.decode(data[i-K]);
    if (tmp >= CODEC::min_log_prob)
      sum -= tmp;
    else
      --invalid;

    tmp = decoder.decode(data[i]);
    if (tmp >= CODEC::min_log_prob)
      sum += tmp;
    else
      ++invalid;

    if (invalid > 0) output.push_back(0);
    else output.push_back(std::exp2(sum));
  }

}




template<typename OT, unsigned char MinInput, unsigned char MaxInput, char MinScore, unsigned int K>
void direct_decode(const std::vector<unsigned char> & data,
                   std::vector<OT>& output) {

  output.clear();

  OT sum = 0;
  int invalid= 0;
  for (int i = 0; i < K; ++i) {
    auto tmp = data[i] - MinInput + MinScore;

    if (tmp > 0)
      sum += std::log2(1.0 - std::exp2(static_cast<OT>(tmp) * log2(10.0) / -10.0));
    else
      ++invalid;
  }
  if (invalid > 0) output.push_back(0);
  else output.push_back(std::exp2(sum));

  for (int i = K; i < data.size(); ++i) {
    auto tmp = data[i-K] - MinInput + MinScore;
    if (tmp > 0)
      sum -= std::log2(1.0 - std::exp2(static_cast<OT>(tmp) * log2(10.0) / -10.0));
    else
      --invalid;

    tmp = data[i] - MinInput + MinScore;
    if (tmp > 0)
      sum += std::log2(1.0 - std::exp2(static_cast<OT>(tmp) * log2(10.0) / -10.0));
    else
      ++invalid;

    if (invalid > 0) output.push_back(0);
    else output.push_back(std::exp2(sum));
  }
}



template<typename OT, unsigned char MinInput, unsigned char MaxInput, char MinScore, unsigned int K>
void testQualityIterator(const std::string & strdata) {

  std::vector<unsigned char> gold(strdata.begin(), strdata.end());

  using Encoder = bliss::index::QualityScoreCodec<OT, MinInput, MaxInput, MinScore>;


  // test process.
  std::vector<OT> goldDecoded;
  direct_decode<OT, MinInput, MaxInput, MinScore, K>(gold, goldDecoded);



  // test codec decoding.
  std::vector<OT> codecDecoded;
  codec_decode< Encoder, K >(gold, codecDecoded);
  bool same = compare_vectors<OT>(codecDecoded, goldDecoded);

  if (!same) {
    std::cout << "codec decode: result not same" << std::endl;

    std::cout << std::endl<< "GOLD decoded: size: " << goldDecoded.size() << std::endl;
    std::copy(goldDecoded.begin() , goldDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
    std::cout << std::endl<< "codec decoded: size: " << codecDecoded.size() << std::endl;
    std::copy(codecDecoded.begin() , codecDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
  }

  EXPECT_TRUE(same);


  std::vector<OT> iterDecoded;
  iter_decode< Encoder, K >(gold, iterDecoded);
  same = compare_vectors<OT>(iterDecoded, goldDecoded);

  if (!same) {
    std::cout << "iter decode: result not same" << std::endl;

    std::cout << std::endl<< "GOLD decoded: size: " << goldDecoded.size() << std::endl;
    std::copy(goldDecoded.begin() , goldDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
    std::cout << std::endl<< "iter decoded: size: " << iterDecoded.size() << std::endl;
    std::copy(iterDecoded.begin() , iterDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
  }

  EXPECT_TRUE(same);


}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGenerationIteratorTest, TestIllunima18QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testQualityIterator<double, 33, 126, 0, 3 >(strdata);
  testQualityIterator<double, 33, 126, 0, 7 >(strdata);
  testQualityIterator<double, 33, 126, 0, 15>(strdata);
  testQualityIterator<double, 33, 126, 0, 31>(strdata);
}




/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGenerationIteratorTest, TestIllunima13QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testQualityIterator<double, 64, 126, 0, 3 >(strdata);
  testQualityIterator<double, 64, 126, 0, 7 >(strdata);
  testQualityIterator<double, 64, 126, 0, 15>(strdata);
  testQualityIterator<double, 64, 126, 0, 31>(strdata);
}



/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGenerationIteratorTest, TestIllunima15QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "CDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testQualityIterator<double, 67, 126, 3, 3 >(strdata);
  testQualityIterator<double, 67, 126, 3, 7 >(strdata);
  testQualityIterator<double, 67, 126, 3, 15>(strdata);
  testQualityIterator<double, 67, 126, 3, 31>(strdata);

}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGenerationIteratorTest, TestIllunima18QualityScoreFloat)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testQualityIterator<float, 33, 126, 0, 3 >(strdata);
  testQualityIterator<float, 33, 126, 0, 7 >(strdata);
  testQualityIterator<float, 33, 126, 0, 15>(strdata);
  testQualityIterator<float, 33, 126, 0, 31>(strdata);

}


//
//
///**
// * Test k-mer generation with 2 bits for each character
// */
//TEST(QualityScoreGenerationIteratorTest, TestIllunima13QualityScoreFloat)
//{
//  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
//  std::string strdata = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
//
//  testQualityIterator<float, 64, 126, 0, 3 >(strdata);
//  testQualityIterator<float, 64, 126, 0, 7 >(strdata);
//  testQualityIterator<float, 64, 126, 0, 15>(strdata);
//  testQualityIterator<float, 64, 126, 0, 31>(strdata);
//
//}
//
//
//
///**
// * Test k-mer generation with 2 bits for each character
// */
//TEST(QualityScoreGenerationIteratorTest, TestIllunima15QualityScoreFloat)
//{
//  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
//  std::string strdata = "CDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
//
//  testQualityIterator<float, 67, 126, 3, 3 >(strdata);
//  testQualityIterator<float, 67, 126, 3, 7 >(strdata);
//  testQualityIterator<float, 67, 126, 3, 15>(strdata);
//  testQualityIterator<float, 67, 126, 3, 31>(strdata);
//
//}
//
