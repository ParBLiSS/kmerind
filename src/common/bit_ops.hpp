/**
 * @file    bit_ops.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares and defines constexpr functions for common bit operations
 *          like shifting and masking.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_BIT_OPS_H
#define BLISS_COMMON_BIT_OPS_H

#include <limits>
#include <type_traits>

#include "common/base_types.hpp"

template <typename T>
constexpr T getLeastSignificantBitsMask(const BitSizeType nBits)
{
  // unsigned integer type:  fine.
  // max for signed integer has sign bit (msb) == 0, but digits does not count sign bit.
  // this does not make sense for floating or any other type, so assert.
  static_assert(std::is_integral<T>::value && !std::is_signed<T>::value, "getBitMask only accepts primitive, unsigned integral type.");

  // if T is signed, digits would not include the sign bit, so shift is 1 less than required.
  //return static_cast<T>(std::numeric_limits<T>::max() >> (std::numeric_limits<T>::digits - nBits));

  // no problem with T being signed, if nBits > 1.  may be ever so slightly faster?  conditional because clang optimized seems to have trouble with bitshifting entire word.
  // negation convert Tmax << nBits to int or long int by zero padding, (e.g. 0x0000FFF0 for 16 bit, nBits = 4), then negate.
  // (0xFFFF000F).  this results in an int larger than T.   that and inlining results in 'integer conversion resulted in truncation'.
  // use static_cast to truncate it explicitly.
  return (static_cast<size_t>(nBits) == sizeof(T) * 8UL) ? std::numeric_limits<T>::max() : static_cast<T>(~(std::numeric_limits<T>::max() << nBits));
  //return static_cast<T>((static_cast<T>(0x1) << nBits) - static_cast<T>(1));
}

template <typename T>
constexpr T getMostSignificantBitsMask(const BitSizeType nBits)
{
  // unsigned integer type:  fine.
  // max for signed integer has sign bit (msb) == 0, but digits does not count sign bit.
  // this does not make sense for floating or any other type, so assert.
  static_assert(std::is_integral<T>::value && !std::is_signed<T>::value, "getBitMask only accepts primitive, unsigned integral type.");

  //conditional because clang optimized seems to have trouble with bitshifting entire word.
  return (static_cast<size_t>(nBits) == sizeof(T) * 8UL) ? std::numeric_limits<T>::max() : (std::numeric_limits<T>::max() << (std::numeric_limits<T>::digits - nBits));
  //return static_cast<T>((static_cast<T>(0x1) << nBits) - static_cast<T>(1));
}

template <typename T>
constexpr T getBitsMask(const BitSizeType nBits, const BitSizeType leastSignficantBitsOffset)
{
  return getLeastSignificantBitsMask<T>(nBits) << leastSignficantBitsOffset;
}


constexpr WordType getWordBitMask(const BitSizeType nBits)
{
  //return static_cast<WordType>((0x1 << nBits) - 1);
  return getLeastSignificantBitsMask<WordType>(nBits);
}

constexpr WordType getWordBitMask(const BitSizeType nBits, const BitSizeType offset)
{
  //return static_cast<WordType>(static_cast<WordType>((0x1 << nBits) - 1) << (offset));
  return getBitsMask<WordType>(nBits, offset);
}

constexpr CharType getCharBitMask(const BitSizeType nBits)
{
  return getLeastSignificantBitsMask<CharType>(nBits);
}

constexpr CharType getCharBitMask(const BitSizeType nBits, const BitSizeType offset)
{
  return getBitsMask<CharType>(nBits, offset);
}





template <typename T>
void copyBits(T& to, const T& from, const BitSizeType bits)
{
  const T mask = getLeastSignificantBitsMask<T>(bits);
  to = ((~mask) & to) | (mask & from);
}

template <typename T, unsigned int bits>
void copyBitsFixed(T& to, const T& from)
{
  constexpr T mask = getLeastSignificantBitsMask<T>(bits);
  to = ((~mask) & to) | (mask & from);
}


// NOTE: for ICPC, constexpr function parameters cannot be a reference.  has to be a copy from a constexpr value defined elsewhere.
// NOTE: for ICPC, constexpr function cannot contain ternary conditional.
// NOTE: for ICPC, not sure if conditional evaluation using value from another conditional evaluation is allowed for constexpr functions.

//#define roundDownToMultiple(value, base) value - (value % base);
template <typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, T>::type roundDownToMultiple(const T value, const T base)
{
  return value - (value % base);
}

//#define roundUpToMultiple(value, base) (value - 1) - ((value - 1) % base) + base;
template <typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, T>::type roundUpToMultiple(const T value, const T base)
{
  return roundDownToMultiple(value - 1, base) + base;
}

//#define intCeil(divident, divisor)   (divident == 0) ? 0 : ((divident - 1) / divisor) + 1;
template <typename T>
inline constexpr typename std::enable_if<std::is_integral<T>::value, T>::type intCeil(const T divident, const T divisor)
{
//  return (divident == 0) ? 0 : ((divident - 1) / divisor) + 1;

	return (divident + divisor - 1) / divisor;

	// TODO: proper fix for checking if divisor is 0.
	// TODO: not sure why the code below does not work.
	//return (divisor == 0) ? static_cast<T>(0) : static_cast<T>(divident + divisor - 1) / divisor;
}

constexpr unsigned Log2(const unsigned int n, const unsigned p = 0) {

  return (n <= 1) ? p : Log2(n >> 1, p + 1);
}

constexpr unsigned ceilLog2(const unsigned n_chars)
{
    return Log2(n_chars - 1) + 1;
}

#endif // BLISS_COMMON_BIT_OPS_H
