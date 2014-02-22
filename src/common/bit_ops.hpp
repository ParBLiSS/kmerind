/**
 * @file    bit_ops.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares and defines constexpr functions for common bit operations
 *          like shifting and masking.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_BIT_OPS_H
#define BLISS_COMMON_BIT_OPS_H

#include <common/base_types.hpp>

constexpr WordType getWordBitMask(const BitSizeType nBits)
{
  return static_cast<WordType>((0x1 << nBits) - 1);
}

constexpr WordType getWordBitMask(const BitSizeType nBits, const BitSizeType offset)
{
  return static_cast<WordType>(static_cast<WordType>((0x1 << nBits) - 1) << (offset));
}

constexpr CharType getCharBitMask(const BitSizeType nBits)
{
  return static_cast<CharType>((0x1 << nBits) - 1);
}

constexpr WordType getCharBitMask(const BitSizeType nBits, const BitSizeType offset)
{
  return static_cast<CharType>(static_cast<CharType>((0x1 << nBits) - 1) << (offset));
}

template <typename T>
constexpr T getBitMask(const BitSizeType nBits)
{
  return static_cast<T>((static_cast<T>(0x1) << nBits) - static_cast<T>(1));
}

template <typename T>
constexpr T roundDownToMultiple(const T& value, const T& base)
{
  return value - (value % base);
}

template <typename T>
constexpr T roundUpToMultiple(const T& value, const T& base)
{
  return roundDownToMultiple(value - 1, base) + base;
}

template <typename T>
constexpr T intCeil(const T& divident, const T& divisor)
{
  return ((divident - 1) / divisor) + 1;
}

constexpr unsigned Log2(unsigned int n, unsigned p = 0) {

  return (n <= 1) ? p : Log2(n >> 1, p + 1);
}

constexpr unsigned ceilLog2(unsigned n_chars)
{
    return Log2(n_chars - 1) + 1;
}

#endif // BLISS_COMMON_BIT_OPS_H
