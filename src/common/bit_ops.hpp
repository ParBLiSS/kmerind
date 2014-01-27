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

inline constexpr WordType getWordBitMask(const BitSizeType nBits)
{
    return static_cast<WordType>((0x1 << nBits) - 1);
}

inline constexpr WordType getWordBitMask(const BitSizeType nBits, const BitSizeType offset)
{
    return static_cast<WordType>(static_cast<WordType>((0x1 << nBits) - 1) << (offset));
}

inline constexpr CharType getCharBitMask(const BitSizeType nBits)
{
    return static_cast<CharType>((0x1 << nBits) - 1);
}

inline constexpr WordType getCharBitMask(const BitSizeType nBits, const BitSizeType offset)
{
    return static_cast<CharType>(static_cast<CharType>((0x1 << nBits) - 1) << (offset));
}

constexpr unsigned Log2(unsigned n, unsigned p = 0) {
    return (n <= 1) ? p : Log2(n >> 1, p + 1);
}

constexpr unsigned bits_per_char_needed(unsigned n_chars)
{
    return Log2(n_chars - 1) + 1;
}

#endif // BLISS_COMMON_BIT_OPS_H
