/**
 * @file    base_types.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares basic typedefs.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_BASE_TYPES_H
#define BLISS_COMMON_BASE_TYPES_H

// C std lib includes:
#include <cstdint>

/// The type for a machine word, the PackedString representation uses this
/// as underlying type.
typedef uint64_t WordType;

/// The character type used by the alphabets
typedef uint8_t CharType;

/// A type for bit sizes (i.e. representing a number of bits)
typedef uint_fast8_t BitSizeType;

/// A type to represent the size of an alphabet
typedef uint_fast16_t AlphabetSizeType;


#endif // BLISS_COMMON_BASE_TYPES_H
