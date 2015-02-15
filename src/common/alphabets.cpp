/**
 * @file    alphabets.cpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Defines all common alphabets, including DNA, DNA5, RNA, RNA5
 *          AA (IUPAC), DNA_IUPAC, and CUSTOM
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#include <common/alphabets.hpp>

// initialize the translation tables in the .cpp file for correct linkage
constexpr uint8_t DNA::FROM_ASCII[256];
constexpr char DNA::TO_ASCII[DNA::SIZE];
constexpr uint8_t DNA::TO_COMPLEMENT[DNA::SIZE];

constexpr uint8_t DNA5::FROM_ASCII[256];
constexpr char DNA5::TO_ASCII[DNA5::SIZE];
constexpr uint8_t DNA5::TO_COMPLEMENT[DNA5::SIZE];

constexpr uint8_t DNA16::FROM_ASCII[256];
constexpr char DNA16::TO_ASCII[DNA16::SIZE];
constexpr uint8_t DNA16::TO_COMPLEMENT[DNA16::SIZE];
