/*
 * bit_test_common.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: tpan
 */

#ifndef BIT_TEST_COMMON_HPP_
#define BIT_TEST_COMMON_HPP_

namespace bliss {
namespace utils {
namespace bit_ops {

namespace test {

template <unsigned char Bits>
struct BitsParam { static constexpr unsigned char bitsPerGroup = Bits; };

template <typename data_type, size_t data_size>
void print(data_type (&w)[data_size]) {
  std::cout << "data type size " << std::dec << sizeof(data_type) << " len " << data_size << ": ";
  for (int64_t k = data_size -1 ; k >= 0; --k) {
    std::cout << std::hex << static_cast<size_t>(w[k]) << " ";
  }
  std::cout << std::endl;
}


template <uint16_t shift, typename data_type, size_t data_size>
bool is_right_shifted(data_type (&in)[data_size], data_type (&out)[data_size])  {

	if (shift == 0) {
		return std::equal(in, in + data_size, out);
	}

	data_type prev = 0;
	uint16_t inv_shift = sizeof(data_type) * 8 - shift;

	bool result = true;

	for (int64_t i = data_size - 1; i >= 0; --i) {
		result &= (out[i] == static_cast<data_type>((in[i] >> shift) | prev));

		if (!result) {
			std::cout << "word " << std::dec << i << " out = " << std::hex <<
					static_cast<size_t>(out[i]) << " in = " <<
					static_cast<size_t>(in[i]) << " gold = " <<
					static_cast<size_t>((in[i] >> shift) | prev) << std::endl;
		}


		prev = in[i] << inv_shift;
	}

	return result;
}

template <uint16_t shift, typename data_type, size_t data_size>
bool is_left_shifted(data_type (&in)[data_size], data_type (&out)[data_size])  {

	if (shift == 0) {
		return std::equal(in, in + data_size, out);
	}

	data_type prev = 0;
	uint16_t inv_shift = sizeof(data_type) * 8 - shift;

	bool result = true;

	for (size_t i = 0; i < data_size; ++i) {
		result &= (out[i] == static_cast<data_type>((in[i] << shift) | prev));

		if (!result) {
			std::cout << "word " << std::dec << i << " out = " << std::hex <<
					static_cast<size_t>(out[i]) << " in = " <<
					static_cast<size_t>(in[i]) << " gold = " <<
					static_cast<size_t>((in[i] << shift) | prev) << std::endl;
		}

		prev = in[i] >> inv_shift;
	}

	return result;
}


}  // ns test
}  // ns bit_ops
}  // ns utils
}  //ns bliss






#endif /* BIT_TEST_COMMON_HPP_ */
