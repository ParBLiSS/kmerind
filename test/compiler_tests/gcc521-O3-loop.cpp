/**
 * @file    gcc521-O3-loop.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2016 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <string>
#include <iostream>


/// static constant for end of line.  note this is same for unicode as well
static constexpr unsigned char eol = '\n';
/// static constant for carriage return.  note this is same for unicode as well
static constexpr unsigned char cr = '\r';

 /**
  * @brief  search for first non-EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
  * @details       iter can point to the previous EOL, or a nonEOL character.
  * @param[in/out] iter    iterator to advance
  * @param[in]     end     position to stop the traversal
  * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
  * @return        iterator at the new position, where the Non EOL char is found, or end.
  */
 template <typename IT,
     typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                          ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                         >::type >
 inline IT findNonEOL(IT& iter, const IT& end, size_t &offset) {
   while ((iter != end) && ((*iter == eol) || (*iter == cr))) {
     ++iter;
     ++offset;
   }
   return iter;
 }

 /**
  * @brief  search for first EOL character in a iterator, returns the stopping position as an iterator.  also records offset.
  * @details       iter can point to a nonEOL char, or an EOL char.
  * @param[in/out] iter    iterator to advance
  * @param[in]     end     position to stop the traversal
  * @param[in/out] offset  the global offset of the sequence record within the file.  used as record id.
  * @return        iterator at the new position, where the EOL char is found, or end.
  */
 template <typename IT,
     typename = typename ::std::enable_if<(::std::is_same<typename ::std::iterator_traits<IT>::value_type, char>::value ||
                                          ::std::is_same<typename ::std::iterator_traits<IT>::value_type, unsigned char>::value)
                                         >::type >
 inline IT findEOL(IT& iter, const IT& end, size_t &offset) {
   while ((iter != end) && ((*iter != eol) && (*iter != cr) ) ) {
     ++iter;
     ++offset;
   }
   return iter;
 }


int main(int argc, char** argv) {

  std::string input("abcdefg\n"
                    "hijklmn\n"
                    "opqrstu\n"
                    "vwxyz12\n");

  size_t len = input.length();

  unsigned char first[4] = {0, 0, 0, 0};
  size_t offsets[4] = {len, len, len, len};

  auto iter = input.begin();
  auto end = input.end();

  size_t i = 0;


  // now find the first non eol
  iter = findNonEOL(iter, end, i);
  if (i == len) return len;
  // assign to first.
  first[0] = *iter;
  offsets[0] = i;

//          std::cout << "1 first chars are " << first[0] << "," << first[1] << "," << first[2] << "," << first[3] << std::endl;


  // lines 2 through 4
  for (int j = 1; j < 4; ++j) {
    iter = findEOL(iter, end, i);
    if (i == len) return len;
    iter = findNonEOL(iter, end, i);
    if (i == len) return len;
    first[j] = *iter;
    offsets[j] = i;
  }

  std::cout << " first chars are " << first[0] << "," << first[1] << "," << first[2] << "," << first[3] << std::endl;
  std::cout << " offsets are " << offsets[0] << "," << offsets[1] << "," << offsets[2] << "," << offsets[3] << std::endl;



}
