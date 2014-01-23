// potentially packed (or not packed) short sequence of alphabet.
// this is primarily for managing packed alphabets
// the word may have padded 0.

// for example, we can fit 5 3 bit letters into 2 bytes, with 1 bit padding, or 8 3 bit letters into 3 bytes.

// template specialization should favor ones that fit nicely into machine word sizes, e.g. 32, 64, or 128 bits.

