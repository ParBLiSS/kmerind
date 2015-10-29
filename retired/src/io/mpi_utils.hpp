/**
 * @file    mpi_utils.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef MPI_UTILS_HPP_
#define MPI_UTILS_HPP_




namespace bliss
{
  namespace io
  {

    typedef int64_t TaggedEpoch;

    /// convenience function to extract tag from TaggedEpoch
    inline int getTagFromTaggedEpoch(const TaggedEpoch & te) {
      return (reinterpret_cast<const int*>(&te))[1];
    }

    /// convenience function to extract epoch from TaggedEpoch
    inline uint64_t getEpochFromTaggedEpoch(const TaggedEpoch & te) {
      return te & 0x00000000FFFFFFFF;
    }

    /// convenience function to extract tag from TaggedEpoch
    inline TaggedEpoch makeTaggedEpoch(const int& tag, const uint64_t& epoch) {
      TaggedEpoch te = static_cast<const TaggedEpoch>(epoch);
      (reinterpret_cast<int *>(&te))[1] = tag;
      return te;
    }

    inline TaggedEpoch updateEpoch(const TaggedEpoch& te, const uint64_t& epoch) {
      return (te & 0xFFFFFFFF00000000) | (epoch & 0x00000000FFFFFFFF);
    }


    /// CONTOL TAG's tagged_epoch
    constexpr TaggedEpoch CONTROL_TAGGED_EPOCH = 0x00000000FFFFFFFF;

    /// CONTOL TAG's epoch
    constexpr uint64_t CONTROL_TAG_EPOCH = 0x00000000FFFFFFFF;

    /// CONTOL TAG
    constexpr int CONTROL_TAG = 0;
  }
}

#endif /* MPI_UTILS_HPP_ */
