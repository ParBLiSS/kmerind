/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    mxx_support.hpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */
#ifndef SRC_DEBRUIJN_MXX_SUPPORT_HPP_
#define SRC_DEBRUIJN_MXX_SUPPORT_HPP_

#include <mxx/datatypes.hpp>

//#include "utils/system_utils.hpp"

#include "debruijn/de_bruijn_node_trait.hpp"

namespace mxx {

  template<typename A, typename T>
    struct datatype_builder<::bliss::de_bruijn::node::edge_counts<A, T> > :
    public datatype_builder<decltype(::bliss::de_bruijn::node::edge_counts<A, T>::counts) > {

      typedef datatype_builder<decltype(::bliss::de_bruijn::node::edge_counts<A, T>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename A, typename T>
    struct datatype_builder<const bliss::de_bruijn::node::edge_counts<A, T> > : 
    public datatype_builder<decltype(bliss::de_bruijn::node::edge_counts<A, T>::counts) > {

      typedef datatype_builder<decltype(bliss::de_bruijn::node::edge_counts<A, T>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename A>
    struct datatype_builder<bliss::de_bruijn::node::edge_exists<A> > : 
    public datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > {

      typedef datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };


  template<typename A>
    struct datatype_builder<const bliss::de_bruijn::node::edge_exists<A> > : 
    public datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > {

      typedef datatype_builder<decltype(bliss::de_bruijn::node::edge_exists<A>::counts) > baseType;

      static MPI_Datatype get_type(){ 
        return baseType::get_type(); 
      }

      static size_t num_basic_elements() {
        return baseType::num_basic_elements();
      }
    };
}  // namespace mxx




#endif /* SRC_IO_MXX_SUPPORT_HPP_ */
