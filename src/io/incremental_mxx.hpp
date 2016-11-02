/*
 * Copyright 2016 Georgia Institute of Technology
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
 * @file    incremental mxx.hpp
 * @ingroup
 * @author  tpan
 * @brief   extends mxx to send/recv incrementally, so as to minimize memory use and allocation
 * @details
 *
 */


namespace imxx
{
// specialized for lower memory, incremental operations
// also  later, non-blocking mxx.
// requirements:  in place in original or fixed-size buffer(s)
// iterative processing

//== get bucket assignment. - grouping completely by processors, or grouping by processors in communication blocks
// parameter: in block per bucket send count (each block is for an all2all operation.  max block size is capped by max int.



//== bucket actually.

//== all2allv using in place all2all as first step, then all2allv + buffer as last part.  pure communication
// to use all2all, the data has to be contiguous.  2 ways to do this:  1. in place.  2. separate buffer
void incremental_all2allv() {
	//  compute send counts
	//  then do in place send/recv
}


//== transform before or after communication

//== communicate, process (one to one), return communication

//== communicate, process (not one to one), return communication


// mpi3 versions?

}
