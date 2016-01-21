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
 * @file    logging.cpp
 * @ingroup utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 *
 * This file defines macros for logging of events of different `severity`.
 *
 * There are 2 sets of macros APIs:
 * 	the first set uses c++ style stream operator <<
 *      - BL_FATAL(msg)
 *      - BL_ERROR(msg)
 *      - BL_WARNING(msg)
 *      - BL_INFO(msg)
 *      - BL_DEBUG(msg)
 *      - BL_TRACE(msg)
 *     Example how the macro has to be called:
 *          DEBUG("current coordinates: x = " << x << ", y = " << y);
 *
 *     Any type that implements the << operator can be printed to the log in this
 *     way.
 * 	the second set uses c style printf syntax.  These are suffixed with F.
 *      - BL_FATALF(format, ...)
 *      - BL_ERRORF(format, ...)
 *      - BL_WARNINGF(format, ...)
 *      - BL_INFOF(format, ...)
 *      - BL_DEBUGF(format, ...)
 *      - BL_TRACEF(format, ...)
 *     Example how the macro has to be called:
 *          DEBUGF("current coordinates: x = %d, y = %d", x, y);
 *
 *     Any type that can be printed by printf can be used in this way.
 *
 * 	note in both cases, end of line char is not needed.
 *
 * 	Benefit of the c++ style macros is its convenience for printing complex types
 * 	Benefit of the c style macros is that they are more thread friendly:
 * 		messages from different threads are less likely to interleave
 * 		and we do not need to declare std::cout or other ostream objects as shared
 * 			in openmp clauses.
 *
 * Logging can be configured with different verbosity levels named after the
 * minimum required severity level, i.e. BL_FATAL, BL_ERROR, BL_WARNING, etc.
 * 	when a verbosity level is specified, messages at the corresponding
 * 	severity level, and any level more severe, are printed, while messages
 * 	from less severe levels are discards.
 * e.g.
 * 		verbosity of INFO included BL_FATAL, BL_ERROR, BL_WARNING, and BL_INFO,
 * 			but not BL_DEBUG or BL_TRACE.
 *
 * Implementation of this selection is via empty preprocessor macros.
 *
 *
 * Prior to using the logging functions, the LOG_INIT(); macro MUST be called.
 *
 * In order to compile the boost log versions, use the following compiler
 * flags:
 *      -DBOOST_LOG_DYN_LINK
 * and the linker flags:
 *      -lboost_log -lboost_thread -lboost_system -lpthread
 *

 */

#include "utils/logging.h"

/*********************************************************************
 *             use a customized boost::log based logger              *
 *********************************************************************/

#if USE_LOGGER == BLISS_LOGGING_BOOST_CUSTOM

namespace bliss
{

namespace log
{

  namespace src = boost::log::sources;

	//BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
	src::severity_logger<severity_level> global_logger;


} // namespace log

} // namespace bliss



#endif


