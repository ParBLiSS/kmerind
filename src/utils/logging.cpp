/**
 * @file    logging.h
 * @ingroup utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 *
 * This file defines macros for logging of events of different `severity`.
 *
 * There are 2 sets of macros APIs:
 * 	the first set uses c++ style stream operator <<
 *      - FATAL(msg)
 *      - ERROR(msg)
 *      - WARNING(msg)
 *      - INFO(msg)
 *      - DEBUG(msg)
 *      - TRACE(msg)
 *     Example how the macro has to be called:
 *          DEBUG("current coordinates: x = " << x << ", y = " << y);
 *
 *     Any type that implements the << operator can be printed to the log in this
 *     way.
 * 	the second set uses c style printf syntax.  These are suffixed with F.
 *      - FATALF(format, ...)
 *      - ERRORF(format, ...)
 *      - WARNINGF(format, ...)
 *      - INFOF(format, ...)
 *      - DEBUGF(format, ...)
 *      - TRACEF(format, ...)
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
 * minimum required severity level, i.e. FATAL, ERROR, WARNING, etc.
 * 	when a verbosity level is specified, messages at the corresponding
 * 	severity level, and any level more severe, are printed, while messages
 * 	from less severe levels are discards.
 * e.g.
 * 		verbosity of INFO included FATAL, ERROR, WARNING, and INFO,
 * 			but not DEBUG or TRACE.
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
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add Licence
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


