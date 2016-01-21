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

// TODO time functions so that we don't have the unused return values.

#ifndef BLISS_LOGGING_H
#define BLISS_LOGGING_H

// cast to void to suppress unused vars when log level is set to minimal.
#define USED_BY_LOGGER_ONLY(x) do { (void)(x); } while (0)


/// LOG ENGINE TYPES
// Disable logging
#define BLISS_LOGGING_NO_LOG         1
// Use std::cerr output for logging (not thread safe)
#define BLISS_LOGGING_CERR           2
// Use boost::log::trivial for logging
#define BLISS_LOGGING_BOOST_TRIVIAL  3
// Use BLISS's customized boost logging wrapper for logging
#define BLISS_LOGGING_BOOST_CUSTOM   4
// using printf
#define BLISS_LOGGING_PRINTF         5

/// logger verbosity.  these are listed in increasing verbosity. each level include all before it.
#define BLISS_LOGGER_VERBOSITY_FATAL   0
#define BLISS_LOGGER_VERBOSITY_ERROR   1
#define BLISS_LOGGER_VERBOSITY_WARNING 2
#define BLISS_LOGGER_VERBOSITY_INFO    3
#define BLISS_LOGGER_VERBOSITY_DEBUG   4
#define BLISS_LOGGER_VERBOSITY_TRACE   5

/// include logger_config.hpp to get the cmake configured logger settings.
#include "bliss-logger_config.hpp"

// set default logger (in case none is specified via compiler flags)
#ifndef USE_LOGGER
// set the default logger to the std::cerr logger
#define USE_LOGGER BLISS_LOGGING_PRINTF
#endif

#ifndef LOGGER_VERBOSITY
// set the default logger to the std::cerr logger
#define LOGGER_VERBOSITY BLISS_LOGGER_VERBOSITY_WARNING
#endif

#include "utils/logging_internal.h"

/********************************************
 *  actual macros, conditioned on verbosity *
 ********************************************/

/// FATAL severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_FATAL
#define BL_FATAL(msg)       PRINT_FATAL(msg)
#define BL_FATALF(msg, ...) PRINT_FATALF(msg, ##__VA_ARGS__)
#else
#define BL_FATAL(msg)       NOPRINT_FATAL(msg)
#define BL_FATALF(msg, ...) NOPRINT_FATALF(msg, ##__VA_ARGS__)
#endif

/// ERROR severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_ERROR
#define BL_ERROR(msg)       PRINT_ERROR(msg)
#define BL_ERRORF(msg, ...) PRINT_ERRORF(msg, ##__VA_ARGS__)
#else
#define BL_ERROR(msg)       NOPRINT_ERROR(msg)
#define BL_ERRORF(msg, ...) NOPRINT_ERRORF(msg, ##__VA_ARGS__)
#endif

/// WARNING severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_WARNING
#define BL_WARNING(msg)       PRINT_WARNING(msg)
#define BL_WARNINGF(msg, ...) PRINT_WARNINGF(msg, ##__VA_ARGS__)
#else
#define BL_WARNING(msg)       NOPRINT_WARNING(msg)
#define BL_WARNINGF(msg, ...) NOPRINT_WARNINGF(msg, ##__VA_ARGS__)
#endif

/// INFO severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_INFO
#define BL_INFO(msg)       PRINT_INFO(msg)
#define BL_INFOF(msg, ...) PRINT_INFOF(msg, ##__VA_ARGS__)
#else
#define BL_INFO(msg)       NOPRINT_INFO(msg)
#define BL_INFOF(msg, ...) NOPRINT_INFOF(msg, ##__VA_ARGS__)
#endif

/// DEBUG severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_DEBUG
#define BL_DEBUG(msg)       PRINT_DEBUG(msg)
#define BL_DEBUGF(msg, ...) PRINT_DEBUGF(msg, ##__VA_ARGS__)
#else
#define BL_DEBUG(msg)       NOPRINT_DEBUG(msg)
#define BL_DEBUGF(msg, ...) NOPRINT_DEBUGF(msg, ##__VA_ARGS__)
#endif

/// TRACE severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_TRACE
#define BL_TRACE(msg)       PRINT_TRACE(msg)
#define BL_TRACEF(msg, ...) PRINT_TRACEF(msg, ##__VA_ARGS__)
#else
#define BL_TRACE(msg)       NOPRINT_TRACE(msg)
#define BL_TRACEF(msg, ...) NOPRINT_TRACEF(msg, ##__VA_ARGS__)
#endif


/// Initialization function for logger
// in case a logger doesn't define an init() function, supply an empty one
#ifndef LOG_INIT
#define LOG_INIT() ;
#endif

#endif // BLISS_LOGGING_H
