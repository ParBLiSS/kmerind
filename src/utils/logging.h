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

#ifndef BLISS_LOGGING_H
#define BLISS_LOGGING_H

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
#include "logger_config.hpp"

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
#define FATAL(msg)       PRINT_FATAL(msg)
#define FATALF(msg, ...) PRINT_FATALF(msg, ##__VA_ARGS__)
#else
#define FATAL(msg)       NOPRINT_FATAL(msg)
#define FATALF(msg, ...) NOPRINT_FATALF(msg, ##__VA_ARGS__)
#endif

/// ERROR severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_ERROR
#define ERROR(msg)       PRINT_ERROR(msg)
#define ERRORF(msg, ...) PRINT_ERRORF(msg, ##__VA_ARGS__)
#else
#define ERROR(msg)       NOPRINT_ERROR(msg)
#define ERRORF(msg, ...) NOPRINT_ERRORF(msg, ##__VA_ARGS__)
#endif

/// WARNING severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_WARNING
#define WARNING(msg)       PRINT_WARNING(msg)
#define WARNINGF(msg, ...) PRINT_WARNINGF(msg, ##__VA_ARGS__)
#else
#define WARNING(msg)       NOPRINT_WARNING(msg)
#define WARNINGF(msg, ...) NOPRINT_WARNINGF(msg, ##__VA_ARGS__)
#endif

/// INFO severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_INFO
#define INFO(msg)       PRINT_INFO(msg)
#define INFOF(msg, ...) PRINT_INFOF(msg, ##__VA_ARGS__)
#else
#define INFO(msg)       NOPRINT_INFO(msg)
#define INFOF(msg, ...) NOPRINT_INFOF(msg, ##__VA_ARGS__)
#endif

/// DEBUG severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_DEBUG
#define DEBUG(msg)       PRINT_DEBUG(msg)
#define DEBUGF(msg, ...) PRINT_DEBUGF(msg, ##__VA_ARGS__)
#else
#define DEBUG(msg)       NOPRINT_DEBUG(msg)
#define DEBUGF(msg, ...) NOPRINT_DEBUGF(msg, ##__VA_ARGS__)
#endif

/// TRACE severity logging functions
#if LOGGER_VERBOSITY >= BLISS_LOGGER_VERBOSITY_TRACE
#define TRACE(msg)       PRINT_TRACE(msg)
#define TRACEF(msg, ...) PRINT_TRACEF(msg, ##__VA_ARGS__)
#else
#define TRACE(msg)       NOPRINT_TRACE(msg)
#define TRACEF(msg, ...) NOPRINT_TRACEF(msg, ##__VA_ARGS__)
#endif

/// Initialization function for logger
// in case a logger doesn't define an init() function, supply an empty one
#ifndef LOG_INIT
#define LOG_INIT() ;
#endif

#endif // BLISS_LOGGING_H
