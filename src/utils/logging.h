/**
 * @file    logging.h
 * @ingroup utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick
 *
 * This file defines macros for logging of events of different `severity`.
 * Each macro takes any number of arguments concatenated via the stream
 * operator <<. NOTE: passing multiple arguments separated by comma (,) DOES
 * NOT WORK.
 * Example how the macro has to be called:
 *      DEBUG("current coordinates: x = " << x << ", y = " << y);
 *
 * Any type that implements the << operator can be printed to the log in this
 * way.
 *
 * The available macros are:
 *      - FATAL(msg)
 *      - ERROR(msg)
 *      - WARNING(msg)
 *      - INFO(msg)
 *      - DEBUG(msg)
 *      - TRACE(msg)
 *
 * Prior to using the logging functions, the LOG_INIT(); macro MUST be called.
 *
 * In order to compile the boost log versions, use the following compiler
 * flags:
 *      -DBOOST_LOG_DYN_LINK
 * and the linker flags:
 *      -lboost_log -lboost_thread -lboost_system -lpthread
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */

#ifndef BLISS_LOGGING_H
#define BLISS_LOGGING_H

/// Disable logging
#define BLISS_LOGGING_NO_LOG         1
/// Use std::cerr output for logging (not thread safe)
#define BLISS_LOGGING_CERR           2
/// Use boost::log::trivial for logging
#define BLISS_LOGGING_BOOST_TRIVIAL  3
/// Use BLISS's customized boost logging wrapper for logging
#define BLISS_LOGGING_BOOST_CUSTOM   4

#include "logger_config.hpp"

// set default logger (in case none is specified via compiler flags)
#ifndef USE_LOGGER
// set the default logger to the std::cerr logger
#define USE_LOGGER BLISS_LOGGING_CERR
#endif


/*********************************************************************
 *                          Disable logging                          *
 *********************************************************************/

#if USE_LOGGER == BLISS_LOGGING_NO_LOG

// empty logging macros (no overhead)
#define FATAL(MSG)
#define ERROR(MSG)
#define WARNING(MSG)
#define INFO(MSG)
#define DEBUG(MSG)
#define TRACE(MSG)



/*********************************************************************
 *                     use std::cerr for logging                     *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_CERR

// simple output via std::cerr
// NOTE: this is not thread safe

#include <iostream>

#define FATAL(MSG)      std::cerr << "[fatal]    " << MSG << std::endl;
#define ERROR(MSG)      std::cerr << "[error]    " << MSG << std::endl;
#define WARNING(MSG)    std::cerr << "[warning]  " << MSG << std::endl;
#define INFO(MSG)       std::cerr << "[info]     " << MSG << std::endl;
#define DEBUG(MSG)      std::cerr << "[debug]    " << MSG << std::endl;
#define TRACE(MSG)      std::cerr << "[trace]    " << MSG << std::endl;



/*********************************************************************
 *                      use boost::log::trivial                      *
 *********************************************************************/

#elif (USE_LOGGER == BLISS_LOGGING_BOOST_TRIVIAL) && defined(USE_BOOST_LOG)

// using boost trival logging
#include <boost/log/trivial.hpp>

#define FATAL(MSG)      BOOST_LOG_TRIVIAL(fatal) << MSG;
#define ERROR(MSG)      BOOST_LOG_TRIVIAL(error) << MSG;
#define WARNING(MSG)    BOOST_LOG_TRIVIAL(warning) << MSG;
#define INFO(MSG)       BOOST_LOG_TRIVIAL(info) << MSG;
#define DEBUG(MSG)      BOOST_LOG_TRIVIAL(debug) << MSG;
#define TRACE(MSG)      BOOST_LOG_TRIVIAL(trace) << MSG;


/*********************************************************************
 *             use a customized boost::log based logger              *
 *********************************************************************/

#elif (USE_LOGGER == BLISS_LOGGING_BOOST_CUSTOM) && defined(USE_BOOST_LOG)

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

namespace bliss
{

namespace log
{

  namespace src = boost::log::sources;


// We define our own severity levels
enum severity_level
{
    trace,
    debug,
    info,
    warning,
    error,
    fatal
};

BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
extern src::severity_logger<severity_level> global_logger;


// The operator puts a human-friendly representation of the severity level to
// the stream
std::ostream& operator<< (std::ostream& strm, severity_level level);


void log_init_text_sink();


// init boost logging
void init();


// macro to add file and line
//#define _LOG_MSG(MSG) MSG << "\t(" << __FILE__ << ":" << __LINE__ << ")"
#define _LOG_MSG(MSG) __FILE__ << ":" << __LINE__ << ":\t" << MSG


#define FATAL(MSG)      BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::fatal) << _LOG_MSG(MSG);
#define ERROR(MSG)      BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::error) << _LOG_MSG(MSG);
#define WARNING(MSG)    BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::warning) << _LOG_MSG(MSG);
#define INFO(MSG)       BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::info) << _LOG_MSG(MSG);
#define DEBUG(MSG)      BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::debug) << _LOG_MSG(MSG);
#define TRACE(MSG)      BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::trace) << _LOG_MSG(MSG);

#define LOG_INIT() bliss::log::init()

} // namespace log

} // namespace bliss



#else
#error "Invalid value of USE_LOGGER: logger not found."
#endif

// in case a logger doesn't define an init() function, supply an empty one
#ifndef LOG_INIT
#define LOG_INIT() ;
#endif

#endif // BLISS_LOGGING_H
