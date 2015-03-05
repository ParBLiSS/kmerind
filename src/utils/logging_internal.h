/**
 * @file    logging-internal.h
 * @ingroup utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 *
 * This file defines internal macros for logging of events of different `severity`.
 * This is meant to be use for code organization.
 *
 * Printing and non printing versions of logging macros are defined here
 * so that we can dynamically define logger macros at different VERBOSITY levels.
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add Licence
 */

#ifndef BLISS_LOGGING_INTERNAL_H
#define BLISS_LOGGING_INTERNAL_H


// MACRO DEFINITIONS  HERE REQUIRE THE "do {...} while (false)" construct for them to play nicely with if ... else ... statements.
// see http://stackoverflow.com/questions/154136/do-while-and-if-else-statements-in-c-c-macros for explanation.
// essentially, this avoids
//  if X
//     Z;    #define Z {a(); b()}
//  else
//     c();
// the construct allows the macros to be used as if they are functions.
// while (false) will be optimized away by compiler.
// notice the macro definition DOES NOT HAVE A SEMICOLON AT END.


/// non-printing version of functions
#define NOPRINT_FATAL(msg) do { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); } while (false)
#define NOPRINT_ERROR(msg) do {} while (false)
#define NOPRINT_WARNING(msg) do {} while (false)
#define NOPRINT_INFO(msg) do {} while (false)
#define NOPRINT_DEBUG(msg) do {} while (false)
#define NOPRINT_TRACE(msg) do {} while (false)

#define NOPRINT_FATALF(msg, ...) do { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); } while (false)
#define NOPRINT_ERRORF(msg, ...) do {} while (false)
#define NOPRINT_WARNINGF(msg, ...) do {} while (false)
#define NOPRINT_INFOF(msg, ...) do {} while (false)
#define NOPRINT_DEBUGF(msg, ...) do {} while (false)
#define NOPRINT_TRACEF(msg, ...) do {} while (false)



/*********************************************************************
 *                          Disable logging                          *
 *********************************************************************/

#if USE_LOGGER == BLISS_LOGGING_NO_LOG

// empty logging macros (no overhead)
#define PRINT_FATAL(msg) do { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); } while (false)
#define PRINT_ERROR(msg) do {} while (false)
#define PRINT_WARNING(msg) do {} while (false)
#define PRINT_INFO(msg) do {} while (false)
#define PRINT_DEBUG(msg) do {} while (false)
#define PRINT_TRACE(msg) do {} while (false)



/*********************************************************************
 *                     use std::cerr for logging                     *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_CERR

// simple output via std::cerr
// NOTE: this is not thread safe

#include <iostream>

#define PRINT_FATAL(msg)      do { std::cerr << "[fatal] " << msg << std::endl << std::flush; exit(-1); } while (false)
#define PRINT_ERROR(msg)      do { std::cerr << "[error] " << msg << std::endl << std::flush; } while (false)
#define PRINT_WARNING(msg)    do { std::cerr << "[warn ] " << msg << std::endl << std::flush; } while (false)
#define PRINT_INFO(msg)       do { std::cerr << "[info ] " << msg << std::endl << std::flush; } while (false)
#define PRINT_DEBUG(msg)      do { std::cerr << "[debug] " << msg << std::endl << std::flush; } while (false)
#define PRINT_TRACE(msg)      do { std::cerr << "[trace] " << msg << std::endl << std::flush; } while (false)


/*********************************************************************
 *                     use printf for logging                     *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_PRINTF

// simple output via std::cerr
// NOTE: this is not thread safe
#include <sstream>

#define PRINT_FATAL(msg)    do { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); } while (false)
#define PRINT_ERROR(msg)    do { std::stringstream ss; ss << msg; printf("[error] %s\n", ss.str().c_str());  } while (false)
#define PRINT_WARNING(msg)  do { std::stringstream ss; ss << msg; printf("[warn ] %s\n", ss.str().c_str());  } while (false)
#define PRINT_INFO(msg)     do { std::stringstream ss; ss << msg; printf("[info ] %s\n", ss.str().c_str());  } while (false)
#define PRINT_DEBUG(msg)    do { std::stringstream ss; ss << msg; printf("[debug] %s\n", ss.str().c_str());  } while (false)
#define PRINT_TRACE(msg)    do { std::stringstream ss; ss << msg; printf("[trace] %s\n", ss.str().c_str());  } while (false)




/*********************************************************************
 *                      use boost::log::trivial                      *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_BOOST_TRIVIAL

// using boost trival logging
#include <boost/log/trivial.hpp>

#define PRINT_FATAL(msg)      do { BOOST_LOG_TRIVIAL(fatal) << msg;  exit(-1);  } while (false)
#define PRINT_ERROR(msg)      do { BOOST_LOG_TRIVIAL(error) << msg;   } while (false)
#define PRINT_WARNING(msg)    do { BOOST_LOG_TRIVIAL(warning) << msg; } while (false)
#define PRINT_INFO(msg)       do { BOOST_LOG_TRIVIAL(info) << msg;    } while (false)
#define PRINT_DEBUG(msg)      do { BOOST_LOG_TRIVIAL(debug) << msg;   } while (false)
#define PRINT_TRACE(msg)      do { BOOST_LOG_TRIVIAL(trace) << msg;   } while (false)


/*********************************************************************
 *             use a customized boost::log based logger              *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_BOOST_CUSTOM

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

  // namespace aliases for boost::log
  namespace logging = boost::log;
  namespace src = boost::log::sources;
  namespace expr = boost::log::expressions;
  namespace sinks = boost::log::sinks;
  namespace attrs = boost::log::attributes;
  namespace keywords = boost::log::keywords;



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
std::ostream& operator<< (std::ostream& strm, severity_level level)
{
    static const char* strings[] =
    {
        "trace",
        "debug",
        "info ",
        "warn ",
        "error",
        "fatal"
    };

    if (static_cast< std::size_t >(level) < sizeof(strings) / sizeof(*strings))
        strm << strings[level];
    else
        strm << static_cast< int >(level);

    return strm;
}


void log_init_text_sink()
{
    boost::shared_ptr< logging::core > core = logging::core::get();

    // Create a backend and attach a couple of streams to it
    boost::shared_ptr< sinks::text_ostream_backend > backend =
        boost::make_shared< sinks::text_ostream_backend >();
    backend->add_stream(
        boost::shared_ptr< std::ostream >(&std::clog, logging::empty_deleter()));
    //backend->add_stream(
    //  boost::shared_ptr< std::ostream >(new std::ofstream("sample.log")));

    // Enable auto-flushing after each log record written
    backend->auto_flush(true);

    // Wrap it into the frontend and register in the core.
    // The backend requires synchronization in the frontend.
    typedef sinks::synchronous_sink< sinks::text_ostream_backend > sink_t;
    boost::shared_ptr< sink_t > sink(new sink_t(backend));

    // Setup the common formatter for all sinks
    boost::log::formatter fmt = boost::log::expressions::stream
        << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S.%f")
        << ": <" << severity << ">\t"
        // the message includes the file, line and function information
        << boost::log::expressions::smessage;
    sink->set_formatter(fmt);

    core->add_sink(sink);
}


// init boost logging
void init()
{
    log_init_text_sink();
    logging::add_common_attributes();
    printf("USING BOOST_LOG_CUSTOM_LOGGER\n");
}


// macro to add file and line
#define _LOG_MSG(msg) __FILE__ << ":" << __LINE__ << ":\t" << msg ;


#define PRINT_FATAL(msg)      do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::fatal) << _LOG_MSG(msg); exit(-1);   } while (false)
#define PRINT_ERROR(msg)      do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::error) << _LOG_MSG(msg);   } while (false)
#define PRINT_WARNING(msg)    do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::warning) << _LOG_MSG(msg); } while (false)
#define PRINT_INFO(msg)       do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::info) << _LOG_MSG(msg);    } while (false)
#define PRINT_DEBUG(msg)      do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::debug) << _LOG_MSG(msg);   } while (false)
#define PRINT_TRACE(msg)      do { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::trace) << _LOG_MSG(msg);   } while (false)

#define LOG_INIT() bliss::log::init()



} // namespace log

} // namespace bliss

#else
#error "Invalid value of USE_LOGGER: logger not found."
#endif







/*******************************
 *  printf style debug macros  *
 *******************************/
#if USE_LOGGER == BLISS_LOGGING_NO_LOG

// empty logging macros (no overhead)
#define PRINT_FATALF(msg, ...) do { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); } while (false)
#define PRINT_ERRORF(msg, ...) do {} while (false)
#define PRINT_WARNINGF(msg, ...) do {} while (false)
#define PRINT_INFOF(msg, ...) do {} while (false)
#define PRINT_DEBUGF(msg, ...) do {} while (false)
#define PRINT_TRACEF(msg, ...) do {} while (false)

#elif USE_LOGGER == BLISS_LOGGING_PRINTF

// C automatically concatenate c string literals, so "blah "msg"\n" where msg is "BLAH", becomes "blah BLAH\n"
// C++11 requires space in between:  "blah " msg "\n".
#define PRINT_FATALF(msg, ...)   do { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); } while (false)
#define PRINT_ERRORF(msg, ...)   do { printf("[error] " msg "\n", ##__VA_ARGS__); } while (false)
#define PRINT_WARNINGF(msg, ...) do { printf("[warn ] " msg "\n", ##__VA_ARGS__); } while (false)
#define PRINT_INFOF(msg, ...)    do { printf("[info ] " msg "\n", ##__VA_ARGS__); } while (false)
#define PRINT_DEBUGF(msg, ...)   do { printf("[debug] " msg "\n", ##__VA_ARGS__); } while (false)
#define PRINT_TRACEF(msg, ...)   do { printf("[trace] " msg "\n", ##__VA_ARGS__); } while (false)


#else
#define BLISS_SPRINTF_BUFFER_SIZE 256

#define PRINT_FATALF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_FATAL(buffer);\
        exit(-1); \
        } while (false)

#define PRINT_ERRORF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_ERROR(buffer);\
        } while (false)

#define PRINT_WARNINGF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_WARNING(buffer);\
        } while (false)

#define PRINT_INFOF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_INFO(buffer);\
        } while (false)

#define PRINT_DEBUGF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_DEBUG(buffer);\
        } while (false);

#define PRINT_TRACEF(msg, ...) do {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_TRACE(buffer);\
        } while (false)
#endif




// in case a logger doesn't define an init() function, supply an empty one
#ifndef LOG_INIT
#define LOG_INIT() ;
#endif

#endif // BLISS_LOGGING_INTERNAL_H
