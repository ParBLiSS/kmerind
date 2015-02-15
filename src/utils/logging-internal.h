/**
 * @file    logging-internal.h
 * @ingroup bliss::utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick
 * @author  Tony Pan
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



/// non-printing version of functions
#define NOPRINT_FATAL(msg) { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); }
#define NOPRINT_ERROR(msg)
#define NOPRINT_WARNING(msg)
#define NOPRINT_INFO(msg)
#define NOPRINT_DEBUG(msg)
#define NOPRINT_TRACE(msg)

#define NOPRINT_FATALF(msg, ...) { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); }
#define NOPRINT_ERRORF(msg, ...)
#define NOPRINT_WARNINGF(msg, ...)
#define NOPRINT_INFOF(msg, ...)
#define NOPRINT_DEBUGF(msg, ...)
#define NOPRINT_TRACEF(msg, ...)



/*********************************************************************
 *                          Disable logging                          *
 *********************************************************************/

#if USE_LOGGER == BLISS_LOGGING_NO_LOG

// empty logging macros (no overhead)
#define PRINT_FATAL(msg) { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); }
#define PRINT_ERROR(msg)
#define PRINT_WARNING(msg)
#define PRINT_INFO(msg)
#define PRINT_DEBUG(msg)
#define PRINT_TRACE(msg)



/*********************************************************************
 *                     use std::cerr for logging                     *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_CERR

// simple output via std::cerr
// NOTE: this is not thread safe

#include <iostream>

#define PRINT_FATAL(msg)      { std::cerr << "[fatal] " << msg << std::endl << std::flush; exit(-1); }
#define PRINT_ERROR(msg)      { std::cerr << "[error] " << msg << std::endl << std::flush; }
#define PRINT_WARNING(msg)    { std::cerr << "[warn ] " << msg << std::endl << std::flush; }
#define PRINT_INFO(msg)       { std::cerr << "[info ] " << msg << std::endl << std::flush; }
#define PRINT_DEBUG(msg)      { std::cerr << "[debug] " << msg << std::endl << std::flush; }
#define PRINT_TRACE(msg)      { std::cerr << "[trace] " << msg << std::endl << std::flush; }


/*********************************************************************
 *                     use printf for logging                     *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_PRINTF

// simple output via std::cerr
// NOTE: this is not thread safe
#include <sstream>

#define PRINT_FATAL(msg)    { std::stringstream ss; ss << msg; printf("[fatal] %s\n", ss.str().c_str());  exit(-1); }
#define PRINT_ERROR(msg)    { std::stringstream ss; ss << msg; printf("[error] %s\n", ss.str().c_str());  }
#define PRINT_WARNING(msg)  { std::stringstream ss; ss << msg; printf("[warn ] %s\n", ss.str().c_str());  }
#define PRINT_INFO(msg)     { std::stringstream ss; ss << msg; printf("[info ] %s\n", ss.str().c_str());  }
#define PRINT_DEBUG(msg)    { std::stringstream ss; ss << msg; printf("[debug] %s\n", ss.str().c_str());  }
#define PRINT_TRACE(msg)    { std::stringstream ss; ss << msg; printf("[trace] %s\n", ss.str().c_str());  }




/*********************************************************************
 *                      use boost::log::trivial                      *
 *********************************************************************/

#elif USE_LOGGER == BLISS_LOGGING_BOOST_TRIVIAL

// using boost trival logging
#include <boost/log/trivial.hpp>

#define PRINT_FATAL(msg)      { BOOST_LOG_TRIVIAL(fatal) << msg;  exit(-1);  }
#define PRINT_ERROR(msg)      { BOOST_LOG_TRIVIAL(error) << msg;   }
#define PRINT_WARNING(msg)    { BOOST_LOG_TRIVIAL(warning) << msg; }
#define PRINT_INFO(msg)       { BOOST_LOG_TRIVIAL(info) << msg;    }
#define PRINT_DEBUG(msg)      { BOOST_LOG_TRIVIAL(debug) << msg;   }
#define PRINT_TRACE(msg)      { BOOST_LOG_TRIVIAL(trace) << msg;   }


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


#define PRINT_FATAL(msg)      { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::fatal) << _LOG_MSG(msg); exit(-1);   }
#define PRINT_ERROR(msg)      { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::error) << _LOG_MSG(msg);   }
#define PRINT_WARNING(msg)    { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::warning) << _LOG_MSG(msg); }
#define PRINT_INFO(msg)       { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::info) << _LOG_MSG(msg);    }
#define PRINT_DEBUG(msg)      { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::debug) << _LOG_MSG(msg);   }
#define PRINT_TRACE(msg)      { BOOST_LOG_SEV(bliss::log::global_logger, bliss::log::trace) << _LOG_MSG(msg);   }

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
#define PRINT_FATALF(msg, ...) { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); }
#define PRINT_ERRORF(msg, ...)
#define PRINT_WARNINGF(msg, ...)
#define PRINT_INFOF(msg, ...)
#define PRINT_DEBUGF(msg, ...)
#define PRINT_TRACEF(msg, ...)

#elif USE_LOGGER == BLISS_LOGGING_PRINTF

// C automatically concatenate c string literals, so "blah "msg"\n" where msg is "BLAH", becomes "blah BLAH\n"
// C++11 requires space in between:  "blah " msg "\n".
#define PRINT_FATALF(msg, ...)   { printf("[fatal] " msg "\n", ##__VA_ARGS__); exit(-1); }
#define PRINT_ERRORF(msg, ...)   { printf("[error] " msg "\n", ##__VA_ARGS__); }
#define PRINT_WARNINGF(msg, ...) { printf("[warn ] " msg "\n", ##__VA_ARGS__); }
#define PRINT_INFOF(msg, ...)    { printf("[info ] " msg "\n", ##__VA_ARGS__); }
#define PRINT_DEBUGF(msg, ...)   { printf("[debug] " msg "\n", ##__VA_ARGS__); }
#define PRINT_TRACEF(msg, ...)   { printf("[trace] " msg "\n", ##__VA_ARGS__); }


#else
#define BLISS_SPRINTF_BUFFER_SIZE 256

#define PRINT_FATALF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_FATAL(buffer);\
        exit(-1); \
        }

#define PRINT_ERRORF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_ERROR(buffer);\
        }

#define PRINT_WARNINGF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_WARNING(buffer);\
        }

#define PRINT_INFOF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_INFO(buffer);\
        }

#define PRINT_DEBUGF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_DEBUG(buffer);\
        }

#define PRINT_TRACEF(msg, ...) {\
        char buffer[BLISS_SPRINTF_BUFFER_SIZE]; \
        snprintf(buffer, BLISS_SPRINTF_BUFFER_SIZE, msg, ##__VA_ARGS__);\
        PRINT_TRACE(buffer);\
        }
#endif




// in case a logger doesn't define an init() function, supply an empty one
#ifndef LOG_INIT
#define LOG_INIT() ;
#endif

#endif // BLISS_LOGGING_INTERNAL_H
