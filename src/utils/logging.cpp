/**
 * @file    logging.h
 * @ingroup utils
 * @brief   Defines macros for logging.
 * @author  Patrick Flick, Tony Pan
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

#include "utils/logging.h"


/*********************************************************************
 *             use a customized boost::log based logger              *
 *********************************************************************/

#if USE_LOGGER == BLISS_LOGGING_BOOST_CUSTOM

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

//BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
src::severity_logger<severity_level> global_logger;


// The operator puts a human-friendly representation of the severity level to
// the stream
std::ostream& operator<< (std::ostream& strm, severity_level level)
{
    static const char* strings[] =
    {
        "trace",
        "debug",
        "info",
        "warning",
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


} // namespace log

} // namespace bliss



#endif


