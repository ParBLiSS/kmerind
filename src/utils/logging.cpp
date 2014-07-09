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

  namespace src = boost::log::sources;

//BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
src::severity_logger<severity_level> global_logger;


} // namespace log

} // namespace bliss



#endif


