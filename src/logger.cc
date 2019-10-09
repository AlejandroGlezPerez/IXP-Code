#include "logger.h"

#include <fstream>
#include <ostream>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/date_time.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

BOOST_LOG_ATTRIBUTE_KEYWORD(line_number, "LineNumber", std::size_t)
BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "Timestamp", boost::posix_time::ptime)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity",
                            boost::log::trivial::severity_level)

BOOST_LOG_GLOBAL_LOGGER_INIT(logger, boost::log::sources::severity_logger_mt) {
    boost::log::sources::severity_logger_mt<
        boost::log::trivial::severity_level> logger;

    // Add logger attributes
    logger.add_attribute("LineNumber",
                         boost::log::attributes::counter<std::size_t>(1));
    logger.add_attribute("Timestamp", boost::log::attributes::local_clock());

    // Create text sink to write to standard output
    typedef boost::log::sinks::synchronous_sink<
        boost::log::sinks::text_ostream_backend> text_sink_t;
    boost::shared_ptr<text_sink_t> console_sink =
        boost::make_shared<text_sink_t>();
    console_sink->locked_backend()->add_stream(boost::shared_ptr<std::ostream>(
        &std::clog, boost::null_deleter()));
  
    // Create file sink to write to files
    typedef boost::log::sinks::synchronous_sink<
        boost::log::sinks::text_file_backend> file_sink_t;
    boost::shared_ptr<file_sink_t> file_sink = boost::make_shared<file_sink_t>(
        boost::log::keywords::file_name = "log/ixp_%Y.%m.%d_%H.%M.%S_%N.log",
        boost::log::keywords::target_file_name = "log/ixp_current.log",
        boost::log::keywords::enable_final_rotation = false);

    // Format log messages and hide logging info from main screen
    boost::log::formatter console_format = boost::log::expressions::stream
        << boost::log::expressions::smessage;
    console_sink->set_formatter(console_format);
  
    boost::log::formatter file_format = boost::log::expressions::stream
        << std::setw(7) << std::setfill('0') << line_number << " | "
        << boost::log::expressions::format_date_time(timestamp, "%H:%M:%S.%f")
        << " [" << boost::log::trivial::severity << "]"
        << " - " << boost::log::expressions::smessage;
    file_sink->set_formatter(file_format);

    // Filter
    console_sink->set_filter(severity >= boost::log::trivial::debug);
    file_sink->set_filter(severity >= boost::log::trivial::info);

    // Add sinks to logger
    boost::log::core::get()->add_sink(console_sink);
    boost::log::core::get()->add_sink(file_sink);

    return logger;
}
