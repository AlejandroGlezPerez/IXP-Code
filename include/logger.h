//  logger.h
//  IXP Code

#ifndef LOGGER_H_
#define LOGGER_H_

#include <cstddef>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/trivial.hpp>

BOOST_LOG_GLOBAL_LOGGER(logger, boost::log::sources::severity_logger_mt<
    boost::log::trivial::severity_level>)
#define LOG(SEVERITY) BOOST_LOG_SEV(logger::get(), \
                                    boost::log::trivial::SEVERITY)

#endif // LOGGER_H_
