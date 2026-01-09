/**
 * @file Logger.h
 * @brief Logging infrastructure for OpenMagnetics library
 * 
 * This file provides a flexible logging system with configurable
 * verbosity levels, output destinations, and module-based filtering.
 */

#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <mutex>
#include <chrono>
#include <iomanip>
#include <functional>
#include <memory>
#include <vector>

namespace OpenMagnetics {

/**
 * @brief Log severity levels
 */
enum class LogLevel : uint8_t {
    TRACE = 0,    ///< Detailed trace information
    DEBUG = 1,    ///< Debug information
    INFO = 2,     ///< General information
    WARNING = 3,  ///< Warning messages
    ERROR = 4,    ///< Error messages
    CRITICAL = 5, ///< Critical errors
    OFF = 6       ///< Disable all logging
};

/**
 * @brief Convert log level to string
 */
inline std::string to_string(LogLevel level) {
    switch (level) {
        case LogLevel::TRACE: return "TRACE";
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::INFO: return "INFO";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::CRITICAL: return "CRITICAL";
        case LogLevel::OFF: return "OFF";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Log output sink interface
 */
class LogSink {
public:
    virtual ~LogSink() = default;
    virtual void write(LogLevel level, const std::string& module, 
                       const std::string& message, 
                       const std::string& timestamp) = 0;
    virtual void flush() = 0;
};

/**
 * @brief Console log sink (stdout/stderr)
 */
class ConsoleSink : public LogSink {
public:
    explicit ConsoleSink(bool useColors = true) : _useColors(useColors) {}
    
    void write(LogLevel level, const std::string& module, 
               const std::string& message, 
               const std::string& timestamp) override {
        std::ostream& out = (level >= LogLevel::ERROR) ? std::cerr : std::cout;
        
        if (_useColors) {
            out << getColorCode(level);
        }
        
        out << "[" << timestamp << "] "
            << "[" << to_string(level) << "] ";
        
        if (!module.empty()) {
            out << "[" << module << "] ";
        }
        
        out << message;
        
        if (_useColors) {
            out << "\033[0m"; // Reset color
        }
        
        out << std::endl;
    }
    
    void flush() override {
        std::cout.flush();
        std::cerr.flush();
    }
    
private:
    bool _useColors;
    
    static std::string getColorCode(LogLevel level) {
        switch (level) {
            case LogLevel::TRACE: return "\033[90m";    // Gray
            case LogLevel::DEBUG: return "\033[36m";    // Cyan
            case LogLevel::INFO: return "\033[32m";     // Green
            case LogLevel::WARNING: return "\033[33m";  // Yellow
            case LogLevel::ERROR: return "\033[31m";    // Red
            case LogLevel::CRITICAL: return "\033[35m"; // Magenta
            default: return "";
        }
    }
};

/**
 * @brief File log sink
 */
class FileSink : public LogSink {
public:
    explicit FileSink(const std::string& filename) 
        : _file(filename, std::ios::app) {}
    
    ~FileSink() override {
        if (_file.is_open()) {
            _file.close();
        }
    }
    
    void write(LogLevel level, const std::string& module, 
               const std::string& message, 
               const std::string& timestamp) override {
        if (_file.is_open()) {
            _file << "[" << timestamp << "] "
                  << "[" << to_string(level) << "] ";
            
            if (!module.empty()) {
                _file << "[" << module << "] ";
            }
            
            _file << message << std::endl;
        }
    }
    
    void flush() override {
        if (_file.is_open()) {
            _file.flush();
        }
    }
    
private:
    std::ofstream _file;
};

/**
 * @brief String buffer sink (for capturing logs in tests)
 */
class StringSink : public LogSink {
public:
    void write(LogLevel level, const std::string& module, 
               const std::string& message, 
               const std::string& timestamp) override {
        std::lock_guard<std::mutex> lock(_mutex);
        _buffer << "[" << timestamp << "] "
                << "[" << to_string(level) << "] ";
        
        if (!module.empty()) {
            _buffer << "[" << module << "] ";
        }
        
        _buffer << message << "\n";
    }
    
    void flush() override {}
    
    std::string getContents() const {
        std::lock_guard<std::mutex> lock(_mutex);
        return _buffer.str();
    }
    
    void clear() {
        std::lock_guard<std::mutex> lock(_mutex);
        _buffer.str("");
        _buffer.clear();
    }
    
private:
    mutable std::mutex _mutex;
    std::ostringstream _buffer;
};

/**
 * @brief Main logger class (singleton)
 */
class Logger {
public:
    /**
     * @brief Get the singleton logger instance
     */
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }
    
    // Prevent copying
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    
    /**
     * @brief Set the minimum log level
     * @param level Messages below this level will be ignored
     */
    void setLevel(LogLevel level) {
        std::lock_guard<std::mutex> lock(_mutex);
        _level = level;
    }
    
    /**
     * @brief Get the current log level
     */
    LogLevel getLevel() const {
        return _level;
    }
    
    /**
     * @brief Add a log sink
     * @param sink Shared pointer to the sink
     */
    void addSink(std::shared_ptr<LogSink> sink) {
        std::lock_guard<std::mutex> lock(_mutex);
        _sinks.push_back(sink);
    }
    
    /**
     * @brief Clear all sinks
     */
    void clearSinks() {
        std::lock_guard<std::mutex> lock(_mutex);
        _sinks.clear();
    }
    
    /**
     * @brief Log a message
     * @param level Log level
     * @param module Module name (optional)
     * @param message The message to log
     */
    void log(LogLevel level, const std::string& module, const std::string& message) {
        if (level < _level) {
            return;
        }
        
        std::lock_guard<std::mutex> lock(_mutex);
        
        auto timestamp = getTimestamp();
        
        for (auto& sink : _sinks) {
            sink->write(level, module, message, timestamp);
        }
    }
    
    /**
     * @brief Flush all sinks
     */
    void flush() {
        std::lock_guard<std::mutex> lock(_mutex);
        for (auto& sink : _sinks) {
            sink->flush();
        }
    }
    
    // Convenience methods
    void trace(const std::string& message, const std::string& module = "") {
        log(LogLevel::TRACE, module, message);
    }
    
    void debug(const std::string& message, const std::string& module = "") {
        log(LogLevel::DEBUG, module, message);
    }
    
    void info(const std::string& message, const std::string& module = "") {
        log(LogLevel::INFO, module, message);
    }
    
    void warning(const std::string& message, const std::string& module = "") {
        log(LogLevel::WARNING, module, message);
    }
    
    void error(const std::string& message, const std::string& module = "") {
        log(LogLevel::ERROR, module, message);
    }
    
    void critical(const std::string& message, const std::string& module = "") {
        log(LogLevel::CRITICAL, module, message);
    }
    
private:
    Logger() {
        // Default: console sink at ERROR level
        _sinks.push_back(std::make_shared<ConsoleSink>());
        _level = LogLevel::ERROR;
    }
    
    static std::string getTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
            << '.' << std::setfill('0') << std::setw(3) << ms.count();
        return oss.str();
    }
    
    std::mutex _mutex;
    LogLevel _level;
    std::vector<std::shared_ptr<LogSink>> _sinks;
};

// ============================================================================
// Logging Macros
// ============================================================================

#define OM_LOG(level, message) \
    OpenMagnetics::Logger::getInstance().log(level, "", message)

#define OM_LOG_MODULE(level, module, message) \
    OpenMagnetics::Logger::getInstance().log(level, module, message)

#define OM_TRACE(message) \
    OpenMagnetics::Logger::getInstance().trace(message)

#define OM_DEBUG(message) \
    OpenMagnetics::Logger::getInstance().debug(message)

#define OM_INFO(message) \
    OpenMagnetics::Logger::getInstance().info(message)

#define OM_WARNING(message) \
    OpenMagnetics::Logger::getInstance().warning(message)

#define OM_ERROR(message) \
    OpenMagnetics::Logger::getInstance().error(message)

#define OM_CRITICAL(message) \
    OpenMagnetics::Logger::getInstance().critical(message)

// Module-specific logging
#define OM_TRACE_M(module, message) \
    OpenMagnetics::Logger::getInstance().trace(message, module)

#define OM_DEBUG_M(module, message) \
    OpenMagnetics::Logger::getInstance().debug(message, module)

#define OM_INFO_M(module, message) \
    OpenMagnetics::Logger::getInstance().info(message, module)

#define OM_WARNING_M(module, message) \
    OpenMagnetics::Logger::getInstance().warning(message, module)

#define OM_ERROR_M(module, message) \
    OpenMagnetics::Logger::getInstance().error(message, module)

#define OM_CRITICAL_M(module, message) \
    OpenMagnetics::Logger::getInstance().critical(message, module)

} // namespace OpenMagnetics
