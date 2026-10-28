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
#include <atomic>
#include <chrono>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <memory>
#include <vector>
#include <optional>

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
    virtual void write(LogLevel level, const std::string& moduleOfOrigin, 
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
    
    void write(LogLevel level, const std::string& moduleOfOrigin, 
               const std::string& message, 
               const std::string& timestamp) override {
        std::ostream& out = (level >= LogLevel::ERROR) ? std::cerr : std::cout;
        
        if (_useColors) {
            out << getColorCode(level);
        }
        
        out << "[" << timestamp << "] "
            << "[" << to_string(level) << "] ";
        
        if (!moduleOfOrigin.empty()) {
            out << "[" << moduleOfOrigin << "] ";
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
    
    void write(LogLevel level, const std::string& moduleOfOrigin, 
               const std::string& message, 
               const std::string& timestamp) override {
        if (_file.is_open()) {
            _file << "[" << timestamp << "] "
                  << "[" << to_string(level) << "] ";
            
            if (!moduleOfOrigin.empty()) {
                _file << "[" << moduleOfOrigin << "] ";
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
    void write(LogLevel level, const std::string& moduleOfOrigin, 
               const std::string& message, 
               const std::string& timestamp) override {
        std::lock_guard<std::mutex> lock(_mutex);
        _buffer << "[" << timestamp << "] "
                << "[" << to_string(level) << "] ";
        
        if (!moduleOfOrigin.empty()) {
            _buffer << "[" << moduleOfOrigin << "] ";
        }
        
        _buffer << message << "\n";
    }
    
    // No-op: StringSink buffers in memory and doesn't need flushing to external storage
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
 * @brief Parse a level name ("TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "OFF").
 * Throws std::invalid_argument for anything else: an unknown level is a caller bug, not a default.
 */
inline LogLevel log_level_from_string(const std::string& name) {
    if (name == "TRACE") return LogLevel::TRACE;
    if (name == "DEBUG") return LogLevel::DEBUG;
    if (name == "INFO") return LogLevel::INFO;
    if (name == "WARNING") return LogLevel::WARNING;
    if (name == "ERROR") return LogLevel::ERROR;
    if (name == "CRITICAL") return LogLevel::CRITICAL;
    if (name == "OFF") return LogLevel::OFF;
    throw std::invalid_argument("Unknown log level '" + name +
                                "'; expected TRACE, DEBUG, INFO, WARNING, ERROR, CRITICAL or OFF");
}

/**
 * @brief The name a log entry's module carries in every sink and in the collector.
 *
 * The collected records reach the web's Engine log, where the module is the "source" column a
 * design engineer reads, so each entry names the part of the design (or the tool) it is about,
 * not the MKF class that wrote it. Code passes either an MKF identifier listed here, which is
 * translated, or one of the names it translates to, which is kept. An empty module or one that
 * is not listed throws std::invalid_argument, whatever the entry's level: there is no blank or
 * generic source. A new module is added here, with a name chosen by what its entries are about.
 */
inline const std::string& user_facing_log_module(const std::string& moduleOfOrigin, const std::string& message) {
    static const std::unordered_map<std::string, std::string> userFacingNames = [] {
        std::unordered_map<std::string, std::string> names = {
            {"CoreAdviser", "Core adviser"},
            {"MagneticAdviser", "Magnetic adviser"},
            {"CoilAdviser", "Coil adviser"},
            {"WireAdviser", "Wire adviser"},
            {"MagneticFilter", "Adviser filters"},
            {"MagneticFilterTemperature", "Adviser filters"},
            {"AdviserScoring", "Adviser scoring"},
            {"CoreCrossReferencer", "Core cross-reference"},
            {"CoreMaterialCrossReferencer", "Core material cross-reference"},
            {"Coil", "Winding"},
            {"WindingLayout", "Winding layout"},
            {"Bobbin", "Bobbins"},
            {"CoreCatalogue", "Core catalogue"},
            {"Core", "Core"},
            {"Inputs", "Operating points"},
            {"MaterialDataExtrapolation", "Material data extrapolation"},
            {"AmplitudePermeability", "Material data"},
            {"ComplexPermeability", "Material data"},
            {"ThermalResistance", "Material data"},
            {"Temperature", "Thermal model"},
            {"InsulationIEC60664Model", "Insulation"},
            {"InsulationIEC61558Model", "Insulation"},
            {"Painter", "Drawing"},
            {"Logger", "Engine log"},
        };
        std::vector<std::pair<std::string, std::string>> kept;
        for (const auto& [identifier, name] : names) {
            kept.emplace_back(name, name);
        }
        names.insert(kept.begin(), kept.end());
        return names;
    }();
    if (moduleOfOrigin.empty()) {
        throw std::invalid_argument("A log entry has no module; every entry must name what it is about "
                                    "(e.g. \"Winding\", \"Operating points\"). Entry: " + message);
    }
    auto found = userFacingNames.find(moduleOfOrigin);
    if (found == userFacingNames.end()) {
        throw std::invalid_argument("Log module '" + moduleOfOrigin + "' has no user-facing name; add it to "
                                    "user_facing_log_module in src/support/Logger.h. Entry: " + message);
    }
    return found->second;
}

/**
 * @brief One structured record kept by the Logger's collector (see Logger::enableCollector).
 * Records with the same level, module and message are merged; `count` says how many times it was logged.
 * `detail` is the part of an entry that varies between repeats of the same message (the value asked, e.g.
 * "asked at 1.2 MHz"): it is not part of the merge key, so repeats stay one record, and it holds the
 * detail of the latest repeat. Empty when the entry has none.
 */
struct LogRecord {
    LogLevel level;
    std::string moduleOfOrigin;
    std::string message;
    size_t count = 1;
    std::string detail;
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
        _sinks.push_back({sink, std::nullopt});
    }

    /**
     * @brief Add a log sink with its own minimum level
     *
     * The sink receives every entry at or above `minimumLevel`, whatever setLevel() says, and
     * setLevel() never changes what it receives. It lets one reader (read_log) see WARNING entries
     * without lowering the level of every other sink for the rest of the process (ABT #1643).
     */
    void addSink(std::shared_ptr<LogSink> sink, LogLevel minimumLevel) {
        std::lock_guard<std::mutex> lock(_mutex);
        _sinks.push_back({sink, minimumLevel});
        if (minimumLevel < _lowestOwnSinkLevel) {
            _lowestOwnSinkLevel = minimumLevel;
        }
    }
    
    /**
     * @brief Clear all sinks
     */
    void clearSinks() {
        std::lock_guard<std::mutex> lock(_mutex);
        _sinks.clear();
        _lowestOwnSinkLevel = LogLevel::OFF;
    }
    
    /**
     * @brief Log a message
     * @param level Log level
     * @param moduleIdentifier What the entry is about (see user_facing_log_module); required
     * @param message The message to log
     */
    void log(LogLevel level, const std::string& moduleIdentifier, const std::string& message) {
        log(level, moduleIdentifier, message, std::string());
    }

    /**
     * @brief Log a message with a detail that varies between its repeats (see LogRecord::detail)
     *
     * The sinks write "message (detail)"; the collector keeps message and detail apart and merges
     * repeats of the message whatever their detail.
     */
    void log(LogLevel level, const std::string& moduleIdentifier, const std::string& message, const std::string& detail) {
        // Checked before the level filter, so a missing or unlisted module is refused even for an
        // entry nobody would have seen.
        const std::string& moduleOfOrigin = user_facing_log_module(moduleIdentifier, message);
        // The sinks see what passes the logger level, exactly as before the collector existed, except
        // a sink added with its own level, which sees what passes that one; the collector, when
        // enabled, sees what passes its own level. None changes the others.
        const LogLevel loggerLevel = _level;
        const bool toSinks = !(level < loggerLevel) || (level < LogLevel::OFF && !(level < _lowestOwnSinkLevel));
        const bool toCollector = _collectorEnabled && level < LogLevel::OFF && !(level < _collectorLevel);
        if (!toSinks && !toCollector) {
            return;
        }

        std::lock_guard<std::mutex> lock(_mutex);

        if (toSinks) {
            std::string timestamp;
            for (auto& entry : _sinks) {
                const LogLevel threshold = entry.ownLevel ? entry.ownLevel.value() : loggerLevel;
                if (level < threshold) {
                    continue;
                }
                if (timestamp.empty()) {
                    timestamp = getTimestamp();
                }
                entry.sink->write(level, moduleOfOrigin, detail.empty() ? message : message + " (" + detail + ")", timestamp);
            }
        }
        if (toCollector) {
            collect(level, moduleOfOrigin, message, detail);
        }
    }

    /**
     * @brief Start keeping structured records (level, module, message) for a caller to drain.
     *
     * Process-wide, like the logger itself, so records logged on adviser worker threads are kept
     * too. Independent of setLevel(): the console sink keeps its own level (ERROR by default) while
     * the collector keeps everything at or above `minimumLevel`. Enabling an already-enabled
     * collector only changes its level; what it holds is kept until drained.
     * Identical records are merged (LogRecord::count); after kMaximumCollectedRecords distinct
     * records further ones are counted, not kept, and drainCollected() reports how many.
     */
    void enableCollector(LogLevel minimumLevel = LogLevel::WARNING) {
        std::lock_guard<std::mutex> lock(_mutex);
        _collectorLevel = minimumLevel;
        _collectorEnabled = true;
    }

    /**
     * @brief Stop collecting and discard whatever was collected and not drained.
     */
    void disableCollector() {
        std::lock_guard<std::mutex> lock(_mutex);
        _collectorEnabled = false;
        _collected.clear();
        _collectedIndex.clear();
        _collectedDropped = 0;
    }

    bool isCollectorEnabled() const {
        return _collectorEnabled;
    }

    LogLevel getCollectorLevel() const {
        return _collectorLevel;
    }

    /**
     * @brief Hand back every record collected since the last drain, in first-logged order, and
     * empty the collector. When records were dropped past the cap, a final WARNING from module
     * "Logger" says how many.
     */
    std::vector<LogRecord> drainCollected() {
        std::lock_guard<std::mutex> lock(_mutex);
        std::vector<LogRecord> drained;
        drained.swap(_collected);
        _collectedIndex.clear();
        if (_collectedDropped > 0) {
            drained.push_back({LogLevel::WARNING, user_facing_log_module("Logger", ""),
                               std::to_string(_collectedDropped) + " further log records were dropped after the first " +
                               std::to_string(kMaximumCollectedRecords) + " distinct ones", 1, ""});
            _collectedDropped = 0;
        }
        return drained;
    }

    static constexpr size_t kMaximumCollectedRecords = 1000;
    
    /**
     * @brief Flush all sinks
     */
    void flush() {
        std::lock_guard<std::mutex> lock(_mutex);
        for (auto& entry : _sinks) {
            entry.sink->flush();
        }
    }
    
    // Convenience methods
    void trace(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::TRACE, moduleOfOrigin, message);
    }
    
    void debug(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::DEBUG, moduleOfOrigin, message);
    }
    
    void info(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::INFO, moduleOfOrigin, message);
    }
    
    void warning(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::WARNING, moduleOfOrigin, message);
    }

    void warning(const std::string& message, const std::string& moduleOfOrigin, const std::string& detail) {
        log(LogLevel::WARNING, moduleOfOrigin, message, detail);
    }
    
    void error(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::ERROR, moduleOfOrigin, message);
    }
    
    void critical(const std::string& message, const std::string& moduleOfOrigin) {
        log(LogLevel::CRITICAL, moduleOfOrigin, message);
    }
    
private:
    Logger() {
        // Default: console sink at ERROR level
        _sinks.push_back({std::make_shared<ConsoleSink>(), std::nullopt});
        _level = LogLevel::ERROR;
    }
    
    static std::string getTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        
        std::tm tm_buf{};
#ifdef _WIN32
        localtime_s(&tm_buf, &time);
#else
        localtime_r(&time, &tm_buf);
#endif
        std::ostringstream oss;
        oss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S")
            << '.' << std::setfill('0') << std::setw(3) << ms.count();
        return oss.str();
    }
    
    // Caller holds _mutex.
    void collect(LogLevel level, const std::string& moduleOfOrigin, const std::string& message, const std::string& detail) {
        std::string key;
        key.reserve(moduleOfOrigin.size() + message.size() + 3);
        key += static_cast<char>('0' + static_cast<int>(level));
        key += '\x1f';
        key += moduleOfOrigin;
        key += '\x1f';
        key += message;
        auto found = _collectedIndex.find(key);
        if (found != _collectedIndex.end()) {
            _collected[found->second].count++;
            _collected[found->second].detail = detail;
            return;
        }
        if (_collected.size() >= kMaximumCollectedRecords) {
            _collectedDropped++;
            return;
        }
        _collectedIndex.emplace(std::move(key), _collected.size());
        _collected.push_back({level, moduleOfOrigin, message, 1, detail});
    }

    std::mutex _mutex;
    std::atomic<LogLevel> _level;
    struct SinkEntry {
        std::shared_ptr<LogSink> sink;
        std::optional<LogLevel> ownLevel;  // empty: follows setLevel()
    };
    std::vector<SinkEntry> _sinks;
    std::atomic<LogLevel> _lowestOwnSinkLevel{LogLevel::OFF};
    std::atomic<bool> _collectorEnabled{false};
    std::atomic<LogLevel> _collectorLevel{LogLevel::WARNING};
    std::vector<LogRecord> _collected;
    std::unordered_map<std::string, size_t> _collectedIndex;
    size_t _collectedDropped = 0;
};

// ============================================================================
// Logging Macros
// ============================================================================

#define OM_LOG_MODULE(level, moduleName, message) \
    OpenMagnetics::Logger::getInstance().log(level, moduleName, message)

// Module-specific logging
#define OM_TRACE_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().trace(message, moduleName)

#define OM_DEBUG_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().debug(message, moduleName)

#define OM_INFO_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().info(message, moduleName)

#define OM_WARNING_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().warning(message, moduleName)

#define OM_ERROR_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().error(message, moduleName)

#define OM_CRITICAL_M(moduleName, message) \
    OpenMagnetics::Logger::getInstance().critical(message, moduleName)

} // namespace OpenMagnetics
