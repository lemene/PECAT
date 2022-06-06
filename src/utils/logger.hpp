#pragma once

#include <stdio.h>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <mutex>

namespace fsa {

class Dumper {
public:

    class Stream {
    public:
        Stream(const std::string &fname="");
        ~Stream();
        void operator() (const char* const format, ...);
        void Write(const char* const format, va_list arglist);
    protected:
        FILE * file_{ nullptr };
        std::mutex mutex_;
    };

    Stream& operator [](const std::string &name);

    void SetDirectory(const std::string &dir) { directory_ = dir; }
    void SetLevel(int level) { level_ = level; }

    bool IsWorking() const { return level_ > 0; }
    std::string MakeFileName(const std::string &name) { return directory_ + "/debug_" + name + ".log"; }
protected:
    std::mutex mutex_;
    std::string directory_ { "." };
    int level_{ 0 };
    std::map<std::string, std::unique_ptr<Stream>> streams_;
    Stream empty_;
};

extern Dumper DUMPER;

class Logger {
public:
    enum Level {
        L_DEBUG=0,
        L_INFO,
        L_WARNING,
        L_ERROR
    };
    Logger();

    class Stream {
    public:
        Stream(Logger& logger, Level level) :logger_(logger), level_(level){}
        void operator() (const char* const format, ...);
    protected:
        Logger & logger_;
        Level level_;
    };

    void SetFileName(const std::string &fname);
    void SetLevel(Level level) { level_ = level; }
    void Log(Level level, const char* format, va_list arglist);
protected:
    std::mutex mutex_;
    std::string filename_;
    Level level_{ L_INFO };
    std::vector<std::string> levelname_{ "DEBUG", "INFO", "WARNING", "ERROR" };
    FILE *file_;
};

#define LOG(s) Logger::Stream(LOGGER, Logger::L_##s)
#define SET_LOG_LEVEL(s) LOGGER.SetLevel(Logger::L_##s)

extern Logger LOGGER;

extern bool print_rubbish;
#define DEBUG_printf if (print_rubbish) printf

#define DEBUG_local_printf if (local_print_rubbish) printf

} // namespace fsa {

