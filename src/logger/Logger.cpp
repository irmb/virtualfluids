#include "Logger.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/daily_file_sink.h>
#include <spdlog/sinks/basic_file_sink.h>

namespace vf::logging
{

    std::string Logger::logPath = {"logs/"};

    void Logger::initalizeLogger() 
    {
        updateDefaultLogger();

        // setting default log level to trace
        // levels: trace < debug < info < warn < error < critical
        spdlog::set_level(spdlog::level::trace);

        // setting the log pattern
        // formatting is documented here: https://github.com/gabime/spdlog/wiki/3.-Custom-formatting
        spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

        // according to the flush policy https://github.com/gabime/spdlog/wiki/7.-Flush-policy
        spdlog::flush_on(spdlog::level::info);
    }


    void Logger::changeLogPath(const std::string& path)
    {
        logPath = path;

        updateDefaultLogger();
    }


    void Logger::updateDefaultLogger()
    {
        // initialize stdout sink with colored output
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

        // initialize daily file sink
        // files will be written into "logs" folder relative to pwd. A new files is created at 0:00 o'clock.
        auto daily_file_sink = std::make_shared<spdlog::sinks::daily_file_sink_mt>(logPath + "daily.txt", 0, 0);

        // initialize last run file sink
        auto last_run_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logPath + "last_run.txt", true);

        // creating default logger with console and file sink
        spdlog::set_default_logger(std::make_shared<spdlog::logger>("default", spdlog::sinks_init_list({console_sink, daily_file_sink, last_run_file_sink})));
    }

}
