#pragma once

// The timer, the machine name and the CSV writer the sweeps share.

#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fsst::search::prefilter::scan::bench {

constexpr int PASSES = 5;

// Seconds of the fastest of PASSES runs after one untimed. A run the clock
// cannot resolve is repeated inside the timed region until it can.
template <typename Run>
double best(Run run) {
    constexpr double FLOOR = 1e-6;
    run();
    size_t reps = 1;
    while (true) {
        double fastest = 1e300;
        for (int pass = 0; pass < PASSES; ++pass) {
            auto start = std::chrono::steady_clock::now();
            for (size_t r = 0; r < reps; ++r) run();
            double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (seconds < fastest) fastest = seconds;
        }
        if (fastest > FLOOR || reps == (size_t{1} << 20)) return fastest / static_cast<double>(reps);
        reps *= 2;
    }
}

inline std::string shell(const char* command) {
    std::string out;
    if (FILE* pipe = popen(command, "r")) {
        char buf[256];
        while (fgets(buf, sizeof buf, pipe)) out += buf;
        pclose(pipe);
    }
    return out;
}

inline std::string labelled(const std::string& text, const std::string& label) {
    size_t at = text.find(label);
    if (at == std::string::npos) return "";
    size_t colon = text.find(':', at);
    size_t end = text.find('\n', colon);
    return colon == std::string::npos ? "" : text.substr(colon + 1, end - colon - 1);
}

inline std::string collapse(std::string s) {
    for (char& c : s)
        if (c == ',' || c == '\n' || c == '\t') c = ' ';
    std::string out;
    for (char c : s) {
        if (c == ' ' && (out.empty() || out.back() == ' ')) continue;
        out.push_back(c);
    }
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out;
}

// "<arch> <cpu model> <isa>", the key a fit groups rows by.
inline std::string machine() {
    std::string model = labelled(shell("cat /proc/cpuinfo 2>/dev/null"), "model name");
    if (model.empty()) model = labelled(shell("lscpu 2>/dev/null"), "Model name");
    if (model.empty()) model = shell("sysctl -n machdep.cpu.brand_string 2>/dev/null");
    if (collapse(model).empty()) model = "unknown cpu";
#if defined(__aarch64__)
    const char* arch = "aarch64";
    const char* isa = "neon";
#else
    const char* arch = "x86_64";
    const char* isa = "scalar";
#endif
    return collapse(std::string(arch) + " " + model + " " + isa);
}

inline std::filesystem::path output_dir() { return std::filesystem::path(BENCH_DIR) / "output"; }

inline std::string timestamp() {
    std::time_t now = std::time(nullptr);
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y-%m-%d_%H-%M-%S", std::localtime(&now));
    return buf;
}

// `output/<prefix>_<stamp>.csv` from a header line and one line per row.
inline std::filesystem::path write_csv(const std::string& prefix, const std::string& header,
                                       const std::vector<std::string>& rows) {
    std::filesystem::create_directories(output_dir());
    std::filesystem::path path = output_dir() / (prefix + "_" + timestamp() + ".csv");
    std::ofstream out(path);
    out << header << '\n';
    for (const std::string& row : rows) out << row << '\n';
    return path;
}

}  // namespace fsst::search::prefilter::scan::bench
