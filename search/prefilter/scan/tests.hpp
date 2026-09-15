#pragma once

// Minimal harness for the scan tests: a failure counter and a generator.

#include <cstdint>
#include <cstdio>

inline int g_failures = 0;

#define CHECK(cond)                                                                   \
    do {                                                                              \
        if (!(cond)) {                                                                \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            ++g_failures;                                                             \
        }                                                                             \
    } while (0)

#define CHECK_MSG(cond, ...)                                                          \
    do {                                                                              \
        if (!(cond)) {                                                                \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s: ", __FILE__, __LINE__, #cond); \
            std::fprintf(stderr, __VA_ARGS__);                                        \
            std::fprintf(stderr, "\n");                                               \
            ++g_failures;                                                             \
        }                                                                             \
    } while (0)

struct Rng {
    uint64_t s;
    uint64_t next() {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        return s * 0x2545F4914F6CDD1DULL;
    }
    uint64_t below(uint64_t n) { return next() % n; }
};

inline int finish(const char* suite) {
    if (g_failures == 0)
        std::printf("%s: ok\n", suite);
    else
        std::printf("%s: %d failures\n", suite, g_failures);
    return g_failures == 0 ? 0 : 1;
}
