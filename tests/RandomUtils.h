#pragma once

#include <random>
#include <string>

namespace OpenMagnetics::TestUtils {

/**
 * @brief Thread-safe random number generator for tests
 * 
 * Replaces srand/rand with C++ <random> facilities as recommended by SonarQube
 */
class RandomGenerator {
public:
    static RandomGenerator& getInstance() {
        static thread_local RandomGenerator instance;
        return instance;
    }

    // Generate random integer in range [min, max]
    int randomInt(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(_engine);
    }

    // Generate random int64_t in range [min, max]
    int64_t randomInt64(int64_t min, int64_t max) {
        std::uniform_int_distribution<int64_t> dist(min, max);
        return dist(_engine);
    }

    // Generate random size_t in range [min, max]
    size_t randomSize(size_t min, size_t max) {
        std::uniform_int_distribution<size_t> dist(min, max);
        return dist(_engine);
    }

    // Generate random double in range [min, max]
    double randomDouble(double min, double max) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(_engine);
    }

    // Generate random bool
    bool randomBool() {
        std::uniform_int_distribution<int> dist(0, 1);
        return dist(_engine) == 1;
    }

    // Get underlying engine for custom distributions
    std::mt19937& engine() { return _engine; }

    // Reseed the generator (call with std::random_device{}() to opt back into fuzzing).
    void seed(unsigned int s) { _engine.seed(s); }

private:
    // Fixed default seed so tests are REPRODUCIBLE: the generator used to seed from
    // std::random_device, making every random-input test nondeterministic run-to-run. Combined
    // with the tests that can only pass (catch-and-continue / if(size>0) checks), that produced
    // the worst case — nondeterministic AND unable to fail. A caller that genuinely wants
    // fuzzing can still opt in via seed(std::random_device{}()). (abt #116)
    static constexpr unsigned int kDefaultSeed = 42;
    RandomGenerator() : _engine(kDefaultSeed) {}
    std::mt19937 _engine;
};

// Convenience functions for quick access
inline int randomInt(int min, int max) {
    return RandomGenerator::getInstance().randomInt(min, max);
}

inline int64_t randomInt64(int64_t min, int64_t max) {
    return RandomGenerator::getInstance().randomInt64(min, max);
}

inline size_t randomSize(size_t min, size_t max) {
    return RandomGenerator::getInstance().randomSize(min, max);
}

inline double randomDouble(double min, double max) {
    return RandomGenerator::getInstance().randomDouble(min, max);
}

inline bool randomBool() {
    return RandomGenerator::getInstance().randomBool();
}

} // namespace OpenMagnetics::TestUtils
