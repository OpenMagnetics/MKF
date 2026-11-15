#pragma once

#include <functional>
#include <string>

namespace OpenMagnetics {

// One top-level call (an adviser's advise, magnetic_autocomplete, mas_autocomplete) logs a warning about an excitation
// once, however many times its solves evaluate that excitation (ABT #2160: the core adviser's turns solve evaluated the
// magnetizing current of one operating point ~370,000 times and logged the same volt-second imbalance warning on every
// evaluation). The entry points open a scope; nested scopes share the outermost one, which ends with its call, so
// nothing carries over to the next call. With no scope open every warning is logged, as before.
// A notice that summarises many uses (material data used outside its tables, ABT #2410: one simulation asked a loss
// factor 4,128 times outside its data) is gathered while the scope is open and logged once when it closes, through
// at_close().
class LogOnceScope {
    public:
        LogOnceScope();
        ~LogOnceScope();
        LogOnceScope(const LogOnceScope&) = delete;
        LogOnceScope& operator=(const LogOnceScope&) = delete;
        // True when the warning `key` (what it is about and what it says) should be logged now: no scope is open on
        // this thread, or it is the first time this scope sees the key.
        static bool first_time(const std::string& key);
        // True when a scope is open on this thread.
        static bool is_open();
        // Runs `action` on this thread when the outermost open scope closes, after the scope is gone (so what the action
        // logs is logged directly). Actions run in the order they were added. Throws std::logic_error when no scope is open.
        static void at_close(std::function<void()> action);
    private:
        bool _owner = false;
};

} // namespace OpenMagnetics
