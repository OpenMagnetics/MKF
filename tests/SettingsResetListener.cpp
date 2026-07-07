#include "support/Settings.h"

#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

namespace {

// ABT #116: the global OpenMagnetics::Settings singleton leaked state between TEST_CASEs
// (TestCoil alone had ~464 set_* calls against ~158 resets), producing order-dependent
// outcomes — the proven "fails standalone, passes in-suite" class. Resetting the singleton
// to its defaults before every test case gives each test a known starting state; tests that
// need non-default settings must set them themselves (and now demonstrably do or fail).
class SettingsResetListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void testCaseStarting(Catch::TestCaseInfo const&) override {
        OpenMagnetics::Settings::GetInstance().reset();
    }
};

}  // namespace

CATCH_REGISTER_LISTENER(SettingsResetListener)
