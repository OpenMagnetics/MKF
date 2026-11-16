// ABT #2274 (Alf, 2026-11-07, rule 3: "a misfit that neither a lane on its own side nor an allowed face move
// cures is reported"). PRIVATE to Coil.cpp and CoilPins.cpp; not part of MKF's interface.
//
// Coil::terminal_lane_plan is static: it sees the routes, their wires and their turn paths, but not the
// bobbin's winding window. Its bundle-order choice (a winding's entrance and exit bundles on one half: the
// entrance-first or the exit-first order) must still prefer the order whose leads keep their far-half
// terminal turns inside the window, which only the lane settle (Coil.cpp) can judge. The settle therefore
// lends terminal_lane_plan its judge for the duration of one call, through this scope; with no judge in
// force (any other caller) the choice is the one spec S0/S1 always made.
#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <set>
#include <string>

namespace OpenMagnetics::lane_plan_detail {

// Whether the terminal lead of routes[index] (the routes terminal_lane_plan was given) leaves the winding
// window with its lane at `lane`.
using LaneWindowJudge = std::function<bool(size_t index, double lane)>;

// The judge in force on this thread, or nullptr (defined in CoilPins.cpp).
const LaneWindowJudge*& lane_window_judge();

class ScopedLaneWindowJudge {
  public:
    explicit ScopedLaneWindowJudge(const LaneWindowJudge& judge) : _previous(lane_window_judge()) { lane_window_judge() = &judge; }
    ~ScopedLaneWindowJudge() { lane_window_judge() = _previous; }
    ScopedLaneWindowJudge(const ScopedLaneWindowJudge&) = delete;
    ScopedLaneWindowJudge& operator=(const ScopedLaneWindowJudge&) = delete;

  private:
    const LaneWindowJudge* _previous;
};

// ABT #2274 (Alf, 2026-11-08 TRY: "ideally the spiral ends on the connection"): the lane x at which the turn of
// routes[index] reaches its lead's row (the slot step is zero there), or nothing when the lead has no slot step.
// Only the settle (Coil.cpp) knows which leads take a slot step and their rows; it lends this judge the same way.
using LaneStepZero = std::function<std::optional<double>(size_t index)>;
const LaneStepZero*& lane_step_zero();

class ScopedLaneStepZero {
  public:
    explicit ScopedLaneStepZero(const LaneStepZero& judge) : _previous(lane_step_zero()) { lane_step_zero() = &judge; }
    ~ScopedLaneStepZero() { lane_step_zero() = _previous; }
    ScopedLaneStepZero(const ScopedLaneStepZero&) = delete;
    ScopedLaneStepZero& operator=(const ScopedLaneStepZero&) = delete;

  private:
    const LaneStepZero* _previous;
};


// ABT #2274 (Alf, 2026-11-08, after agent 21: "Check it, then move the lane"): the least distance from the crossing
// (u = |x|, outwards on the lead's own half) at which the lead of routes[index] may stand, or nothing. The settle sets
// it for a lead whose level row clashes (judge_level_rows) and for every lead of the same end standing outwards of it,
// so the fixed order holds; the plan then takes the lowest clear lane at or beyond it. Lent the same way.
using LaneFloor = std::function<std::optional<double>(size_t index, int half)>;
const LaneFloor*& lane_floor();

class ScopedLaneFloor {
  public:
    explicit ScopedLaneFloor(const LaneFloor& floor) : _previous(lane_floor()) { lane_floor() = &floor; }
    ~ScopedLaneFloor() { lane_floor() = _previous; }
    ScopedLaneFloor(const ScopedLaneFloor&) = delete;
    ScopedLaneFloor& operator=(const ScopedLaneFloor&) = delete;

  private:
    const LaneFloor* _previous;
};

// ABT #2274 (agent 27): how a terminal lead (an index into terminal_lane_plan's routes) at a lane is DRAWN where it leaves
// its turn (Coil's plan_slot_step_ramps): empty when its route is drawn as stated (a stub from its turn's height at the lane
// to its route's row); else its crossing row (rowY) and either its whole route moved to its turn's height at the lane (a
// level row, rowEnd false) or its route kept at the row with its turn ending on that row at the lane (a row end, rowEnd
// true, from fromX at fromY). Read by the tie judge (score_tie_pair) so it judges the copper that is drawn.
struct LaneLevelLead {
    double rowY = 0;
    bool rowEnd = false;
    double fromX = 0;
    // ABT #2274 (decision 2): on a ROUND turn the re-sloped wrap runs a whole revolution: from the path position fromAlong
    // (an exit: +L, the previous crossing; an entrance's row start: -L, the next crossing) at height fromY to its lane on
    // its row. NaN on a flat face.
    double fromAlong = std::numeric_limits<double>::quiet_NaN();
    double fromY = std::numeric_limits<double>::quiet_NaN();
};
using LaneLevelLeads = std::function<std::optional<LaneLevelLead>(size_t index, double lane)>;
const LaneLevelLeads*& lane_level_leads();
// ABT #2274 (Alf, 2026-11-09): the section-25 order notes already printed in this wind ([section25-order] once per pair
// per wind); Coil::wind clears it.
std::set<std::string>& section25_printed();

class ScopedLaneLevelLeads {
  public:
    explicit ScopedLaneLevelLeads(const LaneLevelLeads& leads) : _previous(lane_level_leads()) { lane_level_leads() = &leads; }
    ~ScopedLaneLevelLeads() { lane_level_leads() = _previous; }
    ScopedLaneLevelLeads(const ScopedLaneLevelLeads&) = delete;
    ScopedLaneLevelLeads& operator=(const ScopedLaneLevelLeads&) = delete;

  private:
    const LaneLevelLeads* _previous;
};

}  // namespace OpenMagnetics::lane_plan_detail
