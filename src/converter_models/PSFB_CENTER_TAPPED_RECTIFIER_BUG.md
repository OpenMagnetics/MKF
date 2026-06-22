# ABT — PSFB center-tapped rectifier deck delivers ~half the designed Vout

- **Component:** `converter_models/PhaseShiftedFullBridge.cpp` (`Psfb::generate_ngspice_circuit`)
- **Severity:** High — silent-wrong-physics. Design and simulation disagree by ~2×; no error is raised.
- **Status:** Open (surfaced 2026-06-22 while porting PSFB into OpenConverters/Kirchhoff; needs MKF-owner decision before fix).
- **Scope of bug:** `BRectifierType::CENTER_TAPPED` only. `FULL_BRIDGE` is **correct** and `CURRENT_DOUBLER` is **not** affected by this specific defect (uses two inductors, different node structure — verify separately).

---

## 1. Summary

For `rectifierType = CENTER_TAPPED`, `Psfb::process_design_requirements()` sizes the turns ratio for a
**full-wave** center-tapped secondary (`Vo = Vin·Deff/n − Vd`), but `generate_ngspice_circuit()` emits a
netlist that is electrically a **half-wave** rectifier. The simulated output is therefore ≈ half the
designed/target output. A 48 V → 12 V design simulates at **5.68 V**.

The full-bridge rectifier path designs and simulates consistently (11.84 V for the same operating point),
which is what isolated the defect.

---

## 2. Root cause

A real center-tapped secondary is **two** half-windings sharing a center node:

```
   sec_a ──[½ winding]── CT ──[½ winding]── sec_b
                          │
                         (return / output ground)
   D1: sec_a → out        D2: sec_b → out
   On the +half D1 conducts (sec_a high); on the −half D2 conducts (sec_b high).
   => full-wave: output pulses on BOTH half-cycles.
```

The deck instead emits **one** full secondary winding `L_sec_o1 sec_a → sec_b` and fakes the center tap by
pinning it to **one end** (`sec_b`) through a 1 µΩ resistor. The relevant block
(`generate_ngspice_circuit`, the `else` / center-tapped branch, ~lines 1009–1029):

```spice
* Center-tapped rectifier (output 1)
D_r1_o1 rec_a_o1 out_rect_o1 DIDEAL     ; anode sec_a
D_r2_o1 rec_b_o1 out_rect_o1 DIDEAL     ; anode sec_b
Vct_o1  out_gnd_o1 sec_ct_o1 0          ; "center tap" node...
Rct_o1  sec_ct_o1 sec_b_o1 1u           ; ...tied to sec_b (= ONE END of the winding!)
L_out_o1 out_rect_o1 out_node_o1 <Lo>
C_out_o1 out_node_o1 out_gnd_o1 <Cout>
R_load_o1 out_node_o1 out_gnd_o1 <Rload>
Vout_sense_o1 out_gnd_o1 0 0
```

Because the "center tap" (`sec_ct_o1`) is `sec_b` (an **end** of the winding, not the midpoint), `D_r2`'s
anode and the output return are at essentially the same node. On the half-cycle where `sec_b` is the high
end, `D_r2` cannot deliver power (it sees ≈ the return potential); only `D_r1` (from `sec_a`) ever pushes
current into the output inductor. The result is single-diode / half-wave rectification → ~half the
volt-seconds → ~half Vo. The in-source comment on those lines even admits the simplification
("treat sec_b as the CT reference; the full L_sec winding here is shared").

`compute_turns_ratio(... CENTER_TAPPED)` returns `Vin·Deff/(Vo+Vd)` (full-wave assumption), so the design
expects the full-wave output the deck never produces.

---

## 3. Reproduction (self-contained, standalone — no Kirchhoff needed)

Requires MKF built at `$MKF_ROOT/build` (`libMKF.so` present), same as the converter PtP tests.

`/tmp/repro_psfb_ct.cpp`:

```cpp
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include "converter_models/PhaseShiftedFullBridge.h"
using namespace MAS; using namespace OpenMagnetics;

static double twm(const std::vector<double>& t, const std::vector<double>& v){
    if(t.size()<2) return 0; double a=0;
    for(size_t i=1;i<t.size();++i) a+=0.5*(v[i]+v[i-1])*(t[i]-t[i-1]);
    return a/(t.back()-t.front());
}
static void run(BRectifierType rt, const char* tag){
    const double Vin=48, Vout=12, Iout=2, Fs=100e3, phi=126;   // D_cmd = 126/180 = 0.7
    Psfb p;
    DimensionWithTolerance iv; iv.set_nominal(Vin); iv.set_minimum(Vin*0.95); iv.set_maximum(Vin*1.05);
    p.set_input_voltage(iv);
    p.set_rectifier_type(rt);
    PsfbOperatingPoint op;
    op.set_output_voltages({Vout}); op.set_output_currents({Iout});
    op.set_switching_frequency(Fs); op.set_phase_shift(phi); op.set_ambient_temperature(25);
    p.set_operating_points(std::vector<PsfbOperatingPoint>{op});

    auto dr = p.process_design_requirements();
    std::vector<double> tr; for(auto& t: dr.get_turns_ratios()) tr.push_back(t.get_nominal().value());
    double Lm = dr.get_magnetizing_inductance().get_nominal().value();
    p.set_num_steady_state_periods(2400);   // Cout=100u, Rload=6 -> RC=600us; 2400@100kHz=24ms >> 10*RC
    p.set_num_periods_to_extract(1);
    auto w = p.simulate_and_extract_topology_waveforms(tr, Lm, 1).at(0);
    auto& vo = w.get_output_voltages().at(0);
    std::cout << tag << ": designed n=" << tr[0] << "  -> simulated Vout="
              << twm(vo.get_time().value(), vo.get_data()) << " V (target 12 V)\n";
}
int main(){ run(BRectifierType::CENTER_TAPPED,"CENTER_TAPPED"); run(BRectifierType::FULL_BRIDGE,"FULL_BRIDGE"); }
```

Build + run:

```bash
MKF=/home/alf/OpenMagnetics/MKF ; B=$MKF/build
c++ -std=gnu++23 -O2 \
  -I"$B/_deps/json-src/include/nlohmann" -I"$B/_deps/json-src/include" \
  -I"$B/_deps/magic-enum-src/include/magic_enum" -I"$B/_deps/spline-src/src" \
  -I"$B/_deps/svg-src/src" -I"$B/_deps/eigen-src" \
  -I"$B/MAS" -I"$B/CAS" -I"$B/generated" \
  -I"$MKF/src" -I"$MKF/src/advisers" -I"$MKF/src/constructive_models" \
  -I"$MKF/src/converter_models" -I"$MKF/src/physical_models" \
  -I"$MKF/src/processors" -I"$MKF/src/support" \
  -I"$B/_cmrc/include" -I"$B/_deps/rapidfuzz-src/rapidfuzz/.." -I"$MKF/tests" \
  /tmp/repro_psfb_ct.cpp -o /tmp/repro_psfb_ct -L"$B" -lMKF -Wl,-rpath,"$B"
/tmp/repro_psfb_ct
```

### Observed (the bug)

```
CENTER_TAPPED: designed n=2.59  -> simulated Vout=5.68 V (target 12 V)   <-- ~half
FULL_BRIDGE:   designed n=2.45  -> simulated Vout=11.84 V (target 12 V)  <-- correct
```

Full per-run design numbers for reference:
- CENTER_TAPPED: n=2.59, Lm=1069.24 µH, Lr=2 µH, Lo=62.57 µH, Deff=0.6872 → **Vout 5.68 V, Iout 0.95 A, η 0.66**
- FULL_BRIDGE:   n=2.45, Lm=1009.92 µH, Lr=2 µH, Lo=62.72 µH, Deff=0.6864 → **Vout 11.84 V, Iout 1.97 A, η 0.67**

### Expected (after fix)

`CENTER_TAPPED` should also simulate Vout ≈ 12 V (within the same tolerance the FULL_BRIDGE path meets),
since `compute_turns_ratio(CENTER_TAPPED)` already sizes n for full-wave output.

---

## 4. Proposed fix (for review — do not apply without MKF-owner sign-off)

Emit a **true** center-tapped secondary: two half-windings coupled to the primary with the center node as
the rectifier return. Concretely, replace the single `L_sec_o1 sec_a → sec_b` + fake-CT block with two
half-windings `L_sec_o1a sec_a → sec_ct` and `L_sec_o1b sec_ct → sec_b`, each with its own `K` to the
primary (and to each other), and route both rectifier diodes into the output inductor with `sec_ct` as the
output return:

```spice
L_sec_o1a sec_a_o1  sec_ct_o1  <Ls_half>
L_sec_o1b sec_ct_o1 sec_b_o1   <Ls_half>
K..a L_pri L_sec_o1a 0.9999
K..b L_pri L_sec_o1b 0.9999
K..ab L_sec_o1a L_sec_o1b 0.9999
D_r1_o1 sec_a_o1 out_rect_o1 DIDEAL
D_r2_o1 sec_b_o1 out_rect_o1 DIDEAL
L_out_o1 out_rect_o1 out_node_o1 <Lo>
C_out_o1 out_node_o1 sec_ct_o1  <Cout>
R_load_o1 out_node_o1 sec_ct_o1 <Rload>
Vout_sense_o1 sec_ct_o1 0 0     ; CT = output return = node 0
```

`Ls_half` is the per-half-winding inductance. For a center-tapped secondary, `compute_turns_ratio` defines
n per **half**-secondary (`Vo = Vin·Deff/n`), so each half scales as `Lm/n²` (same as the current single
winding — the half winding carries n turns, the two halves are the full CT secondary). Verify the
half-winding inductance and the K-matrix against `LeakageInductance` / the magnetic the design produces;
sanity-check the emitted network delivers full-wave (output pulses on both half-cycles) at the DC limit.

**Cross-checks the fix must pass** (per the numeric-correctness guardrails):
- The repro above prints CENTER_TAPPED Vout ≈ 12 V.
- The existing PSFB-related PtP / circuit-simulator-exporter test suites still pass.
- Whatever the magnetic-view secondary excitations (`process_operating_point_for_input_voltage`,
  CENTER_TAPPED branch with its two `Secondary Na/Nb` half-windings) assume must stay consistent with the
  new two-half-winding deck (they already model two half-windings — the **deck** was the outlier).

**Alternative (smaller, if two-winding emission is too invasive):** keep one winding but make it a real CT
by exposing the **midpoint** as a node (split `L_sec` into two series-coupled halves) rather than pinning
`sec_b`. Either way the defining property to restore is: **power delivered on both half-cycles.**

---

## 5. Notes / context

- Found while porting PSFB into `OpenConverters/Kirchhoff` (10th MKF-equivalence topology). The Kirchhoff
  port deliberately uses **FULL_BRIDGE** rectification (the correct deck) and matches MKF to −0.9 % Vout, so
  Kirchhoff is **not** blocked by this — but anyone using MKF's PSFB with the default/CENTER_TAPPED setting
  gets silently-wrong (~half) output.
- `CENTER_TAPPED` is the **default** (`get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED)`), so
  this is the path most callers hit unless they explicitly select FULL_BRIDGE.
- Guardrail note: do **not** "fix" this by changing `compute_turns_ratio` to halve n (that would make the
  broken half-wave deck read 12 V while emitting a physically wrong circuit and wrong winding stresses).
  The **deck** is the bug; fix the deck.
