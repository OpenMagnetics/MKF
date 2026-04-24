#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "Defaults.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace OpenMagnetics {

// Forward declaration for free function defined in CircuitSimulatorInterface.cpp.
std::string to_string(double d, size_t precision);

// ============================================================================
// NL5 circuit simulator exporter
// ============================================================================

// ---- NL5 Component builder ----
struct Nl5Cmp {
    std::string type;      // e.g. "R_R", "L_L", "C_C", "W_W", "W_T2"
    int id;
    std::string name;
    int node0, node1;
    int node2 = -1, node3 = -1; // for 4-terminal (K, G)
    // Value parameter: for R=r, L=l, C=c, K=k
    std::string valParam;  // parameter name ("r","l","c","k")
    std::string valStr;    // value as string
    // Extra parameters
    std::map<std::string, std::string> extra;
    // Position
    int px = 0, py = 0, angle = 0;
    bool flip = false;
};

// ---- Engineering notation formatter for NL5 _txt attributes ----
static std::string nl5_to_eng(double v, int sigfigs = 4) {
    if (v == 0.0) return "0";
    bool neg = v < 0;
    double av = std::abs(v);
    std::string suffix; double scaled;
    if      (av >= 1e12) { scaled = av/1e12; suffix = "T"; }
    else if (av >= 1e9)  { scaled = av/1e9;  suffix = "G"; }
    else if (av >= 1e6)  { scaled = av/1e6;  suffix = "Meg"; }
    else if (av >= 1e3)  { scaled = av/1e3;  suffix = "k"; }
    else if (av >= 1.0)  { scaled = av;       suffix = ""; }
    else if (av >= 1e-3) { scaled = av*1e3;   suffix = "m"; }
    else if (av >= 1e-6) { scaled = av*1e6;   suffix = "u"; }
    else if (av >= 1e-9) { scaled = av*1e9;   suffix = "n"; }
    else if (av >= 1e-12){ scaled = av*1e12;  suffix = "p"; }
    else                 { scaled = av*1e15;  suffix = "f"; }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*g", sigfigs, scaled);
    return (neg ? "-" : "") + std::string(buf) + suffix;
}

static std::string nl5_cmp_to_xml(const Nl5Cmp& cmp) {
    std::ostringstream ss;
    ss << "        <Cmp type=\"" << cmp.type << "\" id=\"" << cmp.id
       << "\" name=\"" << cmp.name << "\" descr=\"\" view=\"0\""
       << " node0=\"" << cmp.node0 << "\" node1=\"" << cmp.node1 << "\"";
    if (cmp.node2 >= 0) ss << " node2=\"" << cmp.node2 << "\"";
    if (cmp.node3 >= 0) ss << " node3=\"" << cmp.node3 << "\"";
    // model attribute: use explicit extra["model"] if present, otherwise derive from valParam
    if (cmp.extra.count("model")) {
        ss << " model=\"" << cmp.extra.at("model") << "\"";
    } else if (!cmp.valParam.empty()) {
        std::string mdl(1, static_cast<char>(std::toupper(static_cast<unsigned char>(cmp.valParam[0]))));
        ss << " model=\"" << mdl << "\"";
    }
    if (!cmp.valParam.empty()) {
        ss << " " << cmp.valParam << "=\"" << cmp.valStr << "\"";
        double dval = 0.0;
        try { dval = std::stod(cmp.valStr); } catch (...) {}
        ss << " " << cmp.valParam << "_txt=\"" << nl5_to_eng(dval) << "\"";
    }
    for (auto& [k, v] : cmp.extra) {
        if (k == "model") continue;  // already emitted above
        ss << " " << k << "=\"" << v << "\"";
    }
    ss << ">\n";
    ss << "          <Pos x=\"" << cmp.px << "\" y=\"" << cmp.py
       << "\" angle=\"" << cmp.angle << "\" flip=\"" << (cmp.flip ? "1" : "0") << "\" />\n";
    ss << "        </Cmp>\n";
    return ss.str();
}

static std::string nl5_label_to_xml(int id, const std::string& name, int node0, int px, int py) {
    std::ostringstream ss;
    ss << "        <Cmp type=\"label\" model=\"Label\" id=\"" << id
       << "\" name=\"" << name << "\" node0=\"" << node0 << "\">\n";
    ss << "          <Pos x=\"" << px << "\" y=\"" << py << "\" angle=\"0\" flip=\"0\" />\n";
    ss << "        </Cmp>\n";
    return ss.str();
}

// ============================================================================
// CircuitSimulatorExporterNl5Model::export_magnetic_as_subcircuit
// ============================================================================
// Uses NL5 Winding (W) components as documented in NL5 Reference:
// - Each winding is a W_Winding component with turns ratio
// - All winding "core" pins connect to a common core node
// - Magnetizing inductance connects from core node to ground
// - Leakage is modeled via series inductors
// - DC resistance and AC losses via series R/L/C networks
// ============================================================================
std::string CircuitSimulatorExporterNl5Model::export_magnetic_as_subcircuit(
    Magnetic magnetic,
    double frequency,
    double temperature,
    std::optional<std::string> filePathOrFile,
    CircuitSimulatorExporterCurveFittingModes mode)
{
    auto coil = magnetic.get_coil();
    size_t numWindings = coil.get_functional_description().size();

    // Resolve mode
    auto resolvedMode = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);

    // Precompute ladder coefficients (if needed)
    std::vector<std::vector<double>> acCoeffs;
    std::vector<FractionalPoleNetwork> fracNets;
    std::optional<FractionalPoleNetwork> coreFracNet;

    if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracNets = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
        catch (...) {}
    } else if (resolvedMode != CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
        acCoeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
            magnetic, temperature, CircuitSimulatorExporterCurveFittingModes::LADDER);
    }

    auto coreLossTopology_nl5 = static_cast<CoreLossTopology>(Settings::GetInstance().get_circuit_simulator_core_loss_topology());
    auto coreCoeffs = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature, coreLossTopology_nl5);

    double magnetizingInductance = resolve_dimensional_values(
        MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic)
            .get_magnetizing_inductance());
    auto leakageInductances = LeakageInductance()
        .calculate_leakage_inductance(magnetic, Defaults().measurementFrequency)
        .get_leakage_inductance_per_winding();

    // Get primary turns for reference
    double primaryTurns = static_cast<double>(coil.get_functional_description()[0].get_number_turns());

    // Check if saturation modeling is enabled
    auto& settings = Settings::GetInstance();
    bool includeSaturation = settings.get_circuit_simulator_include_saturation();
    SaturationParameters satParams;
    if (includeSaturation) {
        satParams = get_saturation_parameters(magnetic, temperature);
        // Saturation only supported for single winding (inductor) - coupled windings can't use saturating inductors
        if (satParams.valid && numWindings > 1) {
            // For transformers, fall back to linear model with warning
            // (NL5 XML comment not easily added, so just disable saturation)
            includeSaturation = false;
        }
    }

    // ---- NL5 XML generation ----
    // Node numbering scheme:
    //   0          = ground
    //   1..N       = winding positive terminal Pi+  (external)
    //   N+1..2N    = winding negative terminal Pi-  (external)
    //   50         = core node (common for all windings via ideal transformer)
    //   100+i      = internal node after Rdc (before ladder/fracpole), winding i
    //   200+i      = internal node after ladder/fracpole (before leakage L), winding i
    //   250+i      = internal node after leakage (connects to winding), winding i
    //   300+i*20+k = internal nodes within fracpole/ladder stages for winding i, stage k
    //   900+k      = internal nodes for core loss network

    int nextId = 1;   // component id counter
    std::string cmps; // accumulated <Cmp> XML

    // Core node where all windings connect via ideal transformers
    const int node_core = 50;

    auto add_R = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="R_R"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="r"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    auto add_L = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="L_L"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="l"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    auto add_C = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="C_C"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="c"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    // NL5 saturating inductor using L_Lsat model
    // Parameters: l=initial inductance, Lsat=saturated inductance, Iknee=knee current
    auto add_Lsat = [&](const std::string& name, int n0, int n1, double L0, double Lsat, double Iknee,
                        int px, int py, int angle=0) {
        Nl5Cmp c; c.type="L_Lsat"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="l"; c.valStr=to_string(L0, 12);
        // Lsat is the saturated (minimum) inductance
        c.extra["Lsat"] = to_string(Lsat, 12);
        c.extra["Lsat_txt"] = nl5_to_eng(Lsat);
        // Iknee is the current at knee point (where saturation starts)
        c.extra["Iknee"] = to_string(Iknee, 6);
        c.extra["Iknee_txt"] = nl5_to_eng(Iknee);
        c.extra["model"] = "Lsat";
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    // NL5 Winding component (W_W): ideal transformer with turns ratio
    // node0 = winding +, node1 = winding -, node2 = core
    auto add_Winding = [&](const std::string& name, int n_plus, int n_minus, int n_core,
                           double turns, int px, int py) {
        Nl5Cmp c; c.type="W_W"; c.id=nextId++; c.name=name;
        c.node0=n_plus; c.node1=n_minus; c.node2=n_core;
        c.valParam="w"; c.valStr=std::to_string(static_cast<int>(std::round(turns)));
        c.extra["model"] = "Winding";
        c.px=px; c.py=py;
        cmps += nl5_cmp_to_xml(c);
    };

    // Layout: each winding is 80 units wide, stacked vertically
    for (size_t wi = 0; wi < numWindings; ++wi) {
        int base_y = static_cast<int>(wi) * 80;
        int is = static_cast<int>(wi) + 1;
        std::string ws = std::to_string(is);

        int node_Pplus  = static_cast<int>(wi) + 1;           // Pi+  external
        int node_Pminus = static_cast<int>(numWindings) + is;  // Pi-  external
        int node_after_rdc  = 100 + is;                        // after Rdc
        int node_after_leak = 250 + is;                        // after leakage inductance

        double Rdc = WindingLosses::calculate_effective_resistance_of_winding(magnetic, wi, 0.1, temperature);
        double numTurns = static_cast<double>(coil.get_functional_description()[wi].get_number_turns());

        // Get leakage inductance for this winding (referred to primary)
        double leakageL = 0.0;
        if (wi > 0 && wi <= leakageInductances.size()) {
            leakageL = resolve_dimensional_values(leakageInductances[wi-1]);
            // Scale leakage to this winding's turns
            double turnsRatio = numTurns / primaryTurns;
            leakageL *= turnsRatio * turnsRatio;
        }

        // Add port labels (external terminals for SubCir pin mapping)
        // Use "p"/"n" instead of "+"/"-" because '+' is not URL-safe and gets
        // mangled when NL5 sends component names through URL requests.
        cmps += nl5_label_to_xml(nextId++, "P"+ws+"p", node_Pplus, 0, base_y);
        cmps += nl5_label_to_xml(nextId++, "P"+ws+"n", node_Pminus, 180, base_y+40);

        // 1. Rdc (DC resistance)
        add_R("Rdc"+ws, node_Pplus, node_after_rdc, Rdc, 20, base_y);

        int node_last = node_after_rdc;

        // 2. AC loss network (FRACPOLE or LADDER)
        if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE && wi < fracNets.size()) {
            // FRACPOLE: R/C stages
            const auto& net = fracNets[wi];
            for (size_t k = 0; k < net.stages.size(); ++k) {
                int node_stage = 300 + is * 20 + static_cast<int>(k);
                int px_R = 40 + static_cast<int>(k) * 16;
                // Series R from last node to stage node
                add_R("Rfp"+ws+"_"+std::to_string(k), node_last, node_stage,
                      net.stages[k].R * net.opts.coef, px_R, base_y);
                // Shunt C from stage node to ground
                add_C("Cfp"+ws+"_"+std::to_string(k), node_stage, 0,
                      net.stages[k].C / net.opts.coef, px_R+4, base_y+20, 90);
                node_last = node_stage;
            }
            if (net.Cinf > 0.0) {
                add_C("Cfpinf"+ws, node_last, 0, net.Cinf / net.opts.coef,
                      40 + static_cast<int>(net.stages.size())*16, base_y+20, 90);
            }
        } else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::LADDER
                   && !acCoeffs.empty() && wi < acCoeffs.size()) {
            // LADDER: R/L stages
            const auto& c = acCoeffs[wi];
            size_t stages = c.size() / 2;
            for (size_t k = 0; k < stages; ++k) {
                int node_stage = 300 + is * 20 + static_cast<int>(k);
                int px = 40 + static_cast<int>(k) * 16;
                add_L("Ll"+ws+"_"+std::to_string(k), node_last, node_stage, c[2*k], px, base_y);
                add_R("Rl"+ws+"_"+std::to_string(k), node_stage, 0, c[2*k+1], px+4, base_y+20, 90);
                node_last = node_stage;
            }
            // node_last already points to the last ladder stage node
        } else {
            // ANALYTICAL or fallback: no AC network, node_last stays at node_after_rdc
        }

        // 3. Leakage inductance (series with winding)
        if (leakageL > 0.0 && wi > 0) {
            add_L("Lleak"+ws, node_last, node_after_leak, leakageL, 100, base_y);
            node_last = node_after_leak;
        } else {
            // No leakage for primary or if leakage is zero
            node_after_leak = node_last;
        }

        // 4. Winding (ideal transformer to core)
        // NL5 Winding: node0=+, node1=-, node2=core
        add_Winding("W"+ws, node_after_leak, node_Pminus, node_core, numTurns, 130, base_y);
    }

    // 5. Magnetizing inductance (from core to ground)
    // This represents the magnetic energy storage in the core
    if (includeSaturation && satParams.valid) {
        // Use saturating inductor model
        const double mu0 = 4e-7 * M_PI;
        double Isat = satParams.Isat();
        // Calculate gap contribution for saturated inductance
        double effectiveMuR = satParams.Lmag / (mu0 * satParams.primaryTurns * satParams.primaryTurns * satParams.Ae / satParams.le);
        double gapFactor = std::min(0.99, 1.0 / effectiveMuR);
        double Lsat = satParams.Lmag * gapFactor;  // Saturated inductance (just the gap)
        // Ensure reasonable minimum values
        if (Isat < 1e-6) Isat = 1e-6;
        if (Lsat < 1e-12) Lsat = magnetizingInductance * 0.01;
        // NL5 Lsat model: L0 is initial inductance, Lsat is saturated inductance, Iknee is knee current
        add_Lsat("Lmag", node_core, 0, magnetizingInductance, Lsat, Isat,
                 160, static_cast<int>(numWindings) * 40);
    } else {
        add_L("Lmag", node_core, 0, magnetizingInductance, 160, static_cast<int>(numWindings) * 40);
    }

    // 6. Core loss network in parallel with Lmag
    if (!coreCoeffs.empty()) {
        if (coreLossTopology_nl5 == CoreLossTopology::ROSANO) {
            // Rosano: three stages in series from node_core to ground
            // Stage 1 (R): R1 series, Stage 2 (RL): L1||R2, Stage 3 (RLC): L2||(R3+C1)
            // coefficients: [R1, R2, L1, R3, L2, C1]
            int y_base = static_cast<int>(numWindings) * 40;
            double R1 = coreCoeffs[0];
            double R2 = coreCoeffs[1];
            double L1 = coreCoeffs[2];
            double R3 = coreCoeffs[3];
            double L2 = coreCoeffs[4];
            double C1 = coreCoeffs[5];

            // Stage 1: R1 from core to node 901
            if (R1 > 0) {
                add_R("Rcore_s1", node_core, 901, R1, 180, y_base);
            }
            // Stage 2: L1||R2 from node 901 to node 902
            if (R2 > 0 && L1 > 0) {
                add_R("Rcore_s2", 901, 902, R2, 180, y_base + 16);
                add_L("Lcore_s2", 901, 902, L1, 184, y_base + 24);
            }
            // Stage 3: L2||(R3+C1) from node 902 to ground
            if (R3 > 0 && L2 > 0 && C1 > 0) {
                add_L("Lcore_s3", 902, 0, L2, 180, y_base + 32);
                add_R("Rcore_s3", 902, 903, R3, 184, y_base + 40);
                add_C("Ccore_s3", 903, 0, C1, 188, y_base + 48);
            }
        } else {
            // Ridley: R-L ladder from node_core to ground
            int node_prev = 0;  // ground
            int stageCount = static_cast<int>(coreCoeffs.size() / 2);
            for (int k = stageCount - 1; k >= 0; --k) {
                double R = coreCoeffs[k * 2];
                double L = coreCoeffs[k * 2 + 1];
                if (R <= 0 || L <= 0 || R > 1e6 || L > 1) {
                    continue;
                }
                int node_stage = 900 + k;
                int y_offset = static_cast<int>(numWindings) * 40 + k * 16;
                add_R("Rcore" + std::to_string(k), node_stage, node_prev, R, 180, y_offset);
                int node_next = (k == 0) ? node_core : (900 + k - 1);
                add_L("Lcore" + std::to_string(k), node_next, node_prev, L, 184, y_offset + 8);
                node_prev = node_stage;
            }
        }
    }

    // ---- Assemble full XML ----
    std::ostringstream xml;
    xml << "<?xml version=\"1.0\"?>\n";
    xml << "<NL5>\n";
    xml << "  <c>Generated by OpenMagnetics - do not edit manually</c>\n";
    xml << "  <c>Component: " << magnetic.get_reference() << "</c>\n";
    xml << "  <Version ver=\"3\" rev=\"18\" core=\"78\" build=\"99\" />\n";
    xml << "  <Doc>\n";
    xml << "    <Cir mode=\"2\">\n";
    xml << "      <Cmps>\n";
    xml << cmps;
    xml << "      </Cmps>\n";
    xml << "      <Scopes />\n";
    xml << "      <Scripts />\n";
    xml << "    </Cir>\n";
    xml << "  </Doc>\n";
    xml << "</NL5>\n";

    std::string result = xml.str();

    if (filePathOrFile.has_value()) {
        std::ofstream outFile(filePathOrFile.value());
        outFile << result;
    }
    return result;
}

} // namespace OpenMagnetics
