#include "processors/CircuitSimulatorInterface.h"
#include "physical_models/WindingLosses.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <array>

namespace OpenMagnetics {

// Base64 encoding for PLECS InitializationCommands
std::string CircuitSimulatorExporterPlecsModel::encode_init_commands(const std::string& plainText) {
    static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string encoded;
    encoded.reserve(((plainText.size() + 2) / 3) * 4);

    for (size_t i = 0; i < plainText.size(); i += 3) {
        uint32_t n = static_cast<uint8_t>(plainText[i]) << 16;
        if (i + 1 < plainText.size()) n |= static_cast<uint8_t>(plainText[i + 1]) << 8;
        if (i + 2 < plainText.size()) n |= static_cast<uint8_t>(plainText[i + 2]);

        encoded += table[(n >> 18) & 0x3F];
        encoded += table[(n >> 12) & 0x3F];
        encoded += (i + 1 < plainText.size()) ? table[(n >> 6) & 0x3F] : '=';
        encoded += (i + 2 < plainText.size()) ? table[n & 0x3F] : '=';
    }

    // Split into 76-char lines with PLECS quoting
    std::string result;
    bool first = true;
    for (size_t i = 0; i < encoded.size(); i += 76) {
        std::string chunk = encoded.substr(i, 76);
        if (first) {
            result += "  InitializationCommands base64 \"" + chunk + "\"\n";
            first = false;
        } else {
            result += "\"" + chunk + "\"\n";
        }
    }
    return result;
}

std::string CircuitSimulatorExporterPlecsModel::emit_header(const std::string& name) {
    std::ostringstream ss;
    ss << "Plecs {\n";
    ss << "  Name          \"" << name << "\"\n";
    ss << "  Version       \"4.7\"\n";
    ss << "  CircuitModel  \"ContStateSpace\"\n";
    ss << "  StartTime     \"0.0\"\n";
    ss << "  TimeSpan      \"tsim\"\n";
    ss << "  Timeout       \"\"\n";
    ss << "  Solver        \"dopri\"\n";
    ss << "  MaxStep       \"1\"\n";
    ss << "  InitStep      \"1\"\n";
    ss << "  FixedStep     \"1e-3\"\n";
    ss << "  Refine        \"1\"\n";
    ss << "  ZCStepSize    \"1e-9\"\n";
    ss << "  RelTol        \"1e-3\"\n";
    ss << "  AbsTol        \"-1\"\n";
    ss << "  TurnOnThreshold \"0\"\n";
    ss << "  SyncFixedStepTasks \"2\"\n";
    ss << "  UseSingleCommonBaseRate \"2\"\n";
    ss << "  LossVariableLimitExceededMsg \"3\"\n";
    ss << "  NegativeSwitchLossMsg \"2\"\n";
    ss << "  DivisionByZeroMsg \"2\"\n";
    ss << "  StiffnessDetectionMsg \"2\"\n";
    ss << "  MaxConsecutiveZCs \"10000\"\n";
    ss << "  AlgebraicLoopWithStateMachineMsg \"2\"\n";
    ss << "  AssertionAction \"1\"\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_footer() {
    std::ostringstream ss;
    ss << "  InitialState  \"1\"\n";
    ss << "  SystemState   \"\"\n";
    ss << "  TaskingMode   \"1\"\n";
    ss << "  TaskConfigurations \"\"\n";
    ss << "  CodeGenParameterInlining \"2\"\n";
    ss << "  CodeGenFloatingPointFormat \"2\"\n";
    ss << "  CodeGenAbsTimeUsageMsg \"3\"\n";
    ss << "  CodeGenBaseName \"\"\n";
    ss << "  CodeGenOutputDir \"\"\n";
    ss << "  CodeGenExtraOpts \"\"\n";
    ss << "  CodeGenTarget \"Generic\"\n";
    ss << "  CodeGenTargetSettings \"\"\n";
    ss << "  ExtendedMatrixPrecision \"1\"\n";
    ss << "  MatrixSignificanceCheck \"2\"\n";
    ss << "  EnableStateSpaceSplitting \"2\"\n";
    ss << "  DisplayStateSpaceSplitting \"1\"\n";
    ss << "  DiscretizationMethod \"2\"\n";
    ss << "  ExternalModeSettings \"\"\n";
    ss << "  AlgebraicLoopMethod \"1\"\n";
    ss << "  AlgebraicLoopTolerance \"1e-6\"\n";
    ss << "  ScriptsDialogGeometry \"\"\n";
    ss << "  ScriptsDialogSplitterPos \"0\"\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_schematic_header() {
    std::ostringstream ss;
    ss << "  Schematic {\n";
    ss << "    Location      [400, 200; 1300, 900]\n";
    ss << "    ZoomFactor    1\n";
    ss << "    SliderPosition [0, 0]\n";
    ss << "    ShowBrowser   off\n";
    ss << "    BrowserWidth  100\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_schematic_footer() {
    return "  }\n}\n";
}

std::string CircuitSimulatorExporterPlecsModel::emit_magnetic_interface(
    const std::string& name, const std::string& turnsVar, int polarity,
    std::vector<int> position, const std::string& direction, bool flipped) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          MagneticInterface\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          off\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     " << direction << "\n";
    ss << "      Flipped       " << (flipped ? "on" : "off") << "\n";
    ss << "      LabelPosition east\n";
    ss << "      Parameter {\n";
    ss << "        Variable      \"n\"\n";
    ss << "        Value         \"" << turnsVar << "\"\n";
    ss << "        Show          off\n";
    ss << "      }\n";
    ss << "      Parameter {\n";
    ss << "        Variable      \"Polarity\"\n";
    ss << "        Value         \"" << polarity << "\"\n";
    ss << "        Show          off\n";
    ss << "      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_p_sat(
    const std::string& name, const std::string& areaVar, const std::string& lengthVar,
    const std::string& muRVar, const std::string& bSatVar,
    std::vector<int> position, const std::string& direction, bool flipped) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          Reference\n";
    ss << "      SrcComponent  \"Components/Magnetic/Components/P_sat\"\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          off\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     " << direction << "\n";
    ss << "      Flipped       " << (flipped ? "on" : "off") << "\n";
    ss << "      LabelPosition east\n";
    ss << "      Frame         [-8, -15; 8, 15]\n";
    ss << "      Parameter {\n        Variable      \"fitting\"\n        Value         \"1\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"A\"\n        Value         \"" << areaVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"l\"\n        Value         \"" << lengthVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"mu_r_unsat\"\n        Value         \"" << muRVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"mu_r_sat\"\n        Value         \"1\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"B_sat\"\n        Value         \"" << bSatVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"F_init\"\n        Value         \"0\"\n        Show          off\n      }\n";
    ss << "      Terminal {\n        Type          MagneticPort\n        Position      [0, -20]\n        Direction     up\n      }\n";
    ss << "      Terminal {\n        Type          MagneticPort\n        Position      [0, 20]\n        Direction     down\n      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_p_air(
    const std::string& name, const std::string& areaVar, const std::string& lengthVar,
    std::vector<int> position, const std::string& direction, bool flipped) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          Reference\n";
    ss << "      SrcComponent  \"Components/Magnetic/Components/P_air\"\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          off\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     " << direction << "\n";
    ss << "      Flipped       " << (flipped ? "on" : "off") << "\n";
    ss << "      LabelPosition east\n";
    ss << "      Frame         [-8, -10; 8, 10]\n";
    ss << "      Parameter {\n        Variable      \"A\"\n        Value         \"" << areaVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"l\"\n        Value         \"" << lengthVar << "\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"F_init\"\n        Value         \"0\"\n        Show          off\n      }\n";
    ss << "      Terminal {\n        Type          MagneticPort\n        Position      [0, -15]\n        Direction     up\n      }\n";
    ss << "      Terminal {\n        Type          MagneticPort\n        Position      [0, 15]\n        Direction     down\n      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_ac_voltage_source(
    const std::string& name, std::vector<int> position) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          ACVoltageSource\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          on\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     down\n";
    ss << "      Flipped       on\n";
    ss << "      LabelPosition east\n";
    ss << "      Parameter {\n        Variable      \"V\"\n        Value         \"V_p\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"w\"\n        Value         \"2*pi*f_sin\"\n        Show          off\n      }\n";
    ss << "      Parameter {\n        Variable      \"phi\"\n        Value         \"0\"\n        Show          off\n      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_resistor(
    const std::string& name, const std::string& valueVar, std::vector<int> position) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          Resistor\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          on\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     up\n";
    ss << "      Flipped       off\n";
    ss << "      LabelPosition east\n";
    ss << "      Parameter {\n        Variable      \"R\"\n        Value         \"" << valueVar << "\"\n        Show          off\n      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_scope(
    const std::string& name, std::vector<int> position) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          Scope\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          on\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     up\n";
    ss << "      Flipped       off\n";
    ss << "      LabelPosition south\n";
    ss << "      Location      [588, 348; 938, 595]\n";
    ss << "      State         \"\"\n";
    ss << "      SavedViews    \"\"\n";
    ss << "      HeaderState   \"\"\n";
    ss << "      PlotPalettes  \"\"\n";
    ss << "      Axes          \"1\"\n";
    ss << "      TimeRange     \"0.0\"\n";
    ss << "      ScrollingMode \"1\"\n";
    ss << "      SingleTimeAxis \"1\"\n";
    ss << "      Open          \"0\"\n";
    ss << "      Ts            \"-1\"\n";
    ss << "      SampleLimit   \"0\"\n";
    ss << "      XAxisLabel    \"Time / s\"\n";
    ss << "      ShowLegend    \"1\"\n";
    ss << "      Axis {\n";
    ss << "        Name          \"\"\n";
    ss << "        AutoScale     1\n";
    ss << "        MinValue      0\n";
    ss << "        MaxValue      1\n";
    ss << "        Signals       {}\n";
    ss << "        SignalTypes   [ ]\n";
    ss << "        Untangle      0\n";
    ss << "        KeepBaseline  off\n";
    ss << "        BaselineValue 0\n";
    ss << "      }\n";
    ss << "      Fourier {\n";
    ss << "        SingleXAxis       on\n";
    ss << "        AxisLabel         \"Frequency / Hz\"\n";
    ss << "        Scaling           0\n";
    ss << "        PhaseDisplay      0\n";
    ss << "        ShowFourierLegend off\n";
    ss << "        Axis {\n";
    ss << "          Name          \"\"\n";
    ss << "          AutoScale     1\n";
    ss << "          MinValue      0\n";
    ss << "          MaxValue      1\n";
    ss << "          Signals       {}\n";
    ss << "          Untangle      0\n";
    ss << "          KeepBaseline  off\n";
    ss << "          BaselineValue 0\n";
    ss << "        }\n";
    ss << "      }\n";
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_probe(
    const std::string& name, std::vector<int> position,
    const std::vector<std::pair<std::string, std::string>>& probes) {
    std::ostringstream ss;
    ss << "    Component {\n";
    ss << "      Type          PlecsProbe\n";
    ss << "      Name          \"" << name << "\"\n";
    ss << "      Show          on\n";
    ss << "      Position      [" << position[0] << ", " << position[1] << "]\n";
    ss << "      Direction     right\n";
    ss << "      Flipped       off\n";
    ss << "      LabelPosition south\n";
    for (const auto& [component, signal] : probes) {
        ss << "      Probe {\n";
        ss << "        Component     \"" << component << "\"\n";
        ss << "        Path          \"\"\n";
        ss << "        Signals       {\"" << signal << "\"}\n";
        ss << "      }\n";
    }
    ss << "    }\n";
    return ss.str();
}

std::string CircuitSimulatorExporterPlecsModel::emit_connection(
    const std::string& type,
    const std::string& srcComponent, int srcTerminal,
    const std::string& dstComponent, int dstTerminal,
    std::vector<std::vector<int>> points) {
    std::ostringstream ss;
    ss << "    Connection {\n";
    ss << "      Type          " << type << "\n";
    ss << "      SrcComponent  \"" << srcComponent << "\"\n";
    ss << "      SrcTerminal   " << srcTerminal << "\n";
    if (!points.empty()) {
        ss << "      Points        ";
        for (size_t i = 0; i < points.size(); ++i) {
            ss << "[" << points[i][0] << ", " << points[i][1] << "]";
            if (i < points.size() - 1) ss << "; ";
        }
        ss << "\n";
    }
    ss << "      DstComponent  \"" << dstComponent << "\"\n";
    ss << "      DstTerminal   " << dstTerminal << "\n";
    ss << "    }\n";
    return ss.str();
}

void CircuitSimulatorExporterPlecsModel::build_physical_init(
    std::ostringstream& init, Core& core,
    const std::vector<Winding>& windings,
    const std::vector<ColumnElement>& columns,
    double temperature) {

    size_t numWindings = windings.size();
    init << "mu_0 = 4*pi*1e-7;\n\n";

    // Winding turns
    for (size_t i = 0; i < numWindings; ++i) {
        std::string varName = (numWindings == 1) ? "n_turns" : "n_w" + std::to_string(i);
        init << varName << " = " << windings[i].get_number_turns() << ";\n";
    }
    init << "\n";

    // Core parameters
    double mu_r = core.get_initial_permeability(temperature);
    init << "mu_r = " << mu_r << ";\n";
    double bSat = core.get_magnetic_flux_density_saturation(temperature, false);
    init << "B_sat = " << bSat << ";\n\n";

    bool isMultiColumn = columns.size() > 1;
    if (isMultiColumn) {
        for (size_t ci = 0; ci < columns.size(); ++ci) {
            auto& col = columns[ci];
            std::string suffix;
            if (col.get_coordinates()[0] == 0) suffix = "_center";
            else if (col.get_coordinates()[0] < 0) suffix = "_left";
            else suffix = "_right";
            init << "A" << suffix << " = " << col.get_area() << ";\n";
            init << "l" << suffix << " = " << col.get_height() << ";\n";
        }
        double lPlate = core.get_effective_length();
        for (auto& col : columns) lPlate -= col.get_height();
        lPlate = std::max(lPlate / 2.0, 1e-6);
        init << "l_plate = " << lPlate << ";\n";
        init << "l_outer_total = l_left + 2*l_plate;\n";
    } else {
        init << "A_e = " << core.get_effective_area() << ";\n";
        init << "l_e = " << core.get_effective_length() << ";\n";
    }
    init << "\n";

    // Air gap parameters per column
    for (size_t ci = 0; ci < columns.size(); ++ci) {
        auto gapsInColumn = core.find_gaps_by_column(columns[ci]);
        for (size_t gi = 0; gi < gapsInColumn.size(); ++gi) {
            if (gapsInColumn[gi].get_length() <= 0) continue;
            std::string gapId = "gap_c" + std::to_string(ci) + "_g" + std::to_string(gi);
            init << "l_" << gapId << " = " << gapsInColumn[gi].get_length() << ";\n";
            double gapArea = gapsInColumn[gi].get_area().has_value()
                ? gapsInColumn[gi].get_area().value()
                : core.get_effective_area();
            init << "A_" << gapId << " = " << gapArea << ";\n";
        }
    }
    init << "\n";
}

void CircuitSimulatorExporterPlecsModel::build_magnetic_schematic(
    std::ostringstream& init, std::ostringstream& schematic,
    Core& core,
    const std::vector<Winding>& windings,
    const std::vector<ColumnElement>& columns,
    bool isMultiColumn, int yBase) {

    const int magX = 440;
    const int pSatX = 405;
    const int ySpacing = 60;
    size_t numWindings = windings.size();

    if (!isMultiColumn) {
        // ============ SINGLE-COLUMN TOPOLOGY ============
        int primaryIndex = 0;
        int secondaryIndex = 0;
        for (size_t i = 0; i < numWindings; ++i) {
            std::string miName = "MagInt_w" + std::to_string(i);
            std::string turnsVar = (numWindings == 1) ? "n_turns" : "n_w" + std::to_string(i);
            if (windings[i].get_isolation_side() == IsolationSide::PRIMARY) {
                int yPos = yBase + primaryIndex * ySpacing;
                schematic << emit_magnetic_interface(miName, turnsVar, 1, {magX, yPos}, "up", false);
                ++primaryIndex;
            } else {
                int yPos = yBase - (secondaryIndex + 1) * ySpacing;
                schematic << emit_magnetic_interface(miName, turnsVar, 1, {magX, yPos}, "down", true);
                ++secondaryIndex;
            }
        }

        schematic << emit_p_sat("P_satCore", "A_e", "l_e", "mu_r", "B_sat",
                                {pSatX, yBase - 50}, "up", false);

        auto gapsInColumn = core.find_gaps_by_column(columns[0]);
        int gapY = yBase + static_cast<int>(numWindings) * ySpacing + 10;
        int gapCount = 0;
        for (size_t gi = 0; gi < gapsInColumn.size(); ++gi) {
            if (gapsInColumn[gi].get_length() <= 0) continue;
            std::string gapId = "gap_c0_g" + std::to_string(gi);
            schematic << emit_p_air("P_air_" + gapId, "A_" + gapId, "l_" + gapId,
                                    {magX, gapY + gapCount * 40}, "up", true);
            gapCount++;
        }
        if (gapCount == 0) {
            init << "l_air_equiv = l_e;\nA_air_equiv = A_e;\n";
            schematic << emit_p_air("P_airEquiv", "A_air_equiv", "l_air_equiv",
                                    {magX, gapY}, "up", true);
            gapCount = 1;
        }

        // Magnetic connections
        for (size_t i = 0; i + 1 < numWindings; ++i) {
            schematic << emit_connection("Magnetic", "MagInt_w" + std::to_string(i), 3,
                                         "MagInt_w" + std::to_string(i + 1), 4);
        }
        std::string topMI = (numWindings > 1) ? "MagInt_w" + std::to_string(numWindings - 1) : "MagInt_w0";
        int topMagY = yBase - 80;
        schematic << emit_connection("Magnetic", "P_satCore", 1, topMI, 3,
                                     {{pSatX, topMagY}, {magX, topMagY}});
        std::string lastGap = (gapCount > 0 && !gapsInColumn.empty() && gapsInColumn[0].get_length() > 0)
            ? "P_air_gap_c0_g" + std::to_string(gapsInColumn.size() - 1)
            : "P_airEquiv";
        schematic << emit_connection("Magnetic", lastGap, 2, "MagInt_w0", 4);
        std::string firstGap = (!gapsInColumn.empty() && gapsInColumn[0].get_length() > 0)
            ? "P_air_gap_c0_g0" : "P_airEquiv";
        int bottomY = gapY + (gapCount - 1) * 40 + 30;
        schematic << emit_connection("Magnetic", firstGap, 1, "P_satCore", 2,
                                     {{magX, bottomY}, {pSatX, bottomY}});
    } else {
        // ============ MULTI-COLUMN (E-CORE) TOPOLOGY ============
        int outerX_L = 310, outerX_R = 570;

        size_t centerColIdx = 0;
        for (size_t ci = 0; ci < columns.size(); ++ci) {
            if (columns[ci].get_coordinates()[0] == 0) { centerColIdx = ci; break; }
        }

        int primaryIndex = 0;
        int secondaryIndex = 0;
        for (size_t i = 0; i < numWindings; ++i) {
            std::string miName = "MagInt_w" + std::to_string(i);
            std::string turnsVar = "n_w" + std::to_string(i);
            if (windings[i].get_isolation_side() == IsolationSide::PRIMARY) {
                int yPos = yBase + primaryIndex * ySpacing;
                schematic << emit_magnetic_interface(miName, turnsVar, 1, {magX, yPos}, "up", false);
                ++primaryIndex;
            } else {
                int yPos = yBase - (secondaryIndex + 1) * ySpacing;
                schematic << emit_magnetic_interface(miName, turnsVar, 1, {magX, yPos}, "down", true);
                ++secondaryIndex;
            }
        }

        schematic << emit_p_sat("P_satCenter", "A_center", "l_center", "mu_r", "B_sat",
                                {magX, yBase + static_cast<int>(numWindings) * ySpacing + 20}, "up", false);

        auto centerGaps = core.find_gaps_by_column(columns[centerColIdx]);
        int centerGapCount = 0;
        int centerGapY = yBase - 20;
        for (size_t gi = 0; gi < centerGaps.size(); ++gi) {
            if (centerGaps[gi].get_length() <= 0) continue;
            std::string gapId = "gap_c" + std::to_string(centerColIdx) + "_g" + std::to_string(gi);
            schematic << emit_p_air("P_air_" + gapId, "A_" + gapId, "l_" + gapId,
                                    {magX, centerGapY - centerGapCount * 40}, "up", true);
            centerGapCount++;
        }
        if (centerGapCount == 0) {
            init << "l_no_gap = 1e-6;\nA_no_gap = A_center;\n";
            schematic << emit_p_air("P_air_no_gap", "A_no_gap", "l_no_gap",
                                    {magX, centerGapY}, "up", true);
            centerGapCount = 1;
        }

        schematic << emit_p_sat("P_satOuterL", "A_left", "l_outer_total", "mu_r", "B_sat",
                                {outerX_L, yBase + 50}, "up", false);
        schematic << emit_p_sat("P_satOuterR", "A_right", "l_outer_total", "mu_r", "B_sat",
                                {outerX_R, yBase + 50}, "up", false);

        // Determine first/last gap names for connections
        std::string firstCenterGap, lastCenterGap;
        if (centerGapCount > 0 && !centerGaps.empty() && centerGaps[0].get_length() > 0) {
            firstCenterGap = "P_air_gap_c" + std::to_string(centerColIdx) + "_g0";
            size_t lastGapIdx = 0;
            for (size_t gi = 0; gi < centerGaps.size(); ++gi) {
                if (centerGaps[gi].get_length() > 0) lastGapIdx = gi;
            }
            lastCenterGap = "P_air_gap_c" + std::to_string(centerColIdx) + "_g" + std::to_string(lastGapIdx);
        } else {
            firstCenterGap = "P_air_no_gap";
            lastCenterGap = "P_air_no_gap";
        }

        // Chain gaps in series if multiple
        for (int g = 0; g + 1 < centerGapCount; ++g) {
            std::string gapA = "P_air_gap_c" + std::to_string(centerColIdx) + "_g" + std::to_string(g);
            std::string gapB = "P_air_gap_c" + std::to_string(centerColIdx) + "_g" + std::to_string(g + 1);
            schematic << emit_connection("Magnetic", gapA, 1, gapB, 2);
        }

        // Magnetic connections
        for (size_t i = 0; i + 1 < numWindings; ++i) {
            schematic << emit_connection("Magnetic", "MagInt_w" + std::to_string(i), 3,
                                         "MagInt_w" + std::to_string(i + 1), 4);
        }
        std::string topMI = (numWindings > 1) ? "MagInt_w" + std::to_string(numWindings - 1) : "MagInt_w0";
        schematic << emit_connection("Magnetic", lastCenterGap, 2, topMI, 3);
        schematic << emit_connection("Magnetic", "P_satCenter", 1, "MagInt_w0", 4);
        int topNodeY = yBase - 50;
        schematic << emit_connection("Magnetic", firstCenterGap, 1, "P_satOuterL", 1,
                                     {{magX, topNodeY}, {outerX_L, topNodeY}});
        schematic << emit_connection("Magnetic", firstCenterGap, 1, "P_satOuterR", 1,
                                     {{magX, topNodeY}, {outerX_R, topNodeY}});
        int bottomNodeY = yBase + static_cast<int>(numWindings) * ySpacing + 60;
        schematic << emit_connection("Magnetic", "P_satCenter", 2, "P_satOuterL", 2,
                                     {{magX, bottomNodeY}, {outerX_L, bottomNodeY}});
        schematic << emit_connection("Magnetic", "P_satCenter", 2, "P_satOuterR", 2,
                                     {{magX, bottomNodeY}, {outerX_R, bottomNodeY}});
    }
}

void CircuitSimulatorExporterPlecsModel::build_electrical_schematic(
    std::ostringstream& schematic,
    const std::vector<Winding>& windings,
    bool isMultiColumn, int yBase) {

    const int magX = 440;
    const int vSrcX = 205;
    const int loadX = 635;
    const int scopeX = 720;
    const int probeX = 660;
    const int ySpacing = 60;
    size_t numWindings = windings.size();

    schematic << emit_ac_voltage_source("V_ac", {vSrcX, yBase});

    // Load resistor per secondary winding
    int secCount = 0;
    for (size_t i = 0; i < numWindings; ++i) {
        if (windings[i].get_isolation_side() == IsolationSide::PRIMARY) continue;
        std::string rName = (secCount == 0) ? "R_load" : "R_load_" + std::to_string(secCount);
        int secY = yBase - static_cast<int>(i) * ySpacing + 30;
        schematic << emit_resistor(rName, "RL", {loadX, secY});
        secCount++;
    }

    // Scope + Probe
    int scopeY = isMultiColumn ? yBase - 160 : yBase - 120;
    schematic << emit_scope("Scope", {scopeX, scopeY});
    std::vector<std::pair<std::string, std::string>> probeSignals = {{"V_ac", "Source voltage"}};
    secCount = 0;
    for (size_t i = 0; i < numWindings; ++i) {
        if (windings[i].get_isolation_side() == IsolationSide::PRIMARY) continue;
        std::string rName = (secCount == 0) ? "R_load" : "R_load_" + std::to_string(secCount);
        probeSignals.push_back({rName, "Resistor voltage " + std::to_string(secCount)});
        secCount++;
    }
    schematic << emit_probe("Probe", {probeX, scopeY}, probeSignals);

    // Wire connections
    secCount = 0;
    for (size_t i = 0; i < numWindings; ++i) {
        std::string miName = "MagInt_w" + std::to_string(i);
        if (windings[i].get_isolation_side() == IsolationSide::PRIMARY) {
            schematic << emit_connection("Wire", miName, 1, "V_ac", 1,
                                         {{magX - 75, yBase - 15}, {vSrcX, yBase - 15}});
            schematic << emit_connection("Wire", miName, 2, "V_ac", 2,
                                         {{magX - 75, yBase + 15}, {vSrcX, yBase + 15}});
        } else {
            std::string rName = (secCount == 0) ? "R_load" : "R_load_" + std::to_string(secCount);
            int secY = yBase - static_cast<int>(i) * ySpacing;
            schematic << emit_connection("Wire", miName, 1, rName, 1,
                                         {{magX + 65, secY - 15}, {loadX, secY - 15}});
            schematic << emit_connection("Wire", miName, 2, rName, 2,
                                         {{magX + 65, secY + 15}, {loadX, secY + 15}});
            secCount++;
        }
    }

    schematic << emit_connection("Signal", "Probe", 1, "Scope", 1);
}

std::string CircuitSimulatorExporterPlecsModel::assemble_plecs_file(
    const std::string& name, const std::string& initStr, const std::string& schematicStr) {
    std::string result;
    result += emit_header(name);
    result += encode_init_commands(initStr);
    result += emit_footer();
    result += schematicStr;
    return result;
}

std::string CircuitSimulatorExporterPlecsModel::export_magnetic_as_subcircuit(
    Magnetic magnetic, double frequency, double temperature,
    [[maybe_unused]] std::optional<std::string> filePathOrFile,
    [[maybe_unused]] CircuitSimulatorExporterCurveFittingModes mode) {

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    auto columns = core.get_columns();
    auto windings = coil.get_functional_description();
    bool isMultiColumn = columns.size() > 1;

    std::string name = magnetic.get_reference();
    if (name.empty()) name = "magnetic";

    // Build InitializationCommands
    std::ostringstream init;
    init << std::scientific << std::setprecision(_precision);
    init << "%% PLECS model generated by OpenMagnetics\n";
    init << "%% Component: " << name << "\n\n";
    init << "%% Simulation\n";
    init << "tsim = 0.01;\n\n";
    init << "%% Electrical\n";
    init << "V_p = 1;\n";
    init << "f_sin = " << frequency << ";\n";
    init << "RL = 1e3;\n\n";

    build_physical_init(init, core, windings, columns, temperature);

    // DC resistance per winding
    for (size_t i = 0; i < windings.size(); ++i) {
        double dcR = WindingLosses::calculate_effective_resistance_of_winding(magnetic, i, 0.1, temperature);
        init << "R_dc_w" << i << " = " << dcR << ";\n";
    }

    // Build schematic
    std::ostringstream schematic;
    schematic << emit_schematic_header();

    int yBase = isMultiColumn ? 390 : 290;
    build_magnetic_schematic(init, schematic, core, windings, columns, isMultiColumn, yBase);
    build_electrical_schematic(schematic, windings, isMultiColumn, yBase);

    schematic << emit_schematic_footer();

    return assemble_plecs_file(name, init.str(), schematic.str());
}

std::string CircuitSimulatorExporterPlecsModel::export_magnetic_as_symbol(
    Magnetic magnetic, [[maybe_unused]] std::optional<std::string> filePathOrFile) {

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    auto columns = core.get_columns();
    auto windings = coil.get_functional_description();
    bool isMultiColumn = columns.size() > 1;

    std::string name = magnetic.get_reference();
    if (name.empty()) name = "magnetic";

    // Build InitializationCommands (no simulation/electrical params)
    std::ostringstream init;
    init << std::scientific << std::setprecision(_precision);
    init << "%% PLECS subcircuit symbol generated by OpenMagnetics\n";
    init << "%% Component: " << name << "\n\n";

    build_physical_init(init, core, windings, columns, defaults.ambientTemperature);

    // Build schematic — magnetic circuit only
    std::ostringstream schematic;
    schematic << emit_schematic_header();

    int yBase = isMultiColumn ? 390 : 290;
    build_magnetic_schematic(init, schematic, core, windings, columns, isMultiColumn, yBase);

    schematic << emit_schematic_footer();

    return assemble_plecs_file(name, init.str(), schematic.str());
}

} // namespace OpenMagnetics
