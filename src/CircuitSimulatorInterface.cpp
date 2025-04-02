#include "Utils.h"
#include "Defaults.h"
#include "CircuitSimulatorInterface.h"
#include "WindingOhmicLosses.h"
#include "LeakageInductance.h"
#include "MagnetizingInductance.h"
#include "WindingLosses.h"
#include "Sweeper.h"
#include "Settings.h"
#include <filesystem>
#include <ctime>
#include "levmar.h"
#include <float.h>

namespace OpenMagnetics {

template<typename T> std::vector<T> convolution_valid(std::vector<T> const &f, std::vector<T> const &g) {
    int const nf = f.size();
    int const ng = g.size();
    std::vector<T> const &min_v = (nf < ng)? f : g;
    std::vector<T> const &max_v = (nf < ng)? g : f;
    int const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
    std::vector<T> out(n, T());
    for(auto i(0); i < n; ++i) {
        for(int j(min_v.size() - 1), k(i); j >= 0; --j) {
            out[i] += min_v[j] * max_v[k];
            ++k;
        }
    }
    return out;
}

CircuitSimulatorExporter::CircuitSimulatorExporter(CircuitSimulatorExporterModels program) {
    _model = CircuitSimulatorExporterModel::factory(program);
}

CircuitSimulatorExporter::CircuitSimulatorExporter() {
    _model = CircuitSimulatorExporterModel::factory(CircuitSimulatorExporterModels::SIMBA);
}

double CircuitSimulatorExporter::analytical_model(double x[], double frequency) {
    return x[0] + x[1] * sqrt(frequency) + x[2] * frequency; 
}

void CircuitSimulatorExporter::analytical_func(double *p, double *x, int m, int n, void *data) {
    double* frequencies = static_cast <double*> (data);
    for(int i=0; i<n; ++i) {
        x[i]=analytical_model(p, frequencies[i]);
    }
}

std::complex<double> parallel(std::complex<double> z0, std::complex<double> z1) {
    return 1.0 / (1.0 / z0 + 1.0 / z1);
}

double CircuitSimulatorExporter::ladder_model(double x[], double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    auto R0 = std::complex<double>(dcResistance);
    auto R1 = std::complex<double>(x[0], 0);
    auto L1 = std::complex<double>(0, w * x[1]);
    auto R2 = std::complex<double>(x[2], 0);
    auto L2 = std::complex<double>(0, w * x[3]);
    auto R3 = std::complex<double>(x[4], 0);
    auto L3 = std::complex<double>(0, w * x[5]);
    auto R4 = std::complex<double>(x[6], 0);
    auto L4 = std::complex<double>(0, w * x[7]);
    auto R5 = std::complex<double>(x[8], 0);
    auto L5 = std::complex<double>(0, w * x[9]);

    for(int i=0; i<10; ++i) {
        if (x[i] < 0) {
            return 0;
        }
    }

    return (R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3 + parallel(L4, R4 + parallel(L5, R5)))))).real();
}

void CircuitSimulatorExporter::ladder_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast <double*> (data);
    double dcResistance = dcResistanceAndFrequencies[0];

    for(int i=0; i<n; ++i)
        x[i]=ladder_model(p, dcResistanceAndFrequencies[i + 1], dcResistance);
}


double CircuitSimulatorExporter::core_ladder_model(double x[], double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    auto R0 = std::complex<double>(dcResistance);
    auto R1 = std::complex<double>(x[0], 0);
    auto L1 = std::complex<double>(0, w * x[1]);
    auto R2 = std::complex<double>(x[2], 0);
    auto L2 = std::complex<double>(0, w * x[3]);
    auto R3 = std::complex<double>(x[4], 0);
    auto L3 = std::complex<double>(0, w * x[5]);

    for(int i=0; i<6; ++i) {
        if (x[i] < 0) {
            return 0;
        }
    }

    return (R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3)))).real();
}

void CircuitSimulatorExporter::core_ladder_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast <double*> (data);
    double dcResistance = dcResistanceAndFrequencies[0];

    for(int i=0; i<n; ++i)
        x[i]=core_ladder_model(p, dcResistanceAndFrequencies[i + 1], dcResistance);
}

std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_ladder(MagneticWrapper magnetic) {
    const size_t numberUnknowns = 10;

    const size_t numberElements = 100;
    const size_t numberElementsPlusOne = 101;
    size_t loopIterations = 5;
    size_t maxIterations = 10000;
    double startingFrequency = 0.1;
    double endingFrequency = 1000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex);
        auto frequenciesVector = windingAcResistanceData.get_x_points();

        auto valuePoints = windingAcResistanceData.get_y_points();
        double acResistances[numberElements];
        for (size_t index = 0; index < valuePoints.size(); ++index) {
            acResistances[index] = valuePoints[index];
        }

        double bestError = DBL_MAX;
        double initialState = 10;
        std::vector<double> bestAcResistanceCoefficientsThisWinding;
        for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
            double coefficients[numberUnknowns];
            for (size_t index = 0; index < numberUnknowns; ++index) {
                coefficients[index] = initialState;
            }
            double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

            double lmInitMu = 1e-03;
            double lmStopThresh = 1e-25;
            double lmDiffDelta = 1e-19;

            opts[0]= lmInitMu;
            opts[1]= lmStopThresh;
            opts[2]= lmStopThresh;
            opts[3]= lmStopThresh;
            opts[4]= lmDiffDelta;

            double dcResistanceAndfrequencies[numberElementsPlusOne];
            dcResistanceAndfrequencies[0] = acResistances[0];
            // coefficients[1] = 1e-7;
            for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                dcResistanceAndfrequencies[index + 1] = frequenciesVector[index];
            }

            dlevmar_dif(CircuitSimulatorExporter::ladder_func, coefficients, acResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&dcResistanceAndfrequencies));

            double errorAverage = 0;
            for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                double modeledAcResistance = CircuitSimulatorExporter::ladder_model(coefficients, frequenciesVector[index], acResistances[0]);
                double error = fabs(valuePoints[index] - modeledAcResistance) / valuePoints[index];
                errorAverage += error;
            }

            errorAverage /= frequenciesVector.size();
            initialState /= 10;

            if (errorAverage < bestError) {
                bestError = errorAverage;
                bestAcResistanceCoefficientsThisWinding.clear();
                for (auto coefficient : coefficients) {
                    bestAcResistanceCoefficientsThisWinding.push_back(coefficient);
                }
            }
        }
        acResistanceCoefficientsPerWinding.push_back(bestAcResistanceCoefficientsThisWinding);

    }
    return acResistanceCoefficientsPerWinding;
}

std::vector<double> CircuitSimulatorExporter::calculate_core_resistance_coefficients(MagneticWrapper magnetic) {
    const size_t numberUnknowns = 6;

    const size_t numberElements = 100;
    const size_t numberElementsPlusOne = 101;
    size_t loopIterations = 5;
    size_t maxIterations = 10000;
    double startingFrequency = 0.1;
    double endingFrequency = 1000000;
    auto coil = magnetic.get_coil();

    std::vector<double> coreResistanceCoefficients;


    Curve2D coreResistanceData = Sweeper().sweep_core_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, 25);
    auto frequenciesVector = coreResistanceData.get_x_points();

    auto valuePoints = coreResistanceData.get_y_points();
    double coreResistances[numberElements];
    for (size_t index = 0; index < valuePoints.size(); ++index) {
        coreResistances[index] = valuePoints[index];
    }

    double bestError = DBL_MAX;
    double initialState = 10;
    std::vector<double> bestCoreResistanceCoefficients;
    for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
        double coefficients[numberUnknowns];
        for (size_t index = 0; index < numberUnknowns; ++index) {
            coefficients[index] = initialState;
        }
        double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

        double lmInitMu = 1e-03;
        double lmStopThresh = 1e-25;
        double lmDiffDelta = 1e-19;

        opts[0]= lmInitMu;
        opts[1]= lmStopThresh;
        opts[2]= lmStopThresh;
        opts[3]= lmStopThresh;
        opts[4]= lmDiffDelta;

        double dcResistanceAndfrequencies[numberElementsPlusOne];
        dcResistanceAndfrequencies[0] = coreResistances[0];
        // coefficients[1] = 1e-7;
        for (size_t index = 0; index < frequenciesVector.size(); ++index) {
            dcResistanceAndfrequencies[index + 1] = frequenciesVector[index];
        }

        dlevmar_dif(CircuitSimulatorExporter::core_ladder_func, coefficients, coreResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&dcResistanceAndfrequencies));

        double errorAverage = 0;
        for (size_t index = 0; index < frequenciesVector.size(); ++index) {
            double modeledAcResistance = CircuitSimulatorExporter::ladder_model(coefficients, frequenciesVector[index], coreResistances[0]);
            double error = fabs(valuePoints[index] - modeledAcResistance) / valuePoints[index];
            errorAverage += error;
        }

        errorAverage /= frequenciesVector.size();
        initialState /= 10;

        if (errorAverage < bestError) {
            bestError = errorAverage;
            bestCoreResistanceCoefficients.clear();
            for (auto coefficient : coefficients) {
                bestCoreResistanceCoefficients.push_back(coefficient);
            }
        }
    }
    return bestCoreResistanceCoefficients;
}


std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_analytical(MagneticWrapper magnetic) {
    const size_t numberUnknowns = 4;
    const size_t numberElements = 100;

    size_t maxIterations = 10000;
    double startingFrequency = 0.1;
    double endingFrequency = 1000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex);
        auto frequenciesVector = windingAcResistanceData.get_x_points();

        auto points = windingAcResistanceData.get_y_points();
        double acResistances[numberElements];
        for (size_t index = 0; index < points.size(); ++index) {
            acResistances[index] = points[index];
        }

        double coefficients[numberUnknowns];
        for (size_t index = 0; index < numberUnknowns; ++index) {
            coefficients[index] = 1;
        }
        double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

        double lmInitMu = 1e-03;
        double lmStopThresh = 1e-25;
        double lmDiffDelta = 1e-19;

        opts[0]= lmInitMu;
        opts[1]= lmStopThresh;
        opts[2]= lmStopThresh;
        opts[3]= lmStopThresh;
        opts[4]= lmDiffDelta;

        double frequencies[numberElements];
        for (size_t index = 0; index < frequenciesVector.size(); ++index) {
            frequencies[index] = frequenciesVector[index];
        }

        dlevmar_dif(CircuitSimulatorExporter::analytical_func, coefficients, acResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&frequencies));

        acResistanceCoefficientsPerWinding.push_back(std::vector<double>());
        for (auto coefficient : coefficients) {
            acResistanceCoefficientsPerWinding.back().push_back(coefficient);
        }

    }
    return acResistanceCoefficientsPerWinding;
}


std::vector<std::vector<double>> CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(MagneticWrapper magnetic, CircuitSimulatorExporterCurveFittingModes mode) {
    if (mode == CircuitSimulatorExporterCurveFittingModes::LADDER) {
        return calculate_ac_resistance_coefficients_per_winding_ladder(magnetic);
    }
    else {
        return calculate_ac_resistance_coefficients_per_winding_analytical(magnetic);
    }
}

CircuitSimulatorExporterSimbaModel::CircuitSimulatorExporterSimbaModel() {
    std::random_device rd;
    _gen = std::mt19937(rd());

    std::seed_seq sseq{time(0)};
    
    _gen.seed(sseq);
}

std::shared_ptr<CircuitSimulatorExporterModel> CircuitSimulatorExporterModel::factory(CircuitSimulatorExporterModels programName){
    if (programName == CircuitSimulatorExporterModels::SIMBA) {
        return std::make_shared<CircuitSimulatorExporterSimbaModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::NGSPICE) {
        return std::make_shared<CircuitSimulatorExporterNgspiceModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::LTSPICE) {
        return std::make_shared<CircuitSimulatorExporterLtspiceModel>();
    }
    else
        throw std::runtime_error("Unknown Circuit Simulator program, available options are: {SIMBA, NGSPICE, LTSPICE}");
}

std::string CircuitSimulatorExporter::export_magnetic_as_subcircuit(MagneticWrapper magnetic, double frequency, std::optional<std::string> outputFilename, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    auto result = _model->export_magnetic_as_subcircuit(magnetic, frequency, filePathOrFile, mode);
    if (outputFilename) {
        std::ofstream o(outputFilename.value());
        o << std::setw(2) << result << std::endl;
    }
    return result;
}

std::string CircuitSimulatorExporter::export_magnetic_as_symbol(MagneticWrapper magnetic, std::optional<std::string> outputFilename, std::optional<std::string> filePathOrFile) {
    auto result = _model->export_magnetic_as_symbol(magnetic, filePathOrFile);
    if (outputFilename) {
        std::ofstream o(outputFilename.value());
        o << std::setw(2) << result << std::endl;
    }
    return result;
}

std::string CircuitSimulatorExporterSimbaModel::generate_id() {
    // generator for hex numbers from 0 to F
    std::uniform_int_distribution<int> dis(0x0, 0xF);
    
    const char hexChars[16] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
    std::string id = "";
    for (std::size_t i = 0; i < 8; ++i) {
        id += hexChars[dis(_gen)];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[dis(_gen)];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[dis(_gen)];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[dis(_gen)];
    }
    id += "-";
    for (std::size_t i = 0; i < 12; ++i) {
        id += hexChars[dis(_gen)];
    }
    return id;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_device(std::string libraryName, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson;

    deviceJson["LibraryName"] = libraryName;
    deviceJson["Top"] = coordinates[1];
    deviceJson["Left"] = coordinates[0];
    deviceJson["Angle"] = angle;
    deviceJson["HF"] = false;
    deviceJson["VF"] = false;
    deviceJson["Disabled"] = false;
    deviceJson["Name"] = name;
    deviceJson["ID"] = generate_id();
    deviceJson["Parameters"]["Name"] = name;
    deviceJson["EnabledScopes"] = ordered_json::array();

    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_air_gap(std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Air Gap", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = "1";
    deviceJson["Parameters"]["Area"] = std::to_string(area);
    deviceJson["Parameters"]["Length"] = std::to_string(length);
    return deviceJson;
}
ordered_json CircuitSimulatorExporterSimbaModel::create_core(double initialPermeability, std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Linear Core", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = std::to_string(initialPermeability);
    deviceJson["Parameters"]["Area"] = std::to_string(area);
    deviceJson["Parameters"]["Length"] = std::to_string(length);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_winding(size_t numberTurns, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson = create_device("Winding", coordinates, angle, name);

    deviceJson["Parameters"]["NumberOfTurns"] = std::to_string(numberTurns);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_resistor(double resistance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Resistor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = std::to_string(resistance);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_inductor(double inductance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Inductor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = std::to_string(inductance);
    deviceJson["Parameters"]["Iinit"] = "0";
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_pin(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Electrical Pin", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_magnetic_ground(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Magnetic Ground", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_connector(std::vector<int> startingCoordinates, std::vector<int> endingCoordinates, std::string name) {
    ordered_json connectorJson;

    connectorJson["Name"] = name;
    connectorJson["Segments"] = ordered_json::array();

    if (startingCoordinates[0] == endingCoordinates[0] || startingCoordinates[1] == endingCoordinates[1]) {
        ordered_json segment;
        segment["StartX"] = startingCoordinates[0];
        segment["StartY"] = startingCoordinates[1];
        segment["EndX"] = endingCoordinates[0];
        segment["EndY"] = endingCoordinates[1];
        connectorJson["Segments"].push_back(segment);
    }
    else {
        {
            ordered_json segment;
            segment["StartX"] = startingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = startingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
        {
            ordered_json segment;
            segment["StartX"] = endingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = endingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
    }

    return connectorJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::merge_connectors(ordered_json connectors) {
    auto mergeConnectors = connectors;
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t firstIndex = 0; firstIndex < mergeConnectors.size(); ++firstIndex) {
            bool sharedPoint = false;
            ordered_json mergedConnector;
            for (size_t secondIndex = firstIndex + 1; secondIndex < mergeConnectors.size(); ++secondIndex) {
                for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        if ((firstConnectorSegment["StartX"] == secondConnectorSegment["StartX"] &&
                            firstConnectorSegment["StartY"] == secondConnectorSegment["StartY"]) ||
                            (firstConnectorSegment["EndX"] == secondConnectorSegment["EndX"] &&
                            firstConnectorSegment["EndY"] == secondConnectorSegment["EndY"])) {

                            sharedPoint = true;
                            break;
                        }
                    }
                    if (sharedPoint) {
                        break;
                    }
                }

                if (sharedPoint) {
                    mergedConnector["Segments"] = ordered_json::array();
                    mergedConnector["Name"] = "Merge of connector: " + std::string{mergeConnectors[firstIndex]["Name"]} + " with " + std::string{mergeConnectors[secondIndex]["Name"]};
                    for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(firstConnectorSegment);
                    }
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(secondConnectorSegment);
                    }
                    mergeConnectors.erase(secondIndex);
                    mergeConnectors.erase(firstIndex);
                    mergeConnectors.push_back(mergedConnector);
                    merged = true;
                    break;
                }
            }

            if (sharedPoint) {
                break;
            }
        }
    }
    return mergeConnectors;
}

std::string CircuitSimulatorExporterSimbaModel::export_magnetic_as_subcircuit(MagneticWrapper magnetic, double frequency, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    ordered_json simulation;
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    _scale = _modelSize / core.get_width();


    if (filePathOrFile) {
        try {
            if(!std::filesystem::exists(filePathOrFile.value())) {
                throw std::runtime_error("File not found");
            }
            std::ifstream json_file(filePathOrFile.value());
            if(json_file.is_open()) {
                simulation = ordered_json::parse(json_file);
            }
            else {
                throw std::runtime_error("File not found");
            }
        }
        catch(const std::runtime_error& re) {
            std::stringstream json_file(filePathOrFile.value());
            simulation = ordered_json::parse(json_file);
        }
    }

    if (!simulation.contains("Libraries")) {
        simulation["Libraries"] = ordered_json::array();
    }
    if (!simulation.contains("Designs")) {
        simulation["Designs"] = ordered_json::array();
    }
    ordered_json library;
    ordered_json device;
    library["LibraryItemName"] = "OpenMagnetics components";
    // library["Devices"] = ordered_json::array();

    device["LibraryName"] = magnetic.get_reference();
    device["Angle"] = 0;
    device["Disabled"] = false;
    device["Name"] = magnetic.get_reference();
    device["Id"] = generate_id();
    device["Parameters"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Devices"] = ordered_json::array();
    device["SubcircuitDefinition"]["Connectors"] = ordered_json::array();
    device["SubcircuitDefinition"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Id"] = generate_id();
    device["SubcircuitDefinition"]["Variables"] = ordered_json::array();
    device["SubcircuitDefinition"]["VariableFile"] = "";
    device["SubcircuitDefinitionID"] = device["SubcircuitDefinition"]["Id"];

    auto columns = core.get_columns();
    auto coreEffectiveArea = core.get_effective_area();

    auto coreEffectiveLengthMinusColumns = core.get_effective_length();
    if (columns.size() > 1) {
        for (auto column : columns) {
            if (column.get_coordinates()[0] >= 0) {
                coreEffectiveLengthMinusColumns -= column.get_height();
            }
        }
    }

    std::vector<std::vector<int>> columnBottomCoordinates; 
    std::vector<std::vector<int>> columnTopCoordinates; 
    
    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        auto column = columns[columnIndex];
        auto gapsInThisColumn = core.find_gaps_by_column(column);
        ordered_json coreChunkJson;
        std::vector<int> columnCoordinates;
        if (column.get_coordinates()[0] == 0 && column.get_coordinates()[2] != 0) {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[2] * _scale), 0};
        }
        else {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[0] * _scale), 0};
        }
        if (columnIndex == 0) {
            coreChunkJson = create_core(core.get_initial_permeability(), columnCoordinates, coreEffectiveArea, coreEffectiveLengthMinusColumns, 90, "Core winding column and plates " + std::to_string(columnIndex));
        }
        else {
            coreChunkJson = create_core(core.get_initial_permeability(), columnCoordinates, coreEffectiveArea, column.get_height(), 90, "Core lateral column " + std::to_string(columnIndex));
        }
        std::vector<int> columnTopCoordinate = {columnCoordinates[0] + 3, -2};  // Don't ask
        std::vector<int> columnBottomCoordinate = {columnCoordinates[0] + 3, 4};  // Don't ask
        device["SubcircuitDefinition"]["Devices"].push_back(coreChunkJson);
        int currentColumnGapHeight = 6;
        for (size_t gapIndex = 0; gapIndex < gapsInThisColumn.size(); ++gapIndex) {
            auto gap = gapsInThisColumn[gapIndex];
            if (!gap.get_coordinates()) {
                throw std::runtime_error("Gap is not processed");
            }
            std::vector<int> gapCoordinates = {columnCoordinates[0], currentColumnGapHeight};

            auto gapJson = create_air_gap(gapCoordinates, gap.get_area().value(), gap.get_length(), 90, "Column " + std::to_string(columnIndex) + " gap " + std::to_string(gapIndex));
            device["SubcircuitDefinition"]["Devices"].push_back(gapJson);
            currentColumnGapHeight += 6;
            columnBottomCoordinate[1] += 6;
        }
        columnBottomCoordinates.push_back(columnBottomCoordinate);
        columnTopCoordinates.push_back(columnTopCoordinate);
    }

    double leakageInductance = resolve_dimensional_values(LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0]);

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, windingIndex, frequency, Defaults().ambientTemperature);
        auto winding = coil.get_functional_description()[windingIndex];
        std::vector<int> coordinates = {columnTopCoordinates[0][0] - 2, columnTopCoordinates[0][1] - 6};
        ordered_json windingJson;
        ordered_json topPinJson;
        ordered_json bottomPinJson;
        ordered_json acResistorJson;
        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 0, winding.get_name());
        }
        else {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 180, winding.get_name());
        }

        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            coordinates[0] -= 6;
            acResistorJson = create_resistor(effectiveResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");
            if (windingIndex == 0) {
                coordinates[0] -= 6;
                auto leakageInductanceJson = create_inductor(leakageInductance, coordinates, 0, winding.get_name() + " Leakage inductance");
                device["SubcircuitDefinition"]["Devices"].push_back(leakageInductanceJson);
            }
            coordinates[0] -= 2;
            bottomPinJson = create_pin(coordinates, 0, winding.get_name() + " Input");
            if (windingIndex == 0) {
                coordinates[0] += 12;
            }
            else {
                coordinates[0] += 6;
            }
            coordinates[1] += 4;
            topPinJson = create_pin(coordinates, 0, winding.get_name() + " Output");
        }
        else {
            coordinates[0] += 4;
            acResistorJson = create_resistor(effectiveResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");
            coordinates[0] += 6;
            topPinJson = create_pin(coordinates, 180, winding.get_name() + " Input");
            coordinates[1] += 4;
            coordinates[0] -= 6;
            bottomPinJson = create_pin(coordinates, 180, winding.get_name() + " Output");
        }

        auto connectorTopCoordinates = columnTopCoordinates[0];
        auto connectorBottomCoordinates = columnTopCoordinates[0];
        if (windingIndex == 0) {
            connectorTopCoordinates[1] -= 1;
        }
        else {
            connectorTopCoordinates[1] -= 1;
            connectorBottomCoordinates[1] += 1;
        }
        auto connectorJson = create_connector(connectorBottomCoordinates, connectorTopCoordinates, "Connector between winding " + std::to_string(windingIndex) + " and winding " + std::to_string(windingIndex + 1));
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);

        columnTopCoordinates[0][1] -= 6;
        device["SubcircuitDefinition"]["Devices"].push_back(windingJson);
        device["SubcircuitDefinition"]["Devices"].push_back(topPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(bottomPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(acResistorJson);

        if (windingIndex == coil.get_functional_description().size() - 1) {
            std::vector<int> finalConnectorTopCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[0][1] - 5};
            std::vector<int> finalConnectorBottomCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[0][1] + 1};
            auto connectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding " + std::to_string(windingIndex) + " and top");
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            columnTopCoordinates[0][1] -= 5;

        }
    }

    // Magnetic ground
    {
        std::vector<int> columnBottomCoordinatesAux = {0, columnTopCoordinates[0][1]};
        columnBottomCoordinatesAux[0] += 2;
        columnBottomCoordinatesAux[1] -= 2;
        auto groundJson = create_magnetic_ground(columnBottomCoordinatesAux, 180, "Magnetic ground");
        device["SubcircuitDefinition"]["Devices"].push_back(groundJson);
    }

    for (size_t columnIndex = 1; columnIndex < columns.size(); ++columnIndex) {
        auto connectorJson = create_connector(columnTopCoordinates[0], columnTopCoordinates[columnIndex], "Top Connector between column " + std::to_string(0) + " and columm " + std::to_string(columnIndex));
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
    }
    device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);


    if (columns.size() == 1) {
        auto columnBottomCoordinatesAux = columnBottomCoordinates[0];
        columnBottomCoordinatesAux[1] = 0;
        columnBottomCoordinatesAux[0] += _modelSize / 2;
        auto bottomConnectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinatesAux, "Bottom Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(bottomConnectorJson);
        auto topConnectorJson = create_connector(columnTopCoordinates[0], columnBottomCoordinatesAux, "Rop Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(topConnectorJson);
        device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);
    }
    else {
        for (size_t columnIndex = 1; columnIndex < columns.size(); ++columnIndex) {
            auto connectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinates[columnIndex], "Bottom Connector between column " + std::to_string(0) + " and columm " + std::to_string(columnIndex));
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }
    }

    library["Devices"] = ordered_json::array();
    library["Devices"].push_back(device);
    simulation["Libraries"].push_back(library);
    return simulation.dump(2);
}

std::string CircuitSimulatorExporterNgspiceModel::export_magnetic_as_subcircuit(MagneticWrapper magnetic, double frequency, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    parametersString += ".param MagnetizingInductance_Value=" + std::to_string(magnetizingInductance) + "\n";
    parametersString += ".param Permeance=MagnetizingInductance_Value/NumberTurns_1**2\n";
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, Defaults().ambientTemperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index - 1]);
            double couplingCoefficient = sqrt((magnetizingInductance - leakageInductance) / magnetizingInductance);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient);
        }

        std::vector<std::string> c;
        for (auto coeff : acResistanceCoefficientsPerWinding[index]){

            std::ostringstream out;
            out.precision(20);
            out << std::fixed << coeff;
            c.push_back(std::move(out).str());
        }
        if (mode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            throw std::invalid_argument("Analytica mode not supported in NgSpice");
        }
        else {
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                std::string ladderIndexs = std::to_string(ladderIndex);
                circuitString += "Lladder" + is + "_" + ladderIndexs + " P" + is + "+ Node_Lladder_" + is + "_" + ladderIndexs + " " + c[ladderIndex + 1] + "\n";
                if (ladderIndex == 0) {
                    circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_R_Lmag_" + is + " " + c[ladderIndex] + "\n";
                }
                else {
                    circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + std::to_string(ladderIndex - 2) + " " + c[ladderIndex] + "\n";
                }
            }
            circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            circuitString += "Lmag_" + is + " P" + is + "- Node_R_Lmag_" + is + " {NumberTurns_" + is + "**2*Permeance}\n";

        }
        if (index > 0) {
            circuitString += "K Lmag_1 Lmag_" + is + " {CouplingCoefficient_1" + is + "_Value}\n";
        }

        headerString += " P" + is + "+ P" + is + "-";

    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_subcircuit(MagneticWrapper magnetic, double frequency, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    parametersString += ".param MagnetizingInductance_Value=" + std::to_string(magnetizingInductance) + "\n";
    parametersString += ".param Permeance=MagnetizingInductance_Value/NumberTurns_1**2\n";
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, Defaults().ambientTemperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index - 1]);
            double couplingCoefficient = sqrt((magnetizingInductance - leakageInductance) / magnetizingInductance);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient);
        }

        std::vector<std::string> c;
        for (auto coeff : acResistanceCoefficientsPerWinding[index]){

            std::ostringstream out;
            out.precision(20);
            out << std::fixed << coeff;
            c.push_back(std::move(out).str());
        }
        if (mode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            circuitString += "E" + is + " P" + is + "+ Node_R_Lmag_" + is + " P" + is + "+ Node_R_Lmag_" + is + " Laplace = 1 /(" + c[0] + " + " + c[1] + " * sqrt(abs(s)/(2*pi)) + " + c[2] + " * abs(s)/(2*pi))\n";
            circuitString += "Lmag_" + is + " P" + is + "- Node_R_Lmag_" + is + " {NumberTurns_" + is + "**2*Permeance}\n";
        }
        else {
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                std::string ladderIndexs = std::to_string(ladderIndex);
                circuitString += "Lladder" + is + "_" + ladderIndexs + " P" + is + "+ Node_Lladder_" + is + "_" + ladderIndexs + " " + c[ladderIndex + 1] + "\n";
                if (ladderIndex == 0) {
                    circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_R_Lmag_" + is + " " + c[ladderIndex] + "\n";
                }
                else {
                    circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + std::to_string(ladderIndex - 2) + " " + c[ladderIndex] + "\n";
                }
            }
            circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            circuitString += "Lmag_" + is + " P" + is + "- Node_R_Lmag_" + is + " {NumberTurns_" + is + "**2*Permeance}\n";

        }
        if (index > 0) {
            circuitString += "K Lmag_1 Lmag_" + is + " {CouplingCoefficient_1" + is + "_Value}\n";
        }

        headerString += " P" + is + "+ P" + is + "-";

    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_symbol(MagneticWrapper magnetic, std::optional<std::string> filePathOrFile) {
    std::string symbolString = "Version 4\n";
    symbolString += "SymbolType BLOCK\n";

    auto coil = magnetic.get_coil();

    int rectangleSemiWidth = 72;

    int leftSideSize = 16;
    int rightSideSize = 16;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            leftSideSize += 64;
        }
        else {
            rightSideSize += 64;
        }
    }

    int rectangleHeight = std::max(leftSideSize, rightSideSize);


    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(-rectangleHeight / 2 + 8) + " Left 0 " + magnetic.get_reference() + "\n";
    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(rectangleHeight / 2 - 8) + " Left 0 Made with OpenMagnetics\n";


    symbolString += "RECTANGLE Normal " + std::to_string(-rectangleSemiWidth) + " -" + std::to_string(rectangleHeight / 2) + " " + std::to_string(rectangleSemiWidth) + " " + std::to_string(rectangleHeight / 2) + "\n";
    symbolString += "SYMATTR Prefix X\n";
    symbolString += "SYMATTR Value " + fix_filename(magnetic.get_reference()) + "\n";
    symbolString += "SYMATTR ModelFile " + fix_filename(magnetic.get_reference()) + ".cir\n";

    int currentSpiceOrder = 1;
    int currentRectangleLeftSideHeight = -leftSideSize / 2 + 24;
    int currentRectangleRightSideHeight = -rightSideSize / 2 + 24;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else{
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "+\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else {
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "-\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;
    }


    return symbolString;
}

void CircuitSimulationReader::process_line(std::string line, char separator) {
    std::stringstream ss(line);
    std::string token;

    if (_columns.size() == 0) {
        // Getting column names
        while(getline(ss, token, separator)) {
            if (int(token[0]) < 32) {
                continue;
            }
            CircuitSimulationSignal circuitSimulationSignal;
            token.erase(std::remove(token.begin(), token.end(), (char)13), token.end());
            token.erase(std::remove(token.begin(), token.end(), '\"'), token.end());
            circuitSimulationSignal.name = token;
            _columns.push_back(circuitSimulationSignal);
        }
    }
    else {
        size_t currentColumnIndex = 0;
        while(getline(ss, token, separator)) {
            if (int(token[0]) < 32) {
                continue;
            }
            _columns[currentColumnIndex].data.push_back(stod(token));
            currentColumnIndex++;
        }
    }
}

CircuitSimulationReader::CircuitSimulationReader(std::string filePathOrFile, bool forceFile) {
    char separator = '\0';
    std::string line;

    std::filesystem::path path = filePathOrFile;

    if (path.has_parent_path() && !forceFile) {
        if (!std::filesystem::exists(filePathOrFile)) {
            throw std::runtime_error("File not found");
        }
        std::ifstream is(filePathOrFile);
        if (is.is_open()) {
            while(getline(is, line)) {
                if (separator == '\0') {
                    separator = guess_separator(line);
                }

                process_line(line, separator);
            }
            is.close();
        }
        else {
            throw std::runtime_error("File not found");
        }
    }
    else {
        std::stringstream ss(filePathOrFile);
        while(std::getline(ss, line, '\n')){
            if (separator == '\0') {
                separator = guess_separator(line);
            }

            process_line(line, separator);
        }
    }

    _time = find_time(_columns);
}

std::vector<std::string> CircuitSimulationReader::extract_column_names() {
    std::vector<std::string> names;
    for (auto column : _columns) {
        names.push_back(column.name);
    }

    return names;
}

bool CircuitSimulationReader::can_be_time(std::vector<double> data) {
    if (data.size() == 0) {
        throw std::invalid_argument("vector data cannot be empty");
    }
    if (data.size() == 1) {
        return false;
    }

    for (size_t index = 1; index < data.size(); index++) {
        if (data[index - 1] >= data[index]) {
            return false;
        }
    }

    return true;
}

bool CircuitSimulationReader::can_be_voltage(std::vector<double> data, double limit) {
    if (data.size() == 0) {
        throw std::invalid_argument("vector data cannot be empty");
    }
    if (data.size() == 1) {
        return false;
    }

    std::vector<double> distinctValues;
    std::vector<size_t> distinctValuesCount;
    for (size_t index = 0; index < data.size(); index++) {
        bool isDistinct = true;
        for (size_t distincIndex = 0; distincIndex < distinctValues.size(); distincIndex++) {
            double distinctValue = distinctValues[distincIndex];
            double error = fabs(distinctValue - data[index]) / std::max(fabs(data[index]), fabs(distinctValue));
            if (std::isnan(error)) {
                error = fabs(distinctValue - data[index]);
            }
            if (error <= limit) {
                isDistinct = false;
                distinctValuesCount[distincIndex]++;
            }
        }
        if (isDistinct) {
            distinctValues.push_back(data[index]);
            distinctValuesCount.push_back(1);
        }
    }



    size_t significantDistinctValues = 0;
    for (size_t distincIndex = 0; distincIndex < distinctValuesCount.size(); distincIndex++) {
        if (distinctValuesCount[distincIndex] > data.size() * limit) {
            significantDistinctValues++;
        }
    }

    bool result = significantDistinctValues == 2 || significantDistinctValues == 3;
    return result;
}

bool CircuitSimulationReader::can_be_current(std::vector<double> data, double limit) {
    std::vector<double> diffValues;
    for (size_t index = 0; index < data.size(); index++) {
        double diff;
        if (index == 0) {
            diff = data[index] - data[data.size() - 1];
        }
        else {
            diff = data[index] - data[index - 1];
        }
        diffValues.push_back(diff);
    }

    return can_be_voltage(diffValues, limit);
}

char CircuitSimulationReader::guess_separator(std::string line){
    std::vector<char> possibleSeparators = {',', ';', '\t'};

    for (auto separator : possibleSeparators) {
        std::stringstream ss(line);
        std::string token;
        size_t numberColumns = 0;
        while(getline(ss, token, separator)) {
            numberColumns++;
        }

        if (numberColumns >= 2 && numberColumns <= 30) {
            return separator;
        }
    }

    throw std::runtime_error("No column separator found");
}

Waveform CircuitSimulationReader::get_one_period(Waveform waveform, double frequency, bool sample) {
    double period = 1.0 / frequency;
    if (!waveform.get_time()) {
        throw std::runtime_error("Missing time data");
    }

    auto time = waveform.get_time().value();
    auto data = waveform.get_data();

    double periodEnd = time.back();
    double periodStart = periodEnd - period;
    int periodStartIndex = 0;
    int periodStopIndex = -1;

    if (_periodStartIndex && _periodStopIndex) {
        periodStartIndex = _periodStartIndex.value();
        periodStopIndex = _periodStopIndex.value();
    }
    else {
        for (int i = time.size() - 1; i >= 0; --i)
        {
            if (time[i] <= periodStart) {
                periodStartIndex = i;
                break;
            }
        }

        double previousData = data[periodStartIndex];
        for (int i = periodStartIndex - 1; i >= 0; --i)
        {
            if ((data[i] >= 0 && previousData < 0) || (data[i] <= 0 && previousData > 0)) {
                periodStartIndex = i;
                periodStart = time[i];
                break;
            }
            previousData = data[i];
        }


        for (int i = periodStartIndex; i < time.size(); ++i)
        {
            if (time[i] >= periodStart + period) {
                periodStopIndex = i + 1;
                break;
            }
        }
        _periodStartIndex = periodStartIndex;
        _periodStopIndex = periodStopIndex;
    }

    auto periodData = std::vector<double>(data.begin() + periodStartIndex, data.begin() + periodStopIndex);
    auto periodTime = std::vector<double>(time.begin() + periodStartIndex, time.begin() + periodStopIndex);


    double offset = periodTime[0];
    for (size_t index = 0; index < periodTime.size(); index++) {
        periodTime[index] -= offset;
    }


    Waveform newWaveform;
    newWaveform.set_data(periodData);
    newWaveform.set_time(periodTime);

    if (sample) {
        auto sampledWaveform = InputsWrapper::calculate_sampled_waveform(newWaveform, frequency);
        return sampledWaveform;
    }
    else {
        return newWaveform;
    }
}

CircuitSimulationReader::CircuitSimulationSignal CircuitSimulationReader::find_time(std::vector<CircuitSimulationReader::CircuitSimulationSignal> columns) {
    for (auto column : columns) {
        if (can_be_time(column.data)) {
            column.type = DataType::TIME;
            return column;
        }
    }
    throw std::runtime_error("no time column found");
}

Waveform CircuitSimulationReader::extract_waveform(CircuitSimulationReader::CircuitSimulationSignal signal, double frequency, bool sample) {
    auto settings = Settings::GetInstance();
    Waveform waveform;
    waveform.set_data(signal.data);
    waveform.set_time(_time.data);

    auto waveformOnePeriod = get_one_period(waveform, frequency, sample);
    Waveform reconstructedWaveform = waveformOnePeriod;

    if (false) {
        double originalThreshold = settings->get_harmonic_amplitude_threshold();
        while (reconstructedWaveform.get_data().size() > 8192) {
            settings->set_harmonic_amplitude_threshold(settings->get_harmonic_amplitude_threshold() * 2);
            auto harmonics = InputsWrapper::calculate_harmonics_data(waveformOnePeriod, frequency);
            settings->set_inputs_number_points_sampled_waveforms(2 * OpenMagnetics::round_up_size_to_power_of_2(harmonics.get_frequencies().back() / frequency));
            reconstructedWaveform = InputsWrapper::reconstruct_signal(harmonics, frequency);
        }

        settings->set_harmonic_amplitude_threshold(originalThreshold);
        return reconstructedWaveform;
    }
    else {
        return waveformOnePeriod;
    }
}

std::vector<int> get_numbers_in_string(std::string s) {
    std::regex regex(R"(\d+)");   // matches a sequence of digits
    std::vector<int> numbers;
    std::smatch match;
    while (std::regex_search(s, match, regex)) {
        numbers.push_back(std::stoi(match.str()));
        s = match.suffix();
    }

    return numbers;
}

bool CircuitSimulationReader::extract_winding_indexes(size_t numberWindings) {
    bool result = true;
    size_t numberFoundIndexes = 0;
    std::vector<size_t> indexes;
    std::map<std::string, size_t> windingLabels = {
        {"pri", 0},
        {"sec", 1},
        {"aux", 2},
        {"ter", 2},
        {"a", 0},
        {"b", 1},
        {"c", 2},
        {"HV", 0},
        {"LV", 1},
    };

    std::vector<CircuitSimulationSignal> columnsWithIndexes;
    std::vector<CircuitSimulationSignal> columnsWithResetIndexes;

    for (auto column : _columns) {
        column.windingIndex = -1;
        if (!can_be_time(column.data)) {
            auto numbersInColumnName = get_numbers_in_string(column.name);
            if (numbersInColumnName.size() > 0) {
                numberFoundIndexes++;
                indexes.push_back(numbersInColumnName.back());
                column.windingIndex = numbersInColumnName.back();
            }
            else {
                bool found = false;
                for (auto [label, windingIndex] : windingLabels) {
                    if (column.name.find(label) != std::string::npos)  {
                        numberFoundIndexes++;
                        indexes.push_back(windingIndex);
                        column.windingIndex = windingIndex;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    column.windingIndex = 0;
                    indexes.push_back(0);
                }
            }
        }
        columnsWithIndexes.push_back(column);
    }

    std::sort(indexes.begin(),indexes.end());
    indexes.resize(distance(indexes.begin(),std::unique(indexes.begin(),indexes.end())));

    for (size_t index = 0; index < indexes.size(); index++) {
        for (auto column : columnsWithIndexes) {
            if (column.windingIndex == indexes[index]) {
                column.windingIndex = index;
                columnsWithResetIndexes.push_back(column);
            }
        }
    }

    _columns = columnsWithResetIndexes;
    return numberFoundIndexes >= numberWindings;
}

std::vector<double> rolling_window_filter(std::vector<double> data) {

    size_t rollingFactorDividend = 192;
    size_t rollingFactor = data.size() / rollingFactorDividend;

    if (rollingFactor < 1){
        rollingFactor = 1;
    }

    for (size_t i = 0; i < rollingFactor - 1; ++i) {
        data.push_back(data[i]);
    }
    std::vector<double> ones(rollingFactor, 1.0);

    auto cleanData = convolution_valid(data, ones);
    for (size_t i = 0; i < cleanData.size(); ++i) {
        cleanData[i] /= rollingFactor;
    }
    return cleanData;
}

std::optional<CircuitSimulationReader::DataType> CircuitSimulationReader::guess_type_by_name(std::string name) {
    for (auto timeAlias : _timeAliases) {
        if (name.find(timeAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::TIME;
        }
    }
    for (auto currentAlias : _currentAliases) {
        if (name.find(currentAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::CURRENT;
        }
    }
    for (auto magnetizingCurrentAlias : _magnetizingCurrentAliases) {
        if (name.find(magnetizingCurrentAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::MAGNETIZING_CURRENT;
        }
    }
    for (auto voltageAlias : _voltageAliases) {
        if (name.find(voltageAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::VOLTAGE;
        }
    }

    return std::nullopt;
}

bool CircuitSimulationReader::extract_column_types(double frequency) {
    bool result = true;
    std::vector<std::map<std::string, std::string>> columnNameToSignalPerWinding;
    std::vector<CircuitSimulationSignal> columnsWithTypes;

    for (auto column : _columns) {
        if (can_be_time(column.data)) {
            column.type = DataType::TIME;
        }
        else {

            auto waveform = extract_waveform(column, frequency, false);
            auto data = waveform.get_data();
            int timeout = 100;

            auto guessedType = guess_type_by_name(column.name);

            if (guessedType) {
                column.type = guessedType.value();
            }
            else {
                while (timeout > 0) {
                    if (can_be_current(data)) {
                        column.type = DataType::CURRENT;
                        break;
                    }
                    else if (can_be_voltage(data)) {
                        column.type = DataType::VOLTAGE;
                        break;
                    }
                    else {
                        data = rolling_window_filter(data);
                    }
                    timeout--;
                } 
                if (timeout == 0) {
                    column.type = DataType::UNKNOWN;
                    // throw std::runtime_error("Unknown type");
                }
            }
            columnsWithTypes.push_back(column);
        }
    }
    _columns = columnsWithTypes;
    return result;
}

std::vector<std::map<std::string, std::string>> CircuitSimulationReader::extract_map_column_names(size_t numberWindings, double frequency) {
    OperatingPoint operatingPoint;
    std::vector<std::map<std::string, std::string>> columnNameToSignalPerWinding;

    std::vector<OperatingPointExcitation> excitationsPerWinding;

    extract_winding_indexes(numberWindings);
    extract_column_types(frequency);
    for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
        std::map<std::string, std::string> columnNameToSignal;
        columnNameToSignal["time"] = _time.name;
        for (auto column : _columns) {
            if (column.windingIndex == windingIndex && column.type == DataType::CURRENT) {
                columnNameToSignal["current"] = column.name;
            }
            if (column.windingIndex == windingIndex && column.type == DataType::MAGNETIZING_CURRENT) {
                columnNameToSignal["magnetizingCurrent"] = column.name;
            }
            if (column.windingIndex == windingIndex && column.type == DataType::VOLTAGE) {
                columnNameToSignal["voltage"] = column.name;
            }
        }

        columnNameToSignalPerWinding.push_back(columnNameToSignal);
    }

    return columnNameToSignalPerWinding;
}

bool CircuitSimulationReader::assign_column_names(std::vector<std::map<std::string, std::string>> columnNames) {
    std::vector<CircuitSimulationSignal> assignedColumns;
    for (size_t windingIndex = 0; windingIndex < columnNames.size(); windingIndex++) {
        for (auto [columnType, columnName] : columnNames[windingIndex]) {

            for (auto column : _columns) {
                if (column.name == columnName && columnType == "current") {
                    column.type = DataType::CURRENT;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "magnetizingCurrent") {
                    column.type = DataType::MAGNETIZING_CURRENT;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "voltage") {
                    column.type = DataType::VOLTAGE;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "time") {
                    column.type = DataType::TIME;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
            }
        }
    }

    _columns = assignedColumns;

    return true;
}

OperatingPoint CircuitSimulationReader::extract_operating_point(size_t numberWindings, double frequency, std::optional<std::vector<std::map<std::string, std::string>>> mapColumnNames, double ambientTemperature) {
    OperatingPoint operatingPoint;

    std::vector<OperatingPointExcitation> excitationsPerWinding;
    if (!mapColumnNames) {
        extract_winding_indexes(numberWindings);
        extract_column_types(frequency);
    }
    else {
        assign_column_names(mapColumnNames.value());
    }

    for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
        OperatingPointExcitation excitation;
        for (auto column : _columns) {
            excitation.set_frequency(frequency);
            if (column.windingIndex == int(windingIndex) && column.type == DataType::CURRENT) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor current;
                current.set_waveform(waveform);
                excitation.set_current(current);
            }
            if (column.windingIndex == int(windingIndex) && column.type == DataType::MAGNETIZING_CURRENT) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor current;
                current.set_waveform(waveform);
                excitation.set_magnetizing_current(current);
            }
            if (column.windingIndex == int(windingIndex) && column.type == DataType::VOLTAGE) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor voltage;
                voltage.set_waveform(waveform);
                excitation.set_voltage(voltage);
            }
        }
        excitationsPerWinding.push_back(excitation);
    }

    operatingPoint.set_excitations_per_winding(excitationsPerWinding);
    OperatingConditions conditions;
    conditions.set_cooling(std::nullopt);    
    conditions.set_ambient_temperature(ambientTemperature);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}



} // namespace OpenMagnetics

