#pragma once
#include "Constants.h"
#include "InitialPermeability.h"

#include <CoreWrapper.h>
#include <InputsWrapper.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

enum class ReluctanceModels : int {
    ZHANG,
    PARTRIDGE,
    EFFECTIVE_AREA,
    EFFECTIVE_LENGTH,
    MUEHLETHALER,
    STENGLEIN,
    BALAKRISHNAN,
    CLASSIC
};

class ReluctanceModel {
  private:
    double _magneticFluxDensitySaturation = 0.4;
  protected:
  public:
    virtual std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo) = 0;
    double get_ungapped_core_reluctance(CoreWrapper core, OperatingPoint* operatingPoint = nullptr) {
        auto constants = Constants();
        OpenMagnetics::InitialPermeability initialPermeability;

        auto coreMaterial = core.get_functional_description().get_material();

        double initialPermeabilityValue;
        if (operatingPoint != nullptr) {
            double temperature =
                operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
            _magneticFluxDensitySaturation = core.get_magneticFluxDensitySaturation(temperature, true);
            auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
            initialPermeabilityValue = initialPermeability.get_initial_permeability(
                coreMaterial, &temperature, nullptr, &frequency);
        }
        else {
            initialPermeabilityValue =
                initialPermeability.get_initial_permeability(coreMaterial);
            _magneticFluxDensitySaturation = core.get_magneticFluxDensitySaturation(true);
        }
        double absolutePermeability = constants.vacuumPermeability * initialPermeabilityValue;
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
        double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();

        double reluctanceCore = effectiveLength / (absolutePermeability * effectiveArea);
        return reluctanceCore;
    }

    double get_ungapped_core_reluctance(CoreWrapper core, double initialPermeability) {
        auto constants = Constants();
        double absolutePermeability = constants.vacuumPermeability * initialPermeability;
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
        double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();

        double reluctanceCore = effectiveLength / (absolutePermeability * effectiveArea);
        return reluctanceCore;
    }
    double get_gap_maximum_storable_energy(CoreGap gapInfo, double fringingFactor) {
        auto constants = Constants();
        auto gapLength = gapInfo.get_length();
        auto gapArea = *(gapInfo.get_area());

        double energyStoredInGap = 0.5 / constants.vacuumPermeability * gapLength * gapArea * fringingFactor *
                                      pow(_magneticFluxDensitySaturation, 2);

        return energyStoredInGap;
    }
    static std::map<std::string, std::string> get_models_information() {
        std::map<std::string, std::string> information;
        information["Zhang"] =
            R"(Based on "Improved Calculation Method for Inductance Value of the Air-Gap Inductor" by Xinsheng Zhang.)";
        information["Muehlethaler"] =
            R"(Based on "A Novel Approach for 3D Air Gap Reluctance Calculations" by Jonas Mühlethaler.)";
        information["Partridge"] =
            R"(Based on the method described in page 8-11 from "Transformer and Inductor Design Handbook Fourth Edition" by Colonel Wm. T. McLyman.)";
        information["Effective Area"] =
            R"(Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.)";
        information["Effective Length"] =
            R"(Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.)";
        information["Stenglein"] =
            R"(Based on "The Reluctance of Large Air Gaps in Ferrite Cores" by Erika Stenglein.)";
        information["Balakrishnan"] =
            R"(Based on "Air-gap reluctance and inductance calculations for magnetic circuits using a Schwarz-Christoffel transformation" by A. Balakrishnan.)";
        information["Classic"] = R"(Based on the reluctance of a uniform magnetic circuit)";
        return information;
    }
    static std::map<std::string, double> get_models_errors() {
        // This are taken manually from running the tests. Yes, a pain in the ass
        // TODO: Automate it somehow
        std::map<std::string, double> errors;
        errors["Zhang"] = 0.115811;
        errors["Muehlethaler"] = 0.110996;
        errors["Partridge"] = 0.124488;
        errors["Effective Area"] = 0.175055;
        errors["Effective Length"] = 0.175055;
        errors["Stenglein"] = 0.143346;
        errors["Balakrishnan"] = 0.136754;
        errors["Classic"] = 0.283863;
        return errors;
    }
    static std::map<std::string, std::string> get_models_external_links() {
        std::map<std::string, std::string> externalLinks;
        externalLinks["Zhang"] = "https://ieeexplore.ieee.org/document/9332553";
        externalLinks["Muehlethaler"] = "https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/"
                                         "10_A_Novel_Approach_ECCEAsia2011_01.pdf";
        externalLinks["Partridge"] =
            "https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook";
        externalLinks["Effective Area"] =
            "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
        externalLinks["Effective Length"] =
            "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
        externalLinks["Stenglein"] = "https://ieeexplore.ieee.org/document/7695271/";
        externalLinks["Balakrishnan"] = "https://ieeexplore.ieee.org/document/602560";
        externalLinks["Classic"] = "https://en.wikipedia.org/wiki/Magnetic_reluctance";
        return externalLinks;
    }
    static std::map<std::string, std::string> get_models_internal_links() {
        std::map<std::string, std::string> internalLinks;
        internalLinks["Zhang"] = "";
        internalLinks["Muehlethaler"] = "/musings/10_gap_reluctance_and_muehlethaler_method";
        internalLinks["Partridge"] = "";
        internalLinks["Effective Area"] = "";
        internalLinks["Effective Length"] = "";
        internalLinks["Stenglein"] = "/musings/11_inductance_variables_and_stenglein_method";
        internalLinks["Balakrishnan"] = "";
        internalLinks["Classic"] = "";
        return internalLinks;
    }
    double get_core_reluctance(CoreWrapper core, OperatingPoint* operatingPoint = nullptr) {
        auto coreReluctance = get_ungapped_core_reluctance(core, operatingPoint);

        if (std::isnan(coreReluctance)) {
            throw std::runtime_error("Core Reluctance must be a number, not NaN");
        }
        double calculatedReluctance = coreReluctance + get_gapping_reluctance(core);
        if (std::isnan(calculatedReluctance)) {
            throw std::runtime_error("Reluctance must be a number, not NaN");
        }
        return calculatedReluctance;
    }

    double get_core_reluctance(CoreWrapper core, double initialPermeability) {
        auto coreReluctance = get_ungapped_core_reluctance(core, initialPermeability);

        double calculatedReluctance = coreReluctance + get_gapping_reluctance(core);

        return calculatedReluctance;
    }

    double get_gapping_reluctance(CoreWrapper core) {
        double calculatedReluctance = 0;
        double calculatedCentralReluctance = 0;
        double calculatedLateralReluctance = 0;
        auto gapping = core.get_functional_description().get_gapping();
        if (gapping.size() == 0) {
            return 0;
        }

        // We recompute all gaps in case some is missing coordinates
        for (const auto& gap : gapping) {
            if (!gap.get_coordinates()) {
                core.process_gap();
                gapping = core.get_functional_description().get_gapping();
                break;
            }
        }
        for (const auto& gap : gapping) {
            auto gapReluctance = get_gap_reluctance(gap);
            auto gapColumn = core.find_closest_column_by_coordinates(gap.get_coordinates().value());
            if (gapColumn.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                calculatedLateralReluctance += 1 / gapReluctance["reluctance"];
            }
            else {
                calculatedCentralReluctance += gapReluctance["reluctance"];
            }
            if (gapReluctance["fringing_factor"] < 1) {
                std::cout << "fringing_factor " << gapReluctance["fringing_factor"] << std::endl;
            }
        }
        calculatedReluctance = calculatedCentralReluctance + 1 / calculatedLateralReluctance;
        return calculatedReluctance;
    }
    ReluctanceModel() = default;
    virtual ~ReluctanceModel() = default;


    static std::shared_ptr<ReluctanceModel> factory(ReluctanceModels modelName);
};

// Based on Improved Calculation Method for Inductance Value of the Air-Gap Inductor by Xinsheng Zhang
// https://sci-hub.wf/https://ieeexplore.ieee.org/document/9332553
class ReluctanceZhangModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

// Based on A Novel Approach for 3D Air Gap Reluctance Calculations by Jonas Mühlethaler
// https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf
class ReluctanceMuehlethalerModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);

    double get_basic_reluctance(double l, double w, double h) {
        auto constants = Constants();
        return 1 / constants.vacuumPermeability /
               (w / 2 / l + 2 / std::numbers::pi * (1 + log(std::numbers::pi * h / 4 / l)));
    }

    double get_reluctance_type_1(double l, double w, double h) {
        double basicReluctance = get_basic_reluctance(l, w, h);
        return 1 / (1 / (basicReluctance + basicReluctance) + 1 / (basicReluctance + basicReluctance));
    }
};

class ReluctanceEffectiveAreaModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

class ReluctanceEffectiveLengthModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

class ReluctancePartridgeModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

class ReluctanceClassicModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

// Based on Air-gap reluctance and inductance calculations for magnetic circuits using a Schwarz-Christoffel
// transformation by A. Balakrishnan https://sci-hub.wf/https://ieeexplore.ieee.org/document/602560
class ReluctanceBalakrishnanModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
};

// Based on The Reluctance of Large Air Gaps in Ferrite Cores by Erika Stenglein
// https://sci-hub.wf/10.1109/EPE.2016.7695271
class ReluctanceStengleinModel : public ReluctanceModel {
  public:
    std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);

    double u(double rx, double l1) { return 42.7 * rx / l1 - 50.2; }

    double v(double rx, double l1) { return -55.4 * rx / l1 + 71.6; }

    double w(double rx, double l1) { return 0.88 * rx / l1 - 0.80; }

    double alpha(double rx, double l1, double lg) {
        return u(rx, l1) * pow(lg / l1, 2) + v(rx, l1) * lg / l1 + w(rx, l1);
    }
};

} // namespace OpenMagnetics