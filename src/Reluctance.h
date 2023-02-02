#pragma once
#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <filesystem>
#include <streambuf>
#include "Constants.h"
#include "InitialPermeability.h"
#include <magic_enum.hpp>

#include <InputsWrapper.h>
#include <CoreWrapper.h>

namespace OpenMagnetics {

    enum class ReluctanceModels : int { ZHANG, PARTRIDGE, EFFECTIVE_AREA, EFFECTIVE_LENGTH, MUEHLETHALER, STENGLEIN, BALAKRISHNAN, CLASSIC };

    class ReluctanceModel {
        private:
        protected:
        public:
            virtual std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo) = 0;
            double get_ungapped_core_reluctance(CoreWrapper core, OperationPoint* operationPoint=nullptr) {
                auto constants = Constants();
                OpenMagnetics::InitialPermeability initial_permeability;

                double initial_permeability_value;
                if (operationPoint != nullptr) {
                    double temperature = operationPoint->get_conditions().get_ambient_temperature();  // TODO: Use a future calculated temperature
                    auto frequency = operationPoint->get_excitations_per_winding()[0].get_frequency();
                    initial_permeability_value = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                               &temperature,
                                                                                               nullptr,
                                                                                               &frequency);
                }
                else {
                    initial_permeability_value = initial_permeability.get_initial_permeability(core.get_functional_description().get_material());
                }
                double absolute_permeability = constants.vacuum_permeability * initial_permeability_value;
                double effective_area = core.get_processed_description()->get_effective_parameters().get_effective_area();
                double effective_length = core.get_processed_description()->get_effective_parameters().get_effective_length();

                double reluctance_core = effective_length / (absolute_permeability * effective_area);
                return reluctance_core;
            }

            double get_ungapped_core_reluctance(CoreWrapper core, double initial_permeability) {
                auto constants = Constants();
                double absolute_permeability = constants.vacuum_permeability * initial_permeability;
                double effective_area = core.get_processed_description()->get_effective_parameters().get_effective_area();
                double effective_length = core.get_processed_description()->get_effective_parameters().get_effective_length();

                double reluctance_core = effective_length / (absolute_permeability * effective_area);
                return reluctance_core;
            }
            double get_gap_maximum_storable_energy(CoreGap gapInfo, double fringing_factor) {
                auto constants = Constants();
                double magnetic_flux_density_saturation = constants.magnetic_flux_density_saturation; // HARDCODED TODO: replace when materials are implemented
                auto gap_length = gapInfo.get_length();
                auto gap_area = *(gapInfo.get_area());

                double energy_stored_in_gap = 0.5 / constants.vacuum_permeability * gap_length * gap_area * fringing_factor * pow(magnetic_flux_density_saturation, 2);

                return energy_stored_in_gap;
            }
            static std::map<std::string, std::string> get_models_information() {
                std::map<std::string, std::string> information;
                information["Zhang"] = "Based on \"Improved Calculation Method for Inductance Value of the Air-Gap Inductor\" by Xinsheng Zhang.";
                information["Muehlethaler"] = "Based on \"A Novel Approach for 3D Air Gap Reluctance Calculations\" by Jonas Mühlethaler.";
                information["Partridge"] = "Based on the method described in page 8-11 from \"Transformer and Inductor Design Handbook Fourth Edition\" by Colonel Wm. T. McLyman.";
                information["Effective Area"] = "Based on the method described in page 60 from \"High-Frequency Magnetic Components, Second Edition\" by Marian Kazimierczuk.";
                information["Effective Length"] = "Based on the method described in page 60 from \"High-Frequency Magnetic Components, Second Edition\" by Marian Kazimierczuk.";
                information["Stenglein"] = "Based on \"The Reluctance of Large Air Gaps in Ferrite Cores\" by Erika Stenglein.";
                information["Balakrishnan"] = "Based on \"Air-gap reluctance and inductance calculations for magnetic circuits using a Schwarz-Christoffel transformation\" by A. Balakrishnan.";
                information["Classic"] = "Based on the reluctance of a uniform magnetic circuit";
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
            static std::map<std::string, std::string> get_models_internal_links() {
                std::map<std::string, std::string> external_links;
                external_links["Zhang"] = "https://ieeexplore.ieee.org/document/9332553";
                external_links["Muehlethaler"] = "https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf";
                external_links["Partridge"] = "https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook";
                external_links["Effective Area"] = "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
                external_links["Effective Length"] = "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
                external_links["Stenglein"] = "https://ieeexplore.ieee.org/document/7695271/";
                external_links["Balakrishnan"] = "https://ieeexplore.ieee.org/document/602560";
                external_links["Classic"] = "https://en.wikipedia.org/wiki/Magnetic_reluctance";
                return external_links;
            }
            static std::map<std::string, std::string> get_models_external_links() {
                std::map<std::string, std::string> internal_links;
                internal_links["Zhang"] = "";
                internal_links["Muehlethaler"] = "/musings/10";
                internal_links["Partridge"] = "";
                internal_links["Effective Area"] = "";
                internal_links["Effective Length"] = "";
                internal_links["Stenglein"] = "";
                internal_links["Balakrishnan"] = "";
                internal_links["Classic"] = "";
                return internal_links;
            }
            double get_core_reluctance(CoreWrapper core, OperationPoint* operationPoint=nullptr){
                auto coreReluctance = get_ungapped_core_reluctance(core, operationPoint);
                double calculatedReluctance = coreReluctance + get_gapping_reluctance(core);
                return calculatedReluctance;
            }

            double get_core_reluctance(CoreWrapper core, double initialPermeability){
                auto coreReluctance = get_ungapped_core_reluctance(core, initialPermeability);
                double calculatedReluctance = coreReluctance + get_gapping_reluctance(core);
                return calculatedReluctance;
            }

            double get_gapping_reluctance(CoreWrapper core){
                double calculatedReluctance = 0;
                double calculatedCentralReluctance = 0;
                double calculatedLateralReluctance = 0;
                auto gapping = core.get_functional_description().get_gapping();
                for (const auto& gap: gapping) {
                    auto gap_reluctance = get_gap_reluctance(gap);
                    auto gapColumn = core.find_closest_column_by_coordinates(*gap.get_coordinates());
                    if (gapColumn.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                        calculatedLateralReluctance += 1 / gap_reluctance["reluctance"];
                    }
                    else {
                        calculatedCentralReluctance += gap_reluctance["reluctance"];
                    }
                    if (gap_reluctance["fringing_factor"] < 1) {
                        std::cout << "fringing_factor " << gap_reluctance["fringing_factor"] << std::endl;
                    }

                }
                calculatedReluctance = calculatedCentralReluctance + 1 / calculatedLateralReluctance;
                return calculatedReluctance;
            }

            static std::shared_ptr<ReluctanceModel> factory(ReluctanceModels modelName);
    };

    // Based on Improved Calculation Method for Inductance Value of the Air-Gap Inductor by Xinsheng Zhang
    // https://sci-hub.wf/https://ieeexplore.ieee.org/document/9332553
    class ReluctanceZhangModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    // Based on A Novel Approach for 3D Air Gap Reluctance Calculations by Jonas Mühlethaler
    // https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf
    class ReluctanceMuehlethalerModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);

            double get_basic_reluctance(double l, double w, double h) {
                auto constants = Constants();
                return 1 / constants.vacuum_permeability / (w / 2 / l + 2 / std::numbers::pi * (1 + log(std::numbers::pi * h / 4 / l)));
            }

            double get_reluctance_type_1(double l, double w, double h) {
                double basic_reluctance = get_basic_reluctance(l, w, h);
                return 1 / (1 / (basic_reluctance + basic_reluctance) + 1 / (basic_reluctance + basic_reluctance));
            }
    };

    class ReluctanceEffectiveAreaModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    class ReluctanceEffectiveLengthModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    class ReluctancePartridgeModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    class ReluctanceClassicModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    // Based on Air-gap reluctance and inductance calculations for magnetic circuits using a Schwarz-Christoffel transformation by A. Balakrishnan
    // https://sci-hub.wf/https://ieeexplore.ieee.org/document/602560
    class ReluctanceBalakrishnanModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

    // Based on The Reluctance of Large Air Gaps in Ferrite Cores by Erika Stenglein
    // https://sci-hub.wf/10.1109/EPE.2016.7695271
    class ReluctanceStengleinModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);

            double u(double rx, double l1) {
                return 42.7 * rx / l1 - 50.2;
            }

            double v(double rx, double l1) {
                return -55.4 * rx / l1 + 71.6;
            }

            double w(double rx, double l1) {
                return 0.88 * rx / l1 - 0.80;
            }

            double alpha(double rx, double l1, double lg) {
                return u(rx, l1) * pow(lg / l1, 2) + v(rx, l1) * lg / l1 + w(rx, l1);
            }

    };

}