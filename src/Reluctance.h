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
#include <magic_enum.hpp>

#include <Core.h>

namespace OpenMagnetics {

    enum class ReluctanceModels : int { ZHANG, MCLYMAN, EFFECTIVE_AREA, EFFECTIVE_LENGTH, MUEHLETHALER };

    class ReluctanceModel {
        private:
        protected:
        public:
            virtual std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo) = 0;
            double get_core_reluctance(Core core) {
                auto constants = Constants();
                double absolute_permeability = constants.absolute_permeability; // HARDCODED TODO: replace when materials are implemented
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
                information["McLyman"] = "Based on the method described in page 8-11 from \"Transformer and Inductor Design Handbook Fourth Edition\" by Colonel Wm. T. McLyman.";
                information["Effective Area"] = "Based on the method described in page 60 from \"High-Frequency Magnetic Components, Second Edition\" by Marian Kazimierczuk.";
                information["Effective Length"] = "Based on the method described in page 60 from \"High-Frequency Magnetic Components, Second Edition\" by Marian Kazimierczuk.";
                return information;
            }
            static std::map<std::string, double> get_models_errors() {
                // This are taken manually from running the tests. Yes, a pain in the ass
                // TODO: Automate it somehow
                std::map<std::string, double> errors;
                errors["Zhang"] = 0.115811;
                errors["Muehlethaler"] = 0.110996;
                errors["McLyman"] = 0.124488;
                errors["Effective Area"] = 0.175055;
                errors["Effective Length"] = 0.175055;
                return errors;
            }
            static std::map<std::string, std::string> get_models_internal_links() {
                std::map<std::string, std::string> external_links;
                external_links["Zhang"] = "https://ieeexplore.ieee.org/document/9332553";
                external_links["Muehlethaler"] = "https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf";
                external_links["McLyman"] = "https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook";
                external_links["Effective Area"] = "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
                external_links["Effective Length"] = "https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33";
                return external_links;
            }
            static std::map<std::string, std::string> get_models_external_links() {
                std::map<std::string, std::string> internal_links;
                internal_links["Zhang"] = "";
                internal_links["Muehlethaler"] = "";
                internal_links["McLyman"] = "";
                internal_links["Effective Area"] = "";
                internal_links["Effective Length"] = "";
                return internal_links;
            }
            // virtual double get_total_reluctance() = 0;

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

    class ReluctanceMcLymanModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

}