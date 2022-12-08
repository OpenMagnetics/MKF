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

    enum class ReluctanceModels : int { ZHANG, MCLYMAN, EFFECTIVE_AREA, MUEHLETHALER };

    class ReluctanceModel {
        private:
        protected:
        public:
            virtual std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo) = 0;
            double get_core_reluctance(Core core) {
                auto constants = Constants();
                double absolute_permeability = constants.vacuum_permeability * 2000; // HARDCODED TODO: replace when materials are implemented
                double effective_area = core.get_processed_description()->get_effective_parameters().get_effective_area();
                double effective_length = core.get_processed_description()->get_effective_parameters().get_effective_length();

                double reluctance_core = effective_length / (absolute_permeability * effective_area);
                return reluctance_core;
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

    class ReluctanceMcLymanModel: public ReluctanceModel {
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo);
    };

}