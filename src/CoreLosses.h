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

#include <InputsWrapper.h>
#include <CoreWrapper.h>

namespace OpenMagnetics {

    enum class CoreLossesModels : int { STEINMETZ, IGSE, BARG, ROSHEN, ALBACH, NSE, MSE };

    class CoreLossesModel {
        private:
        protected:
        public:
            std::vector<double> _hysteresisMajorLoopTop;
            std::vector<double> _hysteresisMajorLoopBottom;
            std::vector<double> _hysteresisMajorH;
            std::vector<double> _hysteresisMinorLoopTop;
            std::vector<double> _hysteresisMinorLoopBottom;
            std::vector<double> _hysteresisMinorH;
            SteinmetzCoreLossesMethodRangeDatum _steinmetzDatum;
            bool _steinmetzDatumSet = false;
            CoreLossesModel() = default;
            virtual ~CoreLossesModel() = default;
            virtual std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature) = 0;
            static std::shared_ptr<CoreLossesModel> factory(CoreLossesModels modelName);
            static OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum get_steinmetz_coefficients(CoreMaterialDataOrNameUnion material, double frequency);
            bool is_steinmetz_datum_loaded() { return _steinmetzDatumSet; }
            SteinmetzCoreLossesMethodRangeDatum get_steinmetz_datum() { return _steinmetzDatum; }
            void set_steinmetz_datum(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) { _steinmetzDatumSet = true; _steinmetzDatum = steinmetzDatum; }

            static double apply_temperature_coefficients(double volumetricLosses, SteinmetzCoreLossesMethodRangeDatum steinmetzDatum, double temperature) {
                double volumetricLossesWithTemperature = volumetricLosses;
                if (steinmetzDatum.get_ct0() && steinmetzDatum.get_ct1() && steinmetzDatum.get_ct2()) {
                    double ct0 = steinmetzDatum.get_ct0().value();
                    double ct1 = steinmetzDatum.get_ct1().value();
                    double ct2 = steinmetzDatum.get_ct2().value();
                    volumetricLossesWithTemperature *= (ct2 * pow(temperature, 2) - ct1 * temperature + ct0 );
                }
                return volumetricLossesWithTemperature;
            }

            static std::vector<std::string> get_methods(CoreMaterialDataOrNameUnion material){
                OpenMagnetics::CoreMaterial materialData;
                // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid long loading operations)
                if (std::holds_alternative<std::string>(material) && std::get<std::string>(material) != "dummy") {
                    materialData = OpenMagnetics::find_data_by_name<OpenMagnetics::CoreMaterial>(std::get<std::string>(material));
                }
                else {
                    materialData = std::get<OpenMagnetics::CoreMaterial>(material);
                }

                std::vector<std::string> methods;
                auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();
                for(auto& volumetricLossesMethodVariant: volumetricLossesMethodsVariants) {
                    auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
                    for(auto& volumetricLossesMethod: volumetricLossesMethods) {
                        if (std::holds_alternative<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod)) {
                            auto methodData = std::get<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod);
                            methods.push_back(methodData.get_method());
                        }
                    }
                }

                std::vector<std::string> models;
                if (std::count(methods.begin(), methods.end(), "steinmetz")) {
                    models.push_back("Steinmetz");
                    models.push_back("iGSE");
                    models.push_back("Barg");
                    models.push_back("Albach");
                    models.push_back("MSE");
                }
                if (std::count(methods.begin(), methods.end(), "roshen")) {
                    models.push_back("Roshen");
                }
                return models;
            }

            static std::map<std::string, std::string> get_models_information() {
                std::map<std::string, std::string> information;
                information["Steinmetz"] = "Based on \"On the law of hysteresis\" by Charles Proteus Steinmetz";
                information["iGSE"] = "Based on \"Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters\" by Charles R. Sullivan";
                information["Barg"] = "Based on \"Core Loss Calculation of Symmetric Trapezoidal Magnetic Flux Density Waveform\" by Sobhi Barg";
                information["Roshen"] = "Based on \"Ferrite Core Loss for Power Magnetic Components Design\" and \"A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms\" by Waseem Roshen";
                information["Albach"] = "Based on \"Calculating Core Losses in Transformers for Arbitrary Magnetizing Currents A Comparison of Different Approaches\" by Manfred Albach";
                information["NSE"] = "Based on \"Measurement and Loss Model of Ferrites with Non-sinusoidal Waveforms\" by Alex Van den Bossche";
                information["MSE"] = "Based on \"Calculation of Losses in Ferro- and Ferrimagnetic Materials Based on the Modified Steinmetz Equation\" by Jürgen Reinert";
                return information;
            }
            static std::map<std::string, double> get_models_errors() {
                // This are taken manually from running the tests. Yes, a pain in the ass
                // TODO: Automate it somehow
                std::map<std::string, double> errors;
                errors["Steinmetz"] = 0.39353;
                errors["iGSE"] = 0.358237;
                errors["Barg"] = 0.374326;
                errors["Roshen"] = 0.487881;
                errors["Albach"] = 0.357267;
                errors["NSE"] = 0.358237;
                errors["MSE"] = 0.357267;
                return errors;
            }
            static std::map<std::string, std::string> get_models_external_links() {
                std::map<std::string, std::string> external_links;
                external_links["Steinmetz"] = "https://ieeexplore.ieee.org/document/1457110";
                external_links["iGSE"] = "http://inductor.thayerschool.org/papers/IGSE.pdf";
                external_links["Barg"] = "https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf";
                external_links["Roshen"] = "https://ieeexplore.ieee.org/document/4052433";
                external_links["Albach"] = "https://ieeexplore.ieee.org/iel3/3925/11364/00548774.pdf";
                external_links["NSE"] = "http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf";
                external_links["MSE"] = "https://www.mikrocontroller.net/attachment/129490/Modified_Steinmetz.pdf";
                return external_links;
            }
            static std::map<std::string, std::string> get_models_internal_links() {
                std::map<std::string, std::string> internal_links;
                internal_links["Steinmetz"] = "";
                internal_links["iGSE"] = "";
                internal_links["Barg"] = "";
                internal_links["Roshen"] = "/musings/4_roshen_method_core_losses";
                internal_links["Albach"] = "";
                internal_links["NSE"] = "";
                internal_links["MSE"] = "";
                return internal_links;
            }
    };

    // Based on On the law of hysteresis by Charles Proteus Steinmetz
    // https://sci-hub.wf/10.1109/proc.1984.12842
    class CoreLossesSteinmetzModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
    };


    // Based on Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters by Charles R. Sullivan
    // http://inductor.thayerschool.org/papers/IGSE.pdf
    class CoreLossesIGSEModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
            double get_ki(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum);
    };


    // Based on Core Loss Calculation of Symmetric Trapezoidal Magnetic Flux Density Waveform by Sobhi Barg
    // https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf
    class CoreLossesBargModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
    };


    // Based on Ferrite Core Loss for Power Magnetic Components Design and 
    // A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms by Waseem Roshen
    // https://sci-hub.wf/10.1109/20.278656
    // https://sci-hub.wf/10.1109/TPEL.2006.886608
    class CoreLossesRoshenModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
            double get_hysteresis_losses_density(std::map<std::string, double> parameters, OperationPointExcitation excitation);
            double get_eddy_current_losses_density(CoreWrapper core, OperationPointExcitation excitation, double resistivity);
            double get_excess_eddy_current_losses_density(OperationPointExcitation excitation, double resistivity, double alphaTimesN0);
            std::map<std::string, double> get_roshen_parameters(CoreMaterialDataOrNameUnion material, OperationPointExcitation excitation, double temperature);
    };


    // Based on Calculating Core Losses in Transformers for Arbitrary Magnetizing Currents A Comparison of Different Approaches by Manfred Albach
    // https://sci-hub.wf/10.1109/PESC.1996.548774
    class CoreLossesAlbachModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
    };


    // Based on Measurement and Loss Model of Ferrites with Non-sinusoidal Waveforms by Alex Van den Bossche
    // http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf
    class CoreLossesNSEModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
            double get_kn(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum);
    };


    // Based on Calculation of Losses in Ferro- and Ferrimagnetic Materials Based on the Modified Steinmetz Equation by Jürgen Reinert
    // https://sci-hub.wf/10.1109/28.936396
    class CoreLossesMSEModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
    };


    // Based on Improved Calculation of Core Loss With Nonsinusoidal Waveforms by Charles R. Sullivan
    // http://inductor.thayerschool.org/papers/gse.pdf
    class CoreLossesGSEModel: public CoreLossesModel {
        public:
            std::map<std::string, double> get_core_losses(CoreWrapper core, OperationPointExcitation excitation, double temperature);
    };
}