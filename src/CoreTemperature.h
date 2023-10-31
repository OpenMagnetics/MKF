#pragma once
#include "Constants.h"

#include <CoreWrapper.h>
#include <MAS.hpp>
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

enum class CoreTemperatureModels : int {
    KAZIMIERCZUK,
    MANIKTALA,
    TDK,
    DIXON,
    AMIDON
};

class CoreTemperatureModel {
  private:
  protected:
  public:
    virtual ~CoreTemperatureModel() = default;
    
    static std::shared_ptr<CoreTemperatureModel> factory(CoreTemperatureModels modelName);
    virtual TemperatureOutput get_core_temperature(CoreWrapper core,
                                                   double coreLosses,
                                                   double ambientTemperature) = 0;

    static std::map<std::string, std::string> get_models_information() {
        std::map<std::string, std::string> information;
        information["Maniktala"] = R"(Based on "witching Power Supplies A - Z, 2nd edition" by Sanjaya Maniktala)";
        information["Kazimierczuk"] =
            R"(Based on "High-Frequency Magnetic Components 2nd Edition" by Marian Kazimierczuk)";
        information["TDK"] = R"(Based on "Ferrites and Accessories" by TDK)";
        information["Dixon"] = R"(Based on "Design of Flyback Transformers and Filter Inductors" by Lloyd H. Dixon)";
        information["Amidon"] = R"(Based on "Iron Powder Core Loss Characteristics" by Amidon)";
        return information;
    }
    static std::map<std::string, double> get_models_errors() {
        // This are taken manually from running the tests. Yes, a pain in the ass
        // TODO: Automate it somehow
        std::map<std::string, double> errors;
        errors["Kazimierczuk"] = 0.2577;
        errors["Maniktala"] = 0.2476;
        errors["TDK"] = 0.5148;
        errors["Dixon"] = 0.2464;
        errors["Amidon"] = 0.2537;
        return errors;
    }
    static std::map<std::string, std::string> get_models_external_links() {
        std::map<std::string, std::string> external_links;
        external_links["Kazimierczuk"] =
            "https://www.goodreads.com/book/show/18861402-high-frequency-magnetic-components";
        external_links["Maniktala"] = "https://www.goodreads.com/book/show/12042906-switching-power-supplies-a-z";
        external_links["TDK"] = "https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/"
                                "pdf-applicationnotes.pdf";
        external_links["Dixon"] = "https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086";
        external_links["Amidon"] = "https://www.amidoncorp.com/product_images/specifications/1-38.pdf";
        return external_links;
    }
    static std::map<std::string, std::string> get_models_internal_links() {
        std::map<std::string, std::string> internal_links;
        internal_links["Kazimierczuk"] = "";
        internal_links["Maniktala"] = "";
        internal_links["TDK"] = "";
        internal_links["Dixon"] = "";
        internal_links["Amidon"] = "";
        return internal_links;
    }
};

// Based on Switching Power Supplies A - Z, 2nd edition by Sanjaya Maniktala, page 154
class CoreTemperatureManiktalaModel : public CoreTemperatureModel {
  public:
    TemperatureOutput get_core_temperature(CoreWrapper core, double coreLosses, double ambientTemperature);
};

// Based on High-Frequency Magnetic Components 2nd Edition by Marian Kazimierczuk, page 151
class CoreTemperatureKazimierczukModel : public CoreTemperatureModel {
  public:
    TemperatureOutput get_core_temperature(CoreWrapper core, double coreLosses, double ambientTemperature);
};

// Based on Ferrites and Accessories by TDK, page 23
// https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/pdf-applicationnotes.pdf
class CoreTemperatureTdkModel : public CoreTemperatureModel {
  public:
    TemperatureOutput get_core_temperature(CoreWrapper core, double coreLosses, double ambientTemperature);
};

// Based on Design of Flyback Transformers and Filter Inductors by Lloyd H. Dixon, page 2-4
// https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086
class CoreTemperatureDixonModel : public CoreTemperatureModel {
  public:
    TemperatureOutput get_core_temperature(CoreWrapper core, double coreLosses, double ambientTemperature);
};

// Based on Iron Powder Core Loss Characteristics by Amidon
// https://www.amidoncorp.com/product_images/specifications/1-38.pdf
class CoreTemperatureAmidonModel : public CoreTemperatureModel {
  public:
    TemperatureOutput get_core_temperature(CoreWrapper core, double coreLosses, double ambientTemperature);
};

} // namespace OpenMagnetics