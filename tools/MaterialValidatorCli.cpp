// SPDX-License-Identifier: MIT
// mkf_material_validator: the Magnetic Blade Runner over a MAS core-material catalogue.
//
//   mkf_material_validator [--json] [--all] <core_materials.ndjson> [advanced_core_materials.ndjson]
//   mkf_material_validator --embedded      (the catalogue MKF was built with)
//   mkf_material_validator --class-table   (the envelopes and bounds, with their sources)
//
// Prints per-rule counts, the provenance summary and every IMPOSSIBLE and SUSPICIOUS finding
// (--all adds the skipped checks per record; --json prints the whole summary as JSON).
// Exit status: 0 no IMPOSSIBLE record, 1 at least one, 2 usage or I/O error.
#include "support/MaterialValidator.h"

#include <iostream>

using namespace OpenMagnetics;

int main(int argc, char** argv) {
    bool asJson = false, all = false, embedded = false;
    std::vector<std::string> paths;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--json") asJson = true;
        else if (arg == "--all") all = true;
        else if (arg == "--embedded") embedded = true;
        else if (arg == "--class-table") {
            std::cout << MaterialValidator::class_table_json().dump(2) << std::endl;
            return 0;
        }
        else if (arg == "-h" || arg == "--help") {
            std::cout << "usage: mkf_material_validator [--json] [--all] (--embedded | <core_materials.ndjson> [advanced_core_materials.ndjson])\n"
                         "       mkf_material_validator --class-table\n";
            return 0;
        }
        else paths.push_back(arg);
    }
    if (!embedded && (paths.empty() || paths.size() > 2)) {
        std::cerr << "usage: mkf_material_validator [--json] [--all] (--embedded | <core_materials.ndjson> [advanced_core_materials.ndjson])\n";
        return 2;
    }

    std::vector<MaterialVerdict> verdicts;
    try {
        MaterialValidator validator;
        if (embedded) verdicts = validator.validate_embedded_catalogue();
        else verdicts = validator.validate_catalogue(paths[0], paths.size() == 2 ? std::optional<std::string>(paths[1]) : std::nullopt);
    }
    catch (const std::exception& e) {
        std::cerr << "mkf_material_validator: " << e.what() << std::endl;
        return 2;
    }

    json summary = summarize_material_verdicts(verdicts);
    if (all) summary["verdicts"] = verdicts;
    if (asJson) {
        std::cout << summary.dump(2) << std::endl;
    }
    else {
        std::cout << "Magnetic Blade Runner: " << summary["records"] << " materials, " << summary["invalid"] << " IMPOSSIBLE, "
                  << summary["recordsWithSuspicious"] << " with SUSPICIOUS findings, " << summary["unclassified"] << " unclassified\n";
        if (summary.contains("provenance")) std::cout << summary["provenance"].get<std::string>() << "\n";
        std::cout << "\nfindings by rule:\n";
        for (const auto& [code, counts] : summary["findingsByCode"].items()) {
            std::cout << "  " << code;
            for (const auto& [severity, n] : counts.items()) std::cout << "  " << severity << "=" << n;
            std::cout << "\n";
        }
        std::cout << "\nskipped checks by rule:\n";
        for (const auto& [code, n] : summary["skippedByCode"].items()) std::cout << "  " << code << "  " << n << "\n";
        for (const char* severity : {"impossible", "suspicious"}) {
            std::cout << "\n" << severity << ":\n";
            for (const auto& f : summary[severity]) {
                std::cout << "  [" << f["code"].get<std::string>() << "] " << f["reference"].get<std::string>() << ": " << f["message"].get<std::string>() << "\n";
            }
        }
        if (all) {
            std::cout << "\nskipped per record:\n";
            for (const auto& v : verdicts) {
                for (const auto& s : v.skipped) std::cout << "  " << v.reference << ": " << s << "\n";
            }
        }
    }
    return summary["invalid"].get<size_t>() > 0 ? 1 : 0;
}
