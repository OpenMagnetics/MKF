// SPDX-License-Identifier: MIT
// Magnetic Blade Runner: the class table. Every number below is a PUBLISHED datasheet value and
// carries its source (document, page/table). The documents were collected on 2026-09-28; the
// full list with URLs and document dates is kept with the validator's research notes
// (datasheets/ANCHORS.md). Losses are sinusoidal, peak flux density, W/m^3.
//
// Kinds of value: TAB = read from a table; FIG = read off a curve (uncertainty stated);
// FIT = computed from the vendor's own curve-fit equation in the same document (checked against
// a table value of that document); MAX = a guaranteed maximum.
//
// A class is chosen ONLY from the record's MAS `material` / `materialComposition` /
// `application` fields. A record those fields do not place is "unclassified" and every
// class-dependent check reports itself skipped for it.
#include "support/MaterialValidator.h"

#include <algorithm>

namespace OpenMagnetics {

namespace {

MaterialLossAnchor anchor(std::string manufacturer, std::string grade, double kWm3, std::string source) {
    return MaterialLossAnchor{std::move(manufacturer), std::move(grade), kWm3 * 1e3, std::move(source)};
}

std::vector<MaterialClassBounds> build_class_table() {
    std::vector<MaterialClassBounds> table;

    // ---------------------------------------------------------------- MnZn ferrite
    const std::string mnznPermeabilitySource =
        "mu_i 750: Magnetics L (Magnetics 2021 Ferrite Catalog p6), Ferroxcube 3F46 (MDS 2016-03-03 p2), ACME P63 "
        "(ACME catalogue Oct 2024 pdf p38); mu_i 15000: Magnetics M (2021 Ferrite Catalog p6), ACME A151 (ACME "
        "catalogue p10). TAB";
    const std::string mnznSaturationSource =
        "Bs 600 mT at 25 C: ACME P491 (ACME catalogue Oct 2024 p10, 10 kHz, 1200 A/m). TAB";
    const std::string mnznCurieSource =
        "Tc 100 C: ACME N10 (ACME catalogue p10); Tc 300 C: ACME P491 (>=300, ACME catalogue p10). TAB";
    {
        MaterialClassBounds c;
        c.name = "MnZn power ferrite";
        c.anchorQuality = "firm: table values from TDK, TDK-EPCOS, Ferroxcube, ACME, DMEGC, Magnetics at the stated points";
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.200, 100, true,
            anchor("TDK", "PC47", 250,
                   "TDK 'Mn-Zn Ferrite Material characteristics' (April 2026) p3 list and p6 table, TAB; also "
                   "Ferroxcube 3C98 (MDS 2013-06-03 p2) and ACME P48 (ACME catalogue Oct 2024 p12)"),
            anchor("TDK-EPCOS", "N41", 1400,
                   "TDK-EPCOS SIFERRIT N41 datasheet (June 2025) p2, TAB; ACME P491 1390 (ACME catalogue p13)")});
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.100, 100, false,
            anchor("Ferroxcube", "3C96", 40, "Ferroxcube 3C96 MDS (2008-09-01) p2, TAB (typical)"),
            anchor("Magnetics", "F", 110, "Magnetics 2021 Ferrite Catalog pdf p6 (printed p3), TAB")});
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            500e3, 0.050, 100, false,
            anchor("ACME", "P63", 19,
                   "ACME measurement workbook '40 P63.xlsx' (ACME/USIG e-mail 2025-06-02), sheet 5data, row "
                   "'500kHz;50mT', 100 C column, TAB; TDG TP5R's manufacturer points in MAS give 19.4. TDK-EPCOS "
                   "PC200 (datasheet June 2025 p2), Ferroxcube 3F37 and ACME P53 publish about 60"),
            anchor("DMEGC", "DMR55", 200, "DMEGC DMR55 material datasheet p1, TAB (worst MHz-marketed grade)")});
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            1e6, 0.050, 100, false,
            anchor("ACME", "P63", 80, "ACME catalogue Oct 2024 p13 table 4 and P63 page (pdf p38), TAB"),
            anchor("ACME", "P52", 1000, "ACME catalogue Oct 2024 p13 table 4, TAB")});
        c.minimumInitialPermeability = 750;
        c.maximumInitialPermeability = 15000;
        c.initialPermeabilitySource = mnznPermeabilitySource;
        c.maximumSaturationFluxDensity = 0.600;
        c.saturationSource = mnznSaturationSource;
        c.minimumCurieTemperature = 100;
        c.maximumCurieTemperature = 300;
        c.curieSource = mnznCurieSource;
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "MnZn ferrite (signal, EMI or unspecified application)";
        c.lossAnchorGap = "no published volumetric-loss anchor for non-power MnZn grades";
        c.anchorQuality = "static bounds firm (TAB); no loss anchor";
        c.minimumInitialPermeability = 750;
        c.maximumInitialPermeability = 15000;
        c.initialPermeabilitySource = mnznPermeabilitySource;
        c.maximumSaturationFluxDensity = 0.600;
        c.saturationSource = mnznSaturationSource;
        c.minimumCurieTemperature = 100;
        c.maximumCurieTemperature = 300;
        c.curieSource = mnznCurieSource;
        table.push_back(c);
    }

    // ---------------------------------------------------------------- NiZn ferrite (WEAK)
    {
        MaterialClassBounds c;
        c.name = "NiZn ferrite";
        c.weak = true;
        c.anchorQuality = "WEAK: no point common to all makers; best is a curve reading (+-15%), worst one table value";
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            5e6, 0.010, 100, true,
            anchor("Fair-Rite", "67", 160, "Fair-Rite 67 Pv vs B at 100 C (fair-rite.com 67PLB100.jpg, 2020), FIG +-15%"),
            anchor("Proterial", "NL12S", 340, "Proterial ferrite catalogue HJ-B3 (2023) pdf p37, TAB")});
        c.minimumInitialPermeability = 5;
        c.maximumInitialPermeability = 2500;
        c.initialPermeabilitySource =
            "mu_i 5: Encore N4C (Encore Electronics Ni-Zn material table, encores.com.tw); mu_i 2500: ACME K25 "
            "(ACME catalogue p11). TAB";
        c.maximumSaturationFluxDensity = 0.500;
        c.saturationSource = "Bs 500 mT: Proterial NB25S at 8 kA/m (HJ-B3 catalogue). TAB";
        c.minimumCurieTemperature = 80;
        c.maximumCurieTemperature = 500;
        c.curieSource = "Tc 80 C: TAK L28A (>80 C, takferrite.com/product_detail?id=137); Tc 500 C: Fair-Rite 68 (>500, "
                        "fair-rite.com 68 page). TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "MgZn ferrite";
        c.lossAnchorGap = "no published anchor collected for MgZn ferrite";
        c.anchorQuality = "no anchors";
        table.push_back(c);
    }

    // ---------------------------------------------------------------- powders (WEAK)
    const std::string powderQuality =
        "WEAK: vendors do not state the measuring temperature (evaluated at 25 C); some values are vendor-fit "
        "results";
    {
        MaterialClassBounds c;
        c.name = "MPP powder (FeNiMo)";
        c.weak = true;
        c.anchorQuality = powderQuality;
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.100, 25, true,
            anchor("Magnetics", "MPP 60u", 450, "Magnetics Powder Core Catalog (2025) pdf p12 comparison table, TAB"),
            anchor("Magnetics", "MPP 550u", 1150, "Magnetics Powder Core Catalog pdf p111 loss fit, FIT")});
        c.minimumInitialPermeability = 14;
        c.maximumInitialPermeability = 550;
        c.initialPermeabilitySource = "Magnetics Powder Core Catalog, MPP 14-550 mu. TAB";
        c.maximumSaturationFluxDensity = 0.89;
        c.saturationSource = "Micrometals MPP datasheet (2024) Bsat 8.9 kG; Magnetics MPP 0.8 T. TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "FeNi powder (High Flux, Edge)";
        c.weak = true;
        c.anchorQuality = powderQuality;
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.100, 25, true,
            anchor("Magnetics", "Edge 60u", 375, "Magnetics Powder Core Catalog pdf p12 comparison table, TAB"),
            anchor("Magnetics", "High Flux 14u", 1400, "Magnetics Powder Core Catalog pdf p111 loss fit, FIT")});
        c.minimumInitialPermeability = 14;
        c.maximumInitialPermeability = 160;
        c.initialPermeabilitySource = "Magnetics Powder Core Catalog: High Flux 14-160 mu, Edge 14-125 mu. TAB";
        c.maximumSaturationFluxDensity = 1.51;
        c.saturationSource = "Micrometals HiFlux datasheet Bsat 15.1 kG; Magnetics High Flux / Edge 1.5 T. TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "FeSiAl powder";
        c.weak = true;
        c.anchorQuality = powderQuality;
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.100, 25, true,
            anchor("Magnetics", "Kool Mu Ultra 60u", 160, "Magnetics Powder Core Catalog pdf p12 comparison table, TAB"),
            anchor("Micrometals", "Sendust 26u", 1218,
                   "Micrometals Sendust 26u datasheet (2020) loss formula at 100 kHz/100 mT, FIT")});
        c.minimumInitialPermeability = 14;
        c.maximumInitialPermeability = 125;
        c.initialPermeabilitySource = "Magnetics Powder Core Catalog, Kool Mu 14-125 mu. TAB";
        c.maximumSaturationFluxDensity = 1.2;
        c.saturationSource = "Bs 1.2 T: KDM KPH-HT/HP 12,000 G, listed under Sendust (KDM Alloy Powder Core brochure "
                             "2026-06, p5 table); Magnetics Kool Mu 1.0 T (Powder Core Catalog). TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "FeSi powder";
        c.weak = true;
        c.anchorQuality = powderQuality;
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.100, 25, true,
            anchor("Magnetics", "XFlux Ultra 60u", 1000, "Magnetics Powder Core Catalog pdf p111 loss fit, FIT (p12 table 1035)"),
            anchor("Micrometals", "FluxSan 60u", 1615,
                   "Micrometals FluxSan 60u datasheet (2020) loss formula at 100 kHz/100 mT, FIT")});
        c.maximumSaturationFluxDensity = 1.76;
        c.saturationSource = "Micrometals FluxSan datasheet Bsat 17.6 kG; Magnetics XFlux 1.6 T. TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "iron powder (iron, carbonyl iron)";
        c.weak = true;
        c.anchorQuality = powderQuality + "; only Micrometals publishes, at its own point 100 kHz / 14 mT";
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.014, 25, true,
            anchor("Micrometals", "-2 / -14", 18, "Micrometals mix -2 / -14 datasheets (2019), standard point, TAB (nominal)"),
            anchor("Micrometals", "-40", 127, "Micrometals mix -40 datasheet (2019), standard point, TAB (nominal)")});
        c.maximumSaturationFluxDensity = 1.89;
        c.saturationSource = "Micrometals mix -45 datasheet (2019) Bsat 18.9 kG. TAB";
        table.push_back(c);
    }

    // ---------------------------------------------------------------- nanocrystalline
    {
        MaterialClassBounds c;
        c.name = "nanocrystalline";
        c.anchorQuality = "firm: one maker's table (Proterial FINEMET) at a stated point and temperature";
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.200, 25, true,
            anchor("Proterial", "FT-3 14um", 230,
                   "Proterial FINEMET ribbon catalogue PR-EM07, 'Standard Magnetic Characteristics' table, TAB"),
            anchor("Proterial", "FT-8 18um", 470, "Proterial FINEMET ribbon catalogue PR-EM07, same table, TAB")});
        c.maximumSaturationFluxDensity = 1.30;
        c.saturationSource = "Proterial FINEMET PR-EM07: FT-8 Bs 1.30 T. TAB";
        c.minimumCurieTemperature = 570;
        c.curieSource = "Proterial FINEMET PR-EM07: Tc ~570 C. TAB";
        table.push_back(c);
    }

    // ---------------------------------------------------------------- amorphous (WEAK)
    const std::string amorphousSaturationSource =
        "Metglas technical bulletins: 2605HB1M Bs 1.63 T (the highest). TAB";
    {
        MaterialClassBounds c;
        c.name = "Fe-based amorphous (FeSi)";
        c.weak = true;
        c.anchorQuality = "WEAK: best is a curve reading (+-35%, 85-dpi plot), worst an unnamed grade in a comparison table";
        c.lossReferencePoints.push_back(MaterialLossReferencePoint{
            100e3, 0.200, 25, true,
            anchor("Metglas", "2605S3A", 770, "Metglas 2605S3A technical bulletin p2 loss curve, FIG +-35% (temperature unstated)"),
            anchor("Proterial", "Fe-based amorphous 25um", 2200,
                   "Proterial FINEMET ribbon catalogue PR-EM07 comparison table, TAB")});
        c.maximumSaturationFluxDensity = 1.63;
        c.saturationSource = amorphousSaturationSource;
        c.minimumCurieTemperature = 358;
        c.maximumCurieTemperature = 395;
        c.curieSource = "Metglas technical bulletins: 2605S3A 358 C .. 2605SA1 395 C (Fe-based alloys). TAB";
        table.push_back(c);
    }
    {
        MaterialClassBounds c;
        c.name = "amorphous (composition not stated)";
        c.weak = true;
        c.lossAnchorGap = "MAS gives no composition, so Fe-based and Co-based amorphous (whose losses differ ~3x) cannot be told apart";
        c.anchorQuality = "static bounds only";
        c.maximumSaturationFluxDensity = 1.63;
        c.saturationSource = amorphousSaturationSource;
        c.minimumCurieTemperature = 225;
        c.maximumCurieTemperature = 395;
        c.curieSource = "Metglas technical bulletins: 2714A 225 C .. 2605SA1 395 C. TAB";
        table.push_back(c);
    }

    // ---------------------------------------------------------------- electrical steel / tape
    {
        MaterialClassBounds c;
        c.name = "electrical steel / Ni-Fe tape";
        c.lossAnchorGap = "no class envelope: the MAS type spans SiFe steel to 80% Ni tape (Proterial PR-EM07: 1000..8400 kW/m3 at 100 kHz/0.2 T)";
        c.anchorQuality = "no class envelope";
        table.push_back(c);
    }
    return table;
}

bool has_application(const json& record, const std::string& application) {
    if (!record.contains("application") || !record["application"].is_array()) {
        return false;
    }
    for (const auto& a : record["application"]) {
        if (a.is_string() && a.get<std::string>() == application) {
            return true;
        }
    }
    return false;
}

std::string string_field(const json& record, const char* key) {
    if (record.contains(key) && record[key].is_string()) {
        return record[key].get<std::string>();
    }
    return "";
}

const MaterialClassBounds* find_class(const std::string& name) {
    for (const auto& c : MaterialValidator::class_table()) {
        if (c.name == name) {
            return &c;
        }
    }
    return nullptr;
}

} // namespace

const std::vector<MaterialClassBounds>& MaterialValidator::class_table() {
    static const std::vector<MaterialClassBounds> table = build_class_table();
    return table;
}

const MaterialClassBounds* classify_material(const json& record) {
    const std::string type = string_field(record, "material");
    const std::string composition = string_field(record, "materialComposition");

    if (type == "ferrite") {
        if (composition == "MnZn") {
            return has_application(record, "power") ? find_class("MnZn power ferrite")
                                                    : find_class("MnZn ferrite (signal, EMI or unspecified application)");
        }
        if (composition == "NiZn") return find_class("NiZn ferrite");
        if (composition == "MgZn") return find_class("MgZn ferrite");
        return nullptr;
    }
    if (type == "powder") {
        if (composition == "FeNiMo") return find_class("MPP powder (FeNiMo)");
        if (composition == "FeNi") return find_class("FeNi powder (High Flux, Edge)");
        if (composition == "FeSiAl") return find_class("FeSiAl powder");
        if (composition == "FeSi") return find_class("FeSi powder");
        if (composition == "iron" || composition == "carbonylIron") return find_class("iron powder (iron, carbonyl iron)");
        return nullptr;  // proprietary, FeSiCr, FeMo or not stated
    }
    if (type == "nanocrystalline") return find_class("nanocrystalline");
    if (type == "amorphous") {
        if (composition == "FeSi") return find_class("Fe-based amorphous (FeSi)");
        return find_class("amorphous (composition not stated)");
    }
    if (type == "electricalSteel") return find_class("electrical steel / Ni-Fe tape");
    return nullptr;
}

json MaterialValidator::class_table_json() {
    json out = json::array();
    auto anchorJson = [](const MaterialLossAnchor& a) {
        return json{{"manufacturer", a.manufacturer}, {"grade", a.grade}, {"volumetricLosses", a.volumetricLosses}, {"source", a.source}};
    };
    for (const auto& c : class_table()) {
        json j;
        j["name"] = c.name;
        j["weak"] = c.weak;
        j["anchorQuality"] = c.anchorQuality;
        j["lossReferencePoints"] = json::array();
        for (const auto& p : c.lossReferencePoints) {
            j["lossReferencePoints"].push_back(json{{"frequency", p.frequency},
                                                    {"magneticFluxDensityPeak", p.magneticFluxDensityPeak},
                                                    {"temperature", p.temperature},
                                                    {"primary", p.primary},
                                                    {"best", anchorJson(p.best)},
                                                    {"worst", anchorJson(p.worst)}});
        }
        if (!c.lossAnchorGap.empty()) j["lossAnchorGap"] = c.lossAnchorGap;
        if (c.minimumInitialPermeability) j["minimumInitialPermeability"] = *c.minimumInitialPermeability;
        if (c.maximumInitialPermeability) j["maximumInitialPermeability"] = *c.maximumInitialPermeability;
        if (!c.initialPermeabilitySource.empty()) j["initialPermeabilitySource"] = c.initialPermeabilitySource;
        if (c.maximumSaturationFluxDensity) j["maximumSaturationFluxDensity"] = *c.maximumSaturationFluxDensity;
        if (!c.saturationSource.empty()) j["saturationSource"] = c.saturationSource;
        if (c.minimumCurieTemperature) j["minimumCurieTemperature"] = *c.minimumCurieTemperature;
        if (c.maximumCurieTemperature) j["maximumCurieTemperature"] = *c.maximumCurieTemperature;
        if (!c.curieSource.empty()) j["curieSource"] = c.curieSource;
        out.push_back(j);
    }
    return out;
}

} // namespace OpenMagnetics
