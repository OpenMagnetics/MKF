#pragma once

#include "constructive_models/Mas.h"
#include "constructive_models/Magnetic.h"
#include "processors/Inputs.h"
#include <MAS.hpp>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace OpenMagnetics {

/**
 * @brief Verdict of one design-for-manufacturing rule.
 *
 * NOT_EVALUATED and NOT_APPLICABLE are deliberately distinct and neither is ever
 * silently turned into a PASS: a rule that could not read what it needs says so and
 * carries the reason (house rule: no fallbacks, no silent defaults).
 */
enum class ManufacturabilityStatus : int {
    PASS,            ///< The design complies with the rule.
    WARNING,         ///< Manufacturable, but against Wuerth's stated practice.
    FAIL,            ///< The design breaks the rule.
    INFORMATIONAL,   ///< Report-only rule: numbers, no verdict (R17).
    NOT_APPLICABLE,  ///< The rule does not apply to this design (e.g. an SMT rule on a THT part).
    NOT_EVALUATED    ///< The rule could not be evaluated; `reason` says why.
};

/**
 * @brief One finding of the manufacturability report: exactly one per rule.
 */
class ManufacturabilityFinding {
    private:
        std::string ruleId;
        std::string title;
        ManufacturabilityStatus status = ManufacturabilityStatus::NOT_EVALUATED;
        std::string message;
        std::string source;
        std::optional<double> measuredValue;
        std::optional<double> thresholdValue;
        std::optional<std::string> unit;
        std::optional<std::string> scope;
        std::optional<std::string> reason;

    public:
        ManufacturabilityFinding() = default;
        ManufacturabilityFinding(const std::string& ruleId, const std::string& title, const std::string& source) :
            ruleId(ruleId), title(title), source(source) {}

        const std::string& get_rule_id() const { return ruleId; }
        void set_rule_id(const std::string& value) { ruleId = value; }

        const std::string& get_title() const { return title; }
        void set_title(const std::string& value) { title = value; }

        ManufacturabilityStatus get_status() const { return status; }
        void set_status(ManufacturabilityStatus value) { status = value; }

        const std::string& get_message() const { return message; }
        void set_message(const std::string& value) { message = value; }

        const std::string& get_source() const { return source; }
        void set_source(const std::string& value) { source = value; }

        std::optional<double> get_measured_value() const { return measuredValue; }
        void set_measured_value(std::optional<double> value) { measuredValue = value; }

        std::optional<double> get_threshold_value() const { return thresholdValue; }
        void set_threshold_value(std::optional<double> value) { thresholdValue = value; }

        std::optional<std::string> get_unit() const { return unit; }
        void set_unit(std::optional<std::string> value) { unit = value; }

        std::optional<std::string> get_scope() const { return scope; }
        void set_scope(std::optional<std::string> value) { scope = value; }

        std::optional<std::string> get_reason() const { return reason; }
        void set_reason(std::optional<std::string> value) { reason = value; }
};

void to_json(json& j, const ManufacturabilityFinding& x);

/**
 * @brief The design-for-manufacturing report: one finding per rule, R1 to R17.
 *
 * Rules owned by other work packages are present with status NOT_EVALUATED and the
 * owning work package as the reason, so the report is complete and honest.
 */
class ManufacturabilityReport {
    private:
        std::vector<ManufacturabilityFinding> findings;

    public:
        const std::vector<ManufacturabilityFinding>& get_findings() const { return findings; }
        std::vector<ManufacturabilityFinding>& get_mutable_findings() { return findings; }
        void add_finding(const ManufacturabilityFinding& finding) { findings.push_back(finding); }

        /**
         * @brief The finding for one rule id (e.g. "R1").
         * @throws std::runtime_error when the report carries no such rule.
         */
        const ManufacturabilityFinding& get_finding(const std::string& ruleId) const;

        bool has_finding(const std::string& ruleId) const;

        size_t count_by_status(ManufacturabilityStatus status) const;
};

void to_json(json& j, const ManufacturabilityReport& x);

/**
 * @class Manufacturability
 * @brief Design-for-manufacturing rule pack (ABT #1177, WP8).
 *
 * Encodes as checks the public Wuerth Elektronik design-for-manufacturing guidance
 * catalogued in docs/2026-09-12_manufacturing_features_proposal.md. Every number the
 * rules compare against lives in src/data/dfm_rules.json together with the source that
 * states it; nothing is hardcoded here.
 *
 * No rule changes the design. The report is an MKF struct, never MAS output.
 */
class Manufacturability {
    public:
        Manufacturability();

        /// The rule data as loaded from src/data/dfm_rules.json.
        const json& get_rules() const { return _rules; }

        /// The whole report, one finding per rule R1..R17.
        ManufacturabilityReport calculate_report(Mas& mas);
        ManufacturabilityReport calculate_report(Magnetic& magnetic, Inputs& inputs);

        // One method per standalone rule, so each can be tested (and reverted) on its own.
        ManufacturabilityFinding evaluate_r1_layer_parity(Magnetic& magnetic);
        ManufacturabilityFinding evaluate_r5_smt_heavy_wire(Magnetic& magnetic, Inputs& inputs);
        ManufacturabilityFinding evaluate_r6_thermoset_bobbin(Magnetic& magnetic, Inputs& inputs);
        ManufacturabilityFinding evaluate_r7_size_heuristic(Magnetic& magnetic, Inputs& inputs);
        ManufacturabilityFinding evaluate_r8_insulated_wire(Magnetic& magnetic);
        ManufacturabilityFinding evaluate_r9_margin_tape(Magnetic& magnetic);
        ManufacturabilityFinding evaluate_r12_flying_leads(Magnetic& magnetic, Inputs& inputs);
        ManufacturabilityFinding evaluate_r13_manual_termination(Magnetic& magnetic);
        ManufacturabilityFinding evaluate_r15_lateral_gaps(Magnetic& magnetic, Inputs& inputs);
        ManufacturabilityFinding evaluate_r17_production_tests(Magnetic& magnetic, Inputs& inputs);

        /// The stub finding for a rule owned by another work package.
        ManufacturabilityFinding not_evaluated_rule(const std::string& ruleId);

        /**
         * @brief Classification of a bobbin material name.
         *
         * UNKNOWN is returned for a name that is in neither table; the caller must report
         * it as unknown, never assume one of the two.
         */
        enum class BobbinMaterialClass : int { THERMOPLASTIC, THERMOSET, UNKNOWN };
        BobbinMaterialClass classify_bobbin_material(const std::string& materialName) const;

        /// True when the wire is an insulated (extruded) wire: TIW/TEX-E and its 1/2-layer relatives.
        static bool is_insulated_wire(Wire& wire);
        /// True when the wire is a triple insulated wire (three extruded layers).
        static bool is_triple_insulated_wire(Wire& wire);
        /// True when the wire is a fully insulated wire (enamelled, build grade outside 1..3).
        static bool is_fully_insulated_wire(Wire& wire);
        /// True when the wire is plain enamelled magnet wire of a standard build grade.
        static bool is_enamelled_magnet_wire(Wire& wire);

    private:
        json _rules;

        /// The rule block for a rule id, or a throw if the data file has no such rule.
        const json& rule(const std::string& ruleId) const;
        /// A number from a rule block, or a throw naming the missing key.
        double rule_number(const std::string& ruleId, const std::string& key) const;
        std::string rule_string(const std::string& ruleId, const std::string& key) const;
        ManufacturabilityFinding make_finding(const std::string& ruleId) const;
};

} // namespace OpenMagnetics
