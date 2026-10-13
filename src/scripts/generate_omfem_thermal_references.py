#!/usr/bin/env python3
# Usage: python3 src/scripts/generate_omfem_thermal_references.py tests/testData/omfem_thermal_2d_references.json
#
# Generates tests/testData/omfem_thermal_2d_references.json: for every MAS example, run
# OMFEM's 2D pipeline (mesh -> FEM losses -> radiating thermal FEA) and record the losses
# it used plus the temperatures it produced. The MKF test battery then drives MKF's thermal
# network with the SAME stored losses, so the comparison is loss-model-independent and the
# reference cannot silently go stale when loss models change (the lesson of the frozen
# Icepak bands, ABT #461).
import json, subprocess, sys, os, glob, re, time
_MKF_ROOT = __import__("pathlib").Path(__file__).resolve().parents[2]  # ABT #1596: this checkout, not a developer path

OMFEM = os.environ.get("OMFEM_MAS")
if not OMFEM:
    sys.exit("Set OMFEM_MAS to the omfem_mas binary (e.g. <OMFEM checkout>/build/omfem_mas)")
EXAMPLES = sorted(glob.glob(str(_MKF_ROOT / "MAS" / "examples/*.json")))
OUT = sys.argv[1]
TIMEOUT = 1200

# Examples that run OUTSIDE their core material's tabulated/fitted data, where MKF refuses by default.
# 01_1Hz asks N87's core loss at 1 Hz, below the span its Steinmetz coefficients are fitted over (25 kHz to
# 1 MHz). Alf's decision (2026-10-12, ABT #1796): generate it with MKF's explicit opt-in
# Settings::allowMaterialDataExtrapolation, which omfem_mas forwards when it is given
# --allow-material-data-extrapolation (MasSolveOptions::allow_material_data_extrapolation; a run-time
# option, off by default, no environment variable). Every other example runs with it off. The entry
# records that it was on and MKF's own warning per extrapolated quantity, as omfem_mas returns them.
EXTRAPOLATED_EXAMPLES = {"01_simple_inductor_etd34_n87_1Hz.json"}
LOG_DIR = os.environ.get("OMFEM_REF_LOG_DIR")  # optional: keep each example's omfem_mas stderr there

refs = {"_meta": {
    "generator": "gen_omfem_refs.py (ABT #461 / #454)",
    "source": "omfem_mas: 2D mesh -> FEM losses -> solve_thermal_fem (h from MKF natural-convection law per boundary, out-of-plane face sink, open toroid bore with view factor, k_core from MKF resolve_core_material_heat_conductivity, eps=0.9, radiating; ABT #1751)",
    "losses": "OMFEM's own FEM losses, stored here and fed IDENTICALLY to both models by the test",
    "date": time.strftime("%Y-%m-%d"),
}, "examples": {}}

# Resume support: keep already-successful entries, retry failures.
if os.path.exists(OUT):
    try:
        prev = json.load(open(OUT))
        refs["examples"] = {k: v for k, v in prev.get("examples", {}).items() if v.get("status") == "ok"}
        print(f"[gen] resuming: {len(refs['examples'])} entries kept")
    except Exception:
        pass

for path in EXAMPLES:
    if os.path.basename(path) in refs["examples"]:
        continue
    name = os.path.basename(path)
    entry = {"status": "failed"}
    try:
        mas = json.load(open(path))
        fd = mas["magnetic"]["core"]["functionalDescription"]
        shape = fd.get("shape")
        fam = shape.get("family") if isinstance(shape, dict) else None
        entry["shape_family"] = fam or "unknown"
    except Exception as e:
        entry["error"] = f"mas parse: {e}"
    t0 = time.time()
    try:
        extrapolate = name in EXTRAPOLATED_EXAMPLES
        command = [OMFEM] + (["--allow-material-data-extrapolation"] if extrapolate else []) + [path, OUT + ".scratch.msh"]
        r = subprocess.run(command, capture_output=True, text=True, timeout=TIMEOUT)
        if LOG_DIR:
            open(os.path.join(LOG_DIR, name + ".stderr"), "w").write(r.stderr)
        out = r.stdout.strip().splitlines()
        payload = None
        for line in reversed(out):
            line = line.strip()
            if line.startswith("{"):
                payload = json.loads(line); break
        if payload is None:
            entry["error"] = (r.stderr.strip().splitlines() or ["no json output"])[-1][:300]
        elif not payload.get("temperature_source"):
            entry["error"] = "no thermal result in solve output"
        elif extrapolate and not payload.get("material_data_extrapolations"):
            entry["error"] = "extrapolation was switched on but omfem_mas returned no MKF extrapolation warning"
        elif not extrapolate and "material_data_extrapolation_allowed" in payload:
            entry["error"] = "omfem_mas reports material data extrapolation on an example that did not ask for it"
        else:
            if extrapolate:
                entry["material_data_extrapolation"] = {"allowMaterialDataExtrapolation": True,
                                                        "extrapolated": payload["material_data_extrapolations"]}
            entry.update({
                "status": "ok",
                "ambient": payload["ambient_temperature"],
                "p_core": payload.get("P_core_iGSE", 0.0) + payload.get("P_core_eddy", 0.0),
                "p_cu": payload.get("P_cu", 0.0),
                "omfem_t_core": payload.get("temperature_core"),
                "omfem_t_winding": payload.get("temperature_winding"),
                "omfem_t_max": payload.get("temperature_max"),
            })
    except subprocess.TimeoutExpired:
        entry["error"] = f"timeout after {TIMEOUT}s"
    except Exception as e:
        entry["error"] = str(e)[:300]
    entry["seconds"] = round(time.time() - t0, 1)
    refs["examples"][name] = entry
    print(f"[gen] {name}: {entry['status']} ({entry['seconds']}s)" + (f" core={entry.get('omfem_t_core'):.1f} max={entry.get('omfem_t_max'):.1f} Pc={entry.get('p_core'):.3g} Pcu={entry.get('p_cu'):.3g}" if entry["status"]=="ok" else f" -- {entry.get('error','')[:120]}"), flush=True)
    json.dump(refs, open(OUT, "w"), indent=1)

ok = sum(1 for e in refs["examples"].values() if e["status"] == "ok")
print(f"[gen] DONE: {ok}/{len(refs['examples'])} examples produced references")
