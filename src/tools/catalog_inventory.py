"""catalog_inventory.py — config-driven manufacturer core-catalog importer.

A generalization of hephaestus' per-manufacturer *Inventory classes. Instead of a bespoke
~200-line class per manufacturer (each with a fragile hardcoded scrape), a manufacturer is a
single entry in the MANUFACTURERS registry: a parametric source file (xlsx/csv exported from
the vendor's part finder), a column map, a family map, a material resolver, and coating/gap
rules. One generic engine then builds cores_inventory.ndjson records for it.

What the engine does per part (all via PyOpenMagnetics, like the hephaestus base Stocker):
  * SHAPE: resolve to an existing MAS core shape by family + effective Ae/Le (the same idea as
    Stocker.find_shape_closest_effective_parameters); toroids are emitted directly from OD/ID/HT.
  * COATING: map the vendor's coating column to the MAS coating field.
  * GAPPING: none, or gap-from-length, or gap-from-target-AL — the AL case is solved by
    sweeping the gap through calculate_inductance_from_number_turns_and_gapping until the
    computed AL matches (identical to Stocker.get_gapping), then re-verified.
  * NAMING: "<shape> - [<coating> coated - ]<material> - {Ungapped|Gapped X.XXX mm}".

Stale-embedded-DB workaround: the installed PyOpenMagnetics ships a cmrc snapshot that predates
recently-added materials (CORE_MATERIAL_NOT_FOUND). We load materials from the live MAS
core_materials.ndjson and EMBED the full material object inline so MKF never hits the stale DB.

Adding a manufacturer: append a MANUFACTURERS entry. No new class, no scrape.
"""
import json, math, os, pathlib, re, bisect
import pandas
import PyOpenMagnetics

HERE = pathlib.Path(__file__).parent.resolve()
# live MAS data (sibling checkout preferred over the wheel's stale cmrc snapshot)
MAS_DATA = os.environ.get("MAS_DATA_DIR") or str(HERE.parent.parent.parent / "MAS" / "data")


# ----------------------------------------------------------------- material resolution
def load_materials():
    """name -> full material object, and manufacturer -> [names], from live MAS data."""
    by_name, by_mfr = {}, {}
    path = os.path.join(MAS_DATA, "core_materials.ndjson")
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("version"):
                continue
            o = json.loads(line)
            by_name[o["name"]] = o
            by_mfr.setdefault(o.get("manufacturerInfo", {}).get("name"), []).append(o["name"])
    return by_name, by_mfr


# ----------------------------------------------------------------- MKF helpers (inline material)
_shape_cache = {}

def _shape(name):
    if name not in _shape_cache:
        _shape_cache[name] = PyOpenMagnetics.find_core_shape_by_name(name)
    return _shape_cache[name]

def effective_parameters(shape_name):
    """Ae, Le, Ve for a shape (material-independent geometry)."""
    core = {"name": "x", "functionalDescription": {
        "type": "toroidal" if shape_name.split()[0] in ("T", "R") else "twoPieceSet",
        "material": "3C90", "shape": _shape(shape_name), "gapping": [], "numberStacks": 1}}
    d = PyOpenMagnetics.calculate_core_data(core, False)
    d = json.loads(d) if isinstance(d, str) else d
    ep = d["processedDescription"]["effectiveParameters"]
    return ep["effectiveArea"], ep["effectiveLength"], ep["effectiveVolume"]

_family_shapes = {}
def shapes_of_family(family):
    if family not in _family_shapes:
        out = []
        for nm in PyOpenMagnetics.get_core_shape_names(True):
            if nm.split(" ")[0].lower() == family.lower() and "X" not in nm:
                try:
                    out.append((nm,) + effective_parameters(nm))
                except Exception:
                    pass
        _family_shapes[family] = out
    return _family_shapes[family]

def match_shape_by_effective(family, ae, le, tol=0.05, prefer=None):
    """Closest existing shape by |dAe|+|dLe| within tol on each (Stocker.find_shape_closest_*).

    When several shapes match within tolerance, prefer the variant THIS manufacturer already uses
    most in MAS (`prefer` = {shape_name: usage_count}). Sibling shape variants (e.g. pot cores
    'P 11/7' vs 'P 11/7/I') differ only by a feature dimension that barely moves Ae/Le, so a
    pure-distance match can pick the variant MAS does NOT canonically use for that part — creating
    a physical duplicate under an alternate shape name. Magnetics pot cores, for instance, are
    recorded '/I' 162x vs bare 17x; ranking by usage picks '/I'. Distance only breaks usage ties."""
    prefer = prefer or {}
    best, best_key = None, None  # key = (usage_count, -distance); maximise
    for nm, sae, sle, sve in shapes_of_family(family):
        if abs(sae - ae) / ae < tol and abs(sle - le) / le < tol:
            d = abs(sae - ae) / ae + abs(sle - le) / le
            key = (prefer.get(nm, 0), -d)
            if best_key is None or key > best_key:
                best, best_key = nm, key
    return best

_al_coil = {"bobbin": "Dummy", "functionalDescription": [
    {"name": "pri", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Dummy"}]}
_al_op = {"conditions": {"ambientTemperature": 25.0}, "excitationsPerWinding": [
    {"frequency": 10000, "current": {"waveform": {"data": [-0.001, 0.001, -0.001], "time": [0, 5e-5, 1e-4]},
     "processed": {"label": "Triangular", "peakToPeak": 0.002, "offset": 0, "dutyCycle": 0.5}}}]}
_al_models = {"gapReluctance": "ZHANG"}

def _core_with_gap(shape_name, material_obj, gap, num_residual):
    gapping = [{"type": "subtractive", "length": gap}] + [{"type": "residual", "length": 5e-6}] * num_residual
    core = {"name": "x", "functionalDescription": {
        "type": "toroidal" if shape_name.split()[0] in ("T", "R") else "twoPieceSet",
        "material": material_obj, "shape": _shape(shape_name), "gapping": gapping, "numberStacks": 1}}
    proc = PyOpenMagnetics.calculate_core_data(core, False)          # processes gapping -> sets gap areas
    return json.loads(proc) if isinstance(proc, str) else proc

def al_at_gap(shape_name, material_obj, gap, num_residual=2):
    core = _core_with_gap(shape_name, material_obj, gap, num_residual)
    L = PyOpenMagnetics.calculate_inductance_from_number_turns_and_gapping(core, _al_coil, _al_op, _al_models)
    L = L["magnetizingInductance"]["nominal"] if isinstance(L, dict) else L
    return L * 1e9   # nH

def gap_for_al(shape_name, material_obj, target_al, num_residual=2):
    """Bisection in log-gap space; returns (gap_m, achieved_al) or (None, None)."""
    lo, hi = 2e-6, 6e-3
    try:
        al_lo, al_hi = al_at_gap(shape_name, material_obj, lo, num_residual), al_at_gap(shape_name, material_obj, hi, num_residual)
    except Exception:
        return None, None
    if not (al_hi <= target_al <= al_lo):
        return None, None
    mid = a = None
    for _ in range(30):
        mid = (lo * hi) ** 0.5
        a = al_at_gap(shape_name, material_obj, mid, num_residual)
        if abs(a - target_al) / target_al < 0.01:
            return mid, a
        if a > target_al: lo = mid
        else: hi = mid
    return mid, a


# ----------------------------------------------------------------- naming
def core_name(shape, material, gapping, coating=None):
    if not gapping:
        return f"{shape} - {material} - Ungapped" if coating is None else f"{shape} - {coating} coated - {material} - Ungapped"
    nsub = len([g for g in gapping if g["type"] == "subtractive"])
    tag = "Gapped" if nsub == 1 else "Distributed gapped"
    return f"{shape} - {material} - {tag} {gapping[0]['length'] * 1000:.3f} mm"

def adaptive_round(v):
    v = float(v)
    return round(v) if v >= 10 else round(v, 1) if v >= 1 else round(v, 2)

def _fmt(v):
    """Toroid dims keep their real precision (12.7, not 13) so the emitted shape is exact."""
    return ("%g" % round(float(v), 3))

def load_shape_usage_by_cores(manufacturer=None):
    """{shape_name: usage_count} across existing MAS cores, optionally restricted to one
    manufacturer. This is the canonical variant MAS uses for each part; match_shape_by_effective
    ranks Ae/Le matches by it so we don't emit shadow-duplicate cores under a sibling shape name
    (e.g. Magnetics pot cores are '/I' 162x vs bare 17x)."""
    from collections import Counter
    usage = Counter()
    path = os.path.join(MAS_DATA, "cores.ndjson")
    if not os.path.exists(path):
        return usage
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("version"):
                continue
            try:
                o = json.loads(line)
            except Exception:
                continue
            if manufacturer and (o.get("manufacturerInfo") or {}).get("name") != manufacturer:
                continue
            sh = o.get("functionalDescription", {}).get("shape")
            nm = sh if isinstance(sh, str) else (sh.get("name") if isinstance(sh, dict) else None)
            if nm:
                usage[nm] += 1
    return usage


def load_existing_shapes():
    """Existing MAS core shapes: family -> [(name, [A,B,C mm])] for dedup / numeric matching."""
    from collections import defaultdict
    fams = defaultdict(list)
    path = os.path.join(MAS_DATA, "core_shapes.ndjson")
    if not os.path.exists(path):
        return fams
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("version"):
                continue
            o = json.loads(line)
            dims = o.get("dimensions", {})
            def nom(k):
                d = dims.get(k)
                return (d.get("nominal") if isinstance(d, dict) else d) if d else None
            fams[o.get("family")].append((o["name"], [nom("A"), nom("B"), nom("C")]))
    return fams

def toroid_shape_record(od, idd, ht):
    return {"magneticCircuit": "closed", "type": "standard", "family": "t",
            "name": f"T {_fmt(od)}/{_fmt(idd)}/{_fmt(ht)}",
            "dimensions": {"A": {"nominal": round(od / 1000, 6)},
                           "B": {"nominal": round(idd / 1000, 6)},
                           "C": {"nominal": round(ht / 1000, 6)}}}


# ----------------------------------------------------------------- generic build
def build_inventory(mfr_key, limit=None):
    cfg = MANUFACTURERS[mfr_key]
    by_name, by_mfr = load_materials()
    mfr_mats = by_mfr.get(cfg["manufacturer"], [])
    df = cfg["reader"](cfg)
    existing_shapes = load_existing_shapes()
    # usage-weighted, manufacturer-specific: pick the shape variant THIS mfr already uses most
    shape_usage = load_shape_usage_by_cores(cfg["manufacturer"]) or load_shape_usage_by_cores()
    new_shapes = {}   # name -> shape record emitted for toroid sizes not already in MAS
    cores, stats = {}, {"rows": 0, "no_material": 0, "no_shape": 0, "made": 0, "gap_fail": 0, "new_shapes": 0}

    existing_names = {n for lst in existing_shapes.values() for n, _ in lst}

    def clean_dim(d):
        return {k: v for k, v in d.items() if v is not None} if isinstance(d, dict) else d

    def ensure_concentric_shape(name):
        """Resolve a matched or named shape to its CANONICAL name; emit its record if that
        canonical name is not yet in MAS. Returns the canonical name to reference from the core
        (PyOM lists short aliases like 'E 42/15' that canonicalize to 'E 42/21/15'), or None if
        the name is not a known shape."""
        try:
            s = PyOpenMagnetics.find_core_shape_by_name(name)
        except Exception:
            return None
        if not s or not s.get("dimensions"):
            return None
        canonical = s.get("name") or name
        if canonical in existing_names or canonical in new_shapes:
            return canonical
        rec = {k: s[k] for k in ("magneticCircuit", "type", "family", "aliases", "name") if k in s}
        rec.setdefault("name", canonical)
        rec.setdefault("family", canonical.split(" ")[0].lower())
        rec["dimensions"] = {k: clean_dim(v) for k, v in s.get("dimensions", {}).items() if v}
        new_shapes[canonical] = rec
        existing_names.add(canonical)
        stats["new_shapes"] += 1
        return canonical

    def toroid_shape(od, idd, ht):
        """Return an existing toroid shape name matched by numeric dims (0.06 mm tol), else emit one."""
        for name, (a, b, c) in ((n, (d[0], d[1], d[2])) for n, d in existing_shapes.get("t", [])):
            if None in (a, b, c):
                continue
            if abs(a * 1000 - od) < 0.06 and abs(b * 1000 - idd) < 0.06 and abs(c * 1000 - ht) < 0.06:
                return name
        rec = toroid_shape_record(od, idd, ht)
        if rec["name"] not in new_shapes:
            new_shapes[rec["name"]] = rec
            existing_shapes["t"].append((rec["name"], [od / 1000, idd / 1000, ht / 1000]))
            stats["new_shapes"] += 1
        return rec["name"]

    for _, row in df.iterrows():
        stats["rows"] += 1
        if limit and stats["made"] >= limit:
            break
        r = cfg["map"](row, cfg)
        if r is None:
            continue
        material = cfg["resolve_material"](r, mfr_mats)
        if material is None or material not in by_name:
            stats["no_material"] += 1
            continue
        mat_obj = by_name[material]
        family = r["family"]
        # --- shape ---
        if family in ("T", "R"):
            shape = toroid_shape(float(r["od"]), float(r["id"]), float(r["ht"]))
        elif r.get("shape_name"):        # source gives the shape designation directly (e.g. TDK IEC size)
            shape = ensure_concentric_shape(r["shape_name"])
        else:
            shape = None
            if r.get("ae") and r.get("le"):
                shape = match_shape_by_effective(family, r["ae"] / 1e6, r["le"] / 1e3, prefer=shape_usage)
            if shape is not None:
                shape = ensure_concentric_shape(shape)
        if shape is None:
            stats["no_shape"] += 1   # non-toroid size with no existing MAS shape (needs A-F dims -> not fabricated)
            continue
        coating = r.get("coating")
        # --- gapping ---
        gapping = []
        if r.get("gap_mm"):
            gapping = [{"type": "subtractive", "length": round(r["gap_mm"] / 1000, 8)},
                       {"type": "residual", "length": 5e-6}, {"type": "residual", "length": 5e-6}]
        elif r.get("target_al"):
            gap, al_fit = gap_for_al(shape, mat_obj, r["target_al"])
            if gap is None or abs(al_fit - r["target_al"]) / r["target_al"] > 0.04:
                stats["gap_fail"] += 1
                continue
            gapping = [{"type": "subtractive", "length": round(gap, 8)},
                       {"type": "residual", "length": 5e-6}, {"type": "residual", "length": 5e-6}]
        name = core_name(shape, material, gapping, coating)
        if name in cores:
            continue
        fd = {"type": "toroidal" if family in ("T", "R") else "twoPieceSet",
              "material": material, "shape": shape, "gapping": gapping, "numberStacks": 1}
        if coating is not None:
            fd["coating"] = coating
        cores[name] = {"name": name, "manufacturerInfo": {"name": cfg["manufacturer"], "status": r.get("status", "production"),
                       "reference": r.get("part"), "datasheetUrl": cfg.get("datasheet")},
                       "distributorsInfo": [], "functionalDescription": fd}
        stats["made"] += 1
    return list(new_shapes.values()), list(cores.values()), stats


# ----------------------------------------------------------------- manufacturer registry
def _magnetics_reader(cfg):
    return pandas.read_excel(cfg["source"], engine="openpyxl")

def _magnetics_map(row, cfg):
    shape = str(row.get("Shape") or "")
    part = str(row.get("Sales_Part_Number") or "")
    mat_letter = str(row.get("Material") or "").strip()
    if not shape or not mat_letter:
        return None
    r = {"part": part, "mat_letter": mat_letter, "status": "production"}
    if "oroid" in shape:
        r.update(family="T", od=row["OD__A_mm"], id=row["ID__B_mm"], ht=row["HT__C_mm"])
        cmap = {"Coated": "epoxy", "Parylene": "parylene", "Uncoated": None}
        r["coating"] = next((v for k, v in cmap.items() if k in shape), None)
        if "Nylon" in shape or "X (" in shape:   # discontinued
            return None
    else:
        fam = {"E": "E", "ETD": "ETD", "PQ": "PQ", "RM": "RM", "EP": "EP", "ER": "ER", "EFD": "EFD",
               "EC": "EC", "EER": "EER", "Pot Core": "P", "U": "U", "UR": "UR"}.get(shape)
        if fam is None:
            return None
        r.update(family=fam, ae=row.get("Ae_mm2"), le=row.get("Le_mm"))
        gapped = str(row.get("Gapped") or "").lower() == "gapped"
        if gapped:
            gin = row.get("Gap_length_inches")
            if gin and float(gin) > 0:
                r["gap_mm"] = float(gin) * 25.4
            elif row.get("AL"):
                r["target_al"] = float(row["AL"])
            else:
                return None
    return r

def _magnetics_material(r, mfr_mats):
    # ferrite grades are single letters (C,E,F,J,L,M,P,R,T,V,W); accept if that letter is a known material
    return r["mat_letter"] if r["mat_letter"] in mfr_mats else None

def _csv_reader(cfg):
    return pandas.read_csv(cfg["source"])

def _tdk_map(row, cfg):
    """TDK ferrite-core CSV (pulled from product.tdk.com/pdc_api .../list.csv). Shape designation
    is the IEC 'Size' column; gap is given directly; material is 'Material Name'."""
    shape_col = str(row.get("Core Shape") or "")
    size = row.get("Size (IEC62317) A x B x C")
    size = None if (size is None or (isinstance(size, float) and math.isnan(size))) else str(size).strip()
    part = str(row.get("Part No.") or "").split(" (")[0].strip()
    material = str(row.get("Material Name") or "").strip()
    if not part or not material:
        return None
    r = {"part": part, "material_name": material, "status": "production"}
    is_toroid = "Toroid" in shape_col or "(T)" in shape_col or (size and size[0] in "R")
    if is_toroid:
        try:
            r.update(family="T", od=float(row["Dimm. A"]), id=float(row["Dimm. B"]), ht=float(row["Dimm. C"]),
                     coating="epoxy" if "EMC" in shape_col else None)
        except (TypeError, ValueError):
            return None
    else:
        if not size:
            return None
        r.update(family=size.split()[0], shape_name=size)
    gap = row.get("Air Gap / mm")
    try:
        if gap is not None and float(gap) > 0:
            r["gap_mm"] = float(gap)
    except (TypeError, ValueError):
        pass
    return r

def _by_material_name(r, mfr_mats):
    return r["material_name"] if r["material_name"] in mfr_mats else None


def _normalized_map(row, cfg):
    """Map for a pre-normalized vendor CSV produced by an extraction script. Recognised columns:
    part, family, (od,id,ht | shape_name | ae,le), material_name, coating, gap_mm, target_al.
    Extraction (PDF/HTML table parsing) lives in the per-vendor scripts under vendor_cores/ so the
    engine stays source-agnostic; the CSV is the contract."""
    def val(k):
        v = row.get(k)
        if v is None:
            return None
        try:
            if isinstance(v, float) and math.isnan(v):
                return None
        except (TypeError, ValueError):
            pass
        s = str(v).strip()
        return s if s and s.lower() != "nan" else None
    part = val("part"); material = val("material_name"); family = val("family")
    if not part or not material or not family:
        return None
    r = {"part": part, "material_name": material, "family": family, "status": "production"}
    if family in ("T", "R"):
        for k in ("od", "id", "ht"):
            if val(k) is None:
                return None
        r.update(od=float(row["od"]), id=float(row["id"]), ht=float(row["ht"]))
    elif val("shape_name"):
        r["shape_name"] = val("shape_name")
    elif val("ae") and val("le"):
        r.update(ae=float(row["ae"]), le=float(row["le"]))
    else:
        return None
    if val("coating"):
        r["coating"] = val("coating")
    if val("gap_mm"):
        r["gap_mm"] = float(row["gap_mm"])
    elif val("target_al"):
        r["target_al"] = float(row["target_al"])
    return r


def _norm_cfg(manufacturer, env_var, default_file, datasheet):
    """Build a MANUFACTURERS entry for a vendor sourced from a normalized CSV."""
    return {"manufacturer": manufacturer,
            "source": os.environ.get(env_var, str(HERE / "vendor_cores" / default_file)),
            "datasheet": datasheet, "reader": _csv_reader,
            "map": _normalized_map, "resolve_material": _by_material_name}


MANUFACTURERS = {
    "tdk": {
        "manufacturer": "TDK",
        "source": os.environ.get("TDK_CSV", str(HERE / "tdk_cores.csv")),
        "datasheet": "https://product.tdk.com/en/search/ferrite/ferrite/ferrite-core/list",
        "reader": _csv_reader,
        "map": _tdk_map,
        "resolve_material": _by_material_name,
    },
    "magnetics": {
        "manufacturer": "Magnetics",
        "source": os.environ.get("MAGNETICS_XLSX",
                  "/mnt/c/Users/Alfonso/Downloads/Advanced+Search+Database+Ferrites+Updated+w+Hardware.xlsx"),
        "datasheet": "https://www.mag-inc.com/Media/Magnetics/File-Library/Product%20Literature/Ferrite%20Literature/Magnetics-2021-Ferrite-Catalog.pdf",
        "reader": _magnetics_reader,
        "map": _magnetics_map,
        "resolve_material": _magnetics_material,
    },
    # Vendors sourced from a normalized CSV (extraction script -> vendor_cores/<file>.csv).
    "atm": _norm_cfg("AT&M", "ATM_CSV", "atm_cores.csv",
                     "https://www.antai-emarketing.com/transformer-cores-7/"),
    "proterial": _norm_cfg("Proterial", "PROTERIAL_CSV", "proterial_cores.csv",
                           "https://www.india.proterial.com/products/"),
    "metglas": _norm_cfg("Metglas", "METGLAS_CSV", "metglas_cores.csv",
                         "https://metglas.com/"),
    "vac": _norm_cfg("Vacuumschmelze", "VAC_CSV", "vac_cores.csv",
                     "https://vacuumschmelze.com/products/inductive-components-and-cores/amorphous-and-nanocrystalline-cores"),
    "cosmo": _norm_cfg("Cosmo Ferrites", "COSMO_CSV", "cosmo_cores.csv",
                       "https://www.cosmoferrites.com/uploads/Catalogue_2026.pdf"),
    "cosmoconcentric": _norm_cfg("Cosmo Ferrites", "COSMO_CONC_CSV", "cosmo_concentric_cores.csv",
                                 "https://www.cosmoferrites.com/uploads/Catalogue_2026.pdf"),
    "dmegc": _norm_cfg("DMEGC", "DMEGC_CSV", "dmegc_cores.csv",
                       "https://www.dmegc.com.cn/product/9.html?lang=en"),
    "tdg": _norm_cfg("TDG", "TDG_CSV", "tdg_cores.csv",
                     "https://www.tdgcore.com/en/technical-support/data-download"),
    "tdgpowder": _norm_cfg("TDG", "TDG_POWDER_CSV", "tdg_powder_cores.csv",
                           "https://www.tdgcore.com/en/technical-support/data-download"),
    "poco": _norm_cfg("Poco", "POCO_CSV", "poco_cores.csv",
                      "https://en.pocomagnetic.com/products.html"),
    "changsung": _norm_cfg("Chang Sung", "CSC_CSV", "csc_cores.csv",
                           "https://changsung.com/cores?lang=en"),
    "samwha": _norm_cfg("Samwha", "SAMWHA_CSV", "samwha_cores.csv",
                        "http://www.samwha.com/electronics/eng/"),
    "kdm": _norm_cfg("KDM", "KDM_CSV", "kdm_cores.csv",
                     "https://www.kdm-mag.com/download.html"),
    "acme": _norm_cfg("ACME", "ACME_CSV", "acme_cores.csv",
                      "https://www.acme-ferrite.com.tw/"),
    "sinomag": _norm_cfg("Sinomag", "SINOMAG_CSV", "sinomag_cores.csv",
                         "https://sinomagtech.com/en/download/"),
    "nicera": _norm_cfg("Nicera", "NICERA_CSV", "nicera_cores.csv",
                        "https://www.nicera.co.jp/en/"),
    "magnetec": _norm_cfg("Magnetec", "MAGNETEC_CSV", "magnetec_cores.csv",
                          "https://www.magnetec.de/en/product/standard-cores/"),
    # To add a manufacturer, append an entry with its parametric source file + these callables.
}


def build_all(keys=None, limit=None):
    """Build every configured (or the given) manufacturer and merge shapes + cores, deduped by
    name. Drop-in for hephaestus: writes cores_inventory.ndjson / shapes_inventory.ndjson.
    A manufacturer whose source file is missing is skipped with a warning (its data must be
    fetched first, e.g. the TDK CSV via the Playwright-MCP pull)."""
    keys = keys or list(MANUFACTURERS)
    all_shapes, all_cores = {}, {}
    for key in keys:
        src = MANUFACTURERS[key].get("source")
        if src and not os.path.exists(src):
            print(f"[skip] {key}: source not found ({src}) — fetch it first")
            continue
        shapes, cores, stats = build_inventory(key, limit=limit)
        for s in shapes:
            all_shapes.setdefault(s["name"], s)
        for c in cores:
            all_cores.setdefault(c["name"], c)
        print(f"[{key}] {stats}")
    return list(all_shapes.values()), list(all_cores.values())


if __name__ == "__main__":
    import sys
    key = sys.argv[1] if len(sys.argv) > 1 else "magnetics"
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    if key == "all":
        shapes, cores = build_all(limit=lim)
        base = "inventory"
    else:
        shapes, cores, stats = build_inventory(key, limit=lim)
        print("stats:", stats)
        base = key
    shp = HERE / f"{base}_shapes.ndjson"
    with open(shp, "w") as f:
        for s in shapes:
            f.write(json.dumps(s, ensure_ascii=False) + "\n")
    out = HERE / f"{base}_cores.ndjson"
    with open(out, "w") as f:
        for c in cores:
            f.write(json.dumps(c, ensure_ascii=False) + "\n")
    print(f"wrote {len(shapes)} new shapes -> {shp}")
    print(f"wrote {len(cores)} cores -> {out}")
