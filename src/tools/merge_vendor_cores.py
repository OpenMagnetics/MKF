"""merge_vendor_cores.py — merge catalog_inventory vendor output into MAS data, safely.

Appends <key>_cores.ndjson / <key>_shapes.ndjson (produced by catalog_inventory.py) into the live
MAS data/cores.ndjson and data/core_shapes.ndjson, with full guardrails:
  * schema-validate every new core + shape (MAS + PEAS) — abort on any failure
  * referential integrity: every new core's shape+material must resolve — abort on any dangling
  * skip cores/shapes whose name already exists (idempotent; never duplicate a name)
  * dry-run by default; pass --apply to write

Usage:  python3 merge_vendor_cores.py <key> [<key> ...] [--apply]
"""
import json, os, sys, glob, pathlib
from referencing import Registry, Resource
from jsonschema import Draft202012Validator

HERE = pathlib.Path(__file__).parent.resolve()
MAS = pathlib.Path(os.environ.get("MAS_DATA_DIR") or (HERE.parent.parent.parent / "MAS" / "data"))
SCHEMAS = MAS.parent / "schemas"
PEAS = HERE.parent.parent / "PEAS" / "schemas"


def _registry():
    reg = Registry()
    for root in (SCHEMAS, PEAS):
        for p in glob.glob(str(root) + "/**/*.json", recursive=True):
            try:
                sc = json.load(open(p))
            except Exception:
                continue
            if sc.get("$id"):
                reg = reg.with_resource(sc["$id"], Resource.from_contents(sc))
    return reg


def _load_names(path, key="name"):
    out = {}
    if not os.path.exists(path):
        return out
    for l in open(path):
        s = l.strip()
        if not s or s.startswith("version"):
            continue
        try:
            o = json.loads(s)
        except Exception:
            continue
        if key in o:
            out[o[key]] = o
    return out


def main(keys, apply):
    reg = _registry()
    core_v = Draft202012Validator(json.load(open(SCHEMAS / "magnetic" / "core.json")), registry=reg)
    shape_v = Draft202012Validator(json.load(open(SCHEMAS / "magnetic" / "core" / "shape.json")), registry=reg)

    existing_cores = _load_names(MAS / "cores.ndjson")
    existing_shapes = _load_names(MAS / "core_shapes.ndjson")
    materials = set(_load_names(MAS / "core_materials.ndjson"))

    add_cores, add_shapes = {}, {}
    for key in keys:
        cores = [json.loads(l) for l in open(HERE / f"{key}_cores.ndjson")] if os.path.exists(HERE / f"{key}_cores.ndjson") else []
        shapes = [json.loads(l) for l in open(HERE / f"{key}_shapes.ndjson")] if os.path.exists(HERE / f"{key}_shapes.ndjson") else []
        for s in shapes:
            if s["name"] not in existing_shapes and s["name"] not in add_shapes:
                add_shapes[s["name"]] = s
        for c in cores:
            if c["name"] not in existing_cores and c["name"] not in add_cores:
                add_cores[c["name"]] = c
        print(f"[{key}] cores={len(cores)} shapes={len(shapes)}")

    # validate
    errs = []
    for s in add_shapes.values():
        if next(shape_v.iter_errors(s), None):
            errs.append(("shape-schema", s["name"]))
    all_shape_names = set(existing_shapes) | set(add_shapes)
    for c in add_cores.values():
        if next(core_v.iter_errors(c), None):
            errs.append(("core-schema", c["name"]))
        sh = c["functionalDescription"]["shape"]
        sh = sh if isinstance(sh, str) else sh.get("name")
        mat = c["functionalDescription"]["material"]
        mat = mat if isinstance(mat, str) else mat.get("name")
        if sh not in all_shape_names:
            errs.append(("dangling-shape", f'{c["name"]} -> {sh}'))
        if mat not in materials:
            errs.append(("dangling-material", f'{c["name"]} -> {mat}'))

    print(f"\nNEW shapes: {len(add_shapes)}   NEW cores: {len(add_cores)}   errors: {len(errs)}")
    for e in errs[:20]:
        print("   ERR", e)
    if errs:
        print("ABORT — fix errors before merging."); return 1
    if not apply:
        print("dry-run OK. Re-run with --apply to write."); return 0

    with open(MAS / "core_shapes.ndjson", "a") as f:
        for s in add_shapes.values():
            f.write(json.dumps(s, ensure_ascii=False) + "\n")
    with open(MAS / "cores.ndjson", "a") as f:
        for c in add_cores.values():
            f.write(json.dumps(c, ensure_ascii=False) + "\n")
    print(f"APPLIED: +{len(add_shapes)} shapes, +{len(add_cores)} cores")
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    apply = "--apply" in args
    keys = [a for a in args if not a.startswith("--")]
    sys.exit(main(keys, apply))
