# catalog_inventory — config-driven core-catalog importer

Generalizes hephaestus' per-manufacturer `*Inventory` classes into **one engine + a registry**.
Adding a manufacturer is a `MANUFACTURERS` entry (a parametric source file + a few callables),
not a bespoke ~200-line class with a hardcoded, expiring scrape.

## What the engine does per part
- **Shape** — resolved to an existing MAS/PyOpenMagnetics shape by effective **Ae/Le**
  (`match_shape_by_effective`) *or* directly by an IEC size **name** the source provides
  (`shape_name`). Toroids are emitted from OD/ID/HT at full precision. Any referenced shape not
  yet in MAS is emitted to the shapes output (canonicalized), so the result is self-consistent.
- **Coating** — mapped to the MAS coating field (`epoxy` / `parylene` / …).
- **Gapping** — three modes: none, **gap-from-length**, or **gap-from-target-AL** (the AL is hit
  by sweeping the gap through `calculate_inductance_from_number_turns_and_gapping`, then
  re-verified within 4 %). Canonical 1-subtractive + 2-residual @ 5 µm topology.
- **Naming** — `"<shape> - [<coating> coated - ]<material> - {Ungapped|Gapped X.XXX mm}"`.

## Stale-cmrc workaround
The installed `PyOpenMagnetics` wheel ships a cmrc material snapshot that predates recently
added grades (`CORE_MATERIAL_NOT_FOUND`). The engine loads materials from the live
`MAS/core_materials.ndjson` and **embeds the full material object inline**, so MKF never hits the
stale DB. Point it elsewhere with `MAS_DATA_DIR`.

## Run
```
python3 catalog_inventory.py <manufacturer|all> [limit]
```
Writes `<key>_shapes.ndjson` + `<key>_cores.ndjson` (or `inventory_*` for `all`). A manufacturer
whose source file is missing is skipped with a warning — fetch its parametric list first.

## Add a manufacturer
Append to `MANUFACTURERS`:
```python
"tdg": {
    "manufacturer": "TDG",                 # must match manufacturerInfo.name in core_materials
    "source": "tdg_cores.csv",             # parametric core list (see "Sourcing" below)
    "datasheet": "https://…",
    "reader": _csv_reader,                  # -> DataFrame
    "map": _tdg_map,                        # row -> {part, family, (od/id/ht | shape_name | ae+le),
                                            #         material_name|mat_letter, coating?, gap_mm?, target_al?}
    "resolve_material": _by_material_name,  # (row, mfr_material_names) -> MAS material name | None
},
```
The engine handles the rest (shape resolution, gapping, naming, emission, validation).

## Sourcing the parametric list (the only per-vendor work left)
Two proven patterns:
1. **Local vendor export** — a part-finder xlsx/csv (e.g. Magnetics' Advanced-Search DB). Robust;
   no scraping. Wire `_csv_reader` / `pandas.read_excel`.
2. **Live parametric API via the Playwright-MCP trick** — navigate to the vendor's results page
   to establish the session, then in-page `fetch` the JSON/CSV API (passes Akamai/Cloudflare that
   block raw curl). Used for TDK's `list.csv`. This replaces hephaestus' stale hardcoded cookies.

### Coverage status
- **Existing-cores classes (hephaestus):** Ferroxcube, TDK, Magnetics, Fair-Rite, Micrometals.
- **catalog_inventory configs:** Magnetics (xlsx), TDK (list.csv via MCP) — both verified
  (TDK: 14/63 generated cores exact-name-match existing MAS cores).
- **Still need a parametric list** (materials present in MAS, cores absent): ACME, AT&M, Chang
  Sung, Cosmo, DMEGC, Gaotune, JFE, KDM, Magnetec, Metglas, Nicera, Poco, Proterial, Samwha,
  Sinomag, TDG, VAC, Würth. Each is a config entry once its core list is sourced (pattern 1 or 2).
  Most are Chinese ferrite/powder makers whose catalogs are PDF/site-based — the sourcing is a
  per-vendor discovery task, not an engine change.

## Stocked side (hermes)
`hermes.add_distributor_by_mpn(mpn, …)` matches a distributor product's Manufacturer Part Number
straight against the inventory's `manufacturerInfo.reference`. It is wired into `DigikeyStocker`
as the first path, so **any manufacturer covered here becomes stockable from any distributor with
no bespoke `process_<mfr>_product` parser**; the per-manufacturer parsers remain as a fallback.
