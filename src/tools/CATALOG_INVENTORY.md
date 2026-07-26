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

### Coverage status (2026-07-26: 17 manufacturers, ~18,085 cores)
- **Imported via catalog_inventory:** Magnetics (xlsx), TDK (list.csv via MCP), and the 2026-07
  batch — AT&M, Cosmo, DMEGC, Chang Sung, Magnetec, POCO, Sinomag, TDG (ferrite + metal-powder),
  ACME, KDM, Proterial/Metglas (AMCC/F3CC cut C-cores). Plus the hephaestus classes (Ferroxcube,
  Fair-Rite, Micrometals). Extraction per vendor done by subagents into normalized CSVs; validated
  + merged centrally via `merge_vendor_cores.py`.
- **Import rules that held:** powder/nanocryst materials encode permeability in the MAS name
  (toroid×material); use BARE magnetic dims, never case/finished (VAC deferred — publishes only
  case dims); concentric shapes matched by Ae/Le, proprietary sizes with no MKF match are dropped
  (never fabricated); MKF-missing families skipped per tickets #263-#277; cut C-cores map by name
  to MKF "C <AMCC-designation>" shapes.
- **Deferred / still open:** VAC (case-dims only), Samwha (its MAS PL grades are ferrite, need the
  ferrite catalog not the powder one), Proterial Microlite distributed-gap toroids (per-permeability
  variant materials not in MAS), Würth (no bare power cores), Gaotune (no parametric data),
  Shandong Jianuo (no public catalog), Höganäs (SMC material, not cores). Non-toroid special shapes
  and MKF-missing families await the MKF engine tickets.

## Stocked side (hermes)
`hermes.add_distributor_by_mpn(mpn, …)` matches a distributor product's Manufacturer Part Number
straight against the inventory's `manufacturerInfo.reference`. It is wired into `DigikeyStocker`
as the first path, so **any manufacturer covered here becomes stockable from any distributor with
no bespoke `process_<mfr>_product` parser**; the per-manufacturer parsers remain as a fallback.
