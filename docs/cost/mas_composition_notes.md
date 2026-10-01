# MAS core materials without `materialComposition`: proposal notes

Input: `MAS/data/core_materials.ndjson` (in `.mkf-cost`). It has 70 records with no `materialComposition`: TDG 47, Würth Elektronik 11, Huoh Yow 5, Guangdong Ferrum 3, Gaotune 2, Shincore 2.
Allowed values: the enum in `schemas/magnetic/core/material.json` (MnZn, NiZn, MgZn, FeSiAl, FeSi, FeSiCr, FeNi, FeNiMo, FeMo, carbonylIron, iron, proprietary). No schema change is proposed, and no repo or MAS file was edited.
Everything was fetched between 2026-09-30 23:00 and 23:45 UTC. The raw pages and PDFs are in `~/OpenMagnetics/cost_evidence/mas_composition_evidence/`. Each CSV row names its evidence file(s); the first one holds the quoted text.

## Result

| status | count |
|---|---|
| confirmed | 40 |
| needs_review | 18 |
| blank | 12 |

| manufacturer | confirmed | needs_review | blank |
|---|---|---|---|
| TDG (47) | 35 | 8 | 4 |
| Würth Elektronik (11) | 2 | 8 | 1 |
| Guangdong Ferrum (3) | 3 | 0 | 0 |
| Shincore (2) | 0 | 2 (no enum value fits) | 0 |
| Huoh Yow (5) | 0 | 0 | 5 |
| Gaotune (2) | 0 | 0 | 2 |

## Per manufacturer

**TDG.** The metal powder core catalogue PDF (TDG download page, created 2024-05-29) heads each series page with its alloy, e.g. `材料/Material:铁硅铝系列 TMS SERIES`. Its intro reads "Fe-Si-Al (Sendust), Fe-Si, Fe-Ni (High Flux), FeNi composite powder core and FeSi composite powder core".
- Confirmed: TMS and TMSA are FeSiAl. TMF, TMFA, TMFB and TMFC are FeSi. TMH is FeNi. TMM is FeNiMo (铁镍钼). TMHG (超级铁镍, "super Fe-Ni") is FeNi.
- needs_review: TMHA and TMHB (铁镍复合, "Fe-Ni composite"). The catalogue does not say what the Ni-Fe powder is blended with. FeNi is the proposal; `proprietary` is the alternative.
- blank: TMSC 26/60 and TMFD 60/90. They are not in the 2024 catalogue. The 2025-11 automotive brochure (image-only, read visually) names them only in loss and DC-bias tables. The TDG mall has no item for either. Their names suggest FeSiAl and FeSi siblings, but that is not used.

**Würth Elektronik.** Most WE product pages give only a grade ("Material 3 W 800"), not a composition.
- Confirmed, 91101004 (MAS "1E04_620 91101004") → NiZn. Four maker pages agree. WE-CAR-TEC says "Core material: NiZn" and every article is Material 91101004. WE-AEFA says "Ferrite core made of NiZn" and every article is Material 1E04_620. WE-OEFA says "Material: NiZn", and its datasheet gives Material 91101004. WE-STAR-TEC says "Core material: NiZn" with grade 4 W 620, and design-kit datasheet 742711 gives 91101004.
- Confirmed, 91102001 (MAS "3W5000") → MnZn. WE-STAR-TEC-LFS says "Core material: MnZn" with grade "8 W 5000". Design-kit datasheet 742711 lists the same LFS articles as Material 91102001. The WE 2010 catalogue says WE-AFB LFS (8 W 5000) has a "Ferrite core made of MnZn". **Name issue to review:** MAS calls this grade 3W5000, but WE prints it as 8 W 5000.
- needs_review → NiZn: 91101007 (MAS "1E04_1521", i.e. 4 W 1500), 91101009/14/15 (3 W 800), 91101010 (7 W 700), 91101011 (7 W 850), 91101012 (7 W 380) and 91101013 (7 W 800). These grades appear only in WE-AFB, WE-SAFB, WE-TOF and WE-FLAT, whose WE pages and datasheets state no composition. The only evidence is at series level:
  - The WE eiSos 2018 catalogue (EMI suppression chart, pp. 12–13) colours the "WE-SFA | WE-FLAT | WE-TOF / WE-AFB | WE-SAFB | WE-RIB" bar green, and its legend has green = NiZn. This was read from the page image.
  - Mouser's product text says "The WE-AFB, WE-SAFB, and WE-TOF series all feature a ferrite core made of NiZn". Mouser is a distributor, not the maker.

  This needs review because the claim covers whole series rather than grades, it is colour-coded, and WE also sells MgZn cable ferrites (ANP115).
- blank: 91101006 (4 W 810). It appears only in WE-FLATF ("μi = 620 to 810"), which has no composition statement.
- Grade ↔ article mapping for these records comes from each record's own MAS `description` ("grade X, material article N").

**Guangdong Ferrum.** The earlier run found no composition, but it is there. On ferrum.cn/Material/index_2.html ("合金材料") the tabs 铁硅铬/铁硅/铁镍/铁镍钼/铁硅铝 map in order to five tables. The first tab, 铁硅铬材料 (FeSiCr), holds M800, M850, M600 and M350C, with links to exactly the MAS datasheet PDFs. The other tables hold G10–G40. I checked the tab-to-table order in the HTML. These rows are marked confirmed, but the link is structural (a tab), not a label on the row. The datasheet PDFs themselves state no composition.

**Shincore.** The page says "Material ： Fe-based nanocrystalline （ 1K107 ）". No enum value fits. The proposal is left empty and marked needs_review: either leave it unset, use `proprietary`, or raise a schema gap. No schema change is proposed.

**Huoh Yow.** The page gives only the range "The metal materials are FeSiCr、FeSi、FeSiAl & FeNi". The property table for MPT21, MSB40, MSS70, MST70 and MSS100 has no material row. A separate table does give composition, but only for other grades (MFN120 FeNi, MSL125 FeSiAl, MFS125 FeSi). The per-grade PDFs are image-only. I rendered them and read them visually: they show properties only. All five are blank.

**Gaotune.** The datasheet page is JS-rendered with two images, "Amorphous material" and "Nanocrystalline material". They give properties only, with no alloy and no AF/AN grade name. Any "Fe-based" statement would not map to the enum anyway. Both are blank.

## Checks against the earlier run's claims
- All TDG series claims were re-verified from the re-fetched PDF. TMSC and TMFD are still not in it.
- WE: 4 W 620 → NiZn and 8 W 5000 → MnZn were re-verified. I also added direct article-number links (91101004, 91102001) through WE-CAR-TEC, WE-OEFA and design kit 742711. The earlier "no statement" for other WE grades still holds for maker text specific to a grade.
- Huoh Yow and Gaotune: the earlier findings hold.
- Shincore: the earlier finding holds.
- Guangdong Ferrum: **corrected.** The earlier run found no composition; FeSiCr was found on the alloy-materials page.

## Method notes
- tesseract is not installed, so image-only PDFs and images were rendered with `pdftoppm` and read visually. Quotes from images say so in the CSV notes.
- Each `quoted_text` taken from a text layer was checked by script to appear verbatim in its first evidence file. The WE colour-chart row is described in words, not quoted.
- Keiro search was used for leads only: about 11 queries.
- Some evidence files are kept only to support negative findings, e.g. the WE series pages with no composition and the Huoh Yow PDFs. `mouser_742700790_WE-AFB_7W380.html` was saved through Playwright as a JSON-escaped HTML string.
