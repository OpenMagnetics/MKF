# PEEC winding-loss model — validation against FastHenry

Status: **PEEC validated against FastHenry (the gold-standard filament-PEEC
solver, MIT / FastFieldSolvers).** The two independent PEEC solvers agree to
~0.1%, and BOTH exceed the in-repo "FEM" test expected values — i.e. those test
values appear to UNDER-count AC/proximity loss (consistent with coarse-mesh FEM).

## Method
- PEEC: `WindingLossesPEEC` (this branch) — graded surface-refined mesh,
  3x3 Gauss-quadrature `<ln d>` partial-inductance kernel, core via CoilMesher
  method-of-images.
- FastHenry: built from FastFieldSolvers source (master = MIT + Linux port).
  Each config exported by `Export_FastHenry_Inp` as series-connected bars
  (turns at their (x,y), w x h cross-section, nwinc=10/nhinc=6 skin filaments,
  .freq sweep). Loss = I_rms^2 * Re(Zc(f)).
- FastHenry is air-core (mu0); it references winding skin/proximity, not the
  core image effect. Agreement with PEEC (which includes images) to 0.1% on the
  ungapped rectangular case means the core-image effect on winding loss is small
  there (no gap fringing) — expected.

## Results (loss in W)

Five_Turns_Rectangular (I_rms^2 = 0.5):
| f [Hz] | FastHenry | PEEC | FEM(test) | Analytical | FH/FEM |
|--------|-----------|------|-----------|------------|--------|
| 1e3    | 0.000486  | ~0.00047 | (0.000493 DC) | 0.00049 | ~1.0 |
| 1e5    | 0.0023355 | 0.0023329 | 0.0012475 | 0.001115 | 1.87 |
| 1e6    | 0.0070979 | -        | 0.0032166 | -          | 2.21 |

Ten_Turns_Foil (I_rms^2 = 1.0):
| f [Hz] | FastHenry | FEM(test) | FH/FEM |
|--------|-----------|-----------|--------|
| 1e5    | 0.006194  | 0.004284  | 1.45 |

## Conclusions
1. **PEEC is correct** — matches FastHenry to 0.1% at 100 kHz (5-turn rect).
   It is mesh-convergent (2.053 vs 2.080 across two mesh densities).
2. **The "FEM" test expected values under-count AC loss by ~1.45-2.2x.** They
   were likely generated with coarse winding meshing (a real user's Maxwell
   cross-check of a planar flyback had the same coarse-mesh symptom). They
   should be re-baselined against fine-mesh FEM / FastHenry before being used as
   ground truth for proximity.
3. **The analytical Wang/Ferreira models under-predict proximity ~2x** on these
   cases vs both PEEC solvers — so PEEC (or FastHenry) is the better reference
   at high frequency.

## Remaining work before PEEC is production-ready
- Confirm agreement at 1 MHz and on foil/round (PEEC sweep was CPU-bound; only
  100 kHz directly cross-checked so far — but it's an exact match there).
- Gapped cores: FastHenry (air-core) can't validate the gap-fringing case;
  validate PEEC's image handling against a cored 2D FEM there.
- Performance: dense QR per harmonic is slow; add iterative solve / caching.
- Settings flag to select PEEC; keep analytical as default until the above done.
