# MKF Session Handoff - April 24, 2026

## Status

**Build:** Clean build with latest MAS (185f23e) + schema fixes. Compiles successfully.

**Tests:** Running smoke tests (`[smoke-test]`) still crashes. Need to run with exclusions.

## Completed

### 1. MAS Schema Fixes (in MAS submodule)
- `schemas/magnetic.json`: Added `title: "MagneticDatasheetApplication"` to `$defs.magneticDatasheetApplication`
- `schemas/utils.json`: Added new `$def: "dimension"` with `oneOf` [number, dimensionWithTolerance]
- `schemas/magnetic/bobbin.json`: Replaced inline `oneOf` with `$ref` to `dimension`
- `schemas/magnetic/core/shape.json`: Replaced inline `oneOf` with `$ref` to `dimension`
- `schemas/magnetic.json`: Renamed conflicting `impedanceAtFrequency` title to `DatasheetImpedancePoint`
- Result: MAS.hpp generates correct type names (`Application`, `Dimension`, `ImpedanceAtFrequency`)

### 2. Bug Fixes Applied

#### Bug A: Spline crashes with < 3 points
**Files changed:**
- `src/physical_models/InitialPermeability.h`: Changed `initialPermeabilityTemperatureInterps` from `variant<double, tk::spline>` to `map<string, function<double(double)>>`
- `src/physical_models/InitialPermeability.cpp:620-662`: Sorts points by temperature, handles 2-point case with manual linear interpolation, 1-point as constant
- `src/physical_models/CoreLosses.h`: Changed `lossFactorInterps` from `variant<double, tk::spline>` to `map<string, function<double(double)>>`
- `src/physical_models/CoreLosses.cpp:2375-2418`: Sorts points by frequency, handles 2-point case with manual linear interpolation

**Status:** Fixed. Previously failing tests now pass:
- `Test_CoreAdviserAvailableCores_Only_Toroids_High_Power` ✅
- `Core_Processed_Description_Web_2` ✅
- `Test_CoreMaterialCrossReferencer_All_Core_Materials` ✅

#### Bug B: Coil adviser SIGSEGV - stale turn index cache
**File changed:**
- `src/constructive_models/Coil.cpp:1077-1102`: `get_turn_index_by_name()` now validates cached indices against current turns vector, clears entire cache if stale

**Status:** Partially fixed. Some random tests pass (Random_2, Random_0, Random_1, etc.) but Random_15 and Random_16 still crash with SIGSEGV.

**Next steps needed:** The remaining coil adviser crashes are a deeper memory corruption issue. The GDB trace shows crash in `get_turn_index_by_name` when comparing turn names, suggesting either:
- Dangling reference in `get_turns_description()` (returns optional by value, `.value()` may dangle)
- Vector reallocation during turn placement in `wind_toroidal_additional_turns()`
- Stack corruption from the collision detection loop

## Still Broken

### Coil Adviser Random Tests (11 of 18 still fail)
- `Test_CoilAdviser_Random_15` - SIGSEGV 
- `Test_CoilAdviser_Random_16` - SIGSEGV
- Random_2, 3, 5, 7, 8, 9, 10, 12, 14 also fail (all tagged `[bug]`)

**Crash location:** `src/constructive_models/Coil.cpp:1096` in `get_turn_index_by_name()`
**Stack trace:** `wind_toroidal_additional_turns()` → `get_turn_index_by_name()` → string comparison crash

### Other Known Issues
- Multiple spline assertion failures in various tests when run in certain orders (all pre-existing, tagged `[bug]`)

## Running Tests

Current working command (excludes known failing tests):
```bash
./build_coverage/MKF_tests "[smoke-test]~[bug]" --abort
```

To test individual fixes:
```bash
./build_coverage/MKF_tests "Test_CoreAdviserAvailableCores_Only_Toroids_High_Power"
./build_coverage/MKF_tests "Core_Processed_Description_Web_2"
./build_coverage/MKF_tests "Test_CoilAdviser_Random_2"
```

## Files to Review

1. `src/constructive_models/Coil.cpp:4886-5300` - `wind_toroidal_additional_turns()` - main crash site
2. `src/constructive_models/Coil.cpp:1077-1102` - `get_turn_index_by_name()` - attempted fix
3. `src/physical_models/InitialPermeability.cpp` - spline fix
4. `src/physical_models/CoreLosses.cpp` - spline fix

## MAS Submodule

Currently on `main` branch with uncommitted schema fixes. To commit:
```bash
cd MAS
git add schemas/
git commit -m "fix: resolve quicktype naming conflicts for C++ generation"
```

## Coverage

Coverage build directory: `build_coverage/`
- Configure: `cmake -B build_coverage -G Ninja -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON -DCMAKE_CXX_COMPILER=g++-14 -DINCLUDE_MKF_TESTS=ON -DEMBED_MAS_DATA=ON`
- Build: `ninja -C build_coverage -j4`
- Run: `cd build_coverage && ./MKF_tests "[smoke-test]~[bug]" --abort`
- Collect: `lcov --gcov-tool /usr/bin/gcov-14 --capture --directory . --output-file coverage.info --ignore-errors gcov,source`

## Disk Space

~38GB freed from cleanup. Current usage ~82%.

## Next Priority

1. **Fix remaining coil adviser crashes** - Deep investigation of memory corruption in `wind_toroidal_additional_turns()`
2. **Run full smoke test suite** once crashes are fixed
3. **Generate coverage report** 
4. **Commit MAS schema fixes** to submodule