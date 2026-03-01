# Temperature Test Generator

This script automatically generates C++ test cases for the MKF Temperature model by running Ansyas Icepak simulations on MAS files.

## Overview

The script performs the following workflow:
1. Scans the examples folder for MAS JSON files
2. Checks if a corresponding test already exists in `TestTemperature.cpp`
3. If not, attempts to run an Ansyas Icepak thermal simulation
4. Extracts temperatures for:
   - Core columns (central, lateral)
   - Core yokes (top, bottom)
   - Bobbin column and yokes
   - Each individual turn
5. Generates a C++ test case with the extracted values
6. Appends the test to `TestTemperature.cpp`

## Prerequisites

- Python 3.7+
- Ansys Icepak properly configured (for actual simulations)
- Access to MKF repository in WSL
- MAS files in the examples folder

## Usage

### Basic Usage

```bash
# Run with default examples folder
python3 /home/alf/OpenMagnetics/MKF/src/scripts/generate_temperature_tests.py

# Specify custom examples folder
python3 /home/alf/OpenMagnetics/MKF/src/scripts/generate_temperature_tests.py /path/to/examples

# Use with Windows path
python3 /home/alf/OpenMagnetics/MKF/src/scripts/generate_temperature_tests.py "C:/Users/Alfonso/wuerth/Ansyas/examples"
```

### Options

```bash
# Don't generate tests if Icepak simulation fails
python3 generate_temperature_tests.py --no-placeholders

# Preview what would be generated without making changes
python3 generate_temperature_tests.py --dry-run

# Specify custom test file location
python3 generate_temperature_tests.py --test-file /path/to/TestTemperature.cpp
```

## Current Status

**Important**: The script currently uses placeholder temperature values because:

1. The `cooling.py` module in Ansyas has a stub implementation (`pass` statement)
2. Icepak integration is not fully implemented yet
3. Temperature extraction from simulation results needs implementation

### To Enable Real Icepak Simulations:

1. Complete the implementation in `Ansyas/src/Ansyas/cooling.py`
2. Implement the `run_ansyas_icepak_simulation()` function in this script
3. Update the `extract_temperatures_placeholder()` function to use real data

## Generated Test Structure

Each generated test follows this pattern:

```cpp
TEST_CASE("Temperature_Concentric_ETD_49_2W_35T_Bobbin_Icepak", "[temperature][icepak-validation][smoke-test]") {
    // Icepak validation test generated on 2025-01-15 10:30
    // Source MAS file: concentric_transformer.json
    // Core: ETD 49 (3C97)
    // Windings: 2 with turns: [20, 15]
    
    std::vector<int64_t> numberTurns({20, 15});
    // ... test implementation ...
    
    // Icepak reference thermal resistance validation
    REQUIRE(result.totalThermalResistance > 0.0);
}
```

## Workflow Example

```bash
# 1. Navigate to the script directory
cd /home/alf/OpenMagnetics/MKF/src/scripts

# 2. Run the generator
python3 generate_temperature_tests.py

# 3. Review generated tests
cat /home/alf/OpenMagnetics/MKF/tests/TestTemperature.cpp | tail -100

# 4. Build and run the new tests
cd /home/alf/OpenMagnetics/MKF/build
make -j$(nproc)
./tests/TestTemperature "[icepak-validation]"
```

## Output

The script provides detailed output:

```
======================================================================
Generate Temperature Tests from Ansyas Icepak
======================================================================
Examples folder: C:\Users\Alfonso\wuerth\Ansyas\examples
Test file: /home/alf/OpenMagnetics/MKF/tests/TestTemperature.cpp
======================================================================

Found 5 JSON files

Processing: concentric_transformer.json
  Core: ETD 49 (3C97)
  Windings: 2, Turns: [20, 15]
  Toroidal: False, Bobbin: True
  Test name: Temperature_Concentric_ETD_49_2W_35T_Bobbin_Icepak
  → Test not found - generating new test
  Using placeholder temperature values
  ✓ Test added to TestTemperature.cpp

======================================================================
SUMMARY
======================================================================
Total files processed: 5
Tests generated: 5
Tests skipped (already exist): 0
Errors: 0

Tests added to: /home/alf/OpenMagnetics/MKF/tests/TestTemperature.cpp
Remember to:
1. Review the generated tests
2. Replace placeholder values with actual Icepak results
3. Run the tests to verify they work
```

## Integration with CI/CD

Add to your CI pipeline:

```yaml
# .github/workflows/temperature-tests.yml
- name: Generate Temperature Tests
  run: |
    python3 /home/alf/OpenMagnetics/MKF/src/scripts/generate_temperature_tests.py \
      --no-placeholders  # Only generate if Icepak is available
```

## Future Enhancements

1. **Real Icepak Integration**: Extract actual temperatures from simulation results
2. **Temperature Mapping**: Map Icepak temperature fields to MKF thermal network nodes
3. **Automated Comparison**: Compare MKF results against Icepak reference values
4. **Batch Processing**: Process multiple files in parallel
5. **Test Categorization**: Auto-generate appropriate tags based on component type

## Troubleshooting

### "Could not import Ansyas"
- Ensure Ansyas is installed and in Python path
- Check that you're running from a Python environment with Ansyas accessible

### "Icepak simulation not yet implemented"
- This is expected until `cooling.py` is fully implemented
- Use `--no-placeholders` to skip these files

### "Test file not found"
- The script will create a new TestTemperature.cpp file
- Ensure you have write permissions to the MKF/tests directory

## Contact

For questions about this script, contact the Ansyas/MKF development team.
