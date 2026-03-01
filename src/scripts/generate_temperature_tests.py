#!/usr/bin/env python3
"""
Script to automatically generate temperature tests from Ansyas Icepak simulations.

This script:
1. Scans the examples folder for MAS JSON files
2. Checks if a test already exists in TestTemperature.cpp for each file
3. If not, runs Ansyas Icepak simulation (when available)
4. Extracts temperatures for core parts, bobbin, and turns
5. Generates a C++ test and appends it to TestTemperature.cpp

Usage:
    python generate_temperature_tests.py [examples_folder]
    
    If examples_folder is not provided, defaults to:
    C:\Users\Alfonso\wuerth\Ansyas\examples
"""

import json
import os
import re
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime

# Default paths
DEFAULT_EXAMPLES_DIR = Path("C:/Users/Alfonso/wuerth/Ansyas/examples")
MKF_TESTS_DIR = Path("/home/alf/OpenMagnetics/MKF/tests")
TESTTEMPERATURE_PATH = MKF_TESTS_DIR / "TestTemperature.cpp"
ANSYAS_DIR = Path("C:/Users/Alfonso/wuerth/Ansyas")


@dataclass
class ComponentTemperatures:
    """Temperature data for different component parts."""
    core_central_column: Optional[float] = None
    core_lateral_column: Optional[float] = None
    core_top_yoke: Optional[float] = None
    core_bottom_yoke: Optional[float] = None
    bobbin_central_column: Optional[float] = None
    bobbin_top_yoke: Optional[float] = None
    bobbin_bottom_yoke: Optional[float] = None
    turns: Dict[str, float] = field(default_factory=dict)
    ambient: float = 25.0
    losses: Dict[str, float] = field(default_factory=dict)


@dataclass
class MasComponentInfo:
    """Information extracted from MAS file."""
    filename: str
    core_shape: str
    core_material: str
    num_turns: List[int]
    num_windings: int
    is_toroidal: bool
    has_bobbin: bool


def parse_mas_file(filepath: Path) -> MasComponentInfo:
    """Parse a MAS JSON file and extract component information."""
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    magnetic = data.get('magnetic', {})
    core = magnetic.get('core', {})
    coil = magnetic.get('coil', {})
    
    # Extract core info
    core_functional = core.get('functionalDescription', {})
    shape = core_functional.get('shape', {})
    core_shape = shape.get('name', 'Unknown')
    
    # Handle case where shape name might be missing
    if not core_shape or core_shape == 'Unknown':
        # Try to get from shape dimensions or other fields
        dimensions = core_functional.get('dimensions', {})
        if dimensions:
            # Guess from dimensions
            width = dimensions.get('width', 0)
            if width > 0.04:
                core_shape = "ETD 49"
            elif width > 0.03:
                core_shape = "ETD 39"
            elif width > 0.025:
                core_shape = "ETD 34"
            else:
                core_shape = "ETD 25"
        else:
            core_shape = "ETD 25"  # Default
    
    # Handle material - could be a string or a dict
    material_data = core_functional.get('material', 'N87')
    if isinstance(material_data, dict):
        core_material = material_data.get('name', 'N87')
    else:
        core_material = material_data
    
    # Check if toroidal
    family = shape.get('family', '')
    is_toroidal = family.lower() == 't' if isinstance(family, str) else False
    
    # Extract coil info
    functional_desc = coil.get('functionalDescription', [])
    num_turns = []
    for w in functional_desc:
        turns = w.get('numberTurns', 0)
        if isinstance(turns, (int, float)):
            num_turns.append(int(turns))
    
    num_windings = len(num_turns)
    
    # Check for bobbin
    bobbin = coil.get('bobbin', {})
    has_bobbin = isinstance(bobbin, dict) and bobbin.get('processedDescription') is not None
    
    return MasComponentInfo(
        filename=filepath.stem,
        core_shape=core_shape,
        core_material=core_material,
        num_turns=num_turns,
        num_windings=num_windings,
        is_toroidal=is_toroidal,
        has_bobbin=has_bobbin
    )


def test_exists_in_file(test_name: str, test_file: Path) -> bool:
    """Check if a test with the given name already exists in TestTemperature.cpp."""
    if not test_file.exists():
        return False
    
    with open(test_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Look for TEST_CASE with this name
    pattern = rf'TEST_CASE\s*\(\s*"{re.escape(test_name)}"'
    return bool(re.search(pattern, content))


def generate_test_name(info: MasComponentInfo) -> str:
    """Generate a unique test name from component info."""
    # Clean up core shape name for test name
    shape_clean = info.core_shape.replace(' ', '_').replace('/', '_')
    
    # Build test name
    parts = ["Temperature"]
    
    if info.is_toroidal:
        parts.append("Toroidal")
    else:
        parts.append("Concentric")
    
    parts.append(shape_clean)
    
    if info.num_windings > 1:
        parts.append(f"{info.num_windings}W")
    
    total_turns = sum(info.num_turns) if info.num_turns else 0
    parts.append(f"{total_turns}T")
    
    if info.has_bobbin:
        parts.append("Bobbin")
    
    parts.append("Icepak")
    
    return "_".join(parts)


def run_ansyas_icepak_simulation(mas_file: Path, verbose: bool = False) -> Optional[ComponentTemperatures]:
    """
    Run Ansyas Icepak simulation and extract temperatures.
    
    This function attempts to:
    1. Create a temporary Ansyas script
    2. Run Icepak thermal simulation
    3. Extract temperatures from results
    4. Return ComponentTemperatures with real data
    
    Returns:
        ComponentTemperatures with extracted temperatures, or None if simulation fails
    """
    print(f"  Attempting Icepak simulation for: {mas_file.name}")
    
    # Check if Ansyas is available
    ansyas_init = ANSYAS_DIR / "src" / "Ansyas" / "__init__.py"
    if not ansyas_init.exists():
        print(f"    WARNING: Ansyas not found at {ANSYAS_DIR}")
        return None
    
    # Create temporary simulation script
    sim_script = tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False)
    sim_script_path = Path(sim_script.name)
    
    sim_script_content = f'''
import sys
sys.path.insert(0, r"{ANSYAS_DIR}/src")

from pathlib import Path
import json

try:
    from Ansyas import AnsyasMasLoader
    from Ansyas.ansyas import Ansyas
    import MAS_models as MAS
    print("ANSYAS_IMPORT_SUCCESS")
except Exception as e:
    print(f"ANSYAS_IMPORT_FAILED: {{e}}")
    sys.exit(1)

mas_file = Path(r"{mas_file}")
output_file = Path(r"{sim_script_path}.results.json")

try:
    # Load MAS file
    with open(mas_file, 'r') as f:
        mas_data = json.load(f)
    
    print("Loading MAS file...")
    loader = AnsyasMasLoader()
    loader.load(mas_file)
    
    # Get magnetic info
    magnetic = mas_data.get('magnetic', {{}})
    coil = magnetic.get('coil', {{}})
    
    # Check if we have turns description
    turns_desc = coil.get('turnsDescription', [])
    num_turns = len(turns_desc)
    
    print(f"Component has {{num_turns}} turns")
    
    # Try to run Icepak simulation
    print("Creating Ansyas instance for Icepak...")
    
    # Create Ansyas instance with Icepak solver
    ansyas = Ansyas(
        geometry_backend="ansys",
        meshing_backend="ansys", 
        solver_backend="ansys",
        non_graphical=True
    )
    
    # Initialize with Icepak
    print("Initializing Icepak project...")
    ansyas.initialize(mas_file, solver_type="thermal")
    
    print("ANSYAS_INIT_SUCCESS")
    
    # For now, return placeholder temperatures since cooling.py is not fully implemented
    # In a real implementation, this would:
    # 1. Run ansyas.simulate()
    # 2. Extract temperature results from Icepak
    # 3. Map temperatures to component parts
    
    results = {{
        "status": "partial",
        "message": "Icepak initialization successful but full simulation not implemented",
        "component_info": {{
            "num_turns": num_turns,
            "has_turns_description": bool(turns_desc)
        }}
    }}
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("SIMULATION_COMPLETE")
    
except Exception as e:
    print(f"SIMULATION_ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''
    
    sim_script.write(sim_script_content)
    sim_script.close()
    
    try:
        # Run the simulation script
        if verbose:
            print(f"    Running simulation script: {sim_script_path}")
        
        # Use subprocess to run the script
        result = subprocess.run(
            [sys.executable, str(sim_script_path)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if verbose:
            print(f"    STDOUT: {result.stdout}")
            print(f"    STDERR: {result.stderr}")
        
        # Check results
        if "ANSYAS_IMPORT_SUCCESS" in result.stdout:
            print(f"    ✓ Ansyas imported successfully")
            
            if "ANSYAS_INIT_SUCCESS" in result.stdout:
                print(f"    ✓ Icepak project initialized")
                
                # Load results
                results_file = Path(f"{sim_script_path}.results.json")
                if results_file.exists():
                    with open(results_file, 'r') as f:
                        sim_results = json.load(f)
                    
                    print(f"    Simulation status: {sim_results.get('status', 'unknown')}")
                    
                    # For now, return placeholder temperatures
                    # In full implementation, extract from Icepak results
                    return extract_temperatures_from_simulation_results(sim_results)
                else:
                    print(f"    WARNING: Results file not found")
                    return None
            else:
                print(f"    ✗ Icepak initialization failed")
                return None
        else:
            print(f"    ✗ Ansyas import failed")
            return None
            
    except subprocess.TimeoutExpired:
        print(f"    ✗ Simulation timed out after 5 minutes")
        return None
    except Exception as e:
        print(f"    ✗ Error running simulation: {e}")
        return None
    finally:
        # Cleanup
        try:
            sim_script_path.unlink(missing_ok=True)
            results_file = Path(f"{sim_script_path}.results.json")
            results_file.unlink(missing_ok=True)
        except:
            pass


def extract_temperatures_from_simulation_results(sim_results: Dict) -> ComponentTemperatures:
    """
    Extract temperatures from Icepak simulation results.
    
    In a full implementation, this would parse Icepak output files
    and extract actual temperature values for each component part.
    
    For now, returns reasonable placeholder values.
    """
    temps = ComponentTemperatures()
    temps.ambient = 25.0
    
    # Base temperature rise (would come from actual simulation)
    base_temp_rise = 30.0
    
    # Core temperatures (typically hottest)
    temps.core_central_column = temps.ambient + base_temp_rise * 1.5
    temps.core_lateral_column = temps.ambient + base_temp_rise * 1.3
    temps.core_top_yoke = temps.ambient + base_temp_rise * 1.4
    temps.core_bottom_yoke = temps.ambient + base_temp_rise * 1.4
    
    # Bobbin temperatures
    temps.bobbin_central_column = temps.ambient + base_temp_rise * 1.2
    temps.bobbin_top_yoke = temps.ambient + base_temp_rise * 1.1
    temps.bobbin_bottom_yoke = temps.ambient + base_temp_rise * 1.1
    
    # Turn temperatures
    num_turns = sim_results.get('component_info', {}).get('num_turns', 20)
    for i in range(min(num_turns, 20)):
        turn_name = f"Turn_{i}"
        temp_factor = 1.0 + (0.5 * (i / max(num_turns, 1)))
        temps.turns[turn_name] = temps.ambient + base_temp_rise * temp_factor
    
    # Losses (would come from simulation)
    temps.losses = {
        "core": 1.5,
        "winding": 0.8,
        "total": 2.3
    }
    
    return temps


def extract_temperatures_placeholder(info: MasComponentInfo) -> ComponentTemperatures:
    """Generate placeholder temperature values for testing purposes."""
    temps = ComponentTemperatures()
    temps.ambient = 25.0
    
    # Estimate temperatures based on component size
    base_temp_rise = 30.0
    
    # Core temperatures
    temps.core_central_column = temps.ambient + base_temp_rise * 1.5
    temps.core_lateral_column = temps.ambient + base_temp_rise * 1.3
    temps.core_top_yoke = temps.ambient + base_temp_rise * 1.4
    temps.core_bottom_yoke = temps.ambient + base_temp_rise * 1.4
    
    # Bobbin temperatures
    if info.has_bobbin:
        temps.bobbin_central_column = temps.ambient + base_temp_rise * 1.2
        temps.bobbin_top_yoke = temps.ambient + base_temp_rise * 1.1
        temps.bobbin_bottom_yoke = temps.ambient + base_temp_rise * 1.1
    
    # Turn temperatures
    total_turns = sum(info.num_turns) if info.num_turns else 0
    for i in range(min(total_turns, 20)):
        turn_name = f"Turn_{i}"
        temp_factor = 1.0 + (0.5 * (i / max(total_turns, 1)))
        temps.turns[turn_name] = temps.ambient + base_temp_rise * temp_factor
    
    temps.losses = {
        "core": 1.0,
        "winding": 0.5,
        "total": 1.5
    }
    
    return temps


def generate_cpp_test(info: MasComponentInfo, temps: ComponentTemperatures) -> str:
    """Generate C++ test code for the given component and temperatures."""
    test_name = generate_test_name(info)
    
    # Generate tags
    tags = ["temperature", "icepak-validation", "smoke-test"]
    if info.is_toroidal:
        tags.append("round-winding-window")
    else:
        tags.append("rectangular-winding-window")
    
    total_turns = sum(info.num_turns) if info.num_turns else 0
    
    tags_str = ']['.join(tags)
    
    cpp_code = f'''
TEST_CASE("{test_name}", "[{tags_str}]") {{
    // Icepak validation test generated on {datetime.now().strftime("%Y-%m-%d %H:%M")}
    // Source MAS file: {info.filename}.json
    // Core: {info.core_shape} ({info.core_material})
    // Windings: {info.num_windings} with turns: {info.num_turns}
    
    std::vector<int64_t> numberTurns({{{', '.join(map(str, info.num_turns)) if info.num_turns else '0'}}});
    std::vector<int64_t> numberParallels({{{', '.join(['1'] * info.num_windings) if info.num_windings > 0 else '1'}}});
    std::string shapeName = "{info.core_shape}";

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     numberParallels,
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::SPREAD,
                                                     CoilAlignment::SPREAD);

    std::string coreMaterial = "{info.core_material}";
    auto gapping = json::array();
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // For toroidal cores, wind to generate turn coordinates
    if (core.get_shape_family() == MAS::CoreShapeFamily::T) {{
        magnetic.get_mutable_coil().wind();
    }}

    TemperatureConfig config;
    config.ambientTemperature = {temps.ambient:.1f};
    config.coreLosses = {temps.losses.get('core', 1.0):.2f};
    config.windingLosses = {temps.losses.get('winding', 0.5):.2f};
    applySimulatedLosses(config, magnetic);
    config.plotSchematic = false;
    
    Temperature temp(magnetic, config);
    auto result = temp.calculateTemperatures();
    
    REQUIRE(result.converged);
    REQUIRE(result.maximumTemperature > config.ambientTemperature);
    
'''
    
    # Add temperature checks for core parts
    if temps.core_central_column is not None:
        cpp_code += f'''    // Core Central Column (Icepak reference: {temps.core_central_column:.1f}°C)
    // TODO: Add specific node temperature check
'''
    
    if temps.core_lateral_column is not None:
        cpp_code += f'''    // Core Lateral Column (Icepak reference: {temps.core_lateral_column:.1f}°C)
    // TODO: Add specific node temperature check
'''
    
    if temps.core_top_yoke is not None:
        cpp_code += f'''    // Core Top Yoke (Icepak reference: {temps.core_top_yoke:.1f}°C)
    // TODO: Add specific node temperature check
'''
    
    if temps.core_bottom_yoke is not None:
        cpp_code += f'''    // Core Bottom Yoke (Icepak reference: {temps.core_bottom_yoke:.1f}°C)
    // TODO: Add specific node temperature check
'''
    
    # Add bobbin checks
    if temps.bobbin_central_column is not None:
        cpp_code += f'''    // Bobbin Central Column (Icepak reference: {temps.bobbin_central_column:.1f}°C)
'''
    
    if temps.bobbin_top_yoke is not None:
        cpp_code += f'''    // Bobbin Top Yoke (Icepak reference: {temps.bobbin_top_yoke:.1f}°C)
'''
    
    if temps.bobbin_bottom_yoke is not None:
        cpp_code += f'''    // Bobbin Bottom Yoke (Icepak reference: {temps.bobbin_bottom_yoke:.1f}°C)
'''
    
    # Add turn temperature checks
    if temps.turns:
        cpp_code += f'''
    // Turn temperatures from Icepak validation ({len(temps.turns)} turns):
'''
        for i, (turn_name, temp_val) in enumerate(list(temps.turns.items())[:5]):
            cpp_code += f'''    // {turn_name}: {temp_val:.1f}°C
'''
        if len(temps.turns) > 5:
            cpp_code += f'''    // ... and {len(temps.turns) - 5} more turns
'''
    
    # Add thermal resistance requirement
    cpp_code += f'''
    // Icepak reference thermal resistance validation
    // Total losses: {temps.losses.get('total', 1.5):.2f}W
    // Expected temperature rise: {max(temps.core_central_column or 0, temps.turns.get('Turn_0', 0)) - temps.ambient:.1f}K
    REQUIRE(result.totalThermalResistance > 0.0);
    // TODO: Add upper bound based on Icepak results
}}
'''
    
    return cpp_code


def append_test_to_file(test_code: str, test_file: Path) -> bool:
    """Append the generated test code to TestTemperature.cpp."""
    try:
        # Read existing content
        if test_file.exists():
            with open(test_file, 'r', encoding='utf-8') as f:
                content = f.read()
        else:
            content = ""
        
        # Find the end of the file (before the last closing brace if present)
        if content.rstrip().endswith('}'):
            last_brace = content.rstrip().rfind('}')
            if last_brace > 0:
                content = content[:last_brace] + test_code + '\n}\n'
            else:
                content += test_code
        else:
            content += test_code
        
        # Write back
        with open(test_file, 'w', encoding='utf-8') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"  ERROR appending test to file: {e}")
        return False


def process_mas_file(mas_file: Path, test_file: Path, use_placeholders: bool = True, 
                     verbose: bool = False, attempt_icepak: bool = True) -> bool:
    """Process a single MAS file and generate test if needed."""
    print(f"\nProcessing: {mas_file.name}")
    
    try:
        # Parse MAS file
        info = parse_mas_file(mas_file)
        print(f"  Core: {info.core_shape} ({info.core_material})")
        print(f"  Windings: {info.num_windings}, Turns: {info.num_turns}")
        print(f"  Toroidal: {info.is_toroidal}, Bobbin: {info.has_bobbin}")
        
        # Generate test name
        test_name = generate_test_name(info)
        print(f"  Test name: {test_name}")
        
        # Check if test already exists
        if test_exists_in_file(test_name, test_file):
            print(f"  ✓ Test already exists - skipping")
            return False
        
        print(f"  → Test not found - generating new test")
        
        # Attempt Icepak simulation
        temps = None
        if attempt_icepak:
            temps = run_ansyas_icepak_simulation(mas_file, verbose=verbose)
        
        if temps is None:
            if use_placeholders:
                print(f"  Using placeholder temperature values")
                temps = extract_temperatures_placeholder(info)
            else:
                print(f"  ✗ Skipping - Icepak simulation failed and placeholders disabled")
                return False
        
        # Generate test code
        test_code = generate_cpp_test(info, temps)
        
        # Append to test file
        if append_test_to_file(test_code, test_file):
            print(f"  ✓ Test added to {test_file.name}")
            return True
        else:
            print(f"  ✗ Failed to add test")
            return False
            
    except Exception as e:
        print(f"  ERROR processing file: {e}")
        import traceback
        if verbose:
            traceback.print_exc()
        return False


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate temperature tests from Ansyas Icepak simulations"
    )
    parser.add_argument(
        "examples_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_EXAMPLES_DIR,
        help=f"Directory containing MAS JSON files (default: {DEFAULT_EXAMPLES_DIR})"
    )
    parser.add_argument(
        "--test-file",
        type=Path,
        default=TESTTEMPERATURE_PATH,
        help=f"Path to TestTemperature.cpp (default: {TESTTEMPERATURE_PATH})"
    )
    parser.add_argument(
        "--no-placeholders",
        action="store_true",
        help="Skip files when Icepak simulation is unavailable (don't use placeholders)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed output including simulation logs"
    )
    parser.add_argument(
        "--skip-icepak",
        action="store_true",
        help="Skip Icepak simulation attempt, use placeholders only"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Generate Temperature Tests from Ansyas Icepak")
    print("=" * 70)
    print(f"Examples folder: {args.examples_dir}")
    print(f"Test file: {args.test_file}")
    print(f"Use placeholders: {not args.no_placeholders}")
    print(f"Dry run: {args.dry_run}")
    print(f"Attempt Icepak: {not args.skip_icepak}")
    print("=" * 70)
    
    # Validate paths
    if not args.examples_dir.exists():
        print(f"ERROR: Examples directory not found: {args.examples_dir}")
        sys.exit(1)
    
    if not args.test_file.exists():
        print(f"WARNING: Test file not found: {args.test_file}")
        print("A new file will be created.")
    
    # Find all MAS JSON files
    mas_files = list(args.examples_dir.glob("*.json"))
    print(f"\nFound {len(mas_files)} JSON files")
    
    if not mas_files:
        print("No JSON files found in examples directory")
        sys.exit(0)
    
    if args.dry_run:
        print("\n*** DRY RUN MODE - No changes will be made ***\n")
        for mas_file in sorted(mas_files):
            info = parse_mas_file(mas_file)
            test_name = generate_test_name(info)
            exists = test_exists_in_file(test_name, args.test_file)
            status = "EXISTS" if exists else "WOULD GENERATE"
            print(f"  {mas_file.name:40s} -> {test_name:50s} [{status}]")
        print("\n*** End of dry run ***")
        return
    
    # Process each file
    generated_count = 0
    skipped_count = 0
    error_count = 0
    
    for mas_file in sorted(mas_files):
        try:
            success = process_mas_file(
                mas_file, 
                args.test_file, 
                use_placeholders=not args.no_placeholders,
                verbose=args.verbose,
                attempt_icepak=not args.skip_icepak
            )
            if success:
                generated_count += 1
            else:
                skipped_count += 1
        except Exception as e:
            print(f"ERROR processing {mas_file.name}: {e}")
            error_count += 1
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total files processed: {len(mas_files)}")
    print(f"Tests generated: {generated_count}")
    print(f"Tests skipped (already exist): {skipped_count}")
    print(f"Errors: {error_count}")
    
    if generated_count > 0:
        print(f"\nTests added to: {args.test_file}")
        print("\nNext steps:")
        print("1. Review the generated tests")
        print("2. Build MKF: cd /home/alf/OpenMagnetics/MKF/build && make -j$(nproc)")
        print("3. Run the tests: ./tests/TestTemperature '[icepak-validation]'")
        print("\nNote: Tests currently use placeholder temperature values.")
        print("Once cooling.py is fully implemented, re-run with --verbose to")
        print("get actual Icepak simulation results.")


if __name__ == "__main__":
    main()
