#!/bin/bash
# Demo script for generate_temperature_tests.py

echo "=========================================================================="
echo "Temperature Test Generator - Demo"
echo "=========================================================================="
echo ""

# Set paths
SCRIPT_DIR="/home/alf/OpenMagnetics/MKF/src/scripts"
EXAMPLES_DIR="C:\Users\Alfonso\wuerth\Ansyas\examples"
TEST_FILE="/home/alf/OpenMagnetics/MKF/tests/TestTemperature.cpp"

echo "Script location: $SCRIPT_DIR/generate_temperature_tests.py"
echo "Examples folder: $EXAMPLES_DIR"
echo "Test file: $TEST_FILE"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found"
    exit 1
fi

echo "Python version:"
python3 --version
echo ""

# Show help
echo "=========================================================================="
echo "Available Options:"
echo "=========================================================================="
python3 "$SCRIPT_DIR/generate_temperature_tests.py" --help
echo ""

# Run in dry-run mode first
echo "=========================================================================="
echo "Running in DRY-RUN mode (no changes will be made):"
echo "=========================================================================="
python3 "$SCRIPT_DIR/generate_temperature_tests.py" "$EXAMPLES_DIR" --dry-run
echo ""

# Show actual command to run
echo "=========================================================================="
echo "To actually generate tests, run:"
echo "=========================================================================="
echo "python3 $SCRIPT_DIR/generate_temperature_tests.py $EXAMPLES_DIR"
echo ""
echo "Or with placeholders (recommended until Icepak is fully implemented):"
echo "python3 $SCRIPT_DIR/generate_temperature_tests.py $EXAMPLES_DIR"
echo ""
echo "To skip files when Icepak is unavailable:"
echo "python3 $SCRIPT_DIR/generate_temperature_tests.py $EXAMPLES_DIR --no-placeholders"
echo ""

# Check current test file
echo "=========================================================================="
echo "Current TestTemperature.cpp status:"
echo "=========================================================================="
if [ -f "$TEST_FILE" ]; then
    TEST_COUNT=$(grep -c "TEST_CASE" "$TEST_FILE" 2>/dev/null || echo "0")
    echo "File exists: $TEST_FILE"
    echo "Number of existing tests: $TEST_COUNT"
    echo ""
    echo "Recent tests:"
    grep "TEST_CASE" "$TEST_FILE" | tail -5
else
    echo "Test file not found - will be created on first run"
fi

echo ""
echo "=========================================================================="
echo "Demo complete!"
echo "=========================================================================="
