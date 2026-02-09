#!/usr/bin/env python3
import re

# Read the file
with open('/home/alf/OpenMagnetics/MKF2/tests/TestTemperature.cpp', 'r') as f:
    content = f.read()

# Remove lines with std::cout << "..." statements (debug prints)
# Keep only the export helper function std::cout
lines = content.split('\n')
filtered_lines = []
for line in lines:
    # Skip lines that are std::cout debug prints but keep the schematic export message
    if 'std::cout << "' in line or 'std::cerr << "' in line:
        # Keep the schematic export message
        if 'Thermal circuit schematic exported to' in line:
            filtered_lines.append(line)
        continue
    filtered_lines.append(line)

# Join back
new_content = '\n'.join(filtered_lines)

# Write back
with open('/home/alf/OpenMagnetics/MKF2/tests/TestTemperature.cpp', 'w') as f:
    f.write(new_content)

print("Removed debug prints from TestTemperature.cpp")
