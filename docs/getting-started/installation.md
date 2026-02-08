# Installation

This guide covers building MKF from source on Linux and Windows.

## Prerequisites

### Linux (Ubuntu/Debian)

```bash
# Install build tools
sudo apt-get update
sudo apt-get install -y cmake ninja-build

# Install GCC 14 (C++23 required)
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install -y gcc-14 g++-14
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 10

# Install dependencies
sudo apt install -y gnuplot libopenblas-dev

# Install Node.js and quicktype (for schema generation)
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
npm install -g quicktype
```

### Windows

1. Install [Visual Studio 2022](https://visualstudio.microsoft.com/) with C++ workload
2. Install [CMake](https://cmake.org/download/) (3.18+)
3. Install [Ninja](https://ninja-build.org/) via `choco install ninja`
4. Install [Node.js](https://nodejs.org/) and run `npm install -g quicktype`

## Building from Source

### Clone the Repository

```bash
git clone --recursive https://github.com/OpenMagnetics/MKF.git
cd MKF
```

!!! note "Submodules"
    The `--recursive` flag is important to clone all submodules (MAS schemas, third-party libraries).

### Configure and Build

=== "Linux"

    ```bash
    mkdir build && cd build
    cmake .. -G "Ninja"
    ninja
    ```

=== "Windows"

    ```powershell
    mkdir build
    cd build
    cmake .. -G "Ninja"
    ninja
    ```

### Run Tests

```bash
# Run all tests
./MKF_tests

# Run smoke tests only (fast)
./MKF_tests "[smoke-test]"

# Run tests by category
./MKF_tests "[physical-model]"
./MKF_tests "[adviser]"

# Run a specific test
./MKF_tests "Test_Skin_Depth_Wire_Material_Data_20C"
```

## Python Package

MKF can be installed as a Python package (PyMKF):

```bash
# From the repository root
python3 -m pip install -e . -vvv
```

Then use it in Python:

```python
import PyMKF

# Access the magnetics functionality
core = PyMKF.find_core_by_name("E 42/21/15")
```

## Troubleshooting

### QuickType/MAS Generation Fails

If the schema generation fails during build, generate manually:

```bash
cd build
quicktype -l c++ -s schema ../MAS/schemas/MAS.json \
  -S ../MAS/schemas/magnetic.json \
  -S ../MAS/schemas/magnetic/core.json \
  -S ../MAS/schemas/magnetic/coil.json \
  # ... (see BUILD.md for full command)
  -o MAS/MAS.hpp --namespace MAS \
  --source-style single-source \
  --type-style pascal-case \
  --member-style underscore-case \
  --enumerator-style upper-underscore-case \
  --no-boost
```

### GCC Version Issues

MKF requires C++23 features. Ensure you have GCC 14+ or a compatible compiler:

```bash
g++ --version  # Should show 14.x or higher
```

### Missing Dependencies

If CMake cannot find dependencies, check that all prerequisites are installed:

```bash
# Check for OpenBLAS
ldconfig -p | grep openblas

# Check for gnuplot
which gnuplot
```

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `INCLUDE_MKF_TESTS` | ON | Build test suite |
| `EMBED_MAS_DATA` | ON | Embed database files in binary |
| `CMAKE_BUILD_TYPE` | Release | Build type (Release/Debug) |

Example with options:

```bash
cmake .. -G "Ninja" \
  -DCMAKE_BUILD_TYPE=Release \
  -DINCLUDE_MKF_TESTS=ON \
  -DEMBED_MAS_DATA=ON
```
