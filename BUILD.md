
## Build steps

1. Clone with the schema submodules (MAS, PEAS, CAS, CIAS, TAS, RAS, SAS, cci_coords). The
   build generates its C++ types from these schemas, so a clone without them does not configure:

```
git clone https://github.com/OpenMagnetics/MKF.git
cd MKF
git submodule update --init --recursive
```

2. Install dependencies (Ubuntu; CMake >= 3.18 from Kitware's repository, a signed-by keyring
   instead of the deprecated `apt-key`):

```
sudo install -d -m 0755 /etc/apt/keyrings
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor | sudo tee /etc/apt/keyrings/kitware-archive-keyring.gpg >/dev/null
echo "deb [signed-by=/etc/apt/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" | sudo tee /etc/apt/sources.list.d/kitware.list >/dev/null
sudo apt update
sudo apt install -y cmake ninja-build git

sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install -y gcc-14 g++-14
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 10

# Node.js + quicktype (the MAS/CAS/... C++ types are generated at configure/build time)
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - && sudo apt-get install -y nodejs
npm install -g quicktype

# libngspice: the Kirchhoff converter library (see below) links it in-process
sudo apt install -y libngspice0-dev

sudo apt install -y gnuplot
```

3. Kirchhoff. Natively MKF links the Kirchhoff converter-model library (`MKF_USE_KIRCHHOFF`,
   ON by default, OFF under Emscripten). CMake fetches Kirchhoff and Kelvin from GitHub and builds
   them in an isolated sub-build (`build/kirchhoff-build`), with Kirchhoff's in-process libngspice
   runner enabled, so libngspice must be installed (step 2). To build against local checkouts pass
   `-DKIRCHHOFF_SOURCE_DIR=/path/to/Kirchhoff` and/or `-DKIRCHHOFF_KELVIN_DIR=/path/to/Kelvin`.
   MKF's own `ENABLE_NGSPICE` stays OFF: circuit simulation goes through Kirchhoff.

4. Configure and build (Ninja):

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build MKF_tests
```

   The Kirchhoff sub-build ignores `ninja -j`; bound it with `CMAKE_BUILD_PARALLEL_LEVEL=N` in the
   environment if memory is short.

5. Run the tests. `MKF_tests` is a Catch2 binary; run it directly from the build directory (the
   tests read their fixtures from `../tests/testData`), with Catch2 tag or name filters:

```
cd build
./MKF_tests "[coil]"
```

The Python package is not built from this repository: it lives in
[PyOpenMagnetics](https://github.com/OpenMagnetics/PyOpenMagnetics).

### Fast builds with Ninja (recommended)

If you want the fastest builds (no AddressSanitizer, no gprof), configure a fresh Release build with Ninja:

```bash
# remove any existing build dir to avoid stale flags
rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
# or
ninja -C build -j "$(nproc)"
```

Enable AddressSanitizer for Debug (opt-in):

```bash
rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DENABLE_ASAN=ON
cmake --build build
```

Enable gprof profiling (opt-in):

```bash
rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROFILING=ON
cmake --build build --parallel
```

Quick check that ASan or profiling flags are present after configure:

```bash
grep -R "fsanitize=address\|-pg" build/compile_commands.json || true
```

if quicktype/MAS generation fails in step 4, try to generate it manually:

quicktype -l c++ -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/inputs/topologies/currentTransformer.json -S ../MAS/schemas/inputs/topologies/buck.json -S ../MAS/schemas/inputs/topologies/isolatedBuck.json -S ../MAS/schemas/inputs/topologies/isolatedBuckBoost.json -S ../MAS/schemas/inputs/topologies/boost.json -S ../MAS/schemas/inputs/topologies/pushPull.json -S ../MAS/schemas/inputs/topologies/forward.json -o MAS/MAS.hpp --namespace MAS --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost

quicktype -l python -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/inputs/topologies/currentTransformer.json -S ../MAS/schemas/inputs/topologies/buck.json -S ../MAS/schemas/inputs/topologies/isolatedBuck.json -S ../MAS/schemas/inputs/topologies/isolatedBuckBoost.json -S ../MAS/schemas/inputs/topologies/boost.json -S ../MAS/schemas/inputs/topologies/pushPull.json -S ../MAS/schemas/inputs/topologies/forward.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.py --python-version 3.7 --no-nice-property-names --pydantic-base-model

quicktype -l cs -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/inputs/topologies/currentTransformer.json -S ../MAS/schemas/inputs/topologies/buck.json -S ../MAS/schemas/inputs/topologies/isolatedBuck.json -S ../MAS/schemas/inputs/topologies/isolatedBuckBoost.json -S ../MAS/schemas/inputs/topologies/boost.json -S ../MAS/schemas/inputs/topologies/pushPull.json -S ../MAS/schemas/inputs/topologies/forward.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.cs  --namespace MMM --csharp-version 6  --array-type list 

quicktype -l ts -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/inputs/topologies/currentTransformer.json -S ../MAS/schemas/inputs/topologies/buck.json -S ../MAS/schemas/inputs/topologies/isolatedBuck.json -S ../MAS/schemas/inputs/topologies/isolatedBuckBoost.json -S ../MAS/schemas/inputs/topologies/boost.json -S ../MAS/schemas/inputs/topologies/pushPull.json -S ../MAS/schemas/inputs/topologies/forward.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.ts --nice-property-names --acronym-style camel --converters all-objects --prefer-types

quicktype -l c++ -s schema ../schemas/flyback.json -S ../MAS/schemas/utils.json -o MAS/Flyback.hpp --namespace OpenMagnetics::Topologies --source-style single-source --type-style pascal-
case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost


valgrind --tool=callgrind ./MKF_tests [test]
kcachegrind