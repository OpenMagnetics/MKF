
## Build steps

1. Install dependencies:

```
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt install cmake
sudo apt install ninja-build
nvm install node
or
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - && sudo apt-get install -y nodejs

npm install -g quicktype
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install -y gcc-14 g++-14
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 10
sudo apt install -y gnuplot
sudo apt install libopenblas-dev

if node installation fails:
ln -s /usr/bin/python3 /usr/include/python3

```

2. Create a build directory:

```
mkdir build && cd build
```

3. Configure the CMake project (Using Ninja in this example):

```
cmake .. -G "Ninja"
```

4. Build it:

```
ninja
```

5. Run the tests:

```
./MKF_tests
```
6. Create Python package

```
python3 -m pip install -e ../ -vvv
```


if quicktype/MAS generation fails in step 4, try to generate it manually:

quicktype -l c++ -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -o MAS/MAS.hpp --namespace MAS --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost


quicktype -l python -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.py --python-version 3.7 --no-nice-property-names

quicktype -l cs -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.cs  --namespace MMM --csharp-version 6  --array-type list 

quicktype -l ts -s schema ../MAS/schemas/MAS.json -S ../MAS/schemas/magnetic.json -S ../MAS/schemas/magnetic/core.json -S ../MAS/schemas/magnetic/coil.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/core/gap.json -S ../MAS/schemas/magnetic/core/shape.json -S ../MAS/schemas/magnetic/core/material.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/bobbin.json -S ../MAS/schemas/magnetic/core/piece.json -S ../MAS/schemas/magnetic/core/spacer.json -S ../MAS/schemas/magnetic/wire/basicWire.json -S ../MAS/schemas/magnetic/wire/round.json -S ../MAS/schemas/magnetic/wire/rectangular.json -S ../MAS/schemas/magnetic/wire/foil.json -S ../MAS/schemas/magnetic/wire/planar.json -S ../MAS/schemas/magnetic/wire/litz.json -S ../MAS/schemas/magnetic/wire/material.json -S ../MAS/schemas/magnetic/wire.json -S ../MAS/schemas/utils.json -S ../MAS/schemas/magnetic/insulation/wireCoating.json -S ../MAS/schemas/magnetic/insulation/material.json -S ../MAS/schemas/inputs.json -S ../MAS/schemas/inputs/designRequirements.json -S ../MAS/schemas/inputs/operatingPoint.json -S ../MAS/schemas/inputs/operatingConditions.json -S ../MAS/schemas/inputs/operatingPointExcitation.json -S ../MAS/schemas/inputs/topologies/flyback.json -S ../MAS/schemas/outputs.json -S ../MAS/schemas/outputs/coreLossesOutput.json -o MAS/MAS.ts --nice-property-names --acronym-style camel --converters all-objects --prefer-types

quicktype -l c++ -s schema ../schemas/flyback.json -S ../MAS/schemas/utils.json -o MAS/Flyback.hpp --namespace OpenMagnetics::Topologies --source-style single-source --type-style pascal-
case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost


valgrind --tool=callgrind ./MKF_tests [test]
kcachegrind