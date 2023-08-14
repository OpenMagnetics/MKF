
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
sudo apt install -y gcc-12 g++-12
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 10

sudo apt-get install libboost-all-dev
sudo apt install libeigen3-dev


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
cmake --build .
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

quicktype -l c++ -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/coil.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/magnetic/wire.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/outputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operatingConditions.json -S ../../MAS/schemas/inputs/operatingPointExcitation.json -o MAS/MAS.hpp --namespace OpenMagnetics --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost


quicktype -l python -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/coil.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/magnetic/wire.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operatingConditions.json -S ../../MAS/schemas/inputs/operatingPointExcitation.json -S ../../MAS/schemas/outputs.json -o MAS/MAS.py --python-version 3.7 --just-types --no-nice-property-names

quicktype -l C# -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/coil.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/magnetic/wire.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operatingConditions.json -S ../../MAS/schemas/inputs/operatingPointExcitation.json -S ../../MAS/schemas/outputs.json -o MAS/MAS.cs  --namespace MMM --csharp-version 6  --array-type list 

quicktype -l TypeScript -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/coil.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/magnetic/wire.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operatingConditions.json -S ../../MAS/schemas/inputs/operatingPointExcitation.json -S ../../MAS/schemas/outputs.json -o MAS/MAS.ts --nice-property-names --acronym-style camel --converters all-objects


valgrind --tool=callgrind ./MKF_tests [test]
kcachegrind