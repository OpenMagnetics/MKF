
## Build steps

1. Create a build directory:

    ```
    $ mkdir build && cd build
    ```

2. Configure the CMake project (Using Ninja in this example):

    ```
    $ cmake .. -G "Ninja"
    ```

4. Build it:

    ```
    $ cmake --build .
    ```

5. Run the application:

    ```
    $ ./MKF_tests
    ```
6. Create Python package

    ```
    python3 -m pip install -e ../ -vvv
    ```


wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt install cmake
sudo apt install ninja-build
npm install -g quicktype
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install -y gcc-11 g++-11
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 10

quicktype -l c++ -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/winding.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/core/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operationConditions.json -S ../../MAS/schemas/inputs/operationPointExcitation.json -o MAS/MAS.hpp --namespace OpenMagnetics --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost


quicktype -l python -s schema ../../MAS/schemas/MAS.json -S ../../MAS/schemas/magnetic.json -S ../../MAS/schemas/magnetic/core.json -S ../../MAS/schemas/magnetic/winding.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/core/gap.json -S ../../MAS/schemas/magnetic/core/shape.json -S ../../MAS/schemas/magnetic/core/material.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/core/bobbin.json -S ../../MAS/schemas/magnetic/core/piece.json -S ../../MAS/schemas/magnetic/core/spacer.json -S ../../MAS/schemas/magnetic/wire/solid.json -S ../../MAS/schemas/magnetic/wire/stranded.json -S ../../MAS/schemas/magnetic/wire/material.json -S ../../MAS/schemas/utils.json -S ../../MAS/schemas/magnetic/insulation/wireCoating.json -S ../../MAS/schemas/magnetic/insulation/material.json -S ../../MAS/schemas/inputs.json -S ../../MAS/schemas/inputs/designRequirements.json -S ../../MAS/schemas/inputs/operationConditions.json -S ../../MAS/schemas/inputs/operationPointExcitation.json -o MAS/MAS.py --python-version 3.7 --just-types --no-nice-property-names