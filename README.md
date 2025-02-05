# Magnetics Knowledge Foundation
Repository for holding all the different models, ideally and initially writen in C++


This is the simulation engine of OpenMagnetics, and the core of other modules, like PyMKF or MKFNet.

## List of default models used

### Air gap reluctance
ZHANG - https://sci-hub.wf/https://ieeexplore.ieee.org/document/9332553

### Core Losses
Depending on which models are supported by the material, in the following order:
1) iGSE - http://inductor.thayerschool.org/papers/IGSE.pdf
2) Proprietary (Magnetics or Micrometals)
3) LOSS_FACTOR
4) STEINMETZ
5) ROSHEN - https://sci-hub.wf/10.1109/20.278656 and https://sci-hub.wf/10.1109/TPEL.2006.886608

### Core temperature
Not really used, but for completeness: MANIKTALA - Switching Power Supplies A - Z, 2nd edition by Sanjaya Maniktala, page 154

### Turns Magnetic Field
BINNS_LAWRENSON - https://library.lol/main/78A1F466FA9F3EFBCB6165283FC346B6 Equation 3.34 and 5.4

### Fringing Field Magnetic Field
ROSHEN - https://sci-hub.wf/10.1109/tmag.2007.898908 and https://sci-hub.st/10.1109/tmag.2008.2002302

### Skin effect losses
Divided per wire type:
* Round: ALBACH - Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften and https://sci-hub.wf/10.1109/tpel.2011.2143729
* Litz: ALBACH - Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften and https://sci-hub.wf/10.1109/tpel.2011.2143729
* Rectangular/Planar: KUTKUT - https://sci-hub.wf/10.1109/63.712319
* Foil: KUTKUT - https://sci-hub.wf/10.1109/63.712319

### Proximity effect losses
Divided per wire type:
* Round: LAMMERANER - https://archive.org/details/eddycurrents0000lamm
* Litz: FERREIRA - https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268
* Rectangular/Planar: ALBACH - Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften
* Foil: FERREIRA - https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268

### Capacitance
ALBACH - Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften and https://sci-hub.st/https://ieeexplore.ieee.org/document/793378

### Leakage Inductance
Own model, integrating the magnetic field calculated previously
