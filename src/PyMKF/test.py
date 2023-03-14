import pprint
import PyMKF
# import matplotlib
# from matplotlib import pyplot as plt
# matplotlib.use('TkAgg')

# data = [7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,7.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5]
# plt.plot(data)
# plt.show()
# assert 0

# constants = PyMKF.get_constants()
# pprint.pprint(constants['residualGap'])
# pprint.pprint(constants['minimumNonResidualGap'])
# pprint.pprint(PyMKF.get_constants())


# coreData = {'functionalDescription': {'bobbin': None,
#                            'gapping': [{'area': None,
#                                         'coordinates': None,
#                                         'distanceClosestNormalSurface': None,
#                                         'length': 0.001,
#                                         'sectionDimensions': None,
#                                         'shape': None,
#                                         'type': 'subtractive'},
#                                        {'area': None,
#                                         'coordinates': None,
#                                         'distanceClosestNormalSurface': None,
#                                         'length': 1e-05,
#                                         'sectionDimensions': None,
#                                         'shape': None,
#                                         'type': 'residual'},
#                                        {'area': None,
#                                         'coordinates': None,
#                                         'distanceClosestNormalSurface': None,
#                                         'length': 1e-05,
#                                         'sectionDimensions': None,
#                                         'shape': None,
#                                         'type': 'residual'}],
#                            'material': '3C97',
#                            'name': 'default',
#                            'numberStacks': 1,
#                            'shape': {'aliases': [],
#                                      'dimensions': {'A': 0.0391,
#                                                     'B': 0.0198,
#                                                     'C': 0.0125,
#                                                     'D': 0.0146,
#                                                     'E': 0.030100000000000002,
#                                                     'F': 0.0125,
#                                                     'G': 0.0,
#                                                     'H': 0.0},
#                                      'family': 'etd',
#                                      'familySubtype': '1',
#                                      'name': 'Custom',
#                                      'type': 'custom'},
#                            'type': 'two-piece set'},
#  'geometricalDescription': None,
#  'processedDescription': None}

# core_data = PyMKF.get_core_data(coreData)
# pprint.pprint(core_data)

# coreData = {'functionalDescription': {'bobbin': None,
#                            'gapping': [{'area': 0.000123,
#                                         'coordinates': [0.0, 0.0005, 0.0],
#                                         'distanceClosestNormalSurface': 0.0141,
#                                         'length': 0.001,
#                                         'sectionDimensions': [0.0125, 0.0125],
#                                         'shape': 'round',
#                                         'type': 'subtractive'},
#                                        {'area': 6.2e-05,
#                                         'coordinates': [0.017301, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.014595,
#                                         'length': 1e-05,
#                                         'sectionDimensions': [0.004501, 0.0125],
#                                         'shape': 'irregular',
#                                         'type': 'residual'},
#                                        {'area': 6.2e-05,
#                                         'coordinates': [-0.017301, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.014595,
#                                         'length': 1e-05,
#                                         'sectionDimensions': [0.004501, 0.0125],
#                                         'shape': 'irregular',
#                                         'type': 'residual'}],
#                            'material': '3C97',
#                            'name': 'default',
#                            'numberStacks': 1,
#                            'shape': {'aliases': [],
#                                      'dimensions': {'A': 0.0125,
#                                                     'B': 0.0064,
#                                                     'C': 0.0088,
#                                                     'D': 0.0046,
#                                                     'E': 0.01,
#                                                     'F': 0.0043,
#                                                     'G': 0.0,
#                                                     'H': 0.0,
#                                                     'K': 0.0023},
#                                      'family': 'ep',
#                                      'familySubtype': '1',
#                                      'name': 'Custom',
#                                      'type': 'custom'},
#                            'type': 'two-piece set'},
#  'geometricalDescription': None,
#  'processedDescription': None}

# core_data = PyMKF.get_core_data(coreData)
# pprint.pprint(core_data)


# coreData = {'functionalDescription': {'bobbin': None,
#                            'gapping': [{'area': 1.5e-05,
#                                         'coordinates': [0.0, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.00455,
#                                         'length': 0.0001,
#                                         'sectionDimensions': [0.0043, 0.0043],
#                                         'shape': 'round',
#                                         'type': 'subtractive'},
#                                        {'area': 8.8e-05,
#                                         'coordinates': [0.0, 0.0, 0.00805],
#                                         'distanceClosestNormalSurface': 0.004598,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.058628,
#                                                               0.001501],
#                                         'shape': 'irregular',
#                                         'type': 'residual'}],
#                            'material': '3C97',
#                            'name': 'default',
#                            'numberStacks': 1,
#                            'shape': {'aliases': [],
#                                      'dimensions': {'A': 0.101,
#                                                     'B': 0.076,
#                                                     'C': 0.03,
#                                                     'D': 0.048,
#                                                     'E': 0.044,
#                                                     'G': 0.0,
#                                                     'H': 0.0},
#                                      'family': 'u',
#                                      'familySubtype': '1',
#                                      'name': 'Custom',
#                                      'type': 'custom'},
#                            'type': 'two-piece set'},
#  'geometricalDescription': None,
#  'processedDescription': None}

# core_data = PyMKF.get_core_data(coreData)
# pprint.pprint(core_data)


# coreData = {"functionalDescription":{"name":"default","type":"two-piece set","material":"3C97","shape":{"aliases":[],"dimensions":{"A":0.0595,"B":0.036,"C":0.017,"D":0.0215,"H":0,"G":0},"family":"ur","familySubtype":"1","name":"UR 59/36/17","type":"standard"},"gapping":[{"area":0.000083,"coordinates":[0,0,0],"distanceClosestNormalSurface":0.0015,"length":0.001,"sectionDimensions":[0.0063,0.0145],"shape":"oblong","type":"subtractive"},{"area":0.000043,"coordinates":[0.01145,0,0],"distanceClosestNormalSurface":0.001996,"length":0.00001,"sectionDimensions":[0.002101,0.02],"shape":"rectangular","type":"residual"},{"area":0.000043,"coordinates":[-0.01145,0,0],"distanceClosestNormalSurface":0.001996,"length":0.00001,"sectionDimensions":[0.002101,0.02],"shape":"rectangular","type":"residual"}],"numberStacks":1,"bobbin":None},"geometricalDescription":None,"processedDescription":None}

# pprint.pprint("mierda")
# core_data = PyMKF.get_core_data(coreData)
# pprint.pprint(core_data['processedDescription']['columns'])
# pprint.pprint(len(core_data['functionalDescription']['gapping']))


# coreGap = {'gapping': [{'area': 0.000123,
#               'coordinates': [0, 0.0007, 0],
#               'distanceClosestNormalSurface': 0.01405,
#               'distanceClosestParallelSurface': 0.01405,
#               'length': 0.0011,
#               'sectionDimensions': [0.0125, 0.0125],
#               'shape': 'round',
#               'type': 'subtractive'},
#              {'area': 6.2e-05,
#               'coordinates': [0.017301, 0.00733, 0],
#               'distanceClosestNormalSurface': 0.006801,
#               'distanceClosestParallelSurface': 0.01405,
#               'length': 0.001,
#               'sectionDimensions': [0.004501, 0.0125],
#               'shape': 'irregular',
#               'type': 'subtractive'},
#              {'area': 6.2e-05,
#               'coordinates': [-0.017301, 0, 0],
#               'distanceClosestNormalSurface': 0.014598,
#               'distanceClosestParallelSurface': 0.01405,
#               'length': 5e-06,
#               'sectionDimensions': [0.004501, 0.0125],
#               'shape': 'irregular',
#               'type': 'residual'},
#              {'area': 6.2e-05,
#               'coordinates': [0.017301, 0, 0],
#               'distanceClosestNormalSurface': 0.007245,
#               'distanceClosestParallelSurface': 0.01405,
#               'length': 0.00011,
#               'sectionDimensions': [0.004501, 0.0125],
#               'shape': 'irregular',
#               'type': 'subtractive'},
#              {'area': 6.2e-05,
#               'coordinates': [0.017301, -0.00733, 0],
#               'distanceClosestNormalSurface': 0.00725,
#               'distanceClosestParallelSurface': 0.01405,
#               'length': 0.0001,
#               'sectionDimensions': [0.004501, 0.0125],
#               'shape': 'irregular',
#               'type': 'subtractive'}],
#  'model': 'Stenglein'}

# reluctance_data = PyMKF.get_gap_reluctance(coreGap['gapping'], "STENGLEIN")
# pprint.pprint(reluctance_data)


# coreGap = {
#     'area': 0.00017,
#     'coordinates': [0, 0.0002, 0],
#     'distanceClosestNormalSurface': 0.0104,
#     'height': 0,
#     'length': 0.001,
#     'sectionDimensions': [0.0147, 0.0147],
#     'shape': 'round',
#     'type': 'subtractive'
# }

# reluctance_data = PyMKF.get_gap_reluctance(coreGap, "ZHANG")
# pprint.pprint(reluctance_data)

# models_info = PyMKF.get_gap_reluctance_model_information()
# pprint.pprint(models_info)
# models = {'gapReluctance': 'CLASSIC'}
# core = {'functionalDescription': {'bobbin': None,
#                            'gapping': [{'area': 0.000369,
#                                         'coordinates': [0.0, 0.05, 0.0],
#                                         'distanceClosestNormalSurface': -0.077551,
#                                         'distanceClosestParallelSurface': 0.011524999999999999,
#                                         'length': 0.1,
#                                         'sectionDimensions': [0.02165, 0.02165],
#                                         'shape': 'round',
#                                         'type': 'subtractive'},
#                                        {'area': 0.000184,
#                                         'coordinates': [0.026126, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.022448,
#                                         'distanceClosestParallelSurface': 0.011524999999999999,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.007551,
#                                                               0.02165],
#                                         'shape': 'irregular',
#                                         'type': 'residual'},
#                                        {'area': 0.000184,
#                                         'coordinates': [-0.026126, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.022448,
#                                         'distanceClosestParallelSurface': 0.011524999999999999,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.007551,
#                                                               0.02165],
#                                         'shape': 'irregular',
#                                         'type': 'residual'}],
#                            'material': '3C95',
#                            'name': 'My Core',
#                            'numberStacks': 1,
#                            'shape': {'aliases': ['ETD 54'],
#                                      'dimensions': {'A': 0.0545,
#                                                     'B': 0.0276,
#                                                     'C': 0.0189,
#                                                     'D': 0.0202,
#                                                     'E': 0.0412,
#                                                     'F': 0.0189},
#                                      'family': 'etd',
#                                      'familySubtype': None,
#                                      'magneticCircuit': 'open',
#                                      'name': 'ETD 54/28/19',
#                                      'type': 'standard'},
#                            'type': 'two-piece set'},
#  'geometricalDescription': None,
#  'processedDescription': None}
# winding = {'functionalDescription': [{'isolationSide': 'primary',
#                             'name': 'Primary',
#                             'numberParallels': 1,
#                             'numberTurns': 1,
#                             'wire': 'Dummy'}],
#  'layersDescription': None,
#  'sectionsDescription': None,
#  'turnsDescription': None}
# inputs = {'designRequirements': {'altitude': None,
#                         'cti': None,
#                         'insulationType': None,
#                         'leakageInductance': None,
#                         'magnetizingInductance': {'excludeMaximum': None,
#                                                   'excludeMinimum': None,
#                                                   'maximum': None,
#                                                   'minimum': None,
#                                                   'nominal': 0.004654652816558039},
#                         'name': None,
#                         'operationTemperature': None,
#                         'overvoltageCategory': None,
#                         'pollutionDegree': None,
#                         'turnsRatios': []},
#  'operationPoints': [{'conditions': {'ambientRelativeHumidity': None,
#                                      'ambientTemperature': 25.0,
#                                      'cooling': None,
#                                      'name': None},
#                       'excitationsPerWinding': [{'current': {'harmonics': None,
#                                                              'processed': None,
#                                                              'waveform': {'ancillaryLabel': None,
#                                                                           'data': [41.0,
#                                                                                    51.0,
#                                                                                    41.0],
#                                                                           'numberPeriods': None,
#                                                                           'time': [0.0,
#                                                                                    2.4999999999999998e-06,
#                                                                                    1e-05]}},
#                                                  'frequency': 100000.0,
#                                                  'magneticField': None,
#                                                  'magneticFluxDensity': None,
#                                                  'magnetizingCurrent': None,
#                                                  'name': 'My Operation Point',
#                                                  }],
#                       'name': None}]}
# gappingType = 'GRINDED'

# gapping = PyMKF.get_gapping_from_number_turns_and_inductance(core,
#                                                              winding,
#                                                              inputs,
#                                                                 gappingType,
#                                                                 4,
#                                                                 models)
# print(gapping)
# print(len(gapping))

            # {"frequency", 800000},
            # {"magneticFluxDensityPeak", 0.1},
            # {"magneticFluxDensityDutyCycle", 0.5},
            # {"temperature", 25},
            # {"magneticFluxDensityOffset", 0},
            # {"magneticFieldStrengthDc", 0},
            # {"expectedVolumetricLosses", 1020000.0},

# coeff = PyMKF.get_steinmetz_coefficients("N49", 800000)

# voluemtriclosses = coeff["k"] * 800000**coeff["alpha"] * 0.1**coeff["beta"] * (coeff["ct2"] * 25**2 - coeff["ct1"] * 25 + coeff["ct0"])
# print(coeff)
# print(voluemtriclosses)
# frequencies = [1000000, 500000, 800000, 300000, 400000, 50000, 25000, 100000, 200000]
# materials = ["3F4", "N49", "N49", "3C94", "N27", "N87", "3C90"]
# comb = {}


# std::map<std::string, std::map<int, std::map<std::string, double>>> dynamicCoefficients = {
#     {"3C90",
#         {
#             {25000, {{"alpha", 0.6119032165032708}, {"beta", 2.6414923616949553}, {"ct0", 0.0016607172720359237}, {"ct1", 2.388168471291204e-05}, {"ct2", 1.2273149176903466e-07}, {"k", 25358950.162134618}}},
#             {50000, {{"alpha", 0.6119032165032708}, {"beta", 2.6414923616949553}, {"ct0", 0.0016607172720359237}, {"ct1", 2.388168471291204e-05}, {"ct2", 1.2273149176903466e-07}, {"k", 25358950.162134618}}},
#             {100000, {{"alpha", 1.6922655116583307}, {"beta", 2.6387558212068556}, {"ct0", 0.007035640971373052}, {"ct1", 0.00010466617163345241}, {"ct2", 5.584896170038324e-07}, {"k", 33.12013020258462}}},
#             {200000, {{"alpha", 2.3888138991940306}, {"beta", 2.319170498242314}, {"ct0", 0.0013316374491110425}, {"ct1", 1.2215818567857081e-05}, {"ct2", 6.436835159667792e-08}, {"k", 0.01008697340219455}}},
#             {300000, {{"alpha", 3.207909730760055}, {"beta", 2.208742992181126}, {"ct0", 2.524783110701414e-06}, {"ct1", 1.3817591590708868e-08}, {"ct2", 8.471092014959255e-11}, {"k", 9.16521541023768e-05}}},
#             {400000, {{"alpha", 3.163524246279912}, {"beta", 2.140623673581521}, {"ct0", 8.029803690795726e-07}, {"ct1", 1.5929730782382e-09}, {"ct2", 1.2078361153444373e-11}, {"k", 0.0002309379771197434}}},
#             {500000, {{"alpha", 4.491346273001152}, {"beta", 2.1332489188659687}, {"ct0", 7.90326808296919e-11}, {"ct1", 1.1865327830297114e-13}, {"ct2", 9.007097262275786e-16}, {"k", 6.866630262826199e-08}}},
#             {800000, {{"alpha", 4.491346273001152}, {"beta", 2.1332489188659687}, {"ct0", 7.90326808296919e-11}, {"ct1", 1.1865327830297114e-13}, {"ct2", 9.007097262275786e-16}, {"k", 6.866630262826199e-08}}},
#             {1000000, {{"alpha", 4.491346273001152}, {"beta", 2.1332489188659687}, {"ct0", 7.90326808296919e-11}, {"ct1", 1.1865327830297114e-13}, {"ct2", 9.007097262275786e-16}, {"k", 6.866630262826199e-08}}},
#         },
#     },
# };


# for material in materials:

#     aux = {}
#     for frequency in frequencies:
#         aux[frequency] = PyMKF.get_steinmetz_coefficients(material, frequency)
#     comb[material] = aux

# # import pprint
# # pprint.pprint(comb)

# mapStr = "std::map<std::string, std::map<int, std::map<std::string, double>>> dynamicCoefficients = {\n"

# for material in materials:
#     mapStr += '\t{' + f'"{material}",\n' + '\t\t{\n'
#     for frequency in frequencies:
#         print(str(comb[material][frequency]))
#         aux = '\t{' + str(frequency) + ', {'
#         for key, value in comb[material][frequency].items():
#             aux += '{' + f'"{key}", {value}' + '}, '
#         aux += '}' + '}'
#         print(aux)
#         mapStr += f'\t\t{aux},\n'
#     mapStr += '\t\t},\n\t},\n' 
# mapStr += '};' 
# print(mapStr)

# models = {'coreLosses': 'STEINMETZ', 'gapReluctance': 'BALAKRISHNAN'}
# core = {'functionalDescription': {'gapping': [{'area': 9.8e-05,
#                                         'coordinates': [0.0, 0.0001, 0.0],
#                                         'distanceClosestNormalSurface': 0.011301,
#                                         'distanceClosestParallelSurface': 0.006999999999999999,
#                                         'length': 0.0002,
#                                         'sectionDimensions': [0.0092, 0.01065],
#                                         'shape': 'rectangular',
#                                         'type': 'subtractive'},
#                                        {'area': 4.7e-05,
#                                         'coordinates': [0.0138, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.011498,
#                                         'distanceClosestParallelSurface': 0.006999999999999999,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.004401,
#                                                               0.01065],
#                                         'shape': 'rectangular',
#                                         'type': 'residual'},
#                                        {'area': 4.7e-05,
#                                         'coordinates': [-0.0138, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.011498,
#                                         'distanceClosestParallelSurface': 0.006999999999999999,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.004401,
#                                                               0.01065],
#                                         'shape': 'rectangular',
#                                         'type': 'residual'}],
#                            'material': 'XFlux 19',
#                            'numberStacks': 1,
#                            'shape': {'aliases': [],
#                                      'dimensions': {'A': 0.032,
#                                                     'B': 0.0161,
#                                                     'C': 0.01065,
#                                                     'D': 0.0115,
#                                                     'E': 0.0232,
#                                                     'F': 0.0092,
#                                                     'G': 0.0,
#                                                     'H': 0.0},
#                                      'family': 'e',
#                                      'familySubtype': None,
#                                      'magneticCircuit': 'open',
#                                      'name': 'E 32/16/11',
#                                      'type': 'standard'},
#                            'type': 'two-piece set'},
#  'manufacturerInfo': None,
#  'name': 'My Core'}
# winding = {'bobbin': None,
#  'functionalDescription': [{'isolationSide': 'primary',
#                             'name': 'Primary',
#                             'numberParallels': 1,
#                             'numberTurns': 33,
#                             'wire': 'Dummy'}],
#  'layersDescription': None,
#  'sectionsDescription': None,
#  'turnsDescription': None}
# inputs = {'conditions': {'ambientRelativeHumidity': None,
#                 'ambientTemperature': 37.0,
#                 'cooling': None,
#                 'name': None},
#  'excitationsPerWinding': [{'current': {'harmonics': None,
#                                         'processed': None,
#                                         'waveform': {'ancillaryLabel': None,
#                                                      'data': [-5.0, 5.0, -5.0],
#                                                      'numberPeriods': None,
#                                                      'time': [0.0,
#                                                               2.4999999999999998e-06,
#                                                               1e-05]}},
#                             'frequency': 100000.0,
#                             'magneticFieldStrength': None,
#                             'magneticFluxDensity': None,
#                             'magnetizingCurrent': None,
#                             'name': 'My Operation Point',
#                             'voltage': {'harmonics': None,
#                                         'processed': None,
#                                         'waveform': {'ancillaryLabel': None,
#                                                      'data': [700.5,
#                                                               700.5,
#                                                               -200.4999999999999996,
#                                                               -200.4999999999999996,
#                                                               700.5],
#                                                      'numberPeriods': None,
#                                                      'time': [0.0,
#                                                               2.4999999999999998e-06,
#                                                               2.4999999999999998e-06,
#                                                               1e-05,
#                                                               1e-05]}}}],
#  'name': None}

# core_losses_result = PyMKF.get_core_losses(core,
#                                            winding,
#                                            inputs,
#                                            models)

# print(core_losses_result)
# 'models'
# models = {'coreLosses': 'ROSHEN', 'gapReluctance': 'BALAKRISHNAN'}
# core = {'functionalDescription': {'gapping': [{'area': 5.6e-05,
#                                         'coordinates': [0.0, 0.0005, 0.0],
#                                         'distanceClosestNormalSurface': 0.00205,
#                                         'distanceClosestParallelSurface': 0.004475000000000001,
#                                         'length': 0.001,
#                                         'sectionDimensions': [0.0084, 0.0084],
#                                         'shape': 'round',
#                                         'type': 'subtractive'},
#                                        {'area': 3.3e-05,
#                                         'coordinates': [0.010026, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.003048,
#                                         'distanceClosestParallelSurface': 0.004475000000000001,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.0027, 0.012223],
#                                         'shape': 'irregular',
#                                         'type': 'residual'},
#                                        {'area': 3.3e-05,
#                                         'coordinates': [-0.010026, 0.0, 0.0],
#                                         'distanceClosestNormalSurface': 0.003048,
#                                         'distanceClosestParallelSurface': 0.004475000000000001,
#                                         'length': 5e-06,
#                                         'sectionDimensions': [0.0027, 0.012223],
#                                         'shape': 'irregular',
#                                         'type': 'residual'}],
#                            'material': 'N87',
#                            'numberStacks': 1,
#                            'shape': {'aliases': ['RM 8LP'],
#                                      'dimensions': {'A': 0.02275,
#                                                     'B': 0.00575,
#                                                     'C': 0.0108,
#                                                     'D': 0.00305,
#                                                     'E': 0.01735,
#                                                     'F': 0.0084,
#                                                     'G': 0.009500000000000001,
#                                                     'H': 0.0045,
#                                                     'J': 0.0193,
#                                                     'R': 0.00030000000000000003},
#                                      'family': 'rm',
#                                      'familySubtype': '3',
#                                      'magneticCircuit': 'open',
#                                      'name': 'RM 8/11',
#                                      'type': 'standard'},
#                            'type': 'two-piece set'},
#  'manufacturerInfo': None,
#  'name': 'My Core'}
# winding = {'bobbin': None,
#  'functionalDescription': [{'isolationSide': 'primary',
#                             'name': 'Primary',
#                             'numberParallels': 1,
#                             'numberTurns': 34,
#                             'wire': 'Dummy'}],
#  'layersDescription': None,
#  'sectionsDescription': None,
#  'turnsDescription': None}
# inputs = {'conditions': {'ambientRelativeHumidity': None,
#                 'ambientTemperature': 36.0,
#                 'cooling': None,
#                 'name': None},
#  'excitationsPerWinding': [{'current': {'harmonics': None,
#                                         'processed': None,
#                                         'waveform': {'ancillaryLabel': None,
#                                                      'data': [0.0, 10.0, 0.0],
#                                                      'numberPeriods': None,
#                                                      'time': [0.0,
#                                                               2.4999999999999998e-06,
#                                                               1e-05]}},
#                             'frequency': 100000.0,
#                             'magneticFieldStrength': None,
#                             'magneticFluxDensity': None,
#                             'magnetizingCurrent': None,
#                             'name': 'My Operation Point',
#                             'voltage': None}],
#  'name': None}


# core_losses_result = PyMKF.get_core_losses(core,
#                                            winding,
#                                            inputs,
#                                            models)

# import matplotlib.pyplot as plt
# print(core_losses_result.keys())
# print(core_losses_result['totalLosses'])
# print(core_losses_result['magneticFluxDensityAcPeak'])

# plt.plot(core_losses_result['_hysteresisMajorH'], core_losses_result['_hysteresisMajorLoopBottom'])
# plt.plot(core_losses_result['_hysteresisMajorH'], core_losses_result['_hysteresisMajorLoopTop'])
# plt.show()

# plt.plot(core_losses_result['_hysteresisMinorLoopBottom'])
# plt.plot(core_losses_result['_hysteresisMinorLoopTop'])
# plt.show()


# inputs = {'conditions': {'ambientRelativeHumidity': None,
#                 'ambientTemperature': 37.0,
#                 'cooling': None,
#                 'name': None},
#  'excitationsPerWinding': [{'current': {'harmonics': None,
#                                         'processed': None,
#                                         'waveform': {'ancillaryLabel': None,
#                                                      'data': [0.0, 10.0, 0.0],
#                                                      'numberPeriods': None,
#                                                      'time': [0.0,
#                                                               2.4999999999999998e-06,
#                                                               1e-05]}},
#                             'frequency': 100000.0,
#                             'magneticFieldStrength': None,
#                             'magneticFluxDensity': None,
#                             'magnetizingCurrent': None,
#                             'name': 'My Operation Point',
#                             'voltage': None}],
#  'name': None}


# core_losses_result = PyMKF.get_core_losses(core,
#                                            winding,
#                                            inputs,
#                                            models)

# # print(core_losses_result)
# print(core_losses_result['totalLosses'])
# print(core_losses_result.keys())
# print(core_losses_result['totalLosses'])
# print(core_losses_result['magneticFluxDensityAcPeak'])

# plt.plot(core_losses_result['_hysteresisMajorH'], core_losses_result['_hysteresisMajorLoopBottom'])
# plt.plot(core_losses_result['_hysteresisMajorH'], core_losses_result['_hysteresisMajorLoopTop'])
# plt.show()

# plt.plot(core_losses_result['_hysteresisMinorLoopBottom'])
# plt.plot(core_losses_result['_hysteresisMinorLoopTop'])
# plt.show()
