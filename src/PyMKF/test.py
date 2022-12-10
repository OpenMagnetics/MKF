import pprint
import PyMKF


constants = PyMKF.get_constants()
pprint.pprint(constants['residualGap'])
pprint.pprint(constants['minimumNonResidualGap'])
pprint.pprint(PyMKF.get_constants())


coreData = {
    'functionalDescription': {
        'bobbin': None,
        'gapping': [],
        'material': '3C97',
        'name': 'default',
        'numberStacks': 1,
        'shape': {'aliases': [],
                  'dimensions': {'A': 0.0308,
                                 'B': 0.0264,
                                 'C': 0.0265,
                                 'D': 0.016,
                                 'E': 0.01,
                                 'G': 0.0,
                                 'H': 0.0},
                  'family': 'u',
                  'familySubtype': '1',
                  'name': 'Custom',
                  'type': 'custom'},
        'type': 'two-piece set'
    }
}

core_data = PyMKF.get_core_data(coreData)
pprint.pprint(core_data)


coreGap = {
    'area': 0.000123,
    'coordinates': [0.0, 0.0, 0.0],
    'distanceClosestNormalSurface': 0.0146,
    'length': 0.0001,
    'sectionDimensions': [0.0125, 0.0125],
    'shape': 'round',
    'type': 'additive'
}

reluctance_data = PyMKF.get_gap_reluctance(coreGap, "ZHANG")
pprint.pprint(reluctance_data)


coreGap = {
    'area': 0.00075,
    'coordinates': [0, 0, 0],
    'distanceClosestNormalSurface': 0.0334975,
    'height': 0,
    'length': 0.00013000000000000002,
    'sectionDimensions': [0.025, 0.03],
    'shape': 'rectangular',
    'type': 'subtractive'
}

reluctance_data = PyMKF.get_gap_reluctance(coreGap, "ZHANG")
pprint.pprint(reluctance_data)

models_info = PyMKF.get_gap_reluctance_model_information()
pprint.pprint(models_info)
