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
