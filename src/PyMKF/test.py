import pprint
import PyMKF


constants = PyMKF.Constants()
pprint.pprint(constants.residual_gap)
pprint.pprint(PyMKF.get_constants())
