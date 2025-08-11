import numpy
import json
import ndjson
import math
import PyMKF
import pandas
import pathlib

data_path = pathlib.Path(__file__).parent.resolve().joinpath('../../../MAS/data/advanced_core_materials.ndjson')
advanced_core_materials = pandas.read_json(data_path, lines=True)

data_path = pathlib.Path(__file__).parent.resolve().joinpath('../../../MAS/data/core_materials.ndjson')
core_materials = pandas.read_json(data_path, lines=True)

print(core_materials)

frequencies_ranges = [[[1, 100e3], [100e3, 500e5], [500e5, 1e9]],
                      [[1, 150e3], [150e3, 400e5], [400e5, 1e9]],
                      [[1, 150e3], [150e3, 1e6], [1e6, 1e9]],
                      [[1, 500e3], [500e3, 1e9]],
                      [[1, 1e6], [1e6, 1e9]],
                      [[1, 3e6], [3e6, 1e9]],
                      [[1, 1e9]]]


core_materials = core_materials.to_dict("records")
for core_material_index, core_material in enumerate(core_materials):
    advanced_core_material = advanced_core_materials[advanced_core_materials["name"] == core_material["name"]]

    volumetricLossesData = advanced_core_material['volumetricLosses']
    if volumetricLossesData.empty:
        continue

    volumetricLossesData = volumetricLossesData.iloc[0]
    # print(core_material["name"])
    # print(volumetricLossesData)
    if 'default' not in volumetricLossesData:
        continue

    for method in volumetricLossesData['default']:
        if isinstance(method, list):
            if len(method) == 0:
                continue
            # print(core_material["name"])
            
            best_result = None
            best_result_error = math.inf
            for frequencies_range in frequencies_ranges:
                result = PyMKF.calculate_steinmetz_coefficients_with_error(method, frequencies_range)
                average_error = sum(result["errorPerRange"]) / len(result["errorPerRange"])

                valid = True
                for coefficient_range in result["coefficientsPerRange"]:
                    if coefficient_range["k"] < 0 or coefficient_range["alpha"] < 0 or coefficient_range["beta"] < 0:
                        valid = False
                if not valid:
                    continue

                if average_error < best_result_error:
                    best_result = result
                    best_result_error = average_error
            # print(best_result)
            # print(sum(best_result["errorPerRange"]) / len(best_result["errorPerRange"]))
            for method_index, method_aux in enumerate(core_material['volumetricLosses']['default']):
                if isinstance(method_aux, dict) and method_aux["method"] == "steinmetz":
                    if best_result_error < 0.6:
                        core_materials[core_material_index]['volumetricLosses']['default'][method_index]["ranges"] = best_result["coefficientsPerRange"]
                    break

for core_material_index, core_material in enumerate(core_materials):
    key_to_remove = []
    for key, value in core_material.items():
        if value is numpy.nan:
            key_to_remove.append(key)

    for key in key_to_remove:
        del core_materials[core_material_index][key]

ndjson.dump(core_materials, open(f"{pathlib.Path(__file__).parent.resolve()}/new_core_materials.ndjson", "w"), ensure_ascii=False)
