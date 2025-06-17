import pandas
import ndjson
import pathlib
import numpy
import glob
import os
import re


# spreadsheet_path = "/mnt/c/Users/Alfonso/Downloads/FXC Materials Permeability.xlsx"
current_advanced_materials_path = "/home/alfonso/OpenMagnetics/MAS/data/advanced_core_materials.ndjson"

current_advanced_materials_raw = ndjson.load(open(current_advanced_materials_path, "r"))
# # print(current_advanced_materials_raw[0]["name"])
# # assert 0

# complex_real_permeability_headers_row = 47
# complex_real_permeability_indexes = "A:Q"
# complex_real_permeability_data = pandas.read_excel(
#     spreadsheet_path,
#     sheet_name="summary", 
#     nrows=48,
#     header=complex_real_permeability_headers_row,
#     usecols=complex_real_permeability_indexes
# )

# complex_real_permeability_data = complex_real_permeability_data.dropna()
# complex_real_permeability_data = complex_real_permeability_data.where(pandas.notnull(complex_real_permeability_data), None)
# complex_real_permeability_data = complex_real_permeability_data.replace(numpy.nan, None)
# complex_real_permeability_data = complex_real_permeability_data.replace("NaN", None)
# complex_real_permeability_data = complex_real_permeability_data.rename(columns={"µi'": "frequency"})
# complex_real_permeability_data["frequency"] = complex_real_permeability_data["frequency"] * 1000
# complex_real_permeability_data = complex_real_permeability_data.drop([0], axis=0)
# complex_real_permeability_data = complex_real_permeability_data.set_index("frequency", drop=True)

# complex_imaginary_permeability_headers_row = 96
# complex_imaginary_permeability_indexes = "A:Q"
# complex_imaginary_permeability_data = pandas.read_excel(
#     spreadsheet_path,
#     nrows=48,
#     sheet_name="summary", 
#     header=complex_imaginary_permeability_headers_row,
#     usecols=complex_imaginary_permeability_indexes
# )

# complex_imaginary_permeability_data = complex_imaginary_permeability_data.dropna()
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.where(pandas.notnull(complex_imaginary_permeability_data), None)
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.replace(numpy.nan, None)
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.replace("NaN", None)
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.rename(columns={"µii": "frequency"})
# complex_imaginary_permeability_data.columns = [x.replace("uii ", "") for x in list(complex_imaginary_permeability_data.columns)]
# complex_imaginary_permeability_data["frequency"] = complex_imaginary_permeability_data["frequency"] * 1000
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.drop([0], axis=0)
# complex_imaginary_permeability_data = complex_imaginary_permeability_data.set_index("frequency", drop=True)

# materials = list(complex_imaginary_permeability_data.columns)

# mas_data = []
# mas_advanced_data = []

# for material in materials:
#     data = None
#     for current_advanced_materials_data in current_advanced_materials_raw:
#         if current_advanced_materials_data["name"] == material:
#             data = current_advanced_materials_data
#     if data is None:
#         continue

#     if "permeability" not in data:
#         data["permeability"] = {
#         }

#     if "complex" not in data["permeability"]:
#         data["permeability"]["complex"] = {
#             "real": [],
#             "imaginary": []
#         }

#     data["permeability"]["complex"]["real"] = []
#     data["permeability"]["complex"]["imaginary"] = []

#     for frequency in complex_real_permeability_data[material].index:
#         datum = float(complex_real_permeability_data[material][frequency])
#         data["permeability"]["complex"]["real"].append({
#             "frequency": float(frequency),
#             "value": float(datum)
#         })

#     for frequency in complex_imaginary_permeability_data[material].index:
#         datum = float(complex_imaginary_permeability_data[material][frequency])
#         data["permeability"]["complex"]["imaginary"].append({
#             "frequency": float(frequency),
#             "value": float(datum)
#         })

#     for index in range(len(current_advanced_materials_raw)):
#         if current_advanced_materials_raw[index]["name"] == material:
#             current_advanced_materials_raw[index] = data
#             break

    # print(current_advanced_materials_raw[index]["permeability"])
    # assert 0
# print(mas_advanced_data)

# ndjson.dump(mas_data, open(f"{pathlib.Path(__file__).parent.resolve()}/ferroxcube_updated_core_materials.ndjson", "w"), ensure_ascii=False)

temperatures = [25, 100]

for file in glob.iglob('/home/alfonso/OpenMagnetics/MKF/build/BH loops Ferroxcube/*', recursive=True):

    material = file.split("_")[1]

    for current_advanced_materials_data in current_advanced_materials_raw:
        if current_advanced_materials_data["name"] == material:
            data = current_advanced_materials_data

    data["bhCycle"] = []
    for temperature in temperatures:
        saturation_headers_row = 16
        saturation_indexes = "B:C"
        saturation_data = pandas.read_excel(
            file,
            sheet_name=f"T{temperature}_f10_H250", 
            header=saturation_headers_row,
            usecols=saturation_indexes
        )

        for index, row in saturation_data.iterrows():
            point = {
                "magneticFluxDensity": float(row["B"]) / 1000,
                "magneticField": float(row["H"]),
                "temperature": temperature
            }
            data["bhCycle"].append(point)


ndjson.dump(current_advanced_materials_raw, open(f"{pathlib.Path(__file__).parent.resolve()}/ferroxcube_updated_advanced_core_materials.ndjson", "w"), ensure_ascii=False)
