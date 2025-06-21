import pandas
import ndjson
import pathlib
import numpy
import math


spreadsheet_path = "/mnt/c/Users/Alfonso/Downloads/mmcurvefitcoefficientsall.xlsx"


def autocomplete_materials(data, column):
    current_name = None
    for row_index, row in data.iterrows():
        if row[column] is not None:
            current_name = row[column]
        else:
            data.at[row_index, column] = current_name
    return data


def drop_blocks(data, column):
    current_name = None
    for row_index, row in data.iterrows():
        if row[column] is not None:
            current_name = row[column]
        else:
            data.at[row_index, column] = current_name
    return data[(~data[column].str.contains('B')) & (~data[column].str.contains('OD'))]


def get_material_and_variants(datum):
    if isinstance(datum['Material'], int):
        material = f"Mix {datum['Material']}"
        family = "Mix"
    else:
        material = f"{datum['Material']} {int(datum['Init. Perm. (µ)'])}"
        family = datum['Material']

    if "T" in datum['PartType']:
        variant = "default"
    elif "EQ" in datum['PartType']:
        variant = "EQ"
    elif "E" in datum['PartType']:
        variant = "E"
    elif "PQ" in datum['PartType']:
        variant = "PQ"

    return material, family, int(datum['Init. Perm. (µ)']), variant


def find_material(data, name):
    for row in data:
        if name == row["name"]:
            return row
    else:
        # return None
        assert 0, f"Material {name} not found"


headers_row = 8

permeability_vs_bias_indexes = "A,B,C,F:I"
permeability_vs_bias_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_bias_indexes)
permeability_vs_bias_data = permeability_vs_bias_data.where(pandas.notnull(permeability_vs_bias_data), None)
permeability_vs_bias_data = permeability_vs_bias_data.dropna(subset="a")
permeability_vs_bias_data = autocomplete_materials(data=permeability_vs_bias_data,
                                                   column="Material")
permeability_vs_bias_data = drop_blocks(data=permeability_vs_bias_data,
                                        column="PartType")

permeability_vs_b_peak_indexes = "A,B,C,N:S"
permeability_vs_b_peak_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_b_peak_indexes)
permeability_vs_b_peak_data.columns = ['Material', 'PartType', 'Init. Perm. (µ)', 'a', 'b', 'c', 'd', 'e', 'f']
permeability_vs_b_peak_data = permeability_vs_b_peak_data.where(pandas.notnull(permeability_vs_b_peak_data), None)
permeability_vs_b_peak_data = permeability_vs_b_peak_data.dropna(subset="a")
permeability_vs_b_peak_data = autocomplete_materials(data=permeability_vs_b_peak_data,
                                                     column="Material")
permeability_vs_b_peak_data = drop_blocks(data=permeability_vs_b_peak_data,
                                          column="PartType")

bh_cycle_indexes = "A,B,C,D,E,X:AB"
bh_cycle_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=bh_cycle_indexes)
bh_cycle_data.columns = ['Material', 'PartType', 'Init. Perm. (µ)', 'Bsat(G)', 'Density(g/cm³)', 'a', 'b', 'c', 'd', 'e']
bh_cycle_data = bh_cycle_data.where(pandas.notnull(bh_cycle_data), None)
bh_cycle_data = autocomplete_materials(data=bh_cycle_data, column="Material")
bh_cycle_data = drop_blocks(data=bh_cycle_data, column="PartType")

core_losses_density_indexes = "A,B,C,J:M"
core_losses_density_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=core_losses_density_indexes)
core_losses_density_data = core_losses_density_data.where(pandas.notnull(core_losses_density_data), None)
core_losses_density_data.columns = core_losses_density_data.columns.str.replace('.\\d', '', regex=True)
core_losses_density_data = core_losses_density_data.dropna(subset="a")
core_losses_density_data = autocomplete_materials(data=core_losses_density_data,
                                                  column="Material")
core_losses_density_data = drop_blocks(data=core_losses_density_data,
                                       column="PartType")
permeability_vs_frequency_indexes = "A,B,C,T:W"
permeability_vs_frequency_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_frequency_indexes)
permeability_vs_frequency_data = permeability_vs_frequency_data.where(pandas.notnull(permeability_vs_frequency_data), None)
permeability_vs_frequency_data.columns = permeability_vs_frequency_data.columns.str.replace('.\\d', '', regex=True)
permeability_vs_frequency_data = permeability_vs_frequency_data.dropna(subset="a")
permeability_vs_frequency_data = autocomplete_materials(data=permeability_vs_frequency_data,
                                                        column="Material")
permeability_vs_frequency_data = drop_blocks(data=permeability_vs_frequency_data,
                                             column="PartType")
permeability_vs_temperature_indexes = "A,B,C,AC:AG"
permeability_vs_temperature_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_temperature_indexes)
permeability_vs_temperature_data = permeability_vs_temperature_data.where(pandas.notnull(permeability_vs_temperature_data), None)
permeability_vs_temperature_data.columns = permeability_vs_temperature_data.columns.str.replace('.\\d', '', regex=True)
permeability_vs_temperature_data = permeability_vs_temperature_data.dropna(subset="a")
permeability_vs_temperature_data = autocomplete_materials(data=permeability_vs_temperature_data,
                                                          column="Material")
permeability_vs_temperature_data = drop_blocks(data=permeability_vs_temperature_data,
                                               column="PartType")


oersted_to_ampere_per_meter = 79.5774715459
permeability_vs_bias_data['b'] = permeability_vs_bias_data.apply(lambda row: row['b'] / math.pow(oersted_to_ampere_per_meter, row['c']), axis=1)

mW_cm3_to_W_m3 = 1e-3
gauss_to_tesla = 1e-4
core_losses_density_data['a'] = core_losses_density_data.apply(lambda row: row['a'] * math.pow(gauss_to_tesla, 3) * mW_cm3_to_W_m3, axis=1)
core_losses_density_data['b'] = core_losses_density_data.apply(lambda row: row['b'] * math.pow(gauss_to_tesla, 2.3) * mW_cm3_to_W_m3, axis=1)
core_losses_density_data['c'] = core_losses_density_data.apply(lambda row: row['c'] * math.pow(gauss_to_tesla, 1.65) * mW_cm3_to_W_m3, axis=1)
core_losses_density_data['d'] = core_losses_density_data.apply(lambda row: row['d'] / math.pow(gauss_to_tesla, 2) / mW_cm3_to_W_m3, axis=1)


permeability_vs_b_peak_data['b'] = permeability_vs_b_peak_data.apply(lambda row: row['b'] / math.pow(gauss_to_tesla, row['c']), axis=1)
permeability_vs_b_peak_data['d'] = permeability_vs_b_peak_data.apply(lambda row: row['d'] / math.pow(gauss_to_tesla, row['e']), axis=1)


materials = []
advanced_materials = []
materials_to_ignore = ["Mix 5", "SM 14", "SM 90", "SP 14", "SP 26", "SP 40", "SP 60", "SP 75", "SP 90"]

materials_file = open(f"{pathlib.Path(__file__).parent.parent.parent.resolve()}/MAS/data/core_materials.ndjson", "r")
current_materials = ndjson.load(materials_file)

for row_index, row in permeability_vs_bias_data.iterrows():
    material, family, permeability, variant = get_material_and_variants(row)
    if material in materials_to_ignore:
        continue
    material_data = find_material(current_materials, material)

    if variant == 'default':
        material_data["permeability"] = {
            "initial": {
                "value": permeability,
                "modifiers": {
                }
            }
        }
        material_data["volumetricLosses"] = {}

        materials.append(material_data)

for row_index, row in permeability_vs_bias_data.iterrows():
    material, family, permeability, variant = get_material_and_variants(row)
    if material in materials_to_ignore:
        continue
    if material == "Mix 1":
        print(material)
        print(family)
        print(permeability)
        print(variant)

    material_data = find_material(current_materials, material)

    for elem_index, elem in enumerate(materials):
        if elem["name"] == material:
            current_index = elem_index
            break
    else:
        materials.append(material_data)
        current_index = -1
    if material == "Mix 1":
        print(materials[current_index])

    if variant not in materials[current_index]['permeability']['initial']['modifiers']:
        materials[current_index]['permeability']['initial']['modifiers'][variant] = {
            "method": "micrometals",
            "magneticFieldDcBiasFactor": {"a": row['a'], "b": row['b'], "c": row['c'], "d": row['d']}
        }

    permeability_vs_b_peak_this_point = permeability_vs_b_peak_data[(permeability_vs_b_peak_data['Material'] == row['Material']) & (permeability_vs_b_peak_data['Init. Perm. (µ)'] == row['Init. Perm. (µ)'])]
    if variant == "default":
        permeability_vs_b_peak_this_point = permeability_vs_b_peak_this_point[permeability_vs_b_peak_this_point['PartType'].apply(lambda x: "T" in x)]
    elif variant == "E":
        permeability_vs_b_peak_this_point = permeability_vs_b_peak_this_point[permeability_vs_b_peak_this_point['PartType'].apply(lambda x: x == "E")]
    elif variant == "EQ":
        permeability_vs_b_peak_this_point = permeability_vs_b_peak_this_point[permeability_vs_b_peak_this_point['PartType'].apply(lambda x: x == "EQ")]
    elif variant == "PQ":
        permeability_vs_b_peak_this_point = permeability_vs_b_peak_this_point[permeability_vs_b_peak_this_point['PartType'].apply(lambda x: x == "PQ")]
    else:
        assert 0

    if not permeability_vs_b_peak_this_point.empty:
        assert len(permeability_vs_b_peak_this_point.index) == 1
        permeability_vs_b_peak_this_point = permeability_vs_b_peak_this_point.iloc[0]
        materials[current_index]['permeability']['initial']['modifiers'][variant]["magneticFluxDensityFactor"] = {
            "a": permeability_vs_b_peak_this_point["a"],
            "b": permeability_vs_b_peak_this_point["b"],
            "c": permeability_vs_b_peak_this_point["c"],
            "d": permeability_vs_b_peak_this_point["d"],
            "e": permeability_vs_b_peak_this_point["e"],
            "f": permeability_vs_b_peak_this_point["f"]
        }

    bh_cycle_this_point = bh_cycle_data[(bh_cycle_data['Material'] == row['Material']) & (bh_cycle_data['Init. Perm. (µ)'] == row['Init. Perm. (µ)'])]
    if variant == "default":
        bh_cycle_this_point = bh_cycle_this_point[bh_cycle_this_point['PartType'].apply(lambda x: "T" in x)]
    elif variant == "E":
        bh_cycle_this_point = bh_cycle_this_point[bh_cycle_this_point['PartType'].apply(lambda x: x == "E")]
    elif variant == "EQ":
        bh_cycle_this_point = bh_cycle_this_point[bh_cycle_this_point['PartType'].apply(lambda x: x == "EQ")]
    elif variant == "PQ":
        bh_cycle_this_point = bh_cycle_this_point[bh_cycle_this_point['PartType'].apply(lambda x: x == "PQ")]
    else:
        assert 0

    if not bh_cycle_this_point.empty:
        assert len(bh_cycle_this_point.index) == 1
        bh_cycle_this_point = bh_cycle_this_point.iloc[0]

        advanced_materials.append({
            "name": materials[current_index]["name"],
            "manufacturerInfo": materials[current_index]["manufacturerInfo"],
            "permeability": {"amplitude": []},
            "volumetricLosses": {"default": [[]]},
            "bhCycle": []
        })

        a = bh_cycle_this_point["a"]
        b = bh_cycle_this_point["b"]
        c = bh_cycle_this_point["c"]
        d = bh_cycle_this_point["d"]
        e = bh_cycle_this_point["e"]
        mu = bh_cycle_this_point["Init. Perm. (µ)"]
        prev_B = math.inf
        for H in numpy.arange(1, 10000, 100):
            B = mu / (1 / (H + a * (H**b)) + 1 / (c * H**d) + 1 / e)                
            point = {
                "magneticFluxDensity": round(float(B) / 10000, 5),
                "magneticField": round(float(H) * oersted_to_ampere_per_meter),
                "temperature": 25
            }
            if abs(B - prev_B) / B < 0.001:
                break
            prev_B = B
            advanced_materials[-1]['bhCycle'].append(point)

        materials[current_index]["saturation"] = [{
            "magneticFluxDensity": float(B) / 10000,
            "magneticField": float(H) * oersted_to_ampere_per_meter,
            "temperature": 25
        }]
        materials[current_index]["density"] = bh_cycle_this_point["Density(g/cm³)"] * 1000

    permeability_vs_frequency_this_point = permeability_vs_frequency_data[(permeability_vs_frequency_data['Material'] == row['Material']) & (permeability_vs_frequency_data['Init. Perm. (µ)'] == row['Init. Perm. (µ)'])]
    if variant == "default":
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point[permeability_vs_frequency_this_point['PartType'].apply(lambda x: "T" in x)]
    elif variant == "E":
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point[permeability_vs_frequency_this_point['PartType'].apply(lambda x: x == "E")]
    elif variant == "EQ":
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point[permeability_vs_frequency_this_point['PartType'].apply(lambda x: x == "EQ")]
    elif variant == "PQ":
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point[permeability_vs_frequency_this_point['PartType'].apply(lambda x: x == "PQ")]
    else:
        assert 0

    if not permeability_vs_frequency_this_point.empty:
        assert len(permeability_vs_frequency_this_point.index) == 1
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point.iloc[0]
        materials[current_index]['permeability']['initial']['modifiers'][variant]["frequencyFactor"] = {
            "a": permeability_vs_frequency_this_point["a"],
            "b": permeability_vs_frequency_this_point["b"],
            "c": permeability_vs_frequency_this_point["c"],
            "d": permeability_vs_frequency_this_point["d"]
        }

    permeability_vs_temperature_this_point = permeability_vs_temperature_data[(permeability_vs_temperature_data['Material'] == row['Material']) & (permeability_vs_temperature_data['Init. Perm. (µ)'] == row['Init. Perm. (µ)'])]
    if variant == "default":
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point[permeability_vs_temperature_this_point['PartType'].apply(lambda x: "T" in x)]
    elif variant == "E":
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point[permeability_vs_temperature_this_point['PartType'].apply(lambda x: x == "E")]
    elif variant == "EQ":
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point[permeability_vs_temperature_this_point['PartType'].apply(lambda x: x == "EQ")]
    elif variant == "PQ":
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point[permeability_vs_temperature_this_point['PartType'].apply(lambda x: x == "PQ")]

    if not permeability_vs_temperature_this_point.empty:
        assert len(permeability_vs_temperature_this_point.index) == 1
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point.iloc[0]
        if numpy.isnan(permeability_vs_temperature_this_point["b"]):
            materials[current_index]['permeability']['initial']['modifiers'][variant]["temperatureFactor"] = {
                "a": permeability_vs_temperature_this_point["a"]
            }
        else:
            materials[current_index]['permeability']['initial']['modifiers'][variant]["temperatureFactor"] = {
                "a": permeability_vs_temperature_this_point["a"],
                "b": permeability_vs_temperature_this_point["b"],
                "c": permeability_vs_temperature_this_point["c"],
                "d": permeability_vs_temperature_this_point["d"],
                "e": permeability_vs_temperature_this_point["e"]
            }

    core_loss_data_this_point = core_losses_density_data[(core_losses_density_data['Material'] == row['Material']) & (core_losses_density_data['Init. Perm. (µ)'] == row['Init. Perm. (µ)'])]
    if variant == "default":
        core_loss_data_this_point = core_loss_data_this_point[core_loss_data_this_point['PartType'].apply(lambda x: "T" in x)]
    elif variant == "E":
        core_loss_data_this_point = core_loss_data_this_point[core_loss_data_this_point['PartType'].apply(lambda x: x == "E")]
    elif variant == "EQ":
        core_loss_data_this_point = core_loss_data_this_point[core_loss_data_this_point['PartType'].apply(lambda x: x == "EQ")]
    elif variant == "PQ":
        core_loss_data_this_point = core_loss_data_this_point[core_loss_data_this_point['PartType'].apply(lambda x: x == "PQ")]

    if core_loss_data_this_point.empty:
        continue

    assert len(core_loss_data_this_point.index) == 1
    core_loss_data_this_point = core_loss_data_this_point.iloc[0]

    if variant not in materials[current_index]['volumetricLosses']:
        materials[current_index]['volumetricLosses'][variant] = [
            {
                "method": "micrometals",
                "a": core_loss_data_this_point['a'],
                "b": core_loss_data_this_point['b'],
                "c": core_loss_data_this_point['c'],
                "d": core_loss_data_this_point['d']
            }
        ]

out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/micrometals_materials.ndjson", "w")
ndjson.dump(materials, out_file, ensure_ascii=False)
out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/micrometals_advanced_materials.ndjson", "w")
ndjson.dump(advanced_materials, out_file, ensure_ascii=False)


for micrometals_material in materials:
    for material_index, material in enumerate(current_materials):
        if micrometals_material["name"] == material["name"]:
            current_materials[material_index] = micrometals_material
            break
    else:
        current_materials.append(micrometals_material)

out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/updated_core_materials.ndjson", "w")
ndjson.dump(current_materials, out_file, ensure_ascii=False)
