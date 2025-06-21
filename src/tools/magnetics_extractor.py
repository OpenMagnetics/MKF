import pandas
import ndjson
import pathlib
import math
import numpy


spreadsheet_path = "/mnt/c/Users/Alfonso/Downloads/Curve-Fit-Equation-Tool 11-28-23.xlsx"


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
    return data[~data[column].str.contains('Block')]


def get_material_and_variants(datum):
    datum['Material'] = datum['Material'].replace("\n", " ")
    datum['Material'] = datum['Material'].replace("E-Core", "E-core")
    datum['Material'] = datum['Material'].replace("\"", "")
    if "E-core" in datum['Material']:
        material = datum['Material'].split(" E-core")[0]
        materials = [f"{material} {datum['Permeability']} default"]
        return f"{material} {datum['Permeability']}", material, datum['Permeability'], "E/ER/U"
    elif "EQ" in datum['Material']:
        material = datum['Material'].split(" EQ")[0]
        return f"{material} {datum['Permeability']}", material, datum['Permeability'], "EQ/LP"
    else:
        return f"{datum['Material']} {datum['Permeability']}", datum['Material'], datum['Permeability'], "default"


curieTemperaturePerFamily = {
    'Kool Mµ': 500,
    'Kool Mµ MAX': 500,
    'Kool Mµ Max A': 500,
    'Kool Mµ Max Y': 500,
    'Kool Mµ Hƒ': 500,
    'Kool Mµ Ultra': 500,
    'XFlux': 700,
    'XFlux Ultra': 700,
    'High Flux': 500,
    'High DC Bias XFlux': 500,
    'Edge': 500,
    'High DC Bias Edge': 500,
    'MPP': 460,
    '75-Series': 700,
}
saturationPerFamily = {
    'Kool Mµ': 1,
    'Kool Mµ MAX': 1,
    'Kool Mµ Max A': 1,
    'Kool Mµ Max Y': 1,
    'Kool Mµ Hƒ': 1,
    'Kool Mµ Ultra': 1,
    'XFlux': 1.6,
    'XFlux Ultra': 1.6,
    'High Flux': 1.5,
    'High DC Bias XFlux': 1.5,
    'Edge': 1.5,
    'High DC Bias Edge': 1.5,
    'MPP': 0.8,
    '75-Series': 0.9,
}

headers_row = 7

permeability_vs_bias_indexes = "B:H"
permeability_vs_bias_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_bias_indexes)
permeability_vs_bias_data = permeability_vs_bias_data.where(pandas.notnull(permeability_vs_bias_data), None)
permeability_vs_bias_data = permeability_vs_bias_data.dropna(subset="a")
permeability_vs_bias_data = autocomplete_materials(data=permeability_vs_bias_data,
                                                   column="Material")
permeability_vs_bias_data = drop_blocks(data=permeability_vs_bias_data,
                                        column="Material")

core_losses_density_indexes = "J:O"
core_losses_density_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=core_losses_density_indexes)
core_losses_density_data = core_losses_density_data.where(pandas.notnull(core_losses_density_data), None)
core_losses_density_data.columns = core_losses_density_data.columns.str.replace('.\\d', '', regex=True)
core_losses_density_data = core_losses_density_data.dropna(subset="a")
core_losses_density_data = autocomplete_materials(data=core_losses_density_data,
                                                  column="Material")
core_losses_density_data = drop_blocks(data=core_losses_density_data,
                                       column="Material")

permeability_vs_frequency_indexes = "AA:AI"
permeability_vs_frequency_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_frequency_indexes)
permeability_vs_frequency_data = permeability_vs_frequency_data.where(pandas.notnull(permeability_vs_frequency_data), None)
permeability_vs_frequency_data.columns = permeability_vs_frequency_data.columns.str.replace('.\\d', '', regex=True)
permeability_vs_frequency_data = permeability_vs_frequency_data.dropna(subset="a")
permeability_vs_frequency_data = autocomplete_materials(data=permeability_vs_frequency_data,
                                                        column="Material")
permeability_vs_frequency_data = drop_blocks(data=permeability_vs_frequency_data,
                                             column="Material")

permeability_vs_temperature_indexes = "AK:AS"
permeability_vs_temperature_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=permeability_vs_temperature_indexes)
permeability_vs_temperature_data = permeability_vs_temperature_data.where(pandas.notnull(permeability_vs_temperature_data), None)
permeability_vs_temperature_data.columns = permeability_vs_temperature_data.columns.str.replace('.\\d', '', regex=True)
permeability_vs_temperature_data = permeability_vs_temperature_data.dropna(subset="a")
permeability_vs_temperature_data = autocomplete_materials(data=permeability_vs_temperature_data,
                                                          column="Material")
permeability_vs_temperature_data = drop_blocks(data=permeability_vs_temperature_data,
                                               column="Material")

bh_cycle_indexes = "Q:X"
bh_cycle_data = pandas.read_excel(spreadsheet_path, header=headers_row, usecols=bh_cycle_indexes)
bh_cycle_data = bh_cycle_data.where(pandas.notnull(bh_cycle_data), None)
bh_cycle_data.columns = bh_cycle_data.columns.str.replace('.\\d', '', regex=True)
bh_cycle_data = bh_cycle_data.dropna(subset="a")
bh_cycle_data = autocomplete_materials(data=bh_cycle_data, column="Material")
bh_cycle_data = drop_blocks(data=bh_cycle_data, column="Material")

oersted_to_ampere_per_meter = 79.5774715459
permeability_vs_bias_data['b'] = permeability_vs_bias_data.apply(lambda row: row['b'] / math.pow(oersted_to_ampere_per_meter, row['c']), axis=1)

kilohertz_to_hertz = 1000
core_losses_density_data['a'] = core_losses_density_data.apply(lambda row: 1000 * row['a'] / math.pow(kilohertz_to_hertz, row['c']), axis=1)

megahertz_to_hertz = 1000000
permeability_vs_frequency_data['b'] = permeability_vs_frequency_data.apply(lambda row: row['b'] / math.pow(megahertz_to_hertz, 1), axis=1)
permeability_vs_frequency_data['c'] = permeability_vs_frequency_data.apply(lambda row: row['c'] / math.pow(megahertz_to_hertz, 2), axis=1)
permeability_vs_frequency_data['d'] = permeability_vs_frequency_data.apply(lambda row: row['d'] / math.pow(megahertz_to_hertz, 3), axis=1)
permeability_vs_frequency_data['e'] = permeability_vs_frequency_data.apply(lambda row: row['e'] / math.pow(megahertz_to_hertz, 4), axis=1)

materials = []
advanced_materials = []


for row_index, row in permeability_vs_bias_data.iterrows():
    material, family, permeability, variant = get_material_and_variants(row)
    if variant == 'default':
        core_loss_data_this_point = core_losses_density_data[(core_losses_density_data['Material'] == row['Material']) & (core_losses_density_data['Permeability'] == row['Permeability'])]
        if core_loss_data_this_point.empty:
            continue

        template = {
            "type": "commercial",
            "name": material,
            "family": family,
            "resistivity": [{"value": 100, "temperature": 25}],
            "heatConductivity": {"nominal": 8},
            "curieTemperature": curieTemperaturePerFamily[family],
            "materialComposition": "powder",
            "manufacturerInfo": {"name": "Magnetics"},
            "permeability": {
                "initial": {
                    "value": permeability,
                    "modifiers": {
                    }
                }
            },
            "saturation": [{"magneticFluxDensity": saturationPerFamily[family], "magneticField": 7957, "temperature": 100.0}],
            "volumetricLosses": {
            }
        }
        materials.append(template)
        current_index = -1
    else:
        for elem_index, elem in enumerate(materials):
            if elem["name"] == material:
                current_index = elem_index
                break
        else:
            continue
            # assert 0, "Missing material"

    if variant not in materials[current_index]['permeability']['initial']['modifiers']:
        materials[current_index]['permeability']['initial']['modifiers'][variant] = {
            "method": "magnetics",
            "magneticFieldDcBiasFactor": {"a": row['a'], "b": row['b'], "c": row['c']}
        }

    permeability_vs_frequency_this_point = permeability_vs_frequency_data[(permeability_vs_frequency_data['Material'] == row['Material']) & (permeability_vs_frequency_data['Permeability'] == row['Permeability'])]
    if not permeability_vs_frequency_this_point.empty:
        assert len(permeability_vs_frequency_this_point.index) == 1
        permeability_vs_frequency_this_point = permeability_vs_frequency_this_point.iloc[0]
        materials[current_index]['permeability']['initial']['modifiers'][variant]["frequencyFactor"] = {
            "a": permeability_vs_frequency_this_point["a"],
            "b": permeability_vs_frequency_this_point["b"],
            "c": permeability_vs_frequency_this_point["c"],
            "d": permeability_vs_frequency_this_point["d"],
            "e": permeability_vs_frequency_this_point["e"]
        }

    permeability_vs_temperature_this_point = permeability_vs_temperature_data[(permeability_vs_temperature_data['Material'] == row['Material']) & (permeability_vs_temperature_data['Permeability'] == row['Permeability'])]
    if not permeability_vs_temperature_this_point.empty:
        assert len(permeability_vs_temperature_this_point.index) == 1
        permeability_vs_temperature_this_point = permeability_vs_temperature_this_point.iloc[0]
        materials[current_index]['permeability']['initial']['modifiers'][variant]["temperatureFactor"] = {
            "a": permeability_vs_temperature_this_point["a"],
            "b": permeability_vs_temperature_this_point["b"],
            "c": permeability_vs_temperature_this_point["c"],
            "d": permeability_vs_temperature_this_point["d"],
            "e": permeability_vs_temperature_this_point["e"]
        }

    bh_cycle_this_point = bh_cycle_data[(bh_cycle_data['Material'] == row['Material']) & (bh_cycle_data['Permeability'] == row['Permeability'])]
    if not bh_cycle_this_point.empty:
        # assert len(bh_cycle_this_point.index) == 1
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
        x = bh_cycle_this_point["x"]
        prev_B = math.inf
        for H in numpy.arange(1, 10000, 100):
            B = ((a + b * H + c * H**2) / (1 + d * H + e * H**2))**x                
            point = {
                "magneticFluxDensity": round(float(B), 5),
                "magneticField": round(float(H) * oersted_to_ampere_per_meter),
                "temperature": 25
            }
            if abs(B - prev_B) / B < 0.001:
                break
            prev_B = B
            advanced_materials[-1]['bhCycle'].append(point)

    core_loss_data_this_point = core_losses_density_data[(core_losses_density_data['Material'] == row['Material']) & (core_losses_density_data['Permeability'] == row['Permeability'])]
    if core_loss_data_this_point.empty:
        continue

    assert len(core_loss_data_this_point.index) == 1
    core_loss_data_this_point = core_loss_data_this_point.iloc[0]

    if variant not in materials[current_index]['volumetricLosses']:
        materials[current_index]['volumetricLosses'][variant] = [
            {"method": "magnetics", "a": core_loss_data_this_point['a'], "b": core_loss_data_this_point['b'], "c": core_loss_data_this_point['c']}
        ]

out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/magnetics_materials.ndjson", "w")
ndjson.dump(materials, out_file, ensure_ascii=False)


out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/magnetics_advanced_materials.ndjson", "w")
ndjson.dump(advanced_materials, out_file, ensure_ascii=False)
