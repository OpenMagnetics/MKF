import pandas
import ndjson
import pathlib
import numpy
import re


spreadsheet_path = "/mnt/c/Users/Alfonso/Downloads/Material-Spec_Template_ferrite_cores（SINOMAG).xlsx"


core_loss_headers_row = 2
core_loss_indexes = "B:U"
core_loss_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="Core loss", 
    header=core_loss_headers_row,
    usecols=core_loss_indexes
)
core_loss_data["material"] = core_loss_data["Unnamed: 1"]
core_loss_data["B"] = core_loss_data["Unnamed: 2"] / 1000
core_loss_data["f"] = core_loss_data["Unnamed: 3"] * 1000
core_loss_data = core_loss_data.drop(["Unnamed: 1", "Unnamed: 2", "Unnamed: 3"], axis=1)
core_loss_data = core_loss_data.where(pandas.notnull(core_loss_data), None)
core_loss_data = core_loss_data.replace(numpy.nan, None)
core_loss_data = core_loss_data.replace("NaN", None)
for row_index, row in core_loss_data.iterrows():
    if row["material"] is None and core_loss_data.iloc[row_index - 1]["material"] is not None:
        core_loss_data.loc[row_index, "material"] = core_loss_data.iloc[row_index - 1]["material"]

materials = list(core_loss_data["material"].unique())

initial_permeability_headers_row = 13
initial_permeability_indexes = "B:Q"
initial_permeability_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="μi-T", 
    header=initial_permeability_headers_row,
    usecols=initial_permeability_indexes
)
initial_permeability_data["material"] = initial_permeability_data["Temp./℃"]
initial_permeability_data = initial_permeability_data.drop(["Temp./℃"], axis=1)
initial_permeability_data = initial_permeability_data.dropna()


complex_real_permeability_headers_row = 1
complex_real_permeability_indexes = "B,C,E,G,I,K,M,O,Q"
complex_real_permeability_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="μs‘&μs“（25℃）", 
    header=complex_real_permeability_headers_row,
    usecols=complex_real_permeability_indexes
)
complex_real_permeability_data = complex_real_permeability_data.dropna()
complex_real_permeability_data = complex_real_permeability_data.where(pandas.notnull(complex_real_permeability_data), None)
complex_real_permeability_data = complex_real_permeability_data.replace(numpy.nan, None)
complex_real_permeability_data = complex_real_permeability_data.replace("NaN", None)
complex_real_permeability_data = complex_real_permeability_data.rename(columns={"Material": "frequency"})
complex_real_permeability_data = complex_real_permeability_data.drop([0], axis=0)

complex_imaginary_permeability_headers_row = 1
complex_imaginary_permeability_indexes = "B,D,F,H,J,L,N,P,R"
complex_imaginary_permeability_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="μs‘&μs“（25℃）", 
    header=complex_imaginary_permeability_headers_row,
    usecols=complex_imaginary_permeability_indexes
)
complex_imaginary_permeability_data = complex_imaginary_permeability_data.dropna()
complex_imaginary_permeability_data = complex_imaginary_permeability_data.where(pandas.notnull(complex_imaginary_permeability_data), None)
complex_imaginary_permeability_data = complex_imaginary_permeability_data.replace(numpy.nan, None)
complex_imaginary_permeability_data = complex_imaginary_permeability_data.replace("NaN", None)
complex_imaginary_permeability_data.columns = complex_real_permeability_data.columns
complex_imaginary_permeability_data = complex_imaginary_permeability_data.drop([0], axis=0)

amplitude_permeability_headers_row = 1
amplitude_permeability_indexes = "B:J"
amplitude_permeability_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="μa", 
    header=amplitude_permeability_headers_row,
    usecols=amplitude_permeability_indexes
)
amplitude_permeability_data = amplitude_permeability_data.rename(columns={"Material": "material", "freq(kHz)": "frequency"})
amplitude_permeability_data["frequency"] = amplitude_permeability_data["frequency"] * 1000

saturation_headers_row = 1
saturation_indexes = "B:H"
saturation_data = pandas.read_excel(
    spreadsheet_path,
    sheet_name="BH-data", 
    header=saturation_headers_row,
    usecols=saturation_indexes
)
saturation_data = saturation_data.rename(columns={"Material": "material", "f/kHz": "frequency", "Bm(mT)": "Bm(T)"})
saturation_data["frequency"] = saturation_data["frequency"] * 1000
saturation_data["Bm(T)"] = saturation_data["Bm(T)"] / 1000


materials = list(core_loss_data["material"].unique())

mas_data = []
mas_advanced_data = []

density_per_material = {
    "SMP44": 4800,
    "SMP47": 4900,
    "SMP95": 4900,
    "SMP96": 4900,
    "SMP97": 4800,
    "SMP50": 4800,
    "SMP51": 4800,
    "SMP53": 4800,
}
curie_temperature_per_material = {
    "SMP44": 215,
    "SMP47": 230,
    "SMP95": 230,
    "SMP96": 230,
    "SMP97": 230,
    "SMP50": 215,
    "SMP51": 215,
    "SMP53": 215,
}
resistivity_per_material = {
    "SMP44": 4,
    "SMP47": 4,
    "SMP95": 4,
    "SMP96": 4,
    "SMP97": 4,
    "SMP50": 6.5,
    "SMP51": 6.5,
    "SMP53": 6.5,
}

for material in materials:
    print(material)
    mas_datum = {
        "type": "commercial",
        "name": material,
        "family": "SMP",
        "resistivity": [],
        "remanence": [],
        "coerciveForce": [],
        "density": None,
        "curieTemperature": None,
        "material": "ferrite",
        "materialComposition": "NiZn",
        "manufacturerInfo": {"name": "Sinomag", "datasheetUrl": "https://www.sinomagtech.com/uploadfiles/file/202311/4.pdf"},
        "permeability": {
            "initial": []
        },
        "saturation": [],
        "volumetricLosses": {"default": []}
    }
    mas_advanced_datum = {
        "name": material,
        "manufacturerInfo": {"name": "Sinomag"},
        "permeability": {
            "amplitude": [],
            "complex": {
                "real": [],
                "imaginary": []
            }
        },
        "volumetricLosses": {"default": []},
        "bhCycle": []
    }

    # Core loss
    material_core_loss_data = core_loss_data[core_loss_data['material'] == material]
    core_loss_temperatures = material_core_loss_data.drop(["material", "B", "f"], axis=1).columns.to_list()

    mas_datum["resistivity"] = [{"value": resistivity_per_material[material], "temperature": 25}]
    mas_datum["density"] = density_per_material[material]
    mas_datum["curieTemperature"] = curie_temperature_per_material[material]

    mas_advanced_datum["volumetricLosses"]["default"].append([])
    for row_index, row in material_core_loss_data.iterrows():
        B = row["B"]
        f = row["f"]
        for temperature in core_loss_temperatures:
            if row[temperature] is None:
                continue

            mas_advanced_datum["volumetricLosses"]["default"][0].append(
                {
                    "magneticFluxDensity": {
                        "frequency": f, 
                        "magneticFluxDensity": {
                            "processed": {
                                "label": "Triangular",
                                "peak": B,
                                "offset": 0.0,
                                "dutyCycle": 0.5
                            }
                        }
                    },
                    "temperature": temperature,
                    "value": row[temperature] * 1000,
                    "origin": "manufacturer"
                })

    material_initial_permeability_data = initial_permeability_data[initial_permeability_data['material'] == material]
    initial_permeability_temperatures = material_initial_permeability_data.drop(["material"], axis=1).columns.to_list()
    # print(material_initial_permeability_data)
    for row_index, row in initial_permeability_data.iterrows():
        for temperature in initial_permeability_temperatures:
            mas_datum["permeability"]["initial"].append({
                "temperature": temperature,
                "value": row[temperature]
            })

    material_complex_real_permeability_data = complex_real_permeability_data[["frequency", material]]
    for row_index, row in material_complex_real_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["complex"]["real"].append({
            "frequency": row["frequency"],
            "value": row[material]
        })

    material_complex_imaginary_permeability_data = complex_imaginary_permeability_data[["frequency", material]]
    for row_index, row in material_complex_imaginary_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["complex"]["imaginary"].append({
            "frequency": row["frequency"],
            "value": row[material]
        })

    material_amplitude_permeability_data = amplitude_permeability_data[amplitude_permeability_data['material'] == material]
    for row_index, row in material_amplitude_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["amplitude"].append({
            "value": row["μa"],
            "frequency": row["frequency"],
            "temperature": row["temp/℃"],
            "magneticFieldPeak": row["Hm(A/m)"],
            "magneticFluxDensityPeak": row["Bm(T)"]
        })

    material_saturation_data = saturation_data[saturation_data['material'] == material]
    for row_index, row in material_saturation_data.iterrows():
        mas_datum["saturation"].append({
            "magneticFluxDensity": row["Bm(T)"],
            "magneticField": row["H/Am-1"],
            "temperature": row["Temp/℃"]
        })

    bh_loop_h_headers_row = 2
    bh_loop_h_indexes = "B,D,F,H,J,L,N,P,R"
    bh_loop_h_data = pandas.read_excel(
        spreadsheet_path,
        sheet_name=f"BH curve Data（{material}）", 
        header=bh_loop_h_headers_row,
        usecols=bh_loop_h_indexes
    )
    bh_loop_h_data = bh_loop_h_data.drop([0], axis=0)
    bh_loop_h_data = bh_loop_h_data.rename(columns=lambda x: re.sub('℃', '', x))

    bh_loop_b_headers_row = 2
    bh_loop_b_indexes = "C,E,G,I,K,M,O,Q,S"
    bh_loop_b_data = pandas.read_excel(
        spreadsheet_path,
        sheet_name=f"BH curve Data（{material}）", 
        header=bh_loop_b_headers_row,
        usecols=bh_loop_b_indexes
    )
    bh_loop_b_data = bh_loop_b_data.drop([0], axis=0)
    bh_loop_b_data.columns = bh_loop_h_data.columns

    bh_loop_temperatures = bh_loop_h_data.columns.to_list()
    coercivity_per_temperature = {}
    remanence_per_temperature = {}
    for temperature in bh_loop_temperatures:
        current_closest_h_to_0 = 999999
        current_closest_b_to_0 = 999999
        for row_index, row in bh_loop_h_data.iterrows():

            if abs(bh_loop_h_data.loc[row_index][temperature]) < current_closest_h_to_0:
                current_closest_h_to_0 = abs(bh_loop_h_data.loc[row_index][temperature])
                remanence_per_temperature[int(temperature)] = abs(bh_loop_b_data.loc[row_index][temperature])
            if abs(bh_loop_b_data.loc[row_index][temperature]) < current_closest_b_to_0:
                current_closest_b_to_0 = abs(bh_loop_b_data.loc[row_index][temperature])
                coercivity_per_temperature[int(temperature)] = abs(bh_loop_h_data.loc[row_index][temperature])

            mas_advanced_datum["bhCycle"].append({
                "magneticFluxDensity": bh_loop_b_data.loc[row_index][temperature],
                "magneticField": bh_loop_h_data.loc[row_index][temperature],
                "temperature": temperature
            })

    for temperature, value in coercivity_per_temperature.items():
        mas_datum["coerciveForce"].append({
            "magneticFluxDensity": 0,
            "magneticField": value,
            "temperature": temperature
        })

    for temperature, value in remanence_per_temperature.items():
        mas_datum["remanence"].append({
            "magneticFluxDensity": value,
            "magneticField": 0,
            "temperature": temperature
        })

    mas_data.append(mas_datum)
    mas_advanced_data.append(mas_advanced_datum)

# print(mas_data)
# print(mas_advanced_data)

ndjson.dump(mas_data, open(f"{pathlib.Path(__file__).parent.resolve()}/sinomag_core_materials.ndjson", "w"), ensure_ascii=False)
ndjson.dump(mas_advanced_data, open(f"{pathlib.Path(__file__).parent.resolve()}/sinomag_advanced_core_materials.ndjson", "w"), ensure_ascii=False)
