import pandas
import ndjson
import pathlib
import numpy
import re
import requests
import os
import urllib3
from io import StringIO
from fp.fp import FreeProxy


materials = [15, 20, 31, 43, 44, 46, 51, 52, 61, 67, 68, 73, 75, 76, 77, 78, 79, 80, 95, 96, 97, 98]
# materials = [96]

# Complex real permeability
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/fair-rite_complex_permeability.csv"):
    complex_permeability_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/fair-rite_complex_permeability.csv")
else:
    complex_permeability_data = pandas.DataFrame()
    measurements_dates = ["2015/04", "2015/03", "2019/07", "2018/06", "2024/10"]
    endings = ["Material_for-publish-1.csv", "material-for-publish.csv", "Material-revised-0824_for-publish-1.csv", "Material-Fair-Rite.csv"]

    for material in materials:
        found = False
        for ending in endings:
            for date in measurements_dates:
                try:
                    url = f"https://www.fair-rite.com/wp-content/uploads/{date}/{material}-{ending}"
                    s = requests.get(url).content
                    print(url)
                    if ending.endswith("csv"):
                        s = s.decode("windows-1252")
                        data = pandas.read_csv(StringIO(s), encoding='unicode_escape', skiprows=1)
                    else:
                        assert 0

                    found = True
                    break
                except UnicodeDecodeError:
                    print("error")
                    continue
                except ValueError:
                    print("error")
                    continue
            if found:
                break

        if data.columns[0] != "Frequency(Hz)":
            while data.values[0][0] != "Frequency(Hz)":
                data = data.drop([0]).reset_index(drop=True)

            data.columns = data.values[0]
            data = data.drop([0]).reset_index(drop=True)

            if len(data.columns) > 3:
                data = data.dropna(axis=1)

        complex_permeability_datum = data.rename(columns={
            "Frequency(Hz)": "frequency",
            "m\'": "real",
            "µ\'": "real",
            "æ\'": "real",
            "æs\'": "real",
            "µ\'\'": "imaginary",
            "m\'\'": "imaginary",
            "æ\'\'": "imaginary",
            "æs\'\'": "imaginary",
        })
        complex_permeability_datum["material"] = [material] * len(complex_permeability_datum.index)
        complex_permeability_datum = complex_permeability_datum.dropna()
        print(complex_permeability_datum)
        complex_permeability_data = pandas.concat([complex_permeability_data, complex_permeability_datum], ignore_index=True)

    complex_permeability_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/fair-rite_complex_permeability.csv", index=False)

mas_advanced_data = []

for material in materials:
    # if (material != "N27"):
    #     continue
    mas_advanced_datum = {
        "name": material,
        "manufacturerInfo": {"name": "fair-rite"},
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

    material_complex_permeability_data = complex_permeability_data[complex_permeability_data['material'] == material]
    for row_index, row in material_complex_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["complex"]["real"].append({
            "frequency": row["frequency"],
            "value": row["real"]
        })
        mas_advanced_datum["permeability"]["complex"]["imaginary"].append({
            "frequency": row["frequency"],
            "value": row["imaginary"]
        })

    mas_advanced_data.append(mas_advanced_datum)

# print(mas_data)
# print(mas_advanced_data)

# ndjson.dump(mas_data, open(f"{pathlib.Path(__file__).parent.resolve()}/fair-rite_core_materials.ndjson", "w"), ensure_ascii=False)
ndjson.dump(mas_advanced_data, open(f"{pathlib.Path(__file__).parent.resolve()}/fair-rite_advanced_core_materials.ndjson", "w"), ensure_ascii=False)
