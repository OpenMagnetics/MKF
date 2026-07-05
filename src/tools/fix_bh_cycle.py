import pandas
import ndjson
import os
import matplotlib.pyplot as plt


current_advanced_core_materials_path = "/home/alf/OpenMagnetics/MAS/data/advanced_core_materials.ndjson"

current_advanced_core_materials = pandas.DataFrame(ndjson.load(open(current_advanced_core_materials_path, "r")))
current_advanced_core_materials = current_advanced_core_materials.where(pandas.notnull(current_advanced_core_materials), None)


print(current_advanced_core_materials)

for row_index, row in current_advanced_core_materials.iterrows():
    # plt.plot([x["magneticField"] for x in row["bhCycle"]], [x["magneticFluxDensity"] for x in row["bhCycle"]])
    # plt.show()
    new_bh = []
    if row["manufacturerInfo"]["name"] == "Ferroxcube":
        for bh_point_index, bh_point in enumerate(row["bhCycle"]):
            if bh_point_index < len(row["bhCycle"]) - 1:
                next_bh_point = row["bhCycle"][bh_point_index + 1]
                if bh_point["temperature"] != next_bh_point["temperature"]:
                    continue
                if (bh_point["magneticFluxDensity"] > next_bh_point["magneticFluxDensity"] or abs(bh_point["magneticFluxDensity"]) > 0.6) and abs(bh_point["magneticFluxDensity"]) > 0.05:
                    new_bh_point = {
                        'magneticFluxDensity': bh_point['magneticFluxDensity'] / 1000,
                        'magneticField': bh_point['magneticField'],
                        'temperature': bh_point['temperature']
                    }

                    # plt.plot([x["magneticField"] for x in row["bhCycle"]], [x["magneticFluxDensity"] for x in row["bhCycle"]])
                    # plt.show()
                    # print(bh_point["magneticFluxDensity"])
                    # print(bh_point["temperature"])
                    # print(next_bh_point["magneticFluxDensity"])
                    # print(next_bh_point["temperature"])
                    # print(row["name"])
                    # assert 0
                else:
                    if abs(bh_point["magneticFluxDensity"]) > 0.6:
                        print(bh_point)
                        print(next_bh_point)
                        print(row["name"])
                        assert 0
                    new_bh_point = bh_point
                print(new_bh_point)
                new_bh.append(new_bh_point)

            else:
                new_bh.append(bh_point)

        plt.plot([x["magneticField"] for x in new_bh], [x["magneticFluxDensity"] for x in new_bh])
        plt.show()
        assert 0
