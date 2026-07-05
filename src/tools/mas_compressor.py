import pandas
import json
import pprint
import pathlib
import ndjson


def heavy_compress_magnetic(magnetic: dict) -> dict:

    if isinstance(magnetic["core"]["functionalDescription"]["material"], dict):
        magnetic["core"]["functionalDescription"]["material"] = magnetic["core"]["functionalDescription"]["material"]["name"]
    if isinstance(magnetic["core"]["functionalDescription"]["shape"], dict):
        magnetic["core"]["functionalDescription"]["shape"] = magnetic["core"]["functionalDescription"]["shape"]["name"]

    if "processedDescription" in magnetic["core"]:
        del magnetic["core"]["processedDescription"]
    if "geometricalDescription" in magnetic["core"]:
        del magnetic["core"]["geometricalDescription"]

    for gap_index in range(len(magnetic["core"]["functionalDescription"]["gapping"])):
        magnetic["core"]["functionalDescription"]["gapping"][gap_index] = {
            "length": magnetic["core"]["functionalDescription"]["gapping"][gap_index]["length"],
            "type": magnetic["core"]["functionalDescription"]["gapping"][gap_index]["type"],
        }

    if "layersDescription" in magnetic["coil"]:
        del magnetic["coil"]["layersDescription"]
    if "turnsDescription" in magnetic["coil"]:
        del magnetic["coil"]["turnsDescription"]

    for winding_index in range(len(magnetic["coil"]["functionalDescription"])):
        if isinstance(magnetic["coil"]["functionalDescription"][winding_index]["wire"], dict) and "name" in magnetic["coil"]["functionalDescription"][winding_index]["wire"] and magnetic["coil"]["functionalDescription"][winding_index]["wire"]["name"] is not None:
            magnetic["coil"]["functionalDescription"][winding_index]["wire"] = magnetic["coil"]["functionalDescription"][winding_index]["wire"]["name"]

    return magnetic


def manufacturerInfo(row):
    row["magnetic"]["manufacturerInfo"] = {}
    row["magnetic"]["manufacturerInfo"]["name"] = "Midcom"
    row["magnetic"]["manufacturerInfo"]["reference"] = row["part_number"]
    return row["magnetic"]


data_path = "/home/alf/OpenMagnetics/MKF/output/magnetics.csv"
data = pandas.read_csv(data_path)


data["magnetic"] = data.apply(lambda x: heavy_compress_magnetic(json.loads(x["magnetic"])), axis=1)
# data["magnetic"] = data.apply(lambda x: json.loads(x["magnetic"]), axis=1)
data["magnetic"] = data.apply(lambda x: manufacturerInfo(x), axis=1)

# data.to_csv(data_path, index=False)

ndjson.dump(data["magnetic"].tolist(), open(f"{pathlib.Path(__file__).parent.resolve()}/magnetics.ndjson", "w"), ensure_ascii=False)
# for index, row in data.iterrows():

#     try: 
#         # print(index)
#         magnetic = data.iloc[index]["magnetic"]
#         magnetic = json.loads(magnetic)
#         # magnetic = json.loads(magnetic)
#         data.loc[index, "magnetic"] = heavy_compress_magnetic(magnetic)
#     except ValueError:
#         pprint.pprint(data.iloc[index]["part_number"])
#     #     # magnetic = json.loads(magnetic)
#         assert 0
