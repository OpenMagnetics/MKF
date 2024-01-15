from io import StringIO
import contextlib
import requests, json
import webbrowser
import pathlib
import pprint
import io
import os
import re
import ast
import math
import pandas
import copy
import ndjson
import PyMKF
import urllib
import shutil
from datetime import date
from difflib import SequenceMatcher


class Stocker():
    def __init__(self):
        self.core_data = {}

    def verify_shape_is_in_database(self, family, shape, material, product):
        print("family")
        print(family)
        core_data = {
            "name": "temp",
            "functionalDescription": {
                "type": 'two-piece set' if family != 'T' else 'toroidal',
                "material": material,
                "shape": shape,
                "gapping": [],
                "numberStacks": 1
            }
        }
        try:
            core_data = PyMKF.get_core_data(core_data, False)
            return core_data
        except RuntimeError:
            pprint.pprint(product)
            pprint.pprint(core_data)
            assert 0

    def remove_current_stock(self):
        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.csv"):
            os.remove(f"{pathlib.Path(__file__).parent.resolve()}/inventory.csv")

    def clean_html(self, text):
        clean_text = str(text).replace("\\n", "").replace("\\t", "").replace("<div>", "").replace("</div>", "")
        return clean_text.replace("<span>", "").replace("</span>", "")

    def get_core_name(self, shape, material, gapping, coating=None):
        if len(gapping) == 0:
            if coating is None:
                return f"{shape} - {material} - Ungapped"
            else:
                return f"{shape} - {coating} coated - {material} - Ungapped"
        else:
            number_gaps = len([gap['length'] for gap in gapping if gap['type'] == 'subtractive'])
            if number_gaps == 1:
                return f"{shape} - {material} - Gapped {gapping[0]['length'] * 1000:.3f} mm"
            else:
                return f"{shape} - {material} - Distributed gapped {gapping[0]['length'] * 1000:.3f} mm"

    def find_shape_closest_dimensions(self, family, shape, limit=0.05):
        smallest_distance = math.inf
        smallest_distance_shape = None
        core_data = PyMKF.get_available_core_shapes()

        dimensions = shape.split(" ")[1].split("/")

        print("dimensions")
        print(dimensions)
        try:
            for shape_name in core_data:
                if family == shape_name.split(" ")[0]:
                    aux_dimensions = shape_name.split(" ")[1].split("/")

                    distance = 0
                    for index in range(min(len(dimensions), len(aux_dimensions))):
                        with contextlib.suppress(ValueError):
                            distance += abs(float(dimensions[index]) - float(aux_dimensions[index])) / float(dimensions[index])
                    if smallest_distance > distance:
                        smallest_distance = distance
                        smallest_distance_shape = shape_name
        except IndexError:
            print("shape")
            print(shape)
            print("shape_name")
            print(shape_name)
            print("dimensions")
            print(dimensions)
            print("aux_dimensions")
            print(aux_dimensions)
            assert 0

        print(smallest_distance)
        print(smallest_distance_shape)
        if smallest_distance > limit:
            return None
            # assert 0, f"Unkown shape {shape}"
        return smallest_distance_shape

    def find_shape_closest_effective_parameters(self, family, effective_length, effective_area, effective_volume, limit=0.05):
        smallest_distance = math.inf
        smallest_distance_shape = None
        shape_names = PyMKF.get_available_core_shapes()
        dummyCore = {
            "functionalDescription": {
                "name": "dummy",
                "type": "two-piece set",
                "material": "dummy",
                "shape": None,
                "gapping": [],
                "numberStacks": 1
            }
        }
        try:
            for shape_name in shape_names:
                if family == shape_name.split(" ")[0]:
                    core = copy.deepcopy(dummyCore)
                    if family.lower() in ['t', 'r']:
                        core['functionalDescription']['type'] = "toroidal"
                    if family.lower() in ['ut']:
                        core['functionalDescription']['type'] = "closed shape"
                    core['functionalDescription']['shape'] = shape_name

                    core_data = PyMKF.get_core_data(core, False)

                    if abs(effective_length - core_data['processedDescription']['effectiveParameters']['effectiveLength']) / effective_length < 0.4:
                        if abs(effective_area - core_data['processedDescription']['effectiveParameters']['effectiveArea']) / effective_area < 0.4:
                            if abs(effective_volume - core_data['processedDescription']['effectiveParameters']['effectiveVolume']) / effective_volume < 0.4:

                                distance = abs(effective_length - core_data['processedDescription']['effectiveParameters']['effectiveLength']) / effective_length + abs(effective_area - core_data['processedDescription']['effectiveParameters']['effectiveArea']) / effective_area + abs(effective_volume - core_data['processedDescription']['effectiveParameters']['effectiveVolume']) / effective_volume

                                if smallest_distance > distance:
                                    smallest_distance = distance
                                    smallest_distance_shape = shape_name
        except IndexError:
            print("shape")
            print(shape)
            print("shape_name")
            print(shape_name)
            print("dimensions")
            print(dimensions)
            print("aux_dimensions")
            print(aux_dimensions)
            assert 0

        return smallest_distance_shape

    def add_core(self, family, shape, material, core_type, core_manufacturer, status, manufacturer_part_number, gapping=[], coating=None, datasheet=None, unique=False):
        core_name = self.get_core_name(shape, material, gapping, coating)

        if self.core_data[self.core_data['name'] == core_name].empty:
            core = {
                "name": core_name,
                "manufacturerInfo": {
                    "name": core_manufacturer, 
                    "status": status,
                    "reference": manufacturer_part_number,
                    "datasheetUrl": datasheet
                },
                "distributorsInfo": [],
                "functionalDescription": {
                    "type": core_type,
                    "material": material,
                    "shape": shape,
                    "gapping": gapping,
                    "numberStacks": 1
                }
            }
            if isinstance(shape, str):
                shape = PyMKF.get_core_data(core, False)['functionalDescription']['shape']
                for dimension_key, dimension in shape['dimensions'].items():
                    new_dimension = {}
                    for key, value in dimension.items():
                        if value != None:
                            new_dimension[key] = value
                    shape['dimensions'][dimension_key] = new_dimension
                core['functionalDescription']['shape'] = shape

            print(f"Adding {core_name}")
            core['manufacturer_reference'] = manufacturer_part_number
            if coating is not None:
                core['functionalDescription']['coating'] = coating

            core_series = pandas.Series()
            core_series["name"] = core["name"]
            core_series["manufacturerInfo"] = core["manufacturerInfo"]
            core_series["distributorsInfo"] = core["distributorsInfo"]
            core_series["functionalDescription"] = core["functionalDescription"]
            core_series["manufacturer_reference"] = core["manufacturer_reference"]

            self.core_data = pandas.concat([self.core_data, core_series.to_frame().T], ignore_index=True)
        elif unique:
            assert 0, f"Already inserted core: {core_name}"

        return core_name

    def add_distributor(self, core_name, core_manufacturer, distributor_reference, product_url, quantity, cost, distributor_name=None):
        print("core_name")
        print(core_name)
        row_index = self.core_data[self.core_data['name'] == core_name].index
        core = self.core_data[self.core_data['name'] == core_name].iloc[0]
        print("row_index")
        print(row_index)

        if isinstance(core['distributorsInfo'], str):
            core['distributorsInfo'] = ast.literal_eval(core['distributorsInfo'])

        print(f"Adding distributor for {core_name}")

        # Such is life...
        exceptions = [
            ("TDK", "ELP"),
            ("TDK", "P "),
            ("TDK", "PQ 32/30 - 95 - Distributed gapped 0.710 mm"),
            ("Magnetics", "T 22.1/13.7/7.9 - epoxy coated - MPP 60 - Ungapped"),
            ("TDK", "ER 11/5 - N87 - Ungapped"),
        ]

        for distributor_info in core['distributorsInfo']:
            if (distributor_info["name"] == distributor_name):
                if core_manufacturer == "Fair-Rite":
                    return
                if product_url == distributor_info["link"]:
                    return
                for exception in exceptions:
                    if core_manufacturer == exception[0] and exception[1] in core_name:
                        return
                print(distributor_info["link"])
                print(product_url)
                print("Distributor already in list")
                assert 0

        distributor_info = {
            "name": distributor_name,
            "link": product_url,
            "country": "UK" if distributor_name == "Gateway" else "USA",
            "distributedArea": "International",
            "reference": distributor_reference,
            "quantity": quantity,
            "updatedAt": date.today().strftime("%d/%m/%Y"),
            "cost": cost,
        }
        core['distributorsInfo'].append(distributor_info)

        core_series = pandas.Series()
        core_series["name"] = core["name"]
        core_series["manufacturerInfo"] = core["manufacturerInfo"]
        core_series["distributorsInfo"] = core["distributorsInfo"]
        core_series["functionalDescription"] = core["functionalDescription"]
        core_series["manufacturer_reference"] = core["manufacturer_reference"]

        self.core_data.iloc[row_index] = [core_series]

    def get_gapping(self, core_data, manufacturer_name, al_value, number_gaps=1):
        constants = PyMKF.get_constants()
        inductance = 1
        previous_inductance = 1
        previous_gap_length = constants["residualGap"]
        gap_length = constants["residualGap"]
        while inductance > al_value:
            previous_inductance = inductance
            previous_gap_length = gap_length

            winding_data = {
                "bobbin": "Dummy",
                "functionalDescription":[ 
                    {
                        "name": "MyWinding",
                        "numberTurns": 1,
                        "numberParallels": 1,
                        "isolationSide": "primary",
                        "wire": "Dummy"
                    }
                ]
            }
            operating_point = {
                'conditions': {'ambientTemperature': 25.0},
                'excitationsPerWinding': [
                    {
                        "frequency": 100000,
                        "current": {
                            "waveform": {
                                "data": [-5, 5, -5],
                                "time": [0, 0.0000025, 0.00001],
                            }
                        },
                        "voltage": {
                            "waveform": {
                                "data": [-2.5, 7.5, 7.5, -2.5, -2.5],
                                "time": [0, 0, 0.0000025, 0.0000025, 0.00001],
                            }
                        }
                    }
                ]
            }
            models = {'coreLosses': 'IGSE', 'coreTemperature': 'MANIKTALA', 'gapReluctance': 'ZHANG'}
            core_data['manufacturerInfo'] = {'name': manufacturer_name}
            if number_gaps == 1:
                core_data['functionalDescription']['gapping'] = [
                    {
                        "type": "subtractive",
                        "length": gap_length
                    },
                    {
                        "type": "residual",
                        "length": constants["residualGap"]
                    },
                    {
                        "type": "residual",
                        "length": constants["residualGap"]
                    }
                ]
            else:
                core_data['functionalDescription']['gapping'] = []
                for _ in range(number_gaps):
                    core_data['functionalDescription']['gapping'].append(
                        {
                        "type": "subtractive",
                        "length": gap_length
                        }
                    )
                for _ in range(len(core_data['processedDescription']['columns']) - 1):
                    core_data['functionalDescription']['gapping'].append(
                        {
                        "type": "residual",
                        "length": constants["residualGap"]
                        }
                    )
            try:
                inductance_data = PyMKF.calculate_inductance_from_number_turns_and_gapping(core_data, 
                                                                                           winding_data,
                                                                                           operating_point,
                                                                                           models)
                inductance = inductance_data['magnetizingInductance']['nominal']
            except ValueError:
                print(core_data)
                print(winding_data)
                print(operating_point)
                print(models)
                assert 0
            gap_length += 0.000001
            for index in range(len(core_data['functionalDescription']['gapping'])):
                core_data['functionalDescription']['gapping'][index]['length'] = round(core_data['functionalDescription']['gapping'][index]['length'], 5)

        if abs(previous_inductance - al_value) < abs(inductance - al_value):
            for index in range(len(core_data['functionalDescription']['gapping'])):
                if (core_data['functionalDescription']['gapping'][index]['type'] == 'subtractive'):
                    core_data['functionalDescription']['gapping'][index]['length'] = round(previous_gap_length, 5)

        return core_data['functionalDescription']['gapping']

    def process_gapping(self, core_data, gapping):
        core_data['functionalDescription']['gapping'] = gapping
        core_data = PyMKF.get_core_data(core_data, False)
        return core_data['functionalDescription']['gapping']

    def adaptive_round(self, value):
        if value < 4:
            return round(value, 2)
        elif value < 20:
            return round(value, 1)
        else:
            return round(value)

    def try_get_family(self, text):
        if match := re.match(
            r"([a-z]+)", text, re.I
        ):
            family = match.groups()[0]
            return family
        return None

    def try_get_shape(self, text):
        def fix(text):
            if len(text.split("-")) > 1:
                text = text.split("-")[0]
            return text.replace(" ", "").replace(",", ".").upper().replace("TX", "T").replace("TC", "T").replace("X", "/")

        current_error = 1
        current_shape = None
        for shape_name in PyMKF.get_available_core_shapes():
            error = abs(len(shape_name) - len(text)) / len(text)
            if shape_name in text and error < 0.3:
                print(error)
                if error < current_error:
                    current_error = error
                    current_shape = shape_name

        if current_shape is not None:
            print(current_shape)
            return current_shape

        fixed_text = fix(text)
        current_error = 1
        current_shape = None
        for shape_name in PyMKF.get_available_core_shapes():
            if (len(fixed_text) == 0):
                return None
            error = abs(len(fix(shape_name)) - len(fixed_text)) / len(fixed_text)
            if fix(shape_name) in fixed_text and error < 0.3:
                print(fix(shape_name))
                print(fixed_text)
                print(error)
                if error < current_error:
                    current_error = error
                    current_shape = shape_name
                # if text == "E55/28/25-3C92":
                #     print(abs(len(fix(shape_name)) - len(fixed_text)) / len(fixed_text))
                #     print(fix(shape_name))
                #     print(fixed_text)
                #     assert 0
        if current_shape is not None:
            print(current_shape)
            return current_shape

        print("fixed_text")
        print(fixed_text)
        return None

    def try_get_core_coating(self, text):
        if 'TC' in text:
            return "parylene"
        elif 'TX' in text:
            return "epoxy"
        else:
            return None

    def try_get_material(self, text, manufacturer):  # sourcery skip: use-next
        def fix(text):
            return text.replace(" 0", " ").replace(" ", "").replace(",", ".").replace("Mµ", "Mu").replace("Hƒ", "Hƒ").upper()

        materials = PyMKF.get_available_core_materials(manufacturer)
        materials.reverse()
        for material_name in materials:
            if material_name in text:
                return material_name

        fixed_text = fix(text)
        for material_name in materials:
            if fix(material_name) in fixed_text:
                return material_name
        return None

    def try_get_family_and_shape_statistically(self, text):
        def fix(text):
            return text.replace(" ", "").replace("X", "/").replace(",", ".").upper()

        # aux = text.split("Ferrite Cores & Accessories ")[1]
        fixed_text = fix(text)
        max_ratio = 0
        best_fit = ""

        for shape_name in PyMKF.get_available_core_shapes():
            s = SequenceMatcher(None, fix(shape_name), fixed_text)
            if s.ratio() > max_ratio:
                max_ratio = s.ratio()
                best_fit = shape_name

        print("max_ratio")
        print(max_ratio)
        if max_ratio < 0.9:
            return None, None

        family = best_fit.split(" ")[0]

        return family, best_fit


class DigikeyStocker(Stocker):

    def __init__(self):
        self.client_id = 'cN8i6L6KnNGJB2h3zsQgC7KvWf8AccsC'
        self.client_secret = '8QpIINW6VK9loIeF'
        self.redirect_uri = 'https://localhost:8139/digikey_callback'
        self.credentials_path = f"{pathlib.Path(__file__).parent.resolve()}/credentials.json"
        self.core_data = {}
        self.unfound_shapes = []
        self.ignored_manufacturers = ["VACUUMSCHMELZE",
                                      "Omron Automation and Safety",
                                      "Toshiba Semiconductor and Storage",
                                      "United Chemi-Con",
                                      "Würth Elektronik",
                                      "Laird-Signal Integrity Products",
                                      "Chemi-Con"]

    def get_refresh_token(self):
        # go to the browser and paste:
        # https://sandbox-api.digikey.com/v1/oauth2/authorize?response_type=code&client_id=Cpw48RBF8aYMTiTWeAYLvyUgPo2Db3F7&redirect_uri=https://localhost:8139/digikey_callback
        # https://api.digikey.com/v1/oauth2/authorize?response_type=code&client_id=cN8i6L6KnNGJB2h3zsQgC7KvWf8AccsC&redirect_uri=https://localhost:8139/digikey_callback
        # get auth code from the redirected URL and paste it here:

        auth_code = 'xTYZGMpO'

        token_response = requests.post('https://api.digikey.com/v1/oauth2/token', data={
            'client_id': self.client_id,
            'client_secret': self.client_secret,
            'redirect_uri': self.redirect_uri,
            'code': auth_code,
            'grant_type': 'authorization_code'
        })

        print(token_response.json())
        refresh_token = token_response.json()['refresh_token']

        with open(self.credentials_path, 'w') as outfile:
            json.dump({"refresh_token": refresh_token}, outfile)

    def get_access_token(self):

        with open(self.credentials_path) as json_file:
            data = json.load(json_file)

        refresh_token = data['refresh_token']

        token_response = requests.post('https://api.digikey.com/v1/oauth2/token', data={
            'client_id': self.client_id,
            'client_secret': self.client_secret,
            'refresh_token': refresh_token,
            'grant_type': 'refresh_token'
        })

        print(token_response.json())
        access_token = token_response.json()['access_token']
        refresh_token = token_response.json()['refresh_token']

        with open(self.credentials_path, 'w') as outfile:
            json.dump(token_response.json(), outfile)

        return access_token

    def read_products(self, access_token, offset: int = 0):  # sourcery skip: extract-method
        headers = {
            'accept': 'application/json',
            'Authorization': f"Bearer {access_token}",
            'X-DIGIKEY-Client-Id': self.client_id,
            'Content-Type': 'application/json',
        }

        json_data = {
            'Keywords': '',
            'RecordCount': 50,
            'RecordStartPosition': offset,
            'Filters': {
                'TaxonomyIds': [
                    936,
                ],
                'ManufacturerIds': [
                    # 0,
                ],
                'ParametricFilters': [
                    {
                        'ParameterId': 1989,
                        'ValueId': '0',
                    },
                    # {
                    #     'ParameterId': 1874,
                    #     'ValueId': '350255',
                    # },
                ],
            },
            'Sort': {
                'SortOption': 'SortByDigiKeyPartNumber',
                'Direction': 'Ascending',
                'SortParameterId': 0,
            },
            'RequestedQuantity': 0,
            'SearchOptions': [
                'ManufacturerPartSearch',
            ],
            'ExcludeMarketPlaceProducts': True,
        }

        url = 'https://api.digikey.com/Search/v3/Products/Keyword'

        response = requests.post(url, headers=headers, json=json_data)


        if 'ProductsCount' not in response.json():
            pprint.pprint(f"Daily limit reached, try again in {response.headers['X-RateLimit-Reset']} seconds")
            return {
                'current_offset' : -1,
                'total': -1
            }

        # pprint.pprint(response.json())
        # assert 0
        print(response.json()['ProductsCount'])

        if 'Products' not in response.json():
            print(response.json())

        for product in response.json()['Products']:
            pprint.pprint('=================================================')

            if product['ProductStatus'] == 'Obsolete':
                continue
            if product['QuantityAvailable'] == 0:
                continue
            if product['Manufacturer']['Value'] == 'Ferroxcube':
                # print("Ferroxcube")
                self.process_ferroxcube_product(product)
            elif product['Manufacturer']['Value'] == 'EPCOS - TDK Electronics':
                # continue
                self.process_epcos_tdk_product(product)
            elif product['Manufacturer']['Value'] == 'TDK Corporation':
                self.process_tdk_product(product)
            elif product['Manufacturer']['Value'] == 'Magnetics, a division of Spang & Co.':
                self.process_magnetics_product(product)
            elif product['Manufacturer']['Value'] == 'Fair-Rite Products Corp.':
                self.process_fair_rite_product(product)
            elif product['Manufacturer']['Value'] not in self.ignored_manufacturers:
                pprint.pprint(product['ProductDescription'])
                pprint.pprint(product['ManufacturerPartNumber'])
                pprint.pprint(product['NonStock'])
                pprint.pprint(product['QuantityAvailable'])
                pprint.pprint(product['UnitPrice'])
                pprint.pprint(product['DetailedDescription'])
                pprint.pprint(product['Manufacturer']['Value'])
                pprint.pprint(product['Series']['Value'])
                continue
                # assert 0, f"Unknown manufacturer: {product['Manufacturer']['Value']}"

        print("response.json()['ProductsCount']")
        print(response.json()['ProductsCount'])
        return {
            'current_offset' : offset + 50,
            'total': response.json()['ProductsCount']
        }

    def get_parameter(self, parameters, desired_parameter):
        for parameter in parameters:
            if parameter['Parameter'] == desired_parameter:
                return parameter['Value']

    def try_get_family(self, text):
        if match := re.match(
            r"([a-z]+)", text, re.I
        ):
            family = match.groups()[0]
            print(family)
            return family
        return None

    def process_epcos_tdk_product(self, product):  # sourcery skip: class-extract-method, extract-method, low-code-quality, move-assign
        wrong_products = [
                          "B66291P0000X149",
                          "B66293P0000X187",
                          "B66281P0000X187",
                          "B65877B0000R049",
                          "B64290A0618X087",
                          "B65803A0000R030",
                          "B65803J0040A087",
                          "B65803J0160A087",
                          "B65857C0000R087",
                          "B66281P0000X149",
                          "B66281P0000X197",
                          "B66283P0000X149",
                          "B66295P0000X149",
                          ]
        material_corrections = [("B65879Q0100K097", "N95", "N97"),
                                ("B65879Q0150K097", "N95", "N97")]
        al_value_corrections = [("B65879Q0150K097", "100", "150")]
        exceptions = ["B64290"]
        not_included_families = ["EV", "MHB", "PCH", "PS", "I", "ILP", "QU", "TT", "B", "PTS"]
        typos = [{"typo": "PM 47/28", "corrected_shape": "P 47/28", "corrected_family": "P"},
                 {"typo": "818", "corrected_shape": "RM8"},
                 {"typo": "0", "corrected_shape": "RM12"},
                 {"typo": "97212", "corrected_shape": "TT14"},
                 {"typo": "97217", "corrected_shape": "TT23/18"},
                 {"typo": "97461", "corrected_shape": "P30/19"},
                 {"typo": "P 23/13", "corrected_shape": "P 22/13", "corrected_family": "P"}]
        corrections = ["P 47/28"]
        corrections_family = ["P"]
        not_included_materials = ["H5C2", "PC44"]

        if product['ManufacturerPartNumber'] in wrong_products:
            return

        for parameter_aux in product['Parameters']:
            if parameter_aux['Parameter'] == 'Supplier Device Package':
                if parameter_aux['Value'] == '-' and product['Series']['Value'] == '-':
                   return

        core = self.core_data[self.core_data["manufacturer_reference"].str.contains(product['ManufacturerPartNumber'])]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:
            shape_type = "two-piece set"
            material = None
            gap_length = None
            shape = None
            coating = None
            for parameter in product['Parameters']:
                if parameter['Parameter'] == 'Core Type' and parameter['Value'] == 'Toroid':
                    shape_type = 'toroidal'
                if parameter['Parameter'] == 'Material':
                    material = self.try_get_material(parameter['Value'], 'TDK')

                    for aux in material_corrections:
                        if product['ManufacturerPartNumber'] == aux[0] and material == aux[1]:
                            material = aux[2]
            
                    if material in not_included_materials:
                        return None

            # Some cores have typos in the series value, we correct them with the Id
            for typo in typos:
                if product['Series']['ValueId'] == typo['typo']:
                    product['Series']['Value'] = typo['corrected_shape']

            family = self.try_get_family(product['Series']['Value'])

            if material is None:
                return

            for parameter_aux in product['Parameters']:
                if parameter_aux['Parameter'] == 'Supplier Device Package':
                    if parameter_aux['Value'].startswith("I "):
                        return
                    print("parameter_aux['Value']")
                    print(parameter_aux['Value'])

            if family is not None:

                if family in not_included_families:
                    return

                for parameter_aux in product['Parameters']:
                    if parameter_aux['Parameter'] == 'Supplier Device Package':
                        shape = self.try_get_shape(parameter_aux['Value'])

                if shape is None:
                    try: 
                        shape = product['Series']['Value'].replace(family, f"{family} ")
                        print("shape")
                        print(shape)
                    except ValueError:
                        print(family)
                        print(product['Series']['Value'])
                        print(product['Series']['Value'].replace(family, f"{family} "))
                        assert 0

                    if family in ['T', 'TX', 'TC']:
                        family = 'T'
                        shape = shape.replace('TX', 'T').replace('TC', 'T')

                    if family not in PyMKF.get_available_shape_families() and product['Series']['Value'] not in exceptions:
                        for parameter_aux in product['Parameters']:
                            if parameter_aux['Parameter'] == 'Supplier Device Package':
                                shape = parameter_aux['Value'].replace(" x ", "/")
                                family = shape.split(" ")[0]

                                if family in ['T', 'TX', 'TC']:
                                    family = 'T'
                                    shape = shape.replace('TX', 'T').replace('TC', 'T')

                                if family in not_included_families:
                                    return
                                if family not in PyMKF.get_available_shape_families():
                                    pprint.pprint(product)
                                    print(family)
                                    print(shape)
                                    print("get_available_shape_families")
                                    print(PyMKF.get_available_shape_families())
                                    assert 0

                    # Some cores have typos in the series value, we correct them with the name
                    for typo in typos:
                        if shape == typo['typo']:
                            shape = typo['corrected_shape']
                            family = typo['corrected_family']

                if family in ['T', 'TX', 'TC']:
                    family = 'T'
                    shape = shape.replace('R', 'T')

                # We check that we have the core in the database
                core_data = self.verify_shape_is_in_database(family, shape, material, product)
                coating = self.try_get_core_coating(product['ManufacturerPartNumber'])
                if coating is None:
                    coating = self.try_get_core_coating(product['Series']['Value'])
                if coating is None:
                    for parameter_aux in product['Parameters']:
                        if parameter_aux['Parameter'] == 'Finish':
                            if parameter_aux['Value'] == 'Epoxy':
                                coating = 'epoxy'
                            if parameter_aux['Value'] == 'Parylene':
                                coating = 'parylene'

                core_type = 'toroidal' if family == 'T' else 'two-piece set'
                manufacturer_part_number = product['ManufacturerPartNumber']
                core_manufacturer = 'TDK'
                gapping = []

                status = 'production'
                for parameter in product['Parameters']:
                    if parameter['Parameter'] == 'Gap' and parameter['Value'] not in ['-', 'Ungapped']:
                        number_gaps = 3 if parameter['Value'] == 'Distributed Gapped' else 1
                        al_value = None
                        for parameter_aux in product['Parameters']:
                            if parameter_aux['Parameter'] == 'Inductance Factor (Al)':
                                al_value = parameter_aux['Value'].split(" ")[0]
                                al_value_unit = parameter_aux['Value'].split(" ")[1]

                                for aux in al_value_corrections:
                                    if product['ManufacturerPartNumber'] == aux[0] and al_value == aux[1]:
                                        al_value = aux[2]

                                if al_value_unit == 'nH':
                                    al_value = float(al_value) * 1e-9
                                elif al_value_unit == 'µH':
                                    al_value = float(al_value) * 1e-6
                                else:
                                    print(al_value_unit)
                                    assert 0
                        if family in ["TC", "TX", "T"]:
                            return
                        if al_value is None:
                            return

                        gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)

                    if parameter['Parameter'] == 'Part Status' and parameter['Value'] != 'Active':
                        assert 0, f"Part Status is {parameter['Value']}"

                if family in ['T', 'TX', 'TC']:
                    family = 'T'
                    shape = shape.replace('TX', 'T')
                    shape = shape.replace('TC', 'T')

                core_name = self.add_core(family=family,
                                          shape=shape,
                                          material=material,
                                          core_type=core_type,
                                          core_manufacturer=core_manufacturer,
                                          status=status,
                                          manufacturer_part_number=manufacturer_part_number,
                                          coating=coating,
                                          gapping=gapping)

            else:
                pprint.pprint(product)
                assert 0

        core_manufacturer = 'TDK'
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Digi-Key")

    def process_fair_rite_product(self, product):  # sourcery skip: class-extract-method, extract-method, low-code-quality, move-assign
        exceptions = []
        not_included_families = ['Rod', 'ROD', 'Planar E/I', 'PS', 'Multi-Hole (2)', 'Multi-Hole', 'I', 'Tube']
        typos = []
        shape_type = "two-piece set"
        family = None
        material = None
        gap_length = None
        effective_length = None
        effective_area = None
        effective_volume = None
        coating = None

        core = self.core_data[self.core_data["manufacturer_reference"] == product['ManufacturerPartNumber']]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:

            for parameter in product['Parameters']:
                if parameter['Parameter'] == 'Core Type':
                    if parameter['Value'] in not_included_families:
                        return
                    elif parameter['Value'] == 'Toroid':
                        shape_type = 'toroidal'
                        family = 'T'
                    elif parameter['Value'] == 'EER':
                        family = 'ER'
                    elif parameter['Value'] == 'P (Pot Core)':
                        family = 'P'
                    elif parameter['Value'] in PyMKF.get_available_shape_families():
                        family = parameter['Value']
                    else:
                        print(parameter['Value'])
                        assert 0
                if parameter['Parameter'] == 'Material':
                    material = self.try_get_material(parameter['Value'], 'Fair-Rite')
                if parameter['Parameter'] == 'Effective Length (le) mm':
                    effective_length = float(parameter['Value']) / 1000
                if parameter['Parameter'] == 'Effective Area (Ae) mm²':
                    effective_area = float(parameter['Value']) / 1000000
                if parameter['Parameter'] == 'Effective Magnetic Volume (Ve) mm³':
                    effective_volume = float(parameter['Value']) / 1000000000


            if family is None:
                return
            if effective_length is None:
                assert 0
            if effective_area is None:
                assert 0
            if effective_volume is None:
                assert 0
            if material is None:
                return

            shape = self.find_shape_closest_effective_parameters(family=family,
                                                                 effective_length=effective_length,
                                                                 effective_area=effective_area,
                                                                 effective_volume=effective_volume)

            if shape is None:
                return

            # We check that we have the core in the database
            core_data = self.verify_shape_is_in_database(family, shape, material, product)

            core_type = 'toroidal' if family == 'T' else 'two-piece set'
            core_manufacturer = 'Fair-Rite'
            manufacturer_part_number = product['ManufacturerPartNumber']
            product_url = product['ProductUrl']
            distributor_reference = product['DigiKeyPartNumber']
            quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
            cost = float(product['UnitPrice'])
            gapping = []

            status = 'production'

            core_name = self.add_core(family=family,
                                      shape=shape,
                                      material=material,
                                      core_type=core_type,
                                      core_manufacturer=core_manufacturer,
                                      status=status,
                                      manufacturer_part_number=manufacturer_part_number,
                                      coating=coating,
                                      gapping=gapping)

        core_manufacturer = 'Fair-Rite'
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Digi-Key")

    def process_tdk_product(self, product):  # sourcery skip: class-extract-method, extract-method, low-code-quality, move-assign
        not_included_families = ["EPC", "EER", "EIR", "ELT", "PQI", "EI", "SP", "-"]
        typos = [
                 {"typo": "EI 13/4.4", "corrected_shape": "EL 13/4.4", "corrected_family": "EL"},
                 {"typo": "EI 15.5/5.8", "corrected_shape": "EL 15.5/5.8", "corrected_family": "EL"},
                 {"typo": "EI 14/4.5/9", "corrected_shape": "ER 14/4.5/9", "corrected_family": "ER"},
                 # {"typo": "EIR 14/4.5/9", "corrected_shape": "ER 14/4.5/9", "corrected_family": "ER"},
                 # {"typo": "EI 18/5/12", "corrected_shape": "ER 18/5/12", "corrected_family": "ER"},
                 # {"typo": "EIR 18/5/12", "corrected_shape": "ER 18/5/12", "corrected_family": "ER"},
                 {"typo": "EI 22/5.5/15", "corrected_shape": "ER 22/5.5/15", "corrected_family": "ER"},
                 # {"typo": "EIR 22/5.5/15", "corrected_shape": "ER 22/5.5/15", "corrected_family": "ER"}
                 ]
        corrections = ["P 47/28"]
        corrections_family = ["P"]
        not_included_materials = ["H5C2", "PC44", "PC40", "H5A"]

        core = self.core_data[self.core_data["manufacturer_reference"].str.contains(product['ManufacturerPartNumber'])]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:
            shape_type = "two-piece set"
            material = None
            gap_length = None
            for parameter in product['Parameters']:
                if parameter['Parameter'] == 'Core Type' and parameter['Value'] == 'Toroid':
                    shape_type = 'toroidal'
                if parameter['Parameter'] == 'Material':
                    material = parameter['Value']
                    if material in not_included_materials:
                        return None
            family = None
            coating = None

            for parameter_aux in product['Parameters']:
                if parameter_aux['Parameter'] == 'Supplier Device Package':
                    shape = parameter_aux['Value'].replace(" x ", "/")
                    family = shape.split(" ")[0]

                    # Some cores have typos in the series value, we correct them with the name
                    for typo in typos:
                        if shape == typo['typo']:
                            shape = typo['corrected_shape']
                            family = typo['corrected_family']

                    if family in ['T', 'TX', 'TC']:
                        family = 'T'
                        shape = shape.replace('TX', 'T').replace('TC', 'T')

                    if family in not_included_families:
                        return
                    if family not in PyMKF.get_available_shape_families():
                        pprint.pprint(product)
                        print(family)
                        print(shape)
                        print("get_available_shape_families")
                        print(PyMKF.get_available_shape_families())
                        assert 0

            if material is None:
                return

            if family is not None:

                if family in not_included_families:
                    return

                # Some cores have typos in the series value, we correct them with the name
                for typo in typos:
                    if shape == typo['typo']:
                        shape = typo['corrected_shape']
                        family = typo['corrected_family']


                # We check that we have the core in the database
                core_data = self.verify_shape_is_in_database(family, shape, material, product)

                status = 'production'
                core_type = 'toroidal' if family == 'T' else 'two-piece set'
                core_manufacturer = 'TDK'
                manufacturer_part_number = product['ManufacturerPartNumber']

                gapping = []

                for parameter in product['Parameters']:
                    if parameter['Parameter'] == 'Gap' and parameter['Value'] not in ['-', 'Ungapped']:
                        number_gaps = 3 if parameter['Value'] == 'Distributed Gapped' else 1
                        if family in ["TC", "TX", "T"]:
                            return

                        al_value = None
                        for parameter in product['Parameters']:
                            if parameter['Parameter'] == 'Inductance Factor (Al)':
                                al_value = float(parameter['Value'].split(" nH")[0].strip()) * 1e-9
                        if al_value is None:
                            pprint.pprint(product)
                            assert 0

                        gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)
                    if parameter['Parameter'] == 'Part Status' and parameter['Value'] != 'Active':
                        assert 0, f"Part Status is {parameter['Value']}"

                if family in ['T', 'TX', 'TC']:
                    family = 'T'
                    shape = shape.replace('TX', 'T')
                    shape = shape.replace('TC', 'T')

                core_name = self.add_core(family=family,
                                          shape=shape,
                                          material=material,
                                          core_type=core_type,
                                          core_manufacturer=core_manufacturer,
                                          status=status,
                                          manufacturer_part_number=manufacturer_part_number,
                                          coating=coating,
                                          gapping=gapping)

            else:
                return

        core_manufacturer = 'TDK'
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Digi-Key")

    def process_ferroxcube_product(self, product):


        # sourcery skip: extract-method, low-code-quality, move-assign
        not_included_families = ["I", "PLT", "ROD"]
        typos = [{"typo": "ER 35W/21/11", "corrected_shape": "ER 35/21/11", "corrected_family": "ER"}]
        not_included_materials = ["3D", "3E", "3C11", "4C", "3C81", "3H1", "3E10"]

        # Removing plates for now
        if product['ManufacturerPartNumber'].endswith("-P"):
            return

        core = self.core_data[self.core_data["manufacturer_reference"] == product['ManufacturerPartNumber']]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:

            family = None
            gap_length = None
            if match := re.match(
                r"([a-z]+)", product['ManufacturerPartNumber'], re.I
            ):
                family = match.groups()[0]
                print(family)
            if family is not None:
                if family in not_included_families:
                    return
                try:
                    for not_included_material in not_included_materials:
                        if not_included_material in product['ManufacturerPartNumber']:
                            return None
                    material = self.try_get_material(product['ManufacturerPartNumber'], 'Ferroxcube')
                    shape = self.try_get_shape(product['ManufacturerPartNumber'])

                    if material is None:
                        return None
                        # assert 0, f"Unknown material in {product['ManufacturerPartNumber']}"

                    if shape is None:
                        print(product['ManufacturerPartNumber'])
                        family, shape = self.try_get_family_and_shape_statistically(product['ManufacturerPartNumber'])
                        if shape is None:
                            return

                except ValueError:
                    print(family)
                    print(product['ManufacturerPartNumber'])
                    print(product['ManufacturerPartNumber'].replace(family, f"{family} "))
                    print(product['ManufacturerPartNumber'].replace(family, f"{family} ").split("-"))
                    assert 0

                if family in ['T', 'TX', 'TC']:
                    family = 'T'
                    shape = shape.replace('TX', 'T').replace('TC', 'T')

                if family in not_included_families:
                    return

                # Some cores have typos in the series value, we correct them with the name
                for typo in typos:
                    if shape == typo['typo']:
                        shape = typo['corrected_shape']
                        family = typo['corrected_family']

                core_data = self.verify_shape_is_in_database(family, shape, material, product)

                coating = self.try_get_core_coating(product['ManufacturerPartNumber'])

                status = 'production'
                core_manufacturer = 'Ferroxcube'
                core_type = 'toroidal' if family.lower() in ['t', 'tx', 'tc'] else 'two-piece set'
                manufacturer_part_number = product['ManufacturerPartNumber']
                gapping = []

                try:
                    for parameter in product['Parameters']:
                        if parameter['Parameter'] == 'Gap' and parameter['Value'] not in ['Ungapped']:
                            number_gaps = 3 if parameter['Value'] == 'Distributed Gapped' else 1
                            try:
                                al_value = float(product['ManufacturerPartNumber'].split("-A")[1].split("-")[0]) * 1e-9
                            except IndexError:
                                al_value = float(product['ManufacturerPartNumber'].split("-E")[1].split("-")[0]) * 1e-9

                            gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)

                        if parameter['Parameter'] == 'Part Status' and parameter['Value'] != 'Active':
                            assert 0, f"Part Status is {parameter['Value']}"
                except IndexError:
                    pass

                core_name = self.add_core(family=family,
                                          shape=shape,
                                          material=material,
                                          core_type=core_type,
                                          core_manufacturer=core_manufacturer,
                                          status=status,
                                          manufacturer_part_number=manufacturer_part_number,
                                          coating=coating,
                                          gapping=gapping)

        core_manufacturer = 'Ferroxcube'
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2PC' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Digi-Key")

    def process_magnetics_product(self, product):
        not_included_materials = ["75-Series 125",
                                  "AmoFlux",
                                  "High Flux 19",
                                  "High Flux 90",
                                  "Kool Mu 19",
                                  "MPP 250",
                                  "MPP 367",
                                  "XFlux 14",
                                  "XFlux 50",
                                  "Nanocrystalline"]

        core = self.core_data[self.core_data["manufacturer_reference"] == product['ManufacturerPartNumber']]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:
            # pprint.pprint(product)
            try:
                material = self.try_get_material(self.get_parameter(product['Parameters'], 'Material') + " " + self.get_parameter(product['Parameters'], 'Initial Permeability (µi)'), 'Magnetics')
            except TypeError:
                material = None

            coating = None
            for parameter_aux in product['Parameters']:
                if parameter_aux['Parameter'] == 'Finish':
                    if parameter_aux['Value'] == 'Coated':
                        coating = 'epoxy'

            if self.get_parameter(product['Parameters'], 'Core Type') is None:
                return

            if self.get_parameter(product['Parameters'], 'Core Type') == 'Toroid':
                try:
                    diameter = self.get_parameter(product['Parameters'], 'Diameter')
                    external_diameter = diameter.split('(')[1].split('mm)')[0]
                except IndexError:
                    return
                height = self.get_parameter(product['Parameters'], 'Height').split('(')[1].split('mm)')[0]
                family = 'T'
                pprint.pprint("external_diameter")
                pprint.pprint(external_diameter)
                pprint.pprint("height")
                pprint.pprint(height)
                pprint.pprint("coating")
                pprint.pprint(coating)
                if coating is None:
                    temptative_shape = f"{family} {self.adaptive_round(float(external_diameter))}/{None}/{self.adaptive_round(float(height))}"
                else:
                    temptative_shape = f"{family} {self.adaptive_round(float(external_diameter) - 0.7)}/{None}/{self.adaptive_round(float(height) - 0.7)}"
                pprint.pprint(temptative_shape)
            elif self.get_parameter(product['Parameters'], 'Core Type') in ['I', 'P']:
                return
            else:
                try:
                    family = self.get_parameter(product['Parameters'], 'Core Type')
                    length = self.get_parameter(product['Parameters'], 'Length').split('(')[1].split('mm)')[0]
                    height = self.get_parameter(product['Parameters'], 'Height').split('(')[1].split('mm)')[0]
                    width = self.get_parameter(product['Parameters'], 'Width').split('(')[1].split('mm)')[0]
                    temptative_shape = f"{family} {self.adaptive_round(float(length))}/{self.adaptive_round(float(height))}/{self.adaptive_round(float(width))}"
                except IndexError:
                    # pprint.pprint(product)
                    # assert 0
                    return
            shape = self.find_shape_closest_dimensions(family, temptative_shape, limit=0.1)
            pprint.pprint(family)
            pprint.pprint(material)
            pprint.pprint(shape)

            if material is None:
                material = self.try_get_material(product['ProductDescription'], 'Magnetics')

            if material is None:
                for not_included_material in not_included_materials:
                    if not_included_material in f"{self.get_parameter(product['Parameters'], 'Material')} {self.get_parameter(product['Parameters'], 'Initial Permeability (µi)')}":
                        return None
                # pprint.pprint(product)
                return
                # assert 0, f"Unknown material in {self.get_parameter(product['Parameters'], 'Material')} {self.get_parameter(product['Parameters'], 'Initial Permeability (µi)')}"

            if shape is None:
                return
            core_data = self.verify_shape_is_in_database(family, shape, material, product)

            status = 'production'
            core_manufacturer = 'Magnetics'
            core_type = 'toroidal' if family.lower() in ['t', 'tx', 'tc'] else 'two-piece set'
            manufacturer_part_number = product['ManufacturerPartNumber']
            gapping = []

            core_name = self.add_core(family=family,
                                      shape=shape,
                                      material=material,
                                      core_type=core_type,
                                      core_manufacturer=core_manufacturer,
                                      status=status,
                                      manufacturer_part_number=manufacturer_part_number,
                                      coating=coating,
                                      gapping=gapping)

        core_manufacturer = 'Magnetics'
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Digi-Key")
        # if self.get_parameter(product['Parameters'], 'Core Type') == 'Toroid':
        #     approximated_shape = f"T {}/{}/{}"
        # assert 0

    def get_cores_stock(self):
        # access_token = self.get_refresh_token()
        access_token = self.get_access_token()

        remaining = 1
        current_offset = 0

        if not os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") and os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            shutil.copy(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson")

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        self.core_data = pandas.DataFrame(self.core_data.values())
        self.core_data["manufacturer_reference"] = self.core_data.apply(lambda row: row['manufacturerInfo']['reference'].split(' ')[0], axis=1)

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(access_token, current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            # pprint.pprint(f"self.core_data: {len(self.core_data.values())}")
            # break

        pprint.pprint(self.unfound_shapes)
        # core_data = pandas.DataFrame(self.core_data.values())
        self.core_data = self.core_data.drop(['manufacturer_reference'], axis=1)
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson", "w")
        ndjson.dump(self.core_data.to_dict('records'), out_file)
        out_file.close()


class MouserStocker(Stocker):

    def __init__(self):
        self.apikey = "6c28da81-2a3c-43b4-b03e-fca6847d5370"
        self.core_data = {}
        self.unfound_descriptions = []
        self.unfound_shapes = []
        self.ignored_manufacturers = [
                                        "Proterial",
                                        "Laird Performance Materials",
                                        "Wurth Elektronik",
                                        "United Chemi-Con",
                                        "Vacuumschmelze",
                                        "Chemi-Con",
                                     ]

    def try_get_family(self, text):
        def fix(text):
            return text.replace(" ", "").replace("CORE", "").replace("PLANAR", "").upper()

        aux_split = text.split("Ferrite Cores & Accessories ")
        if len(aux_split) < 2:
            return None

        aux = aux_split[1]
        if not (match := re.match(r"([a-z]+)", aux, re.I)):
            fixed_aux = fix(aux)
            if not (match := re.match(r"([a-z]+)", fixed_aux, re.I)):
                return None
            else:
                family = match.groups()[0]
                return family
            return None

        family = match.groups()[0]

        if family == 'R':
            family = None
        return family

    def process_fair_rite_product(self, product):
        not_included_families = ['ROD', 'Planar E/I', 'PS', 'Multi-Hole (2)', 'Multi-Hole', 'I', 'Tube', 'Bobbin', 'Plate', 'Hold']

        for family in not_included_families:
            if family.lower() in product['Description'].lower().split(' '):
                return

        core = self.core_data[self.core_data["manufacturer_reference"] == product['ManufacturerPartNumber']]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:
            pprint.pprint(product['ProductDetailUrl'])
            import requests

            url = product['ProductDetailUrl']
            headers = {
                'authority': 'www.dickssportinggoods.com',
                'pragma': 'no-cache',
                'cache-control': 'no-cache',
                'sec-ch-ua': '" Not;A Brand";v="99", "Google Chrome";v="91", "Chromium";v="91"',
                'sec-ch-ua-mobile': '?0',
                'upgrade-insecure-requests': '1',
                'user-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.114 Safari/537.36',
                'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
                'sec-fetch-site': 'none',
                'sec-fetch-mode': 'navigate',
                'sec-fetch-user': '?1',
                'sec-fetch-dest': 'document',
                'accept-language': 'es-ES,es;q=0.9,en;q=0.8"',
            }

            try:
                material = self.try_get_material(product['Description'].split('Ferrite Cores & Accessories')[1].split(' ')[1], 'Fair-Rite')
            except IndexError:
                return

            if product['Description'] is None:
                pprint.pprint(product['Description'])
                assert 0

            if material is None:
                return 

            family = self.try_get_family(product['Description'].replace(material, '').replace("Planar EE", "E").split("Core ")[0])
            pprint.pprint(material)
            pprint.pprint(family)

            if family is None:
                assert 0

            length = None
            width = None
            height = None
            shape = None

            for shape_name in PyMKF.get_available_core_shapes():
                if shape_name in product['Description']:
                    shape = shape_name
                    family = shape_name.split(' ')[0]
                    break
                if shape_name.upper().replace(' ', '') in product['Description'].upper().replace(' ', ''):
                    shape = shape_name
                    family = shape_name.split(' ')[0]
                    break

            if shape is None:
                response = requests.get(url, headers=headers)
                try:
                    extra_data = pandas.read_html(StringIO(response.text))[8]
                except ValueError:
                    return
                for datum_index, datum in extra_data.iterrows():
                    if datum['Atributo del producto'] == 'Longitud exterior:':
                        length = float(datum['Valor del atributo'].split(' mm')[0])
                        pprint.pprint(length)
                    if datum['Atributo del producto'] == 'Anchura exterior:':
                        width = float(datum['Valor del atributo'].split(' mm')[0])
                        pprint.pprint(width)
                    if datum['Atributo del producto'] == 'Altura exterior:':
                        height = float(datum['Valor del atributo'].split(' mm')[0])
                        pprint.pprint(height)

                print(length)
                print(width)
                print(height)

                if (length is not None and width is not None and height is not None):

                    temptative_shape = f"{family} {self.adaptive_round(float(length))}/{self.adaptive_round(float(height))}/{self.adaptive_round(float(width))}"

                    shape = self.find_shape_closest_dimensions(family, temptative_shape, limit=0.1)

                if shape is None:
                    return

            print(family)
            print(shape)
            print(material)

            core_data = self.verify_shape_is_in_database(family, shape, material, product['Description'])
            gapping = []
            coating = None

            core_type = 'toroidal' if family == 'T' else 'two-piece set'
            core_manufacturer = 'Fair-Rite'
            status = 'production'
            manufacturer_part_number = product['ManufacturerPartNumber']

            core_name = self.add_core(family=family,
                                      shape=shape,
                                      material=material,
                                      core_type=core_type,
                                      core_manufacturer=core_manufacturer,
                                      status=status,
                                      manufacturer_part_number=manufacturer_part_number,
                                      coating=coating,
                                      gapping=gapping)

        core_manufacturer = 'Fair-Rite'
        product_url = product['ProductDetailUrl']
        distributor_reference = product['MouserPartNumber']
        quantity = int(product['AvailabilityInStock'])
        for price_break in product['PriceBreaks']:
            if price_break['Quantity'] == 1:
                cost = float(price_break['Price'].split("$")[1])
                break
        else:
            assert 0
        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Mouser")

    def process_epcos_tdk_product(self, product):
        not_included_families = ["EV", "MHB", "PCH", "PS", "I", "ILP", "QU", "TT", "B", "PTS"]
        not_included_materials = ["N48", "N22", "T38", "M33", "T65", "T46", "T37", "K10", "K1"]
        typos = [
                    {"typo": "TOROID", "correction": "T"},
                    {"typo": "EER", "correction": "ER"},
                ]
        wrong_products = [
                          "B66417U160K187",
                          ]

        if product['ManufacturerPartNumber'] in wrong_products:
            return

        # print(product['ManufacturerPartNumber'])
        # print(self.core_data["manufacturer_reference"])
        core = self.core_data[self.core_data["manufacturer_reference"].str.contains(product['ManufacturerPartNumber'])]
        if len(core.index) != 1:
            self.unfound_shapes.append(product['ManufacturerPartNumber'])

        if not core.empty:
            core_name = core.iloc[0]['name']
        else:
            def get_core_data(description):
                family = self.try_get_family(description)

                material = self.try_get_material(description, 'TDK')

                if material is None:
                    for not_included_material in not_included_materials:
                        if not_included_material in description:
                            return None
                    # pprint.pprint(product)
                    self.unfound_descriptions.append(description)
                    return None
                    # assert 0, f"Unknown material: {material}"

                print(material)

                cleaned_description = description.split(material)[0]
                cleaned_description = cleaned_description.split('Ferrite Cores & Accessories')[1]
                cleaned_description = cleaned_description.replace('_', '')
                cleaned_description = cleaned_description.strip()
                if cleaned_description.endswith('-'):
                    cleaned_description = cleaned_description[0:-1]

                print(cleaned_description)
                print(family)
                if family in ('TX', 'TC'):
                    family = 'T'
                if family is None:
                    if description.endswith(' T'):
                        family = 'T'
                        cleaned_description = f"{family} {cleaned_description}"
                    elif cleaned_description.startswith('R') and not cleaned_description.startswith('RM'):
                        family = 'T'
                        cleaned_description = cleaned_description.replace('R', 'T')
                    elif cleaned_description.startswith('T'):
                        family = 'T'
                        cleaned_description = cleaned_description.replace('TX', 'T ')
                        cleaned_description = cleaned_description.replace('TC', 'T ')

                shape = None
                if family is None:
                    family, shape = self.try_get_family_and_shape_statistically(cleaned_description)
 
                if family is None:
                    pprint.pprint(product)
                    self.unfound_descriptions.append(cleaned_description)
                    return None
                    # assert 0, "Unknown family"

                if family not in PyMKF.get_available_shape_families() and family not in not_included_families:
                    pprint.pprint(product)
                    self.unfound_descriptions.append(cleaned_description)
                    return None
                    # assert 0, "Unknown family"

                if shape is None:
                    shape = self.try_get_shape(cleaned_description)

                if shape is None:
                    pprint.pprint(product)
                    self.unfound_descriptions.append(cleaned_description)
                    return None
                    # assert 0, "Unknown shape"

                # We check that we have the core in the database
                core_data = self.verify_shape_is_in_database(family, shape, material, product)

                if description.lower().endswith('mm'):
                    try:
                        gap_length = description[-6:-2]
                        gap_length = ast.literal_eval(gap_length.replace(',', '.')) / 1000
                    except SyntaxError:
                        gap_length = description[-5:-2]
                        gap_length = ast.literal_eval(gap_length.replace(',', '.')) / 1000

                    gapping = [
                        {
                            "type": "subtractive",
                            "length": gap_length
                        }]

                    gapping = self.process_gapping(core_data, gapping)

                else:
                    try:


                        if "-DG" in description:
                            number_gaps = 3
                            al_value = float(description.split("-DG")[1].split(" ")[0].strip()) * 1e-9
                        else:
                            number_gaps = 1
                            al_value = float(description.split(material)[1].strip().split(" ")[0]) * 1e-9
                        if al_value > 0:
                            gapping = self.get_gapping(core_data, "TDK", al_value, number_gaps)
                        else:
                            gapping = []
                    except (ValueError, IndexError):
                        try:
                            if product['ManufacturerPartNumber'][7:11] == '0000':
                                gapping = []
                            else:
                                al_value = float(product['ManufacturerPartNumber'][7:11].lstrip('0'))
                                gapping = self.get_gapping(core_data, "TDK", al_value, number_gaps)
                        except (ValueError, IndexError):
                            gapping = []
                    except (RuntimeError):
                        return None

                return family, shape, material, gapping

            exceptions = ["COIL FORMER",
                          "Coil former",
                          "BOBBIN",
                          "PS",
                          "EIR",
                          "YOKE",
                          "MOUNTING",
                          "WAS",
                          "Accessories 7.35X3.6",
                          "Ferrite Core N95 100x100x6",
                          "Accessories 70X14.5",
                          "U-CORE N97U CORE",
                          "DOUBLE-APERT",
                          "Clamp",
                          "COV-",
                          "NOT AVAILABLE",
                          "Accessories Qu",
                          "Accessories I",
                          "Accessories U-CORE N",
                          "Accessories CF-E",
                          "SUMIKON",
                          "P CORE HALF",
                          "CLI",
                          "CLM",
                          "HOUSING",
                          "screw",
                          "ACCS CF-",
                          ]
            for exception in exceptions:
                if exception.upper() in product['Description'].upper():
                    print("exception")
                    print(exception)
                    return

            for typo in typos:
                if typo['typo'] in product['Description']:
                    product['Description'] = product['Description'].replace(typo['typo'], typo['correction'])

            # pprint.pprint(product)
            core_data = get_core_data(product['Description'])
            if core_data is None:
                pprint.pprint(product['Description'])
                return

            family, shape, material, gapping = core_data
            if family != 'T':
                coating = None
            else:
                if product['ManufacturerPartNumber'][6] == 'A':
                    coating = None
                elif product['ManufacturerPartNumber'][6] == 'L':
                    coating = 'epoxy'
                elif product['ManufacturerPartNumber'][6] == 'P':
                    coating = 'parylene'
                else:
                    assert 0, product['ManufacturerPartNumber']


            core_type = 'toroidal' if family == 'T' else 'two-piece set'
            core_manufacturer = 'TDK'
            status = 'production'
            manufacturer_part_number = product['ManufacturerPartNumber']

            core_name = self.add_core(family=family,
                                      shape=shape,
                                      material=material,
                                      core_type=core_type,
                                      core_manufacturer=core_manufacturer,
                                      status=status,
                                      manufacturer_part_number=manufacturer_part_number,
                                      coating=coating,
                                      gapping=gapping)

        core_manufacturer = 'TDK'
        product_url = product['ProductDetailUrl']
        distributor_reference = product['MouserPartNumber']
        quantity = int(product['AvailabilityInStock'])
        for price_break in product['PriceBreaks']:
            if price_break['Quantity'] == 1:
                cost = float(price_break['Price'].split("$")[1])
                break
        else:
            assert 0

        self.add_distributor(core_name=core_name,
                             product_url=product_url,
                             core_manufacturer=core_manufacturer,
                             distributor_reference=distributor_reference,
                             quantity=quantity,
                             cost=cost,
                             distributor_name="Mouser")

        # pprint.pprint(self.core_data)

    def read_products(self, offset: int = 0):  # sourcery skip: extract-method
        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/json',
        }

        params = {
            'apiKey': '6c28da81-2a3c-43b4-b03e-fca6847d5370',
        }

        json_data = {
            'SearchByKeywordRequest': {
                'keyword': 'Ferrite Cores & Accessories',
                'records': 50,
                'startingRecord': offset,
                'searchOptions': '4',
                'searchWithYourSignUpLanguage': 'false',
            },
        }

        response = requests.post(f'https://api.mouser.com/api/v1/search/keyword?apiKey={self.apikey}', params=params, headers=headers, json=json_data)

        if response.json() is not None and 'Parts' not in response.json()['SearchResults']:
            print(response.json())

        if len(response.json()['SearchResults']['Parts']) == 0:
            offset = response.json()['SearchResults']['NumberOfResult']

        for product in response.json()['SearchResults']['Parts']:
            pprint.pprint('=================================================')

            if product['Category'] != "Ferrite Cores & Accessories":
                pprint.pprint(product['Description'])
                print("Wrong category")
                continue

            if product['Manufacturer'] == 'EPCOS / TDK':
                pprint.pprint(product['Description'])
                self.process_epcos_tdk_product(product)
            elif product['Manufacturer'] == 'TDK':
                pprint.pprint(product['Description'])
                self.process_epcos_tdk_product(product)
            elif product['Manufacturer'] == 'Fair-Rite':
                pprint.pprint(product['Description'])
                self.process_fair_rite_product(product)
            elif product['Manufacturer'] in self.ignored_manufacturers:
                print(product['Manufacturer'])
            else:
                continue
                # pprint.pprint(product)
                # assert 0, f"Unknown manufacturer: {product['Manufacturer']}"


        return {
            'current_offset' : offset + 50,
            'total': response.json()['SearchResults']['NumberOfResult']
        }

    def get_cores_stock(self):
        remaining = 1
        current_offset = 0

        if not os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") and os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            shutil.copy(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson")

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row
 
        self.core_data = pandas.DataFrame(self.core_data.values())
        self.core_data["manufacturer_reference"] = self.core_data.apply(lambda row: row['manufacturerInfo']['reference'].split(' ')[0], axis=1)

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            # pprint.pprint(f"self.core_data: {len(self.core_data.values())}")


        # pprint.pprint(self.core_data)
        # core_data = pandas.DataFrame(self.core_data.values())
        self.core_data = self.core_data.drop(['manufacturer_reference'], axis=1)
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson", "w")
        ndjson.dump(self.core_data.to_dict('records'), out_file)
        out_file.close()


class GatewayStocker(Stocker):

    def __init__(self):
        self.apikey = "6c28da81-2a3c-43b4-b03e-fca6847d5370"
        self.core_data = {}
        self.unfound_descriptions = []
        self.unfound_shapes = []
        self.ignored_manufacturers = [
                                        "Proterial",
                                        "Laird Performance Materials",
                                        "Wurth Elektronik",
                                        "United Chemi-Con",
                                        "Vacuumschmelze",
                                     ]

    def read_products(self, offset: int = 0):  # sourcery skip: extract-method
        excluded_words = ['ADJUSTER', 'BOBBIN', 'CLIP']

        try:
            data = pandas.read_excel(f"https://www.shop.gatewaycando.com/gateway_files/stock/gateway_stock_list.xls")
        except urllib.error.HTTPError:
            return {
                'current_offset' : 0,
                'total': 0
            }


        number = 0
        for index, datum in data.iterrows():

            if 'FERROXCUBE' in datum['Franchise']:
                core_manufacturer = 'Ferroxcube'
            if 'MAGNETICS' in datum['Franchise']:
                core_manufacturer = 'Magnetics'
            if 'TDK' in datum['Franchise']:
                core_manufacturer = 'TDK'
            if 'FAIR-RITE' in datum['Franchise']:
                core_manufacturer = 'Fair-Rite'
           
            try:
                invalid = False
                if 'FERROXCUBE' not in datum['Franchise'] and 'MAGNETICS' not in datum['Franchise'] and 'TDK' not in datum['Franchise'] and 'FAIR-RITE' not in datum['Franchise']:
                    continue
                for word in excluded_words:
                    if word in datum['Product'].split(' ') or word in datum['Description'].split(' '):
                        invalid = True
                        break

                if invalid:
                    continue
            except TypeError:
                continue

            manufacturer_part_number = datum['Product'].rstrip('*')
            print(manufacturer_part_number)
            print(manufacturer_part_number)
            print(manufacturer_part_number)
            core = self.core_data[self.core_data["manufacturer_reference"].str.contains(manufacturer_part_number)]
            if len(core.index) != 1:
                self.unfound_shapes.append(manufacturer_part_number)

            core_name = None
            if not core.empty:
                core_name = core.iloc[0]['name']
            else:

                shape = None
                material = None

                try:
                    length = 0
                    for shape_name in PyMKF.get_available_core_shapes():
                        if (shape_name in datum['Description'] or shape_name.upper().replace(' ', '') in datum['Description'].upper().replace('X', '/').replace(' ', '').replace('TOROID', 'T') or
                            shape_name in datum['Product'] or shape_name.upper().replace(' ', '') in datum['Product'].upper().replace('X', '/').replace(' ', '').replace('TOROID', 'T')):

                            if len(shape_name) > length:
                                length = len(shape_name)
                                shape = shape_name
                                if shape[0] == 'R':
                                    shape = shape.replace('R ', 'T ')

                except TypeError:
                    continue

                if shape is not None:
                    family = shape.split(' ')[0]

                    length = 0
                    for material_name in PyMKF.get_available_core_materials(core_manufacturer):
                        if material_name in datum['Description'] or material_name.upper().replace(' ', '') in datum['Description'].upper().replace(' ', ''):
                            material = material_name
                            if len(material_name) > length:
                                length = len(material_name)
                                material = material_name
                                number += 1

                    if material is None:
                        continue

                    core_data = self.verify_shape_is_in_database(family, shape, material, datum['Description'])

                    core_type = 'toroidal' if family == 'T' else 'two-piece set'


                    if ' gapped ' in datum['Description']:
                        number_gaps = 1
                        try:
                            al_value = float(datum['Product'].split("-A")[1].split("-")[0]) * 1e-9
                        except IndexError:
                            try:
                                al_value = float(datum['Product'].split("-G")[1].split("-")[0]) * 1e-9
                            except IndexError:
                                try:
                                    al_value = float(datum['Description'].split(" al ")[1]) * 1e-9
                                except IndexError:
                                    print(datum['Product'])
                                    print(datum['Description'])
                                    assert 0
                        gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)
                    else:
                        gapping = []

                    status = 'production'
                    coating = None
                    core_name = self.add_core(family=family,
                                              shape=shape,
                                              material=material,
                                              core_type=core_type,
                                              core_manufacturer=core_manufacturer,
                                              status=status,
                                              manufacturer_part_number=manufacturer_part_number,
                                              coating=coating,
                                              gapping=gapping)

            product_url = "https://www.shop.gatewaycando.com/magnetics/cores"
            distributor_reference = datum['Product']
            quantity = int(datum['Stock'])
            cost = None
            if core_name is not None:
                self.add_distributor(core_name=core_name,
                                     product_url=product_url,
                                     core_manufacturer=core_manufacturer,
                                     distributor_reference=distributor_reference,
                                     quantity=quantity,
                                     cost=cost,
                                     distributor_name="Gateway")

        return {
            'current_offset' : 0,
            'total': 0
        }

    def get_cores_stock(self):
        remaining = 1
        current_offset = 0

        if not os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") and os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            shutil.copy(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson")

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        self.core_data = pandas.DataFrame(self.core_data.values())
        self.core_data["manufacturer_reference"] = self.core_data.apply(lambda row: row['manufacturerInfo']['reference'].split(' ')[0], axis=1)

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            # pprint.pprint(f"self.core_data: {len(self.core_data.values())}")


        # pprint.pprint(self.core_data)
        # core_data = pandas.DataFrame(self.core_data.values())
        self.core_data = self.core_data.drop(['manufacturer_reference'], axis=1)
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson", "w")
        ndjson.dump(self.core_data.to_dict('records'), out_file)
        out_file.close()



if __name__ == '__main__':  # pragma: no cover
    digikeyStocker = DigikeyStocker()
    digikeyStocker.get_cores_stock()

    mouserStocker = MouserStocker()
    mouserStocker.get_cores_stock()

    gatewayStocker = GatewayStocker()
    gatewayStocker.get_cores_stock()

    cores = {}
    if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson"):
        with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_stock.ndjson") as f:
            previous_data = ndjson.load(f)
            for row in previous_data:
                cores[row['name']] = row
                cores[row['name']]['functionalDescription']['shape'] = cores[row['name']]['functionalDescription']['shape']['name']

    cores = pandas.DataFrame(cores.values())
    out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores.ndjson", "w")
    ndjson.dump(cores.to_dict('records'), out_file)
    out_file.close()

    pprint.pprint(mouserStocker.unfound_descriptions)
