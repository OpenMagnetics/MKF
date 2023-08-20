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

    def get_core_name(self, shape, material, gapping):
        if len(gapping) == 0:
            return f"{shape} - {material} - Ungapped"
        else:
            number_gaps = len([gap['length'] for gap in gapping if gap['type'] == 'subtractive'])
            if number_gaps == 1:
                return f"{shape} - {material} - Gapped {gapping[0]['length'] * 1000} mm"
            else:
                return f"{shape} - {material} - Distributed gapped {gapping[0]['length'] * 1000} mm"

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

    def add_core(self, family, shape, material, core_type, core_manufacturer, status, manufacturer_part_number, product_url, distributor_reference, quantity, cost, gapping=[], coating=None, datasheet=None, unique=False):
        core_name = self.get_core_name(shape, material, gapping)

        if core_name not in self.core_data:
            print(f"Adding {core_name}")
            self.core_data[core_name] = {
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
            if coating is not None:
                self.core_data[core_name]['functionalDescription']['coating'] = coating
        elif unique:
            assert 0, f"Already inserted core: {core_name}"

        if distributor_reference is not None:
            if isinstance(self.core_data[core_name]['distributorsInfo'], str):
                self.core_data[core_name]['distributorsInfo'] = ast.literal_eval(self.core_data[core_name]['distributorsInfo'])

            distributor_info = {
                "name": "Digi-Key",
                "link": product_url,
                "country": "USA",
                "distributedArea": "International",
                "reference": distributor_reference,
                "quantity": quantity,
                "updatedAt": date.today().strftime("%d/%m/%Y"),
                "cost": cost,
            }
            self.core_data[core_name]['distributorsInfo'].append(distributor_info)

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
                                "data": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1]
                            }
                        },
                        "voltage": {
                            "waveform": {
                                "data": [-10, 10, 10, 10, 10, 10, -10, -10, -10, -10]
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

            inductance = PyMKF.calculate_inductance_from_number_turns_and_gapping(core_data, 
                                                                            winding_data,
                                                                            operating_point,
                                                                            models)
            gap_length += 0.000001
            for index in range(len(core_data['functionalDescription']['gapping'])):
                core_data['functionalDescription']['gapping'][index]['length'] = round(core_data['functionalDescription']['gapping'][index]['length'], 5)

        if abs(previous_inductance - al_value) < abs(inductance - al_value):
            for index in range(len(core_data['functionalDescription']['gapping'])):
                if (core_data['functionalDescription']['gapping'][index]['type'] == 'subtractive'):
                    core_data['functionalDescription']['gapping'][index]['length'] = round(previous_gap_length, 5)

        return core_data['functionalDescription']['gapping']

    def adaptive_round(self, value):
        if value < 4:
            return round(value, 2)
        elif value < 10:
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
            return text.replace(" ", "").replace(",", ".").upper().replace("TX", "T").replace("TC", "T").replace("X", "/")

        for shape_name in PyMKF.get_available_core_shapes():
            if shape_name in text:
                return shape_name

        fixed_text = fix(text)
        for shape_name in PyMKF.get_available_core_shapes():
            if fix(shape_name) in fixed_text:
                return shape_name

        print("fixed_text")
        print(fixed_text)
        return None

    def try_get_material(self, text):  # sourcery skip: use-next
        def fix(text):
            return text.replace(" 0", " ").replace(" ", "").replace(",", ".").replace("Mµ", "Mu").upper()

        for material_name in PyMKF.get_available_core_materials():
            if material_name in text:
                return material_name

        fixed_text = fix(text)
        for material_name in PyMKF.get_available_core_materials():
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
        self.ignored_manufacturers = ["VACUUMSCHMELZE",
                                      "Omron Automation and Safety",
                                      "Toshiba Semiconductor and Storage",
                                      "United Chemi-Con",
                                      "Würth Elektronik",
                                      "Laird-Signal Integrity Products"]

    def get_refresh_token(self):
        # go to the browser and paste:
        # https://sandbox-api.digikey.com/v1/oauth2/authorize?response_type=code&client_id=Cpw48RBF8aYMTiTWeAYLvyUgPo2Db3F7&redirect_uri=https://localhost:8139/digikey_callback
        # https://api.digikey.com/v1/oauth2/authorize?response_type=code&client_id=cN8i6L6KnNGJB2h3zsQgC7KvWf8AccsC&redirect_uri=https://localhost:8139/digikey_callback
        # get auth code from the redirected URL and paste it here:

        auth_code = 'BxAnaFf0'

        token_response = requests.post('https://api.digikey.com/v1/oauth2/token', data={
            'client_id': self.client_id,
            'client_secret': self.client_secret,
            'redirect_uri': self.redirect_uri,
            'code': auth_code,
            'grant_type': 'authorization_code'
        })

        refresh_token = token_response.json()['refresh_token']
        print(token_response)

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
                assert 0, f"Unknown manufacturer: {product['Manufacturer']['Value']}"

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

        shape_type = "two-piece set"
        material = None
        gap_length = None
        for parameter in product['Parameters']:
            if parameter['Parameter'] == 'Core Type' and parameter['Value'] == 'Toroid':
                shape_type = 'toroidal'
            if parameter['Parameter'] == 'Material':
                material = self.try_get_material(parameter['Value'])
                if material in not_included_materials:
                    return None

        # Some cores have typos in the series value, we correct them with the Id
        for typo in typos:
            if product['Series']['ValueId'] == typo['typo']:
                product['Series']['Value'] = typo['corrected_shape']

        family = self.try_get_family(product['Series']['Value'])

        if material is None:
            return

        if family is not None:
            try: 
                shape = product['Series']['Value'].replace(family, f"{family} ")
                print("shape")
                print(shape)
            except ValueError:
                print(family)
                print(product['Series']['Value'])
                print(product['Series']['Value'].replace(family, f"{family} "))
                assert 0

            if family in not_included_families:
                return

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


            # We check that we have the core in the database
            core_data = self.verify_shape_is_in_database(family, shape, material, product)

            core_type = 'toroidal' if family == 'T' else 'two-piece set'
            core_manufacturer = 'TDK'
            manufacturer_part_number = product['ManufacturerPartNumber']
            product_url = product['ProductUrl']
            distributor_reference = product['DigiKeyPartNumber']
            quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
            cost = float(product['UnitPrice'])
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

            self.add_core(family=family,
                          shape=shape,
                          material=material,
                          core_type=core_type,
                          core_manufacturer=core_manufacturer,
                          status=status,
                          manufacturer_part_number=manufacturer_part_number,
                          product_url=product_url,
                          distributor_reference=distributor_reference,
                          quantity=quantity,
                          cost=cost,
                          gapping=gapping)

        else:
            pprint.pprint(product)
            assert 0

    def process_fair_rite_product(self, product):  # sourcery skip: class-extract-method, extract-method, low-code-quality, move-assign
        exceptions = []
        not_included_families = ['ROD', 'Planar E/I', 'PS', 'Multi-Hole (2)', 'Multi-Hole', 'I', 'Tube']
        typos = []
        shape_type = "two-piece set"
        family = None
        material = None
        gap_length = None
        effective_length = None
        effective_area = None
        effective_volume = None
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
                material = self.try_get_material(parameter['Value'])
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

        # pprint.pprint(product)
        # pprint.pprint(shape)
        # assert 0
        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping)

    def process_tdk_product(self, product):  # sourcery skip: class-extract-method, extract-method, low-code-quality, move-assign
        not_included_families = ["EPC", "EER", "ELT", "PQI", "EI", "SP", "-"]
        typos = [{"typo": "EI 13/4.4", "corrected_shape": "EL 13/4.4", "corrected_family": "EL"},
                 {"typo": "EI 15.5/5.8", "corrected_shape": "EL 15.5/5.8", "corrected_family": "EL"},
                 {"typo": "EI 14/4.5/9", "corrected_shape": "ER 14/4.5/9", "corrected_family": "ER"},
                 {"typo": "EIR 14/4.5/9", "corrected_shape": "ER 14/4.5/9", "corrected_family": "ER"},
                 {"typo": "EI 18/5/12", "corrected_shape": "ER 18/5/12", "corrected_family": "ER"},
                 {"typo": "EIR 18/5/12", "corrected_shape": "ER 18/5/12", "corrected_family": "ER"},
                 {"typo": "EI 22/5.5/15", "corrected_shape": "ER 22/5.5/15", "corrected_family": "ER"},
                 {"typo": "EIR 22/5.5/15", "corrected_shape": "ER 22/5.5/15", "corrected_family": "ER"}]
        corrections = ["P 47/28"]
        corrections_family = ["P"]
        not_included_materials = ["H5C2", "PC44"]

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
            product_url = product['ProductUrl']
            distributor_reference = product['DigiKeyPartNumber']
            quantity = int(product['QuantityAvailable'] * 2 if '2pcs' in product['ProductDescription'].lower() else product['QuantityAvailable'])
            cost = float(product['UnitPrice'])
            gapping = []

            for parameter in product['Parameters']:
                if parameter['Parameter'] == 'Gap' and parameter['Value'] not in ['-', 'Ungapped']:
                    number_gaps = 3 if parameter['Value'] == 'Distributed Gapped' else 1
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

            self.add_core(family=family,
                          shape=shape,
                          material=material,
                          core_type=core_type,
                          core_manufacturer=core_manufacturer,
                          status=status,
                          manufacturer_part_number=manufacturer_part_number,
                          product_url=product_url,
                          distributor_reference=distributor_reference,
                          quantity=quantity,
                          cost=cost,
                          gapping=gapping)

        else:
            return

    def process_ferroxcube_product(self, product):
        # sourcery skip: extract-method, low-code-quality, move-assign
        not_included_families = ["I", "PLT", "ROD"]
        typos = [{"typo": "ER 35W/21/11", "corrected_shape": "ER 35/21/11", "corrected_family": "ER"}]
        not_included_materials = ["3D", "3E", "3C11", "4C"]

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
                material = self.try_get_material(product['ManufacturerPartNumber'])
                shape = self.try_get_shape(product['ManufacturerPartNumber'])

                if material is None:
                    for not_included_material in not_included_materials:
                        if not_included_material in product['ManufacturerPartNumber']:
                            return None
                    assert 0, f"Unknown material in {product['ManufacturerPartNumber']}"

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

            status = 'production'
            core_manufacturer = 'Ferroxcube'
            core_type = 'toroidal' if family.lower() in ['t', 'tx', 'tc'] else 'two-piece set'
            manufacturer_part_number = product['ManufacturerPartNumber']
            product_url = product['ProductUrl']
            distributor_reference = product['DigiKeyPartNumber']
            quantity = int(product['QuantityAvailable'] * 2 if '2PC' in product['ProductDescription'].lower() else product['QuantityAvailable'])
            cost = float(product['UnitPrice'])
            gapping = []

            for parameter in product['Parameters']:
                if parameter['Parameter'] == 'Gap' and parameter['Value'] not in ['-', 'Ungapped']:
                    number_gaps = 3 if parameter['Value'] == 'Distributed Gapped' else 1
                    try:
                        al_value = float(product['ManufacturerPartNumber'].split("-A")[1].split("-")[0]) * 1e-9
                    except IndexError:
                        try:
                            al_value = float(product['ManufacturerPartNumber'].split("-E")[1].split("-")[0]) * 1e-9
                        except IndexError:
                            pprint.pprint(product)
                            print(product['ManufacturerPartNumber'])
                            assert 0

                    gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)

                if parameter['Parameter'] == 'Part Status' and parameter['Value'] != 'Active':
                    assert 0, f"Part Status is {parameter['Value']}"


            self.add_core(family=family,
                          shape=shape,
                          material=material,
                          core_type=core_type,
                          core_manufacturer=core_manufacturer,
                          status=status,
                          manufacturer_part_number=manufacturer_part_number,
                          product_url=product_url,
                          distributor_reference=distributor_reference,
                          quantity=quantity,
                          cost=cost,
                          gapping=gapping)
        else:
            assert 0

    def process_magnetics_product(self, product):
        not_included_materials = ["75-Series 125",
                                  "AmoFlux",
                                  "High Flux 19",
                                  "High Flux 90",
                                  "Kool Mu 19",
                                  "MPP 250",
                                  "MPP 367",
                                  "XFlux 14",
                                  "XFlux 50"]
        # pprint.pprint(product)
        try:
            material = self.try_get_material(self.get_parameter(product['Parameters'], 'Material') + " " + self.get_parameter(product['Parameters'], 'Initial Permeability (µi)'))
        except TypeError:
            material = None

        if self.get_parameter(product['Parameters'], 'Core Type') is None:
            return

        if self.get_parameter(product['Parameters'], 'Core Type') == 'Toroid':
            external_diameter = self.get_parameter(product['Parameters'], 'Diameter').split('(')[1].split('mm)')[0]
            height = self.get_parameter(product['Parameters'], 'Height').split('(')[1].split('mm)')[0]
            family = 'T'
            pprint.pprint("external_diameter")
            pprint.pprint(external_diameter)
            pprint.pprint("height")
            pprint.pprint(height)
            temptative_shape = f"{family} {self.adaptive_round(float(external_diameter) + 0.7)}/{None}/{self.adaptive_round(float(height) + 0.7)}"
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
                pprint.pprint(product)
                assert 0
        shape = self.find_shape_closest_dimensions(family, temptative_shape, limit=0.1)
        pprint.pprint(family)
        pprint.pprint(material)
        pprint.pprint(shape)

        if material is None:
            material = self.try_get_material(product['ProductDescription'])

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
        product_url = product['ProductUrl']
        distributor_reference = product['DigiKeyPartNumber']
        quantity = int(product['QuantityAvailable'] * 2 if '2PC' in product['ProductDescription'].lower() else product['QuantityAvailable'])
        cost = float(product['UnitPrice'])
        gapping = []

        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping)
        # if self.get_parameter(product['Parameters'], 'Core Type') == 'Toroid':
        #     approximated_shape = f"T {}/{}/{}"
        # assert 0

    def get_cores_stock(self):
        # access_token = self.get_refresh_token()
        access_token = self.get_access_token()

        remaining = 1
        current_offset = 0

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(access_token, current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")
            # break

        # pprint.pprint(self.core_data)
        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()


class MouserStocker(Stocker):

    def __init__(self):
        self.apikey = "6c28da81-2a3c-43b4-b03e-fca6847d5370"
        self.core_data = {}
        self.unfound_descriptions = []
        self.ignored_manufacturers = [
                                        "Proterial",
                                        "Laird Performance Materials",
                                        "Wurth Elektronik",
                                        "United Chemi-Con",
                                        "Vacuumschmelze",
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

        material = self.try_get_material(product['Description'].split('Ferrite Cores & Accessories')[1].split(' ')[1])
        pprint.pprint(product['Description'].split(material)[1])
        pprint.pprint(product['Description'].replace(material, '').replace("Planar EE", "E").split("Core "))
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
            extra_data = pandas.read_html(response.text)[8]
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

        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        core_manufacturer = 'Fair-Rite'
        status = 'production'
        manufacturer_part_number = product['ManufacturerPartNumber']
        product_url = product['ProductDetailUrl']
        distributor_reference = product['MouserPartNumber']
        quantity = int(product['AvailabilityInStock'])
        for price_break in product['PriceBreaks']:
            if price_break['Quantity'] == 1:
                cost = float(price_break['Price'].split("$")[1])
                break
        else:
            assert 0

        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping)



    def process_epcos_tdk_product(self, product):
        not_included_families = ["EV", "MHB", "PCH", "PS", "I", "ILP", "QU", "TT", "B", "PTS"]
        not_included_materials = ["N48", "N22", "T38", "M33", "T65", "T46", "T37", "K10", "K1"]
        typos = [
                    {"typo": "TOROID", "correction": "T"},
                    {"typo": "EER", "correction": "ER"},
                ]

        def get_core_data(description):
            family = self.try_get_family(description)

            shape = None
            if family is None:
                family, shape = self.try_get_family_and_shape_statistically(description)

            if family is None:
                pprint.pprint(product)
                self.unfound_descriptions.append(description)
                return None
                # assert 0, "Unknown family"

            if family not in PyMKF.get_available_shape_families() and family not in not_included_families:
                pprint.pprint(product)
                self.unfound_descriptions.append(description)
                return None
                # assert 0, "Unknown family"

            if shape is None:
                shape = self.try_get_shape(description)

            if shape is None:
                pprint.pprint(product)
                self.unfound_descriptions.append(description)
                return None
                # assert 0, "Unknown shape"

            material = self.try_get_material(description)

            if material is None:
                for not_included_material in not_included_materials:
                    if not_included_material in description:
                        return None
                pprint.pprint(product)
                self.unfound_descriptions.append(description)
                return None
                # assert 0, f"Unknown material: {material}"

            # We check that we have the core in the database
            core_data = self.verify_shape_is_in_database(family, shape, material, product)



            try:
                if "-DG" in description:
                    number_gaps = 3
                    al_value = float(description.split("-DG")[1].strip()) * 1e-9
                else:
                    number_gaps = 1
                    al_value = float(description.split(material)[1].strip().split(" ")[0]) * 1e-9
                if al_value > 0:
                    gapping = self.get_gapping(core_data, "TDK", al_value, number_gaps)
                else:
                    gapping = []
            except ValueError:
                gapping = []

            return family, shape, material, gapping

        exceptions = ["COIL FORMER",
                      "Coil former",
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

        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        core_manufacturer = 'TDK'
        status = 'production'
        manufacturer_part_number = product['ManufacturerPartNumber']
        product_url = product['ProductDetailUrl']
        distributor_reference = product['MouserPartNumber']
        quantity = int(product['AvailabilityInStock'])
        for price_break in product['PriceBreaks']:
            if price_break['Quantity'] == 1:
                cost = float(price_break['Price'].split("$")[1])
                break
        else:
            assert 0

        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping)

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
                pprint.pprint(product)
                assert 0, f"Unknown manufacturer: {product['Manufacturer']}"


        return {
            'current_offset' : offset + 50,
            'total': response.json()['SearchResults']['NumberOfResult']
        }

    def get_cores_stock(self):
        remaining = 1
        current_offset = 0

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")


        # pprint.pprint(self.core_data)
        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()



class GatewayStocker(Stocker):

    def __init__(self):
        self.apikey = "6c28da81-2a3c-43b4-b03e-fca6847d5370"
        self.core_data = {}
        self.unfound_descriptions = []
        self.ignored_manufacturers = [
                                        "Proterial",
                                        "Laird Performance Materials",
                                        "Wurth Elektronik",
                                        "United Chemi-Con",
                                        "Vacuumschmelze",
                                     ]

    def read_products(self, offset: int = 0):  # sourcery skip: extract-method
        excluded_words = ['ADJUSTER', 'BOBBIN', 'CLIP']

        data = pandas.read_excel(f"https://gatewaycando.com/NewStockList/Gateway%20Stock%20File.xlsx")

        number = 0
        for index, datum in data.iterrows():
            try:
                invalid = False
                if 'FERROXCUBE' not in datum['Franchise'] and 'MAGNETICS' not in datum['Franchise'] and 'TDK' not in datum['Franchise'] and 'FAIR-RITE' not in datum['Franchise']:
                    continue
                for word in excluded_words:
                    if word in datum['Part'].split(' ') or word in datum['Description'].split(' '):
                        invalid = True
                        break

                if invalid:
                    continue
            except TypeError:
                continue


            shape = None
            material = None

            try:
                length = 0
                for shape_name in PyMKF.get_available_core_shapes():
                    if (shape_name in datum['Description'] or shape_name.upper().replace(' ', '') in datum['Description'].upper().replace('X', '/').replace(' ', '').replace('TOROID', 'T') or
                        shape_name in datum['Part'] or shape_name.upper().replace(' ', '') in datum['Part'].upper().replace('X', '/').replace(' ', '').replace('TOROID', 'T')):

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
                for material_name in PyMKF.get_available_core_materials():
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
                if 'FERROXCUBE' in datum['Franchise']:
                    core_manufacturer = 'Ferroxcube'
                if 'MAGNETICS' in datum['Franchise']:
                    core_manufacturer = 'Magnetics'
                if 'TDK' in datum['Franchise']:
                    core_manufacturer = 'TDK'
                if 'FAIR-RITE' in datum['Franchise']:
                    core_manufacturer = 'Fair-Rite'
                manufacturer_part_number = datum['Part']
                product_url = datum['www link']
                distributor_reference = datum['Part']
                quantity = int(datum['Qty'])
                cost = None

                if ' gapped ' in datum['Description']:
                    number_gaps = 1
                    try:
                        al_value = float(datum['Part'].split("-A")[1].split("-")[0]) * 1e-9
                    except IndexError:
                        try:
                            al_value = float(datum['Part'].split("-G")[1].split("-")[0]) * 1e-9
                        except IndexError:
                            pprint.pprint(product)
                            print(datum['Part'])
                            assert 0
                    gapping = self.get_gapping(core_data, core_manufacturer, al_value, number_gaps)
                else:
                    gapping = []

                status = 'production'

                self.add_core(family=family,
                              shape=shape,
                              material=material,
                              core_type=core_type,
                              core_manufacturer=core_manufacturer,
                              status=status,
                              manufacturer_part_number=manufacturer_part_number,
                              product_url=product_url,
                              distributor_reference=distributor_reference,
                              quantity=quantity,
                              cost=cost,
                              gapping=gapping)


        return {
            'current_offset' : 0,
            'total': 0
        }

    def get_cores_stock(self):
        remaining = 1
        current_offset = 0

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")


        # pprint.pprint(self.core_data)
        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()


class FerroxcubeInventory(Stocker):
    def remove_current_inventory(self):
        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            os.remove(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson")

    def read_concentric_products(self, offset: int = 0):
        cookies = {
            '_ga': 'GA1.1.382486033.1678131350',
            'PHPSESSID': 'q1c6dif7a9uqvc2agv7ql2mur8',
            '_ga_KDJ2X3F53E': 'GS1.1.1681041728.13.1.1681041845.0.0.0',
        }

        headers = {
            'authority': 'www.ferroxcube.com',
            'accept': 'application/json, text/javascript, */*; q=0.01',
            'accept-language': 'es-ES,es;q=0.9,en;q=0.8',
            'content-type': 'application/x-www-form-urlencoded; charset=UTF-8',
            # 'cookie': '_ga=GA1.1.382486033.1678131350; PHPSESSID=q1c6dif7a9uqvc2agv7ql2mur8; _ga_KDJ2X3F53E=GS1.1.1681041728.13.1.1681041845.0.0.0',
            'dnt': '1',
            'origin': 'https://www.ferroxcube.com',
            'referer': 'https://www.ferroxcube.com/en-global/products_ferroxcube/stepTwo/shape_cores_accessories?s_sel=&series_sel=&material_sel=&material=&part=',
            'sec-ch-ua': '"Chromium";v="110", "Not A(Brand";v="24"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'sec-fetch-dest': 'empty',
            'sec-fetch-mode': 'cors',
            'sec-fetch-site': 'same-origin',
            'sec-gpc': '1',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36',
            'x-requested-with': 'XMLHttpRequest',
        }

        data = {
            'c_sn': '25',
            's_sn': '',
            'series_sn': '',
            'm_no': '',
            'm_value': '',
            'part': '',
            'page': offset,
            'lef_aef_min': '',
            'lef_aef_max': '',
            'lef_min': '',
            'lef_max': '',
            'aef_min': '',
            'aef_max': '',
            'vef_min': '',
            'vef_max': '',
            'aw_min': '',
            'aw_max': '',
            'al_value_min': '',
            'al_value_max': '',
            'gap_value_min': '',
            'gap_value_max': '',
        }

        response = requests.post(
            'https://www.ferroxcube.com/en-global/products_ferroxcube/ajaxFilterData/shape_cores_accessories',
            cookies=cookies,
            headers=headers,
            data=data,
        )

        clean_response = self.clean_html(str(response.json()))

        total_pages = int(re.findall(r'id="page_count" data-page="(.*?)"', clean_response)[0])
        self.process(clean_response)

        return {
            'current_offset' : offset + 1,
            'total': total_pages
        }

    def read_toroidal_products(self, offset: int = 0):
        cookies = {
            '_ga': 'GA1.1.382486033.1678131350',
            'PHPSESSID': '3c5rt12r04uopj0v5cv3c43p1f',
            '_ga_KDJ2X3F53E': 'GS1.1.1681074872.15.1.1681075275.0.0.0',
        }

        headers = {
            'authority': 'www.ferroxcube.com',
            'accept': 'application/json, text/javascript, */*; q=0.01',
            'accept-language': 'es-ES,es;q=0.9,en;q=0.8',
            'content-type': 'application/x-www-form-urlencoded; charset=UTF-8',
            # 'cookie': '_ga=GA1.1.382486033.1678131350; PHPSESSID=3c5rt12r04uopj0v5cv3c43p1f; _ga_KDJ2X3F53E=GS1.1.1681074872.15.1.1681075275.0.0.0',
            'dnt': '1',
            'origin': 'https://www.ferroxcube.com',
            'referer': 'https://www.ferroxcube.com/en-global/products_ferroxcube/stepTwo/toroid?s_sel=&series_sel=&material_sel=3C90&material=&part=',
            'sec-ch-ua': '"Chromium";v="110", "Not A(Brand";v="24"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'sec-fetch-dest': 'empty',
            'sec-fetch-mode': 'cors',
            'sec-fetch-site': 'same-origin',
            'sec-gpc': '1',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36',
            'x-requested-with': 'XMLHttpRequest',
        }

        data = {
            'c_sn': '28',
            's_sn': '',
            'series_sn': '',
            'm_no': '',
            'm_value': '',
            'part': '',
            'page': offset,
            'lef_aef_min': '',
            'lef_aef_max': '',
            'lef_min': '',
            'lef_max': '',
            'aef_min': '',
            'aef_max': '',
            'vef_min': '',
            'vef_max': '',
            'aw_min': '',
            'aw_max': '',
            'al_value_min': '',
            'al_value_max': '',
            'gap_value_min': '',
            'gap_value_max': '',
        }

        response = requests.post(
            'https://www.ferroxcube.com/en-global/products_ferroxcube/ajaxFilterData/toroid',
            cookies=cookies,
            headers=headers,
            data=data,
        )

        clean_response = self.clean_html(str(response.json()))

        total_pages = int(re.findall(r'id="page_count" data-page="(.*?)"', clean_response)[0])
        self.process(clean_response)

        return {
            'current_offset' : offset + 1,
            'total': total_pages
        }

    def process(self, data):
        constants = PyMKF.get_constants()
        not_included_materials = ["3D", "3E", "3C11", "4C", "3H"]
        typos = [{"typo": "ER35W/21/11", "correction": "ER35/21/11"}]


        core_names = re.findall(r'data-title="P/N">(.*?)</td>', data)
        core_datasheet_urls = re.findall(r'media/product/file/Pr_ds(.*?)"', data)
        for index in range(len(core_datasheet_urls)):
            core_datasheet_urls[index] = "https://ferroxcube.com/upload/media/product/file/Pr_ds" + core_datasheet_urls[index]
            core_datasheet_urls[index] = core_datasheet_urls[index].replace(" ", "%20")
        gap_lengths = re.findall(r'µm]">(.*?)</td>', data.replace("≈", ""))
        for index in range(len(gap_lengths)):
            gap_lengths[index] = gap_lengths[index].split("±")[0]
        # assert len(core_names) == len(core_datasheet_urls), data
        assert len(core_names) == len(gap_lengths), data

        assert len(core_names) % 10 == 0, f"len(core_names): {len(core_names)}"

        status = 'production'
        core_manufacturer = "Ferroxcube"
        for core_name_index, core_name in enumerate(core_names):
            if "/R-" in core_name or "-P" in core_name or "PLT" in core_name:
                print(f"Avoiding: {core_name}")
                continue
            if core_name[0] == 'I':
                print(f"Avoiding: {core_name}")
                continue
            if "-E" in core_name or "-A" in core_name:
                check_uniqueness = False
            else:
                check_uniqueness = False
            print(core_name)

            # Some cores have typos in the series value, we correct them with the name
            for typo in typos:
                core_name.replace(typo['typo'], typo['correction'])

            family = self.try_get_family(core_name)
            shape = self.try_get_shape(core_name)
            material = self.try_get_material(core_name)
            if family == 'TX':
                family = 'T'

            print(family)
            print(shape)
            print(material)
            print(gap_lengths[core_name_index])

            if shape is None:
                family, shape = self.try_get_family_and_shape_statistically(core_name)

            if material is None:
                for not_included_material in not_included_materials:
                    if not_included_material in core_name:
                        return None
                assert 0, f"Unknown material in {core_name}"

            if family is None:
                return 
            core_data = self.verify_shape_is_in_database(family, shape, material, core_name)
            core_type = 'toroidal' if family == 'T' else 'two-piece set'
            manufacturer_part_number = core_name
            if len(core_names) == len(core_datasheet_urls):
                print(core_datasheet_urls[core_name_index])
                datasheet = core_datasheet_urls[core_name_index]
            else:
                datasheet = None
            product_url = None
            distributor_reference = None
            quantity = None
            cost = None
            coating = None
            if family == 'T':
                if 'TX' in core_name:
                    coating = "epoxy"

            if gap_lengths[core_name_index] == 'ungapped':
                gapping = []
            else:
                gapping = [
                    {
                        "type": "subtractive",
                        "length": float(gap_lengths[core_name_index]) / 1000000
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

            self.add_core(family=family,
                          shape=shape,
                          material=material,
                          core_type=core_type,
                          core_manufacturer=core_manufacturer,
                          status=status,
                          manufacturer_part_number=manufacturer_part_number,
                          product_url=product_url,
                          distributor_reference=distributor_reference,
                          quantity=quantity,
                          cost=cost,
                          gapping=gapping,
                          coating=coating,
                          datasheet=datasheet,
                          unique=check_uniqueness)


    def get_cores_inventory(self):
        remaining = 1
        current_offset = 1

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_concentric_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_toroidal_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")

        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()


class TdkInventory(Stocker):
    def read_products(self, offset: int = 0):
        cookies = {
            '_gcl_au': '1.1.1345231983.1676319639',
            '_ga': 'GA1.1.852533192.1676319639',
            'wooTracker': 'sSdLALgdnDeE',
            'mhX5qX4': '802dc24a-3f2d-45d6-a0ed-a5a6e2185d38',
            'MHD6Ph5HlS': '1ec7f994-87ed-aaca-aa3e-a724a6a0a00f',
            'OptanonAlertBoxClosed': '2023-02-18T16:25:06.868Z',
            '_ga_W48MFFWN91': 'GS1.1.1679614622.1.1.1679614642.0.0.0',
            '__lt__cid.42694821': '5a1e72da-283d-434a-98c5-945ac91af208',
            '_ga_K1F1PVYW3G': 'GS1.1.1679944616.3.0.1679944617.0.0.0',
            '_ga_PWM9L9Z85D': 'GS1.1.1679944616.3.0.1679944617.0.0.0',
            'ref_data': '%7B%22display_distributor_name%22%3A%22Mouser%20Electronics%22%2C%22short_name%22%3A%22Mouser%20Electronics%22%2C%22net_components_distributor_name%22%3A%22Mouser%20Electronics%20Inc.%22%2C%22host%22%3A%22www.mouser.com%22%2C%22ref_disty%22%3A%22mouser%22%7D',
            '_ga_J72MDETYME': 'GS1.1.1681037908.8.0.1681037911.0.0.0',
            'OptanonConsent': 'isIABGlobal=false&datestamp=Sun+Apr+09+2023+13%3A55%3A22+GMT%2B0200+(hora+de+verano+de+Europa+central)&version=6.6.0&hosts=&consentId=911d4900-ea39-4349-a2c6-005b068c07e0&interactionCount=2&landingPath=NotLandingPage&groups=C001%3A1%2CC002%3A1%2CC003%3A1%2CC004%3A1%2CC006%3A1&AwaitingReconsent=false&geolocation=ES%3BMD',
            '_ga_DEXMV98BYD': 'GS1.1.1681041255.16.1.1681041383.60.0.0',
            '_ga_L03JBTW2NL': 'GS1.1.1681041255.19.1.1681041383.0.0.0',
            'product_search_fw': 'eyJpdiI6IllUbkZ3ck4wZDdlU2xTWUJTUUd3SHc9PSIsInZhbHVlIjoiMGQyME9QQUJuM2ZKa0lQYy83UW5VZWcyT21XOWtXZjB5U29IMUtPOWhEVVdmbzBYcmRqK2tUQTBIY1pyaUJBZHhoVVJDQ3ZzWlliRFBqT3Y4OGRLTVlZT0hpRFgrL2JtK0dhTDA4eEZaSXVMM1VjTTFwZk5kQTNlT2NQT0tldjkiLCJtYWMiOiJlYTU4YTBmY2M0YjczZjhmZGRkOThmNzBhZGE2ODdmYjVjYTQzZWVjNTA0ZTgzNzJhZTk4NWQxM2NhMjc0ZmIwIiwidGFnIjoiIn0%3D',
            'mh8ciXI-3a7HfloWZ': '',
            'mh8oobaIFjn3V': '',
            'mh8oobaJ7gi': '',
            'mh8oobaC7SeJd': '',
            'mhX5qk95': '02d83a15-f22f-4fad-adef-6f02a7662642',
            'ak_bmsc': '9EDD2353602C93C4006855E989CE5EA6~000000000000000000000000000000~YAAQH1UQYBxqeliHAQAAkdDlZhM/lsFidJ39C9LXsGDXfp07dkJdIsPFu3Qx7/0mxxSuuocidaAdTRIvtHA3w9sQD9Cp+DBQndi0viFWWwZc0WlwiA4IxTHwgJWcfn0LbGYQz2aHR7LA0L7le7OumTWqTCQKjqiFrIAXyTlVHY70P6Aaz/MmReJe6tY+zitAYFQ4suqHI9GG67aJzZX4xidxewduzQ9wEoAVPkgDO912URLQ9PL8fHhn3UZnDKUM1NonN7NUoZE1ZLtEjmzKlm0uqPkkGf7dsNhzrKIhuXscAi9PK33iS5Zh2Ua3QmWxStzGk0wzLu6Eiungotl/Xe/mr+5AV8Wev+ZrnQq0nsqFGWhQZOh7cxO8P0/sw4VP0aPd+xMmA6c=',
            'RT': '"z=1&dm=tdk.com&si=0290ba7e-33a6-4ddf-8925-8448e3bd6ce4&ss=lg9clkda&sl=0&tt=0&bcn=%2F%2F684dd313.akstat.io%2F&nu=x9hkqc7&cl=a9kna&ul=a9knr"',
        }

        headers = {
            'authority': 'product.tdk.com',
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
            'accept-language': 'es-ES,es;q=0.9,en;q=0.8',
            # 'cookie': '_gcl_au=1.1.1345231983.1676319639; _ga=GA1.1.852533192.1676319639; wooTracker=sSdLALgdnDeE; mhX5qX4=802dc24a-3f2d-45d6-a0ed-a5a6e2185d38; MHD6Ph5HlS=1ec7f994-87ed-aaca-aa3e-a724a6a0a00f; OptanonAlertBoxClosed=2023-02-18T16:25:06.868Z; _ga_W48MFFWN91=GS1.1.1679614622.1.1.1679614642.0.0.0; __lt__cid.42694821=5a1e72da-283d-434a-98c5-945ac91af208; _ga_K1F1PVYW3G=GS1.1.1679944616.3.0.1679944617.0.0.0; _ga_PWM9L9Z85D=GS1.1.1679944616.3.0.1679944617.0.0.0; ref_data=%7B%22display_distributor_name%22%3A%22Mouser%20Electronics%22%2C%22short_name%22%3A%22Mouser%20Electronics%22%2C%22net_components_distributor_name%22%3A%22Mouser%20Electronics%20Inc.%22%2C%22host%22%3A%22www.mouser.com%22%2C%22ref_disty%22%3A%22mouser%22%7D; _ga_J72MDETYME=GS1.1.1681037908.8.0.1681037911.0.0.0; OptanonConsent=isIABGlobal=false&datestamp=Sun+Apr+09+2023+13%3A55%3A22+GMT%2B0200+(hora+de+verano+de+Europa+central)&version=6.6.0&hosts=&consentId=911d4900-ea39-4349-a2c6-005b068c07e0&interactionCount=2&landingPath=NotLandingPage&groups=C001%3A1%2CC002%3A1%2CC003%3A1%2CC004%3A1%2CC006%3A1&AwaitingReconsent=false&geolocation=ES%3BMD; _ga_DEXMV98BYD=GS1.1.1681041255.16.1.1681041383.60.0.0; _ga_L03JBTW2NL=GS1.1.1681041255.19.1.1681041383.0.0.0; product_search_fw=eyJpdiI6IllUbkZ3ck4wZDdlU2xTWUJTUUd3SHc9PSIsInZhbHVlIjoiMGQyME9QQUJuM2ZKa0lQYy83UW5VZWcyT21XOWtXZjB5U29IMUtPOWhEVVdmbzBYcmRqK2tUQTBIY1pyaUJBZHhoVVJDQ3ZzWlliRFBqT3Y4OGRLTVlZT0hpRFgrL2JtK0dhTDA4eEZaSXVMM1VjTTFwZk5kQTNlT2NQT0tldjkiLCJtYWMiOiJlYTU4YTBmY2M0YjczZjhmZGRkOThmNzBhZGE2ODdmYjVjYTQzZWVjNTA0ZTgzNzJhZTk4NWQxM2NhMjc0ZmIwIiwidGFnIjoiIn0%3D; mh8ciXI-3a7HfloWZ=; mh8oobaIFjn3V=; mh8oobaJ7gi=; mh8oobaC7SeJd=; mhX5qk95=02d83a15-f22f-4fad-adef-6f02a7662642; ak_bmsc=9EDD2353602C93C4006855E989CE5EA6~000000000000000000000000000000~YAAQH1UQYBxqeliHAQAAkdDlZhM/lsFidJ39C9LXsGDXfp07dkJdIsPFu3Qx7/0mxxSuuocidaAdTRIvtHA3w9sQD9Cp+DBQndi0viFWWwZc0WlwiA4IxTHwgJWcfn0LbGYQz2aHR7LA0L7le7OumTWqTCQKjqiFrIAXyTlVHY70P6Aaz/MmReJe6tY+zitAYFQ4suqHI9GG67aJzZX4xidxewduzQ9wEoAVPkgDO912URLQ9PL8fHhn3UZnDKUM1NonN7NUoZE1ZLtEjmzKlm0uqPkkGf7dsNhzrKIhuXscAi9PK33iS5Zh2Ua3QmWxStzGk0wzLu6Eiungotl/Xe/mr+5AV8Wev+ZrnQq0nsqFGWhQZOh7cxO8P0/sw4VP0aPd+xMmA6c=; RT="z=1&dm=tdk.com&si=0290ba7e-33a6-4ddf-8925-8448e3bd6ce4&ss=lg9clkda&sl=0&tt=0&bcn=%2F%2F684dd313.akstat.io%2F&nu=x9hkqc7&cl=a9kna&ul=a9knr"',
            'dnt': '1',
            'referer': 'https://product.tdk.com/en/search/ferrite/ferrite/ferrite-core/list',
            'sec-ch-ua': '"Chromium";v="110", "Not A(Brand";v="24"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'sec-fetch-dest': 'document',
            'sec-fetch-mode': 'navigate',
            'sec-fetch-site': 'same-origin',
            'sec-fetch-user': '?1',
            'sec-gpc': '1',
            'upgrade-insecure-requests': '1',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36',
        }

        params = {
            'ref': 'characteristic',
            '10coreshape_s[1][]': [
                '10',
                '20',
            ],
            '10coreshape_s[2][]': [
                '30',
                '40',
                '50',
                '60',
                '70',
                '80',
                '90',
                '110',
                '120',
            ],
            '10coreshape_s[3][]': [
                '150',
                '160',
                '180',
                '190',
                '200',
                '210',
                '220',
                '230',
                '240',
                '250',
                '260',
                '270',
            ],
            '10coreshape_s[4][]': '280',
            '10coreshape_s[5][]': [
                '290',
                '300',
                '310',
            ],
            '10coreshape_s[6][]': [
                '320',
                '330',
                '340',
            ],
            '10coreshape_s[7][]': [
                '350',
                '370',
            ],
            '10coreshape_s[8][]': '380',
            '10coreshape_s[9][]': '500',
            'psts[]': '0',
            '_l': '100',
            '_p': offset,
            '_c': '10softmtrl-10softmtrl',
            '_d': '0',
        }

        response = requests.get(
            'https://product.tdk.com/pdc_api/en/search/ferrite/ferrite/ferrite-core/list.csv',
            params=params,
            cookies=cookies,
            headers=headers,
        )


        data = pandas.read_csv(io.StringIO(response.text), sep=",")
        data = data.where(pandas.notnull(data), None)
        data = data.rename(columns={"Size (IEC62317) A x B x C": "shape",
                                    "Material Name": "material",
                                    "Air Gap / mm": "gap_length",
                                    "Core Shape": "family"})

        print(data.columns)
        number_products = len(data.index)
        for index, row in data.iterrows():
            self.process(row)
        print(number_products)

        total_pages = offset + 2 if number_products != 0 else offset
        return {
            'current_offset' : offset + 1,
            'total': total_pages
        }

    def process(self, data):
        constants = PyMKF.get_constants()
        not_included_materials = ["HF60", "K1", "K10", "M33", "N22", "N45", "N48", "PC40", "PC90", "T", "PE22", "PC95", "PC50", "PC47"]
        if data['shape'] is None:
            return

        family = self.try_get_family(data['family'].split('(')[1].split(')')[0])
        if family in ["Toroid", "R"]:
            family = "T"
            data['shape'] = data['shape'].replace('R', 'T')
        if family == "PS":
            family = "P"
            data['shape'] = data['shape'].replace('PS', 'P')
        if family == "EE":
            family = "E"
            data['shape'] = data['shape'].replace('EE', 'E')
        shape = self.try_get_shape(data['shape'])
        material = self.try_get_material(data['material'])
        gap_length = float(data['gap_length']) / 1000

        if gap_length is None or math.isnan(gap_length):
            return

        if shape is None:
            family, shape = self.try_get_family_and_shape_statistically(data['shape'])
        if family in ["Toroid", "R"]:
            family = "T"
            data['shape'] = data['shape'].replace('R', 'T')

        if material is None:
            for not_included_material in not_included_materials:
                if not_included_material in data['material']:
                    return None
            assert 0, f"Unknown material in {data['material']}"

        print(family)
        print(shape)
        print(material)
        print(gap_length)
        if family is None:
            return 
        core_data = self.verify_shape_is_in_database(family, shape, material, data)

        status = 'production'
        core_manufacturer = "TDK"
        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        manufacturer_part_number = data['Part No.']
        datasheet = data['Catalog / Data Sheet']
        product_url = None
        distributor_reference = None
        quantity = None
        cost = None
        if gap_length == 0:
            gapping = []
        else:
            gapping = [
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

        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)


    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")

        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()


class MagneticsInventory(Stocker):

    def __init__(self):
        self.core_data = {}
        self.core_effective_parameters = {}
        self.unfound_shapes = {}

    def read_products(self, offset: int = 0):
        cookies = {
            'CMSPreferredCulture': 'en-US',
            'CMSCsrfCookie': 'ksCV/blOkR2v3g2hxK2q0HnHE3iCcewpoBwxoc6u',
            'ASP.NET_SessionId': 'qindzio2lmucl2nvd4ddh4n2',
        }

        headers = {
            'authority': 'www.mag-inc.com',
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
            'accept-language': 'es-ES,es;q=0.9,en;q=0.8',
            'cache-control': 'max-age=0',
            'content-type': 'application/x-www-form-urlencoded',
            # 'cookie': 'CMSPreferredCulture=en-US; CMSCsrfCookie=ksCV/blOkR2v3g2hxK2q0HnHE3iCcewpoBwxoc6u; ASP.NET_SessionId=qindzio2lmucl2nvd4ddh4n2',
            'dnt': '1',
            'origin': 'https://www.mag-inc.com',
            'referer': 'https://www.mag-inc.com/Advanced-Part-Number-Finder',
            'sec-ch-ua': '"Chromium";v="110", "Not A(Brand";v="24"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'sec-fetch-dest': 'document',
            'sec-fetch-mode': 'navigate',
            'sec-fetch-site': 'same-origin',
            'sec-fetch-user': '?1',
            'sec-gpc': '1',
            'upgrade-insecure-requests': '1',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36',
        }

        data = {
            '__CMSCsrfToken': 'soOkB5Je+r3RI4NuIinWi7kh3mf3SbwMNwBqa45koDsF67D6ftcoJCK1688JoAGXM5Ov5ldoKQr/ig/5Sx76jHtMghi4M74drCnibFdafSs=',
            '__EVENTTARGET': '',
            '__EVENTARGUMENT': '',
            '__VIEWSTATE': 'UAqlrfsuukO82cxxFkKAin3rMc0Jrd5rURuuzzPxoqrkKlCAUEwbv6y+xs1jXgGv8V5uttqmNOangY7B+f9aWvBm3Q0uNgnAx33LaKsbYFxSJuB/oyg3HrMVRpB+Oh3S9iTcqbIEjQvdxbg29HxLycCbLngIiQ+HkKHiZgs+S/HlscPQ9JYKiBwbY2CYB0zX6ggiLvpviF1URz77rKGtTb64yr2wqUemSQ19Pr5rHrlH95m0cwPfVSBGdRXt0jR8/+dO6586HeNwTe8fv2rQ5a9pnzW7rxYPIPq9r/G5pHfheXdfgqV8KVvAGe2C1+bgzJ8Xj/OryU6eOPTEWRtc6namEFYUaSRAsn99ZdNsjDEM4z3XzRswhHEVZl8ggB7+dRODkyMzpC3bMsPKmW1lo7H2QP57qYbNlHXVf1sX4Ir6UVKxww9WIIOcEyKtSL4pv/3U98OYjI+O4N573a3t934eIju/8DG+XhV9iN/EPxCHD/iMM7cgdaqvOZ2ZEtzt+LPKchIq2tmQf5rJgiXLXfZxu5jIUQcEUMVwZj8HylSx4jB8lbVdcHCDDFGTRXtM6Fehb0eUCSF9HlA8mKWN9Au4ZGmBuKKyQD2Pl1YLtwnY857gnxIn01IHID+CCACeo7hMImrovTAOETkaHLHPmrkIFOqPHnvXiLBNkEo0dZLQDx018rXx0OaoXcooynlIDLbbHj+dpT9RX4YGhM6wOHIiOZw23cXKAupd7mDJ0Emzl/YmjerWSScPuWQ1inqoC5t/I/DoiBTnEJQGX6aDqpIZBxWcAkiovwxixiroZgtX6bHV0+x4u/q9DaRt168W6a7cpd81Z1FTLUmFiwrDwvb4E1n3qaB4b4Gf0et7e1jJ0M/styw5aMI1my5bStxLgiH6GhCN7ypLsYVOZyg6lcMqlIHhZ8R0H9V9PL0udE1Bcxwh3bi4c/IX0B4TmQRTKJiIEbbTS/5FjA+bD40zDEY0Vprr1G+M+Sl0vf1l8WYlJGfhgHi26oMFLjm5a6v6tdCm5mPRYdc1P4xjkzRbokWLxw7Ex4EKXyR0nblYJbCkRb8LgdbdDE47VBNtuWSw750dbT04mxrOCtRK7vwCgh11FKyhESkd1ZT9ifmbb6gD6Ed7DNXtn+jGDI4Qup8bmGKd1z8zW2+pMDV+gPYfP/mZ6LFjAk2utxHnHMR1t2tZddKbpx9kV7s/lNotYRtWnz1fHiy4TJ/p83eZ1vkhxskIlKF3KzPldz6Geo/0Xk2vpcMHr+SZixMo1wTn/qMecrNgKCBAAG0fL6aqzDpqfMjz4ReoR2dkpoaql72CJquuqx+BIUYBhQMWPdJ2SWiH0eeQufcPISMqVFrm7s5KLD5/mlMgKd38jkmGoSIYYPjbzHa6kY3HPck//tej3/GVQlO2KYdbG2eFq2A8QXWdmQxsvUhC3dVXJmi9bgk0V5MsfPhX2ybpUH/S++5yJkIuS085zYYYQTSyR26xSeDir8WUuF/xlKPJrAHQWGoB30xLGu5gZ3aGzdjk8onoSKDdDn6t+P6wgbJKRQ2gjJIpD3qQAtOi2n7QNPvLkNwpMG9n5tOEkzh2F9zPQEkCv+v04KhshGNS6vdzdEk1HAP7R91O0lgsgY//40cYPL7s5t8kZqrdxpvKC2PyaEJR6CZlTcMy+eB27XVq5LSsW5eC5imqKr53dKS9jg+s2HPPqtfKm0FRO+KeWpua0Z/2klnhw0GDU1zDsxB3MHth9LLg1uST915CJQ7kvlCg09dyNiJrh+Jb+6gSLMBdSTqEZubmNcdiSCI3ZzhQo6RKwjmlvhcaxPSeEbFr/3rSs8OM8qTHlpOV2QDwoB/lMvlDiTEzDAk7eTSCf/GENy6GOn/dfMip3oMj/2GQZtMZCR99ooSaHtmO8fFLD3Xmn2emj8fZIdzJsesE9awdVDb2aOWjK+0AEZ1XqkdK3MWpbQui6cWGbi35bzt6loDbyiO4B2Y5+RjEO1dArnPaP6nAud8cwWGN84bpfJXoKpkLMVEZucSGy1jiO6unkD0R9iIijQ3mXIuSPlb9OA6AzgqpoiW2rS6CvR1h0GYUjLoVD/0yJteX1XrXSWuT37Y5NnnSTsPi/U9O9Ro8uEIudGWCJCLU43n7hQkqwAewLD27ZJTXC0OUm2ZQHf0aDohq7jAJHuxxr3xz0gMnRDtrrG7C4B3sM3o/OKQJwDrTOwKBv5oCzY95P4P4O4XrcrPo10jRk+AqOzUykrRO5v4tnjmky9pVLF7URy9lO+q0VCwJqGktpfKDzPKGZcghGGraRoJUIzpaULc7hIm+LZtVYNUNp3ggo+Q5whNQltn9kWezAcMURyEFEe3GpBQfevlcV1bd79H01T9YkQ/lp7ic9sfZ8C62qtzCCcYUDmu8R3MBtgeZgCltaba4qEQteM/cUukMKYc5Mzho1X3RDZrVTgxfDR36hphVVfFYjCTF7K343F8TSidLM9gD1dx9YkIeFDbH3puLOR2XZKV6eEvLCNgMKrTnPCt0hTEO+0mBhaZLC4+YY/lClZ07Tmp6mZxeOEgFl1bdaSXpn86jvonpVzpTB+XH8/Og9Kj7dTxTxpKbEry1MvN958W/rsIiTALCe96SR5OPpP+grGHkV5mOzcFJSSkAmbChUUs1F1pfrxJoCm35mVUdjwEJT/Z+xwD7csOiArr3JND27KgPES2JG9tA04tghMPK3udTFHeQjVLWW8gIRfRrDqEdqzjWHXuDRFISR1tpx6NXKRNEwUTsEhE+Qt8JQp89ElAOr6ReEjD2HyBAbKX8wEh4XUQ351qa927aPO2s5kdKUDpNA9P9YvHRhLAqZGJKsWsmvdgUpgkdQvVQZCrOL8A0z+08qoouZ7PRqP1ZhhZFpqkdFtYarpmNZ+30QDbmKPvrefZSrHFkouVggh8pzNQFC/zVWxNBBhwmxK8tPtgdW0nxSfpXB4DHmDuF7mLs5sOcBkluOcs/F9wO/TwqtS1s2Vm2tYjR2oTtOLD0zYHxJT8ta1BSMmblsQP7Bwa8G0xLTPw8tHgMECoB1156qcPIofSTAGMWT9peGIEE62PPDA+Nz033oOPn+BOcXFWSKlUEVmw7vha3trzYkTGMG0Lbngw+TQu30YuXH0WmSxFd9jRjGkDTjyC28DqEjRqkz5lOMSeQOGmJ9E9aDNK6jTngluPVqFNx+yN0MRST8cHjzWVkNmiGrdTf86sIZyHVr2H4Mj6T+yOtzHvfrRR/DvsVKhK2lpfOSCZkUiCo/BAAE9RITuBfYQOWOTFGFDpDQeZxRTZFc8yu71cYEu/NVH+QyhZfrjI0r4Z4/QJr0c+9uPEwmPBfN3ZuREhqQl/xByBgIIQKQVinWS2qRkKg1rqsICOFee+DgaaONmSrpUdw1EC0kAjIQ+rZC8v4cA==',
            'lng': 'en-US',
            '__VIEWSTATEGENERATOR': 'A5343185',
            '__SCROLLPOSITIONX': '0',
            '__SCROLLPOSITIONY': '340',
            'p$lt$ctl03$SmartSearchBox$txtWord': '',
            'p$lt$ctl04$MyMagneticsLoginForm$Username': '',
            'p$lt$ctl04$MyMagneticsLoginForm$Password': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PartNumber': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresPerm': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresMaterial': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresShape': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresALInductance': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresODLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresODLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresODLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresIDLegLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresIDLegLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresIDLegLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresHTHeight': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresHTHeightUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresHTHeightPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresLePathLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresLePathLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresLePathLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresAeCrossSection': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresAeCrossSectionUnits': 'mm2',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresAeCrossSectionPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresVeVolume': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresVeVolumeUnits': 'mm3',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresVeVolumePlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$PowderCoresSubmit': 'Search',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerritePerm': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteMaterial': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteShape': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteALInductance': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteODLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteODLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteODLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteIDLegLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteIDLegLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteIDLegLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteHTWidth': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteHTWidthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteHTWidthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteLePathLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteLePathLengthUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteLePathLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteAeCrossSection': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteAeCrossSectionUnits': 'mm2',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteAeCrossSectionPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteVeVolume': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteVeVolumeUnits': 'mm3',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteVeVolumePlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$FerriteWaAc': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreMaterial': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreThickness': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreCaseType': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareOD': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareODUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareODPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareID': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareIDUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareIDPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareHT': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareHTUnits': 'mm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreBareHTPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreLePathLength': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreLePathLengthUnits': 'cm',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreLePathLengthPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreAeCrossSection': '',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreAeCrossSectionUnits': 'cm2',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreAeCrossSectionPlusMinus': '5',
            'p$lt$ctl06$pageplaceholder$p$lt$ctl01$EditableContent$ucEditableText$widget1$ctl00$StripWoundCoreWaAc': '',
        }

        response = requests.post('https://www.mag-inc.com/Advanced-Part-Number-Finder', cookies=cookies, headers=headers, data=data)
        data = pandas.read_html(io.StringIO(response.text), extract_links='body')[0]

        # data = self.clean_html(response.text)
        # pprint.pprint(data)
        for index, row in data.iterrows():
            self.process(row)
        pprint.pprint(self.unfound_shapes)
        pprint.pprint(len(self.unfound_shapes))

        return {
            'current_offset' : 0,
            'total': 0
        }

    def process(self, data):
        not_included_materials = ["75-Series 125",
                                  "AmoFlux",
                                  "High Flux 19",
                                  "High Flux 90",
                                  "Kool Mu 19",
                                  "MPP 250",
                                  "MPP 367",
                                  "XFlux 14",
                                  "XFlux 50"]
        data = data.to_dict()
        for key, value in data.items():
            if key == 'Data Sheet Download':
                if value[1] is None:
                    data[key] = None
                else:
                    data[key] = f"https://www.mag-inc.com{value[1]}"
            else:
                data[key] = value[0]

        data['Material'] = data['Material'].replace("75 Series", "75-Series")
        data['Material'] = data['Material'].replace("Kool Mu HF", "Kool Mu Hf")
        data['Material'] = data['Material'].replace("KOOL MU MAX", "Kool Mu MAX")
        data['Material'] = data['Material'].replace("Xflux", "XFlux")

        family = self.try_get_family(data['Shape'].split(' Core')[0].replace("Toroid", "T").replace("Thinz", "T"))
        material = self.try_get_material(data['Material'] + " " + data['Perm'])

        if family in ['Block', 'Segment']:
            return
        elif family == 'T':
            temptative_shape = f"{family} {self.adaptive_round(float(data['OD / Length mm']))}/{self.adaptive_round(float(data['ID / Leg Length mm']))}/{self.adaptive_round(float(data['HT / Width mm']))}"
        else:
            temptative_shape = f"{family} {self.adaptive_round(float(data['OD / Length mm']))}/{self.adaptive_round(float(data['ID / Leg Length mm']))}/{self.adaptive_round(float(data['HT / Width mm']))}"
        try:
            shape = self.find_shape_closest_dimensions(family, temptative_shape)
        except ZeroDivisionError:
            return 
        pprint.pprint(family)
        pprint.pprint(material)
        pprint.pprint(shape)
        if shape is None:
            if temptative_shape not in self.unfound_shapes:
                self.unfound_shapes[temptative_shape] = (f"{family} {round(float(data['OD / Length mm']) / 1000, 5):.5f}/{round(float(data['ID / Leg Length mm']) / 1000, 5):.5f}/{round(float(data['HT / Width mm']) / 1000, 5):.5f}", data['Data Sheet Download'], data['Part Number'])

        if material is None:
            for not_included_material in not_included_materials:
                if not_included_material in f"{data['Material']} {data['Perm']}":
                    return None
            # assert 0, f"Unknown material in {data['Material']} {data['Perm']}"
            return

        if shape is None:
            return
        print(family)
        print(shape)
        print(material)
        core_data = self.verify_shape_is_in_database(family, shape, material, data)

        status = 'production'
        core_manufacturer = "Magnetics"
        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        manufacturer_part_number = data['Part Number']
        datasheet = data['Data Sheet Download']
        product_url = None
        distributor_reference = None
        quantity = None
        cost = None
        gapping = []

        self.add_core(family=family,
                      shape=shape,
                      material=material,
                      core_type=core_type,
                      core_manufacturer=core_manufacturer,
                      status=status,
                      manufacturer_part_number=manufacturer_part_number,
                      product_url=product_url,
                      distributor_reference=distributor_reference,
                      quantity=quantity,
                      cost=cost,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)


    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            pprint.pprint('==================================================================================================')
            current_status = self.read_products(current_offset)
            current_offset = current_status['current_offset']
            remaining = current_status['total'] - current_status['current_offset']
            pprint.pprint(f"remaining: {remaining}")
            pprint.pprint(f"self.core_data: {len(self.core_data.values())}")

        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file)
        out_file.close()




if __name__ == '__main__':  # pragma: no cover
    ferroxcube_inventory = FerroxcubeInventory()
    ferroxcube_inventory.remove_current_inventory()
    ferroxcube_inventory.get_cores_inventory()

    tdk_inventory = TdkInventory()
    tdk_inventory.get_cores_inventory()

    magnetics_inventory = MagneticsInventory()
    magnetics_inventory.get_cores_inventory()

    digikeyStocker = DigikeyStocker()
    digikeyStocker.get_cores_stock()

    mouserStocker = MouserStocker()
    mouserStocker.get_cores_stock()

    gatewayStocker = GatewayStocker()
    gatewayStocker.get_cores_stock()

    # pprint.pprint(mouserStocker.unfound_descriptions)
