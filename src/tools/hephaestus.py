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
from datetime import date
from difflib import SequenceMatcher


class Stocker():
    def __init__(self):
        self.core_data = {}

    def remove_current_inventory(self):
        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            os.remove(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson")

    def verify_shape_is_in_database(self, family, shape, material, product):
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
            core_data = PyMKF.calculate_core_data(core_data, False)
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
        core_data = PyMKF.get_core_shape_names(True)

        dimensions = shape.split(" ")[1].split("/")

        print("dimensions")
        print(dimensions)
        try:
            for shape_name in core_data:
                if 'X' in shape_name:
                    continue
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
        shape_names = PyMKF.get_core_shape_names(True)
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

                    core_data = PyMKF.calculate_core_data(core, False)

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

    def add_core(self, family, shape, material, core_type, core_manufacturer, status, manufacturer_part_number, product_url, distributor_reference, quantity, cost, coating, gapping=[], datasheet=None, unique=False, distributor_name=None):
        core_name = self.get_core_name(shape, material, gapping, coating)

        print(f"Trying to add {core_name}")

        if core_name not in self.core_data:
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
            # if isinstance(shape, str):
            #     shape = PyMKF.calculate_core_data(core, False)['functionalDescription']['shape']
            #     for dimension_key, dimension in shape['dimensions'].items():
            #         new_dimension = {}
            #         for key, value in dimension.items():
            #             if value != None:
            #                 new_dimension[key] = value
            #         shape['dimensions'][dimension_key] = new_dimension
            #     core['functionalDescription']['shape'] = shape
            if not isinstance(shape, str):
                pprint.pprint(core)
                assert 0
                

            print(f"Adding {core_name}")
            if coating is not None:
                core['functionalDescription']['coating'] = coating
            self.core_data[core_name] = core
        elif unique:
            assert 0, f"Already inserted core: {core_name}"

        if len(gapping) > 0:
            # Just in case, adding core ungapped too
            if not isinstance(shape, str):
                shape = shape["name"]
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
                          gapping=[],
                          coating=coating,
                          datasheet=datasheet,
                          unique=False)

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
        core_data = PyMKF.calculate_core_data(core_data, False)
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
        for shape_name in PyMKF.get_core_shape_names(True):
            error = abs(len(shape_name) - len(text)) / len(text)
            if shape_name in text and error < 0.3:
                if error < current_error:
                    current_error = error
                    current_shape = shape_name

        if current_shape is not None:
            return current_shape

        fixed_text = fix(text)
        current_error = 1
        current_shape = None
        for shape_name in PyMKF.get_core_shape_names(True):
            if (len(fixed_text) == 0):
                return None
            error = abs(len(fix(shape_name)) - len(fixed_text)) / len(fixed_text)
            if fix(shape_name) in fixed_text and error < 0.3:
                if error < current_error:
                    current_error = error
                    current_shape = shape_name
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
            return text.replace(" 0", " ").replace(" ", "").replace(",", ".").replace("Mµ", "Mu").replace("Hƒ", "Hf").upper()

        materials = PyMKF.get_core_material_names_by_manufacturer(manufacturer)

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

        for shape_name in PyMKF.get_core_shape_names(True):
            s = SequenceMatcher(None, fix(shape_name), fixed_text)
            if s.ratio() > max_ratio:
                max_ratio = s.ratio()
                best_fit = shape_name

        print("max_ratio")
        print(max_ratio)
        print("best_fit")
        print(best_fit)
        if max_ratio < 0.9:
            return None, None

        family = best_fit.split(" ")[0]

        return family, best_fit


class FerroxcubeInventory(Stocker):
    def read_concentric_products(self, offset: int = 0):
        cookies = {}

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
            timeout=2
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
            timeout=2
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
            print(f"core_name: {core_name}")
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

            # Some cores have typos in the series value, we correct them with the name
            for typo in typos:
                core_name.replace(typo['typo'], typo['correction'])

            family = self.try_get_family(core_name)
            shape = self.try_get_shape(core_name)
            material = self.try_get_material(core_name, "Ferroxcube")
            if family == 'TX':
                family = 'T'


            if shape is None:
                family, shape = self.try_get_family_and_shape_statistically(core_name)

            print(f"family: {family}")
            print(f"shape: {shape}")
            print(f"material: {material}")
            if material is None:
                for not_included_material in not_included_materials:
                    if not_included_material in core_name:
                        break
                        # return None
                # assert 0, f"Unknown material in {core_name}"

            if material is None:
                continue
            if family is None:
                continue
                # assert 0
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
            coating = self.try_get_core_coating(core_name)

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

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson") as f:
                previous_data = ndjson.load(f)
                for row in previous_data:
                    self.core_data[row['name']] = row

        while remaining > 0:
            try:
                pprint.pprint('==================================================================================================')
                current_status = self.read_concentric_products(current_offset)
                current_offset = current_status['current_offset']
                remaining = current_status['total'] - current_status['current_offset']
                pprint.pprint(f"remaining: {remaining}")
                pprint.pprint(f"self.core_data: {len(self.core_data.values())}")
            except (requests.exceptions.ConnectTimeout, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
                continue

        remaining = 1
        current_offset = 1
        while remaining > 0:
            try:
                pprint.pprint('==================================================================================================')
                current_status = self.read_toroidal_products(current_offset)
                current_offset = current_status['current_offset']
                remaining = current_status['total'] - current_status['current_offset']
                pprint.pprint(f"remaining: {remaining}")
                pprint.pprint(f"self.core_data: {len(self.core_data.values())}")
            except (requests.exceptions.ConnectTimeout, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
                continue

        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file, ensure_ascii=False)
        out_file.close()


class TdkInventory(Stocker):
    def read_products(self, offset: int = 0):
        cookies = {}

        headers = {
            'accept': '*/*',
            'accept-language': 'es-ES,es;q=0.9',
            'content-type': 'application/x-www-form-urlencoded; charset=UTF-8',
            # 'cookie': 'ak_bmsc=7CFA1415F1B17AD986FC6CADA8E5E5E1~000000000000000000000000000000~YAAQXLsUArHr9NqSAQAARNFq8Rk1zMKU5kd3TfdGcYLP3WHh8bcIgLdpnqmSIiVL3Bw72gKPNM4p+iiXT6hExsmGk2XEw/+pHIrX/16ufI9fI1yw0dUOR1bzGYZCTpYgFRPgL2YA0K4K38D+T9GV7ndEad3ZJkdWwTE6YZFDU0nbtoBnJHXJ4TH02Cbad5R0jX3HsMkpNKohfkTHzb2qQbLY9T84Uh3NJDFlbnMrLkgqodNt3xgDH3yTOGeWV96ofgWhxw6cV+E4MezJK4U+ER3S2sfl/83dSvA57QGl4VZIdzl42DvlPIQv+P3KKJBHiObaWM4xym44eH0hPxQ15qyBGqiP/hPjPDjmYB/nP4ZQa37oyDf62ANN65iLvZjrAXoI9xIEG6s=; wooTracker=j3nP9mB8lSGt; _abck=7738C6273CFB64C2B06FF20E4388A16E~0~YAAQXLsUAmjs9NqSAQAAzehq8QzU7lHLX0vTQWpZK2IDk6RAr0zRep62BFnJmWqmIlr2SJ456Fjkd1VWM/fCkJcWE3cczMscQzsw+PrdCVuBc7LawkloM3yzjOwBZn+RbiyYZy44Lr2C3XYqcNCZ4rZ1hukFT3PiUUEocgn5AZ20f3KQg59nxRJD8elDsZKlRS1ybzsuXJ5zOGj+F7XR3pix6e9QyrzmS9zjsZhc6ZDtNgUd7GYvCiLvS5rom4lMFDC/PaGtiovO4L4DCXzUQzNC5ok47dud1srO7liXXk2H3xMbdyWjIULt0+4HQpszArJmPejlmYVMOJdP7LyUIZ4Gaqtu+GXsjLf7ANpHp2KK6NOOerbRSW92yTvo9NFKbBr+3WmUllOCIB5FbA80JyhGZjgzrgyxTcpaZlq/Baq2Ee7Zh3zeLU4EH7ZJBwWUTYTkN30U~-1~-1~-1; ps_ls=0; bm_sz=EAAD71EF804931D7F9B1C5E2EF3A5A14~YAAQXLsUAlf09NqSAQAAI/9r8RlbRu2fzx8H2xuMsUBYW24D1aJ8n+hzgRd+u14jVjp5iL5nADNbgz1kG8emIt18QIn9gAiRNz6IS6Ege40HTc0+VqltPFK3iPOaK7crktk713fbQu6ss8uA55Tt/HIhxfJfEF4juziHN5LYw7i/36YcnqMuS2ONYnKFYDZPtks29/Tt8kx4DK7M1vSy0tl46kxuXHOlolK8tuEZrErwAEPKNnrjP+rMmIPWbvf5EtE5kn6EuBbhQaR03TxUSCxKE0g4QibLqvr0Pxc1bPLwSqmgzWc6FdK5smDmsdi3brSky/mGDbVUXTmciT5PN0xDkIjrZMi+BapwEYSKfDDn+Y/Dgr71HdMlgSVXYhzUrT6xZ5wueJyEvUjCuchuUO5UK5E0YiaWeodF85ihpwTf/blF9kR/N+KjBl0e85PX75nJBzp65QY=~4339766~3421232; product_search_fw=eyJpdiI6InRDU3VJY0xrNGdURmdtYXU2MnBpNHc9PSIsInZhbHVlIjoiVHFROFBpeStnMmlobW44WU93OEU4c0pFVisxUXIvZ2JVZU1rRmV0eGNJNTYrdWpraHBHTzFHamx0RmloK3hFVllvV2c3ZzhZak1lMkNkNDU0clZmR0JaQlRiV0lyc2w2eTdjRUdpSFhibXE3dVczVDI3UE1NNVJGYjU2MThXSUoiLCJtYWMiOiI1NDMxMDY0M2MwNDVjMjg0ZjBmMDc1YjNhYzBhN2ExMjczOWRmNDZkZDBlNTk1OWJlYTJiMTMzMGMxNTI2OTFhIiwidGFnIjoiIn0%3D; bm_sv=0BFC2939F3C4F1993632316C3EC706FA~YAAQXLsUApcI9dqSAQAAQD1u8Rko5RQrLOBwuwhxDIkVNbsDWr0IGR6vpFqsIM0jNlGiOcrtytlcrRxqRKJEt70j2wr8rQCevivajFn/7qjaH2plTpNX2KkG+8/jNJ6LyJS4YFM8YvPPQl9LBQ1/2wVIywHKwV5PBztuKtV4CQAIu9LclQXFIORtRmYnHzF9q7QBpeggA16/XY9MWFo2hK47cUbEULsiPZb/K+QgkF3WikJk8N8HqF+/T/OcFA==~1',
            'origin': 'https://product.tdk.com',
            'priority': 'u=1, i',
            'referer': 'https://product.tdk.com/en/search/ferrite/ferrite/ferrite-core/list',
            'sec-ch-ua': '"Chromium";v="128", "Not;A=Brand";v="24", "Google Chrome";v="128"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'sec-fetch-dest': 'empty',
            'sec-fetch-mode': 'cors',
            'sec-fetch-site': 'same-origin',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36',
            'x-requested-with': 'XMLHttpRequest',
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
            'current_offset': offset + 1,
            'total': total_pages
        }

    def process(self, data):
        constants = PyMKF.get_constants()
        # not_included_materials = ["HF60", "K1", "K10", "M33", "N22", "N45", "N48", "PC40", "PC90", "T", "PE22", "PC95", "PC50", "PC47"]
        not_included_materials = ["HF60", "PC40", "PC90", "T", "PE22", "PC50"]
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
        material = self.try_get_material(data['material'], 'TDK')
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

        if family != 'T':
            coating = None
        else:
            if data['Part No.'][6] == 'A':
                coating = None
            elif data['Part No.'][6] == 'L':
                coating = 'epoxy'
            elif data['Part No.'][6] == 'P':
                coating = 'parylene'
            else:
                assert 0, data['Part No.']


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
                      coating=coating,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)

    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson") as f:
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
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file, ensure_ascii=False)
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
        material = self.try_get_material(data['Material'] + " " + data['Perm'], 'Magnetics')

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

        if family != 'T':
            coating = None
        else:
            coating = 'epoxy'

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
                      coating=coating,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)

    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson") as f:
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
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file, ensure_ascii=False)
        out_file.close()


class FairRiteInventory(Stocker):
    def __init__(self):
        self.core_data = {}
        self.unfound_shapes = []

    def read_products(self, offset: int = 0):
        data = pandas.DataFrame()

        urls = [
            'https://fair-rite.com/product-category/inductive-components/e-cores/',
            'https://fair-rite.com/product-category/inductive-components/pot-cores/',
            'https://fair-rite.com/product-category/inductive-components/efd-cores/',
            'https://fair-rite.com/product-category/inductive-components/toroids/',
            'https://fair-rite.com/product-category/inductive-components/etd-cores/',
            'https://fair-rite.com/product-category/inductive-components/eer-cores/',
            'https://fair-rite.com/product-category/inductive-components/ep-cores/',
            'https://fair-rite.com/product-category/inductive-components/planar-cores/',
            'https://fair-rite.com/product-category/inductive-components/pq-cores/',
            'https://fair-rite.com/product-category/inductive-components/rm-cores/',
        ]

        for url in urls:
            response = requests.get(url)


            aux = pandas.read_html(response.text, extract_links='body')[1]
            columns = [x[0] for x in aux.columns]
            aux.columns = columns
            data = pandas.concat([data, aux])

        data = data.rename(columns={"Generic Size": "shape",
                                    "Part Number": "manufacturer_part_number"})
        data['link'] = data.apply(lambda row: row['manufacturer_part_number'][1], axis=1)
        data['manufacturer_part_number'] = data.apply(lambda row: row['manufacturer_part_number'][0], axis=1)
        data['shape'] = data.apply(lambda row: row['shape'][0] if isinstance(row['shape'], tuple) else row['shape'], axis=1)
        data['A'] = data.apply(lambda row: row['A'][0].split('±')[0].split('Min')[0].split('Max')[0].split('+')[0].split('-')[0].lstrip('0') if isinstance(row['A'], tuple) else row['A'], axis=1)
        data['B'] = data.apply(lambda row: row['B'][0].split('±')[0].split('Min')[0].split('Max')[0].split('+')[0].split('-')[0].lstrip('0') if isinstance(row['B'], tuple) else row['B'], axis=1)
        data['C'] = data.apply(lambda row: row['C'][0].split('±')[0].split('Min')[0].split('Max')[0].split('+')[0].split('-')[0].lstrip('0') if isinstance(row['C'], tuple) else row['C'], axis=1)
        pandas.set_option('display.max_rows', None)
        data = data.where(pandas.notnull(data), None)

        for index, row in data.iterrows():
            self.process(row)
        number_products = len(data.index)
        print(number_products)

        return {
            'current_offset' : number_products,
            'total': number_products
        }

    def try_get_material(self, text):
        material_name = str(text)[2:4]
        return material_name

    def process(self, data):
        constants = PyMKF.get_constants()
        not_included_materials = []
        if data['shape'] is None:
            family = 'T'
            temptative_shape = f"{family} {self.adaptive_round(float(data['A']))}/{self.adaptive_round(float(data['B']))}/{self.adaptive_round(float(data['C']))}"
            print(temptative_shape)
            shape = self.find_shape_closest_dimensions(family, temptative_shape, limit=0.11)

            if shape is None:
                new_shape = '{"magneticCircuit": "closed", "type": "standard", "family": "t", "aliases": [], "name": "' + temptative_shape + '",                                          "dimensions": {"A": {"nominal": '
                new_shape += f"{float(data['A']) / 1000:.5f}"
                new_shape += '}, "B": {"nominal": '
                new_shape += f"{float(data['B']) / 1000:.5f}"
                new_shape += '}, "C": {"nominal": '
                new_shape += f"{float(data['C']) / 1000:.5f}"
                new_shape += '}}}   '
                self.unfound_shapes.append(new_shape)
                return
        else:

            if data['shape'].startswith('I') or data['shape'].startswith('EQI') or data['shape'].startswith('EI') or data['shape'].startswith('UI'):
                return
            print(f"data['shape']: {data['shape']}")

            family = self.try_get_family(data['shape'])
            if family in ["Toroid", "R"]:
                family = "T"
                data['shape'] = data['shape'].replace('R', 'T')
            if family == "PS":
                family = "P"
                data['shape'] = data['shape'].replace('PS', 'P')
            if family == "EE":
                family = "E"
                data['shape'] = data['shape'].replace('EE', 'E')
            if family == "EEQ":
                family = "EQ"
                data['shape'] = data['shape'].replace('EEQ', 'EQ')
            shape = self.try_get_shape(data['shape'])

            # if shape is None:
            #     family, shape = self.try_get_family_and_shape_statistically(data['shape'])
            if family in ["Toroid", "R"]:
                family = "T"
                data['shape'] = data['shape'].replace('R', 'T')

            if family is None:
                family = self.try_get_family(data['shape'].split('(')[0])
            if shape is None:
                shape = self.try_get_shape(data['shape'].split('(')[0])
        gap_length = 0
        material = self.try_get_material(data['manufacturer_part_number'])
        # if material is None:
        #     for not_included_material in not_included_materials:
        #         if not_included_material in data['material']:
        #             return None
        #     assert 0, f"Unknown material in {data['material']}"

        if material in not_included_materials:
            return

        print(f"family: {family}")
        print(f"shape: {shape}")
        print(f"material: {material}")

        if family is None:
            print(data['link'])
            assert 0

        if shape is None:
            print(data['link'])
            assert 0
        core_data = self.verify_shape_is_in_database(family, shape, material, data)

        status = 'production'
        core_manufacturer = "Fair-Rite"
        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        manufacturer_part_number = data['manufacturer_part_number']
        datasheet = data['link']
        product_url = None
        distributor_reference = None
        quantity = None
        cost = None
        coating = None
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
                      coating=coating,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)

    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson") as f:
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

        print(self.unfound_shapes)
        core_data = pandas.DataFrame(self.core_data.values())
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file, ensure_ascii=False)
        out_file.close()


class MicrometalsInventory(Stocker):

    def __init__(self):
        self.core_data = {}
        self.unfound_shapes = {}

    def read_products(self, offset: int = 0):
        response = requests.get("https://www.micrometals.com/api/v1/shapes/vw_plweb_t/products/?units=in&offset=10&limit=")
        t_data = json.loads(response.text)['results']
        print(t_data[0])
        response = requests.get("https://www.micrometals.com/api/v1/shapes/vw_plweb_e/products/?units=in&offset=10&limit=")
        e_data = json.loads(response.text)['results']
        print(e_data[0])
        response = requests.get("https://www.micrometals.com/api/v1/shapes/vw_plweb_pc/products/?units=in&offset=10&limit=")
        p_data = json.loads(response.text)['results']
        print(p_data[0])
        response = requests.get("https://www.micrometals.com/api/v1/shapes/vw_plweb_pq/products/?units=in&offset=10&limit=")
        pq_data = json.loads(response.text)['results']
        print(pq_data[0])
        response = requests.get("https://www.micrometals.com/api/v1/shapes/vw_plweb_eq/products/?units=in&offset=10&limit=")
        eq_data = json.loads(response.text)['results']
        print(eq_data[0])

        for index, row in enumerate(t_data):
            self.process(row, 'T')
        for index, row in enumerate(e_data):
            self.process(row, 'E')
        for index, row in enumerate(p_data):
            self.process(row, 'P')
        for index, row in enumerate(eq_data):
            self.process(row, 'EQ')
        for index, row in enumerate(pq_data):
            self.process(row, 'PQ')
        pprint.pprint(self.unfound_shapes)
        pprint.pprint(len(self.unfound_shapes))

        return {
            'current_offset': 0,
            'total': 0
        }

    def process(self, data, family):
        pprint.pprint("***********************************")
        not_included_materials = ['Mix 0']

        partnumber = data['partnumber']
        if '-' in data['material']:
            data['material'] = data['material'].replace('-', "Mix ")
        else:
            data['material'] = f"{data['material']} {int(data['permeability'])}"
        material = self.try_get_material(data['material'], 'Micrometals')

        if family == 'T':
            temptative_shape = f"{family} {self.adaptive_round(float(data['od']) * 25.4)}/{self.adaptive_round(float(data['id']) * 25.4)}/{self.adaptive_round(float(data['ht']) * 25.4)}"
        elif family in ['E', 'EQ']:
            temptative_shape = f"{family} {self.adaptive_round(float(data['dima']) * 25.4)}/{self.adaptive_round(float(data['dimb']) * 25.4)}/{self.adaptive_round(float(data['dimc']) * 25.4)}"
        elif family in ['PQ']:
            temptative_shape = f"{family} {self.adaptive_round(float(data['dima']) * 25.4)}/{self.adaptive_round(float(data['dimc']) * 2 * 25.4)}"
        elif family in ['P']:
            temptative_shape = f"{family} {self.adaptive_round(float(data['dima']) * 25.4)}/{self.adaptive_round(float(data['dimc']) * 25.4)}"

        try:
            shape = self.find_shape_closest_dimensions(family, temptative_shape)
        except ZeroDivisionError:
            return 

        if shape is None:
            if temptative_shape not in self.unfound_shapes:
                if family == 'T':
                    t_data = {"magneticCircuit": "closed", "type": "standard", "family": "t", "aliases": [], "name": temptative_shape,
                              "dimensions": {
                                  "A": {"nominal": f"{round(float(data['od']) * 25.4 / 1000, 5):.5f}"},
                                  "B": {"nominal": f"{round(float(data['id']) * 25.4 / 1000, 5):.5f}"},
                                  "C": {"nominal": f"{round(float(data['ht']) * 25.4 / 1000, 5):.5f}"}}
                              }

                    self.unfound_shapes[temptative_shape] = (json.dumps(t_data), None, data['partnumber'])
                elif family in ['E']:
                    e_data = {"magneticCircuit": "open", "type": "standard", "family": "e", "aliases": [], "name": temptative_shape,
                              "dimensions": {
                                  "A": {"nominal": f"{round(float(data['dima']) * 25.4 / 1000, 5):.5f}"},
                                  "B": {"nominal": f"{round(float(data['dimb']) * 25.4 / 1000, 5):.5f}"},
                                  "C": {"nominal": f"{round(float(data['dimc']) * 25.4 / 1000, 5):.5f}"},
                                  "D": {"nominal": f"{round(float(data['dimd']) * 25.4 / 1000, 5):.5f}"},
                                  "E": {"nominal": f"{round(float(data['dime']) * 25.4 / 1000, 5):.5f}"},
                                  "F": {"nominal": f"{round(float(data['dimf']) * 25.4 / 1000, 5):.5f}"}}
                              }

                    self.unfound_shapes[temptative_shape] = (json.dumps(e_data), None, data['partnumber'])
                elif family in ['EQ']:
                    self.unfound_shapes[temptative_shape] = (f"{family} {round(float(data['dima']) * 25.4 / 1000, 5):.5f}/{round(float(data['dimb']) * 25.4 / 1000, 5):.5f}/{round(float(data['dimc']) * 25.4 / 1000, 5):.5f}", None, data['partnumber'])
                elif family in ['PQ']:
                    self.unfound_shapes[temptative_shape] = (f"{family} {round(float(data['dima']) * 25.4 / 1000, 5):.5f}/{round(float(data['dimc']) * 2 * 25.4 / 1000, 5):.5f}", None, data['partnumber'])
                elif family in ['P']:
                    self.unfound_shapes[temptative_shape] = (f"{family} {round(float(data['dima']) * 25.4 / 1000, 5):.5f}/{round(float(data['dimc']) * 2 * 25.4 / 1000, 5):.5f}", None, data['partnumber'])

        if material is None:
            for not_included_material in not_included_materials:
                if not_included_material in f"{data['material']}":
                    return None
            pprint.pprint(data)
            # assert 0, f"Unknown material in {data['material']}"
            return

        if shape is None:
            return

        core_data = self.verify_shape_is_in_database(family, shape, material, data)

        status = 'production'
        core_manufacturer = "Micrometals"
        core_type = 'toroidal' if family == 'T' else 'two-piece set'
        manufacturer_part_number = data['partnumber']
        datasheet = f"https://datasheets.micrometals.com/{data['partnumber']}-DataSheet.pdf"
        product_url = data['analyzer_url']
        distributor_reference = data['partnumber']
        quantity = data['qty']
        cost = None
        gapping = []

        if family != 'T':
            coating = None
        else:
            coating = 'parylene'

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
                      coating=coating,
                      gapping=gapping,
                      datasheet=datasheet,
                      unique=False)

    def get_cores_inventory(self):
        remaining = 1
        current_offset = 2

        if os.path.exists(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson"):
            with open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson") as f:
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
        out_file = open(f"{pathlib.Path(__file__).parent.resolve()}/cores_inventory.ndjson", "w")
        ndjson.dump(core_data.to_dict('records'), out_file, ensure_ascii=False)
        out_file.close()


if __name__ == '__main__':  # pragma: no cover
    # ferroxcube_inventory = FerroxcubeInventory()
    # ferroxcube_inventory.remove_current_inventory()
    # ferroxcube_inventory.get_cores_inventory()

    tdk_inventory = TdkInventory()
    tdk_inventory.get_cores_inventory()

    magnetics_inventory = MagneticsInventory()
    magnetics_inventory.get_cores_inventory()

    fair_rite_inventory = FairRiteInventory()
    fair_rite_inventory.get_cores_inventory()

    micrometals_inventory = MicrometalsInventory()
    micrometals_inventory.get_cores_inventory()
