import json
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


# Initial permeability
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_initial_permeability.csv"):
    initial_permeability_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_initial_permeability.csv")
else:

    materials = {
        "C350": 121,
        "C351": 122,
        "K1": 123,
        "K10": 124,
        "M33": 125,
        "M34": 150,
        "N22": 126,
        "N27": 127,
        "N30": 128,
        "N41": 129,
        "N45": 130,
        "N48": 131,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
        "T35": 142,
        "T36": 143,
        "T37": 144,
        "T38": 145,
        "T46": 146,
        "T57": 147,
        "T65": 148,
        "T66": 149,
    }
    initial_permeability_data = pandas.DataFrame()
    for material, material_id in materials.items():
        data = {
            'init_permeab[cbmaterial]': str(material_id),
            'init_permeab[OK]': '',
        }

        response = requests.post(
            'https://tools.tdk-electronics.tdk.com/mdt/index.php/initperm',
            data=data,
        )
        data = pandas.read_html(StringIO(response.text))[0]

        initial_permeability_datum = data.rename(columns={"T [°C]": "temperature", "µi [-]": "value"})
        initial_permeability_datum["material"] = [material] * len(initial_permeability_datum.index)
        initial_permeability_data = pandas.concat([initial_permeability_data, initial_permeability_datum], ignore_index=True)
    print(initial_permeability_data)
    initial_permeability_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_initial_permeability.csv", index=False)

# Core losses
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_volumetric_core_losses.csv"):
    volumetric_core_losses_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_volumetric_core_losses.csv")
else:
    materials = {
        "N27": 127,
        "N41": 129,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
    }
    volumetric_core_losses_data = pandas.DataFrame()
    frequencies = [25, 50, 100, 200, 300, 500, 700, 1000]
    magnetic_flux_densities = [13, 25, 50, 100, 200, 300]

    volumetric_core_losses_data = pandas.DataFrame()
    for material, material_id in materials.items():
        # if (material != "N27"):
        #     continue
        for magnetic_flux_density in magnetic_flux_densities:
            for frequency in frequencies:

                data = {
                    'pv_temp[cbMaterial]': str(material_id),
                    'pv_temp[rgParam]': '0',
                    'pv_temp[cbFirst]': str(frequency),
                    'pv_temp[cbSecond]': str(magnetic_flux_density),
                    'pv_temp[ok]': '',
                }

                proxies = None

                connected = False

                while not connected:
                    try:
                        response = requests.post(
                            'https://tools.tdk-electronics.tdk.com/mdt/index.php/pl_temp',
                            data=data,
                            proxies=proxies,
                        )

                        try:
                            data = pandas.read_html(StringIO(response.text))
                            if len(data) > 0 and not data[0].empty:
                                data = data[0]
                                print(material)
                                print(frequency)
                                print(magnetic_flux_density)
                                print(data)
                                core_losses_datum = data.rename(columns={"T [°C]": "temperature", "Pl / V [kW/m3]": "value", "Pvsin [kW / m3]": "value"})
                                if "α [-]" in core_losses_datum.columns:
                                    core_losses_datum = core_losses_datum.drop(["α [-]"], axis=1)
                                if "β [-]" in core_losses_datum.columns:
                                    core_losses_datum = core_losses_datum.drop(["β [-]"], axis=1)
                                core_losses_datum["value"] = core_losses_datum["value"] * 1000
                                core_losses_datum["magnetic_flux_density"] = [magnetic_flux_density / 1000] * len(core_losses_datum.index)
                                core_losses_datum["frequency"] = [frequency * 1000] * len(core_losses_datum.index)
                                core_losses_datum["material"] = [material] * len(core_losses_datum.index)
                                print(core_losses_datum)
                                volumetric_core_losses_data = pandas.concat([volumetric_core_losses_data, core_losses_datum], ignore_index=True)
                                connected = True

                                print(volumetric_core_losses_data)
                            else:
                                connected = True
                                print(f"No data for {material} at {magnetic_flux_density} mT and {frequency} kHz")

                        except ValueError:
                            print(material)
                            proxy = FreeProxy(https=True, rand=True).get()
                            print(proxy)

                            proxies = {
                                'https': proxy,
                            }

                            pass
                    except urllib3.exceptions.MaxRetryError:
                        proxy = FreeProxy(https=True, rand=True).get()
                        print(proxy)

                        proxies = {
                            'https': proxy,
                        }

    volumetric_core_losses_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_volumetric_core_losses.csv", index=False)


# Amplitude permeability
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_amplitude_permeability.csv"):
    amplitude_permeability_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_amplitude_permeability.csv")
else:

    materials = {
        "N27": 127,
        "N30": 128,
        "N41": 129,
        "N45": 130,
        "N48": 131,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
        "T35": 142,
        "T36": 143,
        "T37": 144,
        "T38": 145,
        "T46": 146,
        "T57": 147,
        "T65": 148,
        "T66": 149,
    }
    temperatures = [25, 100]
    amplitude_permeability_data = pandas.DataFrame()
    for material, material_id in materials.items():
        # if (material != "N27"):
        #     continue
        for temperature in temperatures:
            data = {
                'amplitude_permeability[cbMaterial]': str(material_id),
                'amplitude_permeability[rgTemperature]': str(temperature),
                'amplitude_permeability[ok]': '',
            }

            response = requests.post(
                'https://tools.tdk-electronics.tdk.com/mdt/index.php/amplitpermb',
                data=data,
            )
            data = pandas.read_html(StringIO(response.text))[0]

            if len(data) > 0:
                amplitude_permeability_datum = data.rename(columns={"B [mT]": "magnetic_flux_density", "µa [-]": "value"})
                amplitude_permeability_datum["magnetic_flux_density"] = amplitude_permeability_datum["magnetic_flux_density"] / 1000
                amplitude_permeability_datum["temperature"] = [temperature] * len(amplitude_permeability_datum.index)
                amplitude_permeability_datum["material"] = [material] * len(amplitude_permeability_datum.index)
                amplitude_permeability_data = pandas.concat([amplitude_permeability_data, amplitude_permeability_datum], ignore_index=True)
                print(amplitude_permeability_data)
            else:
                print(f"No amplitude permeability data for {material} at {temperature} C")
    print(amplitude_permeability_data)
    amplitude_permeability_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_amplitude_permeability.csv", index=False)


# Complex real permeability
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_real_permeability.csv"):
    complex_real_permeability_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_real_permeability.csv")
else:

    materials = {
        "C350": 121,
        "C351": 122,
        "K1": 123,
        "K10": 124,
        "M33": 125,
        "M34": 150,
        "N22": 126,
        "N27": 127,
        "N30": 128,
        "N41": 129,
        "N45": 130,
        "N48": 131,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
        "T35": 142,
        "T36": 143,
        "T37": 144,
        "T38": 145,
        "T46": 146,
        "T57": 147,
        "T65": 148,
        "T66": 149,
    }
    complex_real_permeability_data = pandas.DataFrame()
    for material, material_id in materials.items():
        # if (material != "N27"):
        #     continue
        data = {
            'complex_permeability[cbMaterial]': str(material_id),
            'complex_permeability[rgTemperature]': str(0),
            'complex_permeability[ok]': '',
        }

        print(material)
        response = requests.post(
            'https://tools.tdk-electronics.tdk.com/mdt/index.php/complexperm',
            data=data,
        )
        data = pandas.read_html(StringIO(response.text))[0]

        if len(data) > 0:
            complex_real_permeability_datum = data.rename(columns={"f [kHz]": "frequency", "µ` [-]": "value"})
            complex_real_permeability_datum["frequency"] = complex_real_permeability_datum["frequency"] * 1000
            complex_real_permeability_datum["material"] = [material] * len(complex_real_permeability_datum.index)
            complex_real_permeability_data = pandas.concat([complex_real_permeability_data, complex_real_permeability_datum], ignore_index=True)
            print(complex_real_permeability_data)
        else:
            print(f"No complex real permeability data for {material} C")
    print(complex_real_permeability_data)
    complex_real_permeability_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_real_permeability.csv", index=False)

# Complex imaginary permeability
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_imaginary_permeability.csv"):
    complex_imaginary_permeability_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_imaginary_permeability.csv")
else:

    materials = {
        "C350": 121,
        "C351": 122,
        "K1": 123,
        "K10": 124,
        "M33": 125,
        "M34": 150,
        "N22": 126,
        "N27": 127,
        "N30": 128,
        "N41": 129,
        "N45": 130,
        "N48": 131,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
        "T35": 142,
        "T36": 143,
        "T37": 144,
        "T38": 145,
        "T46": 146,
        "T57": 147,
        "T65": 148,
        "T66": 149,
    }
    complex_imaginary_permeability_data = pandas.DataFrame()
    for material, material_id in materials.items():
        # if (material != "N27"):
        #     continue
        data = {
            'complex_permeability[cbMaterial]': str(material_id),
            'complex_permeability[rgTemperature]': str(1),
            'complex_permeability[ok]': '',
        }

        print(material)
        response = requests.post(
            'https://tools.tdk-electronics.tdk.com/mdt/index.php/complexperm',
            data=data,
        )
        data = pandas.read_html(StringIO(response.text))[0]

        if len(data) > 0:
            complex_imaginary_permeability_datum = data.rename(columns={"f [kHz]": "frequency", "µ`` [-]": "value"})
            complex_imaginary_permeability_datum["frequency"] = complex_imaginary_permeability_datum["frequency"] * 1000
            complex_imaginary_permeability_datum["material"] = [material] * len(complex_imaginary_permeability_datum.index)
            complex_imaginary_permeability_data = pandas.concat([complex_imaginary_permeability_data, complex_imaginary_permeability_datum], ignore_index=True)
            print(complex_imaginary_permeability_data)
        else:
            print(f"No complex imaginary permeability data for {material} C")
    print(complex_imaginary_permeability_data)
    complex_imaginary_permeability_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_complex_imaginary_permeability.csv", index=False)

# BH cycle
if os.path.exists(f"{pathlib.Path(__file__).parent}/temp/tdk_bh_cycle.csv"):
    bh_cycle_data = pandas.read_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_bh_cycle.csv")
else:

    materials = {
        "C350": 121,
        "C351": 122,
        "K1": 123,
        "K10": 124,
        "M33": 125,
        "M34": 150,
        "N22": 126,
        "N27": 127,
        "N30": 128,
        "N41": 129,
        "N45": 130,
        "N48": 131,
        "N49": 132,
        "N51": 133,
        "N72": 134,
        "N87": 135,
        "N88": 136,
        "N92": 137,
        "N95": 138,
        "N96": 139,
        "N97": 140,
        "PC200": 141,
        "T35": 142,
        "T36": 143,
        "T37": 144,
        "T38": 145,
        "T46": 146,
        "T57": 147,
        "T65": 148,
        "T66": 149,
    }
    temperatures = [25, 100]
    bh_cycle_data = pandas.DataFrame()
    for material, material_id in materials.items():
        # if (material != "N27"):
        #     continue
        for temperature in temperatures:
            data = {
                'hysteresis[cbMaterial]': str(material_id),
                'hysteresis[rgTemperature]': str(temperature),
                'hysteresis[ok]': '',
            }

            response = requests.post(
                'https://tools.tdk-electronics.tdk.com/mdt/index.php/hysteresis',
                data=data,
            )
            data = pandas.read_html(StringIO(response.text))[0]

            if len(data) > 0:
                print(data)
                bh_cycle_datum = data.rename(columns={"B [mT]": "magnetic_flux_density", "H [A/m]": "magnetic_field_strength"})
                bh_cycle_datum["magnetic_flux_density"] = bh_cycle_datum["magnetic_flux_density"] / 1000
                bh_cycle_datum["temperature"] = [temperature] * len(bh_cycle_datum.index)
                bh_cycle_datum["material"] = [material] * len(bh_cycle_datum.index)
                bh_cycle_data = pandas.concat([bh_cycle_data, bh_cycle_datum], ignore_index=True)
                print(bh_cycle_data)
            else:
                print(f"No BH cycle data for {material} at {temperature} C")
    print(bh_cycle_data)
    bh_cycle_data.to_csv(f"{pathlib.Path(__file__).parent}/temp/tdk_bh_cycle.csv", index=False)


mas_data = []
mas_advanced_data = []

density_per_material = {
    "C350": 2930,
    "C351": 2930,
    "K1": 4650,
    "K10": 4700,
    "M33": 4500,
    "M34": 4800,
    "N22": 4700,
    "N27": 4750,
    "N30": 4800,
    "N41": 4800,
    "N45": 4930,
    "N48": 4700,
    "N49": 4800,
    "N51": 4800,
    "N72": 4800,
    "N87": 4850,
    "N88": 4900,
    "N92": 4850,
    "N95": 4900,
    "N96": 4850,
    "N97": 4920,
    "PC200": 4850,
    "T35": 4900,
    "T36": 4950,
    "T37": 4900,
    "T38": 4950,
    "T46": 5000,
    "T57": 4930,
    "T65": 4930,
    "T66": 5100,
}
curie_temperature_per_material = {
    "C350": 120,
    "C351": 200,
    "K1": 400,
    "K10": 150,
    "M33": 200,
    "M34": 200,
    "N22": 145,
    "N27": 220,
    "N30": 130,
    "N41": 220,
    "N45": 255,
    "N48": 170,
    "N49": 240,
    "N51": 220,
    "N72": 210,
    "N87": 210,
    "N88": 220,
    "N92": 280,
    "N95": 220,
    "N96": 240,
    "N97": 230,
    "PC200": 280,
    "T35": 130,
    "T36": 130,
    "T37": 130,
    "T38": 130,
    "T46": 130,
    "T57": 140,
    "T65": 160,
    "T66": 100,
}
resistivity_per_material = {
    "C350": 500,
    "C351": 500,
    "K1": 100000,
    "K10": 100000,
    "M33": 5,
    "M34": 4,
    "N22": 1,
    "N27": 3,
    "N30": 0.5,
    "N41": 2,
    "N45": 11,
    "N48": 3,
    "N49": 17,
    "N51": 6,
    "N72": 12,
    "N87": 10,
    "N88": 17,
    "N92": 8,
    "N95": 6,
    "N96": 17,
    "N97": 8,
    "PC200": 22,
    "T35": 0.2,
    "T36": 0.2,
    "T37": 0.2,
    "T38": 0.1,
    "T46": 0.01,
    "T57": 3,
    "T65": 0.3,
    "T66": 0.8,
}

family_per_material = {
    "C350": "C",
    "C351": "C",
    "K1": "K",
    "K10": "K",
    "M33": "M",
    "M34": "M",
    "N22": "N",
    "N27": "N",
    "N30": "N",
    "N41": "N",
    "N45": "N",
    "N48": "N",
    "N49": "N",
    "N51": "N",
    "N72": "N",
    "N87": "N",
    "N88": "N",
    "N92": "N",
    "N95": "N",
    "N96": "N",
    "N97": "N",
    "PC200": "PC",
    "T35": "T",
    "T36": "T",
    "T37": "T",
    "T38": "T",
    "T46": "T",
    "T57": "T",
    "T65": "T",
    "T66": "T",
}

material_composition_per_material = {
    "C350": "NiZn",
    "C351": "NiZn",
    "K1": "NiZn",
    "K10": "NiZn",
    "M33": "MnZn",
    "M34": "NiZn",
    "N22": "MnZn",
    "N27": "NiZn",
    "N30": "MnZn",
    "N41": "MnZn",
    "N45": "MnZn",
    "N48": "MnZn",
    "N49": "MnZn",
    "N51": "MnZn",
    "N72": "MnZn",
    "N87": "MnZn",
    "N88": "MnZn",
    "N92": "MnZn",
    "N95": "MnZn",
    "N96": "MnZn",
    "N97": "MnZn",
    "PC200": "MnZn",
    "T35": "MnZn",
    "T36": "MnZn",
    "T37": "MnZn",
    "T38": "MnZn",
    "T46": "MnZn",
    "T57": "MnZn",
    "T65": "MnZn",
    "T66": "MnZn",
}


datasheet_per_material = {
    "C350": "https://www.tdk-electronics.tdk.com/inf/80/db/fer/ferrite_polymer_composites.pdf",
    "C351": "https://www.tdk-electronics.tdk.com/inf/80/db/fer/ferrite_polymer_composites.pdf",
    "K1": "https://www.tdk-electronics.tdk.com/download/528860/df0caad70ee5a0b56cf1d51084fc1fc2/pdf-k1.pdf",
    "K10": "https://www.tdk-electronics.tdk.com/download/528854/d6a20f8b9ec044975113d05211eefff5/pdf-k10.pdf",
    "M33": "https://www.tdk-electronics.tdk.com/download/528878/4f65b2fcd11b8b6b3e99fa2bc44bfd7e/pdf-m33.pdf",
    "M34": "https://www.tdk-electronics.tdk.com/download/528878/4f65b2fcd11b8b6b3e99fa2bc44bfd7e/pdf-m34.pdf",
    "N22": "https://www.tdk-electronics.tdk.com/download/528884/6bb6856ed268a7b9100d14af2d5f3444/pdf-n22.pdf",
    "N27": "https://www.tdk-electronics.tdk.com/download/528850/3a7f957d754f899aec42cd946598c5c4/pdf-n27.pdf",
    "N30": "https://www.tdk-electronics.tdk.com/download/528848/fd7b13cb06de3e2edca6c617ac9de7da/pdf-n30.pdf",
    "N41": "https://www.tdk-electronics.tdk.com/download/528876/07f8fe5e7fb654a26dc779e207a5de15/pdf-n41.pdf",
    "N45": "https://www.tdk-electronics.tdk.com/download/528844/0ca69ad6f73de03e7c30c43b39d0aedf/pdf-n45.pdf",
    "N48": "https://www.tdk-electronics.tdk.com/download/528870/80183d6ef4a48b49c5bee48a6e8b2ac5/pdf-n48.pdf",
    "N49": "https://www.tdk-electronics.tdk.com/download/528856/79b67daa1d54253aa2ae371bb98e9629/pdf-n49.pdf",
    "N51": "https://www.tdk-electronics.tdk.com/download/528840/1b67f37ebca1f89a3d021563b44b4d12/pdf-n51.pdf",
    "N72": "https://www.tdk-electronics.tdk.com/download/528858/fe38d2b52133890d5f7cc6885ff05e75/pdf-n72.pdf",
    "N87": "https://www.tdk-electronics.tdk.com/download/528882/990c299b916e9f3eb7e44ad563b7f0b9/pdf-n87.pdf",
    "N88": "https://www.tdk-electronics.tdk.com/download/528890/cee245ee621a6b10bd5897f906870869/pdf-n88.pdf",
    "N92": "https://www.tdk-electronics.tdk.com/download/528888/4e7e6721e0cb00bfb8893586d22189db/pdf-n92.pdf",
    "N95": "https://www.tdk-electronics.tdk.com/download/528866/d2caba056c7aa06d1e24ece60fb6ed5e/pdf-n95.pdf",
    "N96": "https://www.tdk-electronics.tdk.com/download/528892/62108148dfc602168d141f831c2914e8/pdf-n96.pdf",
    "N97": "https://www.tdk-electronics.tdk.com/download/528886/f7953b4caf14b9c4a15e34b60fb09812/pdf-n97.pdf",
    "PC200": "https://www.tdk-electronics.tdk.com/download/2111336/b9f1a9402fb37cd23fc6c0d495eab3b2/pdf-pc200.pdf",
    "T35": "https://www.tdk-electronics.tdk.com/download/528880/8351bfa4e86ee7129468f0e83b5db861/pdf-t35.pdf",
    "T36": "https://www.tdk-electronics.tdk.com/download/528842/3f54dea2fe7849a8a51f65c2196de15d/pdf-t36.pdf",
    "T37": "https://www.tdk-electronics.tdk.com/download/528874/8251c697a6f71dc2bbf4c1f5362137ed/pdf-t37.pdf",
    "T38": "https://www.tdk-electronics.tdk.com/download/528868/c2cb0f1f5413cd8e83b8ecf7e85152f9/pdf-t38.pdf",
    "T46": "https://www.tdk-electronics.tdk.com/download/528864/23ff89c1652f23879b6710e28ab4ff96/pdf-t46.pdf",
    "T57": "https://www.tdk-electronics.tdk.com/download/528862/91a07e0c1d9d9a5188cd1708d4019f02/pdf-t57.pdf",
    "T65": "https://www.tdk-electronics.tdk.com/download/528846/1598c89da613db594840e493e0ffb9fe/pdf-t65.pdf",
    "T66": "https://www.tdk-electronics.tdk.com/download/528852/528a1f36fa90048c52cd8b25ca078d17/pdf-t66.pdf",
}

materials = [
    "C350",
    "C351",
    "K1",
    "K10",
    "M33",
    "M34",
    "N22",
    "N27",
    "N30",
    "N41",
    "N45",
    "N48",
    "N49",
    "N51",
    "N72",
    "N87",
    "N88",
    "N92",
    "N95",
    "N96",
    "N97",
    "PC200",
    "T35",
    "T36",
    "T37",
    "T38",
    "T46",
    "T57",
    "T65",
    "T66",
]

for material in materials:
    # if (material != "N27"):
    #     continue
    print(material)
    mas_datum = {
        "type": "commercial",
        "name": material,
        "family": None,
        "resistivity": [],
        "remanence": [],
        "coerciveForce": [],
        "density": None,
        "curieTemperature": None,
        "material": "ferrite",
        "materialComposition": None,
        "manufacturerInfo": {"name": "TDK", "datasheetUrl": None},
        "permeability": {
            "initial": []
        },
        "saturation": [],
        "volumetricLosses": {"default": []}
    }
    mas_advanced_datum = {
        "name": material,
        "manufacturerInfo": {"name": "TDK"},
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
    material_volumetric_core_losses_data = volumetric_core_losses_data[volumetric_core_losses_data['material'] == material]

    mas_datum["resistivity"] = [{"value": resistivity_per_material[material], "temperature": 25}]
    mas_datum["family"] = family_per_material[material]
    mas_datum["density"] = density_per_material[material]
    mas_datum["curieTemperature"] = curie_temperature_per_material[material]
    mas_datum["materialComposition"] = material_composition_per_material[material]
    mas_datum["manufacturerInfo"]["datasheetUrl"] = datasheet_per_material[material]

    mas_advanced_datum["volumetricLosses"]["default"].append([])
    for row_index, row in material_volumetric_core_losses_data.iterrows():
        B = row["magnetic_flux_density"]
        f = row["frequency"]
        temperature = row["temperature"]

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
                "value": row["value"],
                "origin": "manufacturer"
            })

    material_initial_permeability_data = initial_permeability_data[initial_permeability_data['material'] == material]
    # print(material_initial_permeability_data)
    for row_index, row in material_initial_permeability_data.iterrows():
        temperature = row["temperature"]
        mas_datum["permeability"]["initial"].append({
            "temperature": temperature,
            "value": row["value"]
        })

    material_complex_real_permeability_data = complex_real_permeability_data[complex_real_permeability_data['material'] == material]
    for row_index, row in material_complex_real_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["complex"]["real"].append({
            "frequency": row["frequency"],
            "value": row["value"]
        })

    material_complex_imaginary_permeability_data = complex_imaginary_permeability_data[complex_imaginary_permeability_data['material'] == material]
    for row_index, row in material_complex_imaginary_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["complex"]["imaginary"].append({
            "frequency": row["frequency"],
            "value": row["value"]
        })

    material_amplitude_permeability_data = amplitude_permeability_data[amplitude_permeability_data['material'] == material]
    for row_index, row in material_amplitude_permeability_data.iterrows():
        mas_advanced_datum["permeability"]["amplitude"].append({
            "value": row["value"],
            "temperature": row["temperature"],
            "magneticFluxDensityPeak": row["magnetic_flux_density"]
        })

    coercivity_per_temperature = {}
    remanence_per_temperature = {}
    current_closest_h_to_0 = {}
    current_closest_b_to_0 = {}
    max_h = {}
    max_b = {}
    material_bh_cycle_data = bh_cycle_data[bh_cycle_data['material'] == material]
    for row_index, row in material_bh_cycle_data.iterrows():
        temperature = row["temperature"]
        current_closest_h_to_0[temperature] = 999999
        current_closest_b_to_0[temperature] = 999999
        max_h[temperature] = 0
        max_b[temperature] = 0

    for row_index, row in material_bh_cycle_data.iterrows():
        magnetic_flux_density = row["magnetic_flux_density"]
        magnetic_field_strength = row["magnetic_field_strength"]
        temperature = row["temperature"]

        if abs(magnetic_flux_density) > max_b[temperature]:
            max_b[temperature] = abs(magnetic_flux_density)
        if abs(magnetic_field_strength) > max_h[temperature]:
            max_h[temperature] = abs(magnetic_field_strength)

        if abs(magnetic_field_strength) < current_closest_h_to_0[temperature]:
            current_closest_h_to_0[temperature] = abs(magnetic_field_strength)
            remanence_per_temperature[int(temperature)] = abs(magnetic_flux_density)

        if abs(magnetic_flux_density) < current_closest_b_to_0[temperature]:
            current_closest_b_to_0[temperature] = abs(magnetic_flux_density)
            coercivity_per_temperature[int(temperature)] = abs(magnetic_field_strength)

            mas_advanced_datum["bhCycle"].append({
                "magneticFluxDensity": magnetic_flux_density,
                "magneticField": magnetic_field_strength,
                "temperature": temperature
            })

    for temperature, value in max_h.items():
        mas_datum["saturation"].append({
            "magneticFluxDensity": max_b[temperature],
            "magneticField": max_h[temperature],
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

    # import pprint
    # pprint.pprint(mas_datum)
    # assert 0

    mas_data.append(mas_datum)
    # Strip empty placeholder blocks before persisting: advanced_core_materials records
    # are name-keyed patches merged onto the base material; an empty block carries no
    # information and corrupts both the MKF merge and merged-schema validation.
    _p = mas_advanced_datum.get("permeability")
    if isinstance(_p, dict):
        if _p.get("amplitude") == []:
            _p.pop("amplitude")
        _c = _p.get("complex")
        if isinstance(_c, dict) and _c.get("real") == [] and _c.get("imaginary") == []:
            _p.pop("complex")
        if not _p:
            mas_advanced_datum.pop("permeability")
    _vl = mas_advanced_datum.get("volumetricLosses")
    if isinstance(_vl, dict):
        _d = [x for x in _vl.get("default", []) if x != []]
        if _d:
            _vl["default"] = _d
        else:
            mas_advanced_datum.pop("volumetricLosses")
    if len(mas_advanced_datum.get("bhCycle", [])) < 4:
        # a B-H "cycle" of fewer than 4 points is an extraction fragment (schema minItems: 4)
        mas_advanced_datum.pop("bhCycle", None)
    _vl2 = mas_advanced_datum.get("volumetricLosses")
    if isinstance(_vl2, dict):
        for _i, _entry in enumerate(_vl2.get("default", [])):
            if isinstance(_entry, list):
                _seen = set(); _out = []
                for _pt in _entry:
                    _k = json.dumps(_pt, sort_keys=True)
                    if _k in _seen:
                        continue  # exact duplicate point (schema uniqueItems)
                    _seen.add(_k); _out.append(_pt)
                _vl2["default"][_i] = _out
    mas_advanced_data.append(mas_advanced_datum)

# print(mas_data)
# print(mas_advanced_data)

ndjson.dump(mas_data, open(f"{pathlib.Path(__file__).parent.resolve()}/tdk_core_materials.ndjson", "w"), ensure_ascii=False)
ndjson.dump(mas_advanced_data, open(f"{pathlib.Path(__file__).parent.resolve()}/tdk_advanced_core_materials.ndjson", "w"), ensure_ascii=False)
