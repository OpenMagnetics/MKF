from scipy.optimize import curve_fit
import numpy
import json
import pandas
import pathlib


def calculate_steinmetz_coefficients(data):
    def steinmetz_fitting_equation(log_f_B_T, log_k, alpha, beta, ct0, ct1, ct2):
        log_f, log_B, T = log_f_B_T
        temp_coeff = ct2 * T**2 - ct1 * T + ct0
        if (temp_coeff < 0).any():
            return log_k + alpha * log_f + beta * log_B
        else:
            return log_k + alpha * log_f + beta * log_B + numpy.log10(temp_coeff)

    [log_k, alpha, beta, ct0, ct1, ct2], pcov = curve_fit(steinmetz_fitting_equation, (numpy.log10(data['frequency']), numpy.log10(data['magneticFluxDensityPeak']), data['temperature']), numpy.log10(data['volumetricLosses']), [0.4, 1.5, 2, 0.001, 1e-06, 1e-08])
    k = 10**log_k
    return {'k': float(k), 'alpha': float(alpha), 'beta': float(beta), 'ct0': float(ct0), 'ct1': float(ct1), 'ct2': float(ct2)}


def calculate_steinmetz_coefficients_no_temp(data):
    def steinmetz_fitting_equation(log_f_B_T, log_k, alpha, beta):
        log_f, log_B, T = log_f_B_T
        return log_k + alpha * log_f + beta * log_B

    [log_k, alpha, beta], pcov = curve_fit(steinmetz_fitting_equation, (numpy.log10(data['frequency']), numpy.log10(data['magneticFluxDensityPeak']), data['temperature']), numpy.log10(data['volumetricLosses']), [0.4, 1.5, 2])
    k = 10**log_k
    return {'k': float(k), 'alpha': float(alpha), 'beta': float(beta)}


def fit(material, frequency):

    data_path = pathlib.Path(__file__).parent.resolve().joinpath('../../../MAS/data/advanced_core_materials.ndjson')
    all_data = pandas.read_json(data_path, lines=True)
    data = all_data[all_data['name'] == material]['volumetricLosses'].iloc[0]
    print(material)
    for method in data['default']:
        if isinstance(method, list):
            if len(method) == 0:
                return {}

            number_rows = 200
            while True:
                try:
                    data = pandas.DataFrame(method)
                    if "MagNet" in data['origin'].to_list():
                        data = data[data['origin'] == "MagNet"]
                    data['frequency'] = data.apply(lambda row: row['magneticFluxDensity']['frequency'], axis=1)
                    data = data.rename(columns={"value": "volumetricLosses"})
                    maxmin_frequency = data['frequency'].max() - data['frequency'].min()

                    if isinstance(frequency, list):
                        data = data[data['frequency'] >= frequency[0]]
                        data = data[data['frequency'] <= frequency[1]]
                        minimumFrequency = frequency[0]
                        maximumFrequency = frequency[1]
                    else:

                        data_greater = data[data['frequency'] >= frequency]
                        data_lesser = data[data['frequency'] <= frequency]
                        data_greater = data_greater.sort_values('frequency', ascending=True).head(number_rows)
                        data_lesser = data_lesser.sort_values('frequency', ascending=False).head(number_rows)
                        minimumFrequency = data['frequency'].min()
                        maximumFrequency = data['frequency'].max()

                        data = pandas.concat([data_greater, data_lesser])

                    if data.empty:
                        return {}

                    data['magneticFluxDensityPeak'] = data.apply(lambda row: row['magneticFluxDensity']['magneticFluxDensity']['processed']['peak'], axis=1)
                    data = data.drop('magneticFluxDensity', axis=1)
                    data = data.drop('origin', axis=1)

                    coefficients = calculate_steinmetz_coefficients(data)
                    # coefficients = calculate_steinmetz_coefficients_no_temp(data)
                    coefficients["minimumFrequency"] = minimumFrequency
                    coefficients["maximumFrequency"] = maximumFrequency

                    # if coefficients['alpha'] < 1:
                    #     number_rows += 100
                    #     print(f"Retrying with number_rows {number_rows}")
                    #     continue
                    # if coefficients['k'] > 10000000:
                    #     number_rows += 100
                    #     print(f"Retrying with number_rows {number_rows}")
                    #     continue
                    return coefficients
                except RuntimeError:
                    number_rows += 100
                    print(f"Retrying with number_rows {number_rows}")

    else:
        return {}


if __name__ == '__main__':  # pragma: no cover
    # coefficients = fit("SMP44", [1, 200000])
    # print(coefficients)
    # coefficients = fit("SMP44", [200000, 600000])
    # print(coefficients)

    # print({"method": "steinmetz", "ranges": [fit("P47", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P45", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P491", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P6", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P61", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P62", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P5", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P51", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("P52", [1, 20000000])]})
    # print({"method": "steinmetz", "ranges": [fit("JSP95", [1, 20000000])]})

    # print({"method": "steinmetz", "ranges": [fit("SMP44", [1, 200000]), fit("SMP44", [200000, 600000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP47", [1, 200000]), fit("SMP47", [200000, 600000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP95", [1, 200000]), fit("SMP95", [200000, 600000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP96", [1, 200000]), fit("SMP96", [200000, 600000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP97", [1, 200000]), fit("SMP97", [200000, 600000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP50", [300000, 1000000]), fit("SMP50", [1000000, 5000000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP51", [300000, 1000000]), fit("SMP51", [1000000, 5000000])]})
    # print({"method": "steinmetz", "ranges": [fit("SMP53", [300000, 1000000]), fit("SMP53", [1000000, 5000000])]})

    # print({"method": "steinmetz", "ranges": [fit("DMR25", [1, 200000]), fit("DMR25", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR24", [1, 200000]), fit("DMR24", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR44", [1, 200000]), fit("DMR44", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR40B", [1, 200000]), fit("DMR40B", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR95", [1, 200000]), fit("DMR95", [200000, 600000000])]})
    print({"method": "steinmetz", "ranges": [fit("DMR40", [1, 200000]), fit("DMR40", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR55", [1, 200000]), fit("DMR55", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR28", [1, 200000]), fit("DMR28", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR96", [1, 200000]), fit("DMR96", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR47", [1, 200000]), fit("DMR47", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR90", [1, 200000]), fit("DMR90", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR71", [1, 200000]), fit("DMR71", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR73", [1, 200000]), fit("DMR73", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR70", [1, 200000]), fit("DMR70", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R5K", [1, 200000]), fit("R5K", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R7K", [1, 200000]), fit("R7K", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R10K", [1, 200000]), fit("R10K", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R12K", [1, 200000]), fit("R12K", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R15K", [1, 200000]), fit("R15K", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R10KZ", [1, 200000]), fit("R10KZ", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR50", [1, 200000]), fit("DMR50", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR50B", [1, 200000]), fit("DMR50B", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR51", [1, 200000]), fit("DMR51", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR52", [1, 200000]), fit("DMR52", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN30B", [1, 200000]), fit("DN30B", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN40B", [1, 200000]), fit("DN40B", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN50B", [1, 200000]), fit("DN50B", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN85H", [1, 200000]), fit("DN85H", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN150H", [1, 200000]), fit("DN150H", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN200L", [1, 200000]), fit("DN200L", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR53", [1, 200000]), fit("DMR53", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR96A", [1, 200000]), fit("DMR96A", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR91", [1, 200000]), fit("DMR91", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R7KC", [1, 200000]), fit("R7KC", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R12KZ", [1, 200000]), fit("R12KZ", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R15KZ", [1, 200000]), fit("R15KZ", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R5KC", [1, 200000]), fit("R5KC", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R5KZ", [1, 200000]), fit("R5KZ", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("R10KC", [1, 200000]), fit("R10KC", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR51W", [1, 200000]), fit("DMR51W", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DMR52W", [1, 200000]), fit("DMR52W", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN15P", [1, 200000]), fit("DN15P", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN20F", [1, 200000]), fit("DN20F", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN2S", [1, 200000]), fit("DN2S", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN33L", [1, 200000]), fit("DN33L", [200000, 600000000])]})
    # print({"method": "steinmetz", "ranges": [fit("DN100H", [1, 200000]), fit("DN100H", [200000, 600000000])]})

    # print({"method": "steinmetz", "ranges": [fit("Finemet", [1, 500000])]})
    # print({"method": "steinmetz", "ranges": [fit("AFCC", [1, 500000])]})
    # print({"method": "steinmetz", "ranges": [fit("ANCC", [1, 500000])]})
    # print({"method": "steinmetz", "ranges": [fit("C351", [1, 200000]), fit("C351", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("K1", [1, 200000]), fit("K1", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("K10", [1, 200000]), fit("K10", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("M33", [1, 200000]), fit("M33", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("M34", [1, 200000]), fit("M34", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N22", [1, 200000]), fit("N22", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N27", [1, 200000]), fit("N27", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N30", [1, 200000]), fit("N30", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N41", [1, 200000]), fit("N41", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N45", [1, 200000]), fit("N45", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N48", [1, 200000]), fit("N48", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N49", [1, 200000]), fit("N49", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N51", [1, 200000]), fit("N51", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N72", [1, 200000]), fit("N72", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N87", [1, 200000]), fit("N87", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N88", [1, 200000]), fit("N88", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N92", [1, 200000]), fit("N92", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N95", [1, 200000]), fit("N95", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N96", [1, 200000]), fit("N96", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("N97", [1, 200000]), fit("N97", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("PC200", [1, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T35", [1, 200000]), fit("T35", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T36", [1, 200000]), fit("T36", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T37", [1, 200000]), fit("T37", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T38", [1, 200000]), fit("T38", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T46", [1, 200000]), fit("T46", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T57", [1, 200000]), fit("T57", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T65", [1, 200000]), fit("T65", [200000, 60000000])]})
    # print({"method": "steinmetz", "ranges": [fit("T66", [1, 200000]), fit("T66", [200000, 60000000])]})
