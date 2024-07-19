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
    return {'k': k, 'alpha': alpha, 'beta': beta, 'ct0': ct0, 'ct1': ct1, 'ct2': ct2}


def fit(material, frequency):

    data_path = pathlib.Path(__file__).parent.resolve().joinpath('../../../MAS/data/advanced_core_materials.ndjson')
    all_data = pandas.read_json(data_path, lines=True)
    data = all_data[all_data['name'] == material]['volumetricLosses'].iloc[0]
    for method in data['default']:
        if isinstance(method, list):
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
                    else:

                        data_greater = data[data['frequency'] >= frequency]
                        data_lesser = data[data['frequency'] <= frequency]
                        data_greater = data_greater.sort_values('frequency', ascending=True).head(number_rows)
                        data_lesser = data_lesser.sort_values('frequency', ascending=False).head(number_rows)

                        data = pandas.concat([data_greater, data_lesser])

                    data['magneticFluxDensityPeak'] = data.apply(lambda row: row['magneticFluxDensity']['magneticFluxDensity']['processed']['peak'], axis=1)
                    data = data.drop('magneticFluxDensity', axis=1)
                    data = data.drop('origin', axis=1)

                    coefficients = calculate_steinmetz_coefficients(data)

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
    coefficients = fit("N49", [1, 300000])
    print(coefficients)
    coefficients = fit("N49", [300000, 600000])
    print(coefficients)
    coefficients = fit("N49", [600000, 800000])
    print(coefficients)
    coefficients = fit("N49", [800000, 1000000])
    print(coefficients)
