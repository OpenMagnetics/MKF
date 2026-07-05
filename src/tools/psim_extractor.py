from scipy.optimize import curve_fit
import numpy
import json
import pandas
import pathlib
import io
import chardet




def convert(filename):
    with open(filename, "rb") as f:
        byte_data = f.read()
        # print(byte_data)
        data = byte_data.decode('MacRoman')
        print(data)
        # print(chardet.detect(byte_data))




if __name__ == '__main__':  # pragma: no cover
    convert("/home/alf/OpenMagnetics/MKF/src/tools/Inductor.dev")
