from scipy.optimize import curve_fit
import numpy
import json
import pandas
import pathlib
import io
import chardet
_MKF_ROOT = __import__("pathlib").Path(__file__).resolve().parents[2]  # ABT #1596: this checkout, not a developer path




def convert(filename):
    with open(filename, "rb") as f:
        byte_data = f.read()
        # print(byte_data)
        data = byte_data.decode('MacRoman')
        print(data)
        # print(chardet.detect(byte_data))




if __name__ == '__main__':  # pragma: no cover
    convert(str(_MKF_ROOT / "src/tools/Inductor.dev"))
