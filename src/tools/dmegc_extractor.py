import pandas
import ndjson
import pathlib
import numpy
import re

density_per_material = {
    "DMR25": 4900,
    "DMR24": 4900,
    "DMR44": 4800,
    "DMR40B": 4800,
    "DMR95": 4900,
    "DMR40": 4800,
    "DMR55": 4800,
    "DMR53": 4800,
    "DMR28": 4900,
    "DMR96": 4800,
    "DMR96A": 4800,
    "DMR47": 4800,
    "DMR90": 4850,
    "DMR91": 4900,
    "DMR71": 4850,
    "DMR73": 4900,
    "DMR70": 4800,
    "R5K": 4850,
    "R7K": 4900,
    "R7KC": 4900,
    "R10K": 4900,
    "R12K": 4900,
    "R12KZ": 4900,
    "R15K": 4900,
    "R15KZ": 4900,
    "R5KC": 4850,
    "R5KZ": 4900,
    "R10KZ": 4900,
    "R10KC": 4900,
    "DMR50": 4800,
    "DMR50B": 4700,
    "DMR51": 4700,
    "DMR51W": 4700,
    "DMR52": 4700,
    "DMR52W": 4700,
    "DN30B": 5200,
    "DN40B": 5200,
    "DN50B": 5200,
    "DN85H": 5100,
    "DN150H": 5100,
    "DN100H": 5100,
    "DN15P": 5200,
    "DN200L": 5200,
    "DN20F": 5200,
    "DN2S": 5000,
    "DN33L": 5200
}
curie_temperature_per_material = {
    "DMR25": 240,
    "DMR24": 280,
    "DMR44": 215,
    "DMR40B": 225,
    "DMR95": 215,
    "DMR40": 215,
    "DMR55": 250,
    "DMR53": 280,
    "DMR28": 300,
    "DMR96": 230,
    "DMR96A": 230,
    "DMR47": 230,
    "DMR90": 280,
    "DMR91": 280,
    "DMR71": 255,
    "DMR73": 160,
    "DMR70": 170,
    "R5K": 140,
    "R7K": 125,
    "R7KC": 170,
    "R10K": 120,
    "R12K": 110,
    "R12KZ": 130,
    "R15K": 110,
    "R15KZ": 130,
    "R5KC": 170,
    "R5KZ": 135,
    "R10KZ": 120,
    "R10KC": 155,
    "DMR50": 240,
    "DMR50B": 240,
    "DMR51": 290,
    "DMR51W": 290,
    "DMR52": 280,
    "DMR52W": 290,
    "DN30B": 260,
    "DN40B": 240,
    "DN50B": 220,
    "DN85H": 150,
    "DN150H": 100,
    "DN100H": 130,
    "DN15P": 245,
    "DN200L": 110,
    "DN20F": 260,
    "DN2S": 220,
    "DN33L": 240
}
remanence_per_material = {
    "R10KZ": [{"magneticFluxDensity": 0.1, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.12, "magneticField": 0, "temperature": 25}],
    "DMR53": [{"magneticFluxDensity": 0.085, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.1, "magneticField": 0, "temperature": 25}],
    "DMR96A": [{"magneticFluxDensity": 0.045, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.055, "magneticField": 0, "temperature": 25}],
    "DMR91": [{"magneticFluxDensity": 0.042, "magneticField": 0, "temperature": 120}, {"magneticFluxDensity": 0.062, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.2, "magneticField": 0, "temperature": 25}],
    "DMR51W": [{"magneticFluxDensity": 0.08, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.12, "magneticField": 0, "temperature": 25}],
    "DMR52W": [{"magneticFluxDensity": 0.05, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.06, "magneticField": 0, "temperature": 25}],
    "DN15P": [{"magneticFluxDensity": 0.218, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.31, "magneticField": 0, "temperature": 25}],
    "DN20F": [{"magneticFluxDensity": 0.24, "magneticField": 0, "temperature": 100}, {"magneticFluxDensity": 0.32, "magneticField": 0, "temperature": 25}],
    "DN2S": [{"magneticFluxDensity": 0.1, "magneticField": 0, "temperature": 25}],
    "DN33L": [{"magneticFluxDensity": 0.35, "magneticField": 0, "temperature": 25}],
}
coercivity_per_material = {
    "R10KZ": [{"magneticFluxDensity": 0, "magneticField": 0.2, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 3, "temperature": 25}],
    "DMR53": [{"magneticFluxDensity": 0, "magneticField": 30, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 35, "temperature": 25}],
    "DMR96A": [{"magneticFluxDensity": 0, "magneticField": 8, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 10, "temperature": 25}],
    "DMR91": [{"magneticFluxDensity": 0, "magneticField": 4.3, "temperature": 120}, {"magneticFluxDensity": 0, "magneticField": 4.2, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 11, "temperature": 25}],
    "DMR51W": [{"magneticFluxDensity": 0, "magneticField": 38, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 45, "temperature": 25}],
    "DMR52W": [{"magneticFluxDensity": 0, "magneticField": 30, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 40, "temperature": 25}],
    "DN15P": [{"magneticFluxDensity": 0, "magneticField": 80, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 106, "temperature": 25}],
    "DN20F": [{"magneticFluxDensity": 0, "magneticField": 40, "temperature": 100}, {"magneticFluxDensity": 0, "magneticField": 60, "temperature": 25}],
    "DN2S": [{"magneticFluxDensity": 0, "magneticField": 630, "temperature": 25}],
    "DN33L": [{"magneticFluxDensity": 0, "magneticField": 65, "temperature": 25}],
}
saturation_per_material = {
    "R10KZ": [{"magneticFluxDensity": 0.260, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.450, "magneticField": 1194, "temperature": 25}],
    "DMR53": [{"magneticFluxDensity": 0.46, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.56, "magneticField": 1194, "temperature": 25}],
    "DMR96A": [{"magneticFluxDensity": 0.43, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.54, "magneticField": 1194, "temperature": 25}],
    "DMR91": [{"magneticFluxDensity": 0.46, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.55, "magneticField": 1194, "temperature": 25}],
    "R7KC": [{"magneticFluxDensity": 0.35, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.49, "magneticField": 1194, "temperature": 25}],
    "R12KZ": [{"magneticFluxDensity": 0.25, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.43, "magneticField": 1194, "temperature": 25}],
    "R15KZ": [{"magneticFluxDensity": 0.26, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.43, "magneticField": 1194, "temperature": 25}],
    "R5KC": [{"magneticFluxDensity": 0.35, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.48, "magneticField": 1194, "temperature": 25}],
    "R5KZ": [{"magneticFluxDensity": 0.27, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.45, "magneticField": 1194, "temperature": 25}],
    "R10KC": [{"magneticFluxDensity": 0.31, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.45, "magneticField": 1194, "temperature": 25}],
    "DMR51W": [{"magneticFluxDensity": 0.43, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.50, "magneticField": 1194, "temperature": 25}],
    "DMR52W": [{"magneticFluxDensity": 0.45, "magneticField": 1194, "temperature": 100}, {"magneticFluxDensity": 0.54, "magneticField": 1194, "temperature": 25}],
    "DN15P": [{"magneticFluxDensity": 0.36, "magneticField": 4000, "temperature": 100}, {"magneticFluxDensity": 0.44, "magneticField": 4000, "temperature": 25}],
    "DN20F": [{"magneticFluxDensity": 0.34, "magneticField": 4000, "temperature": 100}, {"magneticFluxDensity": 0.46, "magneticField": 4000, "temperature": 25}],
    "DN2S": [{"magneticFluxDensity": 0.24, "magneticField": 4000, "temperature": 25}],
    "DN33L": [{"magneticFluxDensity": 0.48, "magneticField": 4000, "temperature": 25}],
}
resistivity_per_material = {
    "DMR25": {"value": 8, "temperature": 25},
    "DMR24": {"value": 8, "temperature": 25},
    "DMR44": {"value": 2, "temperature": 25},
    "DMR40B": {"value": 6.5, "temperature": 25},
    "DMR95": {"value": 6, "temperature": 25},
    "DMR40": {"value": 6.5, "temperature": 25},
    "DMR53": {"value": 6, "temperature": 25},
    "DMR55": {"value": 6, "temperature": 25},
    "DMR28": {"value": 6, "temperature": 25},
    "DMR96": {"value": 6, "temperature": 25},
    "DMR96A": {"value": 6, "temperature": 25},
    "DMR47": {"value": 3.5, "temperature": 25},
    "DMR90": {"value": 6, "temperature": 25},
    "DMR91": {"value": 6, "temperature": 25},
    "DMR71": {"value": 6, "temperature": 25},
    "DMR73": {"value": 6, "temperature": 25},
    "DMR70": {"value": 6, "temperature": 25},
    "R5K": {"value": 0.5, "temperature": 25},
    "R7K": {"value": 0.2, "temperature": 25},
    "R7KC": {"value": 0.2, "temperature": 25},
    "R10K": {"value": 0.15, "temperature": 25},
    "R12K": {"value": 0.15, "temperature": 25},
    "R12KZ": {"value": 0.15, "temperature": 25},
    "R15K": {"value": 0.05, "temperature": 25},
    "R15KZ": {"value": 0.05, "temperature": 25},
    "R5KC": {"value": 0.5, "temperature": 25},
    "R5KZ": {"value": 0.5, "temperature": 25},
    "R10KZ": {"value": 0.15, "temperature": 25},
    "R10KC": {"value": 0.15, "temperature": 25},
    "DMR50": {"value": 6, "temperature": 25},
    "DMR50B": {"value": 6, "temperature": 25},
    "DMR51": {"value": 6, "temperature": 25},
    "DMR51W": {"value": 6, "temperature": 25},
    "DMR52": {"value": 6, "temperature": 25},
    "DMR52W": {"value": 6, "temperature": 25},
    "DN30B": {"value": 1e5, "temperature": 25},
    "DN40B": {"value": 1e5, "temperature": 25},
    "DN50B": {"value": 1e5, "temperature": 25},
    "DN85H": {"value": 1e5, "temperature": 25},
    "DN150H": {"value": 1e5, "temperature": 25},
    "DN100H": {"value": 1e5, "temperature": 25},
    "DN15P": {"value": 1e5, "temperature": 25},
    "DN200L": {"value": 1e6, "temperature": 25},
    "DN20F": {"value": 1e5, "temperature": 25},
    "DN2S": {"value": 1e5, "temperature": 25},
    "DN33L": {"value": 1e5, "temperature": 25}
}
datasheet_per_material = {
    "DMR25": "https://www.dianyuan.com/upload/community/2017/09/05/1504577498-77035.pdf",
    "DMR24": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR44": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR40B": "https://www.dianyuan.com/upload/community/2017/09/05/1504577498-77035.pdf",
    "DMR95": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR40": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR55": "https://ferrite.ru/upload/docs/pdf/DMR55%20Material%20Datasheet.pdf",
    "DMR53": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR28": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR96": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR96A": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR47": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR90": "https://ferrite.ru/upload/docs/pdf/DMR90%20Material%20Datasheet.pdf",
    "DMR91": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR71": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR73": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR70": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "R5K": "https://ferrite.ru/upload/docs/pdf/R%205K%20Material%20Datasheet.pdf",
    "R7K": "https://ferrite.ru/upload/docs/pdf/R%207K%20Material%20Datasheet.pdf",
    "R7KC": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240817/R7KC%20Material%20Characteristics.pdf",
    "R10K": "https://ferrite.ru/upload/docs/pdf/R10K%20Material%20Datasheet.pdf",
    "R12K": "https://ferrite.ru/upload/docs/pdf/R12K%20Material%20Datasheet.pdf",
    "R12KZ": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20230401/R12KZ%20Material%20Characteristics.pdf",
    "R15K": "https://ferrite.ru/upload/docs/pdf/R15K%20Material%20Datasheet.pdf",
    "R15KZ": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20230401/R15KZ%20Material%20Characteristics.pdf",
    "R5KC": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20230401/R5KC%20Material%20Characteristics.pdf",
    "R5KZ": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20230401/R5KZ%20Material%20Characteristics.pdf",
    "R10KZ": "https://www.dianyuan.com/upload/community/2017/09/05/1504577498-77035.pdf",
    "R10KC": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240229/R10KC%20Material%20Characteristics.pdf",
    "DMR50": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR50B": "https://ferrite-group.com/upload/docs/pdf/DMR50B%20Material%20Datasheet.pdf",
    "DMR51": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR51W": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DMR52": "https://www.dianyuan.com/upload/community/2017/09/05/1504577498-77035.pdf",
    "DMR52W": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%94%B0%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN30B": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN40B": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN50B": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN85H": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN150H": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN100H": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN15P": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN200L": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN20F": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN2S": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf",
    "DN33L": "https://dongyangdongci.oss-cn-hangzhou.aliyuncs.com/uploads/20240625/%E9%95%8D%E9%94%8C%E7%9B%AE%E5%BD%95.pdf"
}
composition_per_material = {
    "DMR25": "MnZn",
    "DMR24": "MnZn",
    "DMR44": "MnZn",
    "DMR40B": "MnZn",
    "DMR95": "MnZn",
    "DMR40": "MnZn",
    "DMR55": "MnZn",
    "DMR53": "MnZn",
    "DMR28": "MnZn",
    "DMR96": "MnZn",
    "DMR96A": "MnZn",
    "DMR47": "MnZn",
    "DMR90": "MnZn",
    "DMR91": "MnZn",
    "DMR71": "MnZn",
    "DMR73": "MnZn",
    "DMR70": "MnZn",
    "R5K": "MnZn",
    "R7K": "MnZn",
    "R7KC": "MnZn",
    "R10K": "MnZn",
    "R12K": "MnZn",
    "R12KZ": "MnZn",
    "R15K": "MnZn",
    "R15KZ": "MnZn",
    "R5KC": "MnZn",
    "R5KZ": "MnZn",
    "R10KZ": "MnZn",
    "R10KC": "MnZn",
    "DMR50": "MnZn",
    "DMR50B": "MnZn",
    "DMR51": "MnZn",
    "DMR51W": "MnZn",
    "DMR52": "MnZn",
    "DMR52W": "MnZn",
    "DN30B": "NiZn",
    "DN40B": "NiZn",
    "DN50B": "NiZn",
    "DN85H": "NiZn",
    "DN150H": "NiZn",
    "DN100H": "NiZn",
    "DN15P": "NiZn",
    "DN200L": "NiZn",
    "DN20F": "NiZn",
    "DN2S": "NiZn",
    "DN33L": "NiZn"
}
family_per_material = {
    "DMR25": "DMR",
    "DMR24": "DMR",
    "DMR44": "DMR",
    "DMR40B": "DMR",
    "DMR95": "DMR",
    "DMR40": "DMR",
    "DMR55": "DMR",
    "DMR53": "DMR",
    "DMR28": "DMR",
    "DMR96": "DMR",
    "DMR96A": "DMR",
    "DMR47": "DMR",
    "DMR90": "DMR",
    "DMR91": "DMR",
    "DMR71": "DMR",
    "DMR73": "DMR",
    "DMR70": "DMR",
    "R5K": "R",
    "R7K": "R",
    "R7KC": "R",
    "R10K": "R",
    "R12K": "R",
    "R12KZ": "R",
    "R15K": "R",
    "R15KZ": "R",
    "R5KC": "R",
    "R5KZ": "R",
    "R10KZ": "R",
    "R10KC": "R",
    "DMR50": "DMR",
    "DMR50B": "DMR",
    "DMR51": "DMR",
    "DMR51W": "DMR",
    "DMR52": "DMR",
    "DMR52W": "DMR",
    "DN30B": "DN",
    "DN40B": "DN",
    "DN50B": "DN",
    "DN85H": "DN",
    "DN150H": "DN",
    "DN100H": "DN",
    "DN15P": "DN",
    "DN200L": "DN",
    "DN20F": "DN",
    "DN2S": "DN",
    "DN33L": "DN"
}

steinmetz_per_material = {
    'DMR24': {'method': 'steinmetz', 'ranges': [{'k': 18222.189549990566, 'alpha': 1.2078667042670592, 'beta': 2.631213095612934, 'ct0': 0.003085244732926972, 'ct1': 4.731257601745215e-05, 'ct2': 2.617576127348469e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 2308.205325215401, 'alpha': 1.198069885889299, 'beta': 2.612482857172042, 'ct0': 0.003999469237884526, 'ct1': -4.787642390607381e-05, 'ct2': -6.290513782395236e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR44': {'method': 'steinmetz', 'ranges': [{'k': 672342.7136688298, 'alpha': 0.6045036952560073, 'beta': 1.180460752587825, 'ct0': 0.00526346845291854, 'ct1': 4.786350472717552e-05, 'ct2': 2.3191564579655193e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 2.3469085121989113, 'alpha': 1.9371522264426937, 'beta': 2.375018085259112, 'ct0': 0.0022800504586746786, 'ct1': 2.5762843289932752e-05, 'ct2': 1.2067052204075207e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR95': {'method': 'steinmetz', 'ranges': [{'k': 28.79373330836631, 'alpha': 1.7003980839654613, 'beta': 2.6098910080727107, 'ct0': 0.002742977807908098, 'ct1': 2.43102348747002e-05, 'ct2': 1.7366128167761378e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 0.327586743964156, 'alpha': 2.211731482133195, 'beta': 2.386371618278419, 'ct0': 0.00021986094722014488, 'ct1': 1.117997769692897e-06, 'ct2': 9.271697571054527e-09, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR40': {'method': 'steinmetz', 'ranges': [{'k': 1939.032361765282, 'alpha': 1.2470575389896865, 'beta': 1.6210763411370326, 'ct0': 0.00307821370791562, 'ct1': 4.092994692837942e-05, 'ct2': 2.562809644259293e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 1.3134961169303976, 'alpha': 1.7671808483603777, 'beta': 1.6605738116141604, 'ct0': -0.03440582679244905, 'ct1': -0.001999842970537152, 'ct2': -1.601972490688133e-05, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR55': {'method': 'steinmetz', 'ranges': [{'k': 90006.01819883977, 'alpha': 1.0822820868638068, 'beta': 2.6325925882621104, 'ct0': 0.0046196297112837586, 'ct1': 7.90647303972785e-05, 'ct2': 4.35719654741543e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 1.4490883639950307, 'alpha': 2.0052756814779404, 'beta': 2.5039051983165814, 'ct0': 0.0019308688213245419, 'ct1': 3.192206708985836e-05, 'ct2': 1.8705797769971335e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR28': {'method': 'steinmetz', 'ranges': [{'k': 4577.167505840441, 'alpha': 1.3061601919906327, 'beta': 2.3138888899118335, 'ct0': 0.0027911540958925367, 'ct1': 3.619756508399598e-05, 'ct2': 4.32995714928181e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 114.02582767254282, 'alpha': 1.6004698763622307, 'beta': 2.342435699739984, 'ct0': 0.0030677885800925397, 'ct1': 2.7310194843938795e-05, 'ct2': 4.058881069315225e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR96': {'method': 'steinmetz', 'ranges': [{'k': 39.63738534426162, 'alpha': 1.706614520198055, 'beta': 2.7733080481525985, 'ct0': 0.002171483382186901, 'ct1': 1.2864457021571814e-05, 'ct2': 7.98100363961752e-08, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 2.921644098615229, 'alpha': 1.8992314257280511, 'beta': 2.6578370482807663, 'ct0': 0.002060941824644514, 'ct1': 1.0296910244550779e-05, 'ct2': 7.164801725690915e-08, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR47': {'method': 'steinmetz', 'ranges': [{'k': 1265.7264291943054, 'alpha': 1.3914720217093572, 'beta': 2.4510343446194076, 'ct0': 0.004403080871588602, 'ct1': 6.100860404865853e-05, 'ct2': 2.9735508934239955e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 0.16589227801164724, 'alpha': 1.9879485267646335, 'beta': 2.389358074533084, 'ct0': 0.00490220301585181, 'ct1': -0.000325186948134085, 'ct2': -3.033443188260543e-06, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR90': {'method': 'steinmetz', 'ranges': [{'k': 6062.518317525487, 'alpha': 1.289522129291737, 'beta': 2.463225848408546, 'ct0': 0.003901920460125466, 'ct1': 6.013236097573435e-05, 'ct2': 3.1451514398634154e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 0.5698515444439716, 'alpha': 2.0606428245982857, 'beta': 2.1517912438127005, 'ct0': 0.001129329412475669, 'ct1': 1.177696204906148e-05, 'ct2': 6.260646641931356e-08, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR50': {'method': 'steinmetz', 'ranges': [{'k': 2757436.9211954162, 'alpha': 0.9214787040404904, 'beta': 3.1004578689025353, 'ct0': 0.004091123489122735, 'ct1': 7.62732247059447e-05, 'ct2': 4.2722166323481456e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 2.8300044301996926e-05, 'alpha': 2.8610836012360448, 'beta': 2.876557619967938, 'ct0': 0.0025492088346027075, 'ct1': 4.426013094191998e-05, 'ct2': 2.756948381597386e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR50B': {'method': 'steinmetz', 'ranges': [{'k': 58.202282091181296, 'alpha': 1.7876371187097215, 'beta': 2.917503695440987, 'ct0': 0.003033816311425288, 'ct1': 4.982066369830264e-05, 'ct2': 3.743364861787002e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 0.009345996873610714, 'alpha': 2.4372334045278548, 'beta': 2.7082928806995588, 'ct0': 0.0009331689964818556, 'ct1': 9.423951431326044e-06, 'ct2': 8.193125858274262e-08, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR51': {'method': 'steinmetz', 'ranges': [{'k': 0.01677113713555139, 'alpha': 2.5939500528436947, 'beta': 2.758321187159348, 'ct0': 1.8292305680339463e-05, 'ct1': 6.26588744624651e-08, 'ct2': 4.4883444454514356e-10, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR52': {'method': 'steinmetz', 'ranges': [{'k': 10.695366805658166, 'alpha': 1.373586322968618, 'beta': 1.1294693806286886, 'ct0': 0.0015672109714944295, 'ct1': 1.5114611012671956e-05, 'ct2': 1.0436107289384965e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR96A': {'method': 'steinmetz', 'ranges': [{'k': 354.8207190771065, 'alpha': 1.4604757414198395, 'beta': 2.73396757634519, 'ct0': 0.003376393966725765, 'ct1': 2.4967499544262347e-05, 'ct2': 2.3465202591092188e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 1.8506920796211939, 'alpha': 1.9384031094961522, 'beta': 2.720938268422924, 'ct0': 0.0018130647660418902, 'ct1': 1.8755025506472198e-05, 'ct2': 1.8654276786972097e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR91': {'method': 'steinmetz', 'ranges': [{'k': 2064.8951559988477, 'alpha': 1.4471480975088247, 'beta': 2.5359939450877023, 'ct0': 0.002736459120139975, 'ct1': 4.687512621105782e-05, 'ct2': 2.459247927612495e-07, 'minimumFrequency': 1, 'maximumFrequency': 200000}, {'k': 2.046536774814736, 'alpha': 2.004193358537383, 'beta': 2.4225539279125554, 'ct0': 0.0019258503038307054, 'ct1': 2.9163277457588282e-05, 'ct2': 1.5359233953055897e-07, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR51W': {'method': 'steinmetz', 'ranges': [{'k': 0.011270519819375413, 'alpha': 2.5311150815877412, 'beta': 2.616569883974173, 'ct0': 1.582148018264116e-05, 'ct1': -9.590127894628664e-09, 'ct2': -1.3616551194901908e-11, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DMR52W': {'method': 'steinmetz', 'ranges': [{'k': 0.0036333269853595887, 'alpha': 2.347911810283285, 'beta': 2.2861234329892395, 'ct0': 0.00020815356043378058, 'ct1': 2.6944784888857467e-06, 'ct2': 2.0304347652180706e-08, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
    'DN15P': {'method': 'steinmetz', 'ranges': [{'k': 2.857761847195558e-05, 'alpha': 3.0254109120093093, 'beta': 3.6064442814904973, 'ct0': -0.0002845890155095261, 'ct1': -1.6241003833538758e-05, 'ct2': -6.738746905854904e-08, 'minimumFrequency': 200000, 'maximumFrequency': 600000000}]},
}


def add_material(material):
    mas_datum = {
        "type": "commercial",
        "name": material,
        "family": family_per_material[material],
        "resistivity": [resistivity_per_material[material]],
        "remanence": [],
        "coerciveForce": [],
        "density": density_per_material[material],
        "curieTemperature": curie_temperature_per_material[material],
        "material": "ferrite",
        "materialComposition": composition_per_material[material],
        "manufacturerInfo": {"name": "DMEGC", "datasheetUrl": datasheet_per_material[material]},
        "permeability": {
            "initial": [],
            "complex": {
                "real": [],
                "imaginary": []
            }
        },
        "saturation": [],
        "volumetricLosses": {"default": [{"method": "roshen"}]}
    }
    mas_advanced_datum = {
        "name": material,
        "manufacturerInfo": {"name": "DMEGC"},
        "permeability": {
            "amplitude": []
        },
        "volumetricLosses": {"default": []},
        "bhCycle": []
    }
    mas_materials[material] = mas_datum
    mas_advanced_materials[material] = mas_advanced_datum


permeability_vs_temperature_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/dmegc uT - uT.csv"

data = pandas.read_csv(permeability_vs_temperature_spreadsheet_path)
materials = data.columns[1:]
# print(data)
# print(materials)

mas_materials = {}
mas_advanced_materials = {}

for material in materials:
    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue        

    material_data = data[["temperature", material]].dropna()
    for index, row in material_data.iterrows():
        mas_materials[material]["permeability"]["initial"].append({"magneticFieldDcBias": 0, "temperature": float(row["temperature"]), "value": float(row[material])})


permeability_vs_temperature_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/DMEGC - u vs t.csv"

data = pandas.read_csv(permeability_vs_temperature_spreadsheet_path)
for index in range(len(data.columns)):
    if index % 2 != 0:
        continue
    if index + 1 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1]].dropna()
    material = material_data.columns[0]
    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    for index, row in material_data.iterrows():
        mas_materials[material]["permeability"]["initial"].append({"temperature": float(row.iloc[0]), "value": float(row.iloc[1])})

complex_permeability_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/dmegc u1u2f - u1u2f.csv"

data = pandas.read_csv(complex_permeability_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 3 != 0:
        continue
    if index + 2 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1, index + 2]].dropna()
    material = material_data.columns[1]
    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    for index, row in material_data.iterrows():
        # print(row.iloc[2])
        mas_materials[material]["permeability"]["complex"]["real"].append({"frequency": float(row.iloc[0]), "value": float(row.iloc[1])})
        mas_materials[material]["permeability"]["complex"]["imaginary"].append({"frequency": float(row.iloc[0]), "value": float(row.iloc[2])})

complex_permeability_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/DMEGC - complex.csv"

data = pandas.read_csv(complex_permeability_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 4 != 0:
        continue
    if index + 3 >= len(data.columns):
        break
    imag_material_data = data.iloc[:, [index, index + 1]].dropna()
    real_material_data = data.iloc[:, [index + 2, index + 3]].dropna()
    material = imag_material_data.columns[0]
    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    for index, row in imag_material_data.iterrows():
        mas_materials[material]["permeability"]["complex"]["imaginary"].append({"frequency": float(row.iloc[0]), "value": float(row.iloc[1])})
    for index, row in real_material_data.iterrows():
        mas_materials[material]["permeability"]["complex"]["real"].append({"frequency": float(row.iloc[0]), "value": float(row.iloc[1])})

bh_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/dmegc BH - BH.csv"

data = pandas.read_csv(bh_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 2 != 0:
        continue
    if index + 1 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1]].dropna()
    material = material_data.columns[0].split(".")[0]
    temperature = float(material_data.columns[1].split(".")[0])

    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    current_closest_h_to_0 = 999999
    current_closest_b_to_0 = 999999
    max_b = 0
    max_h = 0
    remanence = None
    coercivity = None
    for index, row in material_data.iterrows():
        h = float(row.iloc[0])
        b = float(row.iloc[1]) / 1000
        if abs(h) < current_closest_h_to_0:
            current_closest_h_to_0 = abs(h)
            remanence = abs(b)
        if abs(b) < current_closest_b_to_0:
            current_closest_b_to_0 = abs(b)
            coercivity = abs(h)
        if abs(b) > max_b:
            max_b = abs(b)
            max_h = abs(h)

        mas_advanced_materials[material]["bhCycle"].append({"magneticFluxDensity": float(row.iloc[1]) / 1000, "magneticField": float(row.iloc[0]), "temperature": temperature})

    mas_materials[material]["saturation"].append({
        "magneticFluxDensity": max_b,
        "magneticField": max_h,
        "temperature": temperature
    })
    mas_materials[material]["coerciveForce"].append({
        "magneticFluxDensity": 0,
        "magneticField": coercivity,
        "temperature": temperature
    })

    mas_materials[material]["remanence"].append({
        "magneticFluxDensity": remanence,
        "magneticField": 0,
        "temperature": temperature
    })

amplitude_permeability_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/dmegc uaB - uaB.csv"

data = pandas.read_csv(amplitude_permeability_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 2 != 0:
        continue
    if index + 1 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1]].dropna()
    material = material_data.columns[0].split(".")[0]
    temperature = float(material_data.columns[1].split(".")[0])

    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    for index, row in material_data.iterrows():
        mas_advanced_materials[material]["permeability"]["amplitude"].append({"value": float(row.iloc[1]), "temperature": temperature, "magneticFluxDensityPeak": float(row.iloc[0]) / 1000})

permeability_vs_bias_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/dmegc uHdc - uHdc.csv"

data = pandas.read_csv(permeability_vs_bias_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 2 != 0:
        continue
    if index + 1 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1]].dropna()
    material = material_data.columns[0].split(".")[0]
    temperature = float(material_data.columns[1].split(".")[0])

    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    for index, row in material_data.iterrows():
        mas_materials[material]["permeability"]["initial"].append({"magneticFieldDcBias": float(row.iloc[0]), "temperature": temperature, "value": float(row.iloc[1])})

losses_vs_temperature_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/DMEGC - loss vs T.csv"

data = pandas.read_csv(losses_vs_temperature_spreadsheet_path)
material_column_ranges = {}
for index in range(len(data.columns)):
    if not data.columns[index].startswith("Unnam"):
        material = data.columns[index]
        material_column_ranges[material] = {
            "start": index
        }
    else:
        material_column_ranges[material]["end"] = index

for material, material_indexes in material_column_ranges.items():
    for index in range(material_indexes["start"], material_indexes["end"] + 1, 2):
        material_data = data.iloc[:, [index, index + 1]].dropna()
        material_data = material_data.rename(columns=material_data.iloc[0])
        material_data = material_data.drop(material_data.index[0])
        frequency = None
        magnetic_flux_density = None
        temperature = None
        for column in material_data.columns:
            if "kHz" in column:
                frequency = float(column.split("kHz")[0]) * 1000
            elif "KHz" in column:
                frequency = float(column.split("KHz")[0]) * 1000
            elif "MHz" in column:
                frequency = float(column.split("MHz")[0]) * 1000000
            elif "z" in column:
                print(column)
                assert 0

            if "mT" in column:
                magnetic_flux_density = float(column.split("mT")[0]) / 1000
            elif "T" in column:
                print(column)
                assert 0
        assert (frequency is not None and magnetic_flux_density is not None) or (frequency is not None and temperature is not None) or (temperature is not None and magnetic_flux_density is not None)
        assert frequency is None or frequency > 1000
        if material not in mas_materials.keys():
            try:
                add_material(material)
            except KeyError:
                print(material)
                assert 0
        if len(mas_advanced_materials[material]["volumetricLosses"]["default"]) == 0:
            mas_advanced_materials[material]["volumetricLosses"]["default"].append([])
        for index, row in material_data.iterrows():
            assert frequency if frequency is not None else float(row.iloc[0]) > 1000, f"Material: {material}, frequency: {frequency if frequency is not None else float(row.iloc[0])}"
            mas_advanced_materials[material]["volumetricLosses"]["default"][0].append(
                {
                    "magneticFluxDensity": {
                        "frequency": frequency if frequency is not None else float(row.iloc[0]),
                        "magneticFluxDensity": {
                            "processed": {
                                "label": "Sinusoidal",
                                "peak": magnetic_flux_density if magnetic_flux_density is not None else float(row.iloc[0]),
                                "offset": 0
                            }
                        }
                    },
                    "temperature": temperature if temperature is not None else float(row.iloc[0]),
                    "value": float(row.iloc[1]),
                    "origin": "manufacturer"
                }
            )

losses_vs_frequency_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/DMEGC - loss vs f.csv"

data = pandas.read_csv(losses_vs_frequency_spreadsheet_path)
material_column_ranges = {}
for index in range(len(data.columns)):
    if not data.columns[index].startswith("Unnam"):
        material = data.columns[index]
        material_column_ranges[material] = {
            "start": index
        }
    else:
        material_column_ranges[material]["end"] = index

for material, material_indexes in material_column_ranges.items():
    for index in range(material_indexes["start"], material_indexes["end"] + 1, 2):
        material_data = data.iloc[:, [index, index + 1]].dropna()
        material_data = material_data.rename(columns=material_data.iloc[0])
        material_data = material_data.drop(material_data.index[0])
        frequency = None
        magnetic_flux_density = None
        temperature = None
        for column in material_data.columns:
            column = str(column)
            if "kHz" in column:
                frequency = float(column.split("kHz")[0]) * 1000
                continue
            elif "KHz" in column:
                frequency = float(column.split("KHz")[0]) * 1000
                continue
            elif "MHz" in column:
                frequency = float(column.split("MHz")[0]) * 1000000
                continue
            elif "z" in column:
                print(column)
                assert 0

            if "mT" in column:
                magnetic_flux_density = float(column.split("mT")[0]) / 1000
                continue
            elif "T" in column:
                print(column)
                assert 0

            temperature = float(column)

        assert (frequency is not None and magnetic_flux_density is not None) or (frequency is not None and temperature is not None) or (temperature is not None and magnetic_flux_density is not None)
        assert frequency is None or frequency > 1000
        if material not in mas_materials.keys():
            try:
                add_material(material)
            except KeyError:
                print(material)
                assert 0
        if len(mas_advanced_materials[material]["volumetricLosses"]["default"]) == 0:
            mas_advanced_materials[material]["volumetricLosses"]["default"].append([])
        for index, row in material_data.iterrows():
            assert frequency if frequency is not None else float(row.iloc[0]) > 1000, f"Material: {material}, frequency: {frequency if frequency is not None else float(row.iloc[0])}"
            mas_advanced_materials[material]["volumetricLosses"]["default"][0].append(
                {
                    "magneticFluxDensity": {
                        "frequency": frequency if frequency is not None else float(row.iloc[0]),
                        "magneticFluxDensity": {
                            "processed": {
                                "label": "Sinusoidal",
                                "peak": magnetic_flux_density if magnetic_flux_density is not None else float(row.iloc[0]),
                                "offset": 0
                            }
                        }
                    },
                    "temperature": temperature if temperature is not None else float(row.iloc[0]),
                    "value": float(row.iloc[1]),
                    "origin": "manufacturer"
                }
            )
permeability_vs_bias_spreadsheet_path = "/home/alf/OpenMagnetics/MKF/src/tools/temp/DMEGC - loss factor.csv"

data = pandas.read_csv(permeability_vs_bias_spreadsheet_path)

for index in range(len(data.columns)):
    if index % 2 != 0:
        continue
    if index + 1 >= len(data.columns):
        break
    material_data = data.iloc[:, [index, index + 1]].dropna()
    material = material_data.columns[0].split(".")[0]

    if material not in mas_materials.keys():
        try:
            add_material(material)
        except KeyError:
            continue  

    mas_materials[material]["volumetricLosses"]["default"].append(
        {
            "method": "lossFactor",
            "factors": []
        })

    for index, row in material_data.iterrows():
        mas_materials[material]["volumetricLosses"]["default"][-1]["factors"].append({"frequency": float(row.iloc[0]), "value": float(row.iloc[1])})


for material in remanence_per_material.keys():
    mas_materials[material]["remanence"].extend(remanence_per_material[material])

for material in coercivity_per_material.keys():
    mas_materials[material]["coerciveForce"].extend(coercivity_per_material[material])

for material in saturation_per_material.keys():
    mas_materials[material]["saturation"].extend(saturation_per_material[material])

for material in steinmetz_per_material.keys():
    mas_materials[material]["volumetricLosses"]["default"].append(steinmetz_per_material[material])

for material in datasheet_per_material.keys():
    if material not in mas_materials.keys():
        print(material)

for material in mas_materials.keys():
    if len(mas_materials[material]["permeability"]["complex"]["imaginary"]) == 0:
        del mas_materials[material]["permeability"]["complex"]
    if len(mas_materials[material]["coerciveForce"]) == 0:
        del mas_materials[material]["coerciveForce"]
    if len(mas_materials[material]["remanence"]) == 0:
        del mas_materials[material]["remanence"]

ndjson.dump(mas_materials.values(), open(f"{pathlib.Path(__file__).parent.resolve()}/dmegc_core_materials.ndjson", "w"), ensure_ascii=False)
ndjson.dump(mas_advanced_materials.values(), open(f"{pathlib.Path(__file__).parent.resolve()}/dmegc_advanced_core_materials.ndjson", "w"), ensure_ascii=False)
