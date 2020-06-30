import json
import numpy as np
import pandas as pd


def to_table(file_name):
    f = open(file_name, 'r')
    j = json.load(f)
    # print(j)

    better_json = {}
    for item in j:
        better_json[item["method_name"]] = [
            float(x)/1000 for x in item["time"]]
        if ("pre_comp" in item):
            better_json["OURS PRE COMPUTATION"] = [
                float(x)/1000 for x in item["pre_comp"]]
    better_json["NUM ROW"] = ['%.2E' % 10**x for x in [2, 3, 4, 5, 6]]
    df = pd.DataFrame.from_dict(better_json)
    rename_dic = {"MKL MULTI THREAD": "MKL MT", "MKL SINGLE THREAD": "MKL ST", "OURS SINGLE THREAD": "OURS ST",
                  "OURS MULTI THREAD": "OURS MT", "EIGEN SINGLE THREAD": "EIGEN ST", "OURS PRE COMPUTATION": "OURS PC"}
    df = df.rename(columns=rename_dic)
    formats = {}
    df = df.set_index("NUM ROW")
    df = df[["EIGEN ST", "MKL ST", "MKL MT", "OURS PC", "OURS ST", "OURS MT"]]
    for item in better_json["NUM ROW"]:
        formats[item] = "{:0.2e}".format
    df = df.T

    la = df.to_latex(index=True, formatters=formats)
    la = la.replace('e-0', 'e-').replace('e+0', 'e')
    print(la)


to_table("result_datas/syrk.json")
to_table("result_datas/syrk_15.json")
