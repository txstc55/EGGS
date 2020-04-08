import os
import argparse
import json

matrix_sizes = [100, 1000, 10000, 100000, 1000000]

if not os.path.isdir("build"):
    print("Creating build folder now")
    os.system("mkdir build")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")
else:
    print("Build exists, rebuilding")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")

os.system("cd build && rm *_result.txt && cd ..")

for i in matrix_sizes:
    for j in range(5):
        print("Matrix size "+str(i))
        os.system("cd build && ./tutorial/803_Numeric_bin -d " +
                  str(j) + " -r "+str(i) + " && cd ..")

os.system("mv build/*_result.txt result_datas/")


numeric_multi_result = {}
numeric_single_result = {}
mkl_multi_result = {}
mkl_single_result = {}

mkl_single_result["method_name"] = "MKL SINGLE THREAD"
mkl_multi_result["method_name"] = "MKL MULTI THREAD"
numeric_single_result["method_name"] = "OURS SINGLE THREAD"
numeric_multi_result["method_name"] = "OURS MULTI THREAD"

mkl_single_result["time"] = []
mkl_multi_result["time"] = []
numeric_single_result["time"] = []
numeric_multi_result["time"] = []

f = open("result_datas/pipe_result.txt", 'r')
count = 0
for line in f:
    time = float(line.split(": ")[-1].split()[0])
    if count % 4 == 0:
        mkl_multi_result["time"].append(time)
    elif count % 4 == 1:
        mkl_single_result["time"].append(time)
    elif count % 4 == 2:
        numeric_multi_result["time"].append(time)
    elif count % 4 == 3:
        numeric_single_result["time"].append(time)
    count += 1

with open("result_datas/pipe.json", 'w')as j:
    json.dump([mkl_multi_result, mkl_single_result, numeric_multi_result, numeric_single_result], j)


mkl_single_result["time"] = []
mkl_multi_result["time"] = []
numeric_single_result["time"] = []
numeric_multi_result["time"] = []
f = open("result_datas/const_result.txt", 'r')
count = 0
for line in f:
    time = float(line.split(": ")[-1].split()[0])
    if count % 4 == 0:
        mkl_multi_result["time"].append(time)
    elif count % 4 == 1:
        mkl_single_result["time"].append(time)
    elif count % 4 == 2:
        numeric_multi_result["time"].append(time)
    elif count % 4 == 3:
        numeric_single_result["time"].append(time)
    count += 1

with open("result_datas/const.json", 'w')as j:
    json.dump([mkl_multi_result, mkl_single_result, numeric_multi_result, numeric_single_result], j)

mkl_single_result["time"] = []
mkl_multi_result["time"] = []
numeric_single_result["time"] = []
numeric_multi_result["time"] = []
f = open("result_datas/sypr_result.txt", 'r')
count = 0
for line in f:
    time = float(line.split(": ")[-1].split()[0])
    if count % 4 == 0:
        mkl_multi_result["time"].append(time)
    elif count % 4 == 1:
        mkl_single_result["time"].append(time)
    elif count % 4 == 2:
        numeric_multi_result["time"].append(time)
    elif count % 4 == 3:
        numeric_single_result["time"].append(time)
    count += 1

with open("result_datas/sypr.json", 'w')as j:
    json.dump([mkl_multi_result, mkl_single_result, numeric_multi_result, numeric_single_result], j)


mkl_single_result["time"] = []
mkl_multi_result["time"] = []
numeric_single_result["time"] = []
numeric_multi_result["time"] = []
f = open("result_datas/syrk_result.txt", 'r')
count = 0
for line in f:
    time = float(line.split(": ")[-1].split()[0])
    if count % 4 == 0:
        mkl_multi_result["time"].append(time)
    elif count % 4 == 1:
        mkl_single_result["time"].append(time)
    elif count % 4 == 2:
        numeric_multi_result["time"].append(time)
    elif count % 4 == 3:
        numeric_single_result["time"].append(time)
    count += 1

with open("result_datas/syrk.json", 'w')as j:
    json.dump([mkl_multi_result, mkl_single_result, numeric_multi_result, numeric_single_result], j)


mkl_single_result["time"] = []
mkl_multi_result["time"] = []
numeric_single_result["time"] = []
numeric_multi_result["time"] = []
f = open("result_datas/spmv_result.txt", 'r')
count = 0
for line in f:
    time = float(line.split(": ")[-1].split()[0])
    if count % 2 == 0:
        mkl_multi_result["time"].append(time)
    elif count % 2 == 1:
        numeric_multi_result["time"].append(time)
    count += 1

with open("result_datas/spmv.json", 'w')as j:
    json.dump([mkl_multi_result, numeric_multi_result], j)