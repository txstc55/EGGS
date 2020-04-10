import os
import argparse
import json

matrix_sizes = [100, 1000, 10000, 100000, 1000000]
numeric_multi_result = {}
numeric_single_result = {}
mkl_multi_result = {}
mkl_single_result = {}
eigen_single_result = {}
mkl_single_result["method_name"] = "MKL SINGLE THREAD"
mkl_multi_result["method_name"] = "MKL MULTI THREAD"
numeric_single_result["method_name"] = "OURS SINGLE THREAD"
numeric_multi_result["method_name"] = "OURS MULTI THREAD"
eigen_single_result["method_name"] = "EIGEN SINGLE THREAD"


if not os.path.isdir("build"):
    print("Creating build folder now")
    os.system("mkdir build")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")
else:
    print("Build exists, rebuilding")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")


os.system("cd build && rm *_result.txt && cd ..")
for i in matrix_sizes:
    for j in range(3):
        print("Matrix size "+str(i))
        os.system("cd build && taskset -c 0-3 ./tutorial/803_Numeric_bin -d " +
                  str(j) + " -r "+str(i) + " && cd ..")
os.system("mv build/*_result.txt result_datas/")


def write_json(input_file, output_file):
    mkl_single_result["time"] = []
    mkl_multi_result["time"] = []
    numeric_single_result["time"] = []
    numeric_multi_result["time"] = []
    eigen_single_result["time"] = []
    numeric_multi_result["pre_comp"] = []

    f = open(input_file, 'r')
    count = 0
    for line in f:
        time = float(line.split(": ")[-1].split()[0])
        if count % 6 == 0:
            numeric_multi_result["pre_comp"].append(time)
        elif count % 6 == 1:
            mkl_multi_result["time"].append(time)
        elif count % 6 == 2:
            mkl_single_result["time"].append(time)
        elif count % 6 == 3:
            numeric_multi_result["time"].append(time)
        elif count % 6 == 4:
            numeric_single_result["time"].append(time)
        elif count % 6 == 5:
            eigen_single_result["time"].append(time)
        count += 1

    with open(output_file, 'w')as j:
        json.dump([mkl_multi_result, mkl_single_result,
                   numeric_multi_result, numeric_single_result, eigen_single_result], j)


write_json("result_datas/pipe_result.txt", "result_datas/pipe.json")
write_json("result_datas/sypr_result.txt", "result_datas/sypr.json")
write_json("result_datas/syrk_result.txt", "result_datas/syrk.json")

# do the test with 15 entries per row
os.system("cd build && rm *_result.txt && cd ..")
for i in matrix_sizes:
    for j in range(3):
        print("Matrix size "+str(i))
        os.system("cd build && taskset -c 0-3 ./tutorial/803_Numeric_bin -d " +
                  str(j) + " -r "+str(i) + " -e 15" + " && cd ..")
os.system("mv build/*_result.txt result_datas/")

write_json("result_datas/pipe_result.txt", "result_datas/pipe_15.json")
write_json("result_datas/sypr_result.txt", "result_datas/sypr_15.json")
write_json("result_datas/syrk_result.txt", "result_datas/syrk_15.json")

os.system("cd build && rm *_result.txt && cd ..")
for i in ["0", "0-1", "0-3", "0-7", "0-15"]:
    os.system("cd build && taskset -c "+i +
              " ./tutorial/803_Numeric_bin -d 1 -r 500000 -e 15 && cd ..")
os.system("mv build/*_result.txt result_datas/")
write_json("result_datas/sypr_result.txt", "result_datas/sypr_thread_diff.json")

# mkl_single_result["time"] = []
# mkl_multi_result["time"] = []
# numeric_single_result["time"] = []
# numeric_multi_result["time"] = []
# eigen_single_result["time"] = []
# numeric_multi_result["pre_comp"] = []
# f = open("result_datas/spmv_result.txt", 'r')
# count = 0
# for line in f:
#     time = float(line.split(": ")[-1].split()[0])
#     if count % 4 == 0:
#         numeric_multi_result["pre_comp"].append(time)
#     elif count % 4 == 1:
#         mkl_multi_result["time"].append(time)
#     elif count % 4 == 2:
#         numeric_multi_result["time"].append(time)
#     elif count % 4 == 3:
#         eigen_single_result["time"].append(time)
#     count += 1

# with open("result_datas/spmv.json", 'w')as j:
#     json.dump([mkl_multi_result, numeric_multi_result], j)


# mkl_single_result["time"] = []
# mkl_multi_result["time"] = []
# numeric_single_result["time"] = []
# numeric_multi_result["time"] = []
# f = open("result_datas/const_result.txt", 'r')
# count = 0
# for line in f:
#     time = float(line.split(": ")[-1].split()[0])
#     if count % 6 == 0:
#         numeric_multi_result["pre_comp"].append(time)
#     elif count % 6 == 1:
#         mkl_multi_result["time"].append(time)
#     elif count % 6 == 2:
#         mkl_single_result["time"].append(time)
#     elif count % 6 == 3:
#         numeric_multi_result["time"].append(time)
#     elif count % 6 == 4:
#         numeric_single_result["time"].append(time)
#     elif count % 6 == 5:
#         eigen_single_result["time"].append(time)
#     count += 1

# with open("result_datas/const.json", 'w')as j:
#     json.dump([mkl_multi_result, mkl_single_result,
#                numeric_multi_result, numeric_single_result, eigen_single_result], j)
