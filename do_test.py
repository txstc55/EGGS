import os
import argparse
import json

parser = argparse.ArgumentParser(description='Process demo type and file name')
parser.add_argument("--d", default=1, help="Demo type, from 1 to 3")
parser.add_argument(
    "--f", default="", help="input mesh name, has to be in tutorial/data folder already")
parser.add_argument(
    "--e", default=0, help="Example, 0 for slim, 1 for cotangent smoothing, 2 for optical flow")
parser.add_argument("--w", default="1000000",
                    help="Example, 0 for slim, 1 for cotangent smoothing")
parser.add_argument("--o", default="",
                    help="output json file name")
args = parser.parse_args()
file = args.f
example = int(args.e)
demo_type = int(args.d)
weight = str(int(args.w))
output = str(args.o)

if (file != ""):
    file = " -f "+file


demo_type = " -d " + str(demo_type)
# check if build exists
if not os.path.isdir("build"):
    print("Creating build folder now")
    os.system("mkdir build")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")
else:
    print("Build exists, rebuilding")
    os.system("cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j && cd ..")

# check if result data folder exists
if not os.path.isdir("result_datas"):
    print("Creating result datas folder now")
    os.system("mkdir result_datas")


eigen_data = {}
eigen_data["name"] = "EIGEN"
eigen_data["COMPUTE"] = 0
eigen_data["SOLVE"] = 0
eigen_data["ASSEMBLE"] = 0

mkl_data = {}
mkl_data["name"] = "MKL"
mkl_data["COMPUTE"] = 0
mkl_data["SOLVE"] = 0
mkl_data["ASSEMBLE"] = 0


numeric1_data = {}
numeric1_data["name"] = "OURS"
numeric1_data["COMPUTE"] = 0
numeric1_data["SOLVE"] = 0
numeric1_data["ASSEMBLE"] = 0

numeric2_data = {}
numeric2_data["name"] = "NUMERIC2"
numeric2_data["COMPUTE"] = 0
numeric2_data["SOLVE"] = 0
numeric2_data["ASSEMBLE"] = 0

if example == 0:
    # do eigen tests
    os.system("cd build && ./tutorial/709_SLIM_bin -m 0" +
              demo_type + file + " && cd ..")
    os.system(
        "mv build/result.txt result_datas/result_eigen.txt && mv build/slim.obj result_datas/slim_eigen.obj")

    # do mkl tests
    os.system("cd build && ./tutorial/709_SLIM_bin -m 2" +
              demo_type + file + " && cd ..")
    os.system(
        "mv build/result.txt result_datas/result_mkl.txt && mv build/slim.obj result_datas/slim_mkl.obj")

    # do numeric 1 tests
    os.system("cd build && ./tutorial/709_SLIM_bin -m 3" +
              demo_type + file + " && cd ..")
    os.system("mv build/result.txt result_datas/result_numeric1.txt && mv build/slim.obj result_datas/slim_numeric1.obj")

    # # do numeric 2 tests
    # os.system("cd build && ./tutorial/709_SLIM_bin -m 4" +
    #           demo_type + file + " && cd ..")
    # os.system("mv build/result.txt result_datas/result_numeric2.txt && mv build/slim.obj result_datas/slim_numeric2.obj")

    f = open("result_datas/result_eigen.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 11
            if count == 0 or count == 1 or count == 2 or count == 3 or count == 5 or count == 6 or count == 8:
                eigen_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 4 or count == 7:
                eigen_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                eigen_data["SOLVE"] += float(line.split(": ")[1])
            count += 1


    f = open("result_datas/result_mkl.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 14
            if (count >= 0 and count <= 8) or count == 11:
                mkl_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 9 or count == 10:
                mkl_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                mkl_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    f = open("result_datas/result_numeric1.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 9
            if (count >= 0 and count <= 5):
                numeric1_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 6:
                numeric1_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                numeric1_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    # f = open("result_datas/result_numeric2.txt")
    # count = 0
    # for line in f:
    #     if not line.startswith("START"):
    #         count = count % 10
    #         if (count >= 0 and count <= 3) or count == 5:
    #             numeric2_data["ASSEMBLE"] += float(line.split(": ")[1])
    #         elif count == 4 or count == 6 or count == 7:
    #             numeric2_data["COMPUTE"] += float(line.split(": ")[1])
    #         else:
    #             numeric2_data["SOLVE"] += float(line.split(": ")[1])
    #         count += 1

    if output == "":
        with open('result_datas/all_result.json', 'w') as j:
            json.dump([eigen_data, mkl_data, numeric1_data], j)
    else:
        with open("result_datas/"+output+".json", 'w')as j:
            json.dump([eigen_data, mkl_data, numeric1_data], j)

elif example == 1:
    # do eigen tests
    os.system("cd build && ./tutorial/205_Laplacian_bin -w " +
              weight+" -m 0" + file + " && cd ..")
    os.system("mv build/result_cot.txt result_datas/result_cot_eigen.txt && mv build/cot_smoothed.obj result_datas/cot_smoothed_eigen.obj")

    # do mkl tests
    os.system("cd build && ./tutorial/205_Laplacian_bin -w " +
              weight+" -m 2" + file + " && cd ..")
    os.system("mv build/result_cot.txt result_datas/result_cot_mkl.txt && mv build/cot_smoothed.obj result_datas/cot_smoothed_mkl.obj")
    # do numeric 1 tests
    os.system("cd build && ./tutorial/205_Laplacian_bin -w " +
              weight+" -m 3" + file + " && cd ..")
    os.system("mv build/result_cot.txt result_datas/result_cot_numeric1.txt && mv build/cot_smoothed.obj result_datas/cot_smoothed_numeric1.obj")

    # # do numeric 2 tests
    # os.system("cd build && ./tutorial/205_Laplacian_bin -w " +
    #           weight+" -m 4" + file + " && cd ..")
    # os.system("mv build/result_cot.txt result_datas/result_cot_numeric2.txt && mv build/cot_smoothed.obj result_datas/cot_smoothed_numeric2.obj")

    f = open("result_datas/result_cot_eigen.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 7
            if count <= 2 or count == 4:
                eigen_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 3:
                eigen_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                eigen_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    f = open("result_datas/result_cot_mkl.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 9
            if count <= 3 or count == 6:
                mkl_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 4 or count == 5:
                mkl_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                mkl_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    f = open("result_datas/result_cot_numeric1.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 6
            if count <= 2:
                numeric1_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 3:
                numeric1_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                numeric1_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    # f = open("result_datas/result_cot_numeric2.txt")
    # count = 0
    # for line in f:
    #     if not line.startswith("START"):
    #         count = count % 7
    #         if count <= 2:
    #             numeric2_data["ASSEMBLE"] += float(line.split(": ")[1])
    #         elif count == 4 or count == 3:
    #             numeric2_data["COMPUTE"] += float(line.split(": ")[1])
    #         else:
    #             numeric2_data["SOLVE"] += float(line.split(": ")[1])
    #         count += 1

    if output == "":
        with open('result_datas/all_result_cot.json', 'w') as j:
            json.dump([eigen_data, mkl_data, numeric1_data], j)
    else:
        with open("result_datas/"+output+".json", 'w')as j:
            json.dump([eigen_data, mkl_data, numeric1_data], j)

elif example == 2:
    # do eigen tests
    os.system("cd build && ./tutorial/801_OpticalFlow_bin -m 0" + demo_type+" && cd ..")
    os.system("mv build/result_opt.txt result_datas/result_opt_eigen.txt")
    # # do numeric multi tests
    # os.system("cd build && ./tutorial/801_OpticalFlow_bin -m 3" + " && cd ..")
    # os.system("mv build/result_opt.txt result_datas/result_opt_numeric1.txt")

    os.system("cd build && ./tutorial/801_OpticalFlow_bin -m 4" + demo_type+" && cd ..")
    os.system("mv build/result_opt.txt result_datas/result_opt_numeric1.txt")

    f = open("result_datas/result_opt_eigen.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 7
            if count == 0 or count == 3 or count == 6:
                eigen_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 1 or count == 2:
                eigen_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                eigen_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    # f = open("result_datas/result_opt_numeric1.txt")
    # count = 0
    # for line in f:
    #     if not line.startswith("START"):
    #         count = count % 6
    #         if count == 0 or count == 2 or count == 5:
    #             numeric1_data["ASSEMBLE"] += float(line.split(": ")[1])
    #         elif count == 1:
    #             numeric1_data["COMPUTE"] += float(line.split(": ")[1])
    #         else:
    #             numeric1_data["SOLVE"] += float(line.split(": ")[1])
    #         count += 1

    f = open("result_datas/result_opt_numeric1.txt")
    count = 0
    for line in f:
        if not line.startswith("START"):
            count = count % 4
            if count == 3:
                numeric1_data["ASSEMBLE"] += float(line.split(": ")[1])
            elif count == 0:
                numeric1_data["COMPUTE"] += float(line.split(": ")[1])
            else:
                numeric1_data["SOLVE"] += float(line.split(": ")[1])
            count += 1

    if output == "":
        with open('result_datas/all_result_opt.json', 'w') as j:
            json.dump([eigen_data, numeric1_data], j)
    else:
        with open("result_datas/"+output+".json", 'w')as j:
            json.dump([eigen_data, numeric1_data], j)

elif example == 3:
    os.system("cd build && ./tutorial/802_CotMatrix_bin"+file+" && cd ..")
    os.system("mv build/result_cot_matrix.txt result_datas/result_cot_matrix.txt")

    f = open("result_datas/result_cot_matrix.txt")
    count = 0
    for line in f:
        if count == 0:
            eigen_data["COMPUTE"] = float(line.split(":")[1])
        elif count == 1:
            numeric1_data["COMPUTE"] = float(line.split(":")[1])
        # elif count == 2:
        #     numeric2_data["COMPUTE"] = float(line.split(":")[1])
        else:
            print(line)
        count += 1

    if output == "":
        with open("result_datas/all_result_cot_matrix.json", 'w') as j:
            json.dump([eigen_data, numeric1_data], j)
    else:
        with open("result_datas/"+output+".json", 'w')as j:
            json.dump([eigen_data, numeric1_data], j)
