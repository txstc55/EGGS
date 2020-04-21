import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from math import ceil, log, log2


file_list = ["pipe", "sypr", "syrk"]
file_list_15 = [x+"_15" for x in file_list]
file_list += file_list_15
file_list_full = ["result_datas/"+x+".json" for x in file_list]

x_range = [2, 3, 4, 5, 6]
x_range = [10**x for x in x_range]

plt.rc('ytick', labelsize=28)    # fontsize of the tick labels
plt.rc('xtick', labelsize=21)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # fontsize of the tick labels


for k in range(len(file_list_full)):
    f = file_list_full[k]
    with open(f, 'r') as j:
        results = json.load(j)
        print(results)
        time_result = {}
        pre_comp_time = []
        for r in results:
            time = r['time']
            time = [float(x) for x in time]
            time_result[r['method_name']] = time
            if (r['method_name'] == "OURS MULTI THREAD"):
                pre_comp_time = r["pre_comp"]

        pre_comp_time = [float(x)/1000 for x in pre_comp_time]

        multi_thread_factor = []
        single_thread_factor = []
        eigen_single_thread_factor = []
        for i in range(5):
            multi_thread_factor.append(
                time_result["MKL MULTI THREAD"][i]/time_result["OURS MULTI THREAD"][i])
        plt.plot(x_range, multi_thread_factor,
                 linewidth=3, label="Eight Threads Speedup", color="#1B9CFC")
        if ("MKL SINGLE THREAD" in time_result):
            for i in range(5):
                single_thread_factor.append(
                    time_result["MKL SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
                eigen_single_thread_factor.append(
                    time_result["EIGEN SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
            plt.plot(x_range, single_thread_factor,
                     linewidth=3, label="Single Thread Speedup", color="#55E6C1")

        plt.legend(fontsize=20)

        # if "_15" in f:
        #     plt.yticks(list(range(0, 45, 5)))
        # else:
        plt.yticks(list(range(0, 19, 4)))
        # plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
        #                     hspace=0, wspace=0)
        ax = plt.axes()
        ax.set_xscale('log')
        # ax.xaxis.set_major_locator(ticker.FixedLocator(x_range))
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + ".pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        # plot the eigen difference
        plt.plot(x_range, eigen_single_thread_factor,
                 linewidth=3, label="Single Thread Speedup Over Eigen", color="#d63031")
        plt.yticks(list(range(0, 80, 10)))
        ax = plt.axes()
        if len(ax.lines) > 1:
            plt.legend()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.yaxis.set_major_locator(ticker.FixedLocator(list(range(0, 80, 10))))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_eigen.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        plt.plot(x_range, pre_comp_time,
                 linewidth=3, label="Pre Computation Time", color="#b8e994")
        # plt.legend(fontsize=16)
        ax = plt.axes()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_pre_comp.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()


# plot the scalability
f = "result_datas/sypr_thread_diff.json"
with open(f, 'r') as j:
    results = json.load(j)
    print(results)
    time_result = {}
    x_range = [1, 2, 4, 8, 16]
    for r in results:
        time = r['time']
        time = [float(x) for x in time]
        if (not "EIGEN SINGLE THREAD" in r['method_name']):
            time_result[r['method_name']] = time

    mkl_speedup = []
    ours_speedup = []
    for i in range(5):
        mkl_speedup.append(
            time_result["MKL SINGLE THREAD"][i]/time_result["MKL MULTI THREAD"][i])
        ours_speedup.append(
            time_result["OURS SINGLE THREAD"][i]/time_result["OURS MULTI THREAD"][i])

    plt.plot(x_range, mkl_speedup, label="MKL", color="#009688", linewidth=2)
    plt.plot(x_range, ours_speedup, label="Ours", color="#FFC107", linewidth=2)
    plt.legend(fontsize=20)
    plt.yticks(list(range(0, 17, 2)), fontsize=20)

    # plt.yticks(np.arange(0, max_y, gap))
    ax = plt.axes()
    # ax.set_xscale('log', basex=2)
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.xaxis.set_major_locator(ticker.FixedLocator([1, 2, 4, 8, 16]))
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig("test_result_graphs/spyr_thread_diff.pdf", bbox_inches='tight',
                pad_inches=0, dpi=200)
    plt.close()

f = "result_datas/sypr_5_thread_diff.json"
with open(f, 'r') as j:
    results = json.load(j)
    print(results)
    time_result = {}
    x_range = [1, 2, 4, 8, 16]
    for r in results:
        time = r['time']
        time = [float(x) for x in time]
        if (not "EIGEN SINGLE THREAD" in r['method_name']):
            time_result[r['method_name']] = time

    mkl_speedup = []
    ours_speedup = []
    for i in range(5):
        mkl_speedup.append(
            time_result["MKL SINGLE THREAD"][i]/time_result["MKL MULTI THREAD"][i])
        ours_speedup.append(
            time_result["OURS SINGLE THREAD"][i]/time_result["OURS MULTI THREAD"][i])

    plt.plot(x_range, mkl_speedup, label="MKL", color="#009688", linewidth=2)
    plt.plot(x_range, ours_speedup, label="Ours", color="#FFC107", linewidth=2)
    plt.legend(fontsize=20)
    plt.yticks(list(range(0, 17, 2)), fontsize=20)

    # plt.yticks(np.arange(0, max_y, gap))
    ax = plt.axes()
    # ax.set_xscale('log', basex=2)
    ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
    ax.xaxis.set_major_locator(ticker.FixedLocator([1, 2, 4, 8, 16]))
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig("test_result_graphs/spyr_5_thread_diff.pdf", bbox_inches='tight',
                pad_inches=0, dpi=200)
    plt.close()


file_list = ["const", "const_2", "const_abc"]
file_list_15 = [x+"_15" for x in file_list]
file_list += file_list_15
file_list_full = ["result_datas/"+x+".json" for x in file_list]

x_range = [3, 4, 5]
x_range = [10**x for x in x_range]
for k in range(len(file_list_full)):
    f = file_list_full[k]
    with open(f, 'r') as j:
        results = json.load(j)
        print(results)
        time_result = {}
        pre_comp_time = []
        for r in results:
            time = r['time']
            time = [float(x) for x in time]
            time_result[r['method_name']] = time
            if (r['method_name'] == "OURS MULTI THREAD"):
                pre_comp_time = r["pre_comp"]

        pre_comp_time = [float(x)/1000 for x in pre_comp_time]

        multi_thread_factor = []
        single_thread_factor = []
        eigen_single_thread_factor = []
        for i in range(3):
            multi_thread_factor.append(
                time_result["MKL MULTI THREAD"][i]/time_result["OURS MULTI THREAD"][i])
        plt.plot(x_range, multi_thread_factor,
                 linewidth=3, label="Eight Threads Speedup", color="#1B9CFC")
        if ("MKL SINGLE THREAD" in time_result):
            for i in range(3):
                single_thread_factor.append(
                    time_result["MKL SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
                eigen_single_thread_factor.append(
                    time_result["EIGEN SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
            plt.plot(x_range, single_thread_factor,
                     linewidth=3, label="Single Thread Speedup", color="#55E6C1")

        plt.legend(fontsize=20)

        # if "_15" in f:
        #     plt.yticks(list(range(0, 45, 5)))
        # else:
        plt.yticks(list(range(0, 5)))
        # plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
        #                     hspace=0, wspace=0)
        ax = plt.axes()
        ax.set_xscale('log')
        # ax.xaxis.set_major_locator(ticker.FixedLocator(x_range))
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + ".pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        # plot the eigen difference
        plt.plot(x_range, eigen_single_thread_factor,
                 linewidth=3, label="Single Thread Speedup Over Eigen", color="#d63031")
        plt.yticks(list(range(0, 30, 5)))
        # plt.legend()
        import math
        localMax = math.ceil(max(eigen_single_thread_factor))
        localMin = math.floor(min(eigen_single_thread_factor))
        plt.yticks(range(localMin, localMax+1,
                         math.ceil((localMax + 1 - localMin)/5)))
        ax = plt.axes()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.yaxis.set_major_locator(ticker.FixedLocator(list(range(0, 30, 5))))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_eigen.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        plt.plot(x_range, pre_comp_time,
                 linewidth=3, label="Pre Computation Time", color="#b8e994")
        # plt.legend()
        ax = plt.axes()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_pre_comp.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()
